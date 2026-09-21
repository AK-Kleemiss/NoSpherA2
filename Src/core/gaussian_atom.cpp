#include "pch.h"
#include "gaussian_atom.h"
#include "crystal_energies.h"
#include "scattering_factors.h"
#include "SALTED_predictor.h"
#include "basis_set.h"
#include "npy.h"
#include <random>

Multipoles Gaussian_Atom::electrical_moments(const int lmax) const
{
	Multipoles M{ lmax, vec((lmax + 1) * (lmax + 1), 0.0) };
	M.q[0] = electrons();
	for (int s = 0; s < tab.n_sh; s++) {
		const int l = tab.sh_l[s];
		if (l == 0 || l > lmax) continue;
		const double I_l = tab.shell_radial_moment(s) * std::sqrt(constants::FOUR_PI / (2 * l + 1));
		for (int m = -l; m <= l; m++) M.q[l * l + l + m] += coefs[tab.coef_off[s] + l + m] * I_l;
	}
	return M;
}

cvec Gaussian_Atom::scattering_factors(const vec2& k_pt) const
{
	cvec sf(k_pt[0].size());
	for (size_t k = 0; k < sf.size(); k++) sf[k] = tab.fourier_atom(k_pt[0][k], k_pt[1][k], k_pt[2][k], coefs.data(), 0);
	return sf;
}

Interaction_Energy Gaussian_Atom::interaction(const Gaussian_Atom& B, const double repulsion_K, const int x_fun) const
{
	WFN a(e_origin::NOT_YET_DEFINED), b(e_origin::NOT_YET_DEFINED);
	a.push_back_atom(at);
	b.push_back_atom(B.at);
	return DensityFitting::interaction_energy(coefs, a, B.coefs, b, repulsion_K, x_fun);
}

namespace {
	//Fills aux with the auxiliary basis the coefficients belong to
	vec fitted_coefficients(const WFN& wavy, const std::filesystem::path& coef_file, WFN& aux, options& opt)
	{
		vec coef;
		if (!coef_file.empty()) {
			err_checkf(!opt.aux_basis.empty(), "No auxiliary basis set specified! Use -ri_fit BEFORE the interaction energy flag", std::cout);
			aux = generate_aux_wfn(wavy, opt.aux_basis);
			std::vector<unsigned long> shape; bool fortran_order;
			npy::LoadArrayFromNumpy(coef_file.string(), shape, fortran_order, coef);
		}
		else if (opt.SALTED) {
			err_checkf(!opt.salted_model_dir.empty(), "No SALTED model directory specified! Use -SALTED <model-dir> BEFORE the interaction energy flag", std::cout);
			SALTEDPredictor SP(wavy, opt);
			if (!SP.basis_set_loaded()) load_basis_into_WFN(SP.wavy, BasisSetLibrary::get_basis_set(SP.get_dfbasis_name()));
			coef = SP.gen_SALTED_densities();
			err_checkf(SP.wavy.get_ncen() == wavy.get_ncen(), "The SALTED model does not cover every atom of " + wavy.get_path().string(), std::cout);
			if (!opt.salted_charge_constraint)
				std::cout << "Hint: -salted_charge_constraint pins the electron count of the prediction; a missing fraction of an electron shifts the whole ESP by q/r" << std::endl;
			aux = SP.wavy;
			aux.set_origin(e_origin::NOT_YET_DEFINED);
		}
		else {
			err_checkf(!opt.aux_basis.empty(), "No auxiliary basis set specified! Use -ri_fit <basis> or -SALTED <model-dir> BEFORE the interaction energy flag", std::cout);
			err_checkf(wavy.get_nmo() > 0, "Only a wavefunction can be fitted; " + wavy.get_path().string() + " needs coefficients or -SALTED <model-dir>", std::cout);
			aux = generate_aux_wfn(wavy, opt.aux_basis);
			coef = DensityFitting::density_fit(wavy, aux, DensityFitting::config_from_options(opt));
		}
		return coef;
	}

	//128 random bits as 32 hex digits
	std::string new_uuid()
	{
		static std::mt19937_64 gen{ std::random_device{}() };
		std::ostringstream s;
		s << std::hex << std::setfill('0') << std::setw(16) << gen() << std::setw(16) << gen();
		return s.str();
	}
}

Gaussian_Molecule::Gaussian_Molecule(const WFN& wavy, options& opt, const std::filesystem::path& coef_file)
	: coefs(fitted_coefficients(wavy, coef_file, aux, opt)), tab(*aux.get_atoms_ptr()), id(new_uuid())
{
	slice();
}

Gaussian_Molecule::Gaussian_Molecule(WFN aux_basis, vec coefficients)
	: aux(std::move(aux_basis)), coefs(std::move(coefficients)), tab(*aux.get_atoms_ptr()), id(new_uuid())
{
	slice();
}

//One Gaussian_Atom per atom with its contiguous run of the coefficients
void Gaussian_Molecule::slice()
{
	err_checkf(static_cast<int>(coefs.size()) == tab.n_coef, "Coefficient count does not match the auxiliary basis", std::cout);
	const std::vector<::atom>& atoms = *aux.get_atoms_ptr();
	ats.reserve(tab.n_at);
	for (int a = 0; a < tab.n_at; a++) {
		const int first = tab.coef_off[tab.sh_start[a]];
		const int last = a + 1 < tab.n_at ? tab.coef_off[tab.sh_start[a + 1]] : tab.n_coef;
		ats.emplace_back(atoms[a], vec(coefs.begin() + first, coefs.begin() + last));
	}
}

double Gaussian_Molecule::electrons() const
{
	const vec pop = populations();
	return std::accumulate(pop.begin(), pop.end(), 0.0);
}

double Gaussian_Molecule::charge() const
{
	double Z = 0.0;
	for (int a = 0; a < tab.n_at; a++) Z += aux.get_atom_charge(a);
	return Z - electrons();
}

std::vector<Multipoles> Gaussian_Molecule::electrical_moments(const int lmax, const std::optional<DensityFitting::CHARGE_SCHEME> partition) const
{
	std::vector<Multipoles> M;
	if (!partition) {
		for (const Gaussian_Atom& a : ats) M.push_back(a.electrical_moments(lmax));
		return M;
	}
	const vec2 rows = DensityFitting::partition_multipole_rows(aux, tab, *partition, lmax, density_batch(*this));
	const int n_mom = (lmax + 1) * (lmax + 1);
	for (int a = 0; a < tab.n_at; a++) {
		Multipoles q{ lmax, vec(n_mom, 0.0) };
		for (int l = 0; l <= lmax; l++)
			for (int m = -l; m <= l; m++) {
				const int k = l * l + l + m;
				q.q[k] = std::inner_product(rows[a * n_mom + k].begin(), rows[a * n_mom + k].end(), coefs.begin(), 0.0) * std::sqrt(constants::FOUR_PI / (2 * l + 1));
			}
		M.push_back(q);
	}
	return M;
}

cvec2 Gaussian_Molecule::scattering_factors(const vec2& k_pt, ivec atoms, ProgressBar* progress) const
{
	if (atoms.empty()) {
		atoms.resize(tab.n_at);
		std::iota(atoms.begin(), atoms.end(), 0);
	}
	cvec2 sf;
	calc_SF_SALTED(k_pt, coefs, tab, atoms, sf, progress);
	return sf;
}

const Interaction_Energy& Gaussian_Molecule::interaction(const Gaussian_Molecule& B, const double repulsion_K, const int x_fun)
{
	auto it = energies.find(B.id);
	if (it == energies.end())
		it = energies.emplace(B.id, DensityFitting::interaction_energy(coefs, aux, B.coefs, B.aux, repulsion_K, x_fun)).first;
	return it->second;
}

Gaussian_Molecule Gaussian_Molecule::transformed(const vec2& R, const vec& t) const
{
	WFN a = aux;
	vec c = coefs;
	crystal_energies::transform(a, c, R, t);
	return Gaussian_Molecule(std::move(a), std::move(c));
}
