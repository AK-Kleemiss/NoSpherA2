#include "pch.h"
#include "gaussian_atom.h"
#include "crystal_energies.h"
#include "scattering_factors.h"
#include "SALTED_predictor.h"
#include "basis_set.h"
#include "npy.h"
#include "spherical_density.h"
#include <random>

Multipoles Gaussian_Atom::electrical_moments(const int lmax) const
{
	Multipoles M{ lmax, vec((size_t)(lmax + 1) * (lmax + 1), 0.0) };
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
	err_checkf(static_cast<int>(coefs.size()) == tab.n_coef,
		"Coefficient count does not match the auxiliary basis: " + std::to_string(coefs.size())
		+ " coefficients for a basis of " + std::to_string(tab.n_coef), std::cout);
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

void Gaussian_Molecule::write_tshb(const std::filesystem::path& fn, const svec& labels, const ivec& atoms, const WFN& wave, const int ECP_mode, const bool ED) const
{
	err_checkf(labels.size() == atoms.size(), "Labels and atoms do not match for " + fn.string(), std::cout);
	std::ofstream f(fn, std::ios::binary);
	err_checkf(f.good(), "Could not open " + fn.string(), std::cout);
	auto u32 = [&f](const uint32_t v) { f.write(reinterpret_cast<const char*>(&v), 4); };
	auto str = [&f, &u32](const std::string& s) { u32(static_cast<uint32_t>(s.size())); f.write(s.data(), s.size()); };
	f.write("TSHB", 4);
	u32(3);
	//ED: smtbx applies the Mott-Bethe conversion of convert_to_ED per symmetry equivalent, with the Z stored per atom
	str(std::string("FRAME: cartesian\nUNITS: angstrom\nSOURCE: NoSpherA2") + (ED ? "\nED: true" : ""));
	u32(static_cast<uint32_t>(atoms.size()));
	const double b = constants::bohr2ang(1.0);
	for (int i = 0; i < atoms.size(); i++) {
		const int a = atoms[i];
		str(labels[i]);
		uint32_t n = 0;
		for (int s = tab.sh_start[a]; s < tab.sh_start[a + 1]; s++) n += tab.pr_start[s + 1] - tab.pr_start[s];
		u32(n);
		for (int s = tab.sh_start[a]; s < tab.sh_start[a + 1]; s++) {
			const int32_t l = tab.sh_l[s];
			//n exp(-a r^2) r^l in bohr is n b^-(l+3) exp(-(a/b^2) R^2) R^l in angstrom and e/A^3
			const double scale = std::pow(b, -(l + 3));
			for (int p = tab.pr_start[s]; p < tab.pr_start[s + 1]; p++) {
				const double alpha = tab.pr_exp[p] / (b * b);
				f.write(reinterpret_cast<const char*>(&l), 4);
				f.write(reinterpret_cast<const char*>(&alpha), 8);
				for (int m = 0; m < 2 * l + 1; m++) {
					const double c = coefs[tab.coef_off[s] + m] * tab.pr_norm[p] * scale;
					f.write(reinterpret_cast<const char*>(&c), 8);
				}
			}
		}
		const int n_ecp = wave.get_atom_ECP_electrons(a);
		std::vector<std::array<double, 3>> sl, ga;
		std::vector<std::array<double, 2>> ex;
		if (n_ecp != 0 && ECP_mode == 0)
			//mode 0 adds n e^(-k/16pi) with k = 4pi stl = 2pi s
			ex.push_back({ double(n_ecp), 1.0 / 8.0 });
		else if (n_ecp != 0) {
			const int Z = wave.get_atom_charge(a);
			//4 pi r rho = T r^N e^(-z r) in bohr is rho = T / (4 pi b^(N+2)) R^(N-1) e^(-(z/b) R) in angstrom
			for (const std::array<double, 3>& t : Thakkar(Z, ECP_mode).core_slater_terms(n_ecp))
				sl.push_back({ t[0] / (constants::FOUR_PI * std::pow(b, t[1] + 2)), t[1] - 1, t[2] / b });
			const Spherical_Gaussian_Density g(Z, ECP_mode);
			for (int k = 0; k < g.size(); k++)
				ga.push_back({ g.coefficient(k) / (b * b * b), g.exponent(k) / (b * b), 0.0 });
		}
		//version 2 records: l = 0 Gaussians c e^(-a R^2), Slater terms c R^n e^(-z R), form factor terms a e^(-b s)
		u32(static_cast<uint32_t>(ga.size()));
		for (const std::array<double, 3>& g : ga) f.write(reinterpret_cast<const char*>(g.data()), 16);
		u32(static_cast<uint32_t>(sl.size()));
		for (const std::array<double, 3>& t : sl) f.write(reinterpret_cast<const char*>(t.data()), 24);
		u32(static_cast<uint32_t>(ex.size()));
		for (const std::array<double, 2>& e : ex) f.write(reinterpret_cast<const char*>(e.data()), 16);
		//version 3: the nuclear charge Mott-Bethe subtracts, the one convert_to_ED uses
		const int32_t Z = wave.get_atom_charge(a);
		f.write(reinterpret_cast<const char*>(&Z), 4);
	}
	err_checkf(f.good(), "Writing " + fn.string() + " failed", std::cout);
}

Gaussian_Molecule Gaussian_Molecule::transformed(const vec2& R, const vec& t) const
{
	WFN a = aux;
	vec c = coefs;
	crystal_energies::transform(a, c, R, t);
	return Gaussian_Molecule(std::move(a), std::move(c));
}
