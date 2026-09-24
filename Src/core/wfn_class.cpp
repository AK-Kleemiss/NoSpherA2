#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "mo_class.h"
#include "cube.h"
#include "constants.h"
#include "fchk.h"
#include "basis_set.h"
#include "nos_math.h"
#include "libCintMain.h"
#include "integrator.h"
#include "cell.h"

#include "occ/OrbitalDefs.h"
long long int WFN::pre[9][5][5][9] = {};
long long int WFN::Afac_pre[9][5][9] = {};
constexpr void WFN::fill_pre()
{
	for (int j = 0; j < 9; j++)
		for (int l = 0; l < 5; l++)
			for (int m = 0; m < 5; m++)
			{
				int imax = std::min(j, l);
				int imin = std::max(0, j - m);
				for (int i = imin; i <= imax; i++)
					pre[j][l][m][i] = constants::ft[j] * constants::ft[l] / constants::ft[l - i] / constants::ft[i] * constants::ft[m] / constants::ft[m - j + i] / constants::ft[j - i];
			}
}

constexpr void WFN::fill_Afac_pre()
{
	for (int l = 0; l < 9; l++)
		for (int r = 0; r <= l / 2; r++)
			for (int s = 0; s <= (l - 2 * r) / 2; s++)
				Afac_pre[l][r][s] = constants::ft[r] * constants::ft[s] * constants::ft[l - 2 * r - 2 * s];
}

void WFN::reset()
{
	ncen = 0;
	nfunc = 0;
	nmo = 0;
	nex = 0;
	charge = 0;
	ECP_m = 0;
	multi = 0;
	origin = e_origin::NOT_YET_DEFINED;
	total_energy = 0.0;
	virial_ratio = 0.0;
	basis_set_name = " ";
	comment = "Test";
	path.clear();
	method.clear();
	MOs.clear();
	coef_primitive_major.clear();
	coef_primitive_major_valid = false;
	prim_exp_group.clear();
	group_exponent.clear();
	center_group_start.clear();
	center_min_exponent.clear();
	centers.clear();
	types.clear();
	exponents.clear();
	UT_DensityMatrix.clear();
	UT_SpinDensityMatrix.clear();
	DM = dMatrix2();
	DM_beta = dMatrix2();
	MO_sph = dMatrix2();
	basis_set = NULL;
	cub.clear();
	atoms.clear();
	fitted.clear();
	modified = false;
	d_f_switch = false;
	distance_switch = false;
	has_ECPs = false;
	isBohr = false;
	is_unrestricted = false;
};

WFN::WFN()
{
	reset();
	fill_pre();
	fill_Afac_pre();
};

WFN::WFN(e_origin given_origin)
{
	reset();
	origin = given_origin;
	fill_pre();
	fill_Afac_pre();
};

WFN::WFN(const std::filesystem::path &filename, const bool &debug)
{
	reset();
	fill_pre();
	fill_Afac_pre();
	read_known_wavefunction_format(filename, std::cout, debug);
};

WFN::WFN(const std::filesystem::path &filename, const int g_charge, const int g_mult, const bool &debug) {
	reset();
	charge = g_charge;
	multi = g_mult;
	fill_pre();
	fill_Afac_pre();
	read_known_wavefunction_format(filename, std::cout, debug);
};

constexpr unsigned int sum_subshells(unsigned int l, bool cartesian = true) {
	if (cartesian)
		return l * (l + 1) * (l + 2) / 6;
	return l * (2 * l + 1);
}

WFN::WFN(const occ::qm::Wavefunction &occ_WF, bool from_file) : WFN()
{
	// Should work for unrestricted now, if something seems weird this is the place to check
	using namespace Eigen;
	using occ::gto::num_subshells;
	set_origin(e_origin::OCC);
	basis_set_name = occ_WF.basis.name();
	has_ECPs = occ_WF.basis.have_ecps();
	charge = occ_WF.charge();
	const occ::Mat3N atom_positions = occ_WF.positions();
	for (long i = 0; i < occ_WF.atoms.size(); i++) {
		push_back_atom(constants::atnr2letter(occ_WF.atoms[i].atomic_number), atom_positions(0, i), atom_positions(1, i), atom_positions(2, i), static_cast<int>(occ_WF.nuclear_charges()(i)), {});
	}
	auto &shells = occ_WF.basis.shells();
	auto &mo = occ_WF.mo;
	const bool unrestricted = (mo.kind == occ::qm::Unrestricted);
	const int n_spin = unrestricted ? 2 : 1;

	for (const auto &shell : shells)
		err_checkf(shell.l <= 10, "Shell with angular momentum l=" + std::to_string(shell.l) + " is not supported by the OCC WFN constructor", std::cout);
	// ── Guard 2: this constructor assumes a spherical basis. ─────────────────
	// It applies A = MappedMatrices[l] (n_cart x n_sph) to transform spherical
	// MO coefficients to Cartesian. A Cartesian basis wavefunction would require
	// a different code path (basis_offset += n_cart, no A transformation).
	err_checkf(occ_WF.basis.is_pure(),
		"The OCC WFN constructor only supports spherical (pure) basis sets. "
		"Cartesian basis sets will produce incorrect MO coefficients.",
		std::cout);

	auto shell2atom = occ_WF.basis.shell_to_atom();
	vec con_coefs;
	VectorXd shellType(shells.size() + 1);
	nex = 0;
	for (const auto shell : shells)
	{
		int l = shell.l;
		const int nprim_shell = static_cast<int>(shell.num_primitives());
		err_checkf(static_cast<int>(shell.exponents.size()) == nprim_shell,
			"Shell exponents.size() (" + std::to_string(shell.exponents.size()) +
			") != num_primitives() (" + std::to_string(nprim_shell) +
			") – cannot build a consistent WFN from this OCC wavefunction.",
			std::cout);
		int n_cart = num_subshells(true, l);
		nex += n_cart * shell.exponents.size();
	}
	auto mo_go = occ::io::conversion::orb::to_gaussian_order(occ_WF.basis, occ_WF.mo);
	//OCC's reordering leaves the phases as libcint has them, while the sph2cart tables are
	//ORCA's, whose |m| = 3, 4, 7, 8 functions carry the opposite sign - the same difference
	//write_nbo corrects for. Without this a def2-TZVP wavefunction integrates to 149.58
	//electrons of 150 and its orbitals lose up to 3% of their norm.
	{
		int row = 0;
		for (const auto &shell : occ_WF.basis.shells()) {
			for (int spin = 0; spin < n_spin; spin++)
				for (int m = 3; m <= shell.l; m++)
					if (m % 4 == 3 || m % 4 == 0) {
						mo_go.C.row(spin * occ_WF.nbf + row + 2 * m - 1) *= -1.0;
						mo_go.C.row(spin * occ_WF.nbf + row + 2 * m) *= -1.0;
					}
			row += 2 * shell.l + 1;
		}
	}
	Vector<int, 10> d_orbital_corr{ 0, 1, 2, 6, 3, 4, 7, 8, 5, 9 };
	auto atom2shell = occ_WF.basis.atom_to_shell();

	unsigned int nprim;
	unsigned int n_cart;
	unsigned int sum_ncart;
	int atom;
	int last_atom = 0;
	int l;
	int k = 0;
	for (int i = 0; i < shells.size(); i++)
	{
		const auto &shell = shells[i];
		l = shell.l;
		shellType(i) = l + 1;
		nprim = shell.num_primitives();
		n_cart = num_subshells(true, l);
		sum_ncart = sum_subshells(l);
		occ::Vec occ_exp = shell.exponents.replicate(n_cart, 1);
		insert_into_exponents(occ_vec_span(occ_exp));
		atom = shell2atom[i];
		if (i != 0) {
			last_atom = shell2atom[i - 1];
		}
		if (last_atom != atom) {
			k = 0;
		}
		// insert_into_centers(std::views::repeat(atom+1, n_cart*nprim));
		for (int j = 0; j < shell.exponents.size(); j++) {
			//l + 1, the convention of every wavefunction basis in this program: WFN's own
			//get_shell_type() hands this very field out and the primitive counters switch over it
			//as 1 = s, 2 = p, ..., so an s shell stored as 0 matches no case and they count
			//nothing. Only the aux bases, which carry origin NOT_YET_DEFINED, store l itself.
			push_back_atom_basis_set(atom, shell.exponents(j), shell.contraction_coefficients(j), shell.l + 1, k);
		}
		k++;
		auto repeated = std::views::iota(0u, n_cart * nprim) | std::views::transform([&](auto) { return atom + 1; });
		insert_into_centers(repeated);

		VectorXi typesVec = ArrayXi::LinSpaced(n_cart, sum_ncart + 1, sum_ncart + n_cart)
			.matrix().transpose().replicate(nprim, 1).reshaped();
		insert_into_types(typesVec);
		// if (l==3) {
		//     Vector<int, 10> reordered = typesVec(d_orbital_corr);
		//     insert_into_types(typesVec);
		// }
		// else
		//     insert_into_types(typesVec);
	}

	// This could also be calculated in the loop above and I tried it, but I noticed when looking at the data that
	// a vector was inserted into coefficients for each MO. I think that probably that would result in more work for the
	// CPU due to cache misses, but I didn't benchmark it.
	unsigned int chunk_size;
	unsigned int n_sph;
	unsigned int n_prim;
	double occ;
	double *coeffs_ptr;
	unsigned int basis_offset;
	unsigned int write_cursor;
	double p;
	double scalar;

	// confac = pow(8 * pow(exp[exp_run], 3) / constants::PI3, 0.25);
	for (int spin = 0; spin < n_spin; spin++) {
		const int row_offset = unrestricted ? spin * occ_WF.nbf : 0;
		const int nocc = unrestricted ? (spin == 0 ? occ_WF.n_alpha() : occ_WF.n_beta()) : occ_WF.n_alpha();
		for (int n = 0; n < occ_WF.nbf; n++)
		{
			basis_offset = 0;
			write_cursor = 0;
			occ = (n < nocc) ? (unrestricted ? 1.0 : 2.0) : 0.0;
			int energy_index = unrestricted ? spin * occ_WF.nbf + n : n;
			push_back_MO(energy_index + 1, occ, mo.energies[energy_index], spin);
			MO& currentMO = MOs.back();
			currentMO.assign_coefficients_size(nex);
			coeffs_ptr = const_cast<double*>(currentMO.get_coefficient_ptr());
			for (const auto& shell : shells) {
				l = shell.l;
				p = (2.0 * l + 3.0) / 4.0;
				scalar = std::pow(2.0, 0.5 * l) / std::pow(constants::PI3, 0.25) / std::sqrt(constants::sph2cart_norm2[l]);
				n_sph = constants::n_spher(l);
				n_prim = shell.num_primitives();
				n_cart = constants::n_cart(l);
				//occ's Gaussian order is x,y,z for p and m = 0,+1,-1,... above, the table columns are m ordered
				MatrixXd A = Map<const Matrix<double, Dynamic, Dynamic, RowMajor>>(constants::sph2cart(l), n_cart, n_sph);
				if (l == 1) A.setIdentity();
				chunk_size = n_cart * n_prim;
				const auto& exp_arr = shell.exponents.array();
				// normalizing the contraction_coeffs
				ArrayXd CC = shell.u_coefficients.array();
				if (!from_file)
					CC = VectorXd::NullaryExpr(CC.size(), [shell](const Index i) {
					return shell.coeff_normalized(0, i);
						});
				const VectorXd& contraction_coeffs = scalar * (2 * exp_arr).pow(p) * CC;
				const auto& MOc = mo_go.C.block(row_offset + basis_offset, n, n_sph, 1);
				Map<MatrixXd> dest_block(coeffs_ptr + write_cursor, n_prim, n_cart);
				// |C>(A|MOc>)^T Has dimensions of (num_subshells(l), m). m: contraction \in R^{(m,1)})
				dest_block = contraction_coeffs * (A * MOc).transpose();

				write_cursor += chunk_size;
				basis_offset += n_sph;
			}
		}
	}
	set_exp_cutoff();

	d_f_switch = false;
	int nmo_ = occ_WF.mo.D.cols();
	DM = dMatrix2(nmo_, nmo_);
	if (occ_WF.mo.kind == occ::qm::Restricted) {
		for (int i = 0; i < nmo_; i++) {
			DM(i, i) = 2 * occ_WF.mo.D(i, i);
			for (int j = i + 1; j < nmo_; j++) {
				DM(i, j) = 2 * occ_WF.mo.D(i, j);
				DM(j, i) = DM(i, j);
			}
		}
	}
	else if (occ_WF.mo.kind == occ::qm::Unrestricted) {
		for (int i = 0; i < nmo_; i++) {
			DM(i, i) = occ_WF.mo.D(i, i) + occ_WF.mo.D(i + nmo_, i);
			for (int j = i + 1; j < nmo_; j++) {
				DM(i, j) = occ_WF.mo.D(i, j) + occ_WF.mo.D(i + nmo_, j);
				DM(j, i) = DM(i, j);
			}
		}
		is_unrestricted = true;
	}
	MO_sph = dMatrix2(occ_WF.mo.C.rows(), occ_WF.mo.C.cols());
	for (long i = 0; i < occ_WF.mo.C.rows(); i++)
		for (long j = 0; j < occ_WF.mo.C.cols(); j++)
			MO_sph(i, j) = occ_WF.mo.C(i, j);
}

// Hokus pokus freestyle modus, hopefully converts a WFN object to the occ wavefunction
// Fully vibe coded without application because WFNs don't include virtual orbitals
void WFN::wfn_to_occ_wavefunction(occ::qm::Wavefunction& occ_wf)
{
	using occ::gto::num_subshells;
	using AOBasisT = std::remove_cvref_t<decltype(occ_wf.basis)>;
	using AtomListT = AOBasisT::AtomList;
	using ShellListT = AOBasisT::ShellList;
	const int n_atoms = get_ncen();
	AtomListT occ_atoms(n_atoms);
	for (int i = 0; i < n_atoms; i++) {
		const atom& at = get_atom(i);
		occ_atoms[i].atomic_number = at.get_charge();
		occ_atoms[i].x = at.get_pos()[0];
		occ_atoms[i].y = at.get_pos()[1];
		occ_atoms[i].z = at.get_pos()[2];
	}
	struct RebuiltShell {
		int atom, l, n_cart, n_prim, nex_start;
		occ::Vec exponents;
	};
	std::vector<RebuiltShell> rebuilt_shells;
	{
		int pos = 0;
		const int nex_total = get_nex();
		while (pos < nex_total) {
			const int type0 = get_type(pos);
			const int atom0 = get_center(pos);
			int l = 0;
			while (type0 > static_cast<int>(sum_subshells(l)) + static_cast<int>(num_subshells(true, l)))
				l++;
			const int n_cart = num_subshells(true, l);
			int n_prim = 0;
			int scan = pos;
			while (scan < nex_total &&
				get_type(scan) == type0 &&
				get_center(scan) == atom0) {
				n_prim++;
				scan++;
			}

			vec exp_v(n_prim);
			for (int i = 0; i < n_prim; i++)
				exp_v[i] = get_exponent(pos + i);

			RebuiltShell rs;
			rs.atom = atom0 - 1;
			rs.l = l;
			rs.n_cart = n_cart;
			rs.n_prim = n_prim;
			rs.nex_start = pos;
			rs.exponents = Eigen::Map<occ::Vec>(exp_v.data(), n_prim);
			rebuilt_shells.push_back(rs);
			pos += n_prim * n_cart;
		}
	}

	const auto& mo0 = get_MO(0);
	const double* mo0_coeffs = mo0.get_coefficient_ptr();
	ShellListT occ_shells;
	occ_shells.reserve(rebuilt_shells.size());
	std::vector<occ::Vec> actual_contraction_coeffs(rebuilt_shells.size());

	for (size_t k = 0; k < rebuilt_shells.size(); k++) {
		const auto& rs = rebuilt_shells[k];
		const int l = rs.l, n_prim = rs.n_prim, n_cart = rs.n_cart;
		int write_cursor = 0;
		for (size_t k2 = 0; k2 < k; k2++)
			write_cursor += rebuilt_shells[k2].n_prim * rebuilt_shells[k2].n_cart;
		Eigen::Map<const Eigen::MatrixXd> dest_block0(mo0_coeffs + write_cursor, n_prim, n_cart);
		int j_ref = 0; double best = 0;
		for (int j = 0; j < n_cart; j++) {
			double m = dest_block0.col(j).cwiseAbs().maxCoeff();
			if (m > best) { best = m; j_ref = j; }
		}
		occ::Vec ratio = dest_block0.col(j_ref);
		const double scalar = std::pow(2.0, 0.5 * l) / std::pow(constants::PI3, 0.25) / std::sqrt(constants::sph2cart_norm2[l]);
		const double p = (2.0 * l + 3.0) / 4.0;
		vec cc_input(n_prim), expo(n_prim);
		for (int i = 0; i < n_prim; i++) {
			expo[i] = rs.exponents(i);
			cc_input[i] = ratio(i) / (scalar * std::pow(2.0 * expo[i], p));
		}
		// ratio carries the sign of the reference MO coefficient; the AO sign is a convention, so make the
		// dominant contraction coefficient positive (what basis set files usually carry)
		int i_ref = 0;
		for (int i = 1; i < n_prim; i++)
			if (std::abs(cc_input[i]) > std::abs(cc_input[i_ref])) i_ref = i;
		if (cc_input[i_ref] < 0)
			for (double& c : cc_input) c = -c;
		std::array<double, 3> pos{ occ_atoms[rs.atom].x,
									occ_atoms[rs.atom].y,
									occ_atoms[rs.atom].z };
		occ::gto::Shell sh(l, expo, { cc_input }, pos);
		sh.kind = occ::gto::Shell::Kind::Spherical;
		// occ's loaders normalise every shell they build; coeff_normalized below assumes it
		sh.incorporate_shell_norm();
		occ::Vec true_cc(n_prim);
		for (int i = 0; i < n_prim; i++)
			true_cc(i) = scalar * std::pow(2.0 * expo[i], p) * sh.coeff_normalized(0, i);

		actual_contraction_coeffs[k] = true_cc;
		occ_shells.push_back(std::move(sh));
	}
	occ_wf.basis = AOBasisT(occ_atoms, occ_shells);
	occ_wf.atoms = occ_atoms;
	occ_wf.basis.set_pure(true);
	const int nbf_sph = static_cast<int>(occ_wf.basis.nbf());
	const int n_mo_total = get_nmo();
	const bool unrestricted = get_is_unrestricted();
	const int n_spin = unrestricted ? 2 : 1;
	occ::Mat C_gaussian_order(nbf_sph * n_spin, n_mo_total / n_spin);
	occ::Vec energies(n_mo_total);
	occ::Vec occupations(n_mo_total);
	int n_alpha = 0, n_beta = 0;
	for (int spin = 0; spin < n_spin; spin++) {
		for (int n = 0; n < n_mo_total / n_spin; n++) {
			const int mo_index = unrestricted ? spin * (n_mo_total / n_spin) + n : n;
			const auto& mo = get_MO(mo_index);
			energies(mo_index) = mo.get_energy();
			occupations(mo_index) = mo.get_occ();
			if (unrestricted) {
				if (spin == 0 && mo.get_occ() > 0.5) n_alpha++;
				if (spin == 1 && mo.get_occ() > 0.5) n_beta++;
			}
			else {
				if (mo.get_occ() > 1.5) { n_alpha++; n_beta++; }
				else if (mo.get_occ() > 0.5) n_alpha++;
			}
			const double* coeffs_ptr = mo.get_coefficient_ptr();
			int write_cursor = 0, sph_offset = spin * nbf_sph;
			for (size_t k = 0; k < rebuilt_shells.size(); k++) {
				const auto& rs = rebuilt_shells[k];
				Eigen::Map<const Eigen::MatrixXd> dest_block(coeffs_ptr + write_cursor, rs.n_prim, rs.n_cart);
				const occ::Vec& cc = actual_contraction_coeffs[k];
				occ::Vec cart_mo_coeffs =
					(dest_block.transpose() * cc) / cc.squaredNorm();
				Eigen::MatrixXd A = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(constants::sph2cart(rs.l), rs.n_cart, constants::n_spher(rs.l));
				if (rs.l == 1) A.setIdentity();
				occ::Vec MOc = A.completeOrthogonalDecomposition().pseudoInverse() * cart_mo_coeffs;
				C_gaussian_order.block(sph_offset, n, MOc.size(), 1) = MOc;
				sph_offset += A.cols();
				write_cursor += rs.n_prim * rs.n_cart;
			}
		}
	}
	occ::qm::MolecularOrbitals mo_go;
	mo_go.kind = unrestricted ? occ::qm::Unrestricted : occ::qm::Restricted;
	mo_go.n_alpha = n_alpha;
	mo_go.n_beta = n_beta;
	mo_go.n_ao = nbf_sph;
	mo_go.C = C_gaussian_order;
	mo_go.energies = energies;
	mo_go.occupation = occupations;
	occ_wf.num_electrons = n_alpha + n_beta;
	occ_wf.mo = occ::io::conversion::orb::from_gaussian_order(occ_wf.basis, mo_go);
	occ_wf.mo.update_occupied_orbitals();
	occ_wf.mo.update_density_matrix();
	occ_wf.nbf = occ_wf.basis.nbf();
}

bool WFN::push_back_atom(const std::string &label, const double &x, const double &y, const double &z, const int &_charge, const atomID &ID)
{
	ncen++;
	if (_charge >= 1)
		atoms.emplace_back(label, ID, ncen, x, y, z, _charge);
	else
	{
		atoms.emplace_back();
		return false;
	}
	return true;
};

bool WFN::push_back_atom(const atom &given)
{
	ncen++;
	atoms.push_back(given);
	return true;
};

bool WFN::erase_atom(const int &nr)
{
	err_checkf(nr < ncen, "unreasonable atom number", std::cout);
	// By index: removeElement matches by value, and the value it is given
	// is a reference into the vector std::remove then shifts under it.
	atoms.erase(atoms.begin() + nr);
	ncen--;
	return true;
};

bool WFN::push_back_MO(const int &nr, const double &occ, const double &ener)
{
	nmo++;
	err_checkf(nr <= nmo, "unreasonable MO number", std::cout);
	MOs.push_back(MO(nr, occ, ener));
	invalidate_coef_cache();
	return true;
};

bool WFN::push_back_MO(const int &nr, const double &occ, const double &ener, const int &oper)
{
	nmo++;
	MOs.push_back(MO(nr, occ, ener, oper));
	invalidate_coef_cache();
	return true;
};

bool WFN::push_back_MO(const MO &given)
{
	nmo++;
	MOs.push_back(given);
	invalidate_coef_cache();
	return true;
};

void WFN::push_back_MO_coef(const int &nr, const double &value)
{
	err_checkf(nr < nmo, "not enough MOs", std::cout);
	MOs[nr].push_back_coef(value);
	invalidate_coef_cache();
};

void WFN::assign_MO_coefs(const int &nr, vec &values)
{
	err_checkf(nr < nmo, "not enough MOs", std::cout);
	MOs[nr].assign_coefs(values);
	invalidate_coef_cache();
};

const double &WFN::get_MO_energy(const int &mo) const
{
	err_checkf(mo < nmo, "not enough MOs", std::cout);
	return MOs[mo].get_energy();
}

const void WFN::clear_MOs()
{
	MOs.clear();
	MOs.shrink_to_fit();
	nmo = 0;
	invalidate_coef_cache();
}

bool WFN::push_back_center(const int &cent)
{
	if (cent <= ncen && cent > 0)
		centers.push_back(cent);
	else
		return false;
	return true;
};

bool WFN::erase_center(const int &g_nr)
{
	centers.erase(centers.begin() + g_nr - 1);
	return true;
};

const std::string WFN::get_centers(const bool &bohr) const
{
	std::string temp;
	for (int i = 0; i < ncen; i++)
	{
		temp.append(atoms[i].get_label());
		temp.append(" ");
		if (bohr)
			temp.append(std::to_string(get_atom_coordinate(i, 0)));
		else
			temp.append(std::to_string(constants::bohr2ang(get_atom_coordinate(i, 0))));
		temp.append(" ");
		if (bohr)
			temp.append(std::to_string(get_atom_coordinate(i, 1)));
		else
			temp.append(std::to_string(constants::bohr2ang(get_atom_coordinate(i, 1))));
		temp.append(" ");
		if (bohr)
			temp.append(std::to_string(get_atom_coordinate(i, 2)));
		else
			temp.append(std::to_string(constants::bohr2ang(get_atom_coordinate(i, 2))));
		temp.append("\n");
	}
	return temp;
};

const void WFN::list_centers() const
{
	for (int i = 0; i < ncen; i++)
	{
		std::cout << atoms[i].get_nr() << " " << atoms[i].get_label() << " "
			<< get_atom_coordinate(i, 0) << " " << get_atom_coordinate(i, 1) << " "
			<< get_atom_coordinate(i, 2) << " " << get_atom_charge(i) << "\n";
	}
};

const MO &WFN::get_MO(const int &n) const
{
	if (n < nmo)
		return MOs[n];
	else
	{
		err_not_impl_f("Wrong MO number", std::cout);
		return MOs[0];
	}
}
const int WFN::get_MO_op_count(const int &op) const
{
	int count = 0;
#pragma omp parallel for reduction(+ : count)
	for (int i = 0; i < nmo; i++)
		if (MOs[i].get_op() == op)
			count++;
	return count;
};

void WFN::delete_MO(const int &nr)
{
	err_checkf(nr < nmo, "not enough MOs", std::cout);
	MOs.erase(MOs.begin() + nr);
	nmo--;
	invalidate_coef_cache();
};

bool WFN::push_back_type(const int &type)
{
	types.push_back(type);
	return true;
};

bool WFN::erase_type(const int &nr)
{
	err_checkf(nr >= 1, "Wrong type to erase!", std::cout);
	types.erase(types.begin() + nr - 1);
	return true;
};

bool WFN::push_back_exponent(const double &e)
{
	exponents.push_back(e);
	return true;
};

bool WFN::erase_exponent(const int &nr)
{
	if (nr < 1)
		return false;
	exponents.erase(exponents.begin() + (nr - 1));
	return true;
};

bool WFN::remove_primitive(const int &nr)
{
	nex--;
	invalidate_coef_cache();
	if (erase_center(nr) && erase_exponent(nr) && erase_type(nr))
	{
		for (int n = 0; n < nmo; n++)
			MOs[n].erase_coef(nr, nex);
		return true;
	}
	else
		return false;
};

bool WFN::add_primitive(const int &cent, const int &type, const double &e, double *values)
{
	if (!push_back_center(cent) || !push_back_type(type) || !push_back_exponent(e))
		return false;
	nex++;
	invalidate_coef_cache();
	for (int n = 0; n < nmo; n++)
		MOs[n].push_back_coef(values[n]);
	return true;
};

void WFN::change_type(const int &nr)
{
	err_checkf(nr < nex, "Wrong input", std::cout);
	bool end = false;
	while (!end)
	{
		std::cout << "Please enter the new type you want to assign: ";
		int new_type = 0;
		std::cin >> new_type;
		if (new_type > ncen || new_type < 0)
		{
			std::cout << "Sorry, wrong input, try again!\n";
			continue;
		}
		types[nr - 1] = new_type;
		end = true;
		set_modified();
	}
};

void WFN::change_exponent(const int &nr)
{
	err_checkf(nr < nex, "Wrong input", std::cout);
	bool end = false;
	while (!end)
	{
		std::cout << "Please enter the new exponent you want to assign: ";
		int new_exp = 0;
		std::cin >> new_exp;
		if (new_exp > ncen || new_exp < 0)
		{
			std::cout << "Sorry, wrong input, try again!\n";
			continue;
		}
		exponents[nr - 1] = new_exp;
		invalidate_coef_cache();
		end = true;
		set_modified();
	}
};

void WFN::change_center(const int &nr)
{
	bool end = false;
	while (!end)
	{
		std::cout << "Please enter the new center you want to assign: ";
		int new_center = 0;
		std::cin >> new_center;
		if (new_center > ncen || new_center < 0)
		{
			std::cout << "Sorry, wrong input, try again!\n";
			continue;
		}
		centers[nr - 1] = new_center;
		invalidate_coef_cache();
		end = true;
		set_modified();
	}
};

bool WFN::set_MO_coef(const int &nr_mo, const int &nr_primitive, const double &value)
{
	err_checkf(nr_mo < MOs.size(), "MO doesn't exist!", std::cout);
	invalidate_coef_cache();
	return MOs[nr_mo].set_coefficient(nr_primitive, value);
};

const void WFN::list_primitives() const
{
	for (int i = 0; i < nex; i++)
	{
		std::cout << i << " center: " << centers[i] << " type: " << types[i] << " exponent: " << exponents[i] << "\n";
	}
};

bool WFN::remove_center(const int &nr)
{
	if (nr < 1 || nr > ncen)
		return false;
	// nr is the 1-based centre: drop its primitives (backwards, remove_primitive is 1-based),
	// then the atom, and renumber the centres behind it
	for (int i = nex - 1; i >= 0; i--)
		if (centers[i] == nr && !remove_primitive(i + 1))
			return false;
	erase_atom(nr - 1);
	for (int i = 0; i < nex; i++)
		if (centers[i] > nr)
			centers[i]--;
	set_modified();
	return true;
}

bool WFN::add_exp(const int &cent, const int &type, const double &e)
{
	if (!push_back_center(cent) || !push_back_type(type) || !push_back_exponent(e))
		return false;
	nex++;
	return true;
};

const double &WFN::get_MO_coef(const int &nr_mo, const int &nr_primitive) const
{
	err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
	return MOs[nr_mo].get_coefficient(nr_primitive);
};

const double &WFN::get_MO_coef_f(const int &nr_mo, const int &nr_primitive) const
{
	err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
	return MOs[nr_mo].get_coefficient_f(nr_primitive);
};

const double *WFN::get_MO_coef_ptr(const int &nr_mo)
{
	err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
	return MOs[nr_mo].get_coefficient_ptr();
};

const int WFN::get_MO_primitive_count(const int &nr_mo) const
{
	err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
	return MOs[nr_mo].get_primitive_count();
};

std::vector<asym_atom> WFN::extract_xyz(const std::string& unit) {
	std::vector<asym_atom> atoms;
	const d3 dummy_coords({0, 0, 0});
	for (int i = 0; i < ncen; i++) {
		atom temp_atom = get_atom(i);
		d3 temp_coords = temp_atom.get_pos();
		if (isBohr == true && unit == "angstrom") {
			temp_coords[0] = constants::bohr2ang(temp_coords[0]);
			temp_coords[1] = constants::bohr2ang(temp_coords[1]);
			temp_coords[2] = constants::bohr2ang(temp_coords[2]);
		}
		else if (isBohr == false && unit == "bohr") {
			temp_coords[0] = constants::ang2bohr(temp_coords[0]);
			temp_coords[1] = constants::ang2bohr(temp_coords[1]);
			temp_coords[2] = constants::ang2bohr(temp_coords[2]);
		}
		atoms.emplace_back(asym_atom(temp_atom.get_label(), temp_atom.get_charge(), temp_coords, dummy_coords, 1.0, cdouble(0.0, 0.0)));
	}
	return atoms;
}
const double WFN::get_maximum_MO_coefficient(bool occu) const {
	double max_coef = 0.0;
	for (int i = 0; i < nmo; i++) {
		if (occu && MOs[i].get_occ() == 0.0)
			continue;
		for (int j = 0; j < nex; j++) {
			if (std::abs(MOs[i].get_coefficients()[j]) > max_coef) {
				max_coef = std::abs(MOs[i].get_coefficients()[j]);
			}
		}
	}
	return max_coef;
};

const double WFN::get_atom_coordinate(const unsigned int &nr, const unsigned int &axis) const
{
	err_checkf(!((int)nr >= ncen || axis > 2), "This input is invalid for get_atom_coordinate!", std::cout);
	return atoms[nr].get_coordinate(axis);
};

const d3 WFN::get_atom_pos(const unsigned int &nr) const
{
	err_checkf(!(nr >= ncen), "This input is invalid for get_atom_coordinate!", std::cout);
	return atoms[nr].get_pos();
};

void WFN::print_primitive(const int &nr) const
{
	std::cout << "center assignement: " << centers[nr] << " type: " << types[nr]
		<< " exponent: " << exponents[nr] << std::endl
		<< "MO coefficients:";
	for (int i = 0; i < nmo; i++)
	{
		std::cout << MOs[i].get_coefficient(nr) << "   ";
		if (i % 5 == 0)
			std::cout << "\n";
	}
};

const int WFN::get_nmo(const bool &only_occ) const
{
	if (!only_occ)
		return nmo;
	else
	{
		int count = 0;
#pragma omp parallel for reduction(+ : count)
		for (int i = 0; i < MOs.size(); i++)
		{
			if (MOs[i].get_occ() != 0.0)
			{
				count++;
			}
		}
		return count;
	}
};

const unsigned int WFN::get_nr_electrons() const
{
	unsigned int count = 0;
	for (int i = 0; i < ncen; i++)
		count += get_atom_charge(i);
	count -= charge;
	return count;
};

const unsigned int WFN::get_nr_ECP_electrons() const
{
	unsigned int count = 0;
	for (int i = 0; i < ncen; i++)
		count += atoms[i].get_ECP_electrons();
	return count;
}

double WFN::count_nr_electrons(void) const
{
	double count = 0;
	for (int i = 0; i < nmo; i++)
		count += MOs[i].get_occ();
	return count;
};

double WFN::count_alpha_electrons(void) const
{
	double count = 0;
	for (int i = 0; i < nmo; i++)
		if (MOs[i].get_spin() == 0)
			count += MOs[i].get_occ();
	return count;
}

double WFN::count_beta_electrons(void) const
{
	double count = 0;
	for (int i = 0; i < nmo; i++)
		if (MOs[i].get_spin() == 1)
			count += MOs[i].get_occ();
	return count;
}

const double WFN::get_atom_basis_set_exponent(const int &nr_atom, const int &nr_prim) const
{
	if (nr_atom < ncen && nr_atom >= 0 && (int)atoms[nr_atom].get_basis_set_size() > nr_prim && nr_prim >= 0)
		return atoms[nr_atom].get_basis_set_exponent(nr_prim);
	else
		return -1;
};

const double WFN::get_atom_basis_set_coefficient(const int &nr_atom, const int &nr_prim) const
{
	if (nr_atom < ncen && nr_atom >= 0 && (int)atoms[nr_atom].get_basis_set_size() > nr_prim && nr_prim >= 0)
		return atoms[nr_atom].get_basis_set_coefficient(nr_prim);
	else
		return -1;
};

bool WFN::change_atom_basis_set_exponent(const int &nr_atom, const int &nr_prim, const double &value)
{
	if (nr_atom < ncen && nr_atom >= 0 && (int)atoms[nr_atom].get_basis_set_size() > nr_prim && nr_prim >= 0)
	{
		atoms[nr_atom].set_basis_set_exponent(nr_prim, value);
		set_modified();
		return true;
	}
	else
		return false;
};

bool WFN::change_atom_basis_set_coefficient(const int &nr_atom, const int &nr_prim, const double &value)
{
	err_checkf(nr_atom < ncen && nr_atom >= 0 && (int)atoms[nr_atom].get_basis_set_size() > nr_prim && nr_prim >= 0, "Wrong input!", std::cout);
	atoms[nr_atom].set_basis_set_coefficient(nr_prim, value);
	set_modified();
	return true;
};

const int WFN::get_atom_primitive_count(const int &nr) const
{
	if (nr < ncen && nr >= 0)
		return (int)atoms[nr].get_basis_set_size();
	else
		return -1;
};

const int WFN::get_basis_set_shell(const unsigned int &nr_atom, const unsigned int &nr_prim) const
{
	if ((int)nr_atom < ncen && atoms[nr_atom].get_basis_set_size() > (int)nr_prim)
	{
		return atoms[nr_atom].get_basis_set_shell(nr_prim);
	}
	else
		return -1;
};

const int WFN::get_atom_primitive_type(const int &nr_atom, const int &nr_prim) const
{
	if (nr_atom < atoms.size() && nr_atom >= 0 && nr_prim < (int)atoms[nr_atom].get_basis_set_size() && nr_prim >= 0)
		return atoms[nr_atom].get_basis_set_type(nr_prim);
	else
		return -1;
};

const int WFN::get_atom_shell_count(const unsigned int &nr) const
{
	if ((int)nr < ncen)
		return (int)atoms[nr].get_shellcount_size();
	else
		return -1;
};

const int WFN::get_atom_shell_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
	if ((int)nr_atom < ncen && (int)nr_shell < atoms[nr_atom].get_shellcount_size())
		return atoms[nr_atom].get_shellcount(nr_shell);
	else
		return -1;
};

const int WFN::get_shell_type(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
	if (static_cast<int>(nr_atom) < ncen && nr_shell < atoms[nr_atom].get_shellcount_size())
	{
		int primitive_counter = 0;
		while (atoms[nr_atom].get_basis_set_shell(primitive_counter) != nr_shell)
			primitive_counter++;
		return atoms[nr_atom].get_basis_set_type(primitive_counter);
	}
	else
		return -1;
};

const int WFN::get_shell_center(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
	if (static_cast<int>(nr_atom) < ncen && nr_shell < atoms[nr_atom].get_shellcount_size())
		return centers[get_shell_start_in_primitives(nr_atom, nr_shell)];
	else
		return -1;
};

const int WFN::get_shell_start(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
	if (static_cast<int>(nr_atom) < ncen && nr_shell < atoms[nr_atom].get_shellcount_size())
	{
		int primitive_counter = 0;
#pragma loop(no_vector)
		for (int s = 0; s < static_cast<int>(nr_shell); s++)
			primitive_counter += atoms[nr_atom].get_shellcount(s);
		return primitive_counter;
	}
	else
		return -1;
};

const int WFN::get_shell_start_in_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
	if (static_cast<int>(nr_atom) < ncen && nr_shell < atoms[nr_atom].get_shellcount_size())
	{
		//The cartesian components of a shell are the gap between consecutive WFN type blocks:
		//1, 3, 6, 10, 15, ... A switch over s/p/d/f used to stand here and added nothing at all for
		//g and above, so every primitive index behind the first g shell of a wavefunction was short
		//by 15 per g shell - which is how Fe.gbw's atom 2 asked for its s shell and was handed a g
		//primitive 540 places later, then wrote past the end of a 1-component buffer.
		const auto cart_components = [](const int shell_type) {
			return (shell_type >= 1 && shell_type < static_cast<int>(std::size(constants::first_type)))
				? constants::first_type[shell_type] - constants::first_type[shell_type - 1]
				: 0;
		};
		int primitive_counter = 0;
		for (unsigned int a = 0; a < nr_atom; a++)
			for (unsigned int s = 0; s < atoms[a].get_shellcount_size(); s++)
				primitive_counter += cart_components(get_shell_type(a, s)) * atoms[a].get_shellcount(s);
		for (unsigned int s = 0; s < nr_shell; s++)
			primitive_counter += cart_components(get_shell_type(nr_atom, s)) * atoms[nr_atom].get_shellcount(s);
		return primitive_counter;
	}
	else
		return -1;
};

const int WFN::get_shell_end(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
	if (static_cast<int>(nr_atom) < ncen && nr_shell < atoms[nr_atom].get_shellcount_size())
	{
		if (nr_shell == atoms[nr_atom].get_shellcount_size() - 1)
			return (int)atoms[nr_atom].get_basis_set_size() - 1;
		int primitive_counter = 0;
		while (atoms[nr_atom].get_basis_set_shell(primitive_counter) != (nr_shell + 1))
			primitive_counter++;
		return primitive_counter - 1;
	}
	else
		return -1;
};

const std::string WFN::get_atom_label(const unsigned int &nr) const
{
	std::string error_return{ '?' };
	if (nr < static_cast<unsigned int>(ncen))
		return atoms[nr].get_label();
	else
		return error_return;
};

const int WFN::get_nr_basis_set_loaded() const
{
	int count = 0;
	for (int a = 0; a < ncen; a++)
		if (atoms[a].get_basis_set_loaded())
			count++;
	return count;
};

const bool WFN::get_atom_basis_set_loaded(const int &nr) const
{
	if (nr < ncen && nr >= 0)
		return atoms[nr].get_basis_set_loaded();
	else
	{
		std::cout << "invalid atom choice in atom_basis_set_loaded!" << std::endl;
		return false;
	}
};

const int WFN::get_atom_charge(const int &nr) const
{
	if (nr < ncen && nr >= 0)
		return atoms[nr].get_charge();
	else
	{
		std::cout << "invalid atom choice in atom_basis_set_loaded!" << std::endl;
		return -1;
	}
};

void WFN::push_back_DM(const double &value)
{
	UT_DensityMatrix.push_back(value);
};

void WFN::resize_DM(const int &size, const double &value)
{
	UT_DensityMatrix.resize(size, value);
};

const double WFN::get_DM(const int &nr) const
{
	if (nr >= 0 && nr < UT_DensityMatrix.size())
		return UT_DensityMatrix[nr];
	else
	{
		std::cout << "Requested nr out of range! Size: " << UT_DensityMatrix.size() << " nr: " << nr << std::endl;
		return -1;
	}
};

bool WFN::set_DM(const int &nr, const double &value)
{
	if (nr >= 0 && nr < UT_DensityMatrix.size())
	{
		UT_DensityMatrix[nr] = value;
		return true;
	}
	else
	{
		std::cout << "invalid arguments for set_DM! Input was: " << nr << ";" << value << std::endl;
		return false;
	}
};

void WFN::push_back_SDM(const double &value)
{
	UT_SpinDensityMatrix.push_back(value);
};

void WFN::resize_SDM(const int &size, const double &value)
{
	UT_SpinDensityMatrix.resize(size, value);
};

const double WFN::get_SDM(const int &nr) const
{
	if (nr >= 0 && nr < UT_SpinDensityMatrix.size())
		return UT_SpinDensityMatrix[nr];
	else
	{
		std::cout << "Requested nr out of range! Size: " << UT_SpinDensityMatrix.size() << " nr: " << nr << std::endl;
		return -1;
	}
};

bool WFN::set_SDM(const int &nr, const double &value)
{
	if (nr >= 0 && nr < UT_SpinDensityMatrix.size())
	{
		UT_SpinDensityMatrix[nr] = value;
		return true;
	}
	else
	{
		std::cout << "invalid arguments for set_SDM! Input was: " << nr << ";" << value << std::endl;
		return false;
	}
};

bool WFN::build_DM(std::string basis_set_path, bool debug) {
	using namespace std;
	int elcount = -get_charge();
	if (debug)
		std::cout << "elcount: " << elcount << std::endl;
	for (int i = 0; i < ncen; i++)
	{
		elcount += get_atom_charge(i);
		elcount -= constants::ECP_core_electrons(constants::ECP_electrons_pTB, get_atom_charge(i));
	}
	if (debug)
		std::cout << "elcount after: " << elcount << std::endl;
	int alpha_els = 0, beta_els = 0, temp_els = elcount;
	while (temp_els > 1)
	{
		alpha_els++;
		beta_els++;
		temp_els -= 2;
		if (debug)
			std::cout << temp_els << "\n";
		err_checkf(alpha_els >= 0 && beta_els >= 0, "Error setting alpha and beta electrons! a or b are negative!", std::cout);
		err_checkf(alpha_els + beta_els <= elcount, "Error setting alpha and beta electrons! Sum a + b > elcount!", std::cout);
		err_checkf(temp_els > -elcount, "Error setting alpha and beta electrons! Ran below -elcount!", std::cout);
	}
	alpha_els += temp_els;
	if (debug)
		std::cout << "al/be els:" << alpha_els << " " << beta_els << std::endl;
	const int mult = get_multi();
	int diff = 0;
	if (mult != 0)
		diff = get_multi() - 1;
	if (debug)
		std::cout << "diff: " << diff << std::endl;
	while (alpha_els - beta_els != diff)
	{
		alpha_els++;
		beta_els--;
		err_checkf(alpha_els >= 0 && beta_els >= 0, "Error setting alpha and beta electrons!", std::cout);
	}
	if (debug)
	{
		std::cout << "alpha, beta, elcount: " << setw(5) << alpha_els << setw(5) << beta_els << setw(5) << elcount << endl;
	}
	if (get_nr_basis_set_loaded() == 0)
	{
		if (debug)
			std::cout << "No basis set loaded, will load a complete basis set now!" << endl;
		err_checkf(BasisSetLibrary::read_basis_set_vanilla(basis_set_path, *this, debug), "ERROR during reading of missing basis set!", std::cout);
	}
	else if (get_nr_basis_set_loaded() < get_ncen())
	{
		std::cout << "Not all atoms have a basis set loaded!\nLaoding the missing atoms..." << flush;
		err_checkf(BasisSetLibrary::read_basis_set_missing(basis_set_path, *this, debug), "ERROR during reading of missing basis set!", std::cout);
	}
	else if (get_nr_basis_set_loaded() > get_ncen())
	{
		err_checkf(false, "# of loaded > # atoms\nSorry, this should not happen... aborting!!!", std::cout);
	}
	// set_modified();
	vec CMO;
	vec CMO_beta;
	if (debug)
	{
		std::cout << "Origin: " << get_origin() << endl;
	}
	if (get_origin() == 2 || get_origin() == 4 || get_origin() == 9 || get_origin() == 8)
	{
		//-----------------------check ordering and order accordingly----------------------
		sort_wfn(check_order(debug), debug);
		//---------------normalize basis set---------------------------------
		if (debug)
			std::cout << "starting to normalize the basis set" << endl;
		vec norm_const;
		//-----------debug output---------------------------------------------------------
		if (debug)
		{
			std::cout << "exemplary output before norm_const of the first atom with all it's properties: " << endl;
			print_atom_long(0);
			std::cout << "ended normalizing the basis set, now for the MO_coeffs" << endl;
			std::cout << "Status report:" << endl;
			std::cout << "size of norm_const: " << norm_const.size() << endl;
			std::cout << "WFN MO counter: " << get_nmo() << endl;
			std::cout << "Number of atoms: " << get_ncen() << endl;
			std::cout << "Primitive count of zero MO: " << get_MO_primitive_count(0) << endl;
			std::cout << "Primitive count of first MO: " << get_MO_primitive_count(1) << endl;
		}

		//-------------------normalize the basis set shell wise into a copy vector---------
		//The shell loop below spells out the cartesian component norms of s, p, d and f by hand and
		//has nothing for g: a g shell used to leave `factor` at the previous shell's value, push no
		//constants at all, and every later shell then read norm_const one shell off. Say so instead.
		//Generalising it needs the component order this file's g types are written in, which is a
		//convention no reader here agrees on yet - and build_DM has no caller outside the tests.
		for (int a = 0; a < get_ncen(); a++)
			for (int s = 0; s < get_atom_shell_count(a); s++)
			{
				const int type = get_shell_type(a, s);
				err_checkf(type >= 1, "The type of shell " + std::to_string(s) + " of atom " +
					std::to_string(a) + " was never read", std::cout);
				if (type > 4)
				{
					std::cout << "build_DM normalises s, p, d and f shells only; shell " << s
						<< " of atom " << a << " is of type " << type
						<< " (l = " << type - 1 << "). Refusing rather than building a density "
						"matrix from normalisation constants that are one shell out of step.\n";
					return false;
				}
			}
		vec2 basis_coefficients(get_ncen());
#pragma omp parallel for
		for (int a = 0; a < get_ncen(); a++)
		{
			for (int p = 0; p < get_atom_primitive_count(a); p++)
			{
				//This was the same s/p/d/f switch that lost a g shell in
				//get_shell_start_in_primitives: a type of 5 matched no case, so the primitive kept
				//its raw exponent as a normalisation constant. constants::axial_prim_norm is the
				//general-l form of exactly these four numbers. The types are validated above,
				//serially - err_checkf exits, and exiting from inside an OpenMP region does not.
				const double temp_c =
					constants::axial_prim_norm(get_atom_primitive_type(a, p) - 1,
						get_atom_basis_set_exponent(a, p)) *
					get_atom_basis_set_coefficient(a, p);
				if (debug)
					std::cout << "temp_c:" << temp_c << "\n";
				basis_coefficients[a].push_back(temp_c);
			}
		}
		for (int a = 0; a < get_ncen(); a++)
		{
			double aiaj = 0.0;
			double factor = 0.0;
			for (int s = 0; s < get_atom_shell_count(a); s++)
			{
				int type_temp = get_shell_type(a, s);
				err_checkf(type_temp != -1, "ERROR in type assignement!!", std::cout);
				if (debug)
				{
					std::cout << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << "\n"
						<< "start: " << get_shell_start(a, s)
						<< " stop: " << get_shell_end(a, s) << "\n"
						<< "factor: ";
				}
				switch (type_temp)
				{
				case 1:
					factor = 0;
					for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
					{
						for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
						{
							aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
							double term = constants::PI3 / pow(aiaj, 3);
							term = pow(term, 0.5);
							factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
						}
					}
					if (factor == 0)
						return false;
					factor = pow(factor, -0.5);
					if (debug)
						std::cout << factor << "\n";
					for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
					{
						if (debug)
						{
							std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i)
								<< " Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << "\n";
						}
						// contraction_coefficients[a][i] = factor * get_atom_basis_set_coefficient(a, i);
						basis_coefficients[a][i] *= factor;
						norm_const.push_back(basis_coefficients[a][i]);
					}
					break;
				case 2:
					factor = 0;
					for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
					{
						for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
						{
							aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
							double term = constants::PI3 / (4 * pow(aiaj, 5));
							term = pow(term, 0.5);
							factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
						}
					}
					if (factor == 0)
						return false;
					factor = pow(factor, -0.5);
					if (debug)
						std::cout << factor << "\n";
					for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
					{
						if (debug)
						{
							std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i)
								<< " Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << "\n";
						}
						// contraction_coefficients[a][i] = factor * get_atom_basis_set_coefficient(a, i);
						basis_coefficients[a][i] *= factor;
						for (int k = 0; k < 3; k++)
							norm_const.push_back(basis_coefficients[a][i]);
					}
					break;
				case 3:
					factor = 0;
					for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
					{
						for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
						{
							aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
							double term = constants::PI3 / (16 * pow(aiaj, 7));
							term = pow(term, 0.5);
							factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
						}
					}
					if (factor == 0)
						return false;
					factor = (pow(factor, -0.5)) / sqrt(3);
					if (debug)
						std::cout << factor << "\n";
					for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
					{
						if (debug)
						{
							std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i)
								<< " Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << "\n";
						}
						// contraction_coefficients[a][i] = factor * get_atom_basis_set_coefficient(a, i);
						basis_coefficients[a][i] *= factor;
						for (int k = 0; k < 3; k++)
							norm_const.push_back(basis_coefficients[a][i]);
						for (int k = 0; k < 3; k++)
							norm_const.push_back(sqrt(3) * basis_coefficients[a][i]);
					}
					break;
				case 4:
					factor = 0;
					for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
					{
						for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
						{
							aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
							double term = constants::PI3 / (64 * pow((aiaj), 9));
							term = pow(term, 0.5);
							factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
						}
					}
					if (factor == 0)
						return false;
					factor = pow(factor, -0.5) / sqrt(15);
					if (debug)
						std::cout << factor << "\n";
					for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
					{
						if (debug)
						{
							std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i)
								<< " Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << "\n";
						}
						// contraction_coefficients[a][i] = factor * get_atom_basis_set_coefficient(a, i);
						basis_coefficients[a][i] *= factor;
						for (int l = 0; l < 3; l++)
							norm_const.push_back(basis_coefficients[a][i]);
						for (int l = 0; l < 6; l++)
							norm_const.push_back(sqrt(5) * basis_coefficients[a][i]);
						norm_const.push_back(sqrt(15) * basis_coefficients[a][i]);
					}
					break;
				}
				if (debug)
					std::cout << "This shell has: " << get_shell_end(a, s) - get_shell_start(a, s) + 1 << " primitives\n";
			}
		}
		//-----------debug output---------------------------------------------------------
		if (debug)
		{
			std::cout << "exemplary output of the first atom with all it's properties: " << endl;
			print_atom_long(0);
			std::cout << "ended normalizing the basis set, now for the norm_cprims" << endl;
			std::cout << "Status report:" << endl;
			std::cout << "size of norm_const: " << norm_const.size() << endl;
			std::cout << "WFN MO counter: " << get_nmo() << endl;
			std::cout << "Number of atoms: " << get_ncen() << endl;
			std::cout << "Primitive count of zero MO: " << get_MO_primitive_count(0) << endl;
			std::cout << "Primitive count of first MO: " << get_MO_primitive_count(1) << endl;
		}
		//---------------------To not mix up anything start normalizing WFN_matrix now--------------------------
		int run = 0;
		vec2 changed_coefs;
		changed_coefs.resize(get_nmo());
		if (debug)
		{
			std::cout << "Opening norm_cprim!" << endl;
			ofstream norm_cprim("norm_prim.debug", ofstream::out);
			for (int m = 0; m < get_nmo(); m++)
			{
				norm_cprim << m << ". MO:\n";
				changed_coefs[m].resize(get_MO_primitive_count(m), 0.0);
				for (int p = 0; p < get_MO_primitive_count(m); p++)
				{
					changed_coefs[m][p] = get_MO_coef(m, p) / norm_const[p];
					if (m == 0)
						std::cout << p << ". primitive; " << m << ". MO "
						<< "norm nonst: " << norm_const[p]
						<< " temp after normalization: " << changed_coefs[m][p] << "\n";
					norm_cprim << " " << changed_coefs[m][p] << "\n";
					run++;
				}
			}
			norm_cprim.flush();
			norm_cprim.close();
			std::cout << "See norm_cprim.debug for the CPRIM vectors" << endl;
			std::cout << "Total count in CPRIM: " << run << endl;
		}
		else
		{
#pragma omp parallel for
			for (int m = 0; m < get_nmo(); m++)
			{
				changed_coefs[m].resize(get_MO_primitive_count(m), 0.0);
				for (int p = 0; p < get_MO_primitive_count(m); p++)
				{
					changed_coefs[m][p] = get_MO_coef(m, p) / norm_const[p];
				}
			}
		}
		//--------------Build CMO of alessandro from the first elements of each shell-------------
		int nao = 0;
		for (int a = 0; a < get_ncen(); a++)
		{
			for (int s = 0; s < get_atom_shell_count(a); s++)
			{
				switch (get_shell_type(a, s))
				{
				case 1:
					nao++;
					break;
				case 2:
					nao += 3;
					break;
				case 3:
					nao += 6;
					break;
				case 4:
					nao += 10;
					break;
				}
			}
		}
		int nshell = 0;
		for (int m = 0; m < get_nmo(); m++)
		{
			int run_2 = 0;
			for (int a = 0; a < get_ncen(); a++)
			{
				for (int s = 0; s < get_atom_shell_count(a); s++)
				{
					// if (debug)std::cout << "Going to load the " << get_shell_start_in_primitives(a, s) << ". value\n"l;
					switch (get_shell_type(a, s))
					{
					case 1:
						CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s)]);
						if (debug && get_atom_shell_primitives(a, s) != 1 && m == 0)
							std::cout << "Pushing back 1 coefficient for S shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives! Shell start is: " << get_shell_start(a, s) << "\n";
						break;
					case 2:
						for (int i = 0; i < 3; i++)
							CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
						if (debug && get_atom_shell_primitives(a, s) != 1 && m == 0)
							std::cout << "Pushing back 3 coefficients for P shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!\n";
						break;
					case 3:
						for (int i = 0; i < 6; i++)
							CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
						if (debug && get_atom_shell_primitives(a, s) != 1 && m == 0)
							std::cout << "Pushing back 6 coefficient for D shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!\n";
						break;
					case 4:
						// this hardcoded piece is due to the order of f-type functions in the fchk
						for (int i = 0; i < 3; i++)
							CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
						CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 6]);
						for (int i = 0; i < 2; i++)
							CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i + 3]);
						for (int i = 0; i < 2; i++)
							CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i + 7]);
						CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 5]);
						CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 9]);
						if (debug && get_atom_shell_primitives(a, s) != 1 && m == 0)
							std::cout << "Pushing back 10 coefficient for F shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!\n";
						break;
					}
					run_2++;
				}
				if (debug && m == 0)
					std::cout << "finished with atom!\n";
			}
			if (debug)
				std::cout << "finished with MO!\n";
			if (nshell != run_2)
				nshell = run_2;
		}
		if (is_unrestricted)
		{
			for (int m = alpha_els; m < alpha_els + beta_els; m++)
			{
				int run_2 = 0;
				for (int a = 0; a < get_ncen(); a++)
				{
					for (int s = 0; s < get_atom_shell_count(a); s++)
					{
						if (debug)
							std::cout << "Going to load the " << get_shell_start_in_primitives(a, s) << ". value\n";
						switch (get_shell_type(a, s))
						{
						case 1:
							CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s)]);
							if (m == 0)
								nao++;
							if (debug && get_atom_shell_primitives(a, s) != 1)
								std::cout << "Pushing back 1 coefficient for S shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives! Shell start is: " << get_shell_start(a, s) << "\n";
							break;
						case 2:
							for (int i = 0; i < 3; i++)
								CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
							if (debug && get_atom_shell_primitives(a, s) != 1)
								std::cout << "Pushing back 3 coefficients for P shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!\n";
							if (m == 0)
								nao += 3;
							break;
						case 3:
							for (int i = 0; i < 6; i++)
								CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
							if (debug && get_atom_shell_primitives(a, s) != 1)
								std::cout << "Pushing back 6 coefficient for D shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!\n";
							if (m == 0)
								nao += 6;
							break;
						case 4:
							// this hardcoded piece is due to the order of f-type functions in the fchk
							for (int i = 0; i < 3; i++)
								CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
							CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 6]);
							for (int i = 0; i < 2; i++)
								CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i + 3]);
							for (int i = 0; i < 2; i++)
								CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i + 7]);
							CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 5]);
							CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 9]);
							if (debug && get_atom_shell_primitives(a, s) != 1)
								std::cout << "Pushing back 10 coefficient for F shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!\n";
							if (m == 0)
								nao += 10;
							break;
						}
						run_2++;
					}
					if (debug)
						std::cout << "finished with atom!\n";
				}
				if (debug)
					std::cout << "finished with MO!\n";
				if (nshell != run_2)
					nshell = run_2;
			}
		}

		if (debug)
		{
			ofstream cmo("cmo.debug", ofstream::out);
			for (int p = 0; p < CMO.size(); p++)
			{
				for (int i = 0; i < 5; i++)
				{
					cmo << scientific << setw(14) << setprecision(7) << CMO[p + i] << " ";
				}
				p += 4;
				cmo << "\n";
			}
			cmo.flush();
			cmo.close();
			std::cout << CMO.size() << " Elements in CMO" << endl;
			std::cout << norm_const.size() << " = nprim" << endl;
			std::cout << nao << " = nao" << endl;
			std::cout << nshell << " = nshell" << endl;
		}
		//------------------ make the DM -----------------------------
		int naotr = nao * (nao + 1) / 2;
		vec kp;
		resize_DM(naotr, 0.0);
		if (is_unrestricted)
			resize_SDM(naotr, 0.0);
		if (debug)
		{
			std::cout << "I made kp!" << endl
				<< nao << " is the maximum for iu" << endl;
			std::cout << "Making DM now!" << endl;
		}
		for (int iu = 0; iu < nao; iu++)
		{
#pragma omp parallel for
			for (int iv = 0; iv <= iu; iv++)
			{
				const int iuv = (iu * (iu + 1) / 2) + iv;
				// if (debug)std::cout << "iu: " << iu << " iv: " << iv << " iuv: " << iuv << " kp(iu): " << iu * (iu + 1) / 2 << endl;
				double temp;
				// if (debug)std::cout << "Working on MO: ";
				for (int m = 0; m < get_nmo(); m++)
				{
					// if (debug && m == 0)std::cout << m << " " << flush;
					// else if (debug && m != get_nmo() - 1)std::cout << "." << flush;
					// elsestd::cout << get_nmo() - 1 << flush;
					if (is_unrestricted)
					{
						if (m < alpha_els)
						{
							temp = get_MO_occ(m) * CMO[iu + (m * nao)] * CMO[iv + (m * nao)];
							err_checkf(set_SDM(iuv, get_SDM(iuv) + temp), "Something went wrong while writing the SDM! iuv=" + to_string(iuv), std::cout);
							err_checkf(set_DM(iuv, get_DM(iuv) + temp), "Something went wrong while writing the DM! iuv=" + to_string(iuv), std::cout);
						}
						else
						{
							temp = get_MO_occ(m) * CMO_beta[iu + ((m - alpha_els) * nao)] * CMO_beta[iv + ((m - alpha_els) * nao)];
							err_checkf(set_SDM(iuv, get_SDM(iuv) - temp), "Something went wrong while writing the SDM! iuv=" + to_string(iuv), std::cout);
							err_checkf(set_DM(iuv, get_DM(iuv) + temp), "Something went wrong while writing the DM! iuv=" + to_string(iuv), std::cout);
						}
					}
					else
					{
						if (get_MO_occ(m) == 0.0)
							continue;
						temp = get_MO_occ(m) * CMO[iu + (m * nao)] * CMO[iv + (m * nao)];
						err_checkf(set_DM(iuv, get_DM(iuv) + temp), "Something went wrong while writing the DM!", std::cout);
					}
					// else if (debug)std::cout << "DM after: " << get_DM(iuv) << endl;
				}
				// if (debug)std::cout << endl;
			}
		}
	}
	else
	{
		std::cout << "Sorry, this origin is not supported yet!" << endl;
		return false;
	}
	return true;
};

int WFN::check_order(const bool &debug) const
{
	for (int i = 0; i < ncen; i++)
	{
		if (!get_atom_basis_set_loaded(i))
		{
			std::cout << "Sorry, consistency check only works if basis set is loaded for all atoms!\n"
				<< "Failing atom: " << i << " " << get_atom_label(i) << "\n";
			return -1;
		}
	}
	//---------------------------check if the wfn is in the right order----------------------
	int order = 0;   // 1=gaussian (P=2222 3333 4444) 2=tonto (234 234 234 234) 3=ORCA (423 423 423 423)
	int f_order = 0; // 1=gaussian (F=11 12 13 17 14 15 18 19 16 20) 2=tonto=ORCA 3=ORCA (11 12 13 14 15 17 16 18 19 20) 4=natural (11 12 13 14 15 16 17 18 19 20)
	int primcounter = 0;
	bool order_found = false;
	for (int a = 0; a < get_ncen(); a++)
	{
		for (int s = 0; s < get_atom_shell_count(a); s++)
		{
			int type = get_shell_type(a, s);
			switch (type)
			{
			case 1: // S Orbital
			{
				for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
				{
					if (types[primcounter] != 1)
					{
						order = -1;
						if (debug)
						{
							std::cout << "This should not happen, the order of your file is not ok for S-types! Checked #:" << primcounter << "\n";
						}
					}
					else
						primcounter++;
				}
				break;
			}
			case 2: // P Orbital
			{
				if (order_found)
				{
					if (order == 1)
					{
						for (int r = 0; r < get_atom_shell_primitives(a, s); r++)
						{
							if (types[primcounter] != 2 || types[primcounter + get_atom_shell_primitives(a, s)] != 3 || types[primcounter + 2 * get_atom_shell_primitives(a, s)] != 4)
							{
								if (debug)
								{
									std::cout << "The found order does not match all type entries! primcounter: " << primcounter << "\n";
								}
							}
							primcounter++;
						}
						primcounter += 2 * get_atom_shell_primitives(a, s);
					}
					else if (order == 2)
					{
						for (int r = 0; r < get_atom_shell_primitives(a, s); r++)
						{
							if (types[primcounter] != 2 || types[primcounter + 1] != 3 || types[primcounter + 2] != 4)
							{
								if (debug)
								{
									std::cout << "The found order does not match all type entries! primcounter: " << primcounter << "\n";
								}
								return -1;
							}
							primcounter += 3;
						}
					}
					else if (order == 3)
					{
						for (int r = 0; r < get_atom_shell_primitives(a, s); r++)
						{
							if (types[primcounter] != 4 || types[primcounter + 1] != 2 || types[primcounter + 2] != 3)
							{
								if (debug)
								{
									std::cout << "The found order does not match all type entries! primcounter: " << primcounter << "\n";
								}
								return -1;
							}
							primcounter += 3;
						}
					}
				}
				else
				{
					if (types[primcounter] == 2)
					{
						if (debug && a == 0)
						{
							std::cout << "Seems to be either tonto or gaussian file...\n";
						}
						if (types[primcounter + 1] == 3)
						{
							order = 2;
							order_found = true;
							if (debug)
							{
								std::cout << "This wfn file is in tonto order!\n";
							}
						}
						else if (types[primcounter + 1] == 2 && get_atom_shell_primitives(a, s) > 1)
						{
							order = 1;
							order_found = true;
							if (debug)
							{
								std::cout << "This wfn file is in gaussian order!\n";
							}
						}
						else
						{
							order = 1;
							order_found = true;
							if (debug)
							{
								std::cout << "Either some error or this shell only has 1 p-primitive and "
									<< "i didn't find any order yet... assuming gaussian\n";
							}
						}
					}
					else if (types[primcounter] == 4)
					{
						if (debug && a == 0)
						{
							std::cout << "Seems as if this file was ORCA ordered...\n";
						}
						order = 3;
						if (types[primcounter + 1] == 2)
						{
							if (debug)
							{
								std::cout << "Yep, that's right! making it permanent now!\n";
							}
							order_found = true;
						}
					}
					else
					{
						std::cout << "I can't recognize this order of the .wfn file...\n";
						return -1;
					}
					s--;
				}
				break;
			}
			case 3: // D Orbital
			{
				if (order_found)
				{
					switch (order)
					{
					case 1:
					{
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
						{
							if (types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] != 5 || types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] != 6 || types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] != 7 || types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] != 8 || types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] != 9 || types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] != 10)
							{
								order = -1;
								if (debug)
								{
									std::cout << "The checked types are 6 from #" << primcounter << " and are:\n"
										<< types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] << "\n";
								}
								std::cout << "Something seems to be wrong in the order of your D-Types...\n";
							}
							else
								primcounter += 6;
						}
					}
					break;
					case 2:
					case 3:
					{
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
						{
							if (types[get_shell_start_in_primitives(a, s) + 0 + 6 * i] != 5 || types[get_shell_start_in_primitives(a, s) + 1 + 6 * i] != 6 || types[get_shell_start_in_primitives(a, s) + 2 + 6 * i] != 7 || types[get_shell_start_in_primitives(a, s) + 3 + 6 * i] != 8 || types[get_shell_start_in_primitives(a, s) + 4 + 6 * i] != 9 || types[get_shell_start_in_primitives(a, s) + 5 + 6 * i] != 10)
							{
								order = -1;
								if (debug)
								{
									std::cout << "The checked types are 6 from #" << primcounter << " and are:\n"
										<< types[get_shell_start_in_primitives(a, s) + 0 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 1 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 2 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 3 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 4 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 5 + 6 * i] << "\n";
								}
								std::cout << "Something seems to be wrong in the order of your D-Types...\n";
							}
							else
								primcounter += 6;
						}
					}
					break;
					}
				}
				else
				{
					std::cout << "That's highly suspicious, no order for P-type but D-type functions? Let me stop before i do something stupid!\n";
					return -1;
				}
			}
			break;
			case 4:
				for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
				{
					if (types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] != 11 || types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] != 12 || types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] != 13 || types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] != 17 || types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] != 14 || types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] != 15 || types[get_shell_start_in_primitives(a, s) + 6 * get_atom_shell_primitives(a, s) + i] != 18 || types[get_shell_start_in_primitives(a, s) + 7 * get_atom_shell_primitives(a, s) + i] != 19 || types[get_shell_start_in_primitives(a, s) + 8 * get_atom_shell_primitives(a, s) + i] != 16 || types[get_shell_start_in_primitives(a, s) + 9 * get_atom_shell_primitives(a, s) + i] != 20)
					{
						if (types[primcounter] == 11 || types[primcounter + 1] == 12 || types[primcounter + 2] == 13 || types[primcounter + 3] == 14 || types[primcounter + 4] == 15 || types[primcounter + 5] == 16 || types[primcounter + 6] == 17 || types[primcounter + 7] == 18 || types[primcounter + 8] == 19 || types[primcounter + 9] == 20)
						{
							f_order = 4;
							if (debug)
							{
								std::cout << "The checked types are 10 from #" << primcounter << " and are:\n"
									<< types[primcounter] << " " << types[primcounter + 1] << " " << types[primcounter + 2] << " "
									<< types[primcounter + 3] << " " << types[primcounter + 4] << " " << types[primcounter + 5] << " "
									<< types[primcounter + 6] << " " << types[primcounter + 7] << " " << types[primcounter + 8] << " "
									<< types[primcounter + 9] << "\n"
									<< "Appears to be already okay...\n";
							}
							primcounter += 10;
						}
						else if (types[primcounter] != 11 || types[primcounter + 1] != 12 || types[primcounter + 2] != 13 || types[primcounter + 3] != 14 || types[primcounter + 4] != 15 || types[primcounter + 5] != 17 || types[primcounter + 6] != 16 || types[primcounter + 7] != 18 || types[primcounter + 8] != 19 || types[primcounter + 9] != 20)
						{
							order = -1;
							if (debug)
							{
								std::cout << "The checked types are 10 from #" << primcounter << " and are:\n"
									<< types[primcounter] << " " << types[primcounter + 1] << " " << types[primcounter + 2] << " "
									<< types[primcounter + 3] << " " << types[primcounter + 4] << " " << types[primcounter + 5] << " "
									<< types[primcounter + 6] << " " << types[primcounter + 7] << " " << types[primcounter + 8] << " "
									<< types[primcounter + 9] << "\n"
									<< "Something seems to be wrong in the order of your F-Types...\n";
							}
						}
						else
						{
							f_order = 3;
							primcounter += 10;
						}
					}
					else
					{
						f_order = 1;
						primcounter += 10;
					}
				}
				break;
			default:
				std::cout << "ERROR in type assignement!!\n";
				return -1;
				break;
			}
		}
	}
	/*if(debug){
	   std::cout << "Going to return " << f_order*10+order << endl;
		Enter();
	}*/
	return (f_order * 10 + order);
};

bool WFN::sort_wfn(const int &g_order, const bool &debug)
{
	set_modified();
	int primcounter = 0;
	int f_order = 0;
	int order = g_order;
	// Sorry for this way of forwarding the order, i think 2 switches would have been more nicely
	while (order >= 10)
	{
		f_order++;
		order -= 10;
	}
	switch (order)
	{
	case 1:
		for (int a = 0; a < get_ncen(); a++)
			for (int s = 0; s < get_atom_shell_count(a); s++)
				switch (get_shell_type(a, s))
				{
				case 1:
					primcounter += get_atom_shell_primitives(a, s);
					break;
				case 2:
				{
					if (get_atom_shell_primitives(a, s) > 1)
					{
						ivec temp_centers;
						temp_centers.resize(3 * get_atom_shell_primitives(a, s));
						ivec temp_types;
						temp_types.resize(3 * get_atom_shell_primitives(a, s));
						vec temp_exponents;
						temp_exponents.resize(3 * get_atom_shell_primitives(a, s));
						vec2 temp_MO_coefficients;
						temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 3);
						for (int i = 0; i < get_atom_shell_primitives(a, s) * 3; i++)
							temp_MO_coefficients[i].resize(nmo);
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
							for (int c = 0; c < 3; c++)
							{
								temp_centers[3 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_types[3 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_exponents[3 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
								for (int m = 0; m < nmo; m++)
									temp_MO_coefficients[3 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
							}
						for (int i = 0; i < 3 * get_atom_shell_primitives(a, s); i++)
						{
							centers[primcounter + i] = temp_centers[i];
							types[primcounter + i] = temp_types[i];
							exponents[primcounter + i] = temp_exponents[i];
							for (int m = 0; m < nmo; m++)
								err_checkf(set_MO_coef(m, primcounter + i, temp_MO_coefficients[i][m]), "Error while assigning new MO coefficient!", std::cout);
						}
					}
					primcounter += get_atom_shell_primitives(a, s) * 3;
					break;
				}
				case 3:
				{
					if (get_atom_shell_primitives(a, s) > 1)
					{
						ivec temp_centers;
						temp_centers.resize(6 * get_atom_shell_primitives(a, s));
						ivec temp_types;
						temp_types.resize(6 * get_atom_shell_primitives(a, s));
						vec temp_exponents;
						temp_exponents.resize(get_atom_shell_primitives(a, s) * 6);
						vec2 temp_MO_coefficients;
						temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 6);
						for (int i = 0; i < get_atom_shell_primitives(a, s) * 6; i++)
							temp_MO_coefficients[i].resize(nmo);
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
							for (int c = 0; c < 6; c++)
							{
								temp_centers[6 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_types[6 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_exponents[6 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
								for (int m = 0; m < nmo; m++)
									temp_MO_coefficients[6 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
							}
						for (int i = 0; i < 6 * get_atom_shell_primitives(a, s); i++)
						{
							centers[primcounter + i] = temp_centers[i];
							types[primcounter + i] = temp_types[i];
							exponents[primcounter + i] = temp_exponents[i];
							for (int m = 0; m < nmo; m++)
								if (!set_MO_coef(m, primcounter + i, temp_MO_coefficients[i][m]))
								{
									std::cout << "Error while assigning new MO coefficient!\n";
									return false;
								}
						}
					}
					primcounter += get_atom_shell_primitives(a, s) * 6;
					break;
				}
				case 4:
				{
					if (get_atom_shell_primitives(a, s) > 1)
					{
						ivec temp_centers;
						temp_centers.resize(10 * get_atom_shell_primitives(a, s));
						ivec temp_types;
						temp_types.resize(10 * get_atom_shell_primitives(a, s));
						vec temp_exponents;
						temp_exponents.resize(get_atom_shell_primitives(a, s) * 10);
						vec2 temp_MO_coefficients;
						temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 10);
						for (int i = 0; i < get_atom_shell_primitives(a, s) * 10; i++)
							temp_MO_coefficients[i].resize(nmo);
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
							for (int c = 0; c < 10; c++)
							{
								temp_centers[10 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_types[10 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_exponents[10 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
								for (int m = 0; m < nmo; m++)
									temp_MO_coefficients[10 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
							}
						for (int i = 0; i < 10 * get_atom_shell_primitives(a, s); i++)
						{
							centers[primcounter + i] = temp_centers[i];
							types[primcounter + i] = temp_types[i];
							exponents[primcounter + i] = temp_exponents[i];
							for (int m = 0; m < nmo; m++)
								if (!set_MO_coef(m, primcounter + i, temp_MO_coefficients[i][m]))
								{
									std::cout << "Error while assigning new MO coefficient!\n";
									return false;
								}
						}
					}
					primcounter += get_atom_shell_primitives(a, s) * 10;
					break;
				}
				}
		break;
	case 2:
		if (debug)
		{
			std::cout << "Nothing to be done here, i like tonto type..." << std::endl;
		}
		break;
	case 3:
		for (int a = 0; a < get_ncen(); a++)
		{
			for (int s = 0; s < get_atom_shell_count(a); s++)
			{
				switch (get_shell_type(a, s))
				{
				case 1:
					primcounter += get_atom_shell_primitives(a, s);
					break;
				case 2:
				{
					int temp_center;
					int temp_type;
					double temp_exponent;
					vec temp_MO_coefficients;
					temp_MO_coefficients.resize(nmo);
					for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
					{
						temp_center = centers[primcounter];
						temp_type = types[primcounter];
						temp_exponent = exponents[primcounter];
						for (int m = 0; m < nmo; m++)
							temp_MO_coefficients[m] = MOs[m].get_coefficient(primcounter);
						for (int j = 0; j < 2; j++)
						{
							centers[primcounter + j] = centers[primcounter + 1 + j];
							types[primcounter + j] = types[primcounter + 1 + j];
							exponents[primcounter + j] = exponents[primcounter + 1 + j];
							for (int m = 0; m < nmo; m++)
							{
								err_checkf(set_MO_coef(m, primcounter + j, MOs[m].get_coefficient(primcounter + 1 + j)), "Error while assigning new MO coefficient!", std::cout);
							}
						}
						centers[primcounter + 2] = temp_center;
						types[primcounter + 2] = temp_type;
						exponents[primcounter + 2] = temp_exponent;
						for (int m = 0; m < nmo; m++)
							if (!set_MO_coef(m, primcounter + 2, temp_MO_coefficients[m]))
							{
								std::cout << "Error while assigning new MO coefficient!\n";
								return false;
							}
						primcounter += 3;
					}
					break;
				}
				case 3:
					primcounter += get_atom_shell_primitives(a, s) * 6;
					break;
				case 4:
					primcounter += get_atom_shell_primitives(a, s) * 10;
					break;
				}
			}
		}
		break;
	default:
		std::cout << "order type: " << f_order << " " << order << " not supported!" << std::endl;
		return false;
		break;
	}
	primcounter = 0;
	switch (f_order)
	{
	case 1:
		for (int a = 0; a < get_ncen(); a++)
			for (int s = 0; s < get_atom_shell_count(a); s++)
				switch (get_shell_type(a, s))
				{
				case 1:
					primcounter += get_atom_shell_primitives(a, s);
					break;
				case 2:
					primcounter += get_atom_shell_primitives(a, s) * 3;
					break;
				case 3:
					primcounter += get_atom_shell_primitives(a, s) * 6;
					break;
				case 4:
				{
					ivec temp_center;
					temp_center.resize(10);
					ivec temp_type;
					temp_type.resize(10);
					vec temp_exponent;
					temp_exponent.resize(10);
					vec2 temp_MO_coefficients;
					temp_MO_coefficients.resize(nmo);
					for (int m = 0; m < nmo; m++)
						temp_MO_coefficients[m].resize(10);
					for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
					{
						for (int j = 0; j < 10; j++)
						{
							temp_center[j] = centers[get_shell_start_in_primitives(a, s) + 10 * i + j];
							temp_type[j] = types[get_shell_start_in_primitives(a, s) + 10 * i + j];
							temp_exponent[j] = exponents[get_shell_start_in_primitives(a, s) + 10 * i + j];
							for (int m = 0; m < nmo; m++)
								temp_MO_coefficients[m][j] = MOs[m].get_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + j);
						}
						// mask[j] is where gaussian's j-th f function (11 12 13 17 14 15 18 19 16 20) lands in natural order
						ivec mask{ 0, 1, 2, 6, 3, 4, 7, 8, 5, 9 };
						for (int j = 0; j < 10; j++)
						{
							centers[get_shell_start_in_primitives(a, s) + 10 * i + mask[j]] = temp_center[j];
							types[get_shell_start_in_primitives(a, s) + 10 * i + mask[j]] = temp_type[j];
							exponents[get_shell_start_in_primitives(a, s) + 10 * i + mask[j]] = temp_exponent[j];
							for (int m = 0; m < nmo; m++)
								set_MO_coef(m, get_shell_start_in_primitives(a, s) + 10 * i + mask[j], temp_MO_coefficients[m][j]);
						}
					}
					primcounter += 10;
				}
				break;
				}
		break;
	case 2:
	case 3:
		for (int a = 0; a < get_ncen(); a++)
			for (int s = 0; s < get_atom_shell_count(a); s++)
				switch (get_shell_type(a, s))
				{
				case 1:
					primcounter += get_atom_shell_primitives(a, s);
					break;
				case 2:
					primcounter += get_atom_shell_primitives(a, s) * 3;
					break;
				case 3:
					primcounter += get_atom_shell_primitives(a, s) * 6;
					break;
				case 4:
				{
					ivec temp_center;
					temp_center.resize(10);
					ivec temp_type;
					temp_type.resize(10);
					vec temp_exponent;
					temp_exponent.resize(10);
					vec2 temp_MO_coefficients;
					temp_MO_coefficients.resize(nmo);
					for (int m = 0; m < nmo; m++)
						temp_MO_coefficients[m].resize(10);
					for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
					{
						for (int j = 0; j < 10; j++)
						{
							temp_center[j] = centers[get_shell_start_in_primitives(a, s) + 10 * i + j];
							temp_type[j] = types[get_shell_start_in_primitives(a, s) + 10 * i + j];
							temp_exponent[j] = exponents[get_shell_start_in_primitives(a, s) + 10 * i + j];
							for (int m = 0; m < nmo; m++)
								temp_MO_coefficients[m][j] = MOs[m].get_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + j);
						}
						ivec mask{ 0, 1, 2, 3, 4, 6, 5, 7, 8, 9 };
						for (int j = 0; j < 10; j++)
						{
							centers[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_center[mask[j]];
							types[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_type[mask[j]];
							exponents[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_exponent[mask[j]];
							for (int m = 0; m < nmo; m++)
								set_MO_coef(m, get_shell_start_in_primitives(a, s) + 10 * i + mask[j], temp_MO_coefficients[m][j]);
						}
					}
					primcounter += 10;
				}
				break;
				}
		break;
	case 4:
		if (debug)
		{
			std::cout << "This is fine, i like them well ordered!" << std::endl;
		}
		break;
	case 0:
		if (debug)
		{
			std::cout << "There seems to be no f-functions!" << std::endl;
		}
		break;
	default:
		std::cout << "f-order type: " << f_order << " not supported!" << std::endl;
		return false;
		break;
	}
	return true;
};

void WFN::set_has_ECPs(const bool &in, const bool &apply_to_atoms, const int &ECP_mode)
{
	has_ECPs = in;
	ECP_m = ECP_mode;
	if (!apply_to_atoms)
		return;
	//1 = def2, 2 = xTB, 3 = pTB. One loop over the atoms for all three, because the only difference
	//between them was which table a table lookup read - and the bound-checked accessor is what keeps
	//an atom heavier than the tables describe from reading past the end of one. A table lookup per
	//atom does not need a thread each, so the three omp loops this replaces are no loss.
	auto core_of = [ECP_mode](const int Z)
		{
			switch (ECP_mode)
			{
			case 1: return constants::ECP_core_electrons(constants::ECP_electrons, Z);
			case 2: return constants::ECP_core_electrons(constants::ECP_electrons_xTB, Z);
			case 3: return constants::ECP_core_electrons(constants::ECP_electrons_pTB, Z);
			default: return 0;
			}
		};
	if (ECP_mode < 1 || ECP_mode > 3)
		return;
	//An atom past the end of the tables gets zero from the accessor, which is the right answer - no
	//core is defined for it - but a silent zero on a run that asked for ECP cores would quietly
	//count its cores as valence. Name the heaviest one instead; the number is in the tables, not here.
	int heaviest_beyond_tables = 0;
	for (int i = 0; i < ncen; i++)
	{
		const int Z = get_atom_charge(i);
		if (Z > constants::heaviest_ECP_element)
			heaviest_beyond_tables = std::max(heaviest_beyond_tables, Z);
		//Zero means the table declares no core for this element, and that is not the same as a
		//declaration of zero: a count read from the file itself stays, as it did before.
		const int core = core_of(Z);
		if (core != 0)
			atoms[i].set_ECP_electrons(core);
	}
	if (heaviest_beyond_tables != 0)
		std::cout << "\nECP cores were asked for, but the tables stop at Z = "
		<< constants::heaviest_ECP_element << " and this structure contains Z = "
		<< heaviest_beyond_tables << ": those atoms are treated as all-electron, which is what "
		"their basis set has to be for the electron count to add up.\n";
};

void WFN::set_ECPs(ivec &nr, ivec &elcount)
{
	has_ECPs = true;
	err_checkf(nr.size() == elcount.size(), "mismatch in size of atoms and ECP electrons!", std::cout);
#pragma omp parallel for
	for (int i = 0; i < ncen; i++)
	{
		for (int j = 0; j < nr.size(); j++)
		{
			std::cout << "checking " << get_atom_charge(i) << " against " << nr[j] << "\n";
			if (get_atom_charge(i) == nr[j])
			{
				std::cout << "Adding " << elcount[j] << " electron to atom " << i << "with charge " << get_atom_charge(i) << "\n";
				atoms[i].set_ECP_electrons(elcount[j]);
				break;
			}
		}
	}
};

WFN::WFN(const WFN &right)
{
	reset();
	*this = right;
};

WFN &WFN::operator=(const WFN &right)
{
	if (this == &right)
		return *this;
	reset();
	ncen = right.ncen;
	nfunc = right.nfunc;
	nmo = right.nmo;
	nex = right.nex;
	charge = right.charge;
	ECP_m = right.ECP_m;
	multi = right.multi;
	origin = right.origin;
	total_energy = right.total_energy;
	virial_ratio = right.virial_ratio;
	basis_set_name = right.basis_set_name;
	comment = right.comment;
	path = right.path;
	method = right.method;
	MOs = right.MOs;
	centers = right.centers;
	types = right.types;
	exponents = right.exponents;
	UT_DensityMatrix = right.UT_DensityMatrix;
	UT_SpinDensityMatrix = right.UT_SpinDensityMatrix;
	DM = right.DM;
	DM_beta = right.DM_beta;
	MO_sph = right.MO_sph;
	basis_set = right.basis_set;
	cub = right.cub;
	fitted = right.fitted;
	atoms = right.atoms;
	modified = right.modified;
	d_f_switch = right.d_f_switch;
	distance_switch = right.distance_switch;
	has_ECPs = right.has_ECPs;
	isBohr = right.isBohr;
	is_unrestricted = right.is_unrestricted;
	fill_pre();
	fill_Afac_pre();
	return *this;
};

int WFN::calculate_charge()
{
	int atomic_charges = 0;
	double mo_charges = 0;
	for (int a = 0; a < ncen; a++)
	{
		int nr = get_atom_charge(a);
		if (nr == 0)
		{
			std::cout << "ERROR: Atomtype misunderstanding!\n";
			return -1000;
		}
		atomic_charges += nr;
	}
	for (int mo = 0; mo < nmo; mo++)
	{
		mo_charges += get_MO_occ(mo);
	}
	return atomic_charges - (int)mo_charges;
};

int WFN::calculate_charge(std::ostream &file)
{
	int atomic_charges = 0;
	double mo_charges = 0;
	for (int a = 0; a < ncen; a++)
	{
		int nr = get_atom_charge(a);
		if (nr == 0)
		{
			file << "ERROR: Atomtype misunderstanding!\n";
			return -1000;
		}
		atomic_charges += nr;
	}
	for (int mo = 0; mo < nmo; mo++)
	{
		mo_charges += get_MO_occ(mo);
	}
	return atomic_charges - (int)mo_charges;
};

bool WFN::guess_multiplicity(std::ostream &file)
{
	if (get_nr_electrons() % 2 == 0)
	{
		file << "With " << get_nr_electrons() << " electrons your system appears to have multiplicity 1." << std::endl;
		assign_multi(1);
	}
	else if (get_nr_electrons() % 2 == 1)
	{
		file << "With " << get_nr_electrons() % 2 << " electron this seems to be open shell, assuming multiplicity 2." << std::endl;
		assign_multi(2);
	}
	else
	{
		file << "This is awkward... i dind't think of this case yet, contact Florian to implement it!" << std::endl;
		return false;
	}
	return true;
};

bool WFN::push_back_cube(const std::string &filepath, const bool &full, const bool &expert)
{
	cub.emplace_back(filepath, full, *this, std::cout, expert);
	return true;
};

void WFN::pop_back_cube()
{
	cub.pop_back();
}

const unsigned int WFN::get_atom_integer_mass(const unsigned int &atomnr) const
{
	if (get_atom_charge(atomnr) > 86)
	{
		std::cout << "Sorry, only implemented until Rn yet, ask Florian for increases!" << std::endl;
		return 0;
	}
	if (get_atom_charge(atomnr) <= 0)
	{
		std::cout << "sorry, something seems wrong with the atoms you requested!" << std::endl;
		return 0;
	}
	return constants::integer_masses[get_atom_charge(atomnr) - 1];
};

const double WFN::get_atom_real_mass(const int &atomnr) const
{
	if (get_atom_charge(atomnr) > 86)
	{
		std::cout << "Sorry, only implemented until Xe yet, ask Florian for increases!" << std::endl;
		return 0;
	}
	if (get_atom_charge(atomnr) <= 0)
	{
		std::cout << "sorry, something seems wrong with the atoms you requested!" << std::endl;
		return 0;
	}
	return constants::real_masses[get_atom_charge(atomnr) - 1];
}

const double &WFN::get_MO_occ(const int &nr) const
{
	return MOs[nr].get_occ();
};

const int &WFN::get_MO_op(const int &nr) const
{
	return MOs[nr].get_op();
};

void WFN::delete_unoccupied_MOs()
{
	for (int i = static_cast<int>(MOs.size()) - 1; i >= 0; i--)
	{
		if (get_MO_occ(i) == 0.0)
		{
			MOs.erase(MOs.begin() + i);
			nmo--;
		}
	}
	invalidate_coef_cache();
};
void WFN::delete_Qs() {
	for (int i = static_cast<int>(atoms.size()) - 1; i >= 0; i--) {
		if (atoms[i].get_charge() == 119) {
			atoms.erase(atoms.begin() + i);
			ncen--;
			// centres are 1-based, i is the 0-based atom index: only atoms behind the dummy shift
			for (int j = 0; j < centers.size(); j++)
				if (centers[j] > i + 1)
					centers[j]--;
		}
	}
}

void WFN::pop_back_MO()
{
	MOs.pop_back();
	nmo--;
	invalidate_coef_cache();
}

bool WFN::delete_basis_set()
{
	for (int a = 0; a < get_ncen(); a++)
	{
		atoms[a].clear_shellcount();
		int nr_prim = get_atom_primitive_count(a);
		for (int p = 0; p < nr_prim; p++)
		{
			bool succes = erase_atom_primitive(a, 0);
			if (!succes)
				return false;
		}
	}
	return true;
};

std::string WFN::get_atom_label(const int &nr) const
{
	return atoms[nr].get_label();
};

int WFN::get_atom_ECP_electrons(const int &nr) const
{
	return atoms[nr].get_ECP_electrons();
};

basis_set_entry WFN::get_atom_basis_set_entry(const int &nr, const int &bs) const
{
	return atoms[nr].get_basis_set_entry(bs);
};

bool WFN::erase_atom_primitive(const unsigned int &nr, const unsigned int &nr_prim)
{
	if ((int)nr < ncen && (int)nr_prim < atoms[nr].get_basis_set_size())
	{
		atoms[nr].erase_basis_set(nr_prim);
		return true;
	}
	else
		return false;
};

std::filesystem::path WFN::get_cube_path(const int &nr) const
{
	return cub[nr].get_path();
};

void WFN::write_cube_file(const int &nr, const std::filesystem::path &filename, const bool &debug) {
	err_checkf(nr < cub.size(), "Wrong cube selected!", std::cout);
	if (cub[nr].get_path() != filename) {
		cub[nr].set_path(filename);
	}
	if (debug) {
		std::cout << "Writing cube file to: " << filename << std::endl;
	}
	cub[nr].write_file(filename);
	if (debug) {
		std::cout << "Cube file written!" << std::endl;
	}
}
void WFN::write_cube_xdgraph(const int &nr, const std::filesystem::path &filename, const bool &debug) {
	err_checkf(nr < cub.size(), "Wrong cube selected!", std::cout);
	if (cub[nr].get_path() != filename) {
		cub[nr].set_path(filename);
	}
	if (debug) {
		std::cout << "Writing cube file to: " << filename << std::endl;
	}
	cub[nr].write_xdgraph(filename);
	if (debug) {
		std::cout << "Cube file written!" << std::endl;
	}
}

void WFN::calc_rho_cube(cube &cube_data) const
{
	_time_point start = get_time();
	const int s1 = cube_data.get_size(0), s2 = cube_data.get_size(1), s3 = cube_data.get_size(2), total_size = s1 * s2 * s3;
	;
	std::cout << "Lets go into the loop! There is " << total_size << " points" << std::endl;
	ProgressBar *progress = new ProgressBar(total_size, 50, "=", " ", "Calculating Rho");

	vec v1{
		cube_data.get_vector(0, 0),
		cube_data.get_vector(1, 0),
		cube_data.get_vector(2, 0) },
		v2{
			cube_data.get_vector(0, 1),
			cube_data.get_vector(1, 1),
			cube_data.get_vector(2, 1) },
			v3{
				cube_data.get_vector(0, 2),
				cube_data.get_vector(1, 2),
				cube_data.get_vector(2, 2) },
				orig{
					cube_data.get_origin(0),
					cube_data.get_origin(1),
					cube_data.get_origin(2) };

#pragma omp parallel
	{
		vec2 d;
		vec phi((*this).get_nmo(), 0.0);
		d.resize((*this).get_ncen());
		for (int i = 0; i < (*this).get_ncen(); i++)
			d[i].resize(16, 0.0);
#pragma omp for schedule(dynamic)
		for (int index = 0; index < total_size; index++)
		{
			int i = index / (s2 * s3);
			int j = (index / s3) % s2;
			int k = index % s3;

			vec PosGrid{
				i * v1[0] + j * v2[0] + k * v3[0] + orig[0],
				i * v1[1] + j * v2[1] + k * v3[1] + orig[1],
				i * v1[2] + j * v2[2] + k * v3[2] + orig[2] };

			cube_data.set_value(i, j, k, (*this).compute_dens(cube_data.get_pos(i, j, k), d, phi));
			progress->update();
		}
	}
	delete (progress);

	using namespace std;
	_time_point end = get_time();
	if (get_sec(start, end) < 60)
		std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
	else if (get_sec(start, end) < 3600)
		std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
	else
		std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
	cube_data.calc_dv();
	std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << cube_data.sum() << std::endl;
};

void WFN::write_rho_cube(const double &radius, const double &increment) const {
	using namespace std;
	properties_options opts;
	opts.radius = radius;
	opts.resolution = increment;
	WFN dummy = (*this);
	readxyzMinMax_fromWFN(dummy, opts);
	cube CubeRho(opts.NbSteps, dummy.get_ncen(), true);
	dummy.delete_unoccupied_MOs();
	CubeRho.give_parent_wfn(dummy);
	std::cout << "Starting work..." << endl;

	for (int i = 0; i < 3; i++)
	{
		CubeRho.set_origin(i, opts.MinMax[i]);
		CubeRho.set_vector(i, i, (opts.MinMax[i + 3] - opts.MinMax[i]) / opts.NbSteps[i]);
	}
	CubeRho.set_comment1("Calculated density using NoSpherA2");
	CubeRho.set_comment2("from " + dummy.get_path().string());
	CubeRho.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho.cube");

	calc_rho_cube(CubeRho);

	std::filesystem::path fn((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho.cube");
	CubeRho.write_file(fn, false);
};
