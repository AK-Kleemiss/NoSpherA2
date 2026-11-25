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

bool WFN::write_nbo(const std::filesystem::path &fileName, const bool &debug, std::ostream* progress_log, const std::string& nbo_keywords)
{
	using namespace std;

	//A FILE47 needs the contracted shells ($BASIS/$CONTRACT) and their AO overlap; a .wfn/.wfx carries primitives only and
	//a guessed archive gives plausible but wrong NBO output. Convert through .molden/.fchk/.gbw instead.
	err_checkf(get_nr_basis_set_loaded() == ncen,
		"Can only write a .47 file when a contracted basis set is present. A primitive-only source"
		" (.wfn/.wfx) does not carry one - use the .gbw, .fchk or .molden of the same calculation.",
		std::cout);
	const auto nbo_start_time = std::chrono::high_resolution_clock::now();
	auto progress_elapsed_seconds = [&]() {
		return get_sec(nbo_start_time, get_time());
	};
	auto progress = [&](const std::string& message) {
		if (progress_log != nullptr) {
			*progress_log << "[FILE47] " << message << " (" << progress_elapsed_seconds() << " s)" << std::endl;
			progress_log->flush();
		}
	};
	auto progress_percent = [&](const std::string& label, const int completed, const int total, int& next_percent) {
		if (progress_log == nullptr || total <= 0)
			return;
		const int percent = static_cast<int>((100LL * completed) / total);
		if (percent >= next_percent || completed == total) {
			*progress_log << "[FILE47] " << label << " " << percent << "% (" << completed << "/" << total << ", "
						  << progress_elapsed_seconds() << " s)" << std::endl;
			progress_log->flush();
			next_percent = percent + 10;
		}
	};

	progress("Starting .47 conversion to " + fileName.string());

	struct NboShell {
		int atom = 0;
		int shell = 0;
		int type = 0;
		int nprim = 0;
		int atom_prim_start = 0;
		int internal_start = 0;
		int nbo_start = 0;
		int cart_components = 0;
		int nbo_components = 0;
	};


	auto nbo_labels = [](const int type) -> ivec {
		switch (type) {
		case 1: return { 1 };
		case 2: return { 103, 101, 102 };
		case 3: return { 255, 252, 253, 254, 251 };
		case 4: return { 351, 352, 353, 354, 355, 356, 357 };
		case 5: return { 451, 452, 453, 454, 455, 456, 457, 458, 459 };
		default: return {};
		}
	};

	auto sph2cart_coefficient = [](const int type, const int cart, const int spher) {
		return constants::sph2cart(type - 1)[cart * constants::n_spher(type - 1) + spher];
	};

	auto project_to_nbo = [&](const int type, const vec& cart_values) {
		const int cart_count = constants::n_cart(type - 1);
		const int nbo_count = constants::n_spher(type - 1);
		vec2 normal(nbo_count, vec(nbo_count, 0.0));
		vec rhs(nbo_count, 0.0);
		for (int c = 0; c < cart_count; c++) {
			for (int i = 0; i < nbo_count; i++) {
				const double transform = sph2cart_coefficient(type, c, i);
				rhs[i] += transform * cart_values[c];
				for (int j = 0; j < nbo_count; j++)
					normal[i][j] += transform * sph2cart_coefficient(type, c, j);
			}
		}
		solve_linear_system(normal, rhs);
		return rhs;
	};

	vector<NboShell> shells;
	int nbo_nao = 0;
	int nbo_nexp = 0;
	int highest_angular = -1;
	for (int a = 0; a < get_ncen(); a++) {
		for (int s = 0; s < get_atom_shell_count(a); s++) {
			const int type = get_shell_type(a, s);
			const int cart_count = constants::n_cart(type - 1);
			const int nbo_count = constants::n_spher(type - 1);
			//FILE47 has h and i codes too; the ceiling is NoSpherA2's (n_cart, n_spher, sph2cart stop at g). Add CH/CI here if a basis ever needs h.
			err_checkf(type <= 5, "Unsupported basis shell in .47 writer: shells beyond g need"
				" constants::sph2cart extended first (FILE47 itself supports h and i)", std::cout);
			NboShell shell;
			shell.atom = a;
			shell.shell = s;
			shell.type = type;
			shell.nprim = get_atom_shell_primitives(a, s);
			shell.atom_prim_start = get_shell_start(a, s);
			shell.internal_start = 0;
			shell.nbo_start = nbo_nao;
			shell.cart_components = cart_count;
			shell.nbo_components = nbo_count;
			shells.push_back(shell);
			nbo_nao += nbo_count;
			nbo_nexp += shell.nprim;
			highest_angular = std::max(highest_angular, type - 1);
		}
	}
	//Int_Params, the gbw reader and OCC group an atom's shells by angular momentum, the loop
	//above follows the wavefunction's shell order. internal_start indexes the spherical
	//overlap, density and coefficient matrices in the grouped order, nbo_start the .47 order.
	{
		int internal_run = 0;
		for (int a = 0; a < get_ncen(); a++) {
			int max_l = 0;
			for (auto& sh : shells)
				if (sh.atom == a) max_l = std::max(max_l, sh.type - 1);
			for (int l = 0; l <= max_l; l++)
				for (auto& sh : shells)
					if (sh.atom == a && sh.type - 1 == l) {
						sh.internal_start = internal_run;
						internal_run += sh.nbo_components;
					}
		}
		err_checkf(internal_run == nbo_nao, "Angular-momentum regrouping lost basis functions in the .47 writer", std::cout);
	}
	progress("Built NBO shell model: atoms=" + std::to_string(get_ncen()) +
		", shells=" + std::to_string(shells.size()) +
		", NBO basis functions=" + std::to_string(nbo_nao) +
		", primitives=" + std::to_string(nbo_nexp));

	progress("Normalizing basis coefficients");
	vec2 basis_coefficients(get_ncen());
#pragma omp parallel for
	for (int a = 0; a < get_ncen(); a++)
	{
		for (int p = 0; p < get_atom_primitive_count(a); p++)
		{
			double temp_c = get_atom_basis_set_exponent(a, p);
			switch (get_atom_primitive_type(a, p))
			{
			case 1:
				temp_c = 8 * pow(temp_c, 3) / constants::PI3;
				break;
			case 2:
				temp_c = 128 * pow(temp_c, 5) / constants::PI3;
				break;
			case 3:
				temp_c = 2048 * pow(temp_c, 7) / (9 * constants::PI3);
				break;
			case 4:
				temp_c = 32768 * pow(temp_c, 9) / (225 * constants::PI3);
				break;
			case -1:
				std::cout << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy...\n";
				break;
			}
			temp_c = pow(temp_c, 0.25) * get_atom_basis_set_coefficient(a, p);
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
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * pow(constants::PI3 / pow(aiaj, 3), 0.5);
					}
				}
				if (factor == 0)
					return false;
				factor = pow(factor, -0.5);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					basis_coefficients[a][i] *= factor;
				}
				break;
			case 2:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * pow(constants::PI3 / (4 * pow(aiaj, 5)), 0.5);
					}
				}
				if (factor == 0)
					return false;
				factor = pow(factor, -0.5);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					basis_coefficients[a][i] *= factor;
				}
				break;
			case 3:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * pow(constants::PI3 / (16 * pow(aiaj, 7)), 0.5);
					}
				}
				if (factor == 0)
					return false;
				factor = (pow(factor, -0.5)) / sqrt(3);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					basis_coefficients[a][i] *= factor;
				}
				break;
			case 4:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * pow(constants::PI3 / (64 * pow((aiaj), 9)), 0.5);
					}
				}
				if (factor == 0)
					return false;
				factor = pow(factor, -0.5) / sqrt(15);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					basis_coefficients[a][i] *= factor;
				}
				break;
			}
			if (debug)
				std::cout << "This shell has: " << get_shell_end(a, s) - get_shell_start(a, s) + 1 << " primitives\n";
		}
	}
	//NBO slot c of a shell (m = 0, +1, -1, +2, -2, ... as ORCA orders it) holds pure function
	//orca_2_pySCF(l, c) of the m = -l..l order libcint, OCC and the gbw reader keep their
	//matrices in; sph_index maps a slot to its row in those matrices.
	ivec sph_index(nbo_nao, 0);
	for (const auto& shell : shells)
		for (int c = 0; c < shell.nbo_components; c++) {
			const auto offset = constants::orca_2_pySCF(shell.type - 1, c);
			err_checkf(offset.has_value(), "Unsupported NBO AO order in .47 writer", std::cout);
			sph_index[shell.nbo_start + c] = shell.internal_start + static_cast<int>(offset.value());
		}

	const int alpha_mos = get_MO_op_count(0);
	const int beta_mos = get_MO_op_count(1);
	vec2 CMO(alpha_mos, vec(nbo_nao, 0.0));
	vec2 CMO_beta(beta_mos, vec(nbo_nao, 0.0));
	//An XCW_fit wavefunction still holds the spherical coefficients OCC converged; rebuilding
	//them from the primitives below needs OCC's normalization convention, which the WFN does
	//not carry. Other origins that cache their coefficients could take this path as well.
	const int spin_blocks = get_is_unrestricted() ? 2 : 1;
	const bool cached_mos = get_origin() == e_origin::XCW_fit
		&& static_cast<int>(MO_sph.extent(0)) == spin_blocks * nbo_nao
		&& static_cast<int>(MO_sph.extent(1)) >= std::max(alpha_mos, beta_mos);
	if (cached_mos) {
		progress("Taking MO coefficients from the cached spherical coefficient matrix");
		for (int m = 0; m < alpha_mos; m++)
			for (int i = 0; i < nbo_nao; i++)
				CMO[m][i] = MO_sph(sph_index[i], m);
		for (int m = 0; m < beta_mos; m++)
			for (int i = 0; i < nbo_nao; i++)
				CMO_beta[m][i] = MO_sph(nbo_nao + sph_index[i], m);
	}
	else
		progress("Reconstructing MO coefficients in NBO AO order: alpha=" + std::to_string(alpha_mos) +
			", beta=" + std::to_string(beta_mos));
	int alpha_run = 0;
	int beta_run = 0;
	//the sph2cart tables return |m| = 3, 4 with ORCA's sign; fchk and OCC coefficients, their DM and the overlap keep libcint's
	const bool libcint_phase = !origin_has_orca_pure_phases(get_origin());
	int mo_progress_next = 10;
	int mo_seen = 0;
	for (int m = 0; m < get_nmo() && !cached_mos; m++)
	{
		if (MOs[m].get_op() != 0 && MOs[m].get_op() != 1)
			continue;
		mo_seen++;
		vec* target = nullptr;
		if (MOs[m].get_op() == 0) {
			target = &CMO[alpha_run++];
		}
		else {
			target = &CMO_beta[beta_run++];
		}
		for (const auto& shell : shells) {
			vec cart_values(shell.cart_components, 0.0);
			const int primitive_start = get_shell_start_in_primitives(shell.atom, shell.shell);
			//Every primitive of a contracted shell is AO coefficient times its stored contraction
			//coefficient, so the largest one is the safest divisor. The gbw reader stores a
			//shell primitive-major (x, y, z of primitive 0, then of primitive 1), the OCC
			//constructor component-major (all primitives of x, then of y); equal consecutive
			//types tell the two apart.
			const int shell_start = get_shell_start(shell.atom, shell.shell);
			int rep = 0;
			for (int p = 1; p < shell.nprim; p++)
				if (std::abs(get_atom_basis_set_coefficient(shell.atom, shell_start + p)) >
					std::abs(get_atom_basis_set_coefficient(shell.atom, shell_start + rep)))
					rep = p;
			const double contraction = get_atom_basis_set_coefficient(shell.atom, shell_start + rep);
			err_checkf(std::abs(contraction) > 1E-14, "Cannot reconstruct pure MO coefficients for .47 output", std::cout);
			const bool component_major = shell.nprim > 1 && shell.cart_components > 1
				&& get_type(primitive_start) == get_type(primitive_start + 1);
			const int stride = component_major ? shell.nprim : 1;
			const int rep_offset = component_major ? rep : rep * shell.cart_components;
			//the components may be permuted as well (gbw stores p as z, x, y): the type says which row
			for (int c = 0; c < shell.cart_components; c++) {
				const int prim = primitive_start + c * stride + rep_offset;
				//A foreign wavefunction decides this index: prim is component c only if the detected primitive order holds
				const int component = get_type(prim) - constants::first_type[shell.type - 1];
				err_checkf(component >= 0 && component < shell.cart_components,
					"Primitive order not understood in the .47 writer: atom " + std::to_string(shell.atom + 1)
					+ ", shell " + std::to_string(shell.shell) + " of type " + std::to_string(shell.type)
					+ ", component " + std::to_string(c) + " -> primitive " + std::to_string(prim)
					+ " of type " + std::to_string(get_type(prim)) + " (component index " + std::to_string(component)
					+ " of " + std::to_string(shell.cart_components) + ")", std::cout);
				cart_values[component] = get_MO_coef(m, prim) / contraction;
			}
			const vec nbo_values = project_to_nbo(shell.type, cart_values);
			for (int c = 0; c < shell.nbo_components; c++)
				(*target)[shell.nbo_start + c] = libcint_phase && orca_pure_sign_flips((c + 1) / 2) ? -nbo_values[c] : nbo_values[c];
		}
		progress_percent("MO coefficient reconstruction", mo_seen, alpha_mos + beta_mos, mo_progress_next);
	}
	progress("Building density matrix");
	int naotr = nbo_nao * (nbo_nao + 1) / 2;
	vec CDM(naotr, 0.0);
	//Occupations and energies per spin, in the row order of CMO / CMO_beta. An open-shell
	//FILE47 carries $DENSITY, $FOCK and $LCAOMO twice - alpha block then beta block - while
	//$OVERLAP stays single; a closed-shell one carries the spin sum once.
	vec occ_spin[2], energy_spin[2];
	for (int m = 0; m < get_nmo(); m++) {
		const int op = MOs[m].get_op();
		if (op != 0 && op != 1)
			continue;
		occ_spin[op].push_back(get_MO_occ(m));
		energy_spin[op].push_back(get_MO_energy(m));
	}
	const bool open_shell = get_is_unrestricted() && beta_mos > 0;
	auto build_mo_density = [&](const int op) {
		const vec2& C = op == 0 ? CMO : CMO_beta;
		const vec& occs = occ_spin[op];
		const int nmo_spin = std::min(static_cast<int>(C.size()), static_cast<int>(occs.size()));
		vec density(naotr, 0.0);
		int density_progress_next = 10;
#pragma omp parallel for schedule(dynamic)
		for (int iu = 0; iu < nbo_nao; iu++) {
			for (int iv = 0; iv <= iu; iv++) {
				const int iuv = (iu * (iu + 1) / 2) + iv;
				for (int m = 0; m < nmo_spin; m++)
					if (occs[m] != 0.0)
						density[iuv] += occs[m] * C[m][iu] * C[m][iv];
			}
#pragma omp critical(nbo_progress)
			progress_percent("Density build", iu + 1, nbo_nao, density_progress_next);
		}
		return density;
	};
	auto build_total_density = [&]() {
		vec density = build_mo_density(0);
		if (beta_mos > 0) {
			const vec beta_density = build_mo_density(1);
			for (int i = 0; i < naotr; i++)
				density[i] += beta_density[i];
		}
		return density;
	};
	vec CDM_alpha, CDM_beta;
	bool density_from_cached_dm = false;
	if (open_shell) {
		progress("Building spin-resolved density matrices for the open-shell FILE47");
		CDM_alpha = build_mo_density(0);
		CDM_beta = build_mo_density(1);
		for (int i = 0; i < naotr; i++)
			CDM[i] = CDM_alpha[i] + CDM_beta[i];
	}
	//A cached DM is the spin sum, so it can only serve the closed-shell layout.
	else if (static_cast<int>(DM.extent(0)) == nbo_nao && static_cast<int>(DM.extent(1)) == nbo_nao) {
		density_from_cached_dm = true;
		for (int iu = 0; iu < nbo_nao; iu++) {
			for (int iv = 0; iv <= iu; iv++) {
				const int iuv = (iu * (iu + 1) / 2) + iv;
				CDM[iuv] = DM(sph_index[iu], sph_index[iv]);
			}
		}
	}
	else {
		CDM = build_total_density();
	}

	vec OVLP_matrix = {};
	//Int_Params reads each shell's angular momentum by origin (OCC and NOT_YET_DEFINED
	//store l, everything else l + 1) and normalises by origin; an XCW_fit wavefunction
	//says which convention it has, so it is handed over as it is.
	Int_Params int_params(*this);
	progress("Computing spherical AO overlap integrals");
	compute2C<Overlap2C_SPH>(int_params, OVLP_matrix);
	err_checkf(static_cast<int>(OVLP_matrix.size()) == nbo_nao * nbo_nao, "Spherical overlap has the wrong size in the .47 writer", std::cout);
	dMatrixRef2 OVLP_sph(OVLP_matrix.data(), nbo_nao, nbo_nao);
	//libcint and OCC use the standard phases; ORCA's f(+-3), g(+-3), g(+-4) have the opposite
	//sign and the gbw reader keeps them, so the overlap takes ORCA's sign there. Other origins
	//with ORCA-like conventions need the same: origin_has_orca_pure_phases names them, and a molden
	//written by orca_2mkl is one.
	vec phase(nbo_nao, 1.0);
	if (origin_has_orca_pure_phases(get_origin()))
		for (const auto& shell : shells)
			for (int c = 0; c < shell.nbo_components; c++)
				if (orca_pure_sign_flips((c + 1) / 2))
					phase[shell.nbo_start + c] = -1.0;
	vec2 OVLP_nbo(nbo_nao, vec(nbo_nao, 0.0));
	for (int i = 0; i < nbo_nao; i++)
		for (int j = 0; j < nbo_nao; j++)
			OVLP_nbo[i][j] = phase[i] * phase[j] * OVLP_sph(sph_index[i], sph_index[j]);
	//A diagonal off 1 means Int_Params normalised this origin with the wrong convention; the
	//density beside it is in a normalised basis, so rescaling S alone would only hide it.
	int unnormalised = 0;
	for (int i = 0; i < nbo_nao; i++)
		if (std::abs(OVLP_nbo[i][i] - 1.0) > 1e-8) unnormalised++;
	if (unnormalised > 0)
		progress(std::to_string(unnormalised) + " of " + std::to_string(nbo_nao)
			+ " AOs are not normalised in the .47 overlap; check the origin's normalisation in Int_Params");
	auto packed_trace_product = [&](const vec& density) {
		double trace = 0.0;
		for (int i = 0; i < nbo_nao; i++) {
			for (int j = 0; j <= i; j++) {
				const int ij = (i * (i + 1) / 2) + j;
				trace += density[ij] * OVLP_nbo[i][j] * (i == j ? 1.0 : 2.0);
			}
		}
		return trace;
	};
	double expected_electrons = 0.0;
	for (int m = 0; m < get_nmo(); m++)
		expected_electrons += get_MO_occ(m);
	double density_electrons = packed_trace_product(CDM);
	if (density_from_cached_dm && std::abs(density_electrons - expected_electrons) > 1.0E-4) {
		progress("Cached GBW density is inconsistent with FILE47 overlap: Tr(P*S)=" +
			std::to_string(density_electrons) + ", expected=" + std::to_string(expected_electrons) +
			". Rebuilding density from NBO-ordered MO coefficients");
		CDM = build_total_density();
		density_from_cached_dm = false;
		density_electrons = packed_trace_product(CDM);
	}
	if (progress_log != nullptr) {
		std::ostringstream density_check;
		density_check << std::fixed << std::setprecision(6) << density_electrons
					  << ", expected=" << expected_electrons;
		*progress_log << "[FILE47] Density electron check Tr(P*S)=" << density_check.str()
					  << " (" << progress_elapsed_seconds() << " s)" << std::endl;
		progress_log->flush();
	}
	//Tr(P*S) decides whether the archive describes the wavefunction at all, and every quantity downstream (NAO populations,
	//NBO occupancies, E2, NRT) is a functional of this density, so it is checked on every path and fatal.
	if (expected_electrons > 0.0 &&
		std::abs(density_electrons - expected_electrons) > 1.0E-4 * std::max(1.0, expected_electrons)) {
		std::ostringstream why;
		why << std::fixed << std::setprecision(6)
			<< "The FILE47 archive does not describe this wavefunction: Tr(P*S)=" << density_electrons
			<< " electrons where the occupied orbitals hold " << expected_electrons << ".";
		if (unnormalised > 0)
			why << " " << unnormalised << " of " << nbo_nao << " AOs are not normalised in the"
				<< " overlap either, so the contracted basis this source was read into carries the"
				<< " wrong normalisation convention (see Int_Params).";
		why << " Every NBO quantity is a functional of this density, so no analysis is written."
			<< " Use the .gbw or .molden of the same calculation.";
		err_checkf(false, why.str(), std::cout);
	}
	auto build_fock = [&](const vec2& C, const vec& energies, const std::string& label) {
		vec2 result;
		if (static_cast<int>(C.size()) != nbo_nao || static_cast<int>(energies.size()) < nbo_nao) {
			progress("Skipping " + label + " Fock matrix: MO count does not match NBO basis size");
			return result;
		}
		progress("Building " + label + " Fock matrix from MO energies with BLAS");
		result = vec2(nbo_nao, vec(nbo_nao, 0.0));
		dMatrix2 overlap(nbo_nao, nbo_nao);
		dMatrix2 cmo(nbo_nao, nbo_nao);
		for (int i = 0; i < nbo_nao; i++) {
			for (int j = 0; j < nbo_nao; j++)
				overlap(i, j) = OVLP_nbo[i][j];
			for (int m = 0; m < nbo_nao; m++)
				cmo(m, i) = C[m][i];
		}
		dMatrix2 eps_cmo(nbo_nao, nbo_nao);
#pragma omp parallel for schedule(dynamic)
		for (int m = 0; m < nbo_nao; m++)
			for (int i = 0; i < nbo_nao; i++)
				eps_cmo(m, i) = energies[m] * cmo(m, i);
		progress("Fock build: C^T * eps * C");
		dMatrix2 ctc = dot<dMatrix2>(cmo, eps_cmo, true, false);
		progress("Fock build: S * (C^T * eps * C)");
		dMatrix2 left = dot<dMatrix2>(overlap, ctc, false, false);
		progress("Fock build: S * (C^T * eps * C) * S");
		dMatrix2 fock = dot<dMatrix2>(left, overlap, false, false);
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nbo_nao; i++)
			for (int j = 0; j < nbo_nao; j++)
				result[i][j] = fock(i, j);
		return result;
	};
	vec2 FOCK_nbo = build_fock(CMO, energy_spin[0], open_shell ? "alpha" : "total");
	vec2 FOCK_beta;
	if (open_shell) {
		FOCK_beta = build_fock(CMO_beta, energy_spin[1], "beta");
		//An open-shell $FOCK is read as two blocks; one alone would be parsed as the alpha
		//block and leave NBO reading the next section as beta, so it is both or neither.
		if (FOCK_nbo.empty() || FOCK_beta.empty()) {
			FOCK_nbo.clear();
			FOCK_beta.clear();
			progress("Skipping $FOCK entirely: an open-shell FILE47 needs both spin blocks");
		}
	}

	ofstream rf(fileName, ios::out);
	if (!rf.is_open())
	{
		std::cout << "Sorry, can't open the file...\n";
		return false;
	}
	progress("Writing FILE47 sections");

	auto write_int_array = [&](const string& name, const ivec& values, const int per_line, const int continuation_indent) {
		rf << name;
		for (int i = 0; i < values.size(); i++) {
			if (i > 0 && i % per_line == 0)
				rf << "\n" << string(continuation_indent, ' ');
			rf << setw(5) << values[i];
		}
		rf << endl;
	};

	auto format_precise_value = [](const double value) {
		if (value == 0.0)
			return string("  0.000000000000E+00");

		const double abs_value = std::abs(value);
		const int exponent = static_cast<int>(std::floor(std::log10(abs_value))) + 1;
		double mantissa = value / std::pow(10.0, exponent);
		int adjusted_exponent = exponent;
		if (std::abs(mantissa) >= 1.0) {
			mantissa /= 10.0;
			adjusted_exponent++;
		}

		stringstream local;
		local << fixed << setprecision(12) << mantissa;
		if (local.str().rfind("1.", 0) == 0 || local.str().rfind("-1.", 0) == 0) {
			mantissa /= 10.0;
			adjusted_exponent++;
			local.str("");
			local.clear();
			local << fixed << setprecision(12) << mantissa;
		}
		string result = local.str() + "E" + (adjusted_exponent >= 0 ? "+" : "-");
		result += (std::abs(adjusted_exponent) < 10 ? "0" : "") + std::to_string(std::abs(adjusted_exponent));
		return string(std::max(0, 20 - static_cast<int>(result.size())), ' ') + result;
	};

	auto write_real_array = [&](const string& name, const vec& values) {
		rf << name;
		for (int i = 0; i < values.size(); i++) {
			if (i > 0 && i % 3 == 0)
				rf << "\n" << string(10, ' ');
			rf << format_precise_value(values[i]);
		}
		rf << endl;
	};

	rf << " $GENNBO NATOMS=" << ncen << " NBAS=" << nbo_nao << (open_shell ? " OPEN" : "") << " UPPER BODM FORMAT=PRECISE $END" << endl;
	rf << " $NBO" << (nbo_keywords.empty() ? "" : " " + nbo_keywords) << " $END" << endl;
	rf << " $COORD" << endl;
	rf << " .47 file generated by NoSpherA2 based on " << path << endl;
	for (int i = 0; i < ncen; i++)
	{
		const int atomic_number = get_atom_charge(i);
		const int nuclear_charge = atomic_number - get_atom_ECP_electrons(i);
		rf << " " << setw(5) << atomic_number << setw(5) << nuclear_charge
			<< fixed << setprecision(6) << setw(15) << constants::bohr2ang(get_atom_coordinate(i, 0))
			<< fixed << setprecision(6) << setw(15) << constants::bohr2ang(get_atom_coordinate(i, 1))
			<< fixed << setprecision(6) << setw(15) << constants::bohr2ang(get_atom_coordinate(i, 2)) << "\n";
	}
	rf << " $END" << endl;

	ivec nbo_centers;
	ivec labels;
	for (const auto& shell : shells) {
		const ivec shell_labels = nbo_labels(shell.type);
		for (int label : shell_labels) {
			nbo_centers.push_back(shell.atom + 1);
			labels.push_back(label);
		}
	}
	rf << " $BASIS" << endl;
	write_int_array("  CENTER =", nbo_centers, 11, 10);
	write_int_array("   LABEL =", labels, 11, 10);
	rf << " $END" << endl;

	ivec ncomp;
	ivec nprim;
	ivec nptr;
	vec exp_values;
	vec2 contract_values(5, vec(nbo_nexp, 0.0));
	int exp_index = 0;
	for (const auto& shell : shells) {
		ncomp.push_back(shell.nbo_components);
		nprim.push_back(shell.nprim);
		nptr.push_back(exp_index + 1);
		for (int p = 0; p < shell.nprim; p++) {
			exp_values.push_back(get_atom_basis_set_exponent(shell.atom, shell.atom_prim_start + p));
			contract_values[shell.type - 1][exp_index] = get_atom_basis_set_coefficient(shell.atom, shell.atom_prim_start + p);
			exp_index++;
		}
	}

	rf << " $CONTRACT" << endl;
	rf << "  NSHELL =" << setw(6) << shells.size() << endl;
	rf << "    NEXP =" << setw(6) << nbo_nexp << endl;
	write_int_array("   NCOMP =", ncomp, 11, 10);
	write_int_array("   NPRIM =", nprim, 11, 10);
	write_int_array("    NPTR =", nptr, 11, 10);
	write_real_array("     EXP =", exp_values);
	if (highest_angular >= 0) write_real_array("      CS =", contract_values[0]);
	if (highest_angular >= 1) write_real_array("      CP =", contract_values[1]);
	if (highest_angular >= 2) write_real_array("      CD =", contract_values[2]);
	if (highest_angular >= 3) write_real_array("      CF =", contract_values[3]);
	if (highest_angular >= 4) write_real_array("      CG =", contract_values[4]);
	rf << " $END" << endl;
	auto write_precise_value = [&](const double value, int& count) {
		rf << format_precise_value(value);
		count++;
		if (count % 4 == 0)
			rf << "\n";
	};
	//NBO reads each spin block with its own Fortran READ, so the beta block starts a fresh record; otherwise TINP fails on
	//$LCAOMO whenever a block length is not a multiple of four.
	auto end_block = [&](int& count) {
		if (count % 4 != 0)
			rf << "\n";
		count = 0;
	};

	rf << " $OVERLAP" << endl;
	int runner = 0;
	for (int i = 0; i < nbo_nao; i++)
		for (int j = 0; j <= i; j++)
			write_precise_value(OVLP_nbo[i][j], runner);
	end_block(runner);
	rf << " $END" << endl;
	rf << " $DENSITY" << endl;
	for (int i = 0; i < naotr; i++)
		write_precise_value(open_shell ? CDM_alpha[i] : CDM[i], runner);
	end_block(runner);
	if (open_shell) {
		for (int i = 0; i < naotr; i++)
			write_precise_value(CDM_beta[i], runner);
		end_block(runner);
	}
	rf << " $END" << endl;
	if (!FOCK_nbo.empty()) {
		auto write_fock_block = [&](const vec2& fock) {
			for (int i = 0; i < nbo_nao; i++)
				for (int j = 0; j <= i; j++)
					write_precise_value(fock[i][j], runner);
			end_block(runner);
		};
		rf << " $FOCK" << endl;
		write_fock_block(FOCK_nbo);
		if (open_shell)
			write_fock_block(FOCK_beta);
		rf << " $END" << endl;
	}
	//$LCAOMO is nbas x nbas per spin block whatever the MO count, so fewer MOs than basis functions are zero-padded
	auto write_lcaomo_block = [&](const vec2& C) {
		for (int mo_counter = 0; mo_counter < nbo_nao; mo_counter++)
		{
			if (debug)
				std::cout << "Writing MO #" << mo_counter + 1 << "...\n";
			for (int i = 0; i < nbo_nao; i++)
				write_precise_value(mo_counter < static_cast<int>(C.size()) ? C[mo_counter][i] : 0.0, runner);
		}
		end_block(runner);
	};
	rf << " $LCAOMO" << endl;
	write_lcaomo_block(CMO);
	if (open_shell)
		write_lcaomo_block(CMO_beta);
	rf << " $END" << endl;
	rf.close();
	progress("Finished .47 conversion");
	return true;
};
