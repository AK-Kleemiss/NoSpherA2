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
		//the shell loop below has cartesian norms for s, p, d and f only: a g shell would leave every later
		//shell's norm_const one shell off, so it is refused.  Generalising needs a g component order no
		//reader here agrees on.
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
				//axial_prim_norm is the general-l form of the s/p/d/f norms.  Types are validated above, serially:
				//err_checkf exits, and exiting from inside an OpenMP region does not.
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
