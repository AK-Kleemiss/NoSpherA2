#include "pch.h"
#include "wfn_class.h"
#include "properties.h"
#ifdef NOSPHERA2_USE_GPU
#include "aux_density_gpu.h"
#include "spherical_density_gpu.h"
#endif
#include "convenience.h"
#include "spherical_density.h"
#include "cube.h"
#include "constants.h"
#include "ibo.h"
#include "nbo.h"
#include "GridManager.h"
#include "isosurface.h"
#include "SALTED_predictor.h"
#include "crystal_energies.h"
#include "b2c.h"
#include "density_source.h"
#include "citations.h"
#include "eli_family.h"

static void Calc_Hirshfeld_atom_2(
	cube &CubeHirsh,
	cube &CubeRho,
	cube &CubeSpherical,
	WFN &wavy,
	int _atom,
	std::ostream &file)
{
	(void)file;
	using namespace std;
	Thakkar atom(wavy.get_atom_charge(_atom));

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < CubeHirsh.get_size(0); i++)
	{
		for (int j = 0; j < CubeHirsh.get_size(1); j++)
			for (int k = 0; k < CubeHirsh.get_size(2); k++)
			{

				const d3 PosGrid{ i * CubeHirsh.get_vector(0, 0) + j * CubeHirsh.get_vector(0, 1) + k * CubeHirsh.get_vector(0, 2) + CubeHirsh.get_origin(0),
										i * CubeHirsh.get_vector(1, 0) + j * CubeHirsh.get_vector(1, 1) + k * CubeHirsh.get_vector(1, 2) + CubeHirsh.get_origin(1),
										i * CubeHirsh.get_vector(2, 0) + j * CubeHirsh.get_vector(2, 1) + k * CubeHirsh.get_vector(2, 2) + CubeHirsh.get_origin(2) };

				// bool skip = true;

				const double dens_choice = atom.get_radial_density(array_length(PosGrid, wavy.get_atom_pos(_atom)));
				const double temp_val = CubeSpherical.get_value(i, j, k);
				if (temp_val != 0)
					CubeHirsh.set_value(i, j, k, (dens_choice / temp_val * CubeRho.get_value(i, j, k)));
			}
	}
};

enum class dipole_types
{
	atom,
	geometry,
	hirshfeld,
	vdW,
	Unknown
};

dipole_types stringTodipole_types(const std::string &str)
{
	static const std::unordered_map<std::string, dipole_types> stringToEnumMap = {
		{"atom", dipole_types::atom},
		{"geometry", dipole_types::geometry},
		{"hirshfeld", dipole_types::hirshfeld},
		{"vdW", dipole_types::vdW} };

	auto it = stringToEnumMap.find(str);
	if (it != stringToEnumMap.end())
	{
		return it->second;
	}
	else
	{
		return dipole_types::Unknown;
	}
}

vec calc_dipole_for_atom(WFN &wavy, const int &i, cube &Hirshfeld_atom, vec &charges, std::string type = "atom")
{
	double mu_x = 0, mu_y = 0, mu_z = 0;
	double scratch = 0;
	const d3 ax = wavy.get_atom_pos(i);
	double dv = Hirshfeld_atom.get_dv();
	// const int c = wavy.get_atom_charge(i);
	double charge = 0;
	vec origin{ 0, 0, 0 };
	vec bound_atoms;
	for (int j = 0; j < wavy.get_ncen(); j++)
	{
		if (i == j)
			continue;
		const double dist = array_length(ax, wavy.get_atom_pos(j));
		const double svdW = constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)];
		if (dist < 1.1 * svdW)
		{
			bound_atoms.push_back(j);
		}
	}
	const double v[9] = { Hirshfeld_atom.get_vector(0, 0), Hirshfeld_atom.get_vector(0, 1), Hirshfeld_atom.get_vector(0, 2),
						 Hirshfeld_atom.get_vector(1, 0), Hirshfeld_atom.get_vector(1, 1), Hirshfeld_atom.get_vector(1, 2),
						 Hirshfeld_atom.get_vector(2, 0), Hirshfeld_atom.get_vector(2, 1), Hirshfeld_atom.get_vector(2, 2) };
	switch (stringTodipole_types(type))
	{
	case dipole_types::atom:
		origin = { ax[0], ax[1], ax[2] };
		break;
	case dipole_types::geometry:
		err_not_impl_f("geometry position not yet implemented", std::cout);
		origin = { 0, 0, 0 };
		break;
	case dipole_types::hirshfeld:
		err_not_impl_f("hirshfeld centers not yet implemented", std::cout);
		break;
	case dipole_types::vdW:
		err_not_impl_f("vdW radius basis not implemented", std::cout);
		break;
	case dipole_types::Unknown:
		err_not_impl_f("Unknown dipole type", std::cout);
		break;
	}
#pragma omp parallel for reduction(+ : mu_x, mu_y, mu_z, charge) private(scratch)
	for (int x = 0; x < Hirshfeld_atom.get_size(0); x++)
	{
		for (int y = 0; y < Hirshfeld_atom.get_size(1); y++)
		{
			for (int z = 0; z < Hirshfeld_atom.get_size(2); z++)
			{
				const d3 PosGrid{
					x * v[0] + y * v[1] + z * v[2] + Hirshfeld_atom.get_origin(0),
					x * v[3] + y * v[4] + z * v[5] + Hirshfeld_atom.get_origin(1),
					x * v[6] + y * v[7] + z * v[8] + Hirshfeld_atom.get_origin(2) };
				scratch = Hirshfeld_atom.get_value(x, y, z) * dv;
				charge += scratch;
				mu_x += (PosGrid[0] - origin[0]) * scratch;
				mu_y += (PosGrid[1] - origin[1]) * scratch;
				mu_z += (PosGrid[2] - origin[2]) * scratch;
			}
		}
	}
	return { mu_x, mu_y, mu_z, charge };
}

void dipole_moments(options &opt, std::ostream &log2)
{
	using namespace std;
	//no header here: run_app_impl has written it to the log before the options were read
	err_checkf(opt.wfn != "", "Error, no wfn file specified!", log2);
	WFN wavy(opt.wfn);
	if (opt.debug)
		log2 << "Starting calculation of dipole moment" << endl;

	if (opt.debug)
		log2 << opt.cif << " " << opt.properties.resolution << " " << opt.properties.radius << endl;
	readxyzMinMax_fromWFN(wavy, opt.properties);
	if (opt.debug)
	{
		log2 << "Resolution: " << opt.properties.resolution << endl;
		log2 << "MinMax:" << endl;
		for (int i = 0; i < 6; i++)
			log2 << setw(14) << scientific << opt.properties.MinMax[i];
		log2 << endl;
		log2 << "Steps:" << endl;
		for (int i = 0; i < 3; i++)
			log2 << setw(14) << scientific << opt.properties.NbSteps[i];
		log2 << endl;
	}
	cube Rho(opt.properties.NbSteps, wavy.get_ncen(), true);
	cube SPHER(opt.properties.NbSteps, wavy.get_ncen(), true);

	Rho.give_parent_wfn(wavy);
	SPHER.give_parent_wfn(wavy);
	vec stepsizes{ (opt.properties.MinMax[3] - opt.properties.MinMax[0]) / opt.properties.NbSteps[0],
				  (opt.properties.MinMax[4] - opt.properties.MinMax[1]) / opt.properties.NbSteps[1],
				  (opt.properties.MinMax[5] - opt.properties.MinMax[2]) / opt.properties.NbSteps[2] };

	for (int i = 0; i < 3; i++)
	{
		Rho.set_origin(i, opt.properties.MinMax[i]);
		SPHER.set_origin(i, opt.properties.MinMax[i]);
		Rho.set_vector(i, i, stepsizes[i]);
		SPHER.set_vector(i, i, stepsizes[i]);
	}
	if (opt.debug)
		log2 << "Origins etc are set up" << endl;
	Rho.set_comment1("Calculated density using NoSpherA2");
	SPHER.set_comment1("Calculated Atomic Hirshfeld deformation density values using NoSpherA2");
	Rho.set_comment2("from " + wavy.get_path().string());
	SPHER.set_comment2("from" + wavy.get_path().string());
	Rho.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
	SPHER.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_spher.cube");
	vector<cube> Hirsh(wavy.get_ncen(), Rho);
	vec charges(wavy.get_ncen(), 0);

	log2 << "Calculating for " << fixed << setprecision(0) << opt.properties.NbSteps[0] * opt.properties.NbSteps[1] * opt.properties.NbSteps[2] << " Gridpoints." << endl;

	log2 << "Calcualting Rho...";
	Calc_Rho(Rho, wavy, opt.properties.radius, log2, false);
	log2 << " ...done!\nCalcualting spherical Rho...";
	Calc_Spherical_Dens(SPHER, wavy, opt.properties.radius, log2, false);
	log2 << " ...done!" << endl;
	vec2 dipole_moments;
	dipole_moments.reserve(wavy.get_ncen());
	for (int i = 0; i < wavy.get_ncen(); i++)
	{
		Hirsh[i].calc_dv();
		Hirsh[i].give_parent_wfn(wavy);
		Hirsh[i].set_zero();
		log2 << "Calcualting Hirshfeld density for atom: " << i << endl;
		Calc_Hirshfeld_atom_2(Hirsh[i], Rho, SPHER, wavy, i, log2);
		charges[i] = Hirsh[i].sum();
		log2 << "..done!" << endl;
	}
	for (int i = 0; i < wavy.get_ncen(); i++)
		dipole_moments.emplace_back(calc_dipole_for_atom(wavy, i, Hirsh[i], charges));
	log2 << " atom   |  dipole moment x,        y,         z" << endl
		<< "======================================" << endl;
	for (int i = 0; i < wavy.get_ncen(); i++)
	{
		log2 << setw(3) << i << " (" << constants::atnr2letter(wavy.get_atom_charge(i)) << ") | " << scientific << setprecision(6) << setw(14) << dipole_moments[i][0] << ", " << setw(14) << dipole_moments[i][1] << ", " << setw(14) << dipole_moments[i][2] << "\n";
	}
	std::cout << "\n\nProperties calculation done!" << std::endl;
}

vec2 dipole_moments(WFN &wavy, cube &SPHER, const properties_options &opts, int threads, std::ostream &log2, bool debug)
{
	using namespace std;
	if (debug)
		log2 << "Starting calculation of dipole moment" << endl;
	cube Rho(opts.NbSteps, wavy.get_ncen(), true);

	Rho.give_parent_wfn(wavy);
	vec stepsizes{ (opts.MinMax[3] - opts.MinMax[0]) / opts.NbSteps[0],
				  (opts.MinMax[4] - opts.MinMax[1]) / opts.NbSteps[1],
				  (opts.MinMax[5] - opts.MinMax[2]) / opts.NbSteps[2] };

	for (int i = 0; i < 3; i++)
	{
		Rho.set_origin(i, opts.MinMax[i]);
		Rho.set_vector(i, i, stepsizes[i]);
	}
	if (debug)
		log2 << "Origins etc are set up" << endl;
	Rho.set_comment1("Calculated density using NoSpherA2");
	Rho.set_comment2("from " + wavy.get_path().string());
	Rho.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
	vector<cube> Hirsh(wavy.get_ncen(), Rho);
	vec charges(wavy.get_ncen(), 0);

	log2 << "Calculating for " << fixed << setprecision(0) << opts.n_grid_points() << " Gridpoints." << endl;

	log2 << "Calcualting Rho...";
	Calc_Rho(Rho, wavy, opts.radius, log2, false);
	log2 << " ...done!\nCalcualting spherical Rho...";
	Calc_Spherical_Dens(SPHER, wavy, opts.radius, log2, false);
	log2 << " ...done!" << endl;
	vec2 dipole_moments;
	for (int i = 0; i < wavy.get_ncen(); i++)
	{
		Hirsh[i].calc_dv();
		Hirsh[i].give_parent_wfn(wavy);
		Hirsh[i].set_zero();
		log2 << "Calcualting Hirshfeld density for atom: " << i << endl;
		Calc_Hirshfeld_atom_2(Hirsh[i], Rho, SPHER, wavy, i, log2);
		charges[i] = Hirsh[i].sum();
		log2 << "..done!" << endl;
	}
	for (int i = 0; i < wavy.get_ncen(); i++)
		dipole_moments.push_back(calc_dipole_for_atom(wavy, i, Hirsh[i], charges));
	log2 << "...done!" << endl;
	log2 << " atom   |    charge    | dipole moment x,        y,         z" << endl
		<< "===================================================" << endl;
	for (int i = 0; i < wavy.get_ncen(); i++)
	{
		log2 << setw(3) << i << " (" << constants::atnr2letter(wavy.get_atom_charge(i)) << ") |" << scientific << setprecision(6) << setw(13) << dipole_moments[i][3] - wavy.get_atom_charge(i) << " | " << scientific << setprecision(6) << setw(14) << dipole_moments[i][0] << ", " << setw(14) << dipole_moments[i][1] << ", " << setw(14) << dipole_moments[i][2] << "\n";
	}
	return dipole_moments;
}
