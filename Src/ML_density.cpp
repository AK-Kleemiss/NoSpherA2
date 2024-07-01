#include "ML_density.h"
using namespace std;


void ML_density::calc_cube(
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	double radius,
	int& exp_coef,
	vec data,
	ostream& file)
{
#ifdef _OPENMP
	if (cpus != -1)
	{
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif

	time_point start = get_time();

	progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
	const int step = (int)max(floor(CubeRho.get_size(0) * 3 / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < CubeRho.get_size(0); i++)
	{
		for (int j = 0; j < CubeRho.get_size(1); j++)
			for (int k = 0; k < CubeRho.get_size(2); k++)
			{

				double PosGrid[3]{
					i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
					i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
					i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2) };

				bool skip = true;
				for (int a = 0; a < wavy.get_ncen(); a++)
					if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
						skip = false;
				if (skip)
					continue;

				CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, wavy.atoms, exp_coef));
			}
		if (i != 0 && i % step == 0)
			progress->write((i + CubeRho.get_size(0)) / double(CubeRho.get_size(0) * 3));
	}
	delete (progress);

	time_point end = get_time();
	if (get_sec(start, end) < 60)
		file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
	else if (get_sec(start, end) < 3600)
		file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
	else
		file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
};


cube ML_density::data_2_Cube(const string& xyz_File, const string& npy_Coeffs, int cpus)
{
	using namespace std;
	vector<unsigned long> shape{};
	bool fortran_order;
	vec data{};
	if (cpus != -1)
	{
#ifdef _OPENMP
		omp_set_num_threads(cpus);
#endif
	}

	// string path{ "water.npy" };
	npy::LoadArrayFromNumpy(npy_Coeffs, shape, fortran_order, data);

	WFN wavy(7);
	wavy.read_xyz(xyz_File, cout);

	int nr_coefs = load_basis_into_WFN(wavy, QZVP_JKfit);// QZVP_JKfit);//def2_SVP_JKFIT
	cout << data.size() << " vs. " << nr_coefs << " coefficients" << endl;


	double MinMax[6]{ 0, 0, 0, 0, 0, 0 };
	int steps[3]{ 0, 0, 0 };
	double radius = 2.5;
	readxyzMinMax_fromWFN(wavy, MinMax, steps, radius, 0.03);
	cube CubeRho(steps[0], steps[1], steps[2], wavy.get_ncen(), true);
	CubeRho.give_parent_wfn(wavy);

	for (int i = 0; i < 3; i++)
	{
		CubeRho.set_origin(i, MinMax[i]);
		CubeRho.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
	}
	CubeRho.set_comment1("Calculated density using NoSpherA2 from ML Data");
	CubeRho.set_comment2("from " + wavy.get_path());

	calc_cube(CubeRho, wavy, cpus, radius, nr_coefs, data, cout);
	cout << "Number of electrons: " << fixed << setprecision(4) << CubeRho.sum() << endl;

	CubeRho.path = get_basename_without_ending(npy_Coeffs) + ".cube";
	CubeRho.write_file(true);

	cout << "Done :)" << endl;
	return CubeRho;
}

void ML_density::calc_diff(const string xyz_File, options& opt) {
	cout << xyz_File + ".xyz" << endl;
	cout << "Calculating Reference:" << endl;
	//cube cube_DFT = data_2_Cube(xyz_File + ".xyz", xyz_File + "_DFT.npy", opt.threads);
	//cout << "Calculating Specific ML Density:" << endl;
	//cube cube_Specific = data_2_Cube(xyz_File + ".xyz", xyz_File + "_ML_Specific.npy", opt.threads);
	cout << "Calculating Diverse ML Density::" << endl;
	cube cube_Diverse = data_2_Cube(xyz_File + ".xyz", xyz_File + "_ML_Diverse.npy", opt.threads);
	////double cube_rrs = cube_ML.rrs(cube_Ref);

	//cout << "Calculating Spherical Density:" << endl;
	//WFN wavy(7);
	////wavy.read_known_wavefunction_format(xyz_File + ".xyz", cout, false);
	//wavy.read_xyz(xyz_File + ".xyz", cout);
	//int nr_coefs = load_basis_into_WFN(wavy, QZVP_JKfit);//def2_SVP_JKFIT

	//double MinMax[6]{ 0, 0, 0, 0, 0, 0 };
	//int steps[3]{ 0, 0, 0 };
	//double radius = 2.5;
	//readxyzMinMax_fromWFN(wavy, MinMax, steps, radius, 0.03);
	//cube cubeSpherical(steps[0], steps[1], steps[2], wavy.get_ncen(), true);
	//cubeSpherical.give_parent_wfn(wavy);
	//for (int i = 0; i < 3; i++)
	//{
	//	cubeSpherical.set_origin(i, MinMax[i]);
	//	cubeSpherical.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
	//}
	//Calc_Spherical_Dens(cubeSpherical, wavy, opt.threads, radius, cout);
	//cout << "Number of electrons: " << fixed << setprecision(4) << cubeSpherical.sum() << endl;
	//cubeSpherical.path = "cube_speherical.cube";
	//cubeSpherical.write_file(true);

	//cout << "Calculating Diff sums:" << endl;
	//cube cube_Diff_Spherical = cube_DFT - cubeSpherical;
	//cube cube_Diff_ML_specific = cube_DFT- cube_Specific;
	//cube cube_Diff_ML_Diverse = cube_DFT - cube_Diverse;
	//cout << "Diff DFT-Spherical:" << fixed << setprecision(4)<< cube_Diff_Spherical.diff_sum() << endl;
	//cout << "Diff DFT-Specific:" << fixed << setprecision(4) << cube_Diff_ML_specific.diff_sum() << endl;
	//cout << "Diff DFT-Diverse:" << fixed << setprecision(4) << cube_Diff_ML_Diverse.diff_sum() << endl;

	//cout << "Writing  Diffs:" << endl;
	//WFN wavy2(7);
	//wavy2.read_xyz(xyz_File + ".xyz", cout, false);
	//cube_Diff_Spherical.give_parent_wfn(wavy2);
	//cube_Diff_ML_specific.give_parent_wfn(wavy2);
	//cube_Diff_ML_Diverse.give_parent_wfn(wavy2);

	//cube_Diff_Spherical.path = "diff_Spherical.cube";
	//cube_Diff_ML_specific.path = "diff_Specific.cube";
	//cube_Diff_ML_Diverse.path = "diff_Diverse.cube";

	//cube_Diff_Spherical.write_file(true);
	//cube_Diff_ML_specific.write_file(true);
	//cube_Diff_ML_Diverse.write_file(true);
}


void ML_density::cubeDiffDaniel() {
	//ofstream log2("NoSpherA2_cube.log", ios::out);
	//auto coutbuf = cout.rdbuf(log2.rdbuf()); // save and redirect
	//log2 << NoSpherA2_message();

	//WFN wavy(0);
	//wavy.read_known_wavefunction_format("HF_def2-SVP.gbw", log2, true);
	//double MinMax[6]{ 0, 0, 0, 0, 0, 0 };
	//int steps[3]{ 0, 0, 0 };
	//vector<vector<double>> cell_matrix;
	//cell_matrix.resize(3);
	//readxyzMinMax_fromCIF("HF_def2-SVP.cif", MinMax, steps, cell_matrix, 0.1, log2, false);
	//wavy2.read_gbw("test2.gbw", cout);
	
}


void ML_density::gbw2DM(string& fn, ostream& file, bool& debug) {
	//UNTESTED AND NOT USED AT THE MOMENT
	ifstream rf(fn.c_str(), ios::binary);
	if (!rf.good())
		err_checkf(false, "Nope!", file);
	string line;

	rf.seekg(24, ios::beg);
	long int MOs_start = 0;
	rf.read((char*)&MOs_start, 8);
	err_checkf(rf.good(), "Error reading MO_start", file);
	err_checkf(MOs_start != 0, "Could not read MO information location from GBW file!", file);
	if (debug)
		file << "I read the pointer of MOs succesfully" << endl;
	rf.seekg(MOs_start, ios::beg);
	int operators = 0, dimension = 0;
	rf.read((char*)&operators, 4);
	err_checkf(rf.good(), "Error reading operators", file);
	rf.read((char*)&dimension, 4);
	err_checkf(rf.good(), "Error reading dimnesion", file);
	size_t coef_nr = size_t(dimension) * size_t(dimension);
	vector<vec> coefficients(operators);
	vector<vec> occupations(operators);
	vector<vec> energies(operators);
	vector<ivec> irreps(operators);
	vector<ivec> cores(operators);
	for (int i = 0; i < operators; i++)
	{
		coefficients[i].resize(coef_nr, 0);
		occupations[i].resize(dimension, 0);
		energies[i].resize(dimension, 0);
		irreps[i].resize(dimension, 0);
		cores[i].resize(dimension, 0);
		if (debug)
			file << "operators: " << operators << " coef_nr: " << coef_nr << " dimension: " << dimension << endl;
		rf.read((char*)coefficients[i].data(), 8 * coef_nr);
		err_checkf(rf.good(), "Error reading coefficients", file);
		if (debug)
			file << "I read the coefficients succesfully" << endl;
		rf.read((char*)occupations[i].data(), 8 * dimension);
		err_checkf(rf.good(), "Error reading occupations", file);
		if (debug)
			file << "I read the occupations succesfully" << endl;
		rf.read((char*)energies[i].data(), 8 * dimension);
		err_checkf(rf.good(), "Error reading energies", file);
		if (debug)
			file << "I read the energies succesfully" << endl;
		rf.read((char*)irreps[i].data(), 4 * dimension);
		err_checkf(rf.good(), "Error reading irreps", file);
		if (debug)
			file << "I read the irreps succesfully" << endl;
		rf.read((char*)cores[i].data(), 4 * dimension);
		err_checkf(rf.good(), "Error reading cores", file);
	}

	int naotr = dimension * (dimension + 1) / 2;
	vector<double> DM(naotr);

	for (int iu = 0; iu < dimension; iu++) {
#pragma omp parallel for
		for (int iv = 0; iv <= iu; iv++) {
			const int iuv = (iu * (iu + 1) / 2) + iv;
			double temp;
			for (int op = 0; op < operators; op++) {
				for (int m = 0; m < dimension; m++) {
					if (occupations[op][m] == 0.0) continue;
					temp = occupations[op][m] * coefficients[op][iu + (m * dimension)] * coefficients[op][iv + (m * dimension)];
					DM[iuv] += temp;
				}
			}
		}
	}
	//npy::write_npy("coeffs_by_black_magic.npy", DM);
}


