#include "ML_density.h"

cube ML::calc_cube(vec data, WFN& dummy, int& exp_coef)
{
	double MinMax[6]{ 0, 0, 0, 0, 0, 0 };
	int steps[3]{ 0, 0, 0 };
	readxyzMinMax_fromWFN(dummy, MinMax, steps, 2.5, 0.03);
	cube CubeRho(steps[0], steps[1], steps[2], dummy.get_ncen(), true);
	CubeRho.give_parent_wfn(dummy);

	for (int i = 0; i < 3; i++)
	{
		CubeRho.set_origin(i, MinMax[i]);
		CubeRho.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
	}
	CubeRho.set_comment1("Calculated density using NoSpherA2 from ML Data");
	CubeRho.set_comment2("from " + dummy.get_path());

	time_point start = get_time();

	progress_bar* progress = new progress_bar{ std::cout, 50u, "Calculating Values" };
	const int step = (int)std::max(floor(CubeRho.get_size(0) / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < CubeRho.get_size(0); i++)
	{
		for (int j = 0; j < CubeRho.get_size(1); j++)
			for (int k = 0; k < CubeRho.get_size(2); k++)
			{
				vec PosGrid{
					i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
					i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
					i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2) };

				CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.atoms, exp_coef));
			}
		if (i != 0 && i % step == 0)
			progress->write(i / double(CubeRho.get_size(0)));
	}
	delete (progress);

	using namespace std;
	time_point end = get_time();
	if (get_sec(start, end) < 60)
		std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
	else if (get_sec(start, end) < 3600)
		std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
	else
		std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
	std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << CubeRho.sum() << std::endl;
	return CubeRho;
};

cube ML::data_2_Cube(const std::string& xyz_File, const std::string& npy_Coeffs)
{
	using namespace std;

	vector<unsigned long> shape{};
	bool fortran_order;
	vec data{};

	// string path{ "water.npy" };
	npy::LoadArrayFromNumpy(npy_Coeffs, shape, fortran_order, data);

	WFN dummy(7);
	dummy.read_xyz(xyz_File, std::cout);

	int nr_coefs = load_basis_into_WFN(dummy, QZVP_JKfit);//def2_SVP_JKFIT
	cout << data.size() << " vs. " << nr_coefs << " coefficients" << endl;
	cube cub = calc_cube(data, dummy, nr_coefs);
	cub.path = get_basename_without_ending(npy_Coeffs) + ".cube";
	cub.write_file(true);

	cout << "Done :)" << endl;
	return cub;
}

void ML::calc_diff(const std::string xyz_File) {
	std::cout << xyz_File + ".xyz" << std::endl;
	cube cube_ML = data_2_Cube(xyz_File + ".xyz", xyz_File + "_pred.npy");
	cube cube_Ref = data_2_Cube(xyz_File + ".xyz", xyz_File + "_ref.npy");
	WFN wavy(7);
	wavy.read_xyz(xyz_File + ".xyz", std::cout, false);
	double cube_rrs = cube_ML.rrs(cube_Ref);
	std::cout << "real space R-Value: " << cube_rrs << std::endl;
	cube cube_Diff = cube_Ref - cube_ML;
	cube_Diff.give_parent_wfn(wavy);
	cube_Diff.path = "diff_Rho.cube";
	cube_Diff.write_file();
}


void ML::cubeDiffDaniel() {
	std::ofstream log2("NoSpherA2_cube.log", std::ios::out);
	auto coutbuf = std::cout.rdbuf(log2.rdbuf()); // save and redirect
	log2 << NoSpherA2_message();

	WFN wavy(0);
	wavy.read_known_wavefunction_format("HF_def2-SVP.gbw", log2, true);
	//double MinMax[6]{ 0, 0, 0, 0, 0, 0 };
	//int steps[3]{ 0, 0, 0 };
	//std::vector<std::vector<double>> cell_matrix;
	//cell_matrix.resize(3);
	//readxyzMinMax_fromCIF("HF_def2-SVP.cif", MinMax, steps, cell_matrix, 0.1, log2, false);
	//wavy2.read_gbw("test2.gbw", std::cout);
	
}

void ML::gbw2DM(std::string& fn, std::ostream& file, bool& debug) {
	using namespace std;
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
	npy::npy_data<double> data;
	data.data = DM;
	data.fortran_order = false;
	data.shape = { unsigned long(naotr) };
	npy::write_npy("coeffs_by_black_magic.npy", data);
}


