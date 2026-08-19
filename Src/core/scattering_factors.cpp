/**
 * @file scattering_factors.cpp
 * @brief This file contains functions related to reading and generating hkl indices and k-points.
 *
 */
#include "pch.h"
#include "tsc_block.h"
#include "tsc_stream.h"
#include "scattering_factors.h"
#include "convenience.h"
#include "cell.h"
#include "wfn_class.h"
#include "spherical_density.h"
#include "AtomGrid.h"
#include "npy.h"
#include "integrator.h"
#include "basis_set.h"
#include "SALTED_utilities.h"
#include "GridManager.h"
#include "cube.h"


#ifdef PEOJECT_NAME
#define FLO_CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "CUDA_utilities.h"
#endif
#ifdef __APPLE__
#include "TargetConditionals.h"
#endif

 /**
  * @brief Reads k-points from a binary file and stores them in provided vectors.
  *
  * This function reads k-points from a file named "kpts.dat". The file is expected to be in binary format.
  * The function first checks if the file exists and can be read. It then reads the number of k-points and
  * stores the k-points and their corresponding hkl values in the provided vectors.
  *
  * @param k_pt A vector of vectors to store the k-points.
  * @param hkl A list to store the hkl values corresponding to the k-points.
  * @param file An output stream to write log messages.
  *
  * @throws runtime_error If the file does not exist or cannot be read.
  */
void read_k_points(vec2& k_pt, hkl_list& hkl, std::ostream& file)
{
	err_checkf(std::filesystem::exists("kpts.dat"), "k-points file does not exist!", file);
	file << "Reading:                                    kpts.dat" << std::flush;
	std::ifstream k_points_file("kpts.dat", std::ios::binary);
	err_checkf(k_points_file.good(), "Error Reading the k-points!", file);
	int nr[1]{ 0 };
	k_points_file.read((char*)&nr, sizeof(nr));
	file << " expecting " << nr[0] << " k points... " << std::flush;
	double temp[1]{ 0.0 };
	int i3emp[1]{ 0 };
	k_pt.resize(3);
	for (int i = 0; i < 3; i++)
	{
		k_pt[i].reserve(nr[0]);
	}
	i3 hkl_;
	for (int run = 0; run < nr[0]; run++)
	{
		for (int i = 0; i < 3; i++)
		{
			k_points_file.read((char*)&temp, sizeof(temp));
			err_checkf(!k_points_file.bad(), "Error reading k-points file!", file);
			k_pt[i].emplace_back(temp[0]);
			k_points_file.read((char*)&i3emp, sizeof(i3emp));
			err_checkf(!k_points_file.bad(), "Error reading hkl values from k-points file!", file);
			hkl_[i] = i3emp[0];
		}
		hkl.emplace(hkl_);
	}
	err_checkf(!k_points_file.bad(), "Error reading k-points file!", file);
	file << " done!" << std::endl
		<< "Size of k_points: " << k_pt[0].size() << std::endl;
	k_points_file.close();
}

/**
 * @brief Saves k-points to a binary file.
 *
 * This function writes k-points and their corresponding hkl values to a file named "kpts.dat". The file is
 * written in binary format. The function first writes the number of k-points, then writes the k-points and
 * their corresponding hkl values.
 *
 * @param k_pt A vector of vectors containing the k-points.
 * @param hkl A list containing the hkl values corresponding to the k-points.
 */
void save_k_points(vec2& k_pt, hkl_list& hkl)
{
	std::ofstream k_points_file("kpts.dat", std::ios::out | std::ios::binary | std::ios::trunc);
	int nr[1] = { (int)k_pt[0].size() };
	k_points_file.write((char*)&nr, sizeof(nr));
	double temp[1]{ 0.0 };
	int i3emp[1]{ 0 };
	hkl_list_it hkl_ = hkl.begin();
	for (int run = 0; run < nr[0]; run++)
	{
		for (int i = 0; i < 3; i++)
		{
			temp[0] = k_pt[i][run];
			k_points_file.write((char*)&temp, sizeof(temp));
			i3emp[0] = (*hkl_)[i];
			k_points_file.write((char*)&i3emp, sizeof(i3emp));
		}
		hkl_ = next(hkl_);
	}
	k_points_file.flush();
	k_points_file.close();
}

/**
 * Generates k-points based on the given parameters.
 *
 * @param read_k_pts Flag indicating whether to read k-points from a file.
 * @param save_k_pts Flag indicating whether to save generated k-points to a file.
 * @param gridsize The size of the grid used for generating k-points.
 * @param unit_cell The unit cell used for generating k-points.
 * @param hkl The list of hkl values used for generating k-points.
 * @param k_pt The vector to store the generated k-points.
 * @param file The output stream to write the generated k-points.
 * @param debug Flag indicating whether to enable debug mode.
 */

void make_k_pts(const bool& read_k_pts,
	const bool& save_k_pts,
	const cell& unit_cell,
	hkl_list& hkl,
	vec2& k_pt,
	std::ostream& file,
	bool debug)
{
	using namespace std;
	const int size = (int)hkl.size();
	if (!read_k_pts)
	{
		file << "Generating k-points..." << flush;
		k_pt.resize(3);
#pragma omp parallel for
		for (int i = 0; i < 3; i++)
			k_pt[i].resize(size, 0.0);

		if (debug)
			file << "K_point_vector is here! size: " << k_pt[0].size() << endl;
		// Create local copy of hkl list for faster access
		const std::vector<i3> hkl_vector(hkl.begin(), hkl.end());

#pragma omp parallel for
		for (int ref = 0; ref < size; ref++)
		{
			const i3 hkl_ = hkl_vector[ref];
			for (int x = 0; x < 3; x++)
			{
				for (int j = 0; j < 3; j++)
				{
					k_pt[x][ref] += unit_cell.get_rcm(x, j) * hkl_[j];
				}
			}
		}
		file << "                            ... done!\nNumber of k-points to evaluate: " << k_pt[0].size() << endl;
		if (save_k_pts)
			save_k_points(k_pt, hkl);
	}
	else
	{
		read_k_points(k_pt, hkl, file);
	}
}

/**
 * Reads the hkl data from the specified file and populates the hkl_list with the data.
 *
 * @param hkl_filename The filename of the hkl file to read.
 * @param hkl The hkl_list to populate with the read data.
 * @param twin_law The vector of twin laws to apply to the hkl data.
 * @param unit_cell The cell object representing the unit cell.
 * @param file The output stream to write debug information to.
 * @param debug Flag indicating whether debug information should be printed.
 */
void read_hkl(const std::filesystem::path& hkl_filename,
	hkl_list& hkl,
	const vec2& twin_law,
	cell& unit_cell,
	std::ostream& file,
	bool debug
)
{
	file << "Reading: " << std::setw(44) << hkl_filename << std::flush;
	i3 hkl_;
	err_checkf(std::filesystem::exists(hkl_filename), "HKL file does not exists!", file);
	std::ifstream hkl_input(hkl_filename, std::ios::in);
	hkl_input.seekg(0, hkl_input.beg);
	std::regex r{ R"([abcdefghijklmnopqrstuvwxyz\(\)ABCDEFGHIJKLMNOPQRSTUVW])" };
	std::string line, temp;
	while (!hkl_input.eof())
	{
		getline(hkl_input, line);
		if (hkl_input.eof())
			break;
		if (line.size() < 2)
			continue;
		std::cmatch result;
		if (regex_search(line.c_str(), result, r))
			continue;
		// if (debug) file << "hkl: ";
		for (int i = 0; i < 3; i++)
		{
			temp = line.substr(4 * size_t(i) + 1, 3);
			temp.erase(remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
			hkl_[i] = stoi(temp);
			// if (debug) file << setw(4) << temp;
		}
		// if (debug) file << endl;
		hkl.emplace(hkl_);
	}
	hkl_list_it found = hkl.find(i3{ 0, 0, 0 });
	if (found != hkl.end())
	{
		if (debug)
			file << "popping back 0 0 0" << std::endl;
		hkl.erase(i3{ 0, 0, 0 });
	}
	hkl_input.close();
	file << " done!\nNr of reflections read from file: " << hkl.size() << std::endl;

	if (debug)
		file << "Number of reflections before twin: " << hkl.size() << std::endl;
	if (twin_law.size() > 0)
	{
		for (const i3& hkl__ : hkl)
			for (int i = 0; i < twin_law.size(); i++)
				hkl.emplace(i3{
					static_cast<int>(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2]),
					static_cast<int>(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2]),
					static_cast<int>(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2]) });
	}
	if (debug)
		file << "Number of reflections after twin: " << hkl.size() << std::endl;

	std::vector<std::vector<ivec>> sym(3);
	for (int i = 0; i < 3; i++)
		sym[i].resize(3);
	sym = unit_cell.get_sym();

	if (debug)
	{
		file << "Read " << sym[0][0].size() << " symmetry elements!" << std::endl;
		for (int i = 0; i < sym[0][0].size(); i++)
		{
			for (int x = 0; x < 3; x++)
			{
				for (int y = 0; y < 3; y++)
					file << std::setw(3) << sym[y][x][i];
				file << std::endl;
			}
			file << std::endl;
		}
	}
	else
		file << "Number of symmetry operations: " << sym[0][0].size() << std::endl;

	i3 tempv;
	hkl_list hkl_enlarged = hkl;
	for (int s = 0; s < sym[0][0].size(); s++)
	{
		if (sym[0][0][s] == 1 && sym[1][1][s] == 1 && sym[2][2][s] == 1 &&
			sym[0][1][s] == 0 && sym[0][2][s] == 0 && sym[1][2][s] == 0 &&
			sym[1][0][s] == 0 && sym[2][0][s] == 0 && sym[2][1][s] == 0)
		{
			continue;
		}
		for (const i3& hkl__ : hkl)
		{
			tempv = { 0, 0, 0 };
			for (int h = 0; h < 3; h++)
			{
				for (int j = 0; j < 3; j++)
					tempv[j] += hkl__[h] * sym[j][h][s];
			}
			hkl_enlarged.emplace(tempv);
		}
	}

	for (const i3& hkl__ : hkl_enlarged)
	{
		tempv = hkl__;
		tempv[0] *= -1;
		tempv[1] *= -1;
		tempv[2] *= -1;
		if (hkl_enlarged.find(tempv) != hkl_enlarged.end())
		{
			hkl_enlarged.erase(tempv);
		}
	}
	if (unit_cell.get_sym().empty()) {
		hkl = hkl_enlarged;
	}
	// Remove 0 0 0 if it exists
	if (hkl.find(i3{ 0, 0, 0 }) != hkl.end())
		hkl.erase(i3{ 0, 0, 0 });
	file << "Nr of reflections to be used: " << hkl.size() << std::endl;
}

hkl_list read_hkl_full(const std::filesystem::path& hkl_filename,
	hkl_list& hkl,
	const vec2& twin_law,
	cell& unit_cell,
	std::ostream& file,
	std::vector<scattering_data>& obs,
	bool debug
)
{
	file << "Reading: " << std::setw(44) << hkl_filename << std::flush;
	i3 hkl_;
	double F_, sigma_;
	int positive_;
	err_checkf(std::filesystem::exists(hkl_filename), "HKL file does not exists!", file);
	std::ifstream hkl_input(hkl_filename, std::ios::in);
	hkl_input.seekg(0, hkl_input.beg);
	std::regex r{ R"([abcdefghijklmnopqrstuvwxyz\(\)ABCDEFGHIJKLMNOPQRSTUVW])" };
	std::string line, temp;
	while (!hkl_input.eof())
	{
		getline(hkl_input, line);
		if (hkl_input.eof())
			break;
		if (line.size() < 2)
			continue;
		std::cmatch result;
		if (regex_search(line.c_str(), result, r))
			continue;
		// if (debug) file << "hkl: ";
		for (int i = 0; i < 3; i++)
		{
			temp = line.substr(4 * size_t(i) + 1, 4);
			temp.erase(remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
			hkl_[i] = stoi(temp);
			// if (debug) file << setw(4) << temp;
		}
		for (int i = 0; i < 2; i++) {
			temp = line;
			temp.erase(0, 12);
			int dot = temp.find_first_of('.');
			std::string temp_F = temp.substr(0, dot + 3);
			std::string temp_sigma = temp.substr(dot + 3, temp.size() - dot - 3);
			temp_sigma.erase(remove_if(temp_sigma.begin(), temp_sigma.end(), ::isspace), temp_sigma.end());
			sigma_ = stof(temp_sigma);
			temp_F.erase(remove_if(temp_F.begin(), temp_F.end(), ::isspace), temp_F.end());
			F_ = stof(temp_F);
			if (F_ < 0) {
				F_ = -std::sqrt(-F_);
				positive_ = -1;
				sigma_ = sigma_ * 0.5 / -F_;
			}
			else {
				positive_ = 1;
				F_ = std::sqrt(F_);
				sigma_ = sigma_ * 0.5 / F_;
			}
		}
		// if (debug) file << endl;
		hkl.emplace(hkl_);
		scattering_data temp_data = { F_, sigma_, positive_ };
		obs.push_back(temp_data);
	}
	hkl_list_it found = hkl.find(i3{ 0, 0, 0 });
	if (found != hkl.end())
	{
		if (debug)
			file << "popping back 0 0 0" << std::endl;
		hkl.erase(i3{ 0, 0, 0 });
	}
	hkl_input.close();
	file << " done!\nNr of reflections read from file: " << hkl.size() << std::endl;

	if (debug)
		file << "Number of reflections before twin: " << hkl.size() << std::endl;
	if (twin_law.size() > 0)
	{
		for (const i3& hkl__ : hkl)
			for (int i = 0; i < twin_law.size(); i++)
				hkl.emplace(i3{
					static_cast<int>(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2]),
					static_cast<int>(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2]),
					static_cast<int>(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2]) });
	}
	if (debug)
		file << "Number of reflections after twin: " << hkl.size() << std::endl;

	std::vector<std::vector<ivec>> sym(3);
	for (int i = 0; i < 3; i++)
		sym[i].resize(3);
	sym = unit_cell.get_sym();

	if (debug)
	{
		file << "Read " << sym[0][0].size() << " symmetry elements!" << std::endl;
		for (int i = 0; i < sym[0][0].size(); i++)
		{
			for (int x = 0; x < 3; x++)
			{
				for (int y = 0; y < 3; y++)
					file << std::setw(3) << sym[y][x][i];
				file << std::endl;
			}
			file << std::endl;
		}
	}
	else
		file << "Number of symmetry operations: " << sym[0][0].size() << std::endl;

	i3 tempv;
	hkl_list hkl_enlarged;
	hkl_enlarged = hkl;
	for (int s = 0; s < sym[0][0].size(); s++)
	{
		if (sym[0][0][s] == 1 && sym[1][1][s] == 1 && sym[2][2][s] == 1 &&
			sym[0][1][s] == 0 && sym[0][2][s] == 0 && sym[1][2][s] == 0 &&
			sym[1][0][s] == 0 && sym[2][0][s] == 0 && sym[2][1][s] == 0)
		{
			continue;
		}
		for (const i3& hkl__ : hkl)
		{
			tempv = { 0, 0, 0 };
			for (int h = 0; h < 3; h++)
			{
				for (int j = 0; j < 3; j++)
					tempv[j] += hkl__[h] * sym[j][h][s];
			}
			hkl_enlarged.emplace(tempv);
		}
	}

	// Remove 0 0 0 if it exists
	if (hkl.find(i3{ 0, 0, 0 }) != hkl.end())
		hkl.erase(i3{ 0, 0, 0 });
	file << "Nr of reflections to be used: " << hkl.size() << std::endl;
	return hkl_enlarged;
}

/**
 * Generates a list of hkl values based on the given parameters.
 *
 * @param dmin The minimum value for d-spacing.
 * @param hkl The output list of hkl values.
 * @param twin_law The list of twin laws.
 * @param unit_cell The unit cell object.
 * @param file The output stream to write the generated hkl values.
 * @param debug A flag indicating whether to enable debug mode.
 */
void generate_hkl(const double& dmin,
	hkl_list& hkl,
	const vec2& twin_law,
	cell& unit_cell,
	std::ostream& file,
	bool debug)
{
	using namespace std;
	file << "Generating hkl indices up to d=: " << fixed << setw(17) << setprecision(2) << dmin << flush;
	i3 hkl_;
	string line, temp;
	/* A reflection with spacing d has |h| <= a/d, and likewise for k and l with
	b and c, in any cell: the bound is the real-space axis length times the
	reciprocal radius. It is tight, so the loops have to reach it. They used to
	stop one short and lean on a 0.01 A widening of dmin to make the difference
	up, which is not the same thing - the slack that buys scales with the axis
	length, so it covers the shortfall on a long axis and not on a short one.
	For a dynamical calculation the list has to be complete rather than nearly
	so: a single index the table misses is a failed refinement, not a slightly
	worse number.
	*/
	const ivec extreme = {
		int(std::floor(unit_cell.get_a() / dmin)),
		int(std::floor(unit_cell.get_b() / dmin)),
		int(std::floor(unit_cell.get_c() / dmin)) };
	if (debug)
		file << "extreme: " << extreme[0] << " " << extreme[1] << " " << extreme[2] << endl;
	for (int h = -extreme[0]; h <= extreme[0]; h++)
	{
		for (int k = -extreme[1]; k <= extreme[1]; k++)
		{
			for (int l = -extreme[2]; l <= extreme[2]; l++)
			{
				hkl_ = { h, k, l };
				hkl.emplace(hkl_);
			}
		}
	}
	file << "... done!\nNr of reflections generated: " << setw(21) << hkl.size() << endl;

	if (debug)
		file << "Number of reflections before twin: " << hkl.size() << endl;
	if (twin_law.size() > 0)
	{
		for (const i3& hkl__ : hkl)
			for (int i = 0; i < twin_law.size(); i++)
				hkl.emplace(i3{
					int(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2]),
					int(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2]),
					int(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2]) });
	}
	if (debug)
		file << "Number of reflections after twin: " << hkl.size() << endl;

	vector<vector<ivec>> sym(3);
	for (int i = 0; i < 3; i++)
		sym[i].resize(3);
	sym = unit_cell.get_sym();

	if (debug)
	{
		file << "Read " << sym[0][0].size() << " symmetry elements!" << endl;
		for (int i = 0; i < sym[0][0].size(); i++)
		{
			for (int x = 0; x < 3; x++)
			{
				for (int y = 0; y < 3; y++)
					file << setw(3) << sym[y][x][i];
				file << endl;
			}
			file << endl;
		}
	}
	else
		file << "Number of symmetry operations: " << setw(19) << sym[0][0].size() << endl;

	i3 tempv;
	hkl_list hkl_enlarged = hkl;
	for (int s = 0; s < sym[0][0].size(); s++)
	{
		if (sym[0][0][s] == 1 && sym[1][1][s] == 1 && sym[2][2][s] == 1 &&
			sym[0][1][s] == 0 && sym[0][2][s] == 0 && sym[1][2][s] == 0 &&
			sym[1][0][s] == 0 && sym[2][0][s] == 0 && sym[2][1][s] == 0)
		{
			continue;
		}
		for (const i3& hkl__ : hkl)
		{
			tempv = { 0, 0, 0 };
			for (int h = 0; h < 3; h++)
			{
				for (int j = 0; j < 3; j++)
					tempv[j] += hkl__[h] * sym[j][h][s];
			}
			hkl_enlarged.emplace(tempv);
		}
	}
	hkl.clear();
	if (debug)
		file << "Number of reflections after sym gen: " << hkl_enlarged.size() << endl;

	for (const i3& hkl__ : hkl_enlarged)
	{
		if (hkl.find(hkl__) != hkl.end())
			continue;
		tempv = hkl__;
		tempv[0] *= -1;
		tempv[1] *= -1;
		tempv[2] *= -1;
		if (hkl.find(tempv) == hkl.end())
		{
			hkl.emplace(hkl__);
		}
	}
	// Remove 0 0 0 if it exists
	if (hkl.find(i3{ 0, 0, 0 }) != hkl.end())
		hkl.erase(i3{ 0, 0, 0 });
	file << "Nr of reflections to be used: " << setw(20) << hkl.size() << endl;
}

void generate_hkl(const ivec2& hkl_min_max,
	hkl_list& hkl,
	const vec2& twin_law,
	cell& unit_cell,
	std::ostream& file,
	bool debug,
	bool ED)
{
	using namespace std;
	i3 hkl_;
	string line, temp;
	int h_max = std::max(abs(hkl_min_max[0][1]), abs(hkl_min_max[0][0])),
		k_max = std::max(abs(hkl_min_max[1][1]), abs(hkl_min_max[1][0])),
		l_max = std::max(abs(hkl_min_max[2][1]), abs(hkl_min_max[2][0]));
	/* Doubling each axis of the measured box is not the same set as the
	reflections out to half the measured spacing, and it is the latter a
	dynamical calculation asks for. A box and a resolution shell are different
	shapes: the shell reaches further along the shorter axes, so this leaves a
	gap there. Olex2 sends -dmin for ED, which takes precedence over the box
	and generates the right set; this stays for callers that pass only a box,
	and says so rather than looking complete.
	*/
	if (ED) {
		h_max *= 2, k_max *= 2, l_max *= 2;
		file << "Warning: an ED table generated from an index box alone may be "
			"short along the shorter axes; pass -dmin to generate to a resolution "
			"instead.\n";
	}
	file << "Generating hkl between ["
		<< setw(2) << -h_max << "," << setw(2) << h_max << "] ; ["
		<< setw(2) << -k_max << "," << setw(2) << k_max << "] ; ["
		<< setw(2) << -l_max << "," << setw(2) << l_max << "] " << flush;
	for (int h = -h_max; h <= h_max; h++)
	{
		for (int k = -k_max; k <= k_max; k++)
		{
			// only need 0 to extreme, since we have no DISP signal
			for (int l = 0; l <= l_max; l++)
			{
				hkl_ = { h, k, l };
				hkl.emplace(hkl_);
			}
		}
	}
	file << "... done!\nNr of reflections generated: " << setw(21) << hkl.size() << endl;

	if (debug)
		file << "Number of reflections before twin: " << hkl.size() << endl;
	if (twin_law.size() > 0)
	{
		for (const i3& hkl__ : hkl)
			for (int i = 0; i < twin_law.size(); i++)
				hkl.emplace(i3{
					int(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2]),
					int(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2]),
					int(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2]) });
	}
	if (debug)
		file << "Number of reflections after twin: " << hkl.size() << endl;

	vector<vector<ivec>> sym(3);
	for (int i = 0; i < 3; i++)
		sym[i].resize(3);
	sym = unit_cell.get_sym();

	if (debug)
	{
		file << "Read " << sym[0][0].size() << " symmetry elements!" << endl;
		for (int i = 0; i < sym[0][0].size(); i++)
		{
			for (int x = 0; x < 3; x++)
			{
				for (int y = 0; y < 3; y++)
					file << setw(3) << sym[y][x][i];
				file << endl;
			}
			file << endl;
		}
	}
	else
		file << "Number of symmetry operations: " << setw(19) << sym[0][0].size() << endl;

	i3 tempv;
	hkl_list hkl_enlarged = hkl;
	for (int s = 0; s < sym[0][0].size(); s++)
	{
		if (sym[0][0][s] == 1 && sym[1][1][s] == 1 && sym[2][2][s] == 1 &&
			sym[0][1][s] == 0 && sym[0][2][s] == 0 && sym[1][2][s] == 0 &&
			sym[1][0][s] == 0 && sym[2][0][s] == 0 && sym[2][1][s] == 0)
		{
			continue;
		}
		for (const i3& hkl__ : hkl)
		{
			tempv = { 0, 0, 0 };
			for (int h = 0; h < 3; h++)
			{
				for (int j = 0; j < 3; j++)
					tempv[j] += hkl__[h] * sym[j][h][s];
			}
			hkl_enlarged.emplace(tempv);
		}
	}
	hkl.clear();
	if (debug)
		file << "Number of reflections after sym gen: " << hkl_enlarged.size() << endl;

	for (const i3& hkl__ : hkl_enlarged)
	{
		if (hkl.find(hkl__) != hkl.end())
			continue;
		tempv = hkl__;
		tempv[0] *= -1;
		tempv[1] *= -1;
		tempv[2] *= -1;
		if (hkl.find(tempv) == hkl.end())
		{
			hkl.emplace(hkl__);
		}
	}
	// Remove 0 0 0 if it exists
	if (hkl.find(i3{ 0, 0, 0 }) != hkl.end())
		hkl.erase(i3{ 0, 0, 0 });
	file << "Nr of reflections to be used: " << setw(20) << hkl.size() << endl;
}

/**
 * Generates fractional hkl values based on the given parameters.
 *
 * @param dmin The minimum value for d-spacing.
 * @param hkl The list to store the generated hkl values.
 * @param twin_law The vector of twin laws.
 * @param unit_cell The cell object representing the unit cell.
 * @param file The output stream to write the generated values.
 * @param stepsize The step size for generating hkl values.
 * @param debug Flag indicating whether to enable debug mode.
 */
void generate_fractional_hkl(const double& dmin,
	hkl_list_d& hkl,
	const vec2& twin_law,
	cell& unit_cell,
	std::ostream& file,
	const d3& stepsize,
	bool debug)
{
	using namespace std;
	file << "Generating hkl indices up to d=: " << fixed << setw(17) << setprecision(2) << dmin << flush;
	d3 hkl_;
	string line, temp;
	const int extreme = 201;
	double dmin_l = 0.9 * dmin;
	const int lim = static_cast<int>(extreme / stepsize[0]);
#pragma omp parallel for private(hkl_)
	for (int h = -lim; h < lim; h++)
	{
		const double _h = h * stepsize[0];
		for (double k = -extreme; k < extreme; k += stepsize[1])
		{
			// only need 0 to extreme, since we have no DISP signal
			for (int l = 0; l < lim; l++)
			{
				hkl_ = { _h, k, l * stepsize[2] };
				if (unit_cell.get_d_of_hkl(hkl_) >= dmin_l)
				{
#pragma omp critical
					hkl.emplace(hkl_);
				}
				else
					break;
			}
		}
	}
	file << "... done!\nNr of reflections generated: " << setw(21) << hkl.size() << endl;

	if (debug)
		file << "Number of reflections before twin: " << hkl.size() << endl;
	if (twin_law.size() > 0)
	{
		for (const d3& hkl__ : hkl)
			for (int i = 0; i < twin_law.size(); i++)
				hkl.emplace(d3{
					(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2]),
					(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2]),
					(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2]) });
	}
	if (debug)
		file << "Number of reflections after twin: " << hkl.size() << endl;

	vector<vector<ivec>> sym(3);
	for (int i = 0; i < 3; i++)
		sym[i].resize(3);
	sym = unit_cell.get_sym();

	if (debug)
	{
		file << "Read " << sym[0][0].size() << " symmetry elements!" << endl;
		for (int i = 0; i < sym[0][0].size(); i++)
		{
			for (int x = 0; x < 3; x++)
			{
				for (int y = 0; y < 3; y++)
					file << setw(3) << sym[y][x][i];
				file << endl;
			}
			file << endl;
		}
	}
	else
		file << "Number of symmetry operations: " << setw(19) << sym[0][0].size() << endl;

	d3 tempv;
	hkl_list_d hkl_enlarged = hkl;
	for (int s = 0; s < sym[0][0].size(); s++)
	{
		if (sym[0][0][s] == 1 && sym[1][1][s] == 1 && sym[2][2][s] == 1 &&
			sym[0][1][s] == 0 && sym[0][2][s] == 0 && sym[1][2][s] == 0 &&
			sym[1][0][s] == 0 && sym[2][0][s] == 0 && sym[2][1][s] == 0)
		{
			continue;
		}
		for (const d3& hkl__ : hkl)
		{
			tempv = { 0, 0, 0 };
			for (int h = 0; h < 3; h++)
			{
				for (int j = 0; j < 3; j++)
					tempv[j] += hkl__[h] * sym[j][h][s];
			}
			hkl_enlarged.emplace(tempv);
		}
	}
	hkl.clear();
	if (debug)
		file << "Number of reflections after sym gen: " << hkl_enlarged.size() << endl;

	for (const d3& hkl__ : hkl_enlarged)
	{
		if (hkl.find(hkl__) != hkl.end())
			continue;
		tempv = hkl__;
		tempv[0] *= -1;
		tempv[1] *= -1;
		tempv[2] *= -1;
		if (hkl.find(tempv) == hkl.end())
		{
			hkl.emplace(hkl__);
		}
	}
	// Remove 0 0 0 if it exists
	if (hkl.find(d3{ 0, 0, 0 }) != hkl.end())
		hkl.erase(d3{ 0, 0, 0 });
	file << "Nr of reflections to be used: " << setw(20) << hkl.size() << endl;
}

/**
 * Reads atoms from a CIF file and performs necessary operations.
 *
 * @param cif_input The input stream for the CIF file.
 * @param input_groups The vector of input groups.
 * @param unit_cell The cell object representing the unit cell.
 * @param wave The WFN object.
 * @param known_atoms The vector of known atoms.
 * @param atom_type_list The vector of atom type list.
 * @param asym_atom_to_type_list The vector of asymmetric atom to type list.
 * @param asym_atom_list The vector of asymmetric atom list.
 * @param needs_grid The vector indicating whether each atom needs a grid.
 * @param file The output stream for the file.
 * @param debug A boolean indicating whether to enable debug mode.
 */
svec read_atoms_from_CIF(std::ifstream& cif_input,
	const ivec& input_groups,
	const cell& unit_cell,
	WFN& wave,
	const svec& known_atoms,
	ivec& atom_type_list,
	ivec& asym_atom_to_type_list,
	ivec& asym_atom_list,
	bvec& needs_grid,
	std::ostream& file,
	const bool debug)
{
	using namespace std;
	if (debug)
		file << "start working on cif" << endl;
	bool atoms_read = false;
	int count_fields = 0;
	int group_field = -1;
	int type_field = -1;
	int position_field[3] = { -1, -1, -1 };
	int label_field = 1000;
	string line;
	cif_input.clear();
	cif_input.seekg(0, cif_input.beg);
	svec labels(wave.get_ncen(), "");
	if (debug && input_groups.size() > 0)
		file << "Group size: " << input_groups.size() << endl;
	else if (debug)
		file << "Starting search loop" << endl;
	while (!cif_input.eof() && !atoms_read)
	{
		getline(cif_input, line);
		// if (debug)
		//     file << "line: " << line << endl;
		if (line.find("loop_") != string::npos)
		{
			count_fields = 0;
			getline(cif_input, line);
			if (debug)
				file << "line in loop field definition: " << trim(line) << endl;
			while (trim(line).find("_") == 0)
			{
				if (line.find("label") != string::npos)
					label_field = count_fields;
				else if (line.find("type_symbol") != string::npos)
					type_field = count_fields;
				else if (line.find("disorder_group") != string::npos)
					group_field = count_fields;
				else if (line.find("fract_x") != string::npos)
					position_field[0] = count_fields;
				else if (line.find("fract_y") != string::npos)
					position_field[1] = count_fields;
				else if (line.find("fract_z") != string::npos)
					position_field[2] = count_fields;
				else if (label_field == 1000)
				{
					if (debug)
						file << "I don't think this is the atom block.. moving on!" << endl;
					break;
				}
				getline(cif_input, line);
				count_fields++;
			}
			if (label_field != 1000) {
				err_checkf(position_field[0] != -1, "No x position found, impossible to continue!", std::cout);
				err_checkf(position_field[1] != -1, "No y position found, impossible to continue!", std::cout);
				err_checkf(position_field[2] != -1, "No z position found, impossible to continue!", std::cout);
				err_checkf(type_field != -1, "No type found, impossible to continue!", std::cout);
			}
			while (trim(line).find("_") > 0 && line.length() > 3)
			{
				atoms_read = true;
				stringstream s(line);
				svec fields;
				fields.resize(count_fields);
				int nr = -1;
				for (int i = 0; i < count_fields; i++)
					s >> fields[i];
				fields[label_field].erase(remove_if(fields[label_field].begin(), fields[label_field].end(), ::isspace), fields[label_field].end());
				fields[type_field].erase(remove_if(fields[type_field].begin(), fields[type_field].end(), ::isspace), fields[type_field].end());
				if (debug)
					file << "label: " << setw(8) << fields[label_field] << " type: " << fields[type_field] << " frac. pos: "
					<< setw(6) << fixed << setprecision(3) << stod(fields[position_field[0]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[0]]) << " "
					<< setw(6) << fixed << setprecision(3) << stod(fields[position_field[1]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[1]]) << " "
					<< setw(6) << fixed << setprecision(3) << stod(fields[position_field[2]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[2]]) << " " << flush;
				vec position = unit_cell.get_coords_cartesian(
					stod(fields[position_field[0]]),
					stod(fields[position_field[1]]),
					stod(fields[position_field[2]]));
				vec precisions = unit_cell.get_coords_cartesian(
					min(0.01, get_decimal_precision_from_CIF_number(fields[position_field[0]])),
					min(0.01, get_decimal_precision_from_CIF_number(fields[position_field[1]])),
					min(0.01, get_decimal_precision_from_CIF_number(fields[position_field[2]])));
				for (int i = 0; i < 3; i++)
				{
					precisions[i] = abs(precisions[i]);
				}
				if (debug)
					file << " cart. pos.: " << setw(8) << position[0] << "+/-" << precisions[0] << " " << setw(8) << position[1] << "+/-" << precisions[1] << " " << setw(8) << position[2] << "+/-" << precisions[2] << endl;

				int group_nr = 0;
				if (group_field != -1 && fields[group_field] != "." && fields[group_field] != "?") {
					group_nr = std::stoi(fields[group_field]);
				}
				// Filter PART before matching this CIF row to a WFN atom.  Disorder
				// positions can be close enough to match an atom in another PART;
				// doing this after the match overwrites that atom's coordinates/ID and
				// corrupts the duplicate check of a later -mtc fragment.
				if (!input_groups.empty() && group_field != -1 &&
					std::find(input_groups.begin(), input_groups.end(), group_nr) == input_groups.end())
				{
					if (debug)
						file << "Wrong part!" << endl;
					getline(cif_input, line);
					continue;
				}
				bool old_atom = false;
				const atomID cif_atom_id(
					stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]),
					group_nr,
					constants::get_Z_from_label(fields[type_field].c_str()) + 1);
				const std::string atom_ID = cif_atom_id.to_hex_string();
#pragma omp parallel for reduction(|| : old_atom)
				for (int run = 0; run < known_atoms.size(); run++)
				{
					if (fields[label_field] == known_atoms[run] || atom_ID == known_atoms[run])
					{
						old_atom = true;
						if (debug)
							file << "I already know this one! " << fields[label_field] << " " << known_atoms[run] << endl;
					}
				}
				if (old_atom)
				{
					getline(cif_input, line);
					continue;
				}
				vec tolerances(3);
				for (int i = 0; i < wave.get_ncen(); i++)
				{
					for (int j = 0; j < 3; j++)
					{
						/* The cap has to be well inside a bond length. It is derived
						from the fractional precision, so it grows with the cell: on a
						55-68 A protein cell it reached 1.3 A per axis and matched
						HB_A:53 to a carbon 0.83 A away, after which the element test
						dropped the atom as a mismatch. The cif and the xyz describe
						the same model, so they agree far better than this.
						*/
						tolerances[j] = 2 * max(min(abs(precisions[j]), 0.1), 0.01);
					}
					if (is_similar_abs(position[0], wave.get_atom_coordinate(i, 0), tolerances[0]) && is_similar_abs(position[1], wave.get_atom_coordinate(i, 1), tolerances[1]) && is_similar_abs(position[2], wave.get_atom_coordinate(i, 2), tolerances[2]))
					{
						wave.set_atom_frac_coords(i, { stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]) });
                        wave.set_atom_group_nr(i, group_nr);
						// Store exactly the identifier used above for the MTC duplicate check.
						wave.set_id_for_atom(i, cif_atom_id);
						string element = constants::atnr2letter(wave.get_atom_charge(i));
						err_checkf(element != "PROBLEM", "Problem identifying atoms!", std::cout);
						string label = fields[label_field];
						string type = fields[type_field];
						transform(element.begin(), element.end(), element.begin(), asciitolower);
						transform(label.begin(), label.end(), label.begin(), asciitolower);
						transform(type.begin(), type.end(), type.begin(), asciitolower);
						if (debug)
						{
							file << "ASYM:  " << setw(8) << element << " charge: " << setw(17) << wave.get_atom_charge(i) << "                             wfn cart. pos: "
								<< fixed << setprecision(3) << setw(16) << wave.get_atom_coordinate(i, 0) << " "
								<< fixed << setprecision(3) << setw(16) << wave.get_atom_coordinate(i, 1) << " "
								<< fixed << setprecision(3) << setw(16) << wave.get_atom_coordinate(i, 2) << flush;
							if (input_groups.size() > 0 && group_field != -1)
							{
								file << " checking disorder group: " << fields[group_field] << " vs. ";
								for (int g = 0; g < input_groups.size(); g++)
									file << input_groups[g] << ",";
							}
						}
						if (label.find(element) == string::npos || label.find(element) > 2)
						{
							if (element != "h")
							{
								if (debug)
								{
									file << "\nElement symbol not found in label, this is a problem!\n checking type...";
									if (type.find(element) == string::npos || label.find(element) > 2)
									{
										file << " ALSO FAILED! WILL IGNORE ATOM!\n";
										continue;
									}
								}
								else
								{
									if (type.find(element) == string::npos || label.find(element) > 2)
									{
										file << "\nAtom " << label << " was not matching by element determined by label reduction or type field, skipping!\n";
										continue;
									}
								}
							}
							else if (label.find("d") == string::npos && label.find("t") == string::npos)
							{
								if (debug)
								{
									file << "\nElement symbol not found in label, this is a problem!\n will check type...";
									if (type.find(element) == string::npos || label.find(element) > 2)
									{
										file << " ALSO FAILED! WILL IGNORE ATOM!\n";
										continue;
									}
								}
								else
								{
									if (type.find(element) == string::npos || label.find(element) > 2)
									{
										file << "\nAtom " << label << " was not matching by element determined by label reduction or type field, skipping!\n";
										continue;
									}
								}
							}
						}

						labels[i] = fields[label_field];
						asym_atom_list.push_back(i);
						needs_grid[i] = true;
						nr = i;
						break;
					}
				}
				if (debug)
					file << " nr= " << nr << endl;
				if (nr != -1)
				{
					bool already_there = false;
					for (int i = 0; i < atom_type_list.size(); i++)
						if (atom_type_list[i] == wave.get_atom_charge(nr))
						{
							already_there = true;
							asym_atom_to_type_list.push_back(i);
							break;
						}
					if (already_there == false && wave.get_atom_charge(nr) != 119)
					{
						asym_atom_to_type_list.push_back((int)atom_type_list.size());
						atom_type_list.push_back(wave.get_atom_charge(nr));
					}
				}
				else if (!old_atom)
				{
					if (debug)
					{
						file << "I did not find this atom! Tolerances were: ";
						for (int j = 0; j < 3; j++)
						{
							file << setw(12) << fixed << setprecision(8) << tolerances[j];
						}
						file << endl;
					}
				}
				getline(cif_input, line);
			}
		}
	}

	// Add missing atom types to be able to calc sphericals correctly
	for (int nr = 0; nr < wave.get_ncen(); nr++)
	{
		bool already_there = false;
		for (int i = 0; i < atom_type_list.size(); i++)
		{
			if (atom_type_list[i] == wave.get_atom_charge(nr))
			{
				already_there = true;
				break;
			}
		}
		if (already_there == false && wave.get_atom_charge(nr) != 119)
		{
			atom_type_list.push_back(wave.get_atom_charge(nr));
		}
	}

	err_checkf(asym_atom_list.size() <= wave.get_ncen(), "More asymmetric unit atoms detected than in the wavefunction! Aborting!", file);
	err_checkf(asym_atom_list.size() != 0, "0 asym atoms is imposible! something is wrong with reading the CIF!", file);

	for (int i = 0; i < atom_type_list.size(); i++)
		err_checkf((atom_type_list[i] <= 113 || atom_type_list[i] == 119) && atom_type_list[i] > 0, "Unreasonable atom type detected: " + toString(atom_type_list[i]) + " (Happens if Atoms were not identified correctly)", file);
	file << " done!" << endl;
	if (debug)
	{
		file << "There are " << atom_type_list.size() << " types of atoms" << endl;
		for (int i = 0; i < atom_type_list.size(); i++)
			file << setw(4) << atom_type_list[i];
		file << endl
			<< "asym_atoms_to_type_list: " << endl;
		for (int i = 0; i < asym_atom_to_type_list.size(); i++)
			file << setw(4) << asym_atom_to_type_list[i];
		file << endl;
		file << "Charges of atoms:" << endl;
		for (int i = 0; i < wave.get_ncen(); i++)
			file << setw(4) << wave.get_atom_charge(i);
		file << endl;
	}
	int size = static_cast<int>(asym_atom_list.size());
	svec labels2;
	for (int i = 0; i < size; i++)
		labels2.emplace_back(labels[asym_atom_list[i]]);
	return labels2;
}

/**
 * Reads atoms from a CIF file and performs necessary operations. Works without wavefunction for XCW routine.
 */
void read_atoms_from_CIF(std::ifstream& cif_input, const cell& unit_cell, int& ncen, bvec& needs_grid, std::vector<asym_atom>& asym_atoms, const bool debug){
	if (!cif_input)
	{
		throw std::runtime_error("Could not open CIF file.");
	}
	std::string line;
	bool in_atom_loop = false;
	bool reading_headers = false;
	std::vector<std::string> headers;
	// Column indices
	int idx_label = -1;
	int idx_type = -1;
	int idx_x = -1;
	int idx_y = -1;
	int idx_z = -1;
	while (std::getline(cif_input, line))
	{
		if (line.empty())
			continue;
		// Trim leading spaces
		size_t first = line.find_first_not_of(" \t");
		if (first == std::string::npos)
			continue;
		std::string trimmed = line.substr(first);
		// Detect loop_
		if (trimmed == "loop_")
		{
			in_atom_loop = false;
			reading_headers = true;
			headers.clear();
			idx_label = -1;
			idx_type = -1;
			idx_x = -1;
			idx_y = -1;
			idx_z = -1;
			continue;
		}
		// Read headers after loop_
		if (reading_headers && !trimmed.empty() && trimmed[0] == '_')
		{
			headers.push_back(trimmed);
			continue;
		}
		// End of headers -> determine whether this is the atom loop
		if (reading_headers)
		{
			reading_headers = false;
			for (size_t i = 0; i < headers.size(); ++i)
			{
				if (headers[i] == "_atom_site_label")
					idx_label = static_cast<int>(i);
				else if (headers[i] == "_atom_site_type_symbol")
					idx_type = static_cast<int>(i);
				else if (headers[i] == "_atom_site_fract_x")
					idx_x = static_cast<int>(i);
				else if (headers[i] == "_atom_site_fract_y")
					idx_y = static_cast<int>(i);
				else if (headers[i] == "_atom_site_fract_z")
					idx_z = static_cast<int>(i);
			}
			// Check if this is the desired atom loop
			if (idx_label >= 0 &&
				idx_type >= 0 &&
				idx_x >= 0 &&
				idx_y >= 0 &&
				idx_z >= 0)
			{
				in_atom_loop = true;
				if (debug)
				{
					std::cout << "Found atom_site loop\n";
					std::cout << "label idx = " << idx_label << "\n";
					std::cout << "type  idx = " << idx_type << "\n";
					std::cout << "x     idx = " << idx_x << "\n";
					std::cout << "y     idx = " << idx_y << "\n";
					std::cout << "z     idx = " << idx_z << "\n";
				}
			}
		}
		// Read atom rows
		if (in_atom_loop)
		{
			// Another loop or new data item ends this loop
			if (trimmed == "loop_" || trimmed[0] == '_')
			{
				in_atom_loop = false;
				// Reprocess line in outer loop logic
				if (trimmed == "loop_")
				{
					reading_headers = true;
					headers.clear();
				}
				continue;
			}
			std::istringstream iss(trimmed);
			std::vector<std::string> tokens;
			std::string token;
			while (iss >> token)
				tokens.push_back(token);
			const int required_cols =
				std::max({ idx_label, idx_type, idx_x, idx_y, idx_z }) + 1;
			if (static_cast<int>(tokens.size()) < required_cols)
			{
				if (debug)
				{
					std::cout << "Skipping malformed atom line:\n";
					std::cout << trimmed << "\n";
				}
				continue;
			}
			asym_atom temp_atom;
			temp_atom.label = tokens[idx_label];
			const std::string& type_str = tokens[idx_type];
			temp_atom.type = constants::get_Z_from_label(type_str.c_str()) + 1;
			const double fx = std::stod(tokens[idx_x]);
			const double fy = std::stod(tokens[idx_y]);
			const double fz = std::stod(tokens[idx_z]);
			temp_atom.frac_pos = { fx, fy, fz };
			auto cart =
				unit_cell.get_coords_cartesian(fx, fy, fz, true);
			temp_atom.pos = { cart[0], cart[1], cart[2] };
			asym_atoms.push_back(temp_atom);
			if (debug)
			{
				std::cout << "Parsed atom:\n";
				std::cout << "  label = " << temp_atom.label << "\n";
				std::cout << "  type  = " << temp_atom.type << "\n";
				std::cout << "  frac  = ("
					<< fx << ", "
					<< fy << ", "
					<< fz << ")\n";
			}
		}
	}
	ncen = static_cast<int>(asym_atoms.size());
	needs_grid.resize(ncen, true);
	if (debug)
	{
		std::cout << "\nTotal atoms parsed: " << ncen << "\n";
	}
}

//svec read_anom_disp_from_CIF(std::ifstream& cif_input,
//    const ivec& input_groups,
//    const cell& unit_cell,
//    const WFN& wave,
//    const svec& known_atoms,
//    ivec& atom_type_list,
//    ivec& asym_atom_to_type_list,
//    ivec& asym_atom_list,
//    bvec& needs_grid,
//    std::ostream& file,
//    bvec& constant_atoms,
//    const bool SALTED,
//    const bool debug)
//{
//    using namespace std;
//    if (debug)
//        file << "start working on cif" << endl;
//    bool atoms_read = false;
//    int count_fields = 0;
//    int group_field = -1;
//    int type_field = -1;
//    int position_field[3] = { -1, -1, -1 };
//    int label_field = 1000;
//    string line;
//    cif_input.clear();
//    cif_input.seekg(0, cif_input.beg);
//    svec labels(wave.get_ncen(), "");
//    if (debug && input_groups.size() > 0)
//        file << "Group size: " << input_groups.size() << endl;
//    else if (debug)
//        file << "Starting search loop" << endl;
//    while (!cif_input.eof() && !atoms_read)
//    {
//        getline(cif_input, line);
//        // if (debug)
//        //     file << "line: " << line << endl;
//        if (line.find("loop_") != string::npos)
//        {
//            count_fields = 0;
//            getline(cif_input, line);
//            if (debug)
//                file << "line in loop field definition: " << trim(line) << endl;
//            while (trim(line).find("_") == 0)
//            {
//                if (line.find("label") != string::npos)
//                    label_field = count_fields;
//                else if (line.find("type_symbol") != string::npos)
//                    type_field = count_fields;
//                else if (line.find("disorder_group") != string::npos)
//                    group_field = count_fields;
//                else if (line.find("fract_x") != string::npos)
//                    position_field[0] = count_fields;
//                else if (line.find("fract_y") != string::npos)
//                    position_field[1] = count_fields;
//                else if (line.find("fract_z") != string::npos)
//                    position_field[2] = count_fields;
//                else if (label_field == 1000)
//                {
//                    if (debug)
//                        file << "I don't think this is the atom block.. moving on!" << endl;
//                    break;
//                }
//                getline(cif_input, line);
//                count_fields++;
//            }
//            if (label_field != 1000) {
//                err_checkf(position_field[0] != -1, "No x position found, impossible to continue!", std::cout);
//                err_checkf(position_field[1] != -1, "No y position found, impossible to continue!", std::cout);
//                err_checkf(position_field[2] != -1, "No z position found, impossible to continue!", std::cout);
//                err_checkf(type_field != -1, "No type found, impossible to continue!", std::cout);
//            }
//            while (trim(line).find("_") > 0 && line.length() > 3)
//            {
//                atoms_read = true;
//                stringstream s(line);
//                svec fields;
//                fields.resize(count_fields);
//                int nr = -1;
//                for (int i = 0; i < count_fields; i++)
//                    s >> fields[i];
//                fields[label_field].erase(remove_if(fields[label_field].begin(), fields[label_field].end(), ::isspace), fields[label_field].end());
//                fields[type_field].erase(remove_if(fields[type_field].begin(), fields[type_field].end(), ::isspace), fields[type_field].end());
//                if (debug)
//                    file << "label: " << setw(8) << fields[label_field] << " type: " << fields[type_field] << " frac. pos: "
//                    << setw(6) << fixed << setprecision(3) << stod(fields[position_field[0]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[0]]) << " "
//                    << setw(6) << fixed << setprecision(3) << stod(fields[position_field[1]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[1]]) << " "
//                    << setw(6) << fixed << setprecision(3) << stod(fields[position_field[2]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[2]]) << " " << flush;
//                vec position = unit_cell.get_coords_cartesian(
//                    stod(fields[position_field[0]]),
//                    stod(fields[position_field[1]]),
//                    stod(fields[position_field[2]]));
//                vec precisions = unit_cell.get_coords_cartesian(
//                    min(0.01, get_decimal_precision_from_CIF_number(fields[position_field[0]])),
//                    min(0.01, get_decimal_precision_from_CIF_number(fields[position_field[1]])),
//                    min(0.01, get_decimal_precision_from_CIF_number(fields[position_field[2]])));
//                for (int i = 0; i < 3; i++)
//                {
//                    precisions[i] = abs(precisions[i]);
//                }
//                if (debug)
//                    file << " cart. pos.: " << setw(8) << position[0] << "+/-" << precisions[0] << " " << setw(8) << position[1] << "+/-" << precisions[1] << " " << setw(8) << position[2] << "+/-" << precisions[2] << endl;
//
//                bool old_atom = false;
//#pragma omp parallel for reduction(|| : old_atom)
//                for (int run = 0; run < known_atoms.size(); run++)
//                {
//                    if (fields[label_field] == known_atoms[run])
//                    {
//                        if (SALTED && (group_field == -1 || fields[group_field].c_str()[0] == '.'))
//                            continue;
//                        old_atom = true;
//                        if (debug)
//                            file << "I already know this one! " << fields[label_field] << " " << known_atoms[run] << endl;
//                    }
//                }
//                if (old_atom)
//                {
//                    getline(cif_input, line);
//                    continue;
//                }
//                vec tolerances(3);
//                for (int i = 0; i < wave.get_ncen(); i++)
//                {
//                    for (int j = 0; j < 3; j++)
//                    {
//                        tolerances[j] = 2 * max(min(abs(precisions[j]), 1.0), 0.01);
//                    }
//                    if (is_similar_abs(position[0], wave.get_atom_coordinate(i, 0), tolerances[0]) && is_similar_abs(position[1], wave.get_atom_coordinate(i, 1), tolerances[1]) && is_similar_abs(position[2], wave.get_atom_coordinate(i, 2), tolerances[2]))
//                    {
//                        string element = constants::atnr2letter(wave.get_atom_charge(i));
//                        err_checkf(element != "PROBLEM", "Problem identifying atoms!", std::cout);
//                        string label = fields[label_field];
//                        string type = fields[type_field];
//                        transform(element.begin(), element.end(), element.begin(), asciitolower);
//                        transform(label.begin(), label.end(), label.begin(), asciitolower);
//                        transform(type.begin(), type.end(), type.begin(), asciitolower);
//                        if (debug)
//                        {
//                            file << "ASYM:  " << setw(8) << element << " charge: " << setw(17) << wave.get_atom_charge(i) << "                             wfn cart. pos: "
//                                << fixed << setprecision(3) << setw(16) << wave.get_atom_coordinate(i, 0) << " "
//                                << fixed << setprecision(3) << setw(16) << wave.get_atom_coordinate(i, 1) << " "
//                                << fixed << setprecision(3) << setw(16) << wave.get_atom_coordinate(i, 2) << flush;
//                            if (input_groups.size() > 0 && group_field != -1)
//                            {
//                                file << " checking disorder group: " << fields[group_field] << " vs. ";
//                                for (int g = 0; g < input_groups.size(); g++)
//                                    file << input_groups[g] << ",";
//                            }
//                        }
//                        if (input_groups.size() > 0)
//                        {
//                            bool yep = false;
//                            for (int g = 0; g < input_groups.size(); g++)
//                            {
//                                if (group_field == -1) {
//                                    yep = true;
//                                    break;
//                                }
//                                if (fields[group_field].c_str()[0] == '.' && input_groups[g] == 0)
//                                {
//                                    if (debug)
//                                        file << "appears to be group 0" << endl;
//                                    yep = true;
//                                    break;
//                                }
//                                else if (stoi(fields[group_field]) == input_groups[g])
//                                    yep = true;
//                            }
//                            if (!yep)
//                            {
//                                if (debug)
//                                    file << "Wrong part!" << endl;
//                                continue;
//                            }
//                        }
//                        if (label.find(element) == string::npos || label.find(element) > 2)
//                        {
//                            if (element != "h")
//                            {
//                                if (debug)
//                                {
//                                    file << "\nElement symbol not found in label, this is a problem!\n checking type...";
//                                    if (type.find(element) == string::npos || label.find(element) > 2)
//                                    {
//                                        file << " ALSO FAILED! WILL IGNORE ATOM!\n";
//                                        continue;
//                                    }
//                                }
//                                else
//                                {
//                                    if (type.find(element) == string::npos || label.find(element) > 2)
//                                    {
//                                        file << "\nAtom " << label << " was not matching by element determined by label reduction or type field, skipping!\n";
//                                        continue;
//                                    }
//                                }
//                            }
//                            else if (label.find("d") == string::npos && label.find("t") == string::npos)
//                            {
//                                if (debug)
//                                {
//                                    file << "\nElement symbol not found in label, this is a problem!\n will check type...";
//                                    if (type.find(element) == string::npos || label.find(element) > 2)
//                                    {
//                                        file << " ALSO FAILED! WILL IGNORE ATOM!\n";
//                                        continue;
//                                    }
//                                }
//                                else
//                                {
//                                    if (type.find(element) == string::npos || label.find(element) > 2)
//                                    {
//                                        file << "\nAtom " << label << " was not matching by element determined by label reduction or type field, skipping!\n";
//                                        continue;
//                                    }
//                                }
//                            }
//                        }
//
//
//                        labels[i] = fields[label_field];
//                        asym_atom_list.push_back(i);
//                        needs_grid[i] = true;
//                        nr = i;
//                        break;
//                    }
//                }
//                if (debug)
//                    file << " nr= " << nr << endl;
//                if (nr != -1)
//                {
//                    bool already_there = false;
//                    for (int i = 0; i < atom_type_list.size(); i++)
//                        if (atom_type_list[i] == wave.get_atom_charge(nr))
//                        {
//                            already_there = true;
//                            asym_atom_to_type_list.push_back(i);
//                            break;
//                        }
//                    if (already_there == false && wave.get_atom_charge(nr) != 119)
//                    {
//                        asym_atom_to_type_list.push_back((int)atom_type_list.size());
//                        atom_type_list.push_back(wave.get_atom_charge(nr));
//                    }
//                }
//                else if (!old_atom)
//                {
//                    if (debug)
//                    {
//                        file << "I did not find this atom! Tolerances were: ";
//                        for (int j = 0; j < 3; j++)
//                        {
//                            file << setw(12) << fixed << setprecision(8) << tolerances[j];
//                        }
//                        file << endl;
//                    }
//                }
//                getline(cif_input, line);
//            }
//        }
//    }
//
//    // Add missing atom types to be able to calc sphericals correctly
//    for (int nr = 0; nr < wave.get_ncen(); nr++)
//    {
//        bool already_there = false;
//        for (int i = 0; i < atom_type_list.size(); i++)
//        {
//            if (atom_type_list[i] == wave.get_atom_charge(nr))
//            {
//                already_there = true;
//                break;
//            }
//        }
//        if (already_there == false && wave.get_atom_charge(nr) != 119)
//        {
//            atom_type_list.push_back(wave.get_atom_charge(nr));
//        }
//    }
//
//    err_checkf(asym_atom_list.size() <= wave.get_ncen(), "More asymmetric unit atoms detected than in the wavefunction! Aborting!", file);
//    err_checkf(asym_atom_list.size() != 0, "0 asym atoms is imposible! something is wrong with reading the CIF!", file);
//
//    for (int i = 0; i < atom_type_list.size(); i++)
//        err_checkf((atom_type_list[i] <= 113 || atom_type_list[i] == 119) && atom_type_list[i] > 0, "Unreasonable atom type detected: " + toString(atom_type_list[i]) + " (Happens if Atoms were not identified correctly)", file);
//    file << " done!" << endl;
//    if (debug)
//    {
//        file << "There are " << atom_type_list.size() << " types of atoms" << endl;
//        for (int i = 0; i < atom_type_list.size(); i++)
//            file << setw(4) << atom_type_list[i];
//        file << endl
//            << "asym_atoms_to_type_list: " << endl;
//        for (int i = 0; i < asym_atom_to_type_list.size(); i++)
//            file << setw(4) << asym_atom_to_type_list[i];
//        file << endl;
//        file << "Charges of atoms:" << endl;
//        for (int i = 0; i < wave.get_ncen(); i++)
//            file << setw(4) << wave.get_atom_charge(i);
//        file << endl;
//    }
//    int size = static_cast<int>(asym_atom_list.size());
//    svec labels2;
//    for (int i = 0; i < size; i++)
//        labels2.emplace_back(labels[asym_atom_list[i]]);
//    return labels2;
//}

constexpr double cutoff(const int& accuracy)
{
	if (accuracy < 3)
		return 1E-10;
	else if (accuracy == 3)
		return 1E-14;
	else
		return 1E-30;
}

// This function yields the fourier bessel transform of the radial integral of a gaussian density function (compare equation 1.2.7.9 in 10.1107/97809553602060000759), assuming that H = 2 \pi S
double fourier_bessel_integral(
	const primitive& p,
	const double& H,
	const int& l)
{
	const double b = p.get_exp();
	return (pow(H, l) * exp(-H * H / (4.0 * b))) / (constants::pow_2[l] * p.get_exp_l_plus_3_2());
}

cdouble sfac_bessel(
	const primitive& p,
	const double* k_point,
	const double* coefs)
{
	const int l = p.get_type();
	switch (l % 4) {
	case 0:
		return cdouble(constants::PI3_2 * fourier_bessel_integral(p, k_point[3], l) * p.get_normalized_coefficient() * constants::spherical_harmonic(l, k_point, coefs), 0);
	case 1:
		return cdouble(0, constants::PI3_2 * fourier_bessel_integral(p, k_point[3], l) * p.get_normalized_coefficient() * constants::spherical_harmonic(l, k_point, coefs));
	case 2:
		return cdouble(-constants::PI3_2 * fourier_bessel_integral(p, k_point[3], l) * p.get_normalized_coefficient() * constants::spherical_harmonic(l, k_point, coefs), 0);
	case 3:
		return cdouble(0, -constants::PI3_2 * fourier_bessel_integral(p, k_point[3], l) * p.get_normalized_coefficient() * constants::spherical_harmonic(l, k_point, coefs));
	default:
		//should normally never happen
		return constants::cnull;
	}
}

//TODO�: This breaks if the aux_basis is contracted... Need to fix that!
void calc_SF_SALTED(const vec2& k_pt,
	const vec& coefs,
	const std::vector<atom>& atom_list,
	const ivec& asym_atom_list,
	cvec2& sf)
{
	const int num_atoms = (int)atom_list.size();
	const int num_asym_atoms = (int)asym_atom_list.size();

	// coefficients per *atom* (full list)
	std::vector<int> atom_ncoefs(num_atoms, 0);

	// global offset for each atom in coefs[]
	std::vector<int> atom_offsets(num_atoms + 1, 0);

#pragma omp parallel for
	for (int iat = 0; iat < num_atoms; ++iat) {
		const atom& a = atom_list[iat];

		int prim = 0;
		int n_this = 0;

		for (int shell = 0; shell < a.get_shellcount_size(); ++shell) {
			int L = a.get_basis_set_type(prim);
			n_this += 2 * L + 1;
			prim += a.get_shellcount(shell);
		}

		atom_ncoefs[iat] = n_this;
	}

	// prefix sum over *all* atoms
	std::partial_sum(atom_ncoefs.begin(),
		atom_ncoefs.end(),
		atom_offsets.begin() + 1);


    ivec coef_offsets(num_asym_atoms, 0);

	for (int ia = 0; ia < num_asym_atoms; ++ia) {
		coef_offsets[ia] = atom_offsets[asym_atom_list[ia]];
	}

	sf.resize(num_asym_atoms);
	ProgressBar pb(k_pt[0].size(), 60, "#", " ", "Generating scattering factors...");

#pragma omp parallel shared(pb, sf)
	{
		// init SF
#pragma omp for
		for (int ia = 0; ia < num_asym_atoms; ++ia) {
			sf[ia].assign(k_pt[0].size(), constants::cnull);
		}

#pragma omp for
		for (int i_kpt = 0; i_kpt < (int)k_pt[0].size(); ++i_kpt)
		{
			double k_pt_local[4] = {
				k_pt[0][i_kpt],
				k_pt[1][i_kpt],
				k_pt[2][i_kpt],
				0.0
			};

			k_pt_local[3] = std::hypot(k_pt_local[0], k_pt_local[1], k_pt_local[2]);

			for (int i = 0; i < 3; i++)
				k_pt_local[i] /= k_pt_local[3];

			for (int ia = 0; ia < num_asym_atoms; ++ia)
			{
				const atom& a = atom_list[asym_atom_list[ia]];

				const basis_set_entry* basis_ptr = &a.get_basis_set_entry(0);
				const int lim = (int)a.get_basis_set_size();

				const double* coef_slice_ptr = coefs.data() + coef_offsets[ia];

				for (int i_basis = 0; i_basis < lim; ++i_basis, ++basis_ptr)
				{
					// IMPORTANT: make basis local, not shared between threads
					const primitive& basis = basis_ptr->get_primitive();
					sf[ia][i_kpt] += sfac_bessel(basis, k_pt_local, coef_slice_ptr);
					coef_slice_ptr += 2 * basis.get_type() + 1;
				}
			}
			pb.update();
		}
	}
}
/**
 * Calculates the scattering factors for a given set of parameters.
 *
 * @param points The number of points.
 * @param k_pt The vector of k points.
 * @param d1 The vector of distances in x.
 * @param d2 The vector of distances in y.
 * @param d3 The vector of distances in z.
 * @param dens The vector of density values.
 * @param sf The vector of scattering factors.
 * @param file The output stream to write the results.
 * @param start The start time of the calculation.
 * @param end1 The end time of the calculation.
 * @param debug Flag indicating whether to enable debug mode.
 * @param no_date Flag indicating whether to exclude the date in the output.
 */
void calc_SF(const int& points,
	vec2& k_pt,
	vec2& d1,
	vec2& d2,
	vec2& d3,
	vec2& dens,
	cvec2& sf,
	std::ostream& file,
	_time_point& start,
	_time_point& end1,
	bool debug,
	bool no_date,
	bool do_XCW)
{
	const long long int imax = static_cast<long long int>(dens.size());
	const long long int smax = static_cast<long long int>(k_pt[0].size());
	sf.reserve(imax * smax);
	sf.resize(imax);
#pragma omp parallel for
	for (int i = 0; i < imax; i++)
		sf[i].resize(smax, constants::cnull);
	using namespace std;
	if (debug)
		file << "Initialized FFs" << std::endl
		<< "asym atom list size: " << imax << " total grid size: " << points << endl;
	end1 = get_time();

	if (!no_date)
	{
		long long int dur = get_sec(start, end1);
		if (dur < 1)
			file << "Time to prepare: " << fixed << setprecision(0) << get_msec(start, end1) << " ms" << endl << endl;
		else
			file << "Time to prepare: " << fixed << setprecision(0) << dur << " s" << endl << endl;
	}
	ProgressBar* progress = nullptr;
	if (!do_XCW) {
		progress = new ProgressBar(imax, 60, "=", " ", "Calculating Scattering Factors");
	}
	long long int pmax, p, s;
	complex<double>* sf_local;
	double work, rho, c, si, * dens_local, re, im, * d1_local, * d2_local, * d3_local;

	// Pre-fetch k_pt data pointers for better cache locality
	const double* k1_data = k_pt[0].data();
	const double* k2_data = k_pt[1].data();
	const double* k3_data = k_pt[2].data();

	for (int i = 0; i < imax; i++)
	{
		pmax = static_cast<long long int>(dens[i].size());
		dens_local = dens[i].data();
		d1_local = d1[i].data();
		d2_local = d2[i].data();
		d3_local = d3[i].data();

#pragma omp parallel for private(work, rho, c, si, re, im, s, p)
		for (s = 0; s < smax; s++)
		{
			re = 0.0, im = 0.0;
			const double& k1_local = k1_data[s];
			const double& k2_local = k2_data[s];
			const double& k3_local = k3_data[s];
			sf_local = sf[i].data();
			// Process loop in blocks of 4 for better instruction-level parallelism
			const long long int pmax_vec = (pmax / 4) * 4;

			// Vectorized main loop processing 4 elements at a time
			for (p = 0; p < pmax_vec; p += 4)
			{
				// Load 4 density values
				const double rho0 = dens_local[p];
				const double rho1 = dens_local[p + 1];
				const double rho2 = dens_local[p + 2];
				const double rho3 = dens_local[p + 3];

				// Calculate work values for 4 points using FMA pattern
				const double work0 = k1_local * d1_local[p] + k2_local * d2_local[p] + k3_local * d3_local[p];
				const double work1 = k1_local * d1_local[p + 1] + k2_local * d2_local[p + 1] + k3_local * d3_local[p + 1];
				const double work2 = k1_local * d1_local[p + 2] + k2_local * d2_local[p + 2] + k3_local * d3_local[p + 2];
				const double work3 = k1_local * d1_local[p + 3] + k2_local * d2_local[p + 3] + k3_local * d3_local[p + 3];

#if (defined(__GNUC__) || defined(__clang__)) && !defined(__APPLE__)
				double si0, c0, si1, c1, si2, c2, si3, c3;
				sincos(work0, &si0, &c0);
				sincos(work1, &si1, &c1);
				sincos(work2, &si2, &c2);
				sincos(work3, &si3, &c3);

				re += rho0 * c0 + rho1 * c1 + rho2 * c2 + rho3 * c3;
				im += rho0 * si0 + rho1 * si1 + rho2 * si2 + rho3 * si3;
#elif defined(__APPLE__)
				double si0, c0, si1, c1, si2, c2, si3, c3;
				__sincos(work0, &si0, &c0);
				__sincos(work1, &si1, &c1);
				__sincos(work2, &si2, &c2);
				__sincos(work3, &si3, &c3);

				re += rho0 * c0 + rho1 * c1 + rho2 * c2 + rho3 * c3;
				im += rho0 * si0 + rho1 * si1 + rho2 * si2 + rho3 * si3;
#else
				const double c0 = cos(work0);
				const double si0 = sin(work0);
				const double c1 = cos(work1);
				const double si1 = sin(work1);
				const double c2 = cos(work2);
				const double si2 = sin(work2);
				const double c3 = cos(work3);
				const double si3 = sin(work3);

				re += rho0 * c0 + rho1 * c1 + rho2 * c2 + rho3 * c3;
				im += rho0 * si0 + rho1 * si1 + rho2 * si2 + rho3 * si3;
#endif
			}

			// Handle remaining elements
			for (p = pmax_vec; p < pmax; p++)
			{
				rho = dens_local[p];
				work = k1_local * d1_local[p] + k2_local * d2_local[p] + k3_local * d3_local[p];
#if (defined(__GNUC__) || defined(__clang__)) && !defined(__APPLE__)
				sincos(work, &si, &c);
				re += rho * c;
				im += rho * si;
#elif defined(__APPLE__)
				__sincos(work, &si, &c);
				re += rho * c;
				im += rho * si;
#else
				c = cos(work);
				si = sin(work);
				re += rho * c;
				im += rho * si;
#endif
			}
			sf_local[s].real(re);
			sf_local[s].imag(im);
		}
		if (!do_XCW) {
			progress->update();
		}
	}
	if (!do_XCW) {
		delete (progress);
	}
}

void calc_SF_CUDA(const int& points,
	vec2& k_pt,
	vec2& d1,
	vec2& d2,
	vec2& d3,
	vec2& dens,
	cvec2& sf,
	std::ostream& file,
	_time_point& start,
	_time_point& end1,
	bool debug,
	bool no_date) {
	const long long int imax = static_cast<long long int>(dens.size());
	const long long int smax = static_cast<long long int>(k_pt[0].size());
	sf.reserve(imax * smax);
	sf.resize(imax);
#pragma omp parallel for
	for (int i = 0; i < imax; i++)
		sf[i].resize(smax, constants::cnull);
	using namespace std;
	if (debug)
		file << "Initialized FFs" << std::endl
		<< "asym atom list size: " << imax << " total grid size: " << points << endl;
	end1 = get_time();

	if (!no_date)
	{
		long long int dur = get_sec(start, end1);
		if (dur < 1)
			file << "Time to prepare: " << fixed << setprecision(0) << get_msec(start, end1) << " ms" << endl
			<< endl;
		else
			file << "Time to prepare: " << fixed << setprecision(0) << dur << " s" << endl
			<< endl;
	}
#ifdef __CUDACC__
	double** gpu_k_pt = NULL,
		** gpu_sf_r = NULL,
		** gpu_sf_i = NULL;
	vector<double> long_kpt;
	long_kpt.resize(3 * k_pt_unique[0].size());
	for (int i = 0; i < k_pt_unique[0].size(); i++)
	{
		long_kpt[3 * i + 0] = k_pt_unique[0][i];
		long_kpt[3 * i + 1] = k_pt_unique[1][i];
		long_kpt[3 * i + 2] = k_pt_unique[2][i];
	}
	gpu_k_pt = (double**)malloc(sizeof(double*));
	gpu_sf_r = (double**)malloc(sizeof(double*));
	gpu_sf_i = (double**)malloc(sizeof(double*));
	cudaMalloc((void**)&gpu_k_pt[0], sizeof(double) * k_pt_unique[0].size() * 3);
	cudaMalloc((void**)&gpu_sf_r[0][i], sizeof(double) * asym_atom_list.size() * k_pt_unique[0].size());
	cudaMalloc((void**)&gpu_sf_i[0][i], sizeof(double) * asym_atom_list.size() * k_pt_unique[0].size());
	cudaMemcpy(gpu_k_pt[0], long_kpt.data(), sizeof(double) * k_pt_unique[0].size() * 3, cudaMemcpyHostToDevice);

	dim3 blocks(asym_atom_list.size(), k_pt_unique[0].size());
	gpu_make_sf << <blocks, 1 >> > (
		gpu_sf_r[0],
		gpu_sf_i[0],
		gpu_k_pt[0],

		);
#endif
}

/**
 * Adds the ECP (Effective Core Potential) contribution to the scattering factors.
 *
 * @param asym_atom_list The list of asymmetric atoms.
 * @param wave The WFN (Wave Function) object.
 * @param sf The scattering factors.
 * @param cell The cell object.
 * @param hkl The hkl_list object.
 * @param file The output stream for writing the results.
 * @param mode The mode of operation. 0 = Gaussian tight core function, 1,2,3 = Thakkar core density based on the ECP type used
 * @param debug Flag indicating whether to enable debug mode.
 */
static void add_ECP_contribution(const ivec& asym_atom_list,
	const WFN& wave,
	cvec2& sf,
	const cell& cell,
	hkl_list& hkl,
	std::ostream& file,
	const int& mode,
	const bool debug)
{
	using namespace std;
	double k = 1.0;
	hkl_list_it it = hkl.begin();
	err_checkf(mode >= 0, "Invalid mode for ECP contribution!", file);
	if (mode == 0)
	{ // Using a gaussian tight core function
		if (debug)
		{
			file << "Using a gaussian tight core function" << endl;
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				if (wave.get_atom_ECP_electrons(asym_atom_list[i]) != 0)
					file << "Atom nr: " << wave.get_atom_charge(asym_atom_list[i]) << " core f000: "
					<< scientific << setw(14) << setprecision(8)
					<< wave.get_atom_ECP_electrons(asym_atom_list[i])
					<< " and at 1 angstrom: " << exp(-pow(constants::bohr2ang(k), 2) / 16.0 / constants::PI) * wave.get_atom_ECP_electrons(asym_atom_list[i]) << endl;
			}
		}
#pragma omp parallel for private(it, k)
		for (int s = 0; s < sf[0].size(); s++)
		{
			it = next(hkl.begin(), s);
			k = constants::FOUR_PI * cell.get_stl_of_hkl(*it);
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				sf[i][s] += wave.get_atom_ECP_electrons(asym_atom_list[i]) * exp(-k / 16.0 / constants::PI);
			}
		}
	}
	else if (mode == 1 || mode == 2 || mode == 3)
	{ // Using a Thakkar core density
		if (debug)
			file << "Using a Thakkar core density" << endl;
		//vector<Thakkar> temp;
		map<int, Thakkar> temp;
		map<int, Spherical_Gaussian_Density> temp_G;
		if (debug) {
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				const int charge = wave.get_atom_charge(asym_atom_list[i]);
				if (temp.find(charge) == temp.end()) {
					temp.emplace(charge, Thakkar{ charge, mode });
					temp_G.emplace(charge, Spherical_Gaussian_Density{ charge, mode });
					if (wave.get_atom_ECP_electrons(asym_atom_list[i]) != 0)
					{
						double k_0001 = temp[charge].get_core_form_factor(0, wave.get_atom_ECP_electrons(asym_atom_list[i]));
						double k_1 = temp[charge].get_core_form_factor(constants::FOUR_PI * constants::bohr2ang(1.0), wave.get_atom_ECP_electrons(asym_atom_list[i]));
						file << "Atom nr: " << charge << " number of ECP electrons: " << wave.get_atom_ECP_electrons(asym_atom_list[i]) << " core f(0) : "
							<< scientific << setw(14) << setprecision(8) << k_0001 << " and at 1 Ang: " << k_1 << endl;
					}
				}
			}
		}
		else {
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				const int charge = wave.get_atom_charge(asym_atom_list[i]);
				if (temp.find(charge) == temp.end()) {
					temp.emplace(charge, Thakkar{ charge, mode });
					temp_G.emplace(charge, Spherical_Gaussian_Density{ charge, mode });
				}
			}
		}

#pragma omp parallel for private(it, k)
		for (int s = 0; s < sf[0].size(); s++)
		{
			int n_el = 0, charge = 0;
			it = next(hkl.begin(), s);
			k = constants::FOUR_PI * constants::bohr2ang(cell.get_stl_of_hkl(*it));
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				n_el = wave.get_atom_ECP_electrons(asym_atom_list[i]);
				charge = wave.get_atom_charge(asym_atom_list[i]);
				if (n_el != 0)
				{
					sf[i][s] += temp_G.at(charge).get_form_factor(k); // This bit will correct for the error of the valence denisty of ECP atoms
					sf[i][s] += temp.at(charge).get_core_form_factor(k, n_el);
				}
			}
		}
	}
	else
	{
		err_not_impl_f("No higher ECP mode than 3 implemented!", file);
	}
}

/**
 * Converts the given asymmetric atom list, wavefunction, unit cell, and hkl list
 * to electron scattering factors using the given X-ray scattering factors.
 * This will be performed in place
 *
 * @param asym_atom_list The list of asymmetric atoms.
 * @param wave The wavefunction.
 * @param sf The scattering factors.
 * @param unit_cell The unit cell.
 * @param hkl The hkl list.
 */
void convert_to_ED(const ivec& asym_atom_list,
	const WFN& wave,
	cvec2& sf,
	const cell& unit_cell,
	const hkl_list& hkl)
{
    const std::vector<i3> hkl_vector(hkl.begin(), hkl.end());
    const int hkl_size = hkl.size();
#pragma omp parallel for shared(hkl_vector)
    for (int s = 0; s < hkl_size; s++)
    {
        const double h2 = pow(unit_cell.get_stl_of_hkl(hkl_vector[s]), 2);
        for (int i = 0; i < asym_atom_list.size(); i++)
            sf[i][s] = cdouble(constants::ED_fact * (wave.get_atom_charge(asym_atom_list[i]) - sf[i][s].real()) / h2, -constants::ED_fact * sf[i][s].imag() / h2);
    }
}


int make_atomic_grids_wrapper(
	const WFN& wave, const bvec& needs_grid, const ivec& asym_atom_list, const cell& unit_cell, const svec& labels, //
	std::vector<_time_point>& time_points, svec& time_descriptions, vec2& d1, vec2& d2, vec2& d3, vec2& dens,
	const options& opt) {

	const int atoms_with_grids = vec_sum(needs_grid);
	err_checkf(atoms_with_grids > 0, "No atoms with grids to generate!", std::cout);
	err_checkf(atoms_with_grids <= wave.get_ncen(), "More atoms with grids than in the wavefunction! Aborting!", std::cout);
	err_checkf(atoms_with_grids == asym_atom_list.size(), "Number of atoms with grids does not match the number of atoms in the CIF file!", std::cout);
	std::cout << "There are:\n"
		<< std::setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
		//<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
		<< std::setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << std::endl;


	std::cout << "\nSelected accuracy: " << opt.accuracy << "\nMaking Integration Grids..." << std::endl;

	GridConfiguration config;
	config.accuracy = opt.accuracy;
	config.partition_type = opt.partition_type;
	config.pbc = opt.pbc;
	config.debug = opt.debug;
	config.all_charges = opt.all_charges;

	GridManager grid_manager(config);

	WFN temp = wave;
	temp.delete_unoccupied_MOs();

	// Setup grids for the molecule
	grid_manager.setup3DGridsForMolecule(temp, asym_atom_list, needs_grid, unit_cell, opt.get_g);
	grid_manager.addTimingInfoToVecs(time_points, time_descriptions);


	// Calculate partitioned charges
	PartitionResults results = grid_manager.calculatePartitionedCharges(temp, unit_cell);
	grid_manager.printChargeTable(labels, temp, asym_atom_list, std::cout, results);
	time_points.push_back(get_time());
	time_descriptions.push_back("calculate charges");

	grid_manager.getDensityVectors(temp, asym_atom_list, d1, d2, d3, dens, opt.get_g);
	time_points.push_back(get_time());
	time_descriptions.push_back("combined density vectors");

	return grid_manager.getTotalGridPoints();
}

itsc_block calculate_scattering_factors_from_cube(
	options& opt,
	WFN& wave,
	const cube& density_cube,
	std::ostream& file)
{
	using namespace std;
	err_checkf(opt.cif != "", "Cube-density SF calculation requires -cif.", file);
	err_checkf(filesystem::exists(opt.cif), "CIF " + opt.cif.string() + " does not exists!", file);
	err_checkf(!opt.groups.empty() && !opt.groups[0].empty(), "No CIF disorder group selected. Use -group 0 for the default group.", file);
	err_checkf(wave.get_ncen() != 0, "Cube file does not contain atom definitions in the header.", file);

	vector<_time_point> time_points;
	vector<string> time_descriptions;
	time_points.push_back(get_time());

	cell unit_cell(opt.cif, file, opt.debug);
	ifstream cif_input(opt.cif.c_str(), ios::in);
	ivec atom_type_list;
	ivec asym_atom_to_type_list;
	ivec asym_atom_list;
	bvec needs_grid(wave.get_ncen(), false);
	svec known_atoms;

	auto labels = read_atoms_from_CIF(cif_input,
		opt.groups[0],
		unit_cell,
		wave,
		known_atoms,
		atom_type_list,
		asym_atom_to_type_list,
		asym_atom_list,
		needs_grid,
		file,
		opt.debug);

	cif_input.close();
	time_points.push_back(get_time());
	time_descriptions.push_back("cif reading");

	vec2 k_pt;
	hkl_list hkl;
	if (opt.m_hkl_list.size() != 0)
	{
		hkl = opt.m_hkl_list;
	}
	else if (opt.read_k_pts == false)
	{
		if (opt.dmin != 99.0)
		{
			if (opt.electron_diffraction)
				generate_hkl(opt.dmin / 2.0, hkl, opt.twin_law, unit_cell, file, opt.debug);
			else
				generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug);
		}
		else if (opt.hkl_min_max[0][0] != -100 && opt.hkl_min_max[2][1] != 100)
		{
			generate_hkl(opt.hkl_min_max, hkl, opt.twin_law, unit_cell, file, opt.debug, opt.electron_diffraction);
		}
		else
		{
			read_hkl(opt.hkl, hkl, opt.twin_law, unit_cell, file, opt.debug);
		}
		opt.m_hkl_list = hkl;
	}

	make_k_pts(
		hkl.size() == 0,
		opt.save_k_pts,
		unit_cell,
		hkl,
		k_pt,
		file,
		opt.debug);

	time_points.push_back(get_time());
	time_descriptions.push_back("k-points preparation");

	GridConfiguration config;
	config.accuracy = opt.accuracy;
	config.partition_type = opt.partition_type;
	config.pbc = opt.pbc;
	config.debug = opt.debug;
	config.all_charges = opt.all_charges;
	config.no_density_eval = true;

	GridManager grid_manager(config);
	grid_manager.setup3DGridsForMolecule(wave, asym_atom_list, needs_grid, unit_cell, false);
	grid_manager.addTimingInfoToVecs(time_points, time_descriptions);

	vec2 d1, d2, d3, dens;
	vec atom_electrons;
	grid_manager.getDensityVectorsFromCube(wave, asym_atom_list, density_cube, d1, d2, d3, dens, atom_electrons);
	time_points.push_back(get_time());
	time_descriptions.push_back("combined density vectors");

	file << "Table of Charges in electrons\n"
		<< setw(10) << "Atom" << setw(12) << "CubePart" << endl;
	for (int i = 0; i < labels.size(); i++)
	{
		const int atom_index = asym_atom_list[i];
		file << setw(10) << labels[i]
			<< fixed << setw(12) << setprecision(3) << wave.get_atom_charge(atom_index) - atom_electrons[i] << endl;
	}
	const double electron_sum = reduce(atom_electrons.begin(), atom_electrons.end(), 0.0);
	file << setprecision(4) << "Total number of partitioned electrons from cube: " << electron_sum << endl;
	time_points.push_back(get_time());
	time_descriptions.push_back("calculate charges");

	cvec2 sf;
	_time_point end1;
	calc_SF(grid_manager.getTotalGridPoints(),
		k_pt,
		d1, d2, d3, dens,
		sf,
		file,
		time_points.front(),
		end1,
		opt.debug,
		opt.no_date);
	time_points.push_back(get_time());
	time_descriptions.push_back("Fourier transform");

	if (opt.electron_diffraction)
	{
		convert_to_ED(asym_atom_list, wave, sf, unit_cell, hkl);
	}

	itsc_block blocky;
	if (opt.label_tsc_output)
	{
		blocky = itsc_block(sf, labels, hkl);
	}
	else
	{
		vector<atomID> ids(asym_atom_list.size());
		for (int atm_idx = 0; atm_idx < asym_atom_list.size(); atm_idx++)
		{
			ids[atm_idx] = wave.get_id_for_atom(asym_atom_list[atm_idx]);
		}
		blocky = itsc_block(sf, ids, hkl);
	}
	time_points.push_back(get_time());
	time_descriptions.push_back("tsc calculation");

	if (!opt.no_date)
	{
		write_timing_to_file(file, time_points, time_descriptions);
	}

	return blocky;
}

/**
 * @brief Calculates scattering factors using the provided calculator and options.
 *
 * This templated function computes scattering factors for a set of atoms and k-points.
 *
 * @tparam tsc_block_type The type used to store the resulting scattering factors block.
 *         Must support assignment and access as used in the function.
 * @tparam calculator_type The type of the calculator object.
 *         Must provide the required interface for calculating scattering factors.
 * @param opt Options for the calculation.
 * @param calculator The calculator object used to compute scattering factors.
 * @param file Output stream for logging messages.
 * @param known_atoms List of known atom symbols.
 * @param nr Number of atoms.
 * @param kpts Pointer to the k-points vector.
 * @return A block of scattering factors of type tsc_block_type.
 */
template <typename tsc_block_type, typename calculator_type>
tsc_block_type calculate_scattering_factors(
	options& opt,
	calculator_type calculator,
	std::ostream& file,
	svec& known_atoms,
	const int& nr,
	vec2* kpts,
	salted_part_prep* prep_out
) {
	using namespace std;
	int nat = 0;
	WFN* wavy = NULL;
	if constexpr (std::is_same_v<calculator_type, std::vector<WFN> &>) {
		wavy = &calculator[nr];
		wavy->delete_Qs();
		err_checkf(wavy->get_ncen() != 0, "No Atoms in the wavefunction, this will not work!!ABORTING!!", file);
		if (!opt.cif_based_combined_tsc_calc)
		{
			err_checkf(filesystem::exists(opt.cif), "CIF " + opt.cif.string() + " does not exists!", file);
		}
		else
		{
			for (int i = 0; i < opt.combined_tsc_calc_cifs.size(); i++)
			{
				err_checkf(filesystem::exists(opt.combined_tsc_calc_cifs[i]), "CIF " + opt.combined_tsc_calc_cifs[i].string() + " does not exists!", file);
			}
		}
		err_checkf(opt.groups[nr].size() >= 1, "Not enough groups specified to work with!", file);
		file << "Number of protons: " << wavy->get_nr_electrons() << endl
			<< "Number of electrons: " << fixed << wavy->count_nr_electrons() << endl;
		if (wavy->get_has_ECPs())
			file << "Number of ECP electrons: " << wavy->get_nr_ECP_electrons() << endl;
		// err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);
		if (opt.debug)
			file << "Working with: " << wavy->get_path() << endl;
		nat = wavy->get_ncen();
	}
	else if constexpr (std::is_same_v<calculator_type, SALTEDPredictor&>) {
		wavy = &calculator.wavy;
		nat = wavy->get_ncen();
	}

	vector<_time_point> time_points;
	vector<string> time_descriptions;
	time_points.push_back(get_time());

	filesystem::path cif;
	if (opt.cif_based_combined_tsc_calc)
		cif = opt.combined_tsc_calc_cifs[nr];
	else
		cif = opt.cif;


	cell unit_cell(cif, file, opt.debug);
	ifstream cif_input(cif.c_str(), std::ios::in);
	ivec atom_type_list;
	ivec asym_atom_to_type_list;
	ivec asym_atom_list;
	bvec needs_grid(nat, false);
	if (opt.debug)
		file << "Reading atoms!!!!" << endl;

	auto labels = read_atoms_from_CIF(cif_input,
		opt.groups[nr],
		unit_cell,
		*wavy,
		known_atoms,
		atom_type_list,
		asym_atom_to_type_list,
		asym_atom_list,
		needs_grid,
		file,
		opt.debug);

	cif_input.close();

	if (opt.debug)
		file << "There are " << atom_type_list.size() << " Types of atoms and " << asym_atom_to_type_list.size() << " atoms in total" << endl;

	time_points.push_back(get_time());
	time_descriptions.push_back("cif reading");

	vec2 k_pt;
	hkl_list hkl;
	if (opt.m_hkl_list.size() != 0)
	{
		hkl = opt.m_hkl_list;
	}
	else if (nr == 0 && opt.read_k_pts == false)
	{
		if (opt.dmin != 99.0)
			if (opt.electron_diffraction)
				generate_hkl(opt.dmin / 2.0, hkl, opt.twin_law, unit_cell, file, opt.debug);
			else
				generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug);
		else if (opt.hkl_min_max[0][0] != -100 && opt.hkl_min_max[2][1] != 100)
			generate_hkl(opt.hkl_min_max, hkl, opt.twin_law, unit_cell, file, opt.debug, opt.electron_diffraction);
		else
			read_hkl(opt.hkl, hkl, opt.twin_law, unit_cell, file, opt.debug);
		opt.m_hkl_list = hkl;
	}
	if (kpts == NULL || kpts->size() == 0)
	{
		make_k_pts(
			nr != 0 && hkl.size() == 0,
			opt.save_k_pts,
			unit_cell,
			hkl,
			k_pt,
			file,
			opt.debug);
		if (kpts != NULL)
		{
			*kpts = k_pt;
		}
	}
	else
	{
		k_pt = *kpts;
	}

	time_points.push_back(get_time());
	time_descriptions.push_back("k-points preparation");
	// Streaming emits the table in reflection blocks instead of holding
	// scatterers x reflections x 16 bytes at once. Only the plain single-part
	// SALTED path is routed this way so far; -mtc still needs every reflection
	// resident because parts are combined by scatterer.
	const bool stream_tsc = opt.tsc_block_size > 0
		&& !opt.needs_Thakkar_fill
		&& !opt.iam_switch
		&& !opt.electron_diffraction
		&& opt.combined_tsc_calc_files.size() <= 1;
	cvec2 sf;
	if (!stream_tsc)
	{
		sf.resize(asym_atom_list.size());
#pragma omp parallel for
		for (int i = 0; i < asym_atom_list.size(); i++)
			sf[i].resize(hkl.size());
	}

    if (opt.iam_switch) {
        if (prep_out != NULL)
        {
            // hand back what a per-block spherical calculation needs and stop
            const std::vector<i3> hv(hkl.begin(), hkl.end());
            prep_out->asym_atom_list = asym_atom_list;
            prep_out->labels = labels;
            prep_out->atom_type_list = atom_type_list;
            prep_out->asym_atom_to_type_list = asym_atom_to_type_list;
            prep_out->hkl_v = hv;
            prep_out->k_of_reflection.resize(hv.size());
            for (size_t s = 0; s < hv.size(); s++)
                prep_out->k_of_reflection[s] =
                    constants::bohr2ang(constants::FOUR_PI * unit_cell.get_stl_of_hkl(hv[s]));
            return tsc_block_type();
        }
        vector<Thakkar> spherical_atoms;
        spherical_atoms.reserve(atom_type_list.size());
        for (int i = 0; i < atom_type_list.size(); i++)
            spherical_atoms.emplace_back(atom_type_list[i]);

        const int imax = (int)asym_atom_list.size();
        const std::vector<i3> hkl_vector(hkl.begin(), hkl.end());
        const int hkl_max = hkl.size();

        if (!opt.electron_diffraction)
        {
#pragma omp parallel for shared(hkl_vector)
            for (int s = 0; s < hkl_max; s++)
            {
                const double k = constants::bohr2ang(constants::FOUR_PI * unit_cell.get_stl_of_hkl(hkl_vector[s]));
                for (int i = 0; i < imax; i++)
                    sf[i][s] = spherical_atoms[asym_atom_to_type_list[i]].get_form_factor(k);
            }
        }
        else
        {
#pragma omp parallel for shared(hkl_vector)
            for (int s = 0; s < hkl_max; s++)
            {
                const double stl = unit_cell.get_stl_of_hkl(hkl_vector[s]);
                const double k = constants::bohr2ang(constants::FOUR_PI * stl);
                const double h2 = pow(stl, 2);
                double sf_x = 0;
                for (int i = 0; i < imax; i++) {
                    sf_x = spherical_atoms[asym_atom_to_type_list[i]].get_form_factor(k);
                    sf[i][s] = cdouble(constants::ED_fact * (atom_type_list[asym_atom_to_type_list[i]] - sf_x) / h2, 0);
                }
            }
        }
    }
    else if constexpr (std::is_same_v<calculator_type, SALTEDPredictor &>)
    {
        // Generation of SALTED density coefficients
        file << "\nGenerating densities... " << endl;
        vec coefs = calculator.gen_SALTED_densities();
        file << setw(13 * 4) << "... done!" << endl;
        time_points.push_back(get_time());
        time_descriptions.push_back("SALTED prediction");
        calculator.shrink_intermediate_vectors();

		err_checkf(labels.size() == asym_atom_list.size(),
			"Inconsistent SALTED atom bookkeeping after disorder filtering!", file);

		vec atom_elecs = calc_atomic_density(calculator.wavy.get_atoms(), coefs);
		file << "Table of Charges in electrons\n"
			<< "       Atom      ML" << endl;

		for (int i = 0; i < labels.size(); i++)
		{
			const int atom_index = asym_atom_list[i];
			file << setw(10) << labels[i]
				<< fixed << setw(10) << setprecision(3) << wavy->get_atom_charge(atom_index) - atom_elecs[atom_index];
			if (opt.debug)
				file << " " << setw(4) << wavy->get_atom_charge(atom_index) << " " << fixed << setw(10) << setprecision(3) << atom_elecs[atom_index];
			file << endl;
		}
		auto el_sum = reduce(atom_elecs.begin(), atom_elecs.end(), 0.0);
		file << setprecision(4) << "Total number of analytical Electrons: " << el_sum << endl;
		time_points.push_back(get_time());
		time_descriptions.push_back("Calculation of Charges");

		if (prep_out != NULL)
		{
			// -mtc streaming: the reflection loop lives outside this call, so hand
			// back everything that does not depend on reflections and stop here.
			prep_out->coefs = std::move(coefs);
			prep_out->asym_atom_list = asym_atom_list;
			prep_out->labels = labels;
			prep_out->atoms = calculator.wavy.get_atoms_ptr();
			prep_out->k_pt = k_pt;
			prep_out->hkl_v.assign(hkl.begin(), hkl.end());
			return tsc_block_type();
		}

		if (stream_tsc)
		{
			// calc_SF_SALTED sizes its output from the k-points handed to it, so a
			// sliced k_pt needs no change on its side; the chunking is a loop around
			// an unmodified kernel.
			ScattererLabels stream_ids;
			if (opt.label_tsc_output)
				for (const auto& l : labels) stream_ids.emplace_back(l);
			else
				for (size_t a = 0; a < asym_atom_list.size(); a++)
					stream_ids.emplace_back(wavy->get_id_for_atom(asym_atom_list[a]));

			const std::vector<i3> hkl_v(hkl.begin(), hkl.end());
			const size_t n_refl = hkl_v.size();
			const size_t bs = std::min(opt.tsc_block_size, n_refl ? n_refl : 1);
			file << "Streaming tsc in blocks of " << bs << " reflections" << endl;
			tsc_stream_writer<int, cdouble> writer(
				"experimental.tscb", stream_ids, std::string(), n_refl, 2);
			size_t block_id = 0;
			for (size_t lo = 0; lo < n_refl; lo += bs)
			{
				const size_t hi = std::min(lo + bs, n_refl);
				vec2 k_slice(3, vec(hi - lo));
				std::vector<std::vector<int>> idx(3, std::vector<int>(hi - lo));
				for (size_t r = lo; r < hi; r++)
					for (int dm = 0; dm < 3; dm++)
					{
						k_slice[dm][r - lo] = k_pt[dm][r];
						idx[dm][r - lo] = hkl_v[r][dm];
					}
				cvec2 chunk;
				calc_SF_SALTED(k_slice, coefs, calculator.wavy.get_atoms(), asym_atom_list, chunk);
				writer.submit(block_id++, std::move(idx), std::move(chunk));
			}
			writer.finish();
			opt.tsc_written_by_stream = true;
		}
		else
		calc_SF_SALTED(
			k_pt,
			coefs,
			calculator.wavy.get_atoms(),
			asym_atom_list,
			sf);
		file << setw(13 * 4) << "... done!\n"
			<< flush;
		time_points.push_back(get_time());
		time_descriptions.push_back("Fourier transform");

		if (calculator.wavy.get_has_ECPs())
		{
			add_ECP_contribution(
				asym_atom_list,
				calculator.wavy,
				sf,
				unit_cell,
				hkl,
				file,
				opt.ECP_mode,
				opt.debug);
		}
	}
	else if constexpr (std::is_same_v<calculator_type, std::vector<WFN> &>) {
		if (opt.partition_type == PartitionType::Hirshfeld ||
			opt.partition_type == PartitionType::Becke ||
			opt.partition_type == PartitionType::TFVC ||
			opt.partition_type == PartitionType::MBIS ||
			opt.partition_type == PartitionType::EMBIS)
		{
			vec2 d1, d2, d3, dens;
			const int points = make_atomic_grids_wrapper(
				*wavy,
				needs_grid,
				asym_atom_list,
				unit_cell,
				labels,
				time_points,
				time_descriptions,
				d1, d2, d3, dens,
				opt);

			_time_point end1;
			calc_SF(points,
				k_pt,
				d1, d2, d3, dens,
				sf,
				file,
				time_points.front(),
				end1,
				opt.debug,
				opt.no_date);

			time_points.push_back(get_time());
			time_descriptions.push_back("Fourier transform");
		}
		else if (opt.partition_type == PartitionType::RI)
		{
			file << "\nGenerating densities... " << endl;
			WFN wavy_aux = generate_aux_wfn(*wavy, opt.aux_basis);

            //TODO: only compute coefs for atoms that are actually in the symmetric unit!
            DensityFitting::CONFIG config;
            config.analyze_quality = opt.debug;
            config.asym_atm_list = asym_atom_list;
            //config.restrain_type = DensityFitting::RESTRAINT_TYPE::SIMPLE_AND_TIK;
            //config.charge_scheme = DensityFitting::CHARGE_SCHEME::HIRSHFELD;
            //if (wavy->get_origin() == e_origin::ptb)
            //    config.restraint_strength = 1.0e-4;

            vec coefs = DensityFitting::density_fit(*wavy, wavy_aux, config);
            file << setw(12 * 4 + 2) << "... done!\n"
                << flush;
            time_points.push_back(get_time());
            time_descriptions.push_back("RI-Fit");

            vec atom_elecs = calc_atomic_density(wavy_aux.get_atoms(), coefs);
            file << "Table of Charges in electrons\n"
                << "       Atom  Charge_RI" << endl;

            for (int i = 0; i < asym_atom_list.size(); i++)
            {
                int a = asym_atom_list[i];
                file << setw(10) << labels[i]
                    << fixed << setw(10) << setprecision(3) << wavy_aux.get_atom_charge(a) - atom_elecs[a];
                if (opt.debug)
                    file << " " << setw(4) << wavy_aux.get_atom_charge(a) << " " << fixed << setw(10) << setprecision(3) << atom_elecs[a];
                file << endl;
            }

            auto el_sum = reduce(atom_elecs.begin(), atom_elecs.end(), 0.0);
            file << setprecision(4) << "Total number of analytical Electrons: " << el_sum << endl;
            time_points.push_back(get_time());
            time_descriptions.push_back("Calculation of Charges");

            if (stream_tsc)
            {
                // calc_SF_SALTED sizes its output from the k-points it is given,
                // so a sliced k_pt needs no change on its side.
                ScattererLabels stream_ids;
                if (opt.label_tsc_output)
                    for (const auto& l : labels) stream_ids.emplace_back(l);
                else
                    for (size_t a = 0; a < asym_atom_list.size(); a++)
                        stream_ids.emplace_back(wavy->get_id_for_atom(asym_atom_list[a]));

                const std::vector<i3> hkl_v(hkl.begin(), hkl.end());
                const size_t n_refl = hkl_v.size();
                const size_t bs = std::min(opt.tsc_block_size, n_refl ? n_refl : 1);
                file << "Streaming tsc in blocks of " << bs << " reflections" << endl;
                tsc_stream_writer<int, cdouble> writer(
                    opt.binary_tsc ? "experimental.tscb" : "experimental.tsc",
                    stream_ids, std::string(), n_refl, 2);
                size_t block_id = 0;
                for (size_t lo = 0; lo < n_refl; lo += bs)
                {
                    const size_t hi = std::min(lo + bs, n_refl);
                    vec2 k_slice(3, vec(hi - lo));
                    std::vector<std::vector<int>> idx(3, std::vector<int>(hi - lo));
                    for (size_t r = lo; r < hi; r++)
                        for (int dm = 0; dm < 3; dm++)
                        {
                            k_slice[dm][r - lo] = k_pt[dm][r];
                            idx[dm][r - lo] = hkl_v[r][dm];
                        }
                    cvec2 chunk;
                    calc_SF_SALTED(k_slice, coefs, wavy_aux.get_atoms(), asym_atom_list, chunk);
                    writer.submit(block_id++, std::move(idx), std::move(chunk));
                }
                writer.finish();
                opt.tsc_written_by_stream = true;
            }
            else
            {
            calc_SF_SALTED(
                k_pt,
                coefs,
                wavy_aux.get_atoms(),
                asym_atom_list,
                sf);
            }
            file << setw(12 * 4 + 2) << "... done!" << endl;
            time_points.push_back(get_time());
            time_descriptions.push_back("Fourier transform");
        }
        else {
            std::cout << "Unknown Partition type, stopping here!" << std::endl;
        }
        if (wavy->get_has_ECPs())
        {
            add_ECP_contribution(
                asym_atom_list,
                *wavy,
                sf,
                unit_cell,
                hkl,
                file,
                opt.ECP_mode,
                opt.debug);
        }
    }

    if (opt.electron_diffraction && !opt.iam_switch)
    {
        convert_to_ED(asym_atom_list,
            *wavy,
            sf,
            unit_cell,
            hkl);
    }

	tsc_block_type blocky;
	if (opt.label_tsc_output)
	{
		blocky = tsc_block_type(sf, labels, hkl);
	}
	else
	{
		std::vector<atomID> IDs(asym_atom_list.size());
		for (int atm_idx = 0; atm_idx < asym_atom_list.size(); atm_idx++) {
			IDs[atm_idx] = wavy->get_id_for_atom(asym_atom_list[atm_idx]);
		}
		blocky = tsc_block_type(sf, IDs, hkl);
	}


    if (opt.needs_Thakkar_fill)
    {
        file << "Performing the remaining calculation of spherical atoms..." << std::endl;
        opt.needs_Thakkar_fill = false;
        // nr is not only an index into this vector: it also picks opt.groups[nr],
        // which is what decides WHICH disorder parts the fill takes atoms from, and
        // the part's cif. Handing over a one-element vector and nr = 0 therefore
        // asked for part 0's atoms no matter which part was being filled, so an
        // atom belonging only to a later part was never restored - 3NIR came out
        // with 1025 scatterers where it has 1026. Build the whole list instead, so
        // the index means the same thing in all three places. (Passing nr with a
        // one-element vector is the variant that reads tempy[1] and crashes.)
        vector<WFN> tempy;
        int fill_nr = 0;
        if (!opt.wfn.empty()) {
            tempy.emplace_back(opt.wfn);
        }
        else {
            for (const auto& part_file : opt.combined_tsc_calc_files)
                tempy.emplace_back(part_file);
            fill_nr = nr;
        }
        opt.m_hkl_list = hkl;
        opt.iam_switch = true; opt.no_date = true;
        tsc_block<int, cdouble> blocky_thakkar = calculate_scattering_factors<itsc_block, std::vector<WFN> &>(opt, tempy, file, labels, fill_nr);
        opt.iam_switch = false; opt.no_date = false;
        blocky.append(std::move(blocky_thakkar), file);
        time_points.push_back(get_time());
        time_descriptions.push_back("Spherical Atoms");
    }


	time_points.push_back(get_time());
	time_descriptions.push_back("tsc calculation");

	if (!opt.no_date)
	{
		write_timing_to_file(file,
			time_points,
			time_descriptions);
	}

	return blocky;
}
// ---------------------------------------------------------------------------
// Streaming a combined (-mtc) table, one block of reflections at a time
// ---------------------------------------------------------------------------
//
// A disordered structure is calculated as several "parts". The ordinary path
// runs one part to completion, then the next, and merges the finished tables:
//
//     for each part:  compute EVERY reflection  ->  append to the table
//
// That needs the whole table in memory, because a part contributes rows and
// every row spans every reflection.
//
// This turns the loops inside out:
//
//     prepare every part once            (no reflections involved yet)
//     for each block of reflections:
//         for each part: compute just this block
//         hand the assembled block to the writer
//
// The inversion is only possible because everything expensive in a part - the
// CIF reading and the SALTED density prediction - does not depend on which
// reflections we ask for. Prepared once and kept, it costs nothing to reuse
// for every block. That is what salted_part_prep holds.
//
// THE SUBTLE PART: ROW ORDER.
//
// The rows of the finished file must appear in exactly the order the merging
// path would have produced, or every value is attributed to the wrong atom -
// a file of the right size, full of plausible, misassigned numbers, which no
// error message would ever reveal.
//
// It works out to a plain concatenation, for a reason worth stating: each part
// is prepared with the identifiers of the parts before it (the growing "known"
// list), exactly as the sequential code does. read_atoms_from_CIF skips any
// atom already in that list, so a part never claims an atom an earlier one
// covered and no duplicates arise. Prepare the parts against an empty list
// instead and they would overlap - silently, and wrongly.
//
// WHAT GOES IN THAT LIST MATTERS, and it is not obvious. The sequential path
// passes result.get_scatterers_string(), which is the atomID hex string unless
// labels were explicitly requested. Atom LABELS are not unique across disorder
// parts - the same label legitimately appears in several - so filling the list
// with labels makes a later part skip atoms it should have kept.
//
// That mistake was made here and cost a table that was 760 KB short on a
// four-part structure: the right size to look plausible, silently missing rows.
// Three-part 1EJG did not reveal it; four-part 3NIR did. If this is ever
// changed, check it against a structure with several parts and compare the
// bytes, not the size.
//
// Atoms of species the SALTED model does not know are erased from every part,
// so they are missing at this stage. The old path computed a whole second
// table of spherical (Thakkar) form factors for them and appended it; here
// they are simply extra rows of each block. A Thakkar factor depends only on
// the element and the reflection, so it needs no prediction and chunks freely.
// ---------------------------------------------------------------------------
bool stream_mtc_salted(options& opt, std::vector<WFN>& wavy, std::ostream& file, vec2* known_kpts)
{
	const size_t n_parts = opt.combined_tsc_calc_files.size();
	if (opt.tsc_block_size == 0 || !opt.SALTED || n_parts < 2
		|| opt.electron_diffraction || opt.iam_switch)
		return false;

	std::vector<std::shared_ptr<SALTEDPredictor>> preds;
	std::vector<salted_part_prep> preps;
	svec known;
	for (size_t i = 0; i < n_parts; i++)
	{
		auto pred = std::make_shared<SALTEDPredictor>(wavy[i], opt);
		if (!pred->basis_set_loaded())
			load_basis_into_WFN(pred->wavy, BasisSetLibrary::get_basis_set(pred->get_dfbasis_name()));
		salted_part_prep prep;
		calculate_scattering_factors<itsc_block, SALTEDPredictor&>(
			opt, *pred, file, known, static_cast<int>(i), known_kpts, &prep);
		// Feed the NEXT part exactly what the sequential path feeds it:
		// result.get_scatterers_string(), which is the atomID hex string unless
		// labels were requested. Labels are NOT unique across disorder parts - the
		// same label appears in several - so passing those makes a later part skip
		// atoms it should keep, and the table silently loses rows.
		for (size_t a = 0; a < prep.asym_atom_list.size(); a++)
		{
			if (opt.label_tsc_output)
				known.push_back(prep.labels[a]);
			else
				known.push_back(pred->wavy.get_id_for_atom(prep.asym_atom_list[a]).to_hex_string());
		}
		preds.push_back(pred);
		preps.push_back(std::move(prep));
	}

	// Atoms the model cannot predict - an unknown species, or nothing inside the
	// descriptor cutoff - were erased from every part, so they are still missing.
	// The old path recomputed a whole second table for them and appended it; here
	// they become extra rows of each block.
	//
	// ONE PASS OVER PART 0 IS NOT ENOUGH. A part file only yields the atoms that
	// belong to the parts it covers, so a missing atom that belongs to part 3 can
	// only come out of part 3's file. Filling from part 0 alone silently dropped
	// it: 3NIR came out with 1025 scatterers where it has 1026, which is exactly
	// the kind of shortfall that looks plausible in a file listing. Walk every
	// part, feeding `known` as we go so no atom is produced twice.
	std::vector<salted_part_prep> spherical(n_parts);
	std::vector<char> have_spherical(n_parts, 0);
	bool any_spherical = false;
	if (opt.needs_Thakkar_fill)
	{
		// Pin the reflections: the whole-table path set m_hkl_list before its
		// recursive fill for exactly this reason. Without it the spherical prep can
		// build its own list and k_of_reflection ends up a different length from the
		// parts, which is read off the end below.
		opt.m_hkl_list.clear();
		for (const auto& h : preps[0].hkl_v) opt.m_hkl_list.emplace(h);
		const bool iam_was = opt.iam_switch;
		opt.iam_switch = true;
		size_t n_filled = 0;
		// The whole list, not one file: the index also selects opt.groups[nr], so a
		// one-element vector with nr = i either reads past the end or asks the wrong
		// part for its atoms. See the note in the sequential path above.
		std::vector<WFN> tempy;
		for (const auto& part_file : opt.combined_tsc_calc_files)
			tempy.emplace_back(part_file);
		for (size_t i = 0; i < n_parts; i++)
		{
			calculate_scattering_factors<itsc_block, std::vector<WFN>&>(
				opt, tempy, file, known, static_cast<int>(i), known_kpts, &spherical[i]);
			have_spherical[i] = spherical[i].asym_atom_list.empty() ? 0 : 1;
			if (!have_spherical[i]) continue;
			any_spherical = true;
			n_filled += spherical[i].asym_atom_list.size();
			err_checkf(spherical[i].k_of_reflection.size() == preps[0].hkl_v.size(),
				"Spherical remainder of part " + std::to_string(i + 1) + " covers " +
				std::to_string(spherical[i].k_of_reflection.size()) +
				" reflections, the parts cover " + std::to_string(preps[0].hkl_v.size()), file);
			for (size_t a = 0; a < spherical[i].labels.size(); a++)
				known.push_back(spherical[i].labels[a]);
		}
		opt.iam_switch = iam_was;
		file << "Spherical remainder: " << n_filled
			 << " atom(s) the model cannot predict" << std::endl;
	}

	ScattererLabels ids;
	for (size_t p = 0; p < preps.size(); p++)
		for (size_t a = 0; a < preps[p].asym_atom_list.size(); a++)
		{
			if (opt.label_tsc_output) ids.emplace_back(preps[p].labels[a]);
			else ids.emplace_back(preds[p]->wavy.get_id_for_atom(preps[p].asym_atom_list[a]));
		}
	for (size_t p = 0; p < n_parts; p++)
		if (have_spherical[p])
			for (size_t a = 0; a < spherical[p].asym_atom_list.size(); a++)
				ids.emplace_back(spherical[p].labels[a]);

	const size_t n_refl = preps[0].hkl_v.size();
	const size_t bs = std::min(opt.tsc_block_size, n_refl ? n_refl : 1);
	file << "Streaming combined tsc: " << ids.size() << " scatterers, "
		<< n_refl << " reflections, blocks of " << bs << std::endl;

	tsc_stream_writer<int, cdouble> writer("experimental.tscb", ids, std::string(), n_refl, 2);
	size_t block_id = 0;
	for (size_t lo = 0; lo < n_refl; lo += bs)
	{
		const size_t hi = std::min(lo + bs, n_refl);
		std::vector<std::vector<int>> idx(3, std::vector<int>(hi - lo));
		for (size_t r = lo; r < hi; r++)
			for (int dm = 0; dm < 3; dm++)
				idx[dm][r - lo] = preps[0].hkl_v[r][dm];

		cvec2 combined;
		combined.reserve(ids.size());
		for (size_t p = 0; p < preps.size(); p++)
		{
			vec2 k_slice(3, vec(hi - lo));
			for (size_t r = lo; r < hi; r++)
				for (int dm = 0; dm < 3; dm++)
					k_slice[dm][r - lo] = preps[p].k_pt[dm][r];
			cvec2 chunk;
			calc_SF_SALTED(k_slice, preps[p].coefs, *preps[p].atoms, preps[p].asym_atom_list, chunk);
			for (auto& row : chunk) combined.push_back(std::move(row));
		}
		// Same order as the id list above: every part's spherical rows, in part order
		for (size_t p = 0; p < n_parts; p++)
		{
			if (!have_spherical[p]) continue;
			std::vector<Thakkar> spheres;
			spheres.reserve(spherical[p].atom_type_list.size());
			for (size_t t = 0; t < spherical[p].atom_type_list.size(); t++)
				spheres.emplace_back(spherical[p].atom_type_list[t]);
			for (size_t a = 0; a < spherical[p].asym_atom_list.size(); a++)
			{
				cvec row(hi - lo);
				for (size_t r = lo; r < hi; r++)
					row[r - lo] = spheres[spherical[p].asym_atom_to_type_list[a]]
						.get_form_factor(spherical[p].k_of_reflection[r]);
				combined.push_back(std::move(row));
			}
		}
		writer.submit(block_id++, std::move(idx), std::move(combined));
	}
	writer.finish();
	opt.tsc_written_by_stream = true;
	return true;
}

template itsc_block calculate_scattering_factors(options& opt,
	std::vector<WFN>& calculator,
	std::ostream& file,
	svec& known_atoms,
	const int& nr,
	vec2* kpts,
	salted_part_prep* prep_out);
template itsc_block calculate_scattering_factors(options& opt,
	SALTEDPredictor& calculator,
	std::ostream& file,
	svec& known_atoms,
	const int& nr,
	vec2* kpts,
	salted_part_prep* prep_out);


/**
 * Calculates the diffuse (that is non integer hkl) scattering factors based on the given options.
 *
 * @param opt The options for calculating the diffuse scattering factors.
 * @param log_file The output stream for logging the calculation process.
 */
void calc_sfac_diffuse(const options& opt, std::ostream& log_file)
{
	using namespace std;
	std::vector<WFN> wavy;
	wavy.emplace_back(e_origin::wfn);
	wavy[0].read_known_wavefunction_format(opt.wfn, std::cout, opt.debug);
	// set number of threads
#ifdef _OPENMP
	if (opt.threads > 0)
		omp_set_num_threads(opt.threads);
#endif
	err_checkf(wavy[0].get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
	err_checkf(exists(opt.cif), "CIF does not exists!", std::cout);
	// err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);

	// time_point start = get_time();
	// time_point end_becke, end_prototypes, end_spherical, end_prune, end_aspherical;
	vector<_time_point> time_points;
	vector<string> time_descriptions;
	time_points.push_back(get_time());

	cell unit_cell(opt.cif, std::cout, opt.debug);
	ifstream cif_input(opt.cif.c_str(), std::ios::in);
	ivec atom_type_list;
	ivec asym_atom_to_type_list;
	ivec asym_atom_list;
	bvec needs_grid(wavy[0].get_ncen(), false);
	svec known_atoms;

	auto labels = read_atoms_from_CIF(cif_input,
		opt.groups[0],
		unit_cell,
		wavy[0],
		known_atoms,
		atom_type_list,
		asym_atom_to_type_list,
		asym_atom_list,
		needs_grid,
		std::cout,
		opt.debug);

	cif_input.close();
	vec2 d1_, d2_, d3_, dens;

	make_atomic_grids_wrapper(
		wavy[0],
		needs_grid,
		asym_atom_list,
		unit_cell,
		labels,
		time_points,
		time_descriptions,
		d1_, d2_, d3_, dens, opt);

	hkl_list_d hkl;
	generate_fractional_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, log_file, opt.sfac_diffuse, opt.debug);

	const long long int size = static_cast<long long int>(hkl.size());
	vec2 k_pt;
	k_pt.reserve(3 * size);
	k_pt.resize(3);
#pragma omp parallel for
	for (int i = 0; i < 3; i++)
		k_pt[i].resize(size, 0.0);

	if (opt.debug)
	{
		log_file << "K_point_vector is here! size: " << k_pt[0].size() << endl;
	}
	int i_ = 0;
	for (const d3& hkl_ : hkl)
	{
		for (int x = 0; x < 3; x++)
		{
			for (int j = 0; j < 3; j++)
			{
				k_pt[x][i_] += unit_cell.get_rcm(x, j) * hkl_[j];
			}
		}
		i_++;
	}

	// below is a strip of Calc_SF without the file IO or progress bar
	cvec2 sf;

	const long long int imax = static_cast<long long int>(dens.size());
	const long long int smax = static_cast<long long int>(k_pt[0].size());
	long long int pmax = static_cast<long long int>(dens[0].size());
	//const long long int step = max(static_cast<long long int>(floor(smax / 20)), 1LL);
	std::cout << "Done with making k_pt " << smax << " " << imax << " " << pmax << endl;
	sf.reserve(imax * smax);
	sf.resize(imax);
#pragma omp parallel for
	for (int i = 0; i < imax; i++)
		sf[i].resize(smax);
	double* dens_local, * d1_local, * d2_local, * d3_local;
	complex<double>* sf_local;
	const double* k1_local = k_pt[0].data();
	const double* k2_local = k_pt[1].data();
	const double* k3_local = k_pt[2].data();
	double work, rho;
	ProgressBar* progress = new ProgressBar(imax * smax, 60, "=", " ", "Calculating Scattering Factors");
	for (long long int i = 0; i < imax; i++)
	{
		pmax = static_cast<long long int>(dens[i].size());
		dens_local = dens[i].data();
		d1_local = d1_[i].data();
		d2_local = d2_[i].data();
		d3_local = d3_[i].data();
		sf_local = sf[i].data();
#pragma omp parallel for private(work, rho)
		for (long long int s = 0; s < smax; s++)
		{
			double re = 0.0, im = 0.0, si, c;
			const double& _k1_local = k1_local[s];
			const double& _k2_local = k2_local[s];
			const double& _k3_local = k3_local[s];
			for (long long int p = pmax - 1; p >= 0; p--)
			{
				rho = dens_local[p];
				work = _k1_local * d1_local[p] + _k2_local * d2_local[p] + _k3_local * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
				if (rho < 0)
				{
					rho = -rho;
					work += M_PI;
				}
#endif
#endif
				c = cos(work);
				si = sin(work);
				re += rho * c;
				im += rho * si;
			}
			sf_local[s].real(re);
			sf_local[s].imag(im);
			progress->update();
		}
	}
	delete (progress);
	tsc_block<double, cdouble> result(sf, labels, hkl);
	result.write_tsc_file_non_integer(opt.cif);
}
