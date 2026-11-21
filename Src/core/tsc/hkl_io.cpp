#include "pch.h"
#include "scattering_factors.h"
#include <map>
#include "cell.h"

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
		const double scale = unit_cell.get_hkl_scale();

#pragma omp parallel for
		for (int ref = 0; ref < size; ref++)
		{
			const i3 hkl_ = hkl_vector[ref];
			for (int x = 0; x < 3; x++)
			{
				for (int j = 0; j < 3; j++)
				{
					k_pt[x][ref] += unit_cell.get_rcm(x, j) * hkl_[j] / scale;
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

//header/comment lines carry a letter or bracket; uppercase X Y Z are left out as before.
//strpbrk instead of the per-line std::regex, which made read_hkl most of rubredoxin's 0.9 s k-point phase (148k lines, i7-7700HQ)
static bool hkl_text_line(const std::string& line)
{
	return std::strpbrk(line.c_str(), "abcdefghijklmnopqrstuvwxyz()ABCDEFGHIJKLMNOPQRSTUVW") != nullptr;
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
	std::string line, temp;
	while (!hkl_input.eof())
	{
		getline_universal(hkl_input, line);
		if (hkl_input.eof())
			break;
		if (line.size() < 2)
			continue;
		if (hkl_text_line(line))
			continue;
		err_checkf(line.size() >= 12, "hkl line too short for h k l: '" + line + "'", file);
		// if (debug) file << "hkl: ";
		for (int i = 0; i < 3; i++)
		{
			temp = line.substr(4 * size_t(i), 4);
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
	err_checkf(!hkl.empty(), "No reflections read from " + hkl_filename.string(), file);
	file << " done!\nNr of reflections read from file: " << hkl.size() << std::endl;

	if (debug)
		file << "Number of reflections before twin: " << hkl.size() << std::endl;
	if (twin_law.size() > 0)
	{
		for (const i3& hkl__ : hkl)
			for (int i = 0; i < twin_law.size(); i++)
				hkl.emplace(i3{
					int(std::lround(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2])),
					int(std::lround(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2])),
					int(std::lround(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2])) });
	}
	if (debug)
		file << "Number of reflections after twin: " << hkl.size() << std::endl;

	ivec3 sym(3);
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
				file << "\n";
			}
			file << "\n";
		}
	}
	else
		file << "Number of symmetry operations: " << sym[0][0].size() << std::endl;

	//the list stays the file's reflections: a symmetry-enlarged, Friedel-pruned copy used to be built
	//here and kept only if get_sym() was empty, which it never is (cell always holds its 3 rows)
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
	scatter_data& obs,
	bool debug
)
{
	file << "Reading: " << std::setw(44) << hkl_filename << std::flush;
	i3 hkl_;
	double F_, abs_F_, F2_, sigma_, sigma2_;
	err_checkf(std::filesystem::exists(hkl_filename), "HKL file does not exists!", file);
	std::ifstream hkl_input(hkl_filename, std::ios::in);
	hkl_input.seekg(0, hkl_input.beg);
	std::string line, temp;
	//hkl is a set, so it hands the reflections back in (h,k,l) order whatever the file order;
	//obs has to follow that order, or F_calc[i] meets the wrong observation
	std::map<i3, scattering_data> obs_by_hkl;
	while (!hkl_input.eof())
	{
		getline_universal(hkl_input, line);
		if (hkl_input.eof())
			break;
		if (line.size() < 2)
			continue;
		if (hkl_text_line(line))
			continue;
		err_checkf(line.size() >= 12, "hkl line too short for h k l: '" + line + "'", file);
		// if (debug) file << "hkl: ";
		for (int i = 0; i < 3; i++)
		{
			temp = line.substr(4 * size_t(i), 4);
			temp.erase(remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
			hkl_[i] = stoi(temp);
			// if (debug) file << setw(4) << temp;
		}
		//F2 runs up to two digits past its decimal point, sigma is the rest: this also reads
		//files whose F2 field is wider than the 3I4,2F8.2 eight characters
		temp = line.substr(12);
		const size_t dot = temp.find_first_of('.');
		err_checkf(dot != std::string::npos, "hkl line without an F2 value: '" + line + "'", file);
		std::string temp_F = temp.substr(0, dot + 3);
		std::string temp_sigma = temp.substr(std::min(dot + 3, temp.size()));
		temp_sigma.erase(remove_if(temp_sigma.begin(), temp_sigma.end(), ::isspace), temp_sigma.end());
		sigma2_ = stof(temp_sigma);
		temp_F.erase(remove_if(temp_F.begin(), temp_F.end(), ::isspace), temp_F.end());
		F_ = stof(temp_F);
		//the two branches differ only in sign; they are NOT std::abs/copysign, which would turn a -0.00 F2 into +0.0 and flip sigma's infinity
		if (F_ < 0) {
			F2_ = -F_;
			abs_F_ = std::sqrt(F2_);
			F_ = -abs_F_;
			sigma_ = sigma2_ * 0.5 / -F_;
		}
		else {
			F2_ = F_;
			abs_F_ = std::sqrt(F2_);
			F_ = abs_F_;
			sigma_ = sigma2_ * 0.5 / F_;
		}
		// if (debug) file << endl;
		hkl.emplace(hkl_);
		obs_by_hkl.try_emplace(hkl_, scattering_data{ F_, abs_F_, F2_, sigma_, sigma2_ });
	}
	hkl_list_it found = hkl.find(i3{ 0, 0, 0 });
	if (found != hkl.end())
	{
		if (debug)
			file << "popping back 0 0 0" << std::endl;
		hkl.erase(i3{ 0, 0, 0 });
	}
	hkl_input.close();
	err_checkf(!hkl.empty(), "No reflections read from " + hkl_filename.string(), file);
	obs.abs_F_obs.clear();
	obs.F_obs.clear();
	obs.F_obs2.clear();
	obs.sigma_obs.clear();
	obs.sigma_obs2.clear();
	std::vector<scattering_data> obs_vector;
	for (const i3& h : hkl) {
		scattering_data temp(obs_by_hkl.at(h));
		obs.abs_F_obs.push_back(temp.abs_F_obs);
		obs.F_obs.push_back(temp.F_obs);
		obs.F_obs2.push_back(temp.F_obs2);
		obs.sigma_obs.push_back(temp.sigma_obs);
		obs.sigma_obs2.push_back(temp.sigma_obs2);
	}
	file << " done!\nNr of reflections read from file: " << hkl.size() << std::endl;

	//The fit stays on the measured reflections (hkl, in step with obs); the twin images R h only join the set
	//the written table covers, so a twinned refinement finds every twin mate it asks for
	obs.hkl_table = hkl;
	for (const i3& hkl__ : hkl)
		for (int i = 0; i < twin_law.size(); i++)
			obs.hkl_table.emplace(i3{
				int(std::lround(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2])),
				int(std::lround(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2])),
				int(std::lround(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2])) });
	if (twin_law.size() > 0)
		file << "Reflections in the table with " << twin_law.size() << " twin law(s): " << obs.hkl_table.size() << std::endl;

	ivec3 sym(3);
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
				file << "\n";
			}
			file << "\n";
		}
	}
	else
		file << "Number of symmetry operations: " << sym[0][0].size() << std::endl;

	i3 tempv;
	hkl_list hkl_enlarged;
	hkl_enlarged = obs.hkl_table;
	for (int s = 0; s < sym[0][0].size(); s++)
	{
		if (sym[0][0][s] == 1 && sym[1][1][s] == 1 && sym[2][2][s] == 1 &&
			sym[0][1][s] == 0 && sym[0][2][s] == 0 && sym[1][2][s] == 0 &&
			sym[1][0][s] == 0 && sym[2][0][s] == 0 && sym[2][1][s] == 0)
		{
			continue;
		}
		for (const i3& hkl__ : obs.hkl_table)
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

	file << "Nr of reflections to be used: " << hkl.size() << std::endl;
	return hkl_enlarged;
}

// Basically a clone of the read_hkl_full function but for reading masks, can be collapsed into a template function in the future
void read_mask(const std::filesystem::path& mask_filename,
	const hkl_list& hkl,
	const hkl_list& hkl_enlarged,
	const cell& unit_cell,
	std::ostream& file,
	vec2& mask_info,
	const int& n_r,
	bool debug
)
{
	// Mask values of the fab file at their position in hkl_enlarged
	cvec mask_enlarged(hkl_enlarged.size(), 0.0);
	bvec in_fab(hkl_enlarged.size(), false);
	const std::vector<i3> enlarged(hkl_enlarged.begin(), hkl_enlarged.end());
	mask_info.resize(n_r, vec(2, 0.0));
	i3 hkl_;
	double real_mask, imag_mask;
	err_checkf(std::filesystem::exists(mask_filename), "fab file does not exists!", file);
	std::ifstream mask_input(mask_filename, std::ios::in);
	mask_input.seekg(0, mask_input.beg);
	std::string line, temp;
	//hkl is a set, so it hands the reflections back in (h,k,l) order whatever the file order;
	//the mask has to follow that order, or F_calc[i] meets the wrong observation
	std::map<i3, scattering_data> mask_by_hkl;
	while (!mask_input.eof())
	{
		getline_universal(mask_input, line);
		if (mask_input.eof())
			break;
		if (line.size() < 2)
			continue;
		if (hkl_text_line(line))
			continue;
		err_checkf(line.size() >= 12, "mask line too short for h k l: '" + line + "'", file);
		std::istringstream iss(line);
		for (int i = 0; i < 3; i++)
		{
			// Values are separated by whitespace, so we can use std::istringstream to parse the line
			iss >> temp;
			hkl_[i] = stoi(temp);
		}
		iss >> temp;
		real_mask = stod(temp);
		iss >> temp;
		imag_mask = stod(temp);
		//a binary search in the sorted copy: a set iterator distance is linear in hkl_enlarged
		const auto it = std::lower_bound(enlarged.begin(), enlarged.end(), hkl_);
		if (it != enlarged.end() && *it == hkl_) {
			const size_t pos = it - enlarged.begin();
			mask_enlarged[pos] = cdouble(real_mask, imag_mask);
			in_fab[pos] = true;
		}
	}
	// The fab file lists one reflection of each set of equivalents: F(h) = exp(2 pi i h.t_s) F(h R_s)
	const vec2 trans = unit_cell.get_trans();
	auto h = hkl.begin();
	const ivec2 lookups = unit_cell.generate_asym_lookup(hkl, hkl_enlarged);
	for (int r = 0; r < n_r; r++, h++) {
		const ivec& lookup = lookups[r];
		for (int s = 0; s < lookup.size(); s++) {
			if (!in_fab[lookup[s]]) continue;
			const cdouble F = std::exp(cdouble(0, constants::TWO_PI * ((*h)[0] * trans[0][s] + (*h)[1] * trans[1][s] + (*h)[2] * trans[2][s]))) * mask_enlarged[lookup[s]];
			mask_info[r][0] = F.real();
			mask_info[r][1] = F.imag();
			break;
		}
	}
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
	bool debug,
	const ivec2& hkl_min_max)
{
	using namespace std;
	const ivec3 sym = unit_cell.get_sym();
	const int n_sym = sym[0][0].size();
	//An index box keeps h when an image h.R, or its Friedel mate, lies in the box: the set smtbx table_based.h
	//resolves for a measured list. A rotation preserves d, so every image of a box index is inside the sphere.
	//sym[j][h][s] holds R^T, so the products below form h.R, the cctbx convention.
	const bool boxed = hkl_min_max.size() == 3;
	if (boxed)
		file << "Keeping only the symmetry images of the index box ["
			<< setw(2) << hkl_min_max[0][0] << "," << setw(2) << hkl_min_max[0][1] << "] ; ["
			<< setw(2) << hkl_min_max[1][0] << "," << setw(2) << hkl_min_max[1][1] << "] ; ["
			<< setw(2) << hkl_min_max[2][0] << "," << setw(2) << hkl_min_max[2][1] << "]\n";
	auto in_box = [&hkl_min_max](const int h, const int k, const int l) {
		return h >= hkl_min_max[0][0] && h <= hkl_min_max[0][1] &&
			k >= hkl_min_max[1][0] && k <= hkl_min_max[1][1] &&
			l >= hkl_min_max[2][0] && l <= hkl_min_max[2][1];
	};
	auto wanted = [&](const int h, const int k, const int l) {
		if (!boxed || in_box(h, k, l) || in_box(-h, -k, -l))
			return true;
		for (int s = 0; s < n_sym; s++)
		{
			const int h_ = h * sym[0][0][s] + k * sym[0][1][s] + l * sym[0][2][s];
			const int k_ = h * sym[1][0][s] + k * sym[1][1][s] + l * sym[1][2][s];
			const int l_ = h * sym[2][0][s] + k * sym[2][1][s] + l * sym[2][2][s];
			if (in_box(h_, k_, l_) || in_box(-h_, -k_, -l_))
				return true;
		}
		return false;
	};
	file << "Generating hkl indices up to d=: " << fixed << setw(17) << setprecision(2) << dmin << flush;
	//The sphere d*^2 <= 1/d_keep^2 (cctbx index_generator) is closed under the point group, so its Friedel half
	//l > 0 | l = 0, k > 0 | k = l = 0, h > 0 is the list; d_keep sits 1e-3 inside dmin to keep a reflection at dmin from a rounded cell.
	const double d_keep = dmin * (1.0 - 1e-3);
	const double s_max = 1.0 / (d_keep * d_keep);
	const array<double, 6> G = unit_cell.get_reciprocal_metric();
	const ivec extreme = {
		int(unit_cell.get_a() / d_keep + 1e-4),
		int(unit_cell.get_b() / d_keep + 1e-4),
		int(unit_cell.get_c() / d_keep + 1e-4) };
	if (debug)
		file << "extreme: " << extreme[0] << " " << extreme[1] << " " << extreme[2] << endl;
	for (int h = -extreme[0]; h <= extreme[0]; h++)
	{
		for (int k = -extreme[1]; k <= extreme[1]; k++)
		{
			//d*^2 = A l^2 + B l + C; the loops run in the order i3 compares, hence the hint
			const double A = G[2];
			const double B = 2.0 * (h * G[4] + k * G[5]);
			const double C = h * h * G[0] + k * k * G[1] + 2.0 * h * k * G[3] - s_max;
			const double disc = B * B - 4.0 * A * C;
			if (disc < 0.0)
				continue;
			const double root = sqrt(disc);
			const int l_lo = max(int(ceil((-B - root) / (2.0 * A))), (k > 0 || (k == 0 && h > 0)) ? 0 : 1);
			const int l_hi = min(int(floor((-B + root) / (2.0 * A))), extreme[2]);
			for (int l = l_lo; l <= l_hi; l++)
				if (wanted(h, k, l))
					hkl.emplace_hint(hkl.end(), i3{ h, k, l });
		}
	}
	file << "... done!\nNr of reflections generated: " << setw(21) << hkl.size() << endl;
	if (debug)
	{
		file << "Read " << sym[0][0].size() << " symmetry elements!" << endl;
		for (int i = 0; i < sym[0][0].size(); i++)
		{
			for (int x = 0; x < 3; x++)
			{
				for (int y = 0; y < 3; y++)
					file << setw(3) << sym[y][x][i];
				file << "\n";
			}
			file << "\n";
		}
	}
	else
		file << "Number of symmetry operations: " << setw(19) << sym[0][0].size() << endl;
	//A twin law is no point-group operation, so its images leave the sphere: expand and reduce
	if (twin_law.size() > 0)
	{
		if (debug)
			file << "Number of reflections before twin: " << hkl.size() << endl;
		hkl_list twinned = hkl;
		for (const i3& hkl__ : hkl)
			for (int i = 0; i < twin_law.size(); i++)
				twinned.emplace(i3{
					int(std::lround(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2])),
					int(std::lround(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2])),
					int(std::lround(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2])) });
		if (debug)
			file << "Number of reflections after twin: " << twinned.size() << endl;
		i3 tempv;
		hkl_list hkl_enlarged = twinned;
		for (int s = 0; s < sym[0][0].size(); s++)
		{
			if (sym[0][0][s] == 1 && sym[1][1][s] == 1 && sym[2][2][s] == 1 &&
				sym[0][1][s] == 0 && sym[0][2][s] == 0 && sym[1][2][s] == 0 &&
				sym[1][0][s] == 0 && sym[2][0][s] == 0 && sym[2][1][s] == 0)
			{
				continue;
			}
			for (const i3& hkl__ : twinned)
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
				hkl.emplace(hkl__);
		}
		hkl.erase(i3{ 0, 0, 0 });
	}
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
	gap there. Olex2 sends -dmin for ED, which generates the sphere; this
	serves callers that pass only a box.
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
					int(std::lround(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2])),
					int(std::lround(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2])),
					int(std::lround(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2])) });
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
				file << "\n";
			}
			file << "\n";
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

void generate_hkl_from_options(const options& opt,
	hkl_list& hkl,
	cell& unit_cell,
	std::ostream& file)
{
	const bool have_dmin = opt.dmin != 99.0;
	const bool have_box = opt.hkl_min_max[0][0] != -100 && opt.hkl_min_max[2][1] != 100;
	if (opt.hkl_scale != 1)
	{
		//scaled indices overflow the 3I4 columns of read_hkl; the list is complete, so no images are added
		err_checkf(opt.hkl != "", "-hkl_scale needs the index list in -hkl", file);
		file << "Reading: " << std::setw(44) << opt.hkl << std::flush;
		std::ifstream in(opt.hkl);
		std::string line;
		i3 h;
		while (getline_universal(in, line))
		{
			std::istringstream row(line);
			if (row >> h[0] >> h[1] >> h[2])
				hkl.emplace(h);
		}
		hkl.erase(i3{ 0, 0, 0 });
		err_checkf(!hkl.empty(), "No reflections read from " + opt.hkl.string(), file);
		file << " done!\nNr of reflections read in units of 1/" << opt.hkl_scale << ": " << hkl.size() << std::endl;
		return;
	}
	if (have_dmin && opt.electron_diffraction)
		generate_hkl(opt.dmin / 2.0 - 0.001, hkl, opt.twin_law, unit_cell, file, opt.debug);
	else if (have_dmin && have_box)
		generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug, opt.hkl_min_max);
	else if (have_dmin)
		generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug);
	else if (have_box)
		generate_hkl(opt.hkl_min_max, hkl, opt.twin_law, unit_cell, file, opt.debug, opt.electron_diffraction);
	else
		read_hkl(opt.hkl, hkl, opt.twin_law, unit_cell, file, opt.debug);
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
	const int lim_k = static_cast<int>(extreme / stepsize[1]);
#pragma omp parallel for private(hkl_)
	for (int h = -lim; h < lim; h++)
	{
		const double _h = h * stepsize[0];
		//k as index * step like h and l: a running sum drifts off the grid the symmetry images land on
		for (int k = -lim_k; k < lim_k; k++)
		{
			// only need 0 to extreme, since we have no DISP signal
			for (int l = 0; l < lim; l++)
			{
				hkl_ = { _h, k * stepsize[1], l * stepsize[2] };
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
				file << "\n";
			}
			file << "\n";
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
			//h - k of two grid values misses the grid value by an ulp and std::set keys on exact doubles
			for (int j = 0; j < 3; j++)
			{
				const double n = std::round(tempv[j] / stepsize[j]);
				if (std::abs(tempv[j] / stepsize[j] - n) < 1e-6)
					tempv[j] = n * stepsize[j];
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
