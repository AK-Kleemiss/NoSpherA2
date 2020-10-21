/*
 * structrue_factors.cpp
 *
 *  Created on: May 27, 2019
 *      Author: florian
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <list>
#include <numeric>
#include <iomanip>
#include <string>
#include <fstream>
#include <regex>
#include <omp.h>
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <omp.h>
#else
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>
#endif
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <vector>
#include <complex>
#include <algorithm>

#include "cell.h"
#include "convenience.h"
#include "structure_factors.h"
#include "wfn_class.h"
#include "spherical_density.h"
using namespace std;
# define M_PI           3.14159265358979323846  /* pi */

#include "AtomGrid.h"
/*
std::vector<std::size_t> sort_permutation(
	const std::vector<double>& vec,
	const std::vector<double>& vec1, 
	const std::vector<double>& vec2)
{
	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j) { return pow(vec[i],2) + pow(vec1[i],2) + pow(vec2[i],2) > pow(vec[j], 2) + pow(vec1[j], 2) + pow(vec2[j], 2); });
	return p;
}

std::vector<std::size_t> sort_permutation(
	const std::vector<double>& vec)
{
	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j) { return vec[i] > vec[j]; });
	return p;
}

template <typename T>
void apply_permutation_in_place(
	std::vector<T>& vec,
	const std::vector<std::size_t>& p)
{
	std::vector<bool> done(vec.size());
	for (std::size_t i = 0; i < vec.size(); ++i)
	{
		if (done[i])
		{
			continue;
		}
		done[i] = true;
		std::size_t prev_j = i;
		std::size_t j = p[i];
		while (i != j)
		{
			std::swap(vec[prev_j], vec[j]);
			done[j] = true;
			prev_j = j;
			j = p[j];
		}
	}
}*/

double compute_dens(
	WFN &wave,
	double * PosGrid,			// [3] array with current position on the grid
	double ** mo_coef,
	int atom = -1
	)
{
	const int nmo = wave.get_nmo(true);
	double* phi = new double[nmo];
	double Rho=0.0;

	for (int i = 0; i < nmo; i++)
		phi[i] = 0.0;

	int iat;
	double d[3];
	int l[3];
	double ex;
	int mo;
	double temp;

	for (int j = 0; j < wave.get_nex(); j++) {
		iat = wave.get_center(j) - 1;
//		if (iat != atom) continue;
		type2vector(wave.get_type(j), l);

		d[0] = PosGrid[0] - wave.atoms[iat].x;
		d[1] = PosGrid[1] - wave.atoms[iat].y;
		d[2] = PosGrid[2] - wave.atoms[iat].z;
		temp = -wave.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < -46.0517) //corresponds to cutoff of ex ~< 1E-20
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0)		continue;
			else if (l[k] == 1)	ex *= d[k];
			else if (l[k] == 2)	ex *= d[k] * d[k];
			else if (l[k] == 3)	ex *= pow(d[k], 3);
			else if (l[k] == 4)	ex *= pow(d[k], 4);
			else if (l[k] == 5)	ex *= pow(d[k], 5);
		}
		for (mo = 0; mo < nmo; mo++)
			phi[mo] += mo_coef[mo][j] * ex;      //build MO values at this point
	}

	for (mo = 0; mo < nmo; mo++)
		Rho += wave.get_MO_occ(mo) * pow(phi[mo], 2);
	delete[](phi);
	return Rho;
}

double linear_interpolate_spherical_density(
	vector <double> & radial_dens,
	vector <double> & spherical_dist,
	const double dist,
	const double lincr,
	const double start
	)
{
	double result = 0;
	if (dist > spherical_dist[spherical_dist.size() - 1])
		return 0;
	else if (dist < spherical_dist[0])
		return radial_dens[0];
	int nr = floor(log(dist / start) / lincr);
	//int end = spherical_dist.size() - 1;
	//for (int i = 0; i < end; i++) {
	//	double prev_dist = spherical_dist[i];
	//	double current_dist = spherical_dist[i + 1];
	//	if (prev_dist < dist && current_dist > dist) {
	//		if (i == 0)
	//			result = radial_dens[0];
	//		else
				result = radial_dens[nr] + (radial_dens[nr + 1] - radial_dens[nr]) / (spherical_dist[nr] - spherical_dist[nr - 1]) * (dist - spherical_dist[nr - 1]);
	//	}
	//}
	if (result < 1E-10) result =  0;
	return result;
}
/*double logarithmic_interpolate_spherical_density(
	vector <double> radial_dens,
	double dist,
	int size,
	double incr
)
{
	double result = 0;
	for (int i = 1; i < size; i++) {
		double prev_dist = 0.000001 * pow(incr, i - 1);
		double current_dist = prev_dist * incr;
		if (prev_dist < dist && current_dist > dist)
			i == 0 ? result = radial_dens[i] : result = radial_dens[i - 1] * exp((dist - prev_dist) * (log(radial_dens[i]) - log(radial_dens[i - 1])) / (current_dist - prev_dist));
	}
	if (result < 1E-15) result = 0;
	return result;
}*/

bool calculate_structure_factors_HF(
	string& hkl_filename,
	string& cif,
	string& asym_cif,
	string& symm,
	WFN& wave,
	bool debug,
	int accuracy,
	ofstream& file,
	vector <int> &input_groups,
	int cpus,
	bool electron_diffraction,
	int pbc
	)
{
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
	}
	int* atom_z = new int[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
	double *x,*y,*z;
	x = new double[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
	y = new double[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
	z = new double[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
	double* alpha_max = new double[wave.get_ncen()];
	int* max_l = new int[wave.get_ncen()];
	int max_l_overall = 0;

#ifdef _WIN64
	time_t start = time(NULL);
	time_t end_becke, end_prototypes, end_spherical, end_prune, end_aspherical, end_tsc;
#else
	struct timeval t1, t2;

	gettimeofday(&t1, 0);
#endif

	if (debug)
		file << "Reading hkl now" << endl;

	vector< vector <int> > hkl;
	string line;
	hkl.resize(3);
	if (!exists(hkl_filename)) {
		file << "HKL file does not exists!" << endl;
		return false;
	}
	ifstream hkl_input(hkl_filename.c_str(), ios::in);
	hkl_input.seekg(0, hkl_input.beg);
	regex r{ R"([abcdefghijklmnopqrstuvwxyz\(\)ABCDEFGHIJKLMNOPQRSTUVW])" };
	while (!hkl_input.eof()) {
		getline(hkl_input, line);
		if (hkl_input.eof())
			break;
		if (line.size() < 2)
			continue;
		cmatch result;
		if (regex_search(line.c_str(), result, r))
			continue;
		//if (debug) file << "hkl: ";
		for (int i = 0; i < 3; i++) {
			string temp = line.substr(4 * size_t(i) + 1, 3);
			temp.erase(remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
			hkl[i].push_back(stoi(temp));
			//if (debug) file << setw(4) << temp;
		}
		//if (debug) file << endl;
		if (hkl[0][hkl[0].size() - 1] == 0 && hkl[1][hkl[0].size() - 1] == 0 && hkl[2][hkl[0].size() - 1] == 0) {
			if (debug) file << "popping back 0 0 0" << endl;
			for (int i = 0; i < 3; i++)
				hkl[i].pop_back();
		}
	}
	hkl_input.close();
	// Remove duplicate reflections
	for (int i = 0; i < hkl[0].size(); i++)
		for (int j = i + 1; j < hkl[0].size(); j++)
			if (hkl[0][i] == hkl[0][j] && hkl[1][i] == hkl[1][j] && hkl[2][i] == hkl[2][j])
				for (int x = 0; x < 3; x++)
					hkl[x].erase(hkl[x].begin() + j);

	if (debug)
		file << "Reflections read! Nr of reflections: " << hkl[0].size() << endl;

	if (debug)
		file << "starting to read cif!" << endl;
	if (!exists(cif)) {
		file << "CIF does not exists!" << endl;
		return false;
	}
	if (!exists(asym_cif)) {
		file << "Asym CIF does not exists!" << endl
			<< asym_cif << endl;
		return false;
	}

	cell unit_cell(cif, file, debug);

	if (debug) {
		file << "RCM done, now labels and asym atoms!" << endl;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j)
				file << setw(10) << fixed << unit_cell.get_rcm(i,j) / 2 / M_PI / 0.529177249 << ' ';
			file << endl;
		}
		file << "CM in bohr:" << endl;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j)
				file << setw(10) << fixed << unit_cell.get_cm(i, j) << ' ';
			file << endl;
		}
	}

	ifstream cif_input(cif.c_str(), std::ios::in);
	ifstream asym_cif_input(asym_cif.c_str(), std::ios::in);
	cif_input.clear();
	cif_input.seekg(0, cif_input.beg);
	vector <string> labels;
	vector <int> atom_type_list;
	vector <int> asym_atom_to_type_list;
	int count_fields = 0;
	int group_field = 0;
	int position_field[3] = { 0,0,0 };
	int label_field = 1000;
	vector <int> asym_atom_list;
	vector <int> all_atom_list;
	vector < bool > is_asym;
	vector < vector <double > > positions;
	positions.resize(wave.get_ncen());
	is_asym.resize(wave.get_ncen());
	for (int i = 0; i < wave.get_ncen(); i++) {
		is_asym[i] = false;
		positions[i].resize(3);
	}
	bool atoms_read = false;
	while (!asym_cif_input.eof() && !atoms_read) {
		getline(asym_cif_input, line);
		//if(debug) file << "line: "<< line << endl;
		if (line.find("loop_") != string::npos) {
			//if(debug) file << "found loop!" << endl;
			while (line.find("_") != string::npos) {
				getline(asym_cif_input, line);
				if (debug) file << "line in loop field definition: " << line << endl;
				if (line.find("label") != string::npos)
					label_field = count_fields;
				else if (line.find("fract_x") != string::npos)
					position_field[0] = count_fields;
				else if (line.find("fract_y") != string::npos)
					position_field[1] = count_fields;
				else if (line.find("fract_z") != string::npos)
					position_field[2] = count_fields;
				else if (label_field == 1000) {
					if (debug) file << "I don't think this is the atom block.. moving on!" << endl;
					break;
				}
				count_fields++;
			}
			while (line.find("_") == string::npos && line.length() > 3) {
				//if(debug) file << "Reading atom!"<< endl;
				atoms_read = true;
				stringstream s(line);
				vector <string> fields;
				fields.resize(count_fields);
				for (int i = 0; i < count_fields; i++)
					s >> fields[i];
				if (debug) file << "label: " << fields[label_field] << " frac_position: " << stod(fields[position_field[0]]) << " " << stod(fields[position_field[1]]) << " " << stod(fields[position_field[2]]) << endl;
				positions[labels.size()] = unit_cell.get_coords_cartesian(stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]));
				bool found_this_one = false;
				if (debug) file << "label: " << fields[label_field] << " position: " << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
				for (int i = 0; i < wave.get_ncen(); i++) {
					if (is_similar(positions[labels.size()][0], wave.atoms[i].x, -1)
						&& is_similar(positions[labels.size()][1], wave.atoms[i].y, -1)
						&& is_similar(positions[labels.size()][2], wave.atoms[i].z, -1)) {
						if (debug) file << "WFN position: " << wave.atoms[i].x << " " << wave.atoms[i].y << " " << wave.atoms[i].z << endl
							<< "Found an atom: " << fields[label_field] << " Corresponding to atom charge " << wave.atoms[i].charge << endl;
						wave.atoms[i].label = fields[label_field];
						all_atom_list.push_back(i);
						found_this_one = true;
						break;
					}
				}
				if (!found_this_one && debug)
					file << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
				labels.push_back(fields[label_field]);
				getline(asym_cif_input, line);
			}
		}
	}
	if (labels.size() != wave.get_ncen()) {
		file << "Number of atoms in labels: " << labels.size() << " and number of atoms in Wavefunction: " << wave.get_ncen() << "!" << endl << "This is not good, i will stop here!" << endl;
		return false;
	}
	atoms_read = false;
	label_field = 1000;
	count_fields = 0;
	while (!cif_input.eof() && !atoms_read) {
		getline(cif_input, line);
		//if(debug) file << "line: "<< line << endl;
		if (line.find("loop_") != string::npos) {
			//if(debug) file << "found loop!" << endl;
			while (line.find("_") != string::npos) {
				getline(cif_input, line);
				if (debug) file << "line in loop field definition: " << line << endl;
				if (line.find("label") != string::npos)
					label_field = count_fields;
				else if (line.find("disorder_group") != string::npos)
					group_field = count_fields;
				else if (label_field == 1000) {
					if (debug) file << "I don't think this is the atom block.. moving on!" << endl;
					break;
				}
				count_fields++;
			}
			while (line.find("_") == string::npos && line.length() > 3) {
				atoms_read = true;
				stringstream s(line);
				vector <string> fields;
				fields.resize(count_fields);
				for (int i = 0; i < count_fields; i++)
					s >> fields[i];
				fields[label_field].erase(remove_if(fields[label_field].begin(), fields[label_field].end(), ::isspace), fields[label_field].end());
				for (int atom = 0; atom < wave.get_ncen(); atom++) {
					if (debug) file << "Comparing atoms: " << fields[label_field] << " / " << labels[atom] << endl;
					if (fields[label_field] == labels[atom]) {
						int nr = -1;
						for (int i = 0; i < wave.get_ncen(); i++) {
							if (is_similar(positions[atom][0], wave.atoms[i].x, -1)
								&& is_similar(positions[atom][1], wave.atoms[i].y, -1)
								&& is_similar(positions[atom][2], wave.atoms[i].z, -1)) {
								if (debug) {
									file << "Found an asymmetric atom: " << fields[label_field] << " Corresponding to atom charge " << wave.atoms[i].charge << endl
										<< "checking disorder group: " << fields[group_field] << " vs. ";
									for (int g = 0; g < input_groups.size(); g++)
										file << input_groups[g] << "+";
									file << endl;
								}
								bool yep = false;
								for (int g = 0; g < input_groups.size(); g++)
									if (stod(fields[group_field]) == input_groups[g])
										yep = true;
								if (!yep && input_groups.size()>0) continue;
								nr = i;
								asym_atom_list.push_back(i);
								break;
							}
						}
						is_asym[atom] = true;
						bool already_there = false;
						for (int i = 0; i < atom_type_list.size(); i++)
							if (atom_type_list[i] == wave.atoms[nr].charge) {
								already_there = true;
								asym_atom_to_type_list.push_back(i);
							}
						if (already_there == false) {
							asym_atom_to_type_list.push_back(atom_type_list.size());
							atom_type_list.push_back(wave.atoms[nr].charge);
						}
						break;
					}
				}
				getline(cif_input, line);
			}
		}
	}

	if (debug) {
		file << "There are " << atom_type_list.size() << " types of atoms" << endl;
		for (int i = 0; i < atom_type_list.size(); i++)
			file << setw(4) << atom_type_list[i];
		file << endl << "asym_atoms_to_type_list: " << endl;
		for (int i = 0; i < asym_atom_to_type_list.size(); i++)
			file << setw(4) << asym_atom_to_type_list[i];
		file << endl;
		file << "Mapping of asym atoms:" << endl;
		for (int i = 0; i < wave.get_ncen(); i++)
			file << setw(4) << wave.atoms[all_atom_list[i]].charge;
		file << endl;
		for (int i = 0; i < wave.get_ncen(); i++)
			file << setw(4) << is_asym[i];
		file << endl;
	}

	bool symm_read = false;
	//Still need to read the sym matrices
	if (symm == "") {
		if (debug) file << "No Symmetry file specified, tyring to read from the CIF!" << endl;
		symm_read = true;
	}
	else if (!exists(symm) && symm_read == true)
		return false;

	vector < vector < vector <int> > > sym;
	sym.resize(3);
	for (int i = 0; i < 3; i++)
		sym[i].resize(3);
	if (!symm_read) {
		ifstream symm_input(symm.c_str(), ios::in);
		string liny;
		int temp_int;
		while (!symm_input.eof()) {
			getline(symm_input, liny);
			stringstream streamy(liny);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) {
					streamy >> temp_int;
					sym[i][j].push_back(temp_int);
				}
		}
	}
	else {
		sym = unit_cell.get_sym();
	}

	if (debug) {
		file << "Read " << sym[0][0].size() << " symmetry elements! Size of sym: " << sym[0][0].size() << endl;
		for (int i = 0; i < sym[0][0].size(); i++) {
			for (int x = 0; x < 3; x++) {
				for (int y = 0; y < 3; y++)
					file << setw(3) << sym[y][x][i];
				file << endl;
			}
			file << endl;
		}
	}

	cif_input.close();

	if (debug)
		file << "There are " << atom_type_list.size() << " Types of atoms and " << asym_atom_to_type_list.size() << " atoms in total" << endl;

	if (asym_atom_list.size() == 0) {
		file << "0 asym atoms is imposible! something is wrong with reading the CIF!" << endl;
		return false;
	}

	if (debug)
		file << "made it post CIF, now make grids!" << endl;

	double*** grid = new double** [asym_atom_list.size()];
	int* num_points = new int[asym_atom_list.size()];
	for (int i = 0; i < asym_atom_list.size(); i++)
		grid[i] = new double* [6];
	// GRID COORDINATES for [a][c][p] a = atom [0,ncen],
	// c = coordinate [0=x, 1=y, 2=z, 3=atomic becke weight and molecular becke weight, 5=total spherical density],
	// p = point in this grid

#pragma omp parallel for
	for (int i = 0; i < wave.get_ncen(); i++) {
		atom_z[i] = wave.atoms[i].charge;
		x[i] = wave.atoms[i].x;
		y[i] = wave.atoms[i].y;
		z[i] = wave.atoms[i].z;
        //if(debug)
        //    file << "xyz= 000 position: " << x[i] << " " << y[i] << " " << z[i] << " Charge: " << atom_z[i] << endl;
		if (pbc != 0) {
			int j = 0;
			for (int pbc_x = -pbc; pbc_x < pbc + 1; pbc_x++)
				for (int pbc_y = -pbc; pbc_y < pbc + 1; pbc_y++)
					for (int pbc_z = -pbc; pbc_z < pbc + 1; pbc_z++) {
						if (pbc_x == 0 && pbc_y == 0 && pbc_z == 0)
                            continue;
                        else{
						    j++;
    						atom_z[i + j * wave.get_ncen()] = wave.atoms[i].charge;
    						x[i + j * wave.get_ncen()] = wave.atoms[i].x + pbc_x * unit_cell.get_cm(0, 0) + pbc_y * unit_cell.get_cm(0, 1) + pbc_z * unit_cell.get_cm(0, 2);
							y[i + j * wave.get_ncen()] = wave.atoms[i].y + pbc_x * unit_cell.get_cm(1, 0) + pbc_y * unit_cell.get_cm(1, 1) + pbc_z * unit_cell.get_cm(1, 2);
							z[i + j * wave.get_ncen()] = wave.atoms[i].z + pbc_x * unit_cell.get_cm(2, 0) + pbc_y * unit_cell.get_cm(2, 1) + pbc_z * unit_cell.get_cm(2, 2);
                            if(debug) 
                                file << "xyz= " << pbc_x << pbc_y << pbc_z << " j = " << j << " position: " << x[i + j * wave.get_ncen()] << " " << y[i + j * wave.get_ncen()] << " " << z[i + j * wave.get_ncen()] << " Charge: " << atom_z[i + j * wave.get_ncen()] << endl;
                        }
					}
		}
		alpha_max[i] = 0.0;
		max_l[i] = 0;
		for (int b = 0; b < wave.get_nex(); b++) {
			if (wave.get_center(b) != i + 1)
				continue;
			if (wave.get_exponent(b) > alpha_max[i])
				alpha_max[i] = wave.get_exponent(b);
			if (wave.get_type(b) > max_l[i]) {
				int l = wave.get_type(b);
				if (l == 1)
					l = 1;
				else if (l >= 2 && l <= 4)
					l = 2;
				else if (l >= 5 && l <= 10)
					l = 3;
				else if (l >= 11 && l <= 20)
					l = 4;
				else if (l >= 21 && l <= 35)
					l = 5;
				max_l[i] = l;
				if (l > max_l_overall)
					max_l_overall = l;
			}
		}
	}

	if (debug)
		file << "Atoms are there!" << endl;

	double** alpha_min = new double* [wave.get_ncen()];
	for (int i = 0; i < wave.get_ncen(); i++)
		alpha_min[i] = new double[max_l_overall];

#pragma omp parallel for
	for (int i = 0; i < wave.get_ncen(); i++) {
		for (int b = 0; b < max_l_overall; b++)
			alpha_min[i][b] = 100000000.0;
	}

#pragma omp parallel for
	for (int i = 0; i < wave.get_ncen(); i++) {
		for (int b = 0; b < wave.get_nex(); b++) {
			if (wave.get_center(b) != i + 1)
				continue;
			int l = wave.get_type(b);
			if (l == 1)
				l = 1;
			else if (l >= 2 && l <= 4)
				l = 2;
			else if (l >= 5 && l <= 10)
				l = 3;
			else if (l >= 11 && l <= 20)
				l = 4;
			else if (l >= 21 && l <= 35)
				l = 5;
			else if (l >= 36 && l <= 56)
				l = 6;
			if (wave.get_exponent(b) < alpha_min[i][l - 1])
				alpha_min[i][l - 1] = wave.get_exponent(b);
		}
		if (debug) {
			file << "alpha_min: ";
			for (int b = 0; b < max_l_overall; b++)
				file << setw(14) << scientific << alpha_min[i][b];
			file << endl;
		}
	}

	if (debug)
		file << "alpha_min is there!" << endl
		<< "Nr of asym atoms: " << asym_atom_list.size() << " Number of all atoms: " << all_atom_list.size() << " Number of atoms in wfn: " << wave.get_ncen() << endl;
	else
		file << "There are:\n" << setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
		//<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
		<< setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << endl;
	// Total grid as a sum of all atomic grids.
	// Dimensions: [c] [p]
	// p = the number of gridpoint
	// c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
	vector < vector <  double > > total_grid;
	// density of spherical atom at each
	// Dimensions: [a] [d]
	// a = atom number in atom type list for which the weight is calcualted
	// d = distance to look at obtained from point_to_distance_map
	vector < vector < double > > spherical_density;
	

	if (debug)
		file << "Making Becke Grid for" << endl;
	else
		file << endl << "Making Becke Grids..." << flush;

	total_grid.resize(7);

	//Make Prototype grids with only single atom weights for all elements
	vector <AtomGrid> Prototype_grids;
	if (debug) file << "max_l_overall: " << max_l_overall << endl;
	for (int i = 0; i < atom_type_list.size(); i++) {
		if (debug) file << "Atom " << i << endl;
		double alpha_max_temp;
		double max_l_temp;
		double* alpha_min_temp = new double[max_l_overall];
		for (int j = 0; j < wave.get_ncen(); j++) {
			if (wave.atoms[j].charge == atom_type_list[i]) {
				if (debug) {
					file << alpha_max[j] << " " << max_l[j] - 1 << " ";
					for (int l = 0; l < max_l_overall; l++)
						file << alpha_min[j][l] << " ";
					file << endl;
				}
				alpha_max_temp = alpha_max[j];
				max_l_temp = max_l[j] - 1;
				for (int l = 0; l <= max_l_temp; l++)
					alpha_min_temp[l] = alpha_min[j][l];
				break;
			}
		}

		if (debug) {
			file << "max_l: " << defaultfloat << max_l_temp << " alpha_max: " << scientific << alpha_max_temp << " alpha_min: ";
			for (int l = 0; l <= max_l_temp; l++)
				file << setw(14) << scientific << alpha_min_temp[l];
			file << " accuracy: " << accuracy << endl;
		}
		int lebedev_high, lebedev_low;
		double radial_acc;
		if (accuracy == 0) {
			lebedev_high = (max_l_temp < 3) ? 0 : 10;
			lebedev_low = (max_l_temp < 3) ? 0 : 10;
			radial_acc = 1e-5;
		}
		else if (accuracy == 1) {
			lebedev_high = (max_l_temp < 3) ? 110 : 146;
			lebedev_low = (max_l_temp < 3) ? 38 : 50;
			radial_acc = 1e-8;
		}
		else if (accuracy == 2) {
			lebedev_high = (max_l_temp < 3) ? 230 : 266;
			lebedev_low = (max_l_temp < 3) ? 110 : 146;
			radial_acc = 1e-10;
		}
		else if (accuracy == 3) {
			lebedev_high = (max_l_temp < 3) ? 350 : 590;
			lebedev_low = (max_l_temp < 3) ? 266 : 350;
			radial_acc = 1e-15;
		}
		else if (accuracy > 3) {
			lebedev_high = (max_l_temp < 3) ? 590 : 5810;
			lebedev_low = (max_l_temp < 3) ? 350 : 4802;
			radial_acc = 1e-20;
		}
		Prototype_grids.push_back(AtomGrid(radial_acc,
			lebedev_low,
			lebedev_high,
			atom_type_list[i],
			alpha_max_temp,
			max_l_temp,
			alpha_min_temp,
			debug,
			file));

	}

#ifdef _WIN64
	end_prototypes = time(NULL);
	if (debug) {
		
		for (int prototype = 0; prototype < Prototype_grids.size(); prototype++)
			file << "Number of gridpoints for atom type " << atom_type_list[prototype] << " :" << Prototype_grids[prototype].get_num_grid_points() << endl;

		//	int diff = end - start;
		if (end_prototypes - start < 1) file << "Time until prototypes are done: <1 s" << endl;
		else if (end_prototypes - start < 60) file << "Time until prototypes are done: " << fixed << setprecision(0) << end_prototypes - start << " s" << endl;
		else if (end_prototypes - start < 3600) file << "Time until prototypes are done: " << fixed << setprecision(0) << floor((end_prototypes - start) / 60) << " m " << (end_prototypes - start) % 60 << " s" << endl;
		else file << "Time until prototypes are done: " << fixed << setprecision(0) << floor((end_prototypes - start) / 3600) << " h " << ((end_prototypes - start) % 3600) / 60 << " m" << endl;
	}
#else
	if (debug) {
		gettimeofday(&t2, 0);

		double time3 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (time3 < 60) printf("Time to prepare: %4.1lf s\n", time3);
		else if (time3 < 3600) printf("Time to prepare: %10.1lf m\n", time3 / 60);
		else printf("Time to prepare: %10.1lf h\n", time3 / 3600);
	}
#endif

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < asym_atom_list.size(); i++) {
		int nr = asym_atom_list[i];
		int type;
		for (int i = 0; i < atom_type_list.size(); i++)
			if (atom_type_list[i] == wave.atoms[nr].charge) 
				type=i;

		num_points[i] = Prototype_grids[type].get_num_grid_points();

		if (debug) file << "Number of points for atom " << i << ": " << num_points[i] << endl;

		for (int n = 0; n < 6; n++)
			grid[i][n] = new double[num_points[i]];
		for (int p = 0; p < num_points[i]; p++)
			grid[i][4][p] = 0.0;

		Prototype_grids[type].get_grid(wave.get_ncen()*pow(pbc*2+1,3),
			nr,
			x,
			y,
			z,
			atom_z,
			grid[i][0],
			grid[i][1],
			grid[i][2],
			grid[i][3],
			grid[i][5]);
	}

	Prototype_grids.clear();
	int points = 0;
	for (int i = 0; i < asym_atom_list.size(); i++)
		points += num_points[i];
	if (debug) file << "Becke Grid exists" << endl;
	else file << "                                 done! Number of gridpoints: " << defaultfloat << points << endl;

#ifdef _WIN64
	if (debug) {
		file << "Taking time..." << endl;
	}
	end_becke = time(NULL);

#else
	if (debug) {
		gettimeofday(&t2, 0);

		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
	}
#endif

	file << "Calculating spherical densities..." << flush;

	vector < vector < double > > radial_density;
	vector < vector < double > > radial_dist;

	radial_density.resize(atom_type_list.size());
	radial_dist.resize(atom_type_list.size());
	spherical_density.resize(asym_atom_list.size());
	for (int i = 0; i < asym_atom_list.size(); i++)
		spherical_density[i].resize(num_points[i]);

	double incr = pow(1.005,max(1,accuracy-1));
	const double lincr = log(incr);
	const double min_dist = 0.0000001;
	//Make radial grids
	for (int i = 0; i < atom_type_list.size(); i++) {
		if (debug) file << "Calculating for atomic number " << atom_type_list[i] << endl;
		double current = 1;
		double dist = min_dist;
		if(accuracy > 3)
			while (current > 1E-10) {
				radial_dist[i].push_back(dist);
				current = get_radial_density(atom_type_list[i], dist);
				if (current == -20)
					return false;
				radial_density[i].push_back(current);
				dist *= incr;
			}
		else
			while (current > 1E-12) {
				radial_dist[i].push_back(dist);
				current = get_radial_density(atom_type_list[i], dist);
				if (current == -20)
					return false;
				radial_density[i].push_back(current);
				dist *= incr;
			}
		if (debug)
			file << "Number of radial density points for atomic number " << atom_type_list[i] << ": " << radial_density[i].size() << endl;
	}
	//if (debug) {
	//	file << "Asym atom list: ";
	//	for (int g = 0; g < asym_atom_list.size(); g++)
	//		file << setw(4) << asym_atom_list[g];
	//	file << endl;
	//	file << "All atom list: ";
	//	for (int g = 0; g < all_atom_list.size(); g++)
	//		file << setw(4) << all_atom_list[g];
	//	file << endl;
	//	file << "spherical radial density:" << endl;
	//	for (int i = 0; i < atom_type_list.size(); i++)
	//	{
	//		file << "-------------------------------" << endl;
	//		file << "atom: " << atom_type_list[i] << endl;
	//		for (int d = 0; d < radial_density[i].size(); d++)
	//			file << "d: " << setw(12) << setprecision(6) << scientific <<  radial_dist[i][d]
	//			<< " rho: " << setw(12) << setprecision(6) << scientific << radial_density[i][d] << endl;
	//		file << endl;
	//		file << "-------------------------------" << endl;
	//	}
	//}
	//apply to the becke grid
	for (int i = 0; i < wave.get_ncen(); i++) {
		int nr = all_atom_list[i];
		int type_list_number = -1;
		for (int j = 0; j < atom_type_list.size(); j++)
			if (wave.atoms[nr].charge == atom_type_list[j])
				type_list_number = j;
		for (int g = 0; g < asym_atom_list.size(); g++) {
#pragma omp parallel for
			for (int p = 0; p < num_points[g]; p++) {
				double temp =
					linear_interpolate_spherical_density(radial_density[type_list_number]
						, radial_dist[type_list_number]
						, sqrt(pow(grid[g][0][p] - wave.atoms[nr].x, 2)
							+ pow(grid[g][1][p] - wave.atoms[nr].y, 2)
							+ pow(grid[g][2][p] - wave.atoms[nr].z, 2))
						, lincr
						, min_dist
					);
				if (all_atom_list[i] == asym_atom_list[g])
					spherical_density[g][p] = temp;
				grid[g][4][p] += temp;
			}
		}
	}
	//fill out with priodic information
	if (pbc != 0) {
		for (int x = -pbc; x < pbc + 1; x++)
			for (int y = -pbc; y < pbc + 1; y++)
				for (int z = -pbc; z < pbc + 1; z++) {
					if (x == 0 && y == 0 && z == 0)
						continue;
					for (int i = 0; i < wave.get_ncen(); i++) {
						int type_list_number = -1;
						int nr = all_atom_list[i];
						for (int j = 0; j < atom_type_list.size(); j++)
							if (wave.atoms[nr].charge == atom_type_list[j])
								type_list_number = j;
						for (int g = 0; g < asym_atom_list.size(); g++)
#pragma omp parallel for
							for (int p = 0; p < num_points[g]; p++)
								grid[g][4][p] += linear_interpolate_spherical_density(
									radial_density[type_list_number],
									radial_dist[type_list_number],
									sqrt(
										  pow(grid[g][0][p] - (wave.atoms[nr].x + x * unit_cell.get_cm(0, 0) + y * unit_cell.get_cm(0, 1) + z * unit_cell.get_cm(0, 2)), 2)
										+ pow(grid[g][1][p] - (wave.atoms[nr].y + x * unit_cell.get_cm(1, 0) + y * unit_cell.get_cm(1, 1) + z * unit_cell.get_cm(1, 2)), 2)
										+ pow(grid[g][2][p] - (wave.atoms[nr].z + x * unit_cell.get_cm(2, 0) + y * unit_cell.get_cm(2, 1) + z * unit_cell.get_cm(2, 2)), 2)
									),
									lincr,
									min_dist
								);
					}
				}
	}

	file << "                    done!" << endl;
	if (debug)
		for (int i = 0; i < asym_atom_list.size(); i++)
			file << endl << "number of points for atom " << i << " " << num_points[i] << " " << spherical_density[i].size() << endl;


	radial_density.clear();
	radial_dist.clear();
		
#ifdef _WIN64
	if (debug) {
		file << "Taking time..." << endl;
	}
	end_spherical = time(NULL);

#else
	if (debug) {
		gettimeofday(&t2, 0);

		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
	}
#endif

	file << "Pruning Grid..." << flush;
	double cutoff;
	if (accuracy < 3)
		cutoff = 1E-10;
	else if (accuracy == 3)
		cutoff = 1E-14;
	else
		cutoff = 1E-30;

	for (int i = 0; i < asym_atom_list.size(); i++) {
		int reduction = 0;
		for (int p = 0; p < num_points[i]; p++) {
			if (grid[i][4][p] != 0.0) {
				if (abs(grid[i][3][p] * spherical_density[i][p - reduction] / grid[i][4][p]) > cutoff) {
					for (int k = 0; k < 5; k++)
						total_grid[k].push_back(grid[i][k][p]);
					total_grid[6].push_back(grid[i][5][p]);
				}
				else {
					spherical_density[i].erase(spherical_density[i].begin() + (p - reduction));
					reduction++;
				}
			}
			else {
				spherical_density[i].erase(spherical_density[i].begin() + (p - reduction) );
				reduction++;
			}

		}
		num_points[i] -= reduction;
		delete[](grid[i]);
		if (debug) file << endl << "number of points for atom " << i << " " << num_points[i] << " " << spherical_density[i].size() << endl;
	}

	if (debug) file << endl;
	points = 0;
	for (int i = 0; i < asym_atom_list.size(); i++)
		points += num_points[i];
	if (debug) file << "sphericals done!" << endl;
	else file << "                                       done! Number of gridpoints: " << defaultfloat << points << endl;

#ifdef _WIN64
	if (debug) {
		file << "Taking time..." << endl;
	}
	end_prune = time(NULL);

#else
	if (debug) {
		gettimeofday(&t2, 0);

		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
	}
#endif

	file << "Calculating aspherical densities..." << flush;
    vector < vector < double > > periodic_grid;

	total_grid[5].resize(total_grid[0].size());

	double** mo_coef;
	mo_coef = new double* [wave.get_nmo()];
	for (int i = 0; i < wave.get_nmo(); i++)
		mo_coef[i] = wave.get_MO_coef_ptr(i);

#pragma omp parallel for
	for (int i = 0; i < total_grid[0].size(); i++) {
		total_grid[5][i] = compute_dens(wave, new double[3]{ 
            total_grid[0][i], 
            total_grid[1][i], 
            total_grid[2][i]},
			mo_coef);
	}
	if (pbc != 0) {
        periodic_grid.resize(pow(pbc*2+1,3));
        int j=0;
        for(int d=0; d<pow(pbc*2+1,3); d++)
            periodic_grid[d].resize(total_grid[5].size());
		for (int x = -pbc; x < pbc+1; x++)
			for (int y = -pbc; y < pbc+1; y++)
				for (int z = -pbc; z < pbc+1; z++) {
					if (x == 0 && y == 0 && z == 0)
						continue;
#pragma omp parallel for
					for (int i = 0; i < total_grid[0].size(); i++)
						periodic_grid[j][i] = compute_dens(wave, new double[3]{
							total_grid[0][i] + x * unit_cell.get_cm(0,0) + y * unit_cell.get_cm(0,1) + z * unit_cell.get_cm(0,2),
							total_grid[1][i] + x * unit_cell.get_cm(1,0) + y * unit_cell.get_cm(1,1) + z * unit_cell.get_cm(1,2),
							total_grid[2][i] + x * unit_cell.get_cm(2,0) + y * unit_cell.get_cm(2,1) + z * unit_cell.get_cm(2,2) },
							mo_coef);
                    j++;
				}
		if (debug){
			for (int i = 0; i < total_grid[0].size(); i++){
                if (i%1000==0)
                    file << "Old dens: " << total_grid[5][i] << " contributions of neighbour-cells:"; 
                for (int j=0; j<pow(pbc*2+1,3)-1; j++){
                    if (i%1000==0)
                        file << " " << periodic_grid[j][i];
                    total_grid[5][i] += periodic_grid[j][i];
                }
                if (i%1000==0)
				    file << endl;
            }
        }
	}

	for (int i = 0; i < wave.get_nmo(); i++)
		mo_coef[i] = nullptr;
	delete[] mo_coef;

	if (debug) file << endl << "with total number of points: " << total_grid[0].size() << endl;
	else file << "                   done!" << endl;

	file << "Applying hirshfeld weights and integrating charges..." << flush;
	double el_sum_becke = 0.0;
	double el_sum_spherical = 0.0;
	double el_sum_hirshfeld = 0.0;
	// Vector containing integrated numbers of electrons
	// dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
	// dimension 1: atoms of asym_atom_list
	vector < vector <double> > atom_els;
	atom_els.resize(3);
	for (int i = 0; i < asym_atom_list.size(); i++)
		for (int n = 0; n < 3; n++)
			atom_els[n].push_back(0.0);

	//Generate Electron sums
	for (int i = 0; i < asym_atom_list.size(); i++) {
		int start_p = 0;
#pragma loop(no_vector)
		for (int a = 0; a < i; a++)
			start_p += num_points[a];
		for (int p = start_p; p < start_p + num_points[i]; p++) {
			if (abs(total_grid[6][p]) > cutoff) {
				atom_els[0][i] += total_grid[6][p] * total_grid[5][p];
				atom_els[1][i] += total_grid[6][p] * total_grid[4][p];
			}
			atom_els[2][i] += total_grid[5][p] * total_grid[3][p] * spherical_density[i][p - start_p] / total_grid[4][p];
		}
		el_sum_becke += atom_els[0][i];
		el_sum_spherical += atom_els[1][i];
		el_sum_hirshfeld += atom_els[2][i];
	}

	if (debug)
		file << "Becke grid with hirshfeld weights done!" << endl;

#ifdef _WIN64
	if (debug) {
		file << "Taking time..." << endl;
	}
	end_aspherical = time(NULL);

#else
	if (debug) {
		gettimeofday(&t2, 0);

		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
	}
#endif

	if (debug) {
		file << "atom_els[2]: ";
		for (int i = 0; i < asym_atom_list.size(); i++) {
			if (isnan(atom_els[2][i]))
				file << "!!!";
			file << atom_els[2][i] << " ";
			if (isnan(atom_els[2][i]))
				file << "!!!";
		}
		file << endl;
	}

#pragma omp parallel for
	for (int p = 0; p < total_grid[0].size(); p++)
			total_grid[5][p] *= total_grid[3][p];

	file << " done!" << endl;
	file << "Number of points evaluated: " << total_grid[0].size() << " with " << el_sum_becke << " electrons in Becke Grid in total." << endl << endl;

	file << "Table of Charges in electrons" << endl << endl << "Atom       Becke   Spherical Hirshfeld" << endl;

	int counter = 0;
	for (int i = 0; i < wave.get_ncen(); i++) {
		if (is_asym[i]) {
			file << setw(6) << labels[i]
			<< fixed << setw(10) << setprecision(3) << wave.atoms[all_atom_list[i]].charge - atom_els[0][counter]
			<< fixed << setw(10) << setprecision(3) << wave.atoms[all_atom_list[i]].charge - atom_els[1][counter]
			<< fixed << setw(10) << setprecision(3) << wave.atoms[all_atom_list[i]].charge - atom_els[2][counter];
			if (debug) file << " " << setw(4) << wave.atoms[all_atom_list[i]].charge << " " << fixed << setw(10) << setprecision(3) << wave.atoms[all_atom_list[i]].charge - atom_els[0][counter]
				<< fixed << setw(10) << setprecision(3) << atom_els[1][counter]
				<< fixed << setw(10) << setprecision(3) << atom_els[2][counter];
			counter++;
			file << endl;
		}
	}

	file << "Total number of electrons in the wavefunction: " << el_sum_becke << endl << " and Hirshfeld electrons (asym unit): " << el_sum_hirshfeld << endl;

	vector < vector <double> > d1,d2,d3;
	vector < vector <double> > dens;
	points = 0;
	dens.resize(asym_atom_list.size());
	if (debug)
		file << "resized outer dens" << endl;
	d1.resize(asym_atom_list.size());
	d2.resize(asym_atom_list.size());
	d3.resize(asym_atom_list.size());
	if (debug)
		file << "resized outer d1-3" << endl;

#pragma omp parallel for reduction(+:points)
	for (int i = 0; i < asym_atom_list.size(); i++) {
		int start_p = 0;
#pragma loop(no_vector)
		for (int a = 0; a < i; a++)
			start_p += num_points[a];
		for (int p = start_p; p < start_p + num_points[i]; p++) {

			dens[i].push_back(total_grid[5][p] * spherical_density[i][p - start_p] / total_grid[4][p]);
			d1[i].push_back(total_grid[0][p] - wave.atoms[asym_atom_list[i]].x);
			d2[i].push_back(total_grid[1][p] - wave.atoms[asym_atom_list[i]].y);
			d3[i].push_back(total_grid[2][p] - wave.atoms[asym_atom_list[i]].z);

		}
		points += dens[i].size();
	}

	for (int i = 0; i < asym_atom_list.size(); i++)
		spherical_density[i].clear();
	spherical_density.clear();

	for (int grid = 0; grid < total_grid.size(); grid++)
		total_grid[grid].resize(0);
	total_grid.resize(0);

	vector < vector <double> > k_pt;
	k_pt.resize(3);
#pragma omp parallel for
	for (int i = 0; i < 3; i++)
		k_pt[i].resize(sym[0][0].size()* hkl[0].size());

	if (debug)
		file << "K_point_vector is here! size: " << k_pt[0].size() << endl;

	bool shrink = false;
	vector < vector<double> > k_pt_unique;
	vector < vector<int> > hkl_unique;
	if (shrink) {
		k_pt_unique.resize(3);
		hkl_unique.resize(3);
		for (int s = 0; s < sym[0][0].size(); s++)
			for (int ref = 0; ref < hkl[0].size(); ref++)
				for (int h = 0; h < 3; h++)
					hkl_unique[h].push_back(hkl[0][ref] * sym[h][0][s] + hkl[1][ref] * sym[h][1][s] + hkl[2][ref] * sym[h][2][s]);
	}

	for (int s = 0; s < sym[0][0].size(); s++) {
#pragma omp parallel for
		for (int ref = 0; ref < hkl[0].size(); ref++) 
			for (int x = 0; x < 3; x++)
				for (int h = 0; h < 3; h++) {
					double rcm_sym = 0.0;
					for (int j = 0; j < 3; j++)
						rcm_sym += unit_cell.get_rcm(x,j) * sym[j][h][s];
					k_pt[x][ref + s * hkl[0].size()] += rcm_sym * hkl[h][ref];
				}
	}

	
	if (shrink) {
		vector <bool> mask;
		mask.resize(k_pt[0].size());
		mask.assign(k_pt[0].size(), true);
		for (int i = 0; i < k_pt[0].size(); i++)
			for (int j = i + 1; j < k_pt[0].size(); j++) {
				if (!mask[j])
					continue;
				if (k_pt[0][i] == k_pt[0][j] && k_pt[1][i] == k_pt[1][j] && k_pt[2][i] == k_pt[2][j])
					mask[j] = false;
			}
		if (debug)  file << "Mask done!" << endl;
		for (int i = k_pt[0].size() - 1; i >= 0; i--) {
			if (debug) file << "i: " << i << " mask; " << mask[i] << endl;
			if (mask[i])
				for (int h = 0; h < 3; h++)
					k_pt_unique[h].insert(k_pt_unique[h].begin(), k_pt[h][i]);
			else
				for (int h = 0; h < 3; h++)
					hkl_unique[h].erase(hkl_unique[h].begin() + i);
		}
		if (debug)
			file << "K-pt_unique size: " << k_pt_unique[0].size() << endl;
	}

	file << endl << "Number of k points to evaluate: " << k_pt[0].size() << " for " << points << " gridpoints." << endl;

	vector< vector < complex<double> > > sf;
	sf.resize(asym_atom_list.size());
	if (shrink)
#pragma omp parallel for
		for (int i = 0; i < asym_atom_list.size(); i++)
			sf[i].resize(hkl_unique[0].size());
	else
#pragma omp parallel for
		for (int i = 0; i < asym_atom_list.size(); i++)
			sf[i].resize(sym[0][0].size() * hkl[0].size());

	if (debug)
		file << "Initialized FFs" << endl
		<< "asym atom list size: " << asym_atom_list.size() << " total grid size: " << points << endl;

#ifdef _WIN64
	time_t end1 = time(NULL);

	if (end1 - start < 60) file << "Time to prepare: " << fixed << setprecision(0) << end1 - start << " s" << endl;
	else if (end1 - start < 3600) file << "Time to prepare: " << fixed << setprecision(0) << floor((end1 - start) / 60) << " m " << (end1 - start) % 60 << " s" << endl;
	else file << "Time to prepare: " << fixed << setprecision(0) << floor((end1 - start) / 3600) << " h " << ((end1 - start) % 3600) / 60 << " m" << endl;
#else
	gettimeofday(&t2, 0);

	double timea = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

	if (timea < 60) printf("Time to prepare: %4.1lf s\n", timea);
	else if (timea < 3600) printf("Time to prepare: %10.1lf m\n", timea / 60);
	else printf("Time to prepare: %10.1lf h\n", timea / 3600);

#endif

	progress_bar * progress = new progress_bar{ file, 60u, "Calculating scattering factors" };
	const int step = max(floor(asym_atom_list.size() / 20),1.0);
	const int smax = shrink ? k_pt_unique[0].size() : k_pt[0].size();
	const int imax = asym_atom_list.size();
	int pmax;
	double* dens_local, * d1_local, * d2_local, * d3_local;
	complex<double>* sf_local;
	const double* k1_local = shrink ? k_pt_unique[0].data() : k_pt[0].data();
	const double* k2_local = shrink ? k_pt_unique[1].data() : k_pt[1].data();
	const double* k3_local = shrink ? k_pt_unique[2].data() : k_pt[2].data();
	for (int i = 0; i < imax; i++) {
		pmax = dens[i].size();
		dens_local = dens[i].data();
		d1_local = d1[i].data();
		d2_local = d2[i].data();
		d3_local = d3[i].data();
		sf_local = sf[i].data();
#pragma omp parallel for
		for (int s = 0; s < smax; s++) {
			double work;
			double rho;
			for (int p = pmax-1; p >= 0; p--) {
				rho = dens_local[p];
				work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
				sf_local[s] += complex<double>(rho * cos(work) + rho * sin(work));
			}
		}
		if (i != 0 && i % step == 0)
			progress->write(i / double(imax));
	}
	delete(progress);

	if (electron_diffraction) {
		const double fact = 2*M_PI*9.1093837015E-31 * pow(1.602176634E-19, 2) / (pow(6.62607015E-34, 2) * 8.8541878128E-12);
		for (int s = 0; s < sym[0][0].size(); s++) {
#pragma omp parallel for
			for (int p = 0; p < hkl[0].size(); p++) {
				double h2 = pow(hkl[0][p] * sym[0][0][s] + hkl[1][p] * sym[0][1][s] + hkl[2][p] * sym[0][2][s], 2)
					+ pow(hkl[0][p] * sym[1][0][s] + hkl[1][p] * sym[1][1][s] + hkl[2][p] * sym[1][2][s], 2)
					+ pow(hkl[0][p] * sym[2][0][s] + hkl[1][p] * sym[2][1][s] + hkl[2][p] * sym[2][2][s], 2);
				for (int i = 0, imax = asym_atom_list.size(); i < imax; i++)
					sf[i][p + hkl[0].size() * s] = fact * ((double) wave.atoms[asym_atom_list[i]].charge - sf[i][p + hkl[0].size() * s]) / h2;
			}
		}
	}


#ifdef _WIN64
	end_tsc = time(NULL);

#else
	gettimeofday(&t2, 0);

	double time4 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

	if (time4 < 60) printf("Time to calculate: %4.1lf s\n", time4);
	else if (time4 < 3600) printf("Time to calculate: %10.1lf m\n", time4 / 60);
	else printf("Time to calculate: %10.1lf h\n", time4 / 3600);

#endif

	if (debug)
		file << endl << "SFs are made, now just write them!" << endl;
	else
		file << endl << "Writing tsc file..." << endl;

	ofstream tsc_file("experimental.tsc", ios::out);

	tsc_file << "TITLE: " << get_filename_from_path(cif).substr(0,cif.find(".cif")) << endl << "SYMM: ";

	if (shrink)
		tsc_file << "expanded";
	else
		for (int s = 0; s < sym[0][0].size(); s++) {
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) {
					if (i != 0 || j != 0) tsc_file << " ";
					tsc_file << sym[j][i][s];
				}
			if (s != sym[0][0].size() - 1)
				tsc_file << ";";
		}

	tsc_file << endl << "SCATTERERS:";
	for (int i = 0; i < asym_atom_list.size(); i++)
		tsc_file << " " << labels[i];
	tsc_file << endl << "DATA:" << endl;

	if (shrink)
		for (int r = 0; r < hkl_unique[0].size(); r++) {
			for (int h = 0; h < 3; h++)
				tsc_file << hkl_unique[h][r] << " ";
			for (int i = 0; i < asym_atom_list.size(); i++)
				tsc_file << scientific << setprecision(8) << real(sf[i][r]) << ","
				<< scientific << setprecision(8) << imag(sf[i][r]) << " ";
			tsc_file << endl;
		}

	else {
		if (debug) file << "Writing a total of " << sym[0][0].size() << " symmetry generated sets of hkls!" << endl;
		for (int s = 0; s < sym[0][0].size(); s++) {
			//if (debug) file << "writing symmetry: " << s << endl;
			for (int p = 0; p < hkl[0].size(); p++) {
				tsc_file
					<< hkl[0][p] * sym[0][0][s] + hkl[1][p] * sym[0][1][s] + hkl[2][p] * sym[0][2][s] << " "
					<< hkl[0][p] * sym[1][0][s] + hkl[1][p] * sym[1][1][s] + hkl[2][p] * sym[1][2][s] << " "
					<< hkl[0][p] * sym[2][0][s] + hkl[1][p] * sym[2][1][s] + hkl[2][p] * sym[2][2][s] << " ";
				for (int i = 0; i < asym_atom_list.size(); i++)
					tsc_file << scientific << setprecision(8) << real(sf[i][p + hkl[0].size() * s]) << ","
					<< scientific << setprecision(8) << imag(sf[i][p + hkl[0].size() * s]) << " ";
				tsc_file << endl;
			}
		}
	}
	tsc_file.close();

#ifdef _WIN64
	time_t end = time(NULL);

	if (end - start < 60) file << "Total Time: "<< fixed << setprecision(0) << end - start << " s\n";
	else if (end - start < 3600) file << "Total Time: " << fixed << setprecision(0) << floor((end - start)/60) << " m " << (end - start) % 60 << " s\n";
	else file << "Total Time: " << fixed << setprecision(0) << floor((end - start)/3600) << " h " << ((end - start) % 3600)/60 << " m\n";
	file << endl;
	file << "Time Breakdown:" << endl;
	file << " ... for Prototype Grid setup:" << setw(6) << end_prototypes - start << " s" << endl;
	file << " ... for Becke Grid setup:    " << setw(6) << end_becke - end_prototypes << " s" << endl;
	file << " ... for spherical density:   " << setw(6) << end_spherical - end_becke << " s" << endl;
	file << " ... for Grid Pruning:        " << setw(6) << end_prune - end_spherical << " s" << endl;
	file << " ... for aspherical density:  " << setw(6) << end_aspherical - end_prune << " s" << endl;
	file << " ... for final preparation:   " << setw(6) << end1 - end_aspherical << " s" << endl;
	file << " ... for tsc calculation:     " << setw(6) << end - end1 << " s" << endl;
#else
	gettimeofday(&t2, 0);

	double time2 = (1000000.0*(t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

	if (time2<60) printf("Total Time: %4.1lf s\n", time2);
	else if (time2 < 3600) printf("Total Time: %10.1lf m\n", time2 / 60);
	else printf("Total Time: %10.1lf h\n", time2 / 3600);

#endif

	return true;
}

//bool calculate_structure_factors_RI(
//	string& hkl_filename,
//	string& cif,
//	string& asym_cif,
//	string& symm,
//	WFN& wave,
//	bool debug,
//	int accuracy,
//	ofstream& file,
//	int cpus,
//	bool electron_diffraction,
//	int pbc
//)
//{
//	if (cpus != -1) {
//		omp_set_num_threads(cpus);
//		omp_set_dynamic(0);
//	}
//	int* atom_z = new int[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
//	double* x, * y, * z;
//	x = new double[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
//	y = new double[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
//	z = new double[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
//	double* alpha_max = new double[wave.get_ncen()];
//	int* max_l = new int[wave.get_ncen()];
//	int max_l_overall = 0;
//
//#ifdef _WIN64
//	time_t start = time(NULL);
//	time_t end_becke, end_prototypes, end_spherical = 0, end_prune = 0, end_aspherical, end_tsc;
//#else
//	struct timeval t1, t2;
//
//	gettimeofday(&t1, 0);
//#endif
//
//	if (debug)
//		file << "Reading hkl now" << endl;
//
//	vector< vector <int> > hkl;
//	string line;
//	hkl.resize(3);
//	if (!exists(hkl_filename)) {
//		file << "HKL file does not exists!" << endl;
//		return false;
//	}
//	ifstream hkl_input(hkl_filename.c_str(), ios::in);
//	hkl_input.seekg(0, hkl_input.beg);
//	regex r{ R"([abcdefghijklmnopqrstuvwxyz\(\)ABCDEFGHIJKLMNOPQRSTUVW])" };
//	while (!hkl_input.eof()) {
//		getline(hkl_input, line);
//		if (hkl_input.eof())
//			break;
//		if (line.size() < 2)
//			continue;
//		cmatch result;
//		if (regex_search(line.c_str(), result, r))
//			continue;
//		//if (debug) file << "hkl: ";
//		for (int i = 0; i < 3; i++) {
//			string temp = line.substr(4 * size_t(i) + 1, 3);
//			temp.erase(remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
//			hkl[i].push_back(stoi(temp));
//			//if (debug) file << setw(4) << temp;
//		}
//		//if (debug) file << endl;
//		if (hkl[0][hkl[0].size() - 1] == 0 && hkl[1][hkl[0].size() - 1] == 0 && hkl[2][hkl[0].size() - 1] == 0) {
//			if (debug) file << "popping back 0 0 0" << endl;
//			for (int i = 0; i < 3; i++)
//				hkl[i].pop_back();
//		}
//	}
//	hkl_input.close();
//	// Remove duplicate reflections
//	for (int i = 0; i < hkl[0].size(); i++)
//		for (int j = i + 1; j < hkl[0].size(); j++)
//			if (hkl[0][i] == hkl[0][j] && hkl[1][i] == hkl[1][j] && hkl[2][i] == hkl[2][j])
//				for (int x = 0; x < 3; x++)
//					hkl[x].erase(hkl[x].begin() + j);
//
//	if (debug)
//		file << "Reflections read! Nr of reflections: " << hkl[0].size() << endl;
//
//	double rcm[3][3];
//	double cm[3][3];
//	if (debug)
//		file << "starting to read cif!" << endl;
//	if (!exists(cif)) {
//		file << "CIF does not exists!" << endl;
//		return false;
//	}
//	if (!exists(asym_cif)) {
//		file << "Asym CIF does not exists!" << endl
//			<< asym_cif << endl;
//		return false;
//	}
//	ifstream cif_input(cif.c_str(), ios::in);
//	ifstream asym_cif_input(asym_cif.c_str(), ios::in);
//	vector<bool> found;
//	found.resize(7);
//	for (int k = 0; k < 7; k++)
//		found[k] = false;
//	double a = 0.0, b = 0.0, c = 0.0, v = 0.0;
//	double alpha = 0.0, beta = 0.0, gamma = 0.0;
//	vector <string> cell_keywords;
//	cell_keywords.push_back("_cell_length_a");
//	cell_keywords.push_back("_cell_length_b");
//	cell_keywords.push_back("_cell_length_c");
//	cell_keywords.push_back("_cell_angle_alpha");
//	cell_keywords.push_back("_cell_angle_beta");
//	cell_keywords.push_back("_cell_angle_gamma");
//	cell_keywords.push_back("_cell_volume");
//	if (debug)
//		file << "Starting while !.eof()" << endl;
//	while (!cif_input.eof()) {
//		if (debug)
//			file << "While line! " << setw(80) << line
//			<< setw(10) << a << found[0]
//			<< setw(10) << b << found[1]
//			<< setw(10) << c << found[2]
//			<< setw(10) << alpha << found[3]
//			<< setw(10) << beta << found[4]
//			<< setw(10) << gamma << found[5]
//			<< setw(10) << v << found[6] << endl;
//		getline(cif_input, line);
//		for (int k = 0; k < cell_keywords.size(); k++) {
//			if (line.find(cell_keywords[k]) != string::npos) {
//				switch (k) {
//				case 0:
//					a = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 1:
//					b = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 2:
//					c = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 3:
//					alpha = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 4:
//					beta = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 5:
//					gamma = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 6:
//					v = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				default:
//					file << "This is weird... should never get here... aborting!" << endl;
//					return false;
//				}
//				found[k] = true;
//			}
//		}
//		if (found[0] == true && found[1] == true && found[2] == true && found[3] == true && found[4] == true && found[5] == true && found[6] == true)
//			break;
//	}
//	double ca = cos(0.0174532925199432944444444444444 * alpha);
//	double cb = cos(0.0174532925199432944444444444444 * beta);
//	double cg = cos(0.0174532925199432944444444444444 * gamma);
//	double sa = sin(0.0174532925199432944444444444444 * alpha);
//	double sb = sin(0.0174532925199432944444444444444 * beta);
//	double sg = sin(0.0174532925199432944444444444444 * gamma);
//	double V = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg);
//	if (V / v > 1.1 || V / v < 0.9) {
//		file << "Volume computed is more than 10% off, please check!" << endl;
//		return false;
//	}
//	double a_star = b * c * sa / V;
//	double b_star = a * c * sb / V;
//	double c_star = a * b * sg / V;
//
//	if (debug)
//		file << "Making cm and rcm" << endl
//		<< ca << " " << cb << " " << cg << " " << sa << " " << sb << " " << sg << " " << V << endl
//		<< a_star << " " << b_star << " " << c_star << endl;
//
//	cm[0][0] = a / 0.529177249;
//	cm[0][1] = sqrt(abs(a * b * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
//	cm[0][2] = sqrt(abs(a * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));
//
//	cm[1][0] = sqrt(abs(a * b * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
//	cm[1][1] = b / 0.529177249;
//	cm[1][2] = sqrt(abs(b * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));
//
//	cm[2][0] = sqrt(abs(a * c * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
//	cm[2][1] = sqrt(abs(b * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));
//	cm[2][2] = c / 0.529177249;
//
//	rcm[0][0] = 2 * M_PI / a;
//	rcm[0][1] = 0;
//	rcm[0][2] = 0;
//
//	rcm[1][0] = 2 * M_PI * -cg / (a * sg);
//	rcm[1][1] = 2 * M_PI * 1 / (b * sg);
//	rcm[1][2] = 0;
//
//	rcm[2][0] = 2 * M_PI * b * c * (ca * cg - cb) / V / sg;
//	rcm[2][1] = 2 * M_PI * a * c * (cb * cg - ca) / V / sg;
//	rcm[2][2] = 2 * M_PI * a * b * sg / V;
//
//	for (int i = 0; i < 3; i++)
//		for (int j = 0; j < 3; j++)
//			if (abs(rcm[i][j]) < 10e-10) {
//				rcm[i][j] = 0.0;
//				//cm[i][j] = 0.0;
//			}
//			else {
//				rcm[i][j] *= 0.529177249;
//				//cm[i][j] *= 0.529177249;
//			}
//
//	if (debug) {
//		file << "RCM done, now labels and asym atoms!" << endl;
//		for (int i = 0; i < 3; ++i) {
//			for (int j = 0; j < 3; ++j)
//				file << setw(10) << fixed << rcm[i][j] / 2 / M_PI / 0.529177249 << ' ';
//			file << endl;
//		}
//		file << "CM in bohr:" << endl;
//		for (int i = 0; i < 3; ++i) {
//			for (int j = 0; j < 3; ++j)
//				file << setw(10) << fixed << cm[i][j] << ' ';
//			file << endl;
//		}
//	}
//	cif_input.clear();
//	cif_input.seekg(0, cif_input.beg);
//	vector <string> labels;
//	vector <int> atom_type_list;
//	vector <int> asym_atom_to_type_list;
//	int count_fields = 0;
//	int position_field[3] = { 0,0,0 };
//	int label_field = 1000;
//	vector <int> asym_atom_list;
//	vector <int> all_atom_list;
//	vector < bool > is_asym;
//	vector < vector <double > > positions;
//	positions.resize(wave.get_ncen());
//	is_asym.resize(wave.get_ncen());
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		is_asym[i] = false;
//		positions[i].resize(3);
//	}
//	bool atoms_read = false;
//	while (!asym_cif_input.eof() && !atoms_read) {
//		getline(asym_cif_input, line);
//		//if(debug) file << "line: "<< line << endl;
//		if (line.find("loop_") != string::npos) {
//			//if(debug) file << "found loop!" << endl;
//			while (line.find("_") != string::npos) {
//				getline(asym_cif_input, line);
//				if (debug) file << "line in loop field definition: " << line << endl;
//				if (line.find("label") != string::npos)
//					label_field = count_fields;
//				else if (line.find("fract_x") != string::npos)
//					position_field[0] = count_fields;
//				else if (line.find("fract_y") != string::npos)
//					position_field[1] = count_fields;
//				else if (line.find("fract_z") != string::npos)
//					position_field[2] = count_fields;
//				else if (label_field == 1000) {
//					if (debug) file << "I don't think this is the atom block.. moving on!" << endl;
//					break;
//				}
//				count_fields++;
//			}
//			while (line.find("_") == string::npos && line.length() > 3) {
//				//if(debug) file << "Reading atom!"<< endl;
//				atoms_read = true;
//				stringstream s(line);
//				vector <string> fields;
//				fields.resize(count_fields);
//				for (int i = 0; i < count_fields; i++)
//					s >> fields[i];
//				positions[labels.size()][0] = (a * stod(fields[position_field[0]]) + b * cg * stod(fields[position_field[1]]) + c * cb * stod(fields[position_field[2]])) / 0.529177249;
//				positions[labels.size()][1] = (b * sg * stod(fields[position_field[1]]) + c * (ca - cb * cg) / sg * stod(fields[position_field[2]])) / 0.529177249;
//				positions[labels.size()][2] = V / (a * b * sg) * stod(fields[position_field[2]]) / 0.529177249;
//				bool found_this_one = false;
//				if (debug) file << "label: " << fields[label_field] << " position: " << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
//				for (int i = 0; i < wave.get_ncen(); i++) {
//					if (is_similar(positions[labels.size()][0], wave.atoms[i].x, -1)
//						&& is_similar(positions[labels.size()][1], wave.atoms[i].y, -1)
//						&& is_similar(positions[labels.size()][2], wave.atoms[i].z, -1)) {
//						if (debug) file << "WFN position: " << wave.atoms[i].x << " " << wave.atoms[i].y << " " << wave.atoms[i].z << endl
//							<< "Found an atom: " << fields[label_field] << " Corresponding to atom charge " << wave.atoms[i].charge << endl;
//						all_atom_list.push_back(i);
//						found_this_one = true;
//						break;
//					}
//				}
//				if (!found_this_one && debug)
//					file << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
//				labels.push_back(fields[label_field]);
//				getline(asym_cif_input, line);
//			}
//		}
//	}
//	if (labels.size() != wave.get_ncen()) {
//		file << "Number of atoms in labels: " << labels.size() << " and number of atoms in Wavefunction: " << wave.get_ncen() << "!" << endl << "This is not good, i will stop here!" << endl;
//		return false;
//	}
//	atoms_read = false;
//	label_field = 1000;
//	count_fields = 0;
//	while (!cif_input.eof() && !atoms_read) {
//		getline(cif_input, line);
//		//if(debug) file << "line: "<< line << endl;
//		if (line.find("loop_") != string::npos) {
//			//if(debug) file << "found loop!" << endl;
//			while (line.find("_") != string::npos) {
//				getline(cif_input, line);
//				if (debug) file << "line in loop field definition: " << line << endl;
//				if (line.find("label") != string::npos)
//					label_field = count_fields;
//				else if (label_field == 1000) {
//					if (debug) file << "I don't think this is the atom block.. moving on!" << endl;
//					break;
//				}
//				count_fields++;
//			}
//			while (line.find("_") == string::npos && line.length() > 3) {
//				atoms_read = true;
//				stringstream s(line);
//				vector <string> fields;
//				fields.resize(count_fields);
//				for (int i = 0; i < count_fields; i++)
//					s >> fields[i];
//				fields[label_field].erase(remove_if(fields[label_field].begin(), fields[label_field].end(), ::isspace), fields[label_field].end());
//				for (int atom = 0; atom < wave.get_ncen(); atom++) {
//					if (debug) file << "Comparing atoms: " << fields[label_field] << " / " << labels[atom] << endl;
//					if (fields[label_field] == labels[atom]) {
//						int nr = -1;
//						for (int i = 0; i < wave.get_ncen(); i++) {
//							if (is_similar(positions[atom][0], wave.atoms[i].x, -1)
//								&& is_similar(positions[atom][1], wave.atoms[i].y, -1)
//								&& is_similar(positions[atom][2], wave.atoms[i].z, -1)) {
//								if (debug) file << "Found an asymmetric atom: " << fields[label_field] << " Corresponding to atom charge " << wave.atoms[i].charge << endl;
//								nr = i;
//								asym_atom_list.push_back(i);
//								break;
//							}
//						}
//						is_asym[atom] = true;
//						bool already_there = false;
//						for (int i = 0; i < atom_type_list.size(); i++)
//							if (atom_type_list[i] == wave.atoms[nr].charge) {
//								already_there = true;
//								asym_atom_to_type_list.push_back(i);
//							}
//						if (already_there == false) {
//							asym_atom_to_type_list.push_back(atom_type_list.size());
//							atom_type_list.push_back(wave.atoms[nr].charge);
//						}
//						break;
//					}
//				}
//				getline(cif_input, line);
//			}
//		}
//	}
//
//	if (debug) {
//		file << "There are " << atom_type_list.size() << " types of atoms" << endl;
//		for (int i = 0; i < atom_type_list.size(); i++)
//			file << setw(4) << atom_type_list[i];
//		file << endl << "asym_atoms_to_type_list: " << endl;
//		for (int i = 0; i < asym_atom_to_type_list.size(); i++)
//			file << setw(4) << asym_atom_to_type_list[i];
//		file << endl;
//		file << "Mapping of asym atoms:" << endl;
//		for (int i = 0; i < wave.get_ncen(); i++)
//			file << setw(4) << wave.atoms[all_atom_list[i]].charge;
//		file << endl;
//		for (int i = 0; i < wave.get_ncen(); i++)
//			file << setw(4) << is_asym[i];
//		file << endl;
//	}
//
//	bool symm_read = false;
//	//Still need to read the sym matrices
//	if (symm == "") {
//		if (debug) file << "No Symmetry file specified, tyring to read from the CIF!" << endl;
//		symm_read = true;
//	}
//	else if (!exists(symm) && symm_read == true)
//		return false;
//
//	vector < vector < vector <int> > > sym;
//	sym.resize(3);
//	for (int i = 0; i < 3; i++)
//		sym[i].resize(3);
//	if (!symm_read) {
//		ifstream symm_input(symm.c_str(), ios::in);
//		string liny;
//		int temp_int;
//		while (!symm_input.eof()) {
//			getline(symm_input, liny);
//			stringstream streamy(liny);
//			for (int i = 0; i < 3; i++)
//				for (int j = 0; j < 3; j++) {
//					streamy >> temp_int;
//					sym[i][j].push_back(temp_int);
//				}
//		}
//	}
//	else {
//		cif_input.clear();
//		cif_input.seekg(0, cif_input.beg);
//		bool symm_found = false;
//		int operation_field = 200;
//		count_fields = 0;
//		while (!cif_input.eof() && !symm_found) {
//			getline(cif_input, line);
//			if (line.find("loop_") != string::npos) {
//				//if(debug) file << "found loop!" << endl;
//				while (line.find("_") != string::npos) {
//					getline(cif_input, line);
//					if (debug) file << "line in loop field definition: " << line << endl;
//					if (line.find("space_group_symop_operation_xyz") != string::npos)
//						operation_field = count_fields;
//					else if (count_fields > 2 || (operation_field == 200 && count_fields != 0)) {
//						if (debug) file << "I don't think this is the symmetry block.. moving on!" << endl;
//						count_fields = 0;
//						break;
//					}
//					count_fields++;
//				}
//				while (line.find("_") == string::npos && line.length() > 3 && count_fields != 0) {
//					if (debug) file << "Reading operation!" << line << endl;
//					symm_found = true;
//					stringstream s(line);
//					vector <string> fields;
//					fields.resize(count_fields);
//					int sym_from_cif[3][3];
//					for (int x = 0; x < 3; x++)
//						for (int y = 0; y < 3; y++)
//							sym_from_cif[x][y] = 0;
//					for (int i = 0; i < count_fields; i++)
//						s >> fields[i];
//					vector<string> vectors;
//					vectors.resize(3);
//					int column = 0;
//					for (int c = 0; c < fields[operation_field].length(); c++) {
//						if (fields[operation_field][c] != ',')
//							vectors[column].push_back(fields[operation_field][c]);
//						else column++;
//					}
//
//					for (int x = 0; x < 3; x++) {
//						if (vectors[x].find("X") != string::npos || vectors[x].find("x") != string::npos) {
//							char sign = ' ';
//							if (vectors[x].find("X") != string::npos && vectors[x].find("X") != 0)
//								sign = vectors[x].at(vectors[x].find("X") - 1);
//							else if (vectors[x].find("X") == 0)
//								sign = '+';
//							if (vectors[x].find("x") != string::npos && vectors[x].find("x") != 0)
//								sign = vectors[x].at(vectors[x].find("x") - 1);
//							else if (vectors[x].find("x") == 0)
//								sign = '+';
//							if (sign == '-')
//								sym_from_cif[x][0] = -1;
//							if (sign == '+')
//								sym_from_cif[x][0] = 1;
//						}
//						if (vectors[x].find("Y") != string::npos || vectors[x].find("y") != string::npos) {
//							char sign = ' ';
//							if (vectors[x].find("Y") != string::npos && vectors[x].find("Y") != 0)
//								sign = vectors[x].at(vectors[x].find("Y") - 1);
//							else if (vectors[x].find("Y") == 0)
//								sign = '+';
//							if (vectors[x].find("y") != string::npos && vectors[x].find("y") != 0)
//								sign = vectors[x].at(vectors[x].find("y") - 1);
//							else if (vectors[x].find("y") == 0)
//								sign = '+';
//							if (sign == '-')
//								sym_from_cif[x][1] = -1;
//							if (sign == '+')
//								sym_from_cif[x][1] = 1;
//						}
//						if (vectors[x].find("Z") != string::npos || vectors[x].find("z") != string::npos) {
//							char sign = ' ';
//							if (vectors[x].find("Z") != string::npos && vectors[x].find("Z") != 0)
//								sign = vectors[x].at(vectors[x].find("Z") - 1);
//							else if (vectors[x].find("Z") == 0)
//								sign = '+';
//							if (vectors[x].find("z") != string::npos && vectors[x].find("z") != 0)
//								sign = vectors[x].at(vectors[x].find("z") - 1);
//							else if (vectors[x].find("z") == 0)
//								sign = '+';
//							if (sign == '-')
//								sym_from_cif[x][2] = -1;
//							if (sign == '+')
//								sym_from_cif[x][2] = 1;
//						}
//					}
//					if (debug) {
//						file << "Comparing ";
//						for (int x = 0; x < 3; x++)
//							for (int y = 0; y < 3; y++)
//								file << sym_from_cif[x][y] << " ";
//						file << endl;
//					}
//					bool already_known = false;
//					for (int s = 0; s < sym[0][0].size(); s++) {
//						bool identical = true;
//						bool inverse = true;
//						for (int x = 0; x < 3; x++)
//							for (int y = 0; y < 3; y++)
//								if (sym[y][x][s] != sym_from_cif[x][y])
//									identical = false;
//						if (!identical)
//							for (int x = 0; x < 3; x++)
//								for (int y = 0; y < 3; y++)
//									if (sym[y][x][s] != sym_from_cif[x][y] * -1)
//										inverse = false;
//						/*if (debug) {
//							file << "Comparison with ";
//							for (int x = 0; x < 3; x++)
//								for (int y = 0; y < 3; y++)
//									file << sym[x][y][s] << " ";
//							file << "resulted in ";
//							file << identical << " " << inverse << endl;
//						}*/
//						if (identical || inverse) {
//							already_known = true;
//							break;
//						}
//					}
//					if (!already_known) {
//						if (debug) file << "This is a new symmetry operation!" << endl;
//						for (int x = 0; x < 3; x++)
//							for (int y = 0; y < 3; y++)
//								sym[y][x].push_back(sym_from_cif[x][y]);
//					}
//					getline(cif_input, line);
//				}
//			}
//		}
//	}
//
//	if (debug) {
//		file << "Read " << sym[0][0].size() << " symmetry elements! Size of sym: " << sym[0][0].size() << endl;
//		for (int i = 0; i < sym[0][0].size(); i++) {
//			for (int x = 0; x < 3; x++) {
//				for (int y = 0; y < 3; y++)
//					file << setw(3) << sym[y][x][i];
//				file << endl;
//			}
//			file << endl;
//		}
//	}
//
//	cif_input.close();
//
//	if (debug)
//		file << "There are " << atom_type_list.size() << " Types of atoms and " << asym_atom_to_type_list.size() << " atoms in total" << endl;
//
//	if (asym_atom_list.size() == 0) {
//		file << "0 asym atoms is imposible! something is wrong with reading the CIF!" << endl;
//		return false;
//	}
//
//	if (debug)
//		file << "made it post CIF, now make grids!" << endl;
//
//	double*** grid = new double** [asym_atom_list.size()];
//	int* num_points = new int[asym_atom_list.size()];
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		grid[i] = new double* [6];
//	// GRID COORDINATES for [a][c][p] a = atom [0,ncen],
//	// c = coordinate [0=x, 1=y, 2=z, 3=atomic becke weight and molecular becke weight, 5=total spherical density],
//	// p = point in this grid
//
//#pragma omp parallel for
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		atom_z[i] = wave.atoms[i].charge;
//		x[i] = wave.atoms[i].x;
//		y[i] = wave.atoms[i].y;
//		z[i] = wave.atoms[i].z;
//		//if(debug)
//		//    file << "xyz= 000 position: " << x[i] << " " << y[i] << " " << z[i] << " Charge: " << atom_z[i] << endl;
//		if (pbc != 0) {
//			int j = 0;
//			for (int pbc_x = -pbc; pbc_x < pbc + 1; pbc_x++)
//				for (int pbc_y = -pbc; pbc_y < pbc + 1; pbc_y++)
//					for (int pbc_z = -pbc; pbc_z < pbc + 1; pbc_z++) {
//						if (pbc_x == 0 && pbc_y == 0 && pbc_z == 0)
//							continue;
//						else {
//							j++;
//							atom_z[i + j * wave.get_ncen()] = wave.atoms[i].charge;
//							x[i + j * wave.get_ncen()] = wave.atoms[i].x + pbc_x * cm[0][0] + pbc_y * cm[0][1] + pbc_z * cm[0][2];
//							y[i + j * wave.get_ncen()] = wave.atoms[i].y + pbc_x * cm[1][0] + pbc_y * cm[1][1] + pbc_z * cm[1][2];
//							z[i + j * wave.get_ncen()] = wave.atoms[i].z + pbc_x * cm[2][0] + pbc_y * cm[2][1] + pbc_z * cm[2][2];
//							if (debug)
//								file << "xyz= " << pbc_x << pbc_y << pbc_z << " j = " << j << " position: " << x[i + j * wave.get_ncen()] << " " << y[i + j * wave.get_ncen()] << " " << z[i + j * wave.get_ncen()] << " Charge: " << atom_z[i + j * wave.get_ncen()] << endl;
//						}
//					}
//		}
//		alpha_max[i] = 0.0;
//		max_l[i] = 0;
//		for (int b = 0; b < wave.get_nex(); b++) {
//			if (wave.get_center(b) != i + 1)
//				continue;
//			if (wave.get_exponent(b) > alpha_max[i])
//				alpha_max[i] = wave.get_exponent(b);
//			if (wave.get_type(b) > max_l[i]) {
//				int l = wave.get_type(b);
//				if (l == 1)
//					l = 1;
//				else if (l >= 2 && l <= 4)
//					l = 2;
//				else if (l >= 5 && l <= 10)
//					l = 3;
//				else if (l >= 11 && l <= 20)
//					l = 4;
//				else if (l >= 21 && l <= 35)
//					l = 5;
//				max_l[i] = l;
//				if (l > max_l_overall)
//					max_l_overall = l;
//			}
//		}
//	}
//
//	if (debug)
//		file << "Atoms are there!" << endl;
//
//	double** alpha_min = new double* [wave.get_ncen()];
//	for (int i = 0; i < wave.get_ncen(); i++)
//		alpha_min[i] = new double[max_l_overall];
//
//#pragma omp parallel for
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		for (int b = 0; b < max_l_overall; b++)
//			alpha_min[i][b] = 100000000.0;
//	}
//
//#pragma omp parallel for
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		for (int b = 0; b < wave.get_nex(); b++) {
//			if (wave.get_center(b) != i + 1)
//				continue;
//			int l = wave.get_type(b);
//			if (l == 1)
//				l = 1;
//			else if (l >= 2 && l <= 4)
//				l = 2;
//			else if (l >= 5 && l <= 10)
//				l = 3;
//			else if (l >= 11 && l <= 20)
//				l = 4;
//			else if (l >= 21 && l <= 35)
//				l = 5;
//			else if (l >= 36 && l <= 56)
//				l = 6;
//			if (wave.get_exponent(b) < alpha_min[i][l - 1])
//				alpha_min[i][l - 1] = wave.get_exponent(b);
//		}
//		if (debug) {
//			file << "alpha_min: ";
//			for (int b = 0; b < max_l_overall; b++)
//				file << setw(14) << scientific << alpha_min[i][b];
//			file << endl;
//		}
//	}
//
//	if (debug)
//		file << "alpha_min is there!" << endl
//		<< "Nr of asym atoms: " << asym_atom_list.size() << " Number of all atoms: " << all_atom_list.size() << " Number of atoms in wfn: " << wave.get_ncen() << endl;
//	else
//		file << "There are:\n" << setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
//		//<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
//		<< setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << endl;
//	// Total grid as a sum of all atomic grids.
//	// Dimensions: [c] [p]
//	// p = the number of gridpoint
//	// c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
//	vector < vector <  double > > total_grid;
//
//	if (debug)
//		file << "Making Becke Grid for" << endl;
//	else
//		file << endl << "Making Becke Grids..." << flush;
//
//	total_grid.resize(6);
//
//	//Make Prototype grids with only single atom weights for all elements
//	vector <AtomGrid> Prototype_grids;
//	if (debug) file << "max_l_overall: " << max_l_overall << endl;
//	for (int i = 0; i < atom_type_list.size(); i++) {
//		if (debug) file << "Atom " << i << endl;
//		double alpha_max_temp;
//		double max_l_temp;
//		double* alpha_min_temp = new double[max_l_overall];
//		for (int j = 0; j < wave.get_ncen(); j++) {
//			if (wave.atoms[j].charge == atom_type_list[i]) {
//				if (debug) {
//					file << alpha_max[j] << " " << max_l[j] - 1 << " ";
//					for (int l = 0; l < max_l_overall; l++)
//						file << alpha_min[j][l] << " ";
//					file << endl;
//				}
//				alpha_max_temp = alpha_max[j];
//				max_l_temp = max_l[j] - 1;
//				for (int l = 0; l <= max_l_temp; l++)
//					alpha_min_temp[l] = alpha_min[j][l];
//				break;
//			}
//		}
//
//		if (debug) {
//			file << "max_l: " << defaultfloat << max_l_temp << " alpha_max: " << scientific << alpha_max_temp << " alpha_min: ";
//			for (int l = 0; l <= max_l_temp; l++)
//				file << setw(14) << scientific << alpha_min_temp[l];
//			file << " accuracy: " << accuracy << endl;
//		}
//		int lebedev_high, lebedev_low;
//		double radial_acc;
//		if (accuracy == 0) {
//			lebedev_high = (max_l_temp < 3) ? 0 : 10;
//			lebedev_low = (max_l_temp < 3) ? 0 : 10;
//			radial_acc = 1e-5;
//		}
//		else if (accuracy == 1) {
//			lebedev_high = (max_l_temp < 3) ? 110 : 146;
//			lebedev_low = (max_l_temp < 3) ? 38 : 50;
//			radial_acc = 1e-8;
//		}
//		else if (accuracy == 2) {
//			lebedev_high = (max_l_temp < 3) ? 230 : 266;
//			lebedev_low = (max_l_temp < 3) ? 110 : 146;
//			radial_acc = 1e-10;
//		}
//		else if (accuracy == 3) {
//			lebedev_high = (max_l_temp < 3) ? 350 : 590;
//			lebedev_low = (max_l_temp < 3) ? 266 : 350;
//			radial_acc = 1e-15;
//		}
//		else if (accuracy > 3) {
//			lebedev_high = (max_l_temp < 3) ? 590 : 5810;
//			lebedev_low = (max_l_temp < 3) ? 350 : 4802;
//			radial_acc = 1e-20;
//		}
//		Prototype_grids.push_back(AtomGrid(radial_acc,
//			lebedev_low,
//			lebedev_high,
//			atom_type_list[i],
//			alpha_max_temp,
//			max_l_temp,
//			alpha_min_temp,
//			debug,
//			file));
//
//	}
//
//#ifdef _WIN64
//	end_prototypes = time(NULL);
//	if (debug) {
//
//		for (int prototype = 0; prototype < Prototype_grids.size(); prototype++)
//			file << "Number of gridpoints for atom type " << atom_type_list[prototype] << " :" << Prototype_grids[prototype].get_num_grid_points() << endl;
//
//		//	int diff = end - start;
//		if (end_prototypes - start < 1) file << "Time until prototypes are done: <1 s" << endl;
//		else if (end_prototypes - start < 60) file << "Time until prototypes are done: " << fixed << setprecision(0) << end_prototypes - start << " s" << endl;
//		else if (end_prototypes - start < 3600) file << "Time until prototypes are done: " << fixed << setprecision(0) << floor((end_prototypes - start) / 60) << " m " << (end_prototypes - start) % 60 << " s" << endl;
//		else file << "Time until prototypes are done: " << fixed << setprecision(0) << floor((end_prototypes - start) / 3600) << " h " << ((end_prototypes - start) % 3600) / 60 << " m" << endl;
//	}
//#else
//	if (debug) {
//		gettimeofday(&t2, 0);
//
//		double time3 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//		if (time3 < 60) printf("Time to prepare: %4.1lf s\n", time3);
//		else if (time3 < 3600) printf("Time to prepare: %10.1lf m\n", time3 / 60);
//		else printf("Time to prepare: %10.1lf h\n", time3 / 3600);
//	}
//#endif
//
//#pragma omp parallel for schedule(dynamic)
//	for (int i = 0; i < asym_atom_list.size(); i++) {
//		int nr = asym_atom_list[i];
//		int type;
//		for (int i = 0; i < atom_type_list.size(); i++)
//			if (atom_type_list[i] == wave.atoms[nr].charge)
//				type = i;
//
//		num_points[i] = Prototype_grids[type].get_num_grid_points();
//
//		if (debug) file << "Number of points for atom " << i << ": " << num_points[i] << endl;
//
//		for (int n = 0; n < 6; n++)
//			grid[i][n] = new double[num_points[i]];
//		for (int p = 0; p < num_points[i]; p++)
//			grid[i][4][p] = 0.0;
//
//		Prototype_grids[type].get_grid(wave.get_ncen() * pow(pbc * 2 + 1, 3),
//			nr,
//			x,
//			y,
//			z,
//			atom_z,
//			grid[i][0],
//			grid[i][1],
//			grid[i][2],
//			grid[i][3],
//			grid[i][5]);
//	}
//
//	Prototype_grids.clear();
//	int points = 0;
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		points += num_points[i];
//	if (debug) file << "Becke Grid exists" << endl;
//	else file << "                                 done! Number of gridpoints: " << defaultfloat << points << endl;
//
//#ifdef _WIN64
//	if (debug) {
//		file << "Taking time..." << endl;
//	}
//	end_becke = time(NULL);
//
//#else
//	if (debug) {
//		gettimeofday(&t2, 0);
//
//		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
//		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
//		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
//	}
//#endif
//
//	if (debug) {
//		file << "Asym atom list: ";
//		for (int g = 0; g < asym_atom_list.size(); g++)
//			file << setw(4) << asym_atom_list[g];
//		file << endl;
//		file << "All atom list: ";
//		for (int g = 0; g < all_atom_list.size(); g++)
//			file << setw(4) << all_atom_list[g];
//		file << endl;
//	}
//	for (int i = 0; i < asym_atom_list.size(); i++) {
//		int reduction = 0;
//		for (int p = 0; p < num_points[i]; p++) {
//			for (int k = 0; k < 4; k++)
//				total_grid[k].push_back(grid[i][k][p]);
//			total_grid[5].push_back(grid[i][5][p]);
//
//		}
//		num_points[i] -= reduction;
//		delete[](grid[i]);
//		if (debug) file << endl << "number of points for atom " << i << " " << num_points[i]<< endl;
//	}
//
//
//	file << "Calculating aspherical densities..." << flush;
//	vector < vector < double > > periodic_grid;
//
//	total_grid[4].resize(total_grid[0].size());
//#pragma omp parallel for
//	for (int i = 0; i < total_grid[0].size(); i++) {
//		total_grid[4][i] = compute_dens(wave, new double[3]{
//			total_grid[0][i],
//			total_grid[1][i],
//			total_grid[2][i] });
//	}
//	if (pbc != 0) {
//		periodic_grid.resize(pow(pbc * 2 + 1, 3));
//		int j = 0;
//		for (int d = 0; d < pow(pbc * 2 + 1, 3); d++)
//			periodic_grid[d].resize(total_grid[5].size());
//		for (int x = -pbc; x < pbc + 1; x++)
//			for (int y = -pbc; y < pbc + 1; y++)
//				for (int z = -pbc; z < pbc + 1; z++) {
//					if (x == 0 && y == 0 && z == 0)
//						continue;
//#pragma omp parallel for
//					for (int i = 0; i < total_grid[0].size(); i++)
//						periodic_grid[j][i] = compute_dens(wave, new double[3]{
//							total_grid[0][i] + x * cm[0][0] + y * cm[0][1] + z * cm[0][2],
//							total_grid[1][i] + x * cm[1][0] + y * cm[1][1] + z * cm[1][2],
//							total_grid[2][i] + x * cm[2][0] + y * cm[2][1] + z * cm[2][2] });
//					j++;
//				}
//		if (debug) {
//			for (int i = 0; i < total_grid[0].size(); i++) {
//				if (i % 1000 == 0)
//					file << "Old dens: " << total_grid[5][i] << " contributions of neighbour-cells:";
//				for (int j = 0; j < pow(pbc * 2 + 1, 3) - 1; j++) {
//					if (i % 1000 == 0)
//						file << " " << periodic_grid[j][i];
//					total_grid[4][i] += periodic_grid[j][i];
//				}
//				if (i % 1000 == 0)
//					file << endl;
//			}
//		}
//	}
//
//	if (debug) file << endl << "with total number of points: " << total_grid[0].size() << endl;
//	else file << "                   done!" << endl;
//
//	file << "Applying hirshfeld weights and integrating charges..." << flush;
//	double el_sum_becke = 0.0;
//	// Vector containing integrated numbers of electrons
//	// dimension 0: 0=Becek grid integration 1=Summed spherical density 2= hirshfeld weighted density
//	// dimension 1: atoms of asym_atom_list
//	vector < vector <double> > atom_els;
//	atom_els.resize(3);
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		for (int n = 0; n < 3; n++)
//			atom_els[n].push_back(0.0);
//
//	//Generate Electron sums
//	for (int i = 0; i < asym_atom_list.size(); i++) {
//		int start_p = 0;
//#pragma loop(no_vector)
//		for (int a = 0; a < i; a++)
//			start_p += num_points[a];
//		for (int p = start_p; p < start_p + num_points[i]; p++) {
//			if (abs(total_grid[6][p]) > 1E-20) {
//				atom_els[0][i] += total_grid[5][p] * total_grid[4][p];
//			}
//		}
//		el_sum_becke += atom_els[0][i];
//	}
//
//	if (debug)
//		file << "Becke grid done!" << endl;
//
//#ifdef _WIN64
//	if (debug) {
//		file << "Taking time..." << endl;
//	}
//	end_aspherical = time(NULL);
//
//#else
//	if (debug) {
//		gettimeofday(&t2, 0);
//
//		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
//		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
//		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
//	}
//#endif
//
//#pragma omp parallel for
//	for (int p = 0; p < total_grid[0].size(); p++)
//		total_grid[5][p] *= total_grid[3][p];
//
//	file << " done!" << endl;
//	file << "Number of points evaluated: " << total_grid[0].size() << " with " << el_sum_becke << " electrons in Becke Grid in total." << endl << endl;
//
//	file << "Table of Charges in electrons" << endl << endl << "Atom       Becke   Spherical Hirshfeld" << endl;
//
//	int counter = 0;
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		if (is_asym[i]) {
//			file << setw(6) << labels[i]
//				<< fixed << setw(10) << setprecision(3) << wave.atoms[all_atom_list[i]].charge - atom_els[0][counter];
//			if (debug) file << " " << wave.atoms[all_atom_list[i]].charge << " " << fixed << setw(10) << setprecision(3) << wave.atoms[all_atom_list[i]].charge - atom_els[0][counter];
//			counter++;
//			file << endl;
//		}
//	}
//
//	file << "Total number of electrons in the wavefunction: " << el_sum_becke << endl;
//
//	vector < vector <double> > d1, d2, d3;
//	vector < vector <double> > dens;
//	points = 0;
//	dens.resize(asym_atom_list.size());
//	d1.resize(asym_atom_list.size());
//	d2.resize(asym_atom_list.size());
//	d3.resize(asym_atom_list.size());
//	if (debug)
//		file << "resized outer d1-3" << endl;
//
//#pragma omp parallel for reduction(+:points)
//	for (int i = 0; i < asym_atom_list.size(); i++) {
//		int start_p = 0;
//#pragma loop(no_vector)
//		for (int a = 0; a < i; a++)
//			start_p += num_points[a];
//		for (int p = start_p; p < start_p + num_points[i]; p++) {
//
//			dens[i].push_back(total_grid[5][p]);
//			d1[i].push_back(total_grid[0][p] - wave.atoms[asym_atom_list[i]].x);
//			d2[i].push_back(total_grid[1][p] - wave.atoms[asym_atom_list[i]].y);
//			d3[i].push_back(total_grid[2][p] - wave.atoms[asym_atom_list[i]].z);
//			
//		}
//		points += dens[i].size();
//	}
//
//	for (int grid = 0; grid < total_grid.size(); grid++)
//		total_grid[grid].clear();
//	total_grid.clear();
//
//	vector < vector <double> > k_pt;
//	k_pt.resize(3);
//#pragma omp parallel for
//	for (int i = 0; i < 3; i++)
//		k_pt[i].resize(sym[0][0].size() * hkl[0].size());
//
//	if (debug)
//		file << "K_point_vector is here" << endl;
//
//	for (int s = 0; s < sym[0][0].size(); s++) {
//#pragma omp parallel for
//		for (int ref = 0; ref < hkl[0].size(); ref++)
//			for (int x = 0; x < 3; x++)
//				for (int h = 0; h < 3; h++) {
//					double rcm_sym = 0.0;
//					for (int j = 0; j < 3; j++)
//						rcm_sym += rcm[x][j] * sym[j][h][s];
//					k_pt[x][ref + s * hkl[0].size()] += rcm_sym * hkl[h][ref];
//				}
//	}
//
//	file << endl << "Number of k points to evaluate: " << k_pt[0].size() << " for " << points << " gridpoints." << endl;
//
//	vector< vector < complex<double> > > sf;
//	sf.resize(asym_atom_list.size());
//
//#pragma omp parallel for
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		sf[i].resize(sym[0][0].size() * hkl[0].size());
//
//	if (debug)
//		file << "Initialized FFs" << endl
//		<< "asym atom list size: " << asym_atom_list.size() << " total grid size: " << points << endl;
//
//#ifdef _WIN64
//	time_t end1 = time(NULL);
//
//	if (end1 - start < 60) file << "Time to prepare: " << fixed << setprecision(0) << end1 - start << " s" << endl;
//	else if (end1 - start < 3600) file << "Time to prepare: " << fixed << setprecision(0) << floor((end1 - start) / 60) << " m " << (end1 - start) % 60 << " s" << endl;
//	else file << "Time to prepare: " << fixed << setprecision(0) << floor((end1 - start) / 3600) << " h " << ((end1 - start) % 3600) / 60 << " m" << endl;
//#else
//	gettimeofday(&t2, 0);
//
//	double timea = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//	if (timea < 60) printf("Time to prepare: %4.1lf s\n", timea);
//	else if (timea < 3600) printf("Time to prepare: %10.1lf m\n", timea / 60);
//	else printf("Time to prepare: %10.1lf h\n", timea / 3600);
//
//#endif
//
//	progress_bar* progress = new progress_bar{ file, 60u, "Calculating scattering factors" };
//	const int step = max(floor(asym_atom_list.size() / 20), 1.0);
//	int smax = k_pt[0].size();
//	int imax = (int)asym_atom_list.size();
//	for (int i = 0; i < imax; i++) {
//		int pmax = (int)dens[i].size();
//#pragma omp parallel for
//		for (int s = 0; s < smax; s++) {
//			double work;
//			double rho;
//			for (int p = pmax - 1; p >= 0; p--) {
//				rho = dens[i][p];
//				work = k_pt[0][s] * d1[i][p] + k_pt[1][s] * d2[i][p] + k_pt[2][s] * d3[i][p];
//				sf[i][s] += complex<double>(rho * cos(work), rho * sin(work));
//			}
//		}
//		if (i != 0 && i % step == 0)
//			progress->write(i / double(imax));
//	}
//	delete(progress);
//
//	if (electron_diffraction) {
//		double fact = 2 * M_PI * 9.1093837015E-31 * pow(1.602176634E-19, 2) / (pow(6.62607015E-34, 2) * 8.8541878128E-12);
//		for (unsigned int s = 0; s < sym[0][0].size(); s++) {
//#pragma omp parallel for
//			for (int p = 0; p < (int)hkl[0].size(); p++) {
//				double h2 = pow(hkl[0][p] * sym[0][0][s] + hkl[1][p] * sym[0][1][s] + hkl[2][p] * sym[0][2][s], 2)
//					+ pow(hkl[0][p] * sym[1][0][s] + hkl[1][p] * sym[1][1][s] + hkl[2][p] * sym[1][2][s], 2)
//					+ pow(hkl[0][p] * sym[2][0][s] + hkl[1][p] * sym[2][1][s] + hkl[2][p] * sym[2][2][s], 2);
//				for (unsigned int i = 0, imax = asym_atom_list.size(); i < imax; i++)
//					sf[i][p + hkl[0].size() * s] = fact * ((double)wave.atoms[asym_atom_list[i]].charge - sf[i][p + hkl[0].size() * s]) / h2;
//			}
//		}
//	}
//
//
//#ifdef _WIN64
//	end_tsc = time(NULL);
//
//#else
//	gettimeofday(&t2, 0);
//
//	double time4 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//	if (time4 < 60) printf("Time to calculate: %4.1lf s\n", time4);
//	else if (time4 < 3600) printf("Time to calculate: %10.1lf m\n", time4 / 60);
//	else printf("Time to calculate: %10.1lf h\n", time4 / 3600);
//
//#endif
//
//	if (debug)
//		file << endl << "SFs are made, now just write them!" << endl;
//	else
//		file << endl << "Writing tsc file..." << endl;
//
//	ofstream tsc_file("experimental.tsc", ios::out);
//
//	tsc_file << "TITLE: " << get_filename_from_path(cif).substr(0, cif.find(".cif")) << endl << "SYMM: ";
//
//	for (unsigned int s = 0; s < sym[0][0].size(); s++) {
//		for (int i = 0; i < 3; i++)
//			for (int j = 0; j < 3; j++) {
//				if (i != 0 || j != 0) tsc_file << " ";
//				tsc_file << sym[j][i][s];
//			}
//		if (s != sym[0][0].size() - 1)
//			tsc_file << ";";
//	}
//	tsc_file << endl << "SCATTERERS:";
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		tsc_file << " " << labels[i];
//	tsc_file << endl << "DATA:" << endl;
//
//	if (debug) file << "Writing a total of " << sym[0][0].size() << " symmetry generated sets of hkls!" << endl;
//	for (unsigned int s = 0; s < sym[0][0].size(); s++) {
//		//if (debug) file << "writing symmetry: " << s << endl;
//		for (unsigned int p = 0; p < hkl[0].size(); p++) {
//			tsc_file
//				<< hkl[0][p] * sym[0][0][s] + hkl[1][p] * sym[0][1][s] + hkl[2][p] * sym[0][2][s] << " "
//				<< hkl[0][p] * sym[1][0][s] + hkl[1][p] * sym[1][1][s] + hkl[2][p] * sym[1][2][s] << " "
//				<< hkl[0][p] * sym[2][0][s] + hkl[1][p] * sym[2][1][s] + hkl[2][p] * sym[2][2][s] << " ";
//			for (unsigned int i = 0; i < asym_atom_list.size(); i++)
//				tsc_file << scientific << setprecision(8) << real(sf[i][p + hkl[0].size() * s]) << ","
//				<< scientific << setprecision(8) << imag(sf[i][p + hkl[0].size() * s]) << " ";
//			tsc_file << endl;
//		}
//	}
//	tsc_file.close();
//
//#ifdef _WIN64
//	time_t end = time(NULL);
//
//	if (end - start < 60) file << "Total Time: " << fixed << setprecision(0) << end - start << " s\n";
//	else if (end - start < 3600) file << "Total Time: " << fixed << setprecision(0) << floor((end - start) / 60) << " m " << (end - start) % 60 << " s\n";
//	else file << "Total Time: " << fixed << setprecision(0) << floor((end - start) / 3600) << " h " << ((end - start) % 3600) / 60 << " m\n";
//	file << endl;
//	file << "Time Breakdown:" << endl;
//	file << " ... for Prototype Grid setup:" << setw(6) << end_prototypes - start << " s" << endl;
//	file << " ... for Becke Grid setup:    " << setw(6) << end_becke - end_prototypes << " s" << endl;
//	file << " ... for spherical density:   " << setw(6) << end_spherical - end_becke << " s" << endl;
//	file << " ... for Grid Pruning:        " << setw(6) << end_prune - end_spherical << " s" << endl;
//	file << " ... for aspherical density:  " << setw(6) << end_aspherical - end_prune << " s" << endl;
//	file << " ... for final preparation:   " << setw(6) << end1 - end_aspherical << " s" << endl;
//	file << " ... for tsc calculation:     " << setw(6) << end - end1 << " s" << endl;
//#else
//	gettimeofday(&t2, 0);
//
//	double time2 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//	if (time2 < 60) printf("Total Time: %4.1lf s\n", time2);
//	else if (time2 < 3600) printf("Total Time: %10.1lf m\n", time2 / 60);
//	else printf("Total Time: %10.1lf h\n", time2 / 3600);
//
//#endif
//
//	return true;
//}
