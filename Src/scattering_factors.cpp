/**
 * @file scattering_factors.cpp
 * @brief This file contains functions related to reading and generating hkl indices and k-points.
 *
 */

#include "tsc_block.h"
#include "scattering_factors.h"
#include "convenience.h"
#include "cell.h"
#include "wfn_class.h"
#include "spherical_density.h"
#include "AtomGrid.h"
#include "npy.h"
using namespace std;

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
void read_k_points(vector<vec> &k_pt, hkl_list &hkl, ostream &file)
{
    err_checkf(exists("kpts.dat"), "k-points file does not exist!", file);
    file << "Reading: kpts.dat" << flush;
    ifstream k_points_file("kpts.dat", ios::binary);
    err_checkf(k_points_file.good(), "Error Reading the k-points!", file);
    int nr[1]{0};
    k_points_file.read((char *)&nr, sizeof(nr));
    file << " expecting " << nr[0] << " k points... " << flush;
    double temp[1]{0.0};
    int hkl_temp[1]{0};
    k_pt.resize(3);
    ivec hkl_(3);
    for (int run = 0; run < nr[0]; run++)
    {
        for (int i = 0; i < 3; i++)
        {
            k_points_file.read((char *)&temp, sizeof(temp));
            k_pt[i].push_back(temp[0]);
            k_points_file.read((char *)&hkl_temp, sizeof(hkl_temp));
            hkl_[i] = hkl_temp[0];
        }
        hkl.emplace(hkl_);
    }
    err_checkf(!k_points_file.bad(), "Error reading k-points file!", file);
    file << " done!" << endl
         << "Size of k_points: " << k_pt[0].size() << endl;
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
void save_k_points(vector<vec> &k_pt, hkl_list &hkl)
{
    ofstream k_points_file("kpts.dat", ios::out | ios::binary | ios::trunc);
    int nr[1] = {(int)k_pt[0].size()};
    k_points_file.write((char *)&nr, sizeof(nr));
    double temp[1]{0.0};
    int hkl_temp[1]{0};
    hkl_list_it hkl_ = hkl.begin();
    for (int run = 0; run < nr[0]; run++)
    {
        for (int i = 0; i < 3; i++)
        {
            temp[0] = k_pt[i][run];
            k_points_file.write((char *)&temp, sizeof(temp));
            hkl_temp[0] = (*hkl_)[i];
            k_points_file.write((char *)&hkl_temp, sizeof(hkl_temp));
        }
        hkl_ = next(hkl_);
    }
    k_points_file.flush();
    k_points_file.close();
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
void read_hkl(const string &hkl_filename,
              hkl_list &hkl,
              const vector<vec> &twin_law,
              cell &unit_cell,
              ostream &file,
              bool debug = false)
{
    file << "Reading: " << setw(44) << hkl_filename << flush;
    ivec hkl_(3);
    err_checkf(exists(hkl_filename), "HKL file does not exists!", file);
    ifstream hkl_input(hkl_filename.c_str(), ios::in);
    hkl_input.seekg(0, hkl_input.beg);
    regex r{R"([abcdefghijklmnopqrstuvwxyz\(\)ABCDEFGHIJKLMNOPQRSTUVW])"};
    string line, temp;
    while (!hkl_input.eof())
    {
        getline(hkl_input, line);
        if (hkl_input.eof())
            break;
        if (line.size() < 2)
            continue;
        cmatch result;
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
    hkl_list_it found = hkl.find(ivec{0, 0, 0});
    if (found != hkl.end())
    {
        if (debug)
            file << "popping back 0 0 0" << endl;
        hkl.erase(ivec{0, 0, 0});
    }
    hkl_input.close();
    file << " done!\nNr of reflections read from file: " << hkl.size() << endl;

    if (debug)
        file << "Number of reflections before twin: " << hkl.size() << endl;
    if (twin_law.size() > 0)
    {
        for (const ivec &hkl__ : hkl)
            for (int i = 0; i < twin_law.size(); i++)
                hkl.emplace(ivec{
                    static_cast<int>(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2]),
                    static_cast<int>(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2]),
                    static_cast<int>(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2])});
    }
    if (debug)
        file << "Number of reflections after twin: " << hkl.size() << endl;

    vector<vector<vector<int>>> sym(3);
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
        file << "Number of symmetry operations: " << sym[0][0].size() << endl;

    ivec tempv(3);
    hkl_list hkl_enlarged = hkl;
    for (int s = 0; s < sym[0][0].size(); s++)
    {
        if (sym[0][0][s] == 1 && sym[1][1][s] == 1 && sym[2][2][s] == 1 &&
            sym[0][1][s] == 0 && sym[0][2][s] == 0 && sym[1][2][s] == 0 &&
            sym[1][0][s] == 0 && sym[2][0][s] == 0 && sym[2][1][s] == 0)
        {
            continue;
        }
        for (const ivec &hkl__ : hkl)
        {
            tempv = {0, 0, 0};
            for (int h = 0; h < 3; h++)
            {
                for (int j = 0; j < 3; j++)
                    tempv[j] += hkl__[h] * sym[j][h][s];
            }
            hkl_enlarged.emplace(tempv);
        }
    }

    for (const ivec &hkl__ : hkl_enlarged)
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
    hkl = hkl_enlarged;
    file << "Nr of reflections to be used: " << hkl.size() << endl;
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
void generate_hkl(const double &dmin,
                  hkl_list &hkl,
                  const vector<vec> &twin_law,
                  cell &unit_cell,
                  ostream &file,
                  bool debug)
{
    file << "Generating hkl indices up to d=: " << fixed << setw(17) << setprecision(2) << dmin << flush;
    ivec hkl_(3);
    string line, temp;
    const int extreme = 201;
    double dmin_l = 0.9 * dmin;
    for (int h = -extreme; h < extreme; h++)
    {
        for (int k = -extreme; k < extreme; k++)
        {
            // only need 0 to extreme, since we have no DISP signal
            for (int l = 0; l < extreme; l++)
            {
                hkl_ = {h, k, l};
                if (unit_cell.get_d_of_hkl(hkl_) >= dmin_l)
                    hkl.emplace(hkl_);
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
        for (const ivec &hkl__ : hkl)
            for (int i = 0; i < twin_law.size(); i++)
                hkl.emplace(ivec{
                    int(twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2]),
                    int(twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2]),
                    int(twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2])});
    }
    if (debug)
        file << "Number of reflections after twin: " << hkl.size() << endl;

    vector<vector<vector<int>>> sym(3);
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

    ivec tempv(3);
    hkl_list hkl_enlarged = hkl;
    for (int s = 0; s < sym[0][0].size(); s++)
    {
        if (sym[0][0][s] == 1 && sym[1][1][s] == 1 && sym[2][2][s] == 1 &&
            sym[0][1][s] == 0 && sym[0][2][s] == 0 && sym[1][2][s] == 0 &&
            sym[1][0][s] == 0 && sym[2][0][s] == 0 && sym[2][1][s] == 0)
        {
            continue;
        }
        for (const ivec &hkl__ : hkl)
        {
            tempv = {0, 0, 0};
            for (int h = 0; h < 3; h++)
            {
                for (int j = 0; j < 3; j++)
                    tempv[j] += hkl__[h] * sym[j][h][s];
            }
            hkl_enlarged.emplace(tempv);
        }
    }

    for (const ivec &hkl__ : hkl_enlarged)
    {
        tempv = hkl__;
        tempv[0] *= -1;
        tempv[1] *= -1;
        tempv[2] *= -1;
        if (hkl.find(tempv) != hkl.end() && hkl.find(hkl__) == hkl.end())
        {
            hkl.emplace(hkl__);
        }
    }
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
void generate_fractional_hkl(const double &dmin,
                             hkl_list_d &hkl,
                             const vector<vec> &twin_law,
                             cell &unit_cell,
                             ostream &file,
                             double stepsize,
                             bool debug)
{
    file << "Generating hkl indices up to d=: " << fixed << setw(17) << setprecision(2) << dmin << flush;
    vec hkl_(3);
    string line, temp;
    const int extreme = 201;
    double dmin_l = 0.9 * dmin;
    const int lim = extreme / stepsize;
#pragma omp parallel for private(hkl_)
    for (int h = -lim; h < lim; h ++)
    {
			double _h = h * stepsize;
        for (double k = -extreme; k < extreme; k += stepsize)
        {
            // only need 0 to extreme, since we have no DISP signal
            for (int l = 0; l < lim; l++)
            {
                hkl_ = {_h, k, l * stepsize};
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
        for (const vec &hkl__ : hkl)
            for (int i = 0; i < twin_law.size(); i++)
                hkl.emplace(vec{
                    (twin_law[i][0] * hkl__[0] + twin_law[i][1] * hkl__[1] + twin_law[i][2] * hkl__[2]),
                    (twin_law[i][3] * hkl__[0] + twin_law[i][4] * hkl__[1] + twin_law[i][5] * hkl__[2]),
                    (twin_law[i][6] * hkl__[0] + twin_law[i][7] * hkl__[1] + twin_law[i][8] * hkl__[2])});
    }
    if (debug)
        file << "Number of reflections after twin: " << hkl.size() << endl;

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
void read_atoms_from_CIF(ifstream &cif_input,
                         const vector<int> &input_groups,
                         const cell &unit_cell,
                         WFN &wave,
                         const vector<string> &known_atoms,
                         vector<int> &atom_type_list,
                         vector<int> &asym_atom_to_type_list,
                         vector<int> &asym_atom_list,
                         vector<bool> &needs_grid,
                         ostream &file,
                         const bool debug)
{
    bool atoms_read = false;
    int count_fields = 0;
    int group_field = 0;
    int type_field = 0;
    int position_field[3] = {0, 0, 0};
    int label_field = 1000;
    string line;
    cif_input.clear();
    cif_input.seekg(0, cif_input.beg);
    if (debug && input_groups.size() > 0)
        file << "Group size: " << input_groups.size() << endl;
    while (!cif_input.eof() && !atoms_read)
    {
        getline(cif_input, line);
        if (debug)
            file << "line: " << line << endl;
        if (line.find("loop_") != string::npos)
        {
            getline(cif_input, line);
            if (debug)
                file << "line in loop field definition: " << trim(line) << endl;
            while (trim(line).find("_") == 0)
            {
                if (debug)
                    file << "line in loop field definition: " << trim(line) << endl;
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
            while (trim(line).find("_") > 0 && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                vector<string> fields;
                fields.resize(count_fields);
                int nr = -1;
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                fields[label_field].erase(remove_if(fields[label_field].begin(), fields[label_field].end(), ::isspace), fields[label_field].end());
                fields[type_field].erase(remove_if(fields[type_field].begin(), fields[type_field].end(), ::isspace), fields[type_field].end());
                if (debug)
                    file << "label: " << setw(8) << fields[label_field] << " type: " << fields[type_field] << " frac. pos: "
                         << fixed << setprecision(3) << stod(fields[position_field[0]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[0]]) << " "
                         << fixed << setprecision(3) << stod(fields[position_field[1]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[1]]) << " "
                         << fixed << setprecision(3) << stod(fields[position_field[2]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[2]]) << " " << flush;
                vector<double> position = unit_cell.get_coords_cartesian(
                    stod(fields[position_field[0]]),
                    stod(fields[position_field[1]]),
                    stod(fields[position_field[2]]));
                vector<double> precisions = unit_cell.get_coords_cartesian(
                    get_decimal_precision_from_CIF_number(fields[position_field[0]]),
                    get_decimal_precision_from_CIF_number(fields[position_field[1]]),
                    get_decimal_precision_from_CIF_number(fields[position_field[2]]));
                for (int i = 0; i < 3; i++)
                {
                    precisions[i] = abs(precisions[i]);
                }
                if (debug)
                    file << " cart. pos.: " << setw(8) << position[0] << "+/-" << precisions[0] << " " << setw(8) << position[1] << "+/-" << precisions[1] << " " << setw(8) << position[2] << "+/-" << precisions[2] << endl;
                bool old_atom = false;
#pragma omp parallel for reduction(|| : old_atom)
                for (int run = 0; run < known_atoms.size(); run++)
                {
                    if (fields[label_field] == known_atoms[run])
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
                        tolerances[j] = 2 * max(min(abs(precisions[j]), 1.0), 0.01);
                    }
                    if (is_similar_abs(position[0], wave.atoms[i].x, tolerances[0]) && is_similar_abs(position[1], wave.atoms[i].y, tolerances[1]) && is_similar_abs(position[2], wave.atoms[i].z, tolerances[2]))
                    {
                        string element = atnr2letter(wave.get_atom_charge(i));
                        err_checkf(element != "PROBLEM", "Problem identifying atoms!", std::cout);
                        string label = fields[label_field];
                        string type = fields[type_field];
                        transform(element.begin(), element.end(), element.begin(), asciitolower);
                        transform(label.begin(), label.end(), label.begin(), asciitolower);
                        transform(type.begin(), type.end(), type.begin(), asciitolower);
                        if (debug)
                        {
                            file << "ASYM:  " << setw(8) << element << " charge: " << setw(17) << wave.get_atom_charge(i) << "                          wfn cart. pos: "
                                 << fixed << setprecision(3) << setw(16) << wave.atoms[i].x << " "
                                 << fixed << setprecision(3) << setw(16) << wave.atoms[i].y << " "
                                 << fixed << setprecision(3) << setw(16) << wave.atoms[i].z << flush;
                            if (input_groups.size() > 0)
                            {
                                file << " checking disorder group: " << fields[group_field] << " vs. ";
                                for (int g = 0; g < input_groups.size(); g++)
                                    file << input_groups[g] << ",";
                            }
                        }
                        if (input_groups.size() > 0)
                        {
                            bool yep = false;
                            for (int g = 0; g < input_groups.size(); g++)
                            {
                                if (fields[group_field].c_str()[0] == '.' && input_groups[g] == 0)
                                {
                                    if (debug)
                                        file << "appears to be group 0" << endl;
                                    yep = true;
                                    break;
                                }
                                else if (stoi(fields[group_field]) == input_groups[g])
                                    yep = true;
                            }
                            if (!yep)
                            {
                                if (debug)
                                    file << "Wrong part!" << endl;
                                continue;
                            }
                        }
                        if (label.find(element) == string::npos)
                        {
                            if (element != "h")
                            {
                                if (debug)
                                {
                                    file << "\nElement symbol not found in label, this is a problem!\n checking type...";
                                    if (type.find(element) == string::npos)
                                    {
                                        file << " ALSO FAILED! WILL IGNORE ATOM!\n";
                                        continue;
                                    }
                                }
                                else
                                {
                                    if (type.find(element) == string::npos)
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
                                    if (type.find(element) == string::npos)
                                    {
                                        file << " ALSO FAILED! WILL IGNORE ATOM!\n";
                                        continue;
                                    }
                                }
                                else
                                {
                                    if (type.find(element) == string::npos)
                                    {
                                        file << "\nAtom " << label << " was not matching by element determined by label reduction or type field, skipping!\n";
                                        continue;
                                    }
                                }
                            }
                        }
                        wave.atoms[i].label = fields[label_field];
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
    file << "... done!" << endl;
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
}

/**
 * Generates Hirshfeld grids based on the specified parameters.
 *
 * @param pbc The periodic boundary condition flag.
 * @param accuracy The accuracy level for grid generation.
 * @param unit_cell The cell object representing the unit cell.
 * @param wave The WFN object containing wavefunction information.
 * @param atom_type_list The list of atom types.
 * @param asym_atom_list The list of asymmetric atoms.
 * @param needs_grid The vector indicating whether each atom needs a grid.
 * @param d1 The vector of grid vectors in the x-direction.
 * @param d2 The vector of grid vectors in the y-direction.
 * @param d3 The vector of grid vectors in the z-direction.
 * @param dens The vector of density vectors.
 * @param file The output stream for writing the grid data.
 * @param start The start time point for timing measurements.
 * @param end_becke The end time point for the Becke grid generation.
 * @param end_prototypes The end time point for the prototype grid generation.
 * @param end_spherical The end time point for the spherical grid generation.
 * @param end_prune The end time point for the pruning step.
 * @param end_aspherical The end time point for the aspherical grid generation.
 * @param debug Flag indicating whether to enable debug mode.
 * @param no_date Flag indicating whether to exclude the date in the output.
 *
 * @return The number of grid points in the final total grid.
 */
int make_hirshfeld_grids(const int &pbc,
                         const int &accuracy,
                         cell &unit_cell,
                         const WFN &wave,
                         const ivec &atom_type_list,
                         const ivec &asym_atom_list,
                         vector<bool> &needs_grid,
                         vector<vec> &d1,
                         vector<vec> &d2,
                         vector<vec> &d3,
                         vector<vec> &dens,
                         ostream &file,
                         time_point &start,
                         time_point &end_becke,
                         time_point &end_prototypes,
                         time_point &end_spherical,
                         time_point &end_prune,
                         time_point &end_aspherical,
                         bool debug,
                         bool no_date)
{
    int atoms_with_grids = 0;
    for (int i = 0; i < needs_grid.size(); i++)
    {
        if (needs_grid[i])
            atoms_with_grids++;
    }
    ivec num_points(atoms_with_grids);
    vector<vector<vec>> grid(atoms_with_grids);
    const int nr_of_atoms = (wave.get_ncen() * (int)pow(pbc * 2 + 1, 3));
    vec x(nr_of_atoms), y(nr_of_atoms), z(nr_of_atoms);
    ivec atom_z(nr_of_atoms);
    vec alpha_max(wave.get_ncen());
    ivec max_l(wave.get_ncen());
    int max_l_overall = 0;
#pragma omp parallel for
    for (int i = 0; i < atoms_with_grids; i++)
        grid[i].resize(6);
        // GRID COORDINATES for [a][c][p] a = atom [0,ncen],
        // c = coordinate [0=x, 1=y, 2=z, 3=atomic becke weight, 4= molecular becke weight, 5=total spherical density],
        // p = point in this grid

        // Accumulate vectors with information about all atoms
#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        atom_z[i] = wave.get_atom_charge(i);
        x[i] = wave.atoms[i].x;
        y[i] = wave.atoms[i].y;
        z[i] = wave.atoms[i].z;
        // if(debug)
        //     file << "xyz= 000 position: " << x[i] << " " << y[i] << " " << z[i] << " Charge: " << atom_z[i] << endl;
        if (pbc != 0)
        {
            int j = 0;
            for (int pbc_x = -pbc; pbc_x < pbc + 1; pbc_x++)
                for (int pbc_y = -pbc; pbc_y < pbc + 1; pbc_y++)
                    for (int pbc_z = -pbc; pbc_z < pbc + 1; pbc_z++)
                    {
                        if (pbc_x == 0 && pbc_y == 0 && pbc_z == 0)
                            continue;
                        else
                        {
                            j++;
                            atom_z[i + j * wave.get_ncen()] = wave.get_atom_charge(i);
                            x[i + j * wave.get_ncen()] = wave.atoms[i].x + pbc_x * unit_cell.get_cm(0, 0) + pbc_y * unit_cell.get_cm(0, 1) + pbc_z * unit_cell.get_cm(0, 2);
                            y[i + j * wave.get_ncen()] = wave.atoms[i].y + pbc_x * unit_cell.get_cm(1, 0) + pbc_y * unit_cell.get_cm(1, 1) + pbc_z * unit_cell.get_cm(1, 2);
                            z[i + j * wave.get_ncen()] = wave.atoms[i].z + pbc_x * unit_cell.get_cm(2, 0) + pbc_y * unit_cell.get_cm(2, 1) + pbc_z * unit_cell.get_cm(2, 2);
                            if (debug)
                                file << "xyz= " << pbc_x << pbc_y << pbc_z << " j = " << j << " position: " << x[i + j * wave.get_ncen()] << " " << y[i + j * wave.get_ncen()] << " " << z[i + j * wave.get_ncen()] << " Charge: " << atom_z[i + j * wave.get_ncen()] << endl;
                        }
                    }
        }
        alpha_max[i] = 0.0;
        max_l[i] = 0;
        for (int b = 0; b < wave.get_nex(); b++)
        {
            if (wave.get_center(b) != i + 1)
                continue;
            if (wave.get_exponent(b) > alpha_max[i])
                alpha_max[i] = wave.get_exponent(b);
            if (wave.get_type(b) > max_l[i])
            {
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
#pragma omp critical
                {
                    if (l > max_l_overall)
                        max_l_overall = l;
                }
            }
        }
    }

    if (debug)
    {
        file << "Atoms are there! max_l:" << setw(5) << max_l_overall << endl;
        for (int i = 0; i < max_l.size(); i++)
            file << "max_l: " << setw(5) << max_l[i] << endl;
    }

    vector<vec> alpha_min(wave.get_ncen());
    for (int i = 0; i < wave.get_ncen(); i++)
        alpha_min[i].resize(max_l_overall, 100000000.0);

#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        for (int b = 0; b < max_l_overall; b++)
            alpha_min[i][b] = 100000000.0;
    }

#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        for (int b = 0; b < wave.get_nex(); b++)
        {
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
    }

    if (debug)
    {
        for (int i = 0; i < wave.get_ncen(); i++)
        {
            file << "alpha_min: ";
            for (int b = 0; b < max_l_overall; b++)
                file << setw(14) << scientific << alpha_min[i][b];
            file << endl;
        }
    }

    if (debug)
        file << "alpha_min is there!" << endl
             << "Nr of asym atoms: " << asym_atom_list.size() << " Number of atoms in wfn: " << wave.get_ncen() << " atoms that needs a grid: " << atoms_with_grids << endl;
    else
        file << "There are:\n"
             << setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
             //<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
             << setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << endl;
#ifdef FLO_CUDA
    int nDevices;                  /*Number of devices available (running time)*/
    cudaGetDeviceCount(&nDevices); /*Get the number of devices*/
    int dev = 0;
    cudaDeviceProp *prop = NULL;
    cudaDeviceProp deviceProp;
    prop = (cudaDeviceProp *)malloc(sizeof(cudaDeviceProp) * nDevices);
    for (int devl = 0; devl < nDevices; devl++)
    { // Make CUDA information available in prop
        cudaGetDeviceProperties(&(deviceProp), devl);
        prop[devl] = deviceProp;
    }
    cudaSetDevice(dev);
    printCUDA(prop, nDevices, file);
    gpuErrchk(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
#endif

    if (no_date)
        file << "\nMaking Becke Grids..." << flush;
    else
    {
        if (debug)
            file << "max_l_overall: " << max_l_overall << endl
                 << "Selected accuracy: " << accuracy << "\nMaking Becke Grid for" << endl;
        else
            file << endl
                 << "Selected accuracy: " << accuracy << "\nMaking Becke Grids..." << flush;
    }

    // Make Prototype grids with only single atom weights for all elements
    vector<AtomGrid> Prototype_grids;

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        if (debug)
            file << "Atom Type " << i << ": " << atom_type_list[i] << endl;
        double alpha_max_temp(0);
        int max_l_temp(0);
        vec alpha_min_temp(max_l_overall);
        for (int j = 0; j < wave.get_ncen(); j++)
        {
            if (wave.get_atom_charge(j) == 119)
            {
                continue;
            }
            if (wave.get_atom_charge(j) == atom_type_list[i])
            {
                if (debug)
                {
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

        if (debug)
        {
            file << "max_l: " << defaultfloat << max_l_temp << " alpha_max: " << scientific << alpha_max_temp << " alpha_min: ";
            for (int l = 0; l <= max_l_temp; l++)
                file << setw(14) << scientific << alpha_min_temp[l];
            file << " accuracy: " << accuracy << endl;
        }
        int lebedev_high, lebedev_low;
        double radial_acc;
        err_checkf(accuracy >= 0, "Negative accuracy is not defined!", file);
        if (accuracy == 0)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                radial_acc = 1e-4;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                radial_acc = 1e-3;
            }
        }
        else if (accuracy == 1)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[7] : constants::lebedev_table[8];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[3] : constants::lebedev_table[4];
                radial_acc = 1e-4;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[6] : constants::lebedev_table[7];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[2] : constants::lebedev_table[3];
                radial_acc = 1e-5;
            }
        }
        else if (accuracy == 2)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[11] : constants::lebedev_table[12];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[7] : constants::lebedev_table[8];
                radial_acc = 1e-5;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[10] : constants::lebedev_table[11];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[6] : constants::lebedev_table[7];
                radial_acc = 1e-6;
            }
        }
        else if (accuracy == 3)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[14] : constants::lebedev_table[16];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[12] : constants::lebedev_table[14];
                radial_acc = 1e-10;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[13] : constants::lebedev_table[15];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[11] : constants::lebedev_table[13];
                radial_acc = 1e-11;
            }
        }
        else if (accuracy == 4)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[19] : constants::lebedev_table[21];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[14] : constants::lebedev_table[17];
                radial_acc = 1e-20;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[18] : constants::lebedev_table[20];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[13] : constants::lebedev_table[16];
                radial_acc = 1e-15;
            }
        }
        else
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[31] : constants::lebedev_table[32];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[29] : constants::lebedev_table[31];
                radial_acc = 1e-20;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[30] : constants::lebedev_table[32];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[28] : constants::lebedev_table[30];
                radial_acc = 1e-15;
            }
        }
        Prototype_grids.push_back(AtomGrid(radial_acc,
                                           lebedev_low,
                                           lebedev_high,
                                           atom_type_list[i],
                                           alpha_max_temp,
                                           max_l_temp,
                                           alpha_min_temp.data(),
                                           file));
    }

    end_prototypes = get_time();
    if (debug)
    {

        for (int prototype = 0; prototype < Prototype_grids.size(); prototype++)
            file << "Number of gridpoints for atom type " << atom_type_list[prototype] << ": " << Prototype_grids[prototype].get_num_grid_points() << endl;

        int dur = get_sec(start, end_prototypes);

        if (dur < 1)
            file << "Time until prototypes are done: " << fixed << setprecision(0) << get_msec(start, end_prototypes) << " ms" << endl;
        else
            file << "Time until prototypes are done: " << fixed << setprecision(0) << dur << " s" << endl;
    }

#ifdef FLO_CUDA
    vector<vector<double>> radial_density;
    vector<vector<double>> radial_dist;

    radial_density.resize(atom_type_list.size());
    radial_dist.resize(atom_type_list.size());
    spherical_density.resize(asym_atom_list.size());
    for (int i = 0; i < asym_atom_list.size(); i++)
        spherical_density[i].resize(num_points[i]);

    const double incr = pow(1.005, max(1, accuracy - 1));
    const double lincr = log(incr);
    const double min_dist = 0.0000001;
    vector<Thakkar> sphericals;
    for (int i = 0; i < atom_type_list.size(); i++)
        sphericals.push_back(Thakkar(atom_type_list[i]));
    // Make radial grids
    if (debug)
    {
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            file << "Calculating for atomic number " << atom_type_list[i] << endl;
            double current = 1;
            double dist = min_dist;
            if (accuracy > 3)
                while (current > 1E-10)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    if (current == -20)
                        return false;
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    if (current == -20)
                        return false;
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            file << "Number of radial density points for atomic number " << atom_type_list[i] << ": " << radial_density[i].size() << endl;
            // for (int j = 0; j < radial_density[i].size(); j++) {
            //	if (radial_density[i][j] < 0.1)
            //		break;
            //	file << scientific << setprecision(8) << radial_density[i][j] << endl;
            // }
        }
    }
    else
    {
#pragma omp parallel for
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            double current = 1;
            double dist = min_dist;
            if (accuracy > 3)
                while (current > 1E-10)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
        }
    }
    sphericals.clear();

    float **gpu_PosAtomsx = NULL,
          **gpu_PosAtomsy = NULL,
          **gpu_PosAtomsz = NULL,
          **gpu_GridRho = NULL,
          **gpu_Gridx = NULL,
          **gpu_Gridy = NULL,
          **gpu_Gridz = NULL,
          **gpu_Gridaw = NULL,
          **gpu_Gridmw = NULL,
          **gpu_exponents = NULL,
          **gpu_coefficients = NULL,
          **gpu_occ = NULL;
    double ***gpu_atomgrid_x = NULL,
           ***gpu_atomgrid_y = NULL,
           ***gpu_atomgrid_z = NULL,
           ***gpu_atomgrid_w = NULL;
    int **gpu_types = NULL,
        **gpu_centers = NULL,
        **gpu_asym_atom_list = NULL,
        **gpu_atom_type_list = NULL,
        **gpu_numpoints = NULL,
        **gpu_atom_z = NULL;
    vector<vector<float>> PosAtoms;
    PosAtoms.resize(3);
    for (int i = 0; i < 3; i++)
        PosAtoms[i].resize(wave.get_ncen());
    for (int a = 0; a < wave.get_ncen(); a++)
    {
        PosAtoms[0][a] = wave.atoms[a].x;
        PosAtoms[1][a] = wave.atoms[a].y;
        PosAtoms[2][a] = wave.atoms[a].z;
    }
    /*Allocation GPU Pointer*/
    gpu_PosAtomsx = (float **)malloc(sizeof(float *));
    gpu_PosAtomsy = (float **)malloc(sizeof(float *));
    gpu_PosAtomsz = (float **)malloc(sizeof(float *));
    gpu_atom_z = (int **)malloc(sizeof(int *));
    gpu_types = (int **)malloc(sizeof(int *));
    gpu_centers = (int **)malloc(sizeof(int *));
    gpu_asym_atom_list = (int **)malloc(sizeof(int *));
    gpu_atom_type_list = (int **)malloc(sizeof(int *));
    gpu_numpoints = (int **)malloc(sizeof(int *));
    gpu_exponents = (float **)malloc(sizeof(float *));
    gpu_occ = (float **)malloc(sizeof(float *));
    gpu_coefficients = (float **)malloc(sizeof(float *));
    gpu_GridRho = (float **)malloc(sizeof(float *));
    gpu_Gridx = (float **)malloc(sizeof(float *));
    gpu_Gridy = (float **)malloc(sizeof(float *));
    gpu_Gridz = (float **)malloc(sizeof(float *));
    gpu_Gridaw = (float **)malloc(sizeof(float *));
    gpu_Gridmw = (float **)malloc(sizeof(float *));
    gpu_atomgrid_x = (double ***)malloc(sizeof(double **));
    gpu_atomgrid_y = (double ***)malloc(sizeof(double **));
    gpu_atomgrid_z = (double ***)malloc(sizeof(double **));
    gpu_atomgrid_w = (double ***)malloc(sizeof(double **));
    gpu_atomgrid_x[0] = (double **)malloc(sizeof(double *) * asym_atom_list.size());
    gpu_atomgrid_y[0] = (double **)malloc(sizeof(double *) * asym_atom_list.size());
    gpu_atomgrid_z[0] = (double **)malloc(sizeof(double *) * asym_atom_list.size());
    gpu_atomgrid_w[0] = (double **)malloc(sizeof(double *) * asym_atom_list.size());

    int nex_temp = wave.get_nex();
    int nmo_temp = wave.get_nmo(true);
    int ncen_temp = wave.get_ncen();
    for (int i = 0; i < wave.get_ncen(); i++)
        atom_z[i] = wave.atoms[i].charge;
    int MaxGrid = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int nr = asym_atom_list[i];
        int type;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.atoms[nr].charge)
                type = j;

        num_points[i] = Prototype_grids[type].get_num_grid_points();
        MaxGrid += num_points[i];
    }
    int numBlocks, blocks, gridSize;
    size_t size;
    gpuErrchk(cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize));
    if (debug)
        file << "\nold Heapsize: " << size / 1024 / 1024 << " MB" << endl;
    float result;
    gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
        &numBlocks,
        &blocks,
        (void *)gpu_calc_dens_per_MO_static1,
        0,
        MaxGrid));
    gridSize = (MaxGrid + blocks - 1) / blocks;
    result = ((sizeof(int) * 10 + sizeof(float) * (6 * ncen_temp + 20)) * blocks * prop[dev].multiProcessorCount) / 1024 / 1024;

    if (debug)
        file << "result: " << fixed << result << " MB for " << gridSize << " blocks with " << blocks << " threads" << endl
             << "sizeof float: " << sizeof(float) << " B" << endl;
    if (result > size / 1024 / 1024)
        gpuErrchk(cudaDeviceSetLimit(cudaLimitMallocHeapSize, result * 1024 * 1024 * 1.1));

    gpuErrchk(cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize));
    if (debug)
        file << "new Heapsize: " << size / 1024 / 1024 << " MB" << endl;

    gridSize = (MaxGrid + blocks - 1) / blocks;

    // Allocate and copy vectors to device for WFN
    if (debug)
        file << "Copying WFN to devices now!" << endl;
    gpuErrchk(cudaMalloc((void **)&gpu_PosAtomsx[0], sizeof(float) * ncen_temp));
    gpuErrchk(cudaMalloc((void **)&gpu_PosAtomsy[0], sizeof(float) * ncen_temp));
    gpuErrchk(cudaMalloc((void **)&gpu_PosAtomsz[0], sizeof(float) * ncen_temp));
    gpuErrchk(cudaMalloc((void **)&gpu_atom_z[0], sizeof(int) * ncen_temp));
    gpuErrchk(cudaMalloc((void **)&gpu_Gridx[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void **)&gpu_Gridy[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void **)&gpu_Gridz[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void **)&gpu_Gridaw[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void **)&gpu_Gridmw[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void **)&gpu_asym_atom_list[0], sizeof(int) * asym_atom_list.size()));
    gpuErrchk(cudaMalloc((void **)&gpu_atom_type_list[0], sizeof(int) * atom_type_list.size()));
    gpuErrchk(cudaMalloc((void **)&gpu_numpoints[0], sizeof(int) * asym_atom_list.size()));
    if (debug)
        file << "Mallocs done!" << endl;
    gpuErrchk(cudaMemcpyToSymbol(gpu_nex, &nex_temp, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_nmo, &nmo_temp, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_ncen, &ncen_temp, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_MaxGrid, &MaxGrid, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_start_radial_dens, &min_dist, sizeof(float)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_log_incr, &lincr, sizeof(float)));
    gpuErrchk(cudaMemcpy(gpu_PosAtomsx[0], PosAtoms[0].data(), sizeof(float) * wave.get_ncen(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_PosAtomsy[0], PosAtoms[1].data(), sizeof(float) * wave.get_ncen(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_PosAtomsz[0], PosAtoms[2].data(), sizeof(float) * wave.get_ncen(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_asym_atom_list[0], asym_atom_list.data(), sizeof(int) * asym_atom_list.size(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_atom_type_list[0], atom_type_list.data(), sizeof(int) * atom_type_list.size(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_numpoints[0], num_points.data(), sizeof(int) * asym_atom_list.size(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaPeekAtLastError());
    file << "All copying done!" << endl;

    vector<cudaStream_t> streams;
    streams.resize(asym_atom_list.size());
    int offset = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
            &numBlocks,
            &blocks,
            (void *)gpu_make_grid,
            0,
            num_points[i]));

        gridSize = (MaxGrid + blocks - 1) / blocks;
        if (debug)
            file << i << ": num points: " << num_points[i] << " blocks: " << gridSize << " threads: " << blocks << " grid > points? " << (blocks * gridSize > num_points[i]) << endl;
        int type;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.atoms[asym_atom_list[i]].charge)
                type = j;
        gpuErrchk(cudaStreamCreate(&streams[i]));
        gpuErrchk(cudaMalloc((void **)&gpu_atomgrid_x[0][i], sizeof(double) * num_points[i]));
        gpuErrchk(cudaMalloc((void **)&gpu_atomgrid_y[0][i], sizeof(double) * num_points[i]));
        gpuErrchk(cudaMalloc((void **)&gpu_atomgrid_z[0][i], sizeof(double) * num_points[i]));
        gpuErrchk(cudaMalloc((void **)&gpu_atomgrid_w[0][i], sizeof(double) * num_points[i]));
        gpuErrchk(cudaMemcpyAsync(gpu_atomgrid_x[0][i], Prototype_grids[type].get_gridx_ptr(), sizeof(double) * num_points[i], cudaMemcpyHostToDevice, streams[i]));
        gpuErrchk(cudaMemcpyAsync(gpu_atomgrid_y[0][i], Prototype_grids[type].get_gridy_ptr(), sizeof(double) * num_points[i], cudaMemcpyHostToDevice, streams[i]));
        gpuErrchk(cudaMemcpyAsync(gpu_atomgrid_z[0][i], Prototype_grids[type].get_gridz_ptr(), sizeof(double) * num_points[i], cudaMemcpyHostToDevice, streams[i]));
        gpuErrchk(cudaMemcpyAsync(gpu_atomgrid_w[0][i], Prototype_grids[type].get_gridw_ptr(), sizeof(double) * num_points[i], cudaMemcpyHostToDevice, streams[i]));

        gpu_make_grid<<<gridSize, blocks, 0, streams[i]>>>(
            i,
            gpu_PosAtomsx[0],
            gpu_PosAtomsy[0],
            gpu_PosAtomsz[0],
            gpu_atomgrid_x[0][i],
            gpu_atomgrid_y[0][i],
            gpu_atomgrid_z[0][i],
            gpu_atomgrid_w[0][i],
            gpu_atom_z[0],
            gpu_asym_atom_list[0],
            gpu_numpoints[0],
            offset,
            gpu_Gridx[0],
            gpu_Gridy[0],
            gpu_Gridz[0],
            gpu_Gridaw[0],
            gpu_Gridmw[0]);

        offset += num_points[i];
    }
    gpuErrchk(cudaDeviceSynchronize());
    gpuErrchk(cudaPeekAtLastError());

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        gpuErrchk(cudaFree(gpu_atomgrid_w[0][i]));
        gpuErrchk(cudaFree(gpu_atomgrid_x[0][i]));
        gpuErrchk(cudaFree(gpu_atomgrid_y[0][i]));
        gpuErrchk(cudaFree(gpu_atomgrid_z[0][i]));
    }
#else
    vector<vector<double>> radial_density(atom_type_list.size());
    vector<vector<double>> radial_dist(atom_type_list.size());
    if (!debug)
    {
        file << " ...  " << flush;
    }
    // get_grid is parallelized, therefore not parallel here
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        // skip atoms, that do not need a grid
        if (!needs_grid[i])
        {
            continue;
        }
        if (debug)
        {
            file << "Making grid for atom " << i << endl;
        }
        int type = 0;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.get_atom_charge(i))
                type = j;

        int grid_number = 0;
        for (int j = 0; j < i; j++)
        {
            if (needs_grid[j])
            {
                grid_number++;
            }
        }
        num_points[grid_number] = Prototype_grids[type].get_num_grid_points();

        for (int n = 0; n < 6; n++)
            grid[grid_number][n].resize(num_points[grid_number], 0.0);

        Prototype_grids[type].get_grid(int(wave.get_ncen() * pow(pbc * 2 + 1, 3)),
                                       i,
                                       &x[0],
                                       &y[0],
                                       &z[0],
                                       &atom_z[0],
                                       grid[grid_number][0].data(),
                                       grid[grid_number][1].data(),
                                       grid[grid_number][2].data(),
                                       grid[grid_number][3].data(),
                                       grid[grid_number][5].data());
    }
    if (debug)
    {
        int run = 0;
        file << "  label | needs_grid | number of gridpoints\n";
        for (int j = 0; j < wave.get_ncen(); j++)
        {
            file << setw(8) << wave.atoms[j].label << setw(13) << needs_grid[j];
            if (needs_grid[j])
            {
                file << setw(7) << num_points[run];
                run++;
            }
            else
            {
                file << setw(6) << "---";
            }
            file << endl;
        }
    }
    Prototype_grids.clear();
#endif

    int points = 0;
    for (int i = 0; i < atoms_with_grids; i++)
        points += num_points[i];
    if (debug)
        file << "Becke Grid exists" << endl;
    else
        file << "                           done! Number of gridpoints: " << defaultfloat << points << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_becke = get_time();

    file << "Calculating spherical densities..." << flush;

    // Total grid as a sum of all atomic grids.
    // Dimensions: [c] [p]
    // p = the number of gridpoint
    // c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
#ifdef FLO_CUDA
    vector<vector<float>> total_grid(7);
    float ***gpu_spherical_density = NULL,
          ***gpu_radial_density = NULL,
          ***gpu_radial_dist = NULL,
          **gpu_Grids = NULL;
    gpu_radial_density = (float ***)malloc(sizeof(float **));
    gpu_radial_dist = (float ***)malloc(sizeof(float **));
    gpu_spherical_density = (float ***)malloc(sizeof(float **));

    gpu_radial_density[0] = (float **)malloc(sizeof(float *) * atom_type_list.size());
    gpu_radial_dist[0] = (float **)malloc(sizeof(float *) * atom_type_list.size());
    gpu_spherical_density[0] = (float **)malloc(sizeof(float *) * asym_atom_list.size());
    gpu_Grids = (float **)malloc(sizeof(float *));

    for (int i = 0; i < asym_atom_list.size(); i++)
        gpuErrchk(cudaMalloc((void **)&(gpu_spherical_density[0][i]), sizeof(float) * num_points[i]));

    gpuErrchk(cudaMalloc((void **)&gpu_Grids[0], sizeof(float) * MaxGrid));

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        gpuErrchk(cudaMalloc((void **)&gpu_radial_density[0][i], sizeof(float) * radial_density[i].size()));
        gpuErrchk(cudaMalloc((void **)&gpu_radial_dist[0][i], sizeof(float) * radial_dist[i].size()));
        gpuErrchk(cudaMemcpy(gpu_radial_density[0][i], radial_density[i].data(), sizeof(float) * radial_density[i].size(), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(gpu_radial_dist[0][i], radial_dist[i].data(), sizeof(float) * radial_dist[i].size(), cudaMemcpyHostToDevice));
    }
    gpuErrchk(cudaDeviceSynchronize());
    gpuErrchk(cudaPeekAtLastError());
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        int nr = all_atom_list[i];
        int type_list_number = -1;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (wave.atoms[nr].charge == atom_type_list[j])
                type_list_number = j;
        offset = 0;
        for (int g = 0; g < asym_atom_list.size(); g++)
        {
            if (g == 0)
            {
                gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
                    &numBlocks,
                    &blocks,
                    (void *)gpu_linear_interpolate_spherical_density,
                    0,
                    num_points[i]));

                gridSize = (MaxGrid + blocks - 1) / blocks;
                // file << i << "/" << g << ": blocks: " << gridSize << " threads: " << blocks << endl;
            }
            bool match = (all_atom_list[i] == asym_atom_list[g]);
            gpu_linear_interpolate_spherical_density<<<gridSize, blocks>>>(
                i,
                gpu_radial_density[0][type_list_number],
                gpu_radial_dist[0][type_list_number],
                radial_density[type_list_number].size(),
                match,
                offset,
                gpu_Gridx[0],
                gpu_Gridy[0],
                gpu_Gridz[0],
                gpu_numpoints[0],
                gpu_spherical_density[0][g],
                gpu_Grids[0],
                gpu_PosAtomsx[0],
                gpu_PosAtomsy[0],
                gpu_PosAtomsz[0]);

            offset += num_points[i];
        }
    }
    gpuErrchk(cudaDeviceSynchronize());
    gpuErrchk(cudaPeekAtLastError());

    if (debug)
    {
        // Copy grid from GPU to print:
        //  Dimensions: [c] [p]
        //  p = the number of gridpoint
        //  c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
        for (int i = 0; i < 7; i++)
            total_grid[i].resize(MaxGrid);
        gpuErrchk(cudaMemcpy(total_grid[0].data(), gpu_Gridx[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[1].data(), gpu_Gridy[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[2].data(), gpu_Gridz[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[3].data(), gpu_Gridaw[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[6].data(), gpu_Gridmw[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[4].data(), gpu_Grids[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));

        gpuErrchk(cudaDeviceSynchronize());
        gpuErrchk(cudaPeekAtLastError());

        ofstream grid0("grid0.file", ios::out);
        for (int i = 0; i < MaxGrid; i++)
        {
            for (int j = 0; j < 7; j++)
                grid0 << setw(16) << scientific << setprecision(8) << total_grid[j][i];
            grid0 << "\n";
        }
        grid0.flush();
        grid0.close();
    }

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        gpuErrchk(cudaFree(gpu_radial_density[0][i]));
        gpuErrchk(cudaFree(gpu_radial_dist[0][i]));
    }

#else
    vector<vector<double>> total_grid(7);
    // density of spherical atom at each
    // Dimensions: [a] [d]
    // a = atom number in atom type list for which the weight is calcualted
    // d = distance to look at obtained from point_to_distance_map
    vector<vector<double>> spherical_density(atoms_with_grids);

    const double incr = pow(1.005, max(1, accuracy - 1));
    const double lincr = log(incr);
    const double min_dist = 0.0000001;
    vector<Thakkar> sphericals;
    for (int i = 0; i < atom_type_list.size(); i++)
        sphericals.push_back(Thakkar(atom_type_list[i]));
    // Make radial grids
    if (debug)
    {
        file << "\nSize of atom_type_list:" << setw(5) << atom_type_list.size() << endl;
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            file << "\nCalculating for atomic number " << atom_type_list[i] << endl;
            double current = 1;
            double dist = min_dist;
            if (accuracy < 3)
                while (current > 1E-10)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    if (current == -20)
                        return false;
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    if (current == -20)
                        return false;
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            file << "Number of radial density points for atomic number " << atom_type_list[i] << ": " << radial_density[i].size() << endl;
        }
    }
    else
    {
#pragma omp parallel for
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            double current = 1;
            double dist = min_dist;
            if (accuracy < 3)
                while (current > 1E-10)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
        }
    }
    // sphericals.clear();
    int type_list_number = -1;
    if (debug)
    {
        file << "Cleared the sphericals!" << endl;
    }
#pragma omp parallel
    {
#pragma omp for
        for (int g = 0; g < atoms_with_grids; g++)
        {
            spherical_density[g].resize(num_points[g]);
        }
        for (int i = 0; i < wave.get_ncen(); i++)
        {
#pragma omp single
            {
                type_list_number = -1;
                // Determine which type in the type list of sphericals to use
                for (int j = 0; j < atom_type_list.size(); j++)
                    if (wave.get_atom_charge(i) == atom_type_list[j])
                        type_list_number = j;
                if (debug && type_list_number != -1)
                {
                    file << type_list_number << " Atom type: " << atom_type_list[type_list_number] << endl;
                }
            }
            if (type_list_number == -1)
            {
#pragma omp single
                {
                    file << "I skipped an atom! make sure this is okay!" << endl;
                }
                continue;
            }
#pragma omp for
            for (int g = 0; g < atoms_with_grids; g++)
            {
                for (int p = 0; p < num_points[g]; p++)
                {
                    double temp =
                        linear_interpolate_spherical_density(radial_density[type_list_number], radial_dist[type_list_number], sqrt(pow(grid[g][0][p] - wave.atoms[i].x, 2) + pow(grid[g][1][p] - wave.atoms[i].y, 2) + pow(grid[g][2][p] - wave.atoms[i].z, 2)), lincr, min_dist);
                    if (i == asym_atom_list[g])
                        spherical_density[g][p] = temp;
                    grid[g][4][p] += temp;
                }
            }
        }
    }
    // fill out with priodic information
    if (pbc != 0)
    {
        for (int _x = -pbc; _x < pbc + 1; _x++)
            for (int _y = -pbc; _y < pbc + 1; _y++)
                for (int _z = -pbc; _z < pbc + 1; _z++)
                {
                    if (_x == 0 && _y == 0 && _z == 0)
                        continue;
                    for (int i = 0; i < wave.get_ncen(); i++)
                    {
                        type_list_number = -1;
                        // int nr = all_atom_list[i];
                        for (int j = 0; j < atom_type_list.size(); j++)
                            if (wave.get_atom_charge(i) == atom_type_list[j])
                                type_list_number = j;
                        for (int g = 0; g < atoms_with_grids; g++)
                        {
#pragma omp parallel for
                            for (int p = 0; p < num_points[g]; p++)
                                grid[g][4][p] += linear_interpolate_spherical_density(
                                    radial_density[type_list_number],
                                    radial_dist[type_list_number],
                                    sqrt(
                                        pow(grid[g][0][p] - (wave.atoms[i].x + _x * unit_cell.get_cm(0, 0) + _y * unit_cell.get_cm(0, 1) + _z * unit_cell.get_cm(0, 2)), 2) + pow(grid[g][1][p] - (wave.atoms[i].y + _x * unit_cell.get_cm(1, 0) + _y * unit_cell.get_cm(1, 1) + _z * unit_cell.get_cm(1, 2)), 2) + pow(grid[g][2][p] - (wave.atoms[i].z + _x * unit_cell.get_cm(2, 0) + _y * unit_cell.get_cm(2, 1) + _z * unit_cell.get_cm(2, 2)), 2)),
                                    lincr,
                                    min_dist);
                        }
                    }
                }
    }
    for (int i = 0; i < atoms_with_grids; i++)
        err_checkf(num_points[i] == spherical_density[i].size(), "mismatch in number of spherical density points! i=" + toString(i), file);
#endif

    file << "                    done!" << endl;

    shrink_vector<vec>(radial_density);
    shrink_vector<vec>(radial_dist);

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_spherical = get_time();

    double _cutoff;
    if (accuracy < 3)
        _cutoff = 1E-10;
    else if (accuracy == 3)
        _cutoff = 1E-14;
    else
        _cutoff = 1E-30;
#ifndef FLO_CUDA
    ivec new_gridsize(atoms_with_grids, 0);
    ivec reductions(atoms_with_grids, 0);
    int final_size = 0;
    bool prune = true;
    if (prune)
    {
        file << "Pruning Grid..." << flush;
#pragma omp parallel for reduction(+ : final_size)
        for (int i = 0; i < atoms_with_grids; i++)
        {
            for (int p = 0; p < num_points[i]; p++)
            {
                if (grid[i][4][p] != 0.0 && abs(grid[i][3][p] * spherical_density[i][p] / grid[i][4][p]) > _cutoff)
                {
                    new_gridsize[i]++;
                }
            }
            reductions[i] = num_points[i] - new_gridsize[i];
            final_size += new_gridsize[i];
        }
        for (int k = 0; k < 7; k++)
            total_grid[k].resize(final_size);
#pragma omp parallel for
        for (int i = 0; i < atoms_with_grids; i++)
        {
            int offset = 0;
            for (int j = 0; j < i; j++)
            {
                offset += new_gridsize[j];
            }
            int reduction = 0;
            for (int p = 0; p < num_points[i]; p++)
            {
                if (grid[i][4][p] != 0.0 && abs(grid[i][3][p] * spherical_density[i][p - reduction] / grid[i][4][p]) > _cutoff)
                {
                    for (int k = 0; k < 5; k++)
                        total_grid[k][p + offset - reduction] = grid[i][k][p];
                    total_grid[6][p + offset - reduction] = grid[i][5][p];
                }
                else
                {
                    spherical_density[i].erase(spherical_density[i].begin() + (p - reduction));
                    reduction++;
                }
            }
            num_points[i] -= reduction;
            shrink_vector<vec>(grid[i]);
        }
    }
    else
    {
        for (int i = 0; i < atoms_with_grids; i++)
        {
            for (int p = 0; p < num_points[i]; p++)
            {
                for (int k = 0; k < 5; k++)
                    total_grid[k].push_back(grid[i][k][p]);
                total_grid[6].push_back(grid[i][5][p]);
            }
            shrink_vector<vec>(grid[i]);
        }
    }
    shrink_vector<vector<vec>>(grid);
    points = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
        points += num_points[i];

    // total_grid[5].resize(total_grid[0].size());
    if (debug)
        file << "sphericals done!" << endl;
    else
        file << "                                       done! Number of gridpoints: " << defaultfloat << points << endl;
#endif

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_prune = get_time();

    file << "Calculating non-spherical densities..." << flush;
    vector<vector<double>> periodic_grid;

#ifdef FLO_CUDA
    // Vector containing integrated numbers of electrons
    // dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
    // dimension 1: atoms of asym_atom_list
    vector<vector<double>> atom_els;
    atom_els.resize(3);
    for (int i = 0; i < asym_atom_list.size(); i++)
        for (int n = 0; n < 3; n++)
            atom_els[n].push_back(0.0);

    vector<float> coef;
    vector<float> ex;
    for (int i = 0; i < nex_temp; i++)
        ex.push_back(wave.get_exponent(i));

    gpuErrchk(cudaMalloc((void **)&gpu_types[0], sizeof(int) * nex_temp));
    gpuErrchk(cudaMalloc((void **)&gpu_centers[0], sizeof(int) * nex_temp));
    gpuErrchk(cudaMalloc((void **)&gpu_exponents[0], sizeof(float) * nex_temp));
    gpuErrchk(cudaMalloc((void **)&gpu_GridRho[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMemcpy(gpu_types[0], wave.get_ptr_types(), sizeof(int) * nex_temp, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_centers[0], wave.get_ptr_centers(), sizeof(int) * nex_temp, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_exponents[0], ex.data(), sizeof(float) * nex_temp, cudaMemcpyHostToDevice));

    const bool per_MO = true;
    if (per_MO)
    {
        for (int mo = 0; mo < wave.get_nmo(false); mo++)
            for (int i = 0; i < nex_temp; i++)
                if (wave.get_MO_occ(mo) != 0)
                    coef.push_back(wave.get_MO_coef(mo, i));
        if (debug)
            file << "Number of coefs: " << coef.size() << endl;
        gpuErrchk(cudaMalloc((void **)&gpu_coefficients[0], sizeof(float) * nex_temp));
        gpuErrchk(cudaMemset(gpu_GridRho[0], 0.0, sizeof(float) * MaxGrid));
        if (nex_temp < 10000)
        {
            file << "Using shared memory kernel" << endl;
            gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
                &numBlocks,
                &blocks,
                (void *)gpu_calc_dens_per_MO_shared,
                0,
                MaxGrid));
        }
        else
        {
            file << "Using static memory kernel" << endl;
            gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
                &numBlocks,
                &blocks,
                (void *)gpu_calc_dens_per_MO_static1,
                0,
                MaxGrid));
        }

        gridSize = (MaxGrid + blocks - 1) / blocks;
        if (debug)
            file << "running " << gridSize << " blocks with " << blocks << " threads" << endl;
        unsigned int nex_offset = 0;

        for (int mo = 0; mo < nmo_temp; mo++)
        {
            gpuErrchk(cudaPeekAtLastError());
            gpuErrchk(cudaDeviceSynchronize());

            if (debug)
                file << "MO " << mo << " starting at coef: " << nex_offset << endl;
            gpuErrchk(cudaMemcpy(gpu_coefficients[0], &coef[nex_offset], sizeof(float) * nex_temp, cudaMemcpyHostToDevice));
            if (nex_temp < 10000)
            {
                gpu_calc_dens_per_MO_shared<<<gridSize, blocks, nex_temp * sizeof(float)>>>(
                    gpu_GridRho[0],
                    gpu_Gridx[0],
                    gpu_Gridy[0],
                    gpu_Gridz[0],
                    gpu_PosAtomsx[0],
                    gpu_PosAtomsy[0],
                    gpu_PosAtomsz[0],
                    gpu_types[0],
                    gpu_centers[0],
                    gpu_exponents[0],
                    gpu_coefficients[0],
                    wave.get_MO_occ(mo));
            }
            else
            {
                gpu_calc_dens_per_MO_static2<<<gridSize, blocks>>>(
                    gpu_GridRho[0],
                    gpu_Gridx[0],
                    gpu_Gridy[0],
                    gpu_Gridz[0],
                    gpu_PosAtomsx[0],
                    gpu_PosAtomsy[0],
                    gpu_PosAtomsz[0],
                    gpu_types[0],
                    gpu_centers[0],
                    gpu_exponents[0],
                    gpu_coefficients[0],
                    wave.get_MO_occ(mo));
            }
            nex_offset += nex_temp;
        }
    }
    else
    {
        for (int i = 0; i < nex_temp; i++)
            for (int mo = 0; mo < wave.get_nmo(false); mo++)
                if (wave.get_MO_occ(mo) != 0)
                    coef.push_back(wave.get_MO_coef(mo, i));
        // seems broken and is slower, especially if L1 is sufficient for coefficients with size nex_temp
        gpuErrchk(cudaMalloc((void **)&gpu_coefficients[0], sizeof(float) * nex_temp * nmo_temp));
        gpuErrchk(cudaMemcpy(gpu_coefficients[0], coef.data(), sizeof(float) * nex_temp * nmo_temp, cudaMemcpyHostToDevice));
        vector<float> occ;
        for (int i = 0; i < wave.get_nmo(false); i++)
        {
            occ.push_back(wave.get_MO_occ(i));
            if (occ[occ.size() - 1] == 0)
                occ.pop_back();
        }
        gpuErrchk(cudaMalloc((void **)&gpu_occ[0], sizeof(float) * nmo_temp));
        gpuErrchk(cudaMemcpy(gpu_occ[0], occ.data(), sizeof(float) * occ.size(), cudaMemcpyHostToDevice));
        gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
            &numBlocks,
            &blocks,
            (void *)gpu_calc_dens,
            0,
            MaxGrid));

        gridSize = (MaxGrid + blocks - 1) / blocks;
        gpu_calc_dens<<<gridSize, blocks>>>(
            gpu_GridRho[0],
            gpu_Gridx[0],
            gpu_Gridy[0],
            gpu_Gridz[0],
            gpu_PosAtomsx[0],
            gpu_PosAtomsy[0],
            gpu_PosAtomsz[0],
            gpu_types[0],
            gpu_centers[0],
            gpu_exponents[0],
            gpu_coefficients[0],
            gpu_occ[0]);

        occ.clear();
        gpuErrchk(cudaFree(gpu_occ[0]));
    }
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    coef.clear();
    ex.clear();

    if (debug)
    {
        // Copy grid from GPU to print:
        //  Dimensions: [c] [p]
        //  p = the number of gridpoint
        //  c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
        gpuErrchk(cudaMemcpy(total_grid[5].data(), gpu_GridRho[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        ofstream grid("grid.file", ios::out);
        for (int i = 0; i < MaxGrid; i++)
        {
            for (int j = 0; j < 7; j++)
                grid << setw(16) << scientific << setprecision(8) << total_grid[j][i];
            grid << "\n";
        }
        grid.flush();
        grid.close();
    }

    gpuErrchk(cudaFree(gpu_types[0]));
    gpuErrchk(cudaFree(gpu_centers[0]));
    gpuErrchk(cudaFree(gpu_exponents[0]));
    gpuErrchk(cudaFree(gpu_coefficients[0]));
    free(gpu_types);
    free(gpu_centers);
    free(gpu_exponents);
    free(gpu_coefficients);
    free(gpu_occ);

    offset = 0;
    double el_sum_becke = 0.0;
    double el_sum_spherical = 0.0;
    double el_sum_hirshfeld = 0.0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        double charges[3];
        double *result = (double *)malloc(sizeof(double));
        gpuErrchk(cudaMalloc((void **)&(result), sizeof(double) * 3));
        gpu_calc_charges<<<1, 1>>>(
            gpu_GridRho[0],
            gpu_Gridaw[0],
            gpu_Gridmw[0],
            gpu_Grids[0],
            gpu_spherical_density[0][i],
            num_points[i],
            offset,
            cutoff,
            result);
        gpuErrchk(cudaMemcpy(&charges[0], result, sizeof(double) * 3, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaDeviceSynchronize());
        gpuErrchk(cudaPeekAtLastError());
        offset += num_points[i];
        el_sum_becke += charges[0];
        el_sum_spherical += charges[1];
        el_sum_hirshfeld += charges[2];
        for (int j = 0; j < 3; j++)
            atom_els[j][i] = charges[j];
    }

    offset = 0;

    file << "Applying weights..." << endl;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        cudaOccupancyMaxPotentialBlockSize(
            &numBlocks,
            &blocks,
            (void *)gpu_apply_weights,
            0,
            num_points[i]);

        gridSize = (MaxGrid + blocks - 1) / blocks;
        // file << i << ": blocks: " << gridSize << " threads: " << blocks << endl;

        gpu_apply_weights<<<gridSize, blocks, 0, streams[i]>>>(
            offset,
            gpu_GridRho[0],
            gpu_Grids[0],
            gpu_spherical_density[0][i],
            gpu_Gridaw[0],
            num_points[i]);

        offset += num_points[i];
    }
    gpuErrchk(cudaDeviceSynchronize());
    gpuErrchk(cudaPeekAtLastError());
    file << "After Applying!" << endl;

    file << endl;
    file << " done!" << endl;
    file << "Number of points evaluated: " << MaxGrid;
#else
    {
        WFN temp = wave;
        temp.delete_unoccupied_MOs();
        const int nr_atoms = (int)total_grid[0].size();
        const int nr_mos = temp.get_nmo(true);
        const int nr_cen = temp.get_ncen();
        if (debug)
        {
            file << endl
                 << "Using " << temp.get_nmo() << " MOs in temporary wavefunction" << endl;
            temp.write_wfn("temp_wavefunction.wfn", false, true);
        }
#pragma omp parallel
        {
            vector<vec> d_temp(16);
            for (int i = 0; i < 16; i++)
            {
                d_temp[i].resize(nr_cen);
            }
            vec phi_temp(nr_mos);
#pragma omp for
            for (int i = 0; i < nr_atoms; i++)
            {
                total_grid[5][i] = temp.compute_dens(
                    total_grid[0][i],
                    total_grid[1][i],
                    total_grid[2][i],
                    d_temp,
                    phi_temp,
                    false);
            }
            for (int i = 0; i < 4; i++)
                shrink_vector<double>(d_temp[i]);
            shrink_vector<vec>(d_temp);
            shrink_vector<double>(phi_temp);
        }
        // if (debug) {
        //	//Copy grid from GPU to print:
        //	// Dimensions: [c] [p]
        //	// p = the number of gridpoint
        //	// c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
        //	ofstream grid("grid.file", ios::out);
        //	for (int i = 0; i < total_grid[0].size(); i++) {
        //		for (int j = 0; j < 7; j++)
        //			grid << setw(16) << scientific << setprecision(8) << total_grid[j][i];
        //		grid << "\n";
        //	}
        //	grid.flush();
        //	grid.close();
        // }
        if (pbc != 0)
        {
            periodic_grid.resize((int)pow(pbc * 2 + 1, 3));
            int j = 0;
            for (int d = 0; d < (int)pow(pbc * 2 + 1, 3); d++)
                periodic_grid[d].resize(total_grid[5].size());
            for (int _x = -pbc; _x < pbc + 1; _x++)
                for (int _y = -pbc; _y < pbc + 1; _y++)
                    for (int _z = -pbc; _z < pbc + 1; _z++)
                    {
                        if (_x == 0 && _y == 0 && _z == 0)
                            continue;
#pragma omp parallel for
                        for (int i = 0; i < total_grid[0].size(); i++)
                        {
                            periodic_grid[j][i] = temp.compute_dens(total_grid[0][i] + _x * unit_cell.get_cm(0, 0) + _y * unit_cell.get_cm(0, 1) + _z * unit_cell.get_cm(0, 2),
                                                                    total_grid[1][i] + _x * unit_cell.get_cm(1, 0) + _y * unit_cell.get_cm(1, 1) + _z * unit_cell.get_cm(1, 2),
                                                                    total_grid[2][i] + _x * unit_cell.get_cm(2, 0) + _y * unit_cell.get_cm(2, 1) + _z * unit_cell.get_cm(2, 2), true);
                        }
                        j++;
                    }
            if (debug)
            {
                for (int i = 0; i < total_grid[0].size(); i++)
                {
                    if (i % 1000 == 0)
                        file << "Old dens: " << total_grid[5][i] << " contributions of neighbour-cells:";
                    for (int _j = 0; _j < pow(pbc * 2 + 1, 3) - 1; _j++)
                    {
                        if (i % 1000 == 0)
                            file << " " << periodic_grid[_j][i];
                        total_grid[5][i] += periodic_grid[_j][i];
                    }
                    if (i % 1000 == 0)
                        file << endl;
                }
            }
        }
    }

    if (debug)
        file << endl
             << "with total number of points: " << total_grid[0].size() << endl;
    else
        file << "                done!" << endl;

    file << "Applying hirshfeld weights and integrating charges..." << flush;
    double el_sum_becke = 0.0;
    //double el_sum_spherical = 0.0;
    double el_sum_hirshfeld = 0.0;
    // Vector containing integrated numbers of electrons
    // dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
    // dimension 1: atoms of asym_atom_list
    vector<vector<double>> atom_els(3);
    for (int n = 0; n < 3; n++)
    {
        atom_els[n].resize(asym_atom_list.size(), 0.0);
    }

    if (debug)
        file << "before loop" << endl;
        // Generate Electron sums
#pragma omp parallel for reduction(+ : el_sum_becke, el_sum_hirshfeld) //, el_sum_spherical)
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        // if (debug) file << "i=" << i << endl;
        int start_p = 0;
        for (int a = 0; a < i; a++)
            start_p += num_points[a];
        for (int p = start_p; p < start_p + num_points[i]; p++)
        {
            if (abs(total_grid[6][p]) > _cutoff)
            {
                atom_els[0][i] += total_grid[6][p] * total_grid[5][p];
                atom_els[1][i] += total_grid[6][p] * total_grid[4][p];
            }
            if (total_grid[4][p] != 0)
            {
                atom_els[2][i] += total_grid[5][p] * total_grid[3][p] * spherical_density[i][p - start_p] / total_grid[4][p];
            }
        }
        el_sum_becke += atom_els[0][i];
        //el_sum_spherical += atom_els[1][i];
        el_sum_hirshfeld += atom_els[2][i];
        if (wave.get_has_ECPs())
        {
            int n = wave.atoms[asym_atom_list[i]].ECP_electrons;
            el_sum_becke += n;
            //el_sum_spherical += n;
            el_sum_hirshfeld += n;
            atom_els[0][i] += n;
            atom_els[2][i] += n;
        }
    }

    if (debug)
    {
        file << "Becke grid with hirshfeld weights done!" << endl;
        file << "atom_els[2]: ";
        for (int i = 0; i < asym_atom_list.size(); i++)
            file << fixed << setw(10) << setprecision(3) << atom_els[2][i] << " ";
        file << endl;
    }

#pragma omp parallel for
    for (int p = 0; p < total_grid[0].size(); p++)
        total_grid[5][p] *= total_grid[3][p];
    file << " done!" << endl;
    file << "Number of points evaluated: " << total_grid[0].size();
#endif

    file << " with " << fixed << setw(10) << setprecision(6) << el_sum_becke << " electrons in Becke Grid in total." << endl
         << endl;

    file << "Table of Charges in electrons" << endl
         << endl
         << "    Atom       Becke   Spherical Hirshfeld" << endl;

    int counter = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int a = asym_atom_list[i];
        file << setw(10) << wave.atoms[a].label
             << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[0][counter]
             << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[1][counter]
             << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[2][counter];
        if (debug)
            file << " " << setw(4) << wave.get_atom_charge(a) << " " << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[0][counter]
                 << fixed << setw(10) << setprecision(3) << atom_els[1][counter]
                 << fixed << setw(10) << setprecision(3) << atom_els[2][counter];
        counter++;
        file << endl;
    }

    file << "Total number of electrons in the wavefunction: " << el_sum_becke << endl
         << " and Hirshfeld electrons (asym unit): " << el_sum_hirshfeld << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_aspherical = get_time();

    dens.resize(asym_atom_list.size());
    if (debug)
        file << "resized outer dens" << endl;
    d1.resize(asym_atom_list.size());
    d2.resize(asym_atom_list.size());
    d3.resize(asym_atom_list.size());
    if (debug)
        file << "resized outer d1-3" << endl;
#ifdef FLO_CUDA
    points = 0;
#pragma omp parallel for
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int type;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.atoms[asym_atom_list[i]].charge)
                type = j;
        vector<float> temp_dens;
        temp_dens.resize(num_points[i]);
        // file << "At atom: " << i;
        gpuErrchk(cudaMemcpy(temp_dens.data(), gpu_spherical_density[0][i], sizeof(float) * num_points[i], cudaMemcpyDeviceToHost));
        for (int p = 0; p < 0 + num_points[i]; p++)
        {
            if (abs(temp_dens[p]) > _cutoff)
            {
                dens[i].push_back(temp_dens[p]);
                d1[i].push_back(Prototype_grids[type].get_gridx(p));
                d2[i].push_back(Prototype_grids[type].get_gridy(p));
                d3[i].push_back(Prototype_grids[type].get_gridz(p));
            }
        }
        // file << " dens size: " << dens[i].size() << " num_points[i]: " << num_points[i] << endl;
    }
    for (int i = 0; i < asym_atom_list.size(); i++)
        points += dens[i].size();

    // gpuErrchk(cudaFree(gpu_Grids[0]));
    // gpuErrchk(cudaFree(gpu_PosAtomsx[0]));
    // gpuErrchk(cudaFree(gpu_PosAtomsy[0]));
    // gpuErrchk(cudaFree(gpu_PosAtomsz[0]));
    // gpuErrchk(cudaFree(gpu_GridRho[0]));
    // gpuErrchk(cudaFree(gpu_Gridx[0]));
    // gpuErrchk(cudaFree(gpu_Gridy[0]));
    // gpuErrchk(cudaFree(gpu_Gridz[0]));
    // gpuErrchk(cudaFree(gpu_Gridaw[0]));
    // gpuErrchk(cudaFree(gpu_Gridmw[0]));
    // gpuErrchk(cudaFree(gpu_asym_atom_list[0]));
    // gpuErrchk(cudaFree(gpu_atom_type_list[0]));
    // gpuErrchk(cudaFree(gpu_numpoints[0]));
    // gpuErrchk(cudaFree(gpu_atom_z[0]));
    // for (int i = 0; i < asym_atom_list.size(); i++) {
    //	gpuErrchk(cudaFree(gpu_spherical_density[0][i]));
    //	gpuErrchk(cudaFree(gpu_atomgrid_x[0][i]));
    //	gpuErrchk(cudaFree(gpu_atomgrid_y[0][i]));
    //	gpuErrchk(cudaFree(gpu_atomgrid_z[0][i]));
    // }
    free(gpu_Grids);
    free(gpu_PosAtomsx);
    free(gpu_PosAtomsy);
    free(gpu_PosAtomsz);
    free(gpu_GridRho);
    free(gpu_Gridx);
    free(gpu_Gridy);
    free(gpu_Gridz);
    free(gpu_Gridaw);
    free(gpu_Gridmw);
    free(gpu_asym_atom_list);
    free(gpu_atom_type_list);
    free(gpu_numpoints);
    free(gpu_atom_z);
    free(gpu_radial_density);
    free(gpu_radial_dist);
    free(gpu_spherical_density);
    free(gpu_atomgrid_x);
    free(gpu_atomgrid_y);
    free(gpu_atomgrid_z);
    free(gpu_atomgrid_w);
    cudaDeviceReset();
    file << "CUDA device resetted!" << endl;
#else
    points = 0;
#pragma omp parallel for
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        dens[i].resize(num_points[i]);
        d1[i].resize(num_points[i]);
        d2[i].resize(num_points[i]);
        d3[i].resize(num_points[i]);
    }
    double upper = 0, diffs = 0, avg = 0, lower = 0;
#pragma omp parallel for reduction(+ : points, upper, avg, diffs, lower)
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int start_p = 0;
        int run = 0;
        double res;
        double diff;
        double dist = 0;
        double densy;
        for (int a = 0; a < i; a++)
            start_p += num_points[a];
        for (int p = start_p; p < start_p + num_points[i]; p++)
        {
            res = total_grid[5][p] * spherical_density[i][p - start_p] / total_grid[4][p];
            if (abs(res) > _cutoff)
            {
                densy = 0;
                dens[i][run] = (res);
                d1[i][run] = (total_grid[0][p] - wave.atoms[asym_atom_list[i]].x);
                d2[i][run] = (total_grid[1][p] - wave.atoms[asym_atom_list[i]].y);
                d3[i][run] = (total_grid[2][p] - wave.atoms[asym_atom_list[i]].z);
                diff = total_grid[5][p] - total_grid[4][p] * total_grid[3][p];
                if (wave.atoms[asym_atom_list[i]].ECP_electrons != 0)
                {
                    dist = sqrt(pow(d1[i][run], 2) + pow(d2[i][run], 2) + pow(d3[i][run], 2));
                    int type_list_number = -1;
                    // Determine which type in the type list of sphericals to use
                    for (int j = 0; j < atom_type_list.size(); j++)
                        if (wave.get_atom_charge(i) == atom_type_list[j])
                            type_list_number = j;
                    densy = sphericals[type_list_number].get_core_density(dist, wave.atoms[asym_atom_list[i]].ECP_electrons) * total_grid[3][p];
                    diff += densy;
                }
                diffs += pow(diff, 2);
                upper += abs(abs(total_grid[5][p]) - abs(total_grid[4][p] * total_grid[3][p]) + densy);
                lower += abs(total_grid[5][p] + densy);
                avg += diff;
                run++;
            }
        }
        points += run;
        shrink_vector<double>(spherical_density[i]);
        dens[i].resize(run);
        d1[i].resize(run);
        d2[i].resize(run);
        d3[i].resize(run);
    }
    if (no_date == false)
    {
        file << "NRMSD value of density =              " << setw(9) << setprecision(4) << fixed << sqrt(diffs / points) / (avg / points);
        file << "\nR value of sph. vs non-sph. density = " << setw(9) << setprecision(4) << fixed << upper / lower * 100 << " %" << endl;
    }
    shrink_vector<vec>(spherical_density);
#pragma omp parallel for
    for (int _grid = 0; _grid < total_grid.size(); _grid++)
        shrink_vector<double>(total_grid[_grid]);
    shrink_vector<vec>(total_grid);
#endif
    return points;
}

/**
 * @brief Generates Hirshfeld grids for the RI method.
 *
 * This function generates Hirshfeld grids for the RI (Resolution of Identity) method.
 * It takes various parameters such as accuracy, wavefunction, coefficient filename,
 * atom type list, asymmetric atom list, needs grid vector, displacement vectors,
 * density vectors, expansion coefficients, output file stream, timing variables,
 * debug flag, and no date flag.
 *
 * @param accuracy The accuracy of the grids.
 * @param wave The wavefunction.
 * @param coef_filename The coefficient filename.
 * @param atom_type_list The list of atom types.
 * @param asym_atom_list The list of asymmetric atoms.
 * @param needs_grid The vector indicating if a grid is needed for each atom.
 * @param d1 The displacement vector d1.
 * @param d2 The displacement vector d2.
 * @param d3 The displacement vector d3.
 * @param dens The density vector.
 * @param exp_coefs The number of expansion coefficients.
 * @param file The output file stream.
 * @param start The start time point.
 * @param end_becke The end time point for Becke grid generation.
 * @param end_prototypes The end time point for prototype generation.
 * @param end_spherical The end time point for spherical grid generation.
 * @param end_prune The end time point for pruning.
 * @param end_aspherical The end time point for aspherical grid generation.
 * @param debug The debug flag.
 * @param no_date The no date flag.
 *
 * @return number of gridpoints in the final total grid
 */
static int make_hirshfeld_grids_RI(
    const int &accuracy,
    const WFN &wave,
    const string coef_filename,
    const vector<int> &atom_type_list,
    const vector<int> &asym_atom_list,
    vector<bool> &needs_grid,
    vector<vec> &d1,
    vector<vec> &d2,
    vector<vec> &d3,
    vector<vec> &dens,
    const int exp_coefs,
    ostream &file,
    time_point &start,
    time_point &end_becke,
    time_point &end_prototypes,
    time_point &end_spherical,
    time_point &end_prune,
    time_point &end_aspherical,
    bool debug,
    bool no_date)
{
    int atoms_with_grids = 0;
    for (int i = 0; i < needs_grid.size(); i++)
    {
        if (needs_grid[i])
            atoms_with_grids++;
    }
    ivec num_points(atoms_with_grids);
    // GRID COORDINATES for [a][c][p] a = atom [0,ncen],
    // c = coordinate [0=x, 1=y, 2=z, 3=atomic becke weight, 4=total spherical density, 5=molecular becke weight],
    // p = point in this grid
    vector<vector<vec>> grid(atoms_with_grids);
    const int nr_of_atoms = (wave.get_ncen());
    vec x(nr_of_atoms), y(nr_of_atoms), z(nr_of_atoms);
    ivec atom_z(nr_of_atoms);
    vec alpha_max(wave.get_ncen());
    ivec max_l(wave.get_ncen());
    int max_l_overall = 0;
#pragma omp parallel for
    for (int i = 0; i < atoms_with_grids; i++)
        grid[i].resize(6);

        // Accumulate vectors with information about all atoms
#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        const atom *ai = &(wave.atoms[i]);
        atom_z[i] = wave.get_atom_charge(i);
        x[i] = ai->x;
        y[i] = ai->y;
        z[i] = ai->z;
        alpha_max[i] = 0.0;
        max_l[i] = 0;
        for (int b = 0; b < ai->basis_set.size(); b++)
        {
            if (ai->basis_set[b].exponent > alpha_max[i])
                alpha_max[i] = ai->basis_set[b].exponent * 2;
            int l = ai->basis_set[b].type;
            if (l > max_l[i])
            {
                max_l[i] = l;
#pragma omp critical
                {
                    if (l > max_l_overall)
                        max_l_overall = l;
                }
            }
        }
    }

    if (debug)
    {
        file << "Atoms are there! max_l:" << setw(5) << max_l_overall << endl;
        for (int i = 0; i < max_l.size(); i++)
            file << "max_l: " << setw(5) << max_l[i] << endl;
    }

    vector<vec> alpha_min(wave.get_ncen());
    for (int i = 0; i < wave.get_ncen(); i++)
        alpha_min[i].resize(max_l_overall, 100000000.0);

#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        for (int b = 0; b < max_l_overall; b++)
            alpha_min[i][b] = 100000000.0;
    }

#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        const atom *ai = &(wave.atoms[i]);
        for (int b = 0; b < ai->basis_set.size(); b++)
        {
            int l = ai->basis_set[b].type;
            if (ai->basis_set[b].exponent < alpha_min[i][l])
                alpha_min[i][l] = ai->basis_set[b].exponent;
        }
    }

    if (debug)
    {
        for (int i = 0; i < wave.get_ncen(); i++)
        {
            file << "alpha_min: ";
            for (int b = 0; b < max_l_overall; b++)
                file << setw(14) << scientific << alpha_min[i][b];
            file << endl;
        }
    }

    if (debug)
        file << "alpha_min is there!" << endl
             << "Nr of asym atoms: " << asym_atom_list.size() << " Number of atoms in wfn: " << wave.get_ncen() << " atoms that needs a grid: " << atoms_with_grids << endl;
    else
        file << "There are:\n"
             << setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
             //<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
             << setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << endl;

    if (no_date)
        file << "\nMaking Becke Grids..." << flush;
    else
    {
        if (debug)
            file << "max_l_overall: " << max_l_overall << endl
                 << "Selected accuracy: " << accuracy << "\nMaking Becke Grid for" << endl;
        else
            file << endl
                 << "Selected accuracy: " << accuracy << "\nMaking Becke Grids..." << flush;
    }

    // Make Prototype grids with only single atom weights for all elements
    vector<AtomGrid> Prototype_grids;

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        if (debug)
            file << "Atom Type " << i << ": " << atom_type_list[i] << endl;
        double alpha_max_temp(0);
        int max_l_temp(0);
        vec alpha_min_temp(max_l_overall);
        for (int j = 0; j < wave.get_ncen(); j++)
        {
            if (wave.get_atom_charge(j) == 119)
            {
                continue;
            }
            if (wave.get_atom_charge(j) == atom_type_list[i])
            {
                if (debug)
                {
                    file << alpha_max[j] << " " << max_l[j] - 1 << " ";
                    for (int l = 0; l < max_l_overall; l++)
                        file << alpha_min[j][l] << " ";
                    file << endl;
                }
                alpha_max_temp = alpha_max[j];
                max_l_temp = max_l[j];
                for (int l = 0; l <= max_l_temp; l++)
                    alpha_min_temp[l] = alpha_min[j][l];
                break;
            }
        }

        if (debug)
        {
            file << "max_l: " << defaultfloat << max_l_temp << " alpha_max: " << scientific << alpha_max_temp << " alpha_min: ";
            for (int l = 0; l <= max_l_temp; l++)
                file << setw(14) << scientific << alpha_min_temp[l];
            file << " accuracy: " << accuracy << endl;
        }
        int lebedev_high, lebedev_low;
        double radial_acc;
        err_checkf(accuracy >= 0, "Negative accuracy is not defined!", file);
        if (accuracy == 0)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                radial_acc = 1e-4;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                radial_acc = 1e-3;
            }
        }
        else if (accuracy == 1)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[7] : constants::lebedev_table[8];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[3] : constants::lebedev_table[4];
                radial_acc = 1e-4;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[6] : constants::lebedev_table[7];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[2] : constants::lebedev_table[3];
                radial_acc = 1e-5;
            }
        }
        else if (accuracy == 2)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[11] : constants::lebedev_table[12];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[7] : constants::lebedev_table[8];
                radial_acc = 1e-5;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[10] : constants::lebedev_table[11];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[6] : constants::lebedev_table[7];
                radial_acc = 1e-6;
            }
        }
        else if (accuracy == 3)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[14] : constants::lebedev_table[16];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[12] : constants::lebedev_table[14];
                radial_acc = 1e-10;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[13] : constants::lebedev_table[15];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[11] : constants::lebedev_table[13];
                radial_acc = 1e-11;
            }
        }
        else if (accuracy == 4)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[19] : constants::lebedev_table[21];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[14] : constants::lebedev_table[17];
                radial_acc = 1e-20;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[18] : constants::lebedev_table[20];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[13] : constants::lebedev_table[16];
                radial_acc = 1e-15;
            }
        }
        else
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[31] : constants::lebedev_table[32];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[29] : constants::lebedev_table[31];
                radial_acc = 1e-20;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[30] : constants::lebedev_table[32];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[28] : constants::lebedev_table[30];
                radial_acc = 1e-15;
            }
        }
        Prototype_grids.push_back(AtomGrid(radial_acc,
                                           lebedev_low,
                                           lebedev_high,
                                           atom_type_list[i],
                                           alpha_max_temp,
                                           max_l_temp,
                                           alpha_min_temp.data(),
                                           file));
    }

    end_prototypes = get_time();
    if (debug)
    {
        int dur = get_sec(start, end_prototypes);
        for (int prototype = 0; prototype < Prototype_grids.size(); prototype++)
            file << "Number of gridpoints for atom type " << atom_type_list[prototype] << ": " << Prototype_grids[prototype].get_num_grid_points() << endl;

        //	int diff = end - start;
        if (dur < 1)
            file << "Time until prototypes are done: " << fixed << setprecision(0) << get_msec(start, end_prototypes) << " ms" << endl;
        else
            file << "Time until prototypes are done: " << fixed << setprecision(0) << dur << " s" << endl;
    }

    vector<vec> radial_density(atom_type_list.size());
    vector<vec> radial_dist(atom_type_list.size());
    if (!debug)
    {
        file << " ...  " << flush;
    }
    // get_grid is parallelized, therefore not parallel here
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        // skip atoms, that do not need a grid
        if (!needs_grid[i])
        {
            continue;
        }
        if (debug)
        {
            file << "Making grid for atom " << i << endl;
        }
        int type = 0;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.get_atom_charge(i))
                type = j;

        int grid_number = 0;
        for (int j = 0; j < i; j++)
        {
            if (needs_grid[j])
            {
                grid_number++;
            }
        }
        num_points[grid_number] = Prototype_grids[type].get_num_grid_points();

        for (int n = 0; n < 6; n++)
            grid[grid_number][n].resize(num_points[grid_number], 0.0);

        Prototype_grids[type].get_grid(int(wave.get_ncen()),
                                       i,
                                       &x[0],
                                       &y[0],
                                       &z[0],
                                       &atom_z[0],
                                       grid[grid_number][0].data(),
                                       grid[grid_number][1].data(),
                                       grid[grid_number][2].data(),
                                       grid[grid_number][3].data(),
                                       grid[grid_number][5].data());
    }
    if (debug)
    {
        int run = 0;
        file << "  label | needs_grid | number of gridpoints\n";
        for (int j = 0; j < wave.get_ncen(); j++)
        {
            file << setw(8) << wave.atoms[j].label << setw(13) << needs_grid[j];
            if (needs_grid[j])
            {
                file << setw(7) << num_points[run];
                run++;
            }
            else
            {
                file << setw(6) << "---";
            }
            file << endl;
        }
    }
    Prototype_grids.clear();

    int points = 0;
    for (int i = 0; i < atoms_with_grids; i++)
        points += num_points[i];
    if (debug)
        file << "Becke Grid exists" << endl;
    else
        file << "                           done! Number of gridpoints: " << defaultfloat << points << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_becke = get_time();

    file << "Calculating spherical densities..." << flush;

    // Total grid as a sum of all atomic grids.
    // Dimensions: [c] [p]
    // p = the number of gridpoint
    // c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
    vector<vec> total_grid(7);
    // density of spherical atom at each
    // Dimensions: [a] [d]
    // a = atom number in atom type list for which the weight is calcualted
    // d = distance to look at obtained from point_to_distance_map
    vector<vec> spherical_density(atoms_with_grids);

    const double incr = pow(1.005, max(1, accuracy - 1));
    const double lincr = log(incr);
    const double min_dist = 0.0000001;
    vector<Thakkar> sphericals;
    for (int i = 0; i < atom_type_list.size(); i++)
        sphericals.push_back(Thakkar(atom_type_list[i]));
    // Make radial grids
    if (debug)
    {
        file << "\nSize of atom_type_list:" << setw(5) << atom_type_list.size() << endl;
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            file << "\nCalculating for atomic number " << atom_type_list[i] << endl;
            double current = 1;
            double dist = min_dist;
            if (accuracy > 3)
                while (current > 1E-10)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    if (current == -20)
                        return false;
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    if (current == -20)
                        return false;
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            file << "Number of radial density points for atomic number " << atom_type_list[i] << ": " << radial_density[i].size() << endl;
        }
    }
    else
    {
#pragma omp parallel for
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            double current = 1;
            double dist = min_dist;
            if (accuracy > 3)
                while (current > 1E-10)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
        }
    }
    sphericals.clear();
    int type_list_number = -1;
    if (debug)
    {
        file << "Cleared the sphericals!" << endl;
    }
#pragma omp parallel
    {
#pragma omp for
        for (int g = 0; g < atoms_with_grids; g++)
        {
            spherical_density[g].resize(num_points[g]);
        }
        for (int i = 0; i < wave.get_ncen(); i++)
        {
#pragma omp single
            {
                type_list_number = -1;
                // Determine which type in the type list of sphericals to use
                for (int j = 0; j < atom_type_list.size(); j++)
                    if (wave.get_atom_charge(i) == atom_type_list[j])
                        type_list_number = j;
                if (debug && type_list_number != -1)
                {
                    file << type_list_number << " Atom type: " << atom_type_list[type_list_number] << endl;
                }
            }
            if (type_list_number == -1)
            {
#pragma omp single
                {
                    file << "I skipped an atom! make sure this is okay!" << endl;
                }
                continue;
            }
#pragma omp for
            for (int g = 0; g < atoms_with_grids; g++)
            {
                for (int p = 0; p < num_points[g]; p++)
                {
                    double temp =
                        linear_interpolate_spherical_density(radial_density[type_list_number], radial_dist[type_list_number], sqrt(pow(grid[g][0][p] - wave.atoms[i].x, 2) + pow(grid[g][1][p] - wave.atoms[i].y, 2) + pow(grid[g][2][p] - wave.atoms[i].z, 2)), lincr, min_dist);
                    if (i == asym_atom_list[g])
                        spherical_density[g][p] = temp;
                    grid[g][4][p] += temp;
                }
            }
        }
    }
    for (int i = 0; i < atoms_with_grids; i++)
        err_checkf(num_points[i] == spherical_density[i].size(), "mismatch in number of spherical density points! i=" + toString(i), file);

    file << "                    done!" << endl;

    shrink_vector<vec>(radial_density);
    shrink_vector<vec>(radial_dist);

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_spherical = get_time();

    double _cutoff;
    if (accuracy < 3)
        _cutoff = 1E-10;
    else if (accuracy == 3)
        _cutoff = 1E-14;
    else
        _cutoff = 1E-30;
#ifndef FLO_CUDA
    ivec new_gridsize(atoms_with_grids, 0);
    ivec reductions(atoms_with_grids, 0);
    int final_size = 0;
    bool prune = true;
    if (prune)
    {
        file << "Pruning Grid..." << flush;
#pragma omp parallel for reduction(+ : final_size)
        for (int i = 0; i < atoms_with_grids; i++)
        {
            for (int p = 0; p < num_points[i]; p++)
            {
                if (grid[i][4][p] != 0.0 && abs(grid[i][3][p] * spherical_density[i][p] / grid[i][4][p]) > _cutoff)
                {
                    new_gridsize[i]++;
                }
            }
            reductions[i] = num_points[i] - new_gridsize[i];
            final_size += new_gridsize[i];
        }
        for (int k = 0; k < 7; k++)
            total_grid[k].resize(final_size);
#pragma omp parallel for
        for (int i = 0; i < atoms_with_grids; i++)
        {
            int offset = 0;
            for (int j = 0; j < i; j++)
            {
                offset += new_gridsize[j];
            }
            int reduction = 0;
            for (int p = 0; p < num_points[i]; p++)
            {
                if (grid[i][4][p] != 0.0 && abs(grid[i][3][p] * spherical_density[i][p - reduction] / grid[i][4][p]) > _cutoff)
                {
                    for (int k = 0; k < 5; k++)
                        total_grid[k][p + offset - reduction] = grid[i][k][p];
                    total_grid[6][p + offset - reduction] = grid[i][5][p];
                }
                else
                {
                    spherical_density[i].erase(spherical_density[i].begin() + (p - reduction));
                    reduction++;
                }
            }
            num_points[i] -= reduction;
            shrink_vector<vec>(grid[i]);
        }
    }
    else
    {
        for (int i = 0; i < atoms_with_grids; i++)
        {
            for (int p = 0; p < num_points[i]; p++)
            {
                for (int k = 0; k < 5; k++)
                    total_grid[k].push_back(grid[i][k][p]);
                total_grid[6].push_back(grid[i][5][p]);
            }
            shrink_vector<vec>(grid[i]);
        }
    }
    shrink_vector<vector<vec>>(grid);
    points = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
        points += num_points[i];

    // total_grid[5].resize(total_grid[0].size());
    if (debug)
        file << "sphericals done!" << endl;
    else
        file << "                                       done! Number of gridpoints: " << defaultfloat << points << endl;
#endif

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_prune = get_time();

    file << "Calculating non-spherical densities..." << flush;
    vector<vector<double>> periodic_grid;

    {
        WFN temp = wave;
        const int nr_pts = (int)total_grid[0].size();
        vector<unsigned long> shape{};
        bool fortran_order;
        vec data{};

        string path{coef_filename};
        npy::LoadArrayFromNumpy(path, shape, fortran_order, data);

#pragma omp parallel for
        for (int i = 0; i < nr_pts; i++)
        {
            total_grid[5][i] = calc_density_ML(
                total_grid[0][i],
                total_grid[1][i],
                total_grid[2][i],
                data,
                temp.atoms,
                exp_coefs);
        }
        shrink_vector<double>(data);
    }

    if (debug)
        file << endl
             << "with total number of points: " << total_grid[0].size() << endl;
    else
        file << "                done!" << endl;

    file << "Applying hirshfeld weights and integrating charges..." << flush;
    double el_sum_becke = 0.0;
    //double el_sum_spherical = 0.0;
    double el_sum_hirshfeld = 0.0;
    // Vector containing integrated numbers of electrons
    // dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
    // dimension 1: atoms of asym_atom_list
    vector<vector<double>> atom_els(3);
    for (int n = 0; n < 3; n++)
    {
        atom_els[n].resize(asym_atom_list.size(), 0.0);
    }

    if (debug)
        file << "before loop" << endl;
        // Generate Electron sums
#pragma omp parallel for reduction(+ : el_sum_becke, el_sum_hirshfeld) // el_sum_spherical)
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        // if (debug) file << "i=" << i << endl;
        int start_p = 0;
        for (int a = 0; a < i; a++)
            start_p += num_points[a];
        for (int p = start_p; p < start_p + num_points[i]; p++)
        {
            if (abs(total_grid[6][p]) > _cutoff)
            {
                atom_els[0][i] += total_grid[6][p] * total_grid[5][p]; // Molecular grid * WFN rho
                atom_els[1][i] += total_grid[6][p] * total_grid[4][p]; // Molecular grid * spheircal rho
            }
            if (total_grid[4][p] != 0)
            {
                // WFN rho * atomic weight * hirshfeld weight
                atom_els[2][i] += total_grid[5][p] * total_grid[3][p] * spherical_density[i][p - start_p] / total_grid[4][p];
            }
        }
        el_sum_becke += atom_els[0][i];
        //el_sum_spherical += atom_els[1][i];
        el_sum_hirshfeld += atom_els[2][i];
        if (wave.get_has_ECPs())
        {
            int n = wave.atoms[asym_atom_list[i]].ECP_electrons;
            el_sum_becke += n;
            //el_sum_spherical += n;
            el_sum_hirshfeld += n;
            atom_els[0][i] += n;
            atom_els[2][i] += n;
        }
    }

    if (debug)
    {
        file << "Becke grid with hirshfeld weights done!" << endl;
        file << "atom_els[2]: ";
        for (int i = 0; i < asym_atom_list.size(); i++)
            file << fixed << setw(10) << setprecision(3) << atom_els[2][i] << " ";
        file << endl;
    }

#pragma omp parallel for
    for (int p = 0; p < total_grid[0].size(); p++)
        total_grid[5][p] *= total_grid[3][p];
    file << " done!" << endl;
    file << "Number of points evaluated: " << total_grid[0].size();

    file << " with " << fixed << setw(10) << setprecision(6) << el_sum_becke << " electrons in Becke Grid in total." << endl
         << endl;

    file << "Table of Charges in electrons" << endl
         << endl
         << "    Atom       Becke   Spherical Hirshfeld" << endl;

    int counter = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int a = asym_atom_list[i];
        file << setw(10) << wave.atoms[a].label
             << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[0][counter]
             << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[1][counter]
             << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[2][counter];
        if (debug)
            file << " " << setw(4) << wave.get_atom_charge(a) << " " << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[0][counter]
                 << fixed << setw(10) << setprecision(3) << atom_els[1][counter]
                 << fixed << setw(10) << setprecision(3) << atom_els[2][counter];
        counter++;
        file << endl;
    }

    file << "Total number of electrons in the wavefunction: " << el_sum_becke << endl
         << " and Hirshfeld electrons (asym unit): " << el_sum_hirshfeld << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_aspherical = get_time();

    dens.resize(asym_atom_list.size());
    if (debug)
        file << "resized outer dens" << endl;
    d1.resize(asym_atom_list.size());
    d2.resize(asym_atom_list.size());
    d3.resize(asym_atom_list.size());
    if (debug)
        file << "resized outer d1-3" << endl;
    points = 0;
#pragma omp parallel for
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        dens[i].resize(num_points[i]);
        d1[i].resize(num_points[i]);
        d2[i].resize(num_points[i]);
        d3[i].resize(num_points[i]);
    }
    double upper = 0, diffs = 0, avg = 0, lower = 0;
#pragma omp parallel for reduction(+ : points, upper, avg, diffs, lower)
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int start_p = 0;
        int run = 0;
        double res;
        double diff;
        double dist = 0;
        double densy;
        for (int a = 0; a < i; a++)
            start_p += num_points[a];
        for (int p = start_p; p < start_p + num_points[i]; p++)
        {
            res = total_grid[5][p] * spherical_density[i][p - start_p] / total_grid[4][p];
            if (abs(res) > _cutoff)
            {
                densy = 0;
                dens[i][run] = (res);
                d1[i][run] = (total_grid[0][p] - wave.atoms[asym_atom_list[i]].x);
                d2[i][run] = (total_grid[1][p] - wave.atoms[asym_atom_list[i]].y);
                d3[i][run] = (total_grid[2][p] - wave.atoms[asym_atom_list[i]].z);
                diff = total_grid[5][p] - total_grid[4][p] * total_grid[3][p];
                if (wave.atoms[asym_atom_list[i]].ECP_electrons != 0)
                {
                    dist = sqrt(pow(d1[i][run], 2) + pow(d2[i][run], 2) + pow(d3[i][run], 2));
                    int type_list_number = -1;
                    // Determine which type in the type list of sphericals to use
                    for (int j = 0; j < atom_type_list.size(); j++)
                        if (wave.get_atom_charge(i) == atom_type_list[j])
                            type_list_number = j;
                    densy = sphericals[type_list_number].get_core_density(dist, wave.atoms[asym_atom_list[i]].ECP_electrons);
                    diff += densy;
                }
                diffs += pow(diff, 2);
                upper += abs(abs(total_grid[5][p]) - abs(total_grid[4][p] * total_grid[3][p]) + densy);
                lower += abs(total_grid[5][p] + densy);
                avg += diff;
                run++;
            }
        }
        points += run;
        shrink_vector<double>(spherical_density[i]);
        dens[i].resize(run);
        d1[i].resize(run);
        d2[i].resize(run);
        d3[i].resize(run);
    }
    if (no_date == false)
    {
        file << "NRMSD value of density =              " << setw(9) << setprecision(4) << fixed << sqrt(diffs / points) / (avg / points);
        file << "\nR value of sph. vs non-sph. density = " << setw(9) << setprecision(4) << fixed << upper / lower * 100 << " %" << endl;
    }
    shrink_vector<vec>(spherical_density);
#pragma omp parallel for
    for (int _grid = 0; _grid < total_grid.size(); _grid++)
        shrink_vector<double>(total_grid[_grid]);
    shrink_vector<vec>(total_grid);
    return points;
}

/**
 * @brief Generates integration grids for calculating scattering factors.
 *
 * This function generates integration grids based on the given parameters and stores them in the provided vectors.
 *
 * @param accuracy The desired accuracy of the integration grids.
 * @param unit_cell The unit cell object containing information about the crystal lattice.
 * @param wave The WFN object containing information about the wavefunction.
 * @param coef_filename The filename of the coefficient file.
 * @param atom_type_list The list of atom types.
 * @param asym_atom_list The list of asymmetric atoms.
 * @param needs_grid The vector indicating whether each atom needs a grid.
 * @param d1 The vector to store the first integration grid.
 * @param d2 The vector to store the second integration grid.
 * @param d3 The vector to store the third integration grid.
 * @param dens The vector to store the density.
 * @param exp_coefs The number of expansion coefficients.
 * @param file The output stream to write the integration grids.
 * @param start The start time point for measuring execution time.
 * @param end_becke The end time point for measuring execution time of the Becke grid generation.
 * @param end_prototypes The end time point for measuring execution time of the prototype grid generation.
 * @param end_spherical The end time point for measuring execution time of the spherical grid generation.
 * @param end_prune The end time point for measuring execution time of the pruning step.
 * @param end_aspherical The end time point for measuring execution time of the aspherical grid generation.
 * @param debug Flag indicating whether to enable debug mode.
 * @param no_date Flag indicating whether to exclude the date from the output.
 *
 * @return The number of points in the integration grids.
 */
static int make_integration_grids(
    const int &accuracy,
    cell &unit_cell,
    const WFN &wave,
    const string coef_filename,
    const vector<int> &atom_type_list,
    const vector<int> &asym_atom_list,
    vector<bool> &needs_grid,
    vector<vec> &d1,
    vector<vec> &d2,
    vector<vec> &d3,
    vector<vec> &dens,
    const int exp_coefs,
    ostream &file,
    time_point &start,
    time_point &end_becke,
    time_point &end_prototypes,
    time_point &end_spherical,
    time_point &end_prune,
    time_point &end_aspherical,
    bool debug,
    bool no_date)
{
    int atoms_with_grids = 0;
    for (int i = 0; i < needs_grid.size(); i++)
    {
        if (needs_grid[i])
            atoms_with_grids++;
    }
    // counts number of points inside each atomic grid
    ivec num_points(atoms_with_grids);
    // GRID COORDINATES for [a][c][p]
    // a = atom [0,ncen],
    // c = coordinate [0=x, 1=y, 2=z, 3=atomic becke weight, 4=molecular weight],
    // p = point in this grid
    vector<vector<vec>> grid(atoms_with_grids);
#pragma omp parallel for
    for (int i = 0; i < atoms_with_grids; i++)
        grid[i].resize(5);

    const int nr_of_atoms = (wave.get_ncen());
    // positions
    vec x(nr_of_atoms), y(nr_of_atoms), z(nr_of_atoms);
    // charges
    ivec atom_z(nr_of_atoms);
    vec alpha_max(wave.get_ncen());
    ivec max_l(wave.get_ncen());
    int max_l_overall = 0;

    // Accumulate vectors with information about all atoms
#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        const atom *ai = &(wave.atoms[i]);
        atom_z[i] = wave.get_atom_charge(i);
        x[i] = ai->x;
        y[i] = ai->y;
        z[i] = ai->z;
        alpha_max[i] = 0.0;
        max_l[i] = 0;
        for (int b = 0; b < ai->basis_set.size(); b++)
        {
            if (ai->basis_set[b].exponent > alpha_max[i])
                alpha_max[i] = ai->basis_set[b].exponent * 2;
            int l = ai->basis_set[b].type;
            if (l > max_l[i])
            {
                max_l[i] = l;
#pragma omp critical
                {
                    if (l > max_l_overall)
                        max_l_overall = l;
                }
            }
        }
    }

    if (debug)
    {
        file << "Atoms are there! max_l:" << setw(5) << max_l_overall << endl;
        for (int i = 0; i < max_l.size(); i++)
            file << "max_l: " << setw(5) << max_l[i] << endl;
    }

    vector<vec> alpha_min(wave.get_ncen());
    for (int i = 0; i < wave.get_ncen(); i++)
        alpha_min[i].resize(max_l_overall, 100000000.0);

#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        for (int b = 0; b < max_l_overall; b++)
            alpha_min[i][b] = 100000000.0;
    }

#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        const atom *ai = &(wave.atoms[i]);
        for (int b = 0; b < ai->basis_set.size(); b++)
        {
            int l = ai->basis_set[b].type;
            if (ai->basis_set[b].exponent < alpha_min[i][l])
                alpha_min[i][l] = ai->basis_set[b].exponent;
        }
    }

    if (debug)
    {
        for (int i = 0; i < wave.get_ncen(); i++)
        {
            file << "alpha_min: ";
            for (int b = 0; b < max_l_overall; b++)
                file << setw(14) << scientific << alpha_min[i][b];
            file << endl;
        }
    }

    if (debug)
        file << "alpha_min is there!" << endl
             << "Nr of asym atoms: " << asym_atom_list.size() << " Number of atoms in wfn: " << wave.get_ncen() << " atoms that needs a grid: " << atoms_with_grids << endl;
    else
        file << "There are:\n"
             << setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
             //<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
             << setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << endl;

    if (no_date)
        file << "\nMaking Becke Grids..." << flush;
    else
    {
        if (debug)
            file << "max_l_overall: " << max_l_overall << endl
                 << "Selected accuracy: " << accuracy << "\nMaking Becke Grid for" << endl;
        else
            file << endl
                 << "Selected accuracy: " << accuracy << "\nMaking Becke Grids..." << flush;
    }

    // Make Prototype grids with only single atom weights for all elements
    vector<AtomGrid> Prototype_grids;

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        if (debug)
            file << "Atom Type " << i << ": " << atom_type_list[i] << endl;
        double alpha_max_temp(0);
        int max_l_temp(0);
        vec alpha_min_temp(max_l_overall);
        for (int j = 0; j < wave.get_ncen(); j++)
        {
            if (wave.get_atom_charge(j) == 119)
            {
                continue;
            }
            if (wave.get_atom_charge(j) == atom_type_list[i])
            {
                if (debug)
                {
                    file << alpha_max[j] << " " << max_l[j] - 1 << " ";
                    for (int l = 0; l < max_l_overall; l++)
                        file << alpha_min[j][l] << " ";
                    file << endl;
                }
                alpha_max_temp = alpha_max[j];
                max_l_temp = max_l[j];
                for (int l = 0; l <= max_l_temp; l++)
                    alpha_min_temp[l] = alpha_min[j][l];
                break;
            }
        }

        if (debug)
        {
            file << "max_l: " << defaultfloat << max_l_temp << " alpha_max: " << scientific << alpha_max_temp << " alpha_min: ";
            for (int l = 0; l <= max_l_temp; l++)
                file << setw(14) << scientific << alpha_min_temp[l];
            file << " accuracy: " << accuracy << endl;
        }
        int lebedev_high, lebedev_low;
        double radial_acc;
        err_checkf(accuracy >= 0, "Negative accuracy is not defined!", file);
        if (accuracy == 0)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                radial_acc = 1e-4;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                radial_acc = 1e-3;
            }
        }
        else if (accuracy == 1)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[7] : constants::lebedev_table[8];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[3] : constants::lebedev_table[4];
                radial_acc = 1e-4;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[6] : constants::lebedev_table[7];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[2] : constants::lebedev_table[3];
                radial_acc = 1e-5;
            }
        }
        else if (accuracy == 2)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[11] : constants::lebedev_table[12];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[7] : constants::lebedev_table[8];
                radial_acc = 1e-5;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[10] : constants::lebedev_table[11];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[6] : constants::lebedev_table[7];
                radial_acc = 1e-6;
            }
        }
        else if (accuracy == 3)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[14] : constants::lebedev_table[16];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[12] : constants::lebedev_table[14];
                radial_acc = 1e-10;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[13] : constants::lebedev_table[15];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[11] : constants::lebedev_table[13];
                radial_acc = 1e-11;
            }
        }
        else if (accuracy == 4)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[19] : constants::lebedev_table[21];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[14] : constants::lebedev_table[17];
                radial_acc = 1e-20;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[18] : constants::lebedev_table[20];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[13] : constants::lebedev_table[16];
                radial_acc = 1e-15;
            }
        }
        else
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[31] : constants::lebedev_table[32];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[29] : constants::lebedev_table[31];
                radial_acc = 1e-20;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[30] : constants::lebedev_table[32];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[28] : constants::lebedev_table[30];
                radial_acc = 1e-15;
            }
        }
        Prototype_grids.push_back(AtomGrid(radial_acc,
                                           lebedev_low,
                                           lebedev_high,
                                           atom_type_list[i],
                                           alpha_max_temp,
                                           max_l_temp,
                                           alpha_min_temp.data(),
                                           file));
    }

    end_prototypes = get_time();
    if (debug)
    {

        for (int prototype = 0; prototype < Prototype_grids.size(); prototype++)
            file << "Number of gridpoints for atom type " << atom_type_list[prototype] << ": " << Prototype_grids[prototype].get_num_grid_points() << endl;

        int dur = get_sec(start, end_prototypes);
        if (dur < 1)
            file << "Time until prototypes are done: " << fixed << setprecision(0) << get_msec(start, end_prototypes) << " ms" << endl;
        else
            file << "Time until prototypes are done: " << fixed << setprecision(0) << dur << " s" << endl;
    }

    if (!debug)
    {
        file << " ...  " << flush;
    }
    // get_grid is parallelized, therefore not parallel here
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        // skip atoms, that do not need a grid
        if (!needs_grid[i])
        {
            continue;
        }
        if (debug)
        {
            file << "Making grid for atom " << i << endl;
        }
        int type = 0;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.get_atom_charge(i))
                type = j;

        int grid_number = 0;
        for (int j = 0; j < i; j++)
        {
            if (needs_grid[j])
            {
                grid_number++;
            }
        }
        num_points[grid_number] = Prototype_grids[type].get_num_grid_points();

        for (int n = 0; n < grid[grid_number].size(); n++)
            grid[grid_number][n].resize(num_points[grid_number], 0.0);

        Prototype_grids[type].get_grid(wave.get_ncen(),
                                       i,
                                       &x[0],
                                       &y[0],
                                       &z[0],
                                       &atom_z[0],
                                       grid[grid_number][0].data(),
                                       grid[grid_number][1].data(),
                                       grid[grid_number][2].data(),
                                       grid[grid_number][3].data(),
                                       grid[grid_number][4].data());
    }
    if (debug)
    {
        int run = 0;
        file << "  label | needs_grid | number of gridpoints\n";
        for (int j = 0; j < wave.get_ncen(); j++)
        {
            file << setw(8) << wave.atoms[j].label << setw(13) << needs_grid[j];
            if (needs_grid[j])
            {
                file << setw(7) << num_points[run];
                run++;
            }
            else
            {
                file << setw(6) << "---";
            }
            file << endl;
        }
    }
    Prototype_grids.clear();

    int points = 0;
    for (int i = 0; i < atoms_with_grids; i++)
        points += num_points[i];
    if (debug)
        file << "Becke Grid exists" << endl;
    else
        file << "                           done! Number of gridpoints: " << defaultfloat << points << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_becke = get_time();
    // Total grid as a sum of all atomic grids.
    // Dimensions: [c] [p]
    // p = the number of gridpoint
    // c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=wavefunction density, 5=MW
    vector<vec> total_grid(6);

    //int type_list_number = -1;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_spherical = get_time();

    double _cutoff;
    if (accuracy < 3)
        _cutoff = 1E-10;
    else if (accuracy == 3)
        _cutoff = 1E-14;
    else
        _cutoff = 1E-30;
    ivec new_gridsize(atoms_with_grids, 0);
    ivec reductions(atoms_with_grids, 0);
    int final_size = 0;
    bool prune = true;
    if (prune)
    {
        file << "Pruning Grid..." << flush;
#pragma omp parallel for reduction(+ : final_size)
        for (int i = 0; i < atoms_with_grids; i++)
        {
            for (int p = 0; p < num_points[i]; p++)
            {
                if (abs(grid[i][4][p]) > _cutoff)
                {
                    new_gridsize[i]++;
                }
            }
            reductions[i] = num_points[i] - new_gridsize[i];
            final_size += new_gridsize[i];
        }
        for (int k = 0; k < total_grid.size(); k++)
            total_grid[k].resize(final_size);
#pragma omp parallel for
        for (int i = 0; i < atoms_with_grids; i++)
        {
            int offset = 0;
            for (int j = 0; j < i; j++)
            {
                offset += new_gridsize[j];
            }
            int reduction = 0;
            for (int p = 0; p < num_points[i]; p++)
            {
                if (abs(grid[i][4][p]) > _cutoff)
                {
                    for (int k = 0; k < 4; k++)
                        total_grid[k][p + offset - reduction] = grid[i][k][p];
                    total_grid[5][p + offset - reduction] = grid[i][4][p];
                }
                else
                {
                    reduction++;
                }
            }
            num_points[i] -= reduction;
            shrink_vector<vec>(grid[i]);
        }
    }
    else
    {
        for (int i = 0; i < atoms_with_grids; i++)
        {
            for (int p = 0; p < num_points[i]; p++)
            {
                for (int k = 0; k < 4; k++)
                    total_grid[k].push_back(grid[i][k][p]);
                total_grid[5].push_back(grid[i][4][p]);
            }
            shrink_vector<vec>(grid[i]);
        }
    }
    shrink_vector<vector<vec>>(grid);
    points = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
        points += num_points[i];

    file << "                                       done! Number of gridpoints: " << defaultfloat << points << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_prune = get_time();

    file << "Calculating non-spherical densities..." << flush;

    WFN temp = wave;
    const int nr_pts = (int)total_grid[0].size();
    vector<unsigned long> shape{};
    bool fortran_order;
    vec data{};

    string path{coef_filename};
    npy::LoadArrayFromNumpy(path, shape, fortran_order, data);

#pragma omp parallel for
    for (int i = 0; i < nr_pts; i++)
    {
        total_grid[4][i] = calc_density_ML(
            total_grid[0][i],
            total_grid[1][i],
            total_grid[2][i],
            data,
            temp.atoms,
            exp_coefs);
    }
    shrink_vector<double>(data);

    if (debug)
        file << endl
             << "with total number of points: " << total_grid[0].size() << endl;
    else
        file << "                done!" << endl;

    file << "Applying weights and integrating charges...          " << flush;
    double el_sum_SALTED = 0.0;
    // Vector containing integrated numbers of electrons
    // dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
    // dimension 1: atoms of asym_atom_list
    vector<vec> atom_els(3);
    for (int n = 0; n < atom_els.size(); n++)
    {
        atom_els[n].resize(asym_atom_list.size(), 0.0);
    }

    if (debug)
        file << "before loop" << endl;
        // Generate Electron sums
#pragma omp parallel for reduction(+ : el_sum_SALTED)
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int start_p = 0;
        for (int a = 0; a < i; a++)
            start_p += num_points[a];
        for (int p = start_p; p < start_p + num_points[i]; p++)
        {
            atom_els[2][i] += total_grid[4][p] * total_grid[5][p];
        }
        el_sum_SALTED += atom_els[2][i];
        if (wave.get_has_ECPs())
        {
            int n = wave.atoms[asym_atom_list[i]].ECP_electrons;
            el_sum_SALTED += n;
            atom_els[2][i] += n;
        }
    }

    if (debug)
    {
        file << "Becke grid with hirshfeld weights done!" << endl;
        file << "atom_els[2]: ";
        for (int i = 0; i < asym_atom_list.size(); i++)
            file << fixed << setw(10) << setprecision(3) << atom_els[2][i] << " ";
        file << endl;
    }

#pragma omp parallel for
    for (int p = 0; p < total_grid[0].size(); p++)
        total_grid[4][p] *= total_grid[5][p];
    file << " done!" << endl;
    file << "Number of points evaluated: " << total_grid[0].size();

    file << " with " << fixed << setw(10) << setprecision(6) << el_sum_SALTED << " electrons in Becke Grid in total." << endl
         << endl;

    file << "Table of Charges in electrons" << endl
         << endl
         << "    Atom      Charge" << endl;

    int counter = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int a = asym_atom_list[i];
        file << setw(10) << wave.atoms[a].label
             << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[2][counter];
        counter++;
        file << endl;
    }

    file << "Total number of electrons in the wavefunction: " << el_sum_SALTED << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_aspherical = get_time();

    dens.resize(asym_atom_list.size());
    if (debug)
        file << "resized outer dens" << endl;
    d1.resize(asym_atom_list.size());
    d2.resize(asym_atom_list.size());
    d3.resize(asym_atom_list.size());
    if (debug)
        file << "resized outer d1-3" << endl;

    points = 0;
#pragma omp parallel
    {
#pragma omp for reduction(+ : points)
        for (int i = 0; i < asym_atom_list.size(); i++)
        {
            dens[i].resize(num_points[i]);
            d1[i].resize(num_points[i]);
            d2[i].resize(num_points[i]);
            d3[i].resize(num_points[i]);
            int start_p = 0;
            int run = 0;
            for (int a = 0; a < i; a++)
                start_p += num_points[a];
            for (int p = start_p; p < start_p + num_points[i]; p++)
            {
                if (abs(total_grid[4][p]) > _cutoff)
                {
                    dens[i][run] = (total_grid[4][p]);
                    d1[i][run] = (total_grid[0][p] - wave.atoms[asym_atom_list[i]].x);
                    d2[i][run] = (total_grid[1][p] - wave.atoms[asym_atom_list[i]].y);
                    d3[i][run] = (total_grid[2][p] - wave.atoms[asym_atom_list[i]].z);
                    run++;
                }
            }
            points += run;
            dens[i].resize(run);
            d1[i].resize(run);
            d2[i].resize(run);
            d3[i].resize(run);
        }
#pragma omp for
        for (int _grid = 0; _grid < total_grid.size(); _grid++)
            shrink_vector<double>(total_grid[_grid]);
    }
    shrink_vector<vec>(total_grid);
    return points;
}

/**
 * @brief Generates integration grids for SALTED method.
 *
 * This function generates integration grids for the SALTED method based on the provided parameters.
 *
 * @param accuracy The accuracy level of the integration grids.
 * @param unit_cell The unit cell object.
 * @param wave The WFN object.
 * @param coef_filename The filename of the coefficient file.
 * @param atom_type_list The list of atom types.
 * @param asym_atom_list The list of asymmetric atoms.
 * @param needs_grid The vector indicating whether each atom needs a grid.
 * @param d1 The vector of grid vectors d1.
 * @param d2 The vector of grid vectors d2.
 * @param d3 The vector of grid vectors d3.
 * @param dens The vector of density vectors.
 * @param exp_coefs The number of expansion coefficients.
 * @param file The output stream for writing results.
 * @param start The start time point.
 * @param end_becke The end time point for Becke grid generation.
 * @param end_prototypes The end time point for prototype grid generation.
 * @param end_prune The end time point for pruning grid points.
 * @param end_aspherical The end time point for aspherical density calculation.
 * @param debug Flag indicating whether to enable debug mode.
 * @param no_date Flag indicating whether to exclude date information in the output.
 *
 * @return The number of points in the integration grids.
 */
static int make_integration_grids_SALTED(
    const int &accuracy,
    cell &unit_cell,
    const WFN &wave,
    const string coef_filename,
    const vector<int> &atom_type_list,
    const vector<int> &asym_atom_list,
    vector<bool> &needs_grid,
    vector<vec> &d1,
    vector<vec> &d2,
    vector<vec> &d3,
    vector<vec> &dens,
    const int exp_coefs,
    ostream &file,
    time_point &start,
    time_point &end_becke,
    time_point &end_prototypes,
    time_point &end_prune,
    time_point &end_aspherical,
    bool debug,
    bool no_date)
{
    int atoms_with_grids = 0;
    for (int i = 0; i < needs_grid.size(); i++)
    {
        if (needs_grid[i])
            atoms_with_grids++;
    }
    // counts number of points inside each atomic grid
    ivec num_points(atoms_with_grids);
    // GRID COORDINATES for [a][c][p]
    // a = atom [0,ncen],
    // c = coordinate [0=x, 1=y, 2=z, 3=atomic weight],
    // p = point in this grid
    vector<vector<vec>> grid(atoms_with_grids);
#pragma omp parallel for
    for (int i = 0; i < atoms_with_grids; i++)
        grid[i].resize(4);

    const int nr_of_atoms = (wave.get_ncen());
    // positions
    vec x(nr_of_atoms), y(nr_of_atoms), z(nr_of_atoms);
    // charges
    ivec atom_z(nr_of_atoms);
    vec alpha_max(wave.get_ncen());
    ivec max_l(wave.get_ncen());
    int max_l_overall = 0;

    // Accumulate vectors with information about all atoms
#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        const atom *ai = &(wave.atoms[i]);
        atom_z[i] = wave.get_atom_charge(i);
        x[i] = ai->x;
        y[i] = ai->y;
        z[i] = ai->z;
        alpha_max[i] = 0.0;
        max_l[i] = 0;
        for (int b = 0; b < ai->basis_set.size(); b++)
        {
            if (ai->basis_set[b].exponent > alpha_max[i])
                alpha_max[i] = ai->basis_set[b].exponent * 2;
            int l = ai->basis_set[b].type;
            if (l > max_l[i])
            {
                max_l[i] = l;
#pragma omp critical
                {
                    if (l > max_l_overall)
                        max_l_overall = l;
                }
            }
        }
    }

    if (debug)
    {
        file << "Atoms are there! max_l:" << setw(5) << max_l_overall << endl;
        for (int i = 0; i < max_l.size(); i++)
            file << "max_l: " << setw(5) << max_l[i] << endl;
    }

    vector<vec> alpha_min(wave.get_ncen());
    for (int i = 0; i < wave.get_ncen(); i++)
        alpha_min[i].resize(max_l_overall, 100000000.0);

#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        for (int b = 0; b < max_l_overall; b++)
            alpha_min[i][b] = 100000000.0;
    }

#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        const atom *ai = &(wave.atoms[i]);
        for (int b = 0; b < ai->basis_set.size(); b++)
        {
            int l = ai->basis_set[b].type;
            if (ai->basis_set[b].exponent < alpha_min[i][l])
                alpha_min[i][l] = ai->basis_set[b].exponent;
        }
    }

    if (debug)
    {
        for (int i = 0; i < wave.get_ncen(); i++)
        {
            file << "alpha_min: ";
            for (int b = 0; b < max_l_overall; b++)
                file << setw(14) << scientific << alpha_min[i][b];
            file << endl;
        }
    }

    if (debug)
        file << "alpha_min is there!" << endl
             << "Nr of asym atoms: " << asym_atom_list.size() << " Number of atoms in wfn: " << wave.get_ncen() << " atoms that needs a grid: " << atoms_with_grids << endl;
    else
        file << "There are:\n"
             << setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
             //<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
             << setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << endl;

    if (no_date)
        file << "\nMaking Becke Grids..." << flush;
    else
    {
        if (debug)
            file << "max_l_overall: " << max_l_overall << endl
                 << "Selected accuracy: " << accuracy << "\nMaking Becke Grid for" << endl;
        else
            file << endl
                 << "Selected accuracy: " << accuracy << "\nMaking Becke Grids..." << flush;
    }

    // Make Prototype grids with only single atom weights for all elements
    vector<AtomGrid> Prototype_grids;

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        if (debug)
            file << "Atom Type " << i << ": " << atom_type_list[i] << endl;
        double alpha_max_temp(0);
        int max_l_temp(0);
        vec alpha_min_temp(max_l_overall);
        for (int j = 0; j < wave.get_ncen(); j++)
        {
            if (wave.get_atom_charge(j) == 119)
            {
                continue;
            }
            if (wave.get_atom_charge(j) == atom_type_list[i])
            {
                if (debug)
                {
                    file << alpha_max[j] << " " << max_l[j] - 1 << " ";
                    for (int l = 0; l < max_l_overall; l++)
                        file << alpha_min[j][l] << " ";
                    file << endl;
                }
                alpha_max_temp = alpha_max[j];
                max_l_temp = max_l[j];
                for (int l = 0; l <= max_l_temp; l++)
                    alpha_min_temp[l] = alpha_min[j][l];
                break;
            }
        }

        if (debug)
        {
            file << "max_l: " << defaultfloat << max_l_temp << " alpha_max: " << scientific << alpha_max_temp << " alpha_min: ";
            for (int l = 0; l <= max_l_temp; l++)
                file << setw(14) << scientific << alpha_min_temp[l];
            file << " accuracy: " << accuracy << endl;
        }
        int lebedev_high, lebedev_low;
        double radial_acc;
        err_checkf(accuracy >= 0, "Negative accuracy is not defined!", file);
        if (accuracy == 0)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                radial_acc = 1e-4;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[0] : constants::lebedev_table[1];
                radial_acc = 1e-3;
            }
        }
        else if (accuracy == 1)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[7] : constants::lebedev_table[8];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[3] : constants::lebedev_table[4];
                radial_acc = 1e-4;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[6] : constants::lebedev_table[7];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[2] : constants::lebedev_table[3];
                radial_acc = 1e-5;
            }
        }
        else if (accuracy == 2)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[11] : constants::lebedev_table[12];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[7] : constants::lebedev_table[8];
                radial_acc = 1e-5;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[10] : constants::lebedev_table[11];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[6] : constants::lebedev_table[7];
                radial_acc = 1e-6;
            }
        }
        else if (accuracy == 3)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[14] : constants::lebedev_table[16];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[12] : constants::lebedev_table[14];
                radial_acc = 1e-10;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[13] : constants::lebedev_table[15];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[11] : constants::lebedev_table[13];
                radial_acc = 1e-11;
            }
        }
        else if (accuracy == 4)
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[19] : constants::lebedev_table[21];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[14] : constants::lebedev_table[17];
                radial_acc = 1e-20;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[18] : constants::lebedev_table[20];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[13] : constants::lebedev_table[16];
                radial_acc = 1e-15;
            }
        }
        else
        {
            if (atom_type_list[i] != 1)
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[31] : constants::lebedev_table[32];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[29] : constants::lebedev_table[31];
                radial_acc = 1e-20;
            }
            else
            {
                lebedev_high = (max_l_temp < 3) ? constants::lebedev_table[30] : constants::lebedev_table[32];
                lebedev_low = (max_l_temp < 3) ? constants::lebedev_table[28] : constants::lebedev_table[30];
                radial_acc = 1e-15;
            }
        }
        Prototype_grids.push_back(AtomGrid(radial_acc,
                                           lebedev_low,
                                           lebedev_high,
                                           atom_type_list[i],
                                           alpha_max_temp,
                                           max_l_temp,
                                           alpha_min_temp.data(),
                                           file));
    }

    end_prototypes = get_time();
    if (debug)
    {

        for (int prototype = 0; prototype < Prototype_grids.size(); prototype++)
            file << "Number of gridpoints for atom type " << atom_type_list[prototype] << ": " << Prototype_grids[prototype].get_num_grid_points() << endl;

        int dur = get_sec(start, end_prototypes);
        if (dur < 1)
            file << "Time until prototypes are done: " << fixed << setprecision(0) << get_msec(start, end_prototypes) << " ms" << endl;
        else
            file << "Time until prototypes are done: " << fixed << setprecision(0) << dur << " s" << endl;
    }
    else
    {
        file << " ...  " << flush;
    }
    // get_grid is parallelized, therefore not parallel here
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        // skip atoms, that do not need a grid
        if (!needs_grid[i])
        {
            continue;
        }
        if (debug)
        {
            file << "Making grid for atom " << i << endl;
        }
        int type = 0;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.get_atom_charge(i))
                type = j;

        int grid_number = 0;
        for (int j = 0; j < i; j++)
        {
            if (needs_grid[j])
            {
                grid_number++;
            }
        }
        num_points[grid_number] = Prototype_grids[type].get_num_grid_points();

        for (int n = 0; n < grid[grid_number].size(); n++)
            grid[grid_number][n].resize(num_points[grid_number], 0.0);

        Prototype_grids[type].get_atomic_grid(
            i,
            &x[0],
            &y[0],
            &z[0],
            grid[grid_number][0].data(),
            grid[grid_number][1].data(),
            grid[grid_number][2].data(),
            grid[grid_number][3].data());
    }
    if (debug)
    {
        int run = 0;
        file << "  label | needs_grid | number of gridpoints\n";
        for (int j = 0; j < wave.get_ncen(); j++)
        {
            file << setw(8) << wave.atoms[j].label << setw(13) << needs_grid[j];
            if (needs_grid[j])
            {
                file << setw(7) << num_points[run];
                run++;
            }
            else
            {
                file << setw(6) << "---";
            }
            file << endl;
        }
    }
    Prototype_grids.clear();

    int points = 0;
    for (int i = 0; i < atoms_with_grids; i++)
        points += num_points[i];
    if (debug)
        file << "Becke Grid exists" << endl;
    else
        file << "                           done! Number of gridpoints: " << defaultfloat << points << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_becke = get_time();

    // Total grid as a sum of all atomic grids.
    // Dimensions: [c] [p]
    // p = the number of gridpoint
    // c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=atomic density
    vector<vec> total_grid(5);

    //int type_list_number = -1;

    double _cutoff;
    if (accuracy < 3)
        _cutoff = 1E-10;
    else if (accuracy == 3)
        _cutoff = 1E-14;
    else
        _cutoff = 1E-30;
    ivec new_gridsize(atoms_with_grids, 0);
    int final_size = 0;
    bool prune = false;
    if (prune)
    {
        file << "Pruning Grid..." << flush;
#pragma omp parallel for reduction(+ : final_size)
        for (int i = 0; i < atoms_with_grids; i++)
        {
            for (int p = 0; p < num_points[i]; p++)
            {
                if (abs(grid[i][3][p]) > _cutoff)
                {
                    new_gridsize[i]++;
                }
            }
            final_size += new_gridsize[i];
        }
        for (int k = 0; k < total_grid.size(); k++)
            total_grid[k].resize(final_size, 0.0);
#pragma omp parallel for
        for (int i = 0; i < atoms_with_grids; i++)
        {
            int offset = 0;
            for (int j = 0; j < i; j++)
            {
                offset += new_gridsize[j];
            }
            int reduction = 0;
            for (int p = 0; p < num_points[i]; p++)
            {
                if (abs(grid[i][3][p]) > _cutoff)
                {
                    for (int k = 0; k < 4; k++)
                        total_grid[k][p + offset - reduction] = grid[i][k][p];
                }
                else
                {
                    reduction++;
                }
            }
            num_points[i] -= reduction;
            shrink_vector<vec>(grid[i]);
        }
    }
    else
    {
        for (int i = 0; i < atoms_with_grids; i++)
        {
            for (int p = 0; p < num_points[i]; p++)
            {
                for (int k = 0; k < 4; k++)
                    total_grid[k].push_back(grid[i][k][p]);
                total_grid[4].push_back(0);
            }
            shrink_vector<vec>(grid[i]);
        }
    }
    shrink_vector<vector<vec>>(grid);
    points = 0;
    for (int i = 0; i < atoms_with_grids; i++)
        points += num_points[i];

    file << "                                       done! Number of gridpoints: " << defaultfloat << points << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_prune = get_time();

    file << "Calculating non-spherical densities..." << flush;

    WFN temp = wave;
    const int nr_pts = (int)total_grid[0].size();
    vector<unsigned long> shape{};
    bool fortran_order;
    vec data{};

    string path{coef_filename};
    npy::LoadArrayFromNumpy(path, shape, fortran_order, data);

    // #pragma omp parallel for
    //   for (int i = 0; i < nr_pts; i++) {
    //     total_grid[4][i] = calc_density_ML(
    //       total_grid[0][i],
    //       total_grid[1][i],
    //       total_grid[2][i],
    //       data,
    //       temp.atoms,
    //       exp_coefs,
    //
    //     );
    //   }

    if (debug)
        file << endl
             << "with total number of points: " << nr_pts << endl;
    else
        file << "                done!" << endl;

    file << "Applying weights and integrating charges...          " << flush;
    double el_sum_SALTED = 0.0;
    // Vector containing integrated numbers of electrons
    // dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
    // dimension 1: atoms of asym_atom_list
    vec atom_els;
    atom_els.resize(asym_atom_list.size(), 0.0);

    if (debug)
        file << "before loop" << endl;
        // Generate Electron sums
#pragma omp parallel for reduction(+ : el_sum_SALTED)
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int start_p = 0;
        for (int a = 0; a < i; a++)
            start_p += num_points[a];
        for (int p = start_p; p < start_p + num_points[i]; p++)
        {
            total_grid[4][p] = calc_density_ML(
                                   total_grid[0][i],
                                   total_grid[1][i],
                                   total_grid[2][i],
                                   data,
                                   temp.atoms,
                                   exp_coefs,
                                   asym_atom_list[i]) *
                               total_grid[3][p];
            atom_els[i] += total_grid[4][p];
        }
        el_sum_SALTED += atom_els[i];
        if (wave.get_has_ECPs())
        {
            int n = wave.atoms[asym_atom_list[i]].ECP_electrons;
            el_sum_SALTED += n;
            atom_els[i] += n;
        }
    }
    shrink_vector<double>(data);

    if (debug)
    {
        file << "Becke grid with hirshfeld weights done!" << endl;
        file << "atom_els[2]: ";
        for (int i = 0; i < asym_atom_list.size(); i++)
            file << fixed << setw(10) << setprecision(3) << atom_els[i] << " ";
        file << endl;
    }

    file << " done!\nNumber of points evaluated: " << total_grid[0].size();

    file << " with " << fixed << setw(10) << setprecision(6) << el_sum_SALTED << " electrons in Becke Grid in total.\n\n";

    file << "Table of Charges in electrons\n\n    Atom      Charge" << endl;

    int counter = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int a = asym_atom_list[i];
        file << setw(10) << wave.atoms[a].label
             << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[counter];
        counter++;
        file << endl;
    }

    file << "Total number of electrons in the wavefunction: " << el_sum_SALTED << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_aspherical = get_time();

    dens.resize(asym_atom_list.size());
    if (debug)
        file << "resized outer dens" << endl;
    d1.resize(asym_atom_list.size());
    d2.resize(asym_atom_list.size());
    d3.resize(asym_atom_list.size());
    if (debug)
        file << "resized outer d1-3" << endl;

    points = 0;
#pragma omp parallel
    {
#pragma omp for reduction(+ : points)
        for (int i = 0; i < asym_atom_list.size(); i++)
        {
            dens[i].resize(num_points[i]);
            d1[i].resize(num_points[i]);
            d2[i].resize(num_points[i]);
            d3[i].resize(num_points[i]);
            int start_p = 0;
            int run = 0;
            for (int a = 0; a < i; a++)
                start_p += num_points[a];
            for (int p = start_p; p < start_p + num_points[i]; p++)
            {
                if (abs(total_grid[4][p]) > _cutoff)
                {
                    dens[i][run] = (total_grid[4][p]);
                    d1[i][run] = (total_grid[0][p] - wave.atoms[asym_atom_list[i]].x);
                    d2[i][run] = (total_grid[1][p] - wave.atoms[asym_atom_list[i]].y);
                    d3[i][run] = (total_grid[2][p] - wave.atoms[asym_atom_list[i]].z);
                    run++;
                }
            }
            points += run;
            dens[i].resize(run);
            d1[i].resize(run);
            d2[i].resize(run);
            d3[i].resize(run);
        }
#pragma omp for
        for (int _grid = 0; _grid < total_grid.size(); _grid++)
            shrink_vector<double>(total_grid[_grid]);
    }
    shrink_vector<vec>(total_grid);
    return points;
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
void make_k_pts(const bool &read_k_pts,
                const bool &save_k_pts,
                const int gridsize,
                const cell &unit_cell,
                hkl_list &hkl,
                vector<vec> &k_pt,
                ostream &file,
                bool debug = false)
{
    const int size = (int)hkl.size();
    if (!read_k_pts)
    {
        k_pt.resize(3);
#pragma omp parallel for
        for (int i = 0; i < 3; i++)
            k_pt[i].resize(size, 0.0);

        if (debug)
            file << "K_point_vector is here! size: " << k_pt[0].size() << endl;

#pragma omp parallel for
        for (int ref = 0; ref < size; ref++)
        {
            hkl_list_it hkl_ = next(hkl.begin(), ref);
            for (int x = 0; x < 3; x++)
            {
                for (int j = 0; j < 3; j++)
                {
                    k_pt[x][ref] += unit_cell.get_rcm(x, j) * (*(hkl_))[j];
                }
            }
        }

        file << endl
             << "Number of k-points to evaluate: ";
        file << k_pt[0].size();
        file << " for " << gridsize << " gridpoints." << endl;
        if (save_k_pts)
            save_k_points(k_pt, hkl);
    }
    else
    {
        read_k_points(k_pt, hkl, file);
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
void calc_SF(const int &points,
             vector<vec> &k_pt,
             vector<vec> &d1,
             vector<vec> &d2,
             vector<vec> &d3,
             vector<vec> &dens,
             vector<cvec> &sf,
             ostream &file,
             time_point &start,
             time_point &end1,
             bool debug,
             bool no_date = false)
{
    const long long int imax = static_cast<long long int>(dens.size());
    const long long int step = std::max(static_cast<long long int>(std::floor(imax / 20)), 1LL);
    const long long int smax = static_cast<long long int>(k_pt[0].size());
    sf.reserve(imax * smax);
    sf.resize(imax);
#pragma omp parallel for
    for (int i = 0; i < imax; i++)
        sf[i].resize(smax, constants::cnull);

    if (debug)
        file << "Initialized FFs" << endl
             << "asym atom list size: " << imax << " total grid size: " << points << endl;
    end1 = get_time();
    if (!no_date)
    {
        int dur = get_sec(start, end1);
        if (dur < 1)
            file << "Time to prepare: " << fixed << setprecision(0) << get_msec(start, end1) << " ms" << endl;
        else
            file << "Time to prepare: " << fixed << setprecision(0) << dur << " s" << endl;
    }

#ifdef FLO_CUDA
    double **gpu_k_pt = NULL,
           **gpu_sf_r = NULL,
           **gpu_sf_i = NULL;
    vector<double> long_kpt;
    long_kpt.resize(3 * k_pt_unique[0].size());
    for (int i = 0; i < k_pt_unique[0].size(); i++)
    {
        long_kpt[3 * i + 0] = k_pt_unique[0][i];
        long_kpt[3 * i + 1] = k_pt_unique[1][i];
        long_kpt[3 * i + 2] = k_pt_unique[2][i];
    }
    gpu_k_pt = (double **)malloc(sizeof(double *));
    gpu_sf_r = (double **)malloc(sizeof(double *));
    gpu_sf_i = (double **)malloc(sizeof(double *));
    cudaMalloc((void **)&gpu_k_pt[0], sizeof(double) * k_pt_unique[0].size() * 3);
    cudaMalloc((void **)&gpu_sf_r[0][i], sizeof(double) * asym_atom_list.size() * k_pt_unique[0].size());
    cudaMalloc((void **)&gpu_sf_i[0][i], sizeof(double) * asym_atom_list.size() * k_pt_unique[0].size());
    cudaMemcpy(gpu_k_pt[0], long_kpt.data(), sizeof(double) * k_pt_unique[0].size() * 3, cudaMemcpyHostToDevice);

    dim3 blocks(asym_atom_list.size(), k_pt_unique[0].size());
    gpu_make_sf<<<blocks, 1>>>(
        gpu_sf_r[0],
        gpu_sf_i[0],
        gpu_k_pt[0],

    );
#else

    progress_bar *progress = new progress_bar{file, 60u, "Calculating scattering factors"};
    long long int pmax;
    double *dens_local, *d1_local, *d2_local, *d3_local;
    complex<double> *sf_local;
    const double *k1_local = k_pt[0].data();
    const double *k2_local = k_pt[1].data();
    const double *k3_local = k_pt[2].data();
    double work, rho;
    for (int i = 0; i < imax; i++)
    {
        pmax = static_cast<long long int>(dens[i].size());
        dens_local = dens[i].data();
        d1_local = d1[i].data();
        d2_local = d2[i].data();
        d3_local = d3[i].data();
        sf_local = sf[i].data();
#pragma omp parallel for private(work, rho)
        for (long long int s = 0; s < smax; s++)
        {
            for (long long int p = pmax - 1; p >= 0; p--)
            {
                rho = dens_local[p];
                work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
                if (rho < 0)
                {
                    rho = -rho;
                    work += M_PI;
                }
#endif
#endif
                sf_local[s] += polar(rho, work);
            }
        }
        if (i != 0 && i % step == 0)
            progress->write(i / static_cast<double>(imax));
    }
    delete (progress);

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
static void add_ECP_contribution(const vector<int> &asym_atom_list,
                                 const WFN &wave,
                                 vector<vector<complex<double>>> &sf,
                                 const cell &cell,
                                 hkl_list &hkl,
                                 ostream &file,
                                 const int &mode,
                                 const bool debug)
{
    double k = 1.0;
    hkl_list_it it = hkl.begin();
    if (mode == 0)
    { // Using a gaussian tight core function
        if (debug)
        {
            file << "Using a gaussian tight core function" << endl;
            for (int i = 0; i < asym_atom_list.size(); i++)
            {
                if (wave.atoms[asym_atom_list[i]].ECP_electrons != 0)
                    file << "Atom nr: " << wave.atoms[asym_atom_list[i]].charge << " core f000: "
                         << scientific << setw(14) << setprecision(8)
                         << wave.atoms[asym_atom_list[i]].ECP_electrons
                         << " and at 1 angstrom: " << exp(-pow(constants::bohr2ang(k), 2) / 16.0 / constants::PI) * wave.atoms[asym_atom_list[i]].ECP_electrons << endl;
            }
        }
#pragma omp parallel for private(it, k)
        for (int s = 0; s < sf[0].size(); s++)
        {
            it = next(hkl.begin(), s);
            k = constants::FOUR_PI * cell.get_stl_of_hkl(*it);
            for (int i = 0; i < asym_atom_list.size(); i++)
            {
                sf[i][s] += wave.atoms[asym_atom_list[i]].ECP_electrons * exp(-k / 16.0 / constants::PI);
            }
        }
    }
    else if (mode == 1 || mode == 2 || mode == 3)
    { // Using a Thakkar core density
        if (debug)
            file << "Using a Thakkar core density" << endl;
        vector<Thakkar> temp;
        for (int i = 0; i < asym_atom_list.size(); i++)
        {
            temp.push_back(Thakkar(wave.atoms[asym_atom_list[i]].charge, mode));
            if (debug && wave.atoms[asym_atom_list[i]].ECP_electrons != 0)
            {
                double k_0001 = temp[i].get_core_form_factor(0, wave.atoms[asym_atom_list[i]].ECP_electrons);
                double k_1 = temp[i].get_core_form_factor(constants::FOUR_PI * constants::bohr2ang(1.0), wave.atoms[asym_atom_list[i]].ECP_electrons);
                file << "Atom nr: " << wave.atoms[asym_atom_list[i]].charge << " number of ECP electrons: " << wave.atoms[asym_atom_list[i]].ECP_electrons << " core f(0) : "
                     << scientific << setw(14) << setprecision(8) << k_0001 << " and at 1 Ang: " << k_1 << endl;
            }
        }

#pragma omp parallel for private(it, k)
        for (int s = 0; s < sf[0].size(); s++)
        {
            it = next(hkl.begin(), s);
            k = constants::FOUR_PI * constants::bohr2ang(cell.get_stl_of_hkl(*it));
            for (int i = 0; i < asym_atom_list.size(); i++)
            {
                if (wave.atoms[asym_atom_list[i]].ECP_electrons != 0)
                    sf[i][s] += temp[i].get_core_form_factor(k, wave.atoms[asym_atom_list[i]].ECP_electrons);
            }
        }
    }
    else
    {
        err_not_impl_f("No higher ECP mode than 1 implemented!", file);
    }
}

/**
 * Converts the given asymmetric atom list, wavefunction, unit cell, and hkl list
 * to electron scattering factors using the given X-ray scattering factors.
 *
 * @param asym_atom_list The list of asymmetric atoms.
 * @param wave The wavefunction.
 * @param sf The scattering factors.
 * @param unit_cell The unit cell.
 * @param hkl The hkl list.
 */
void convert_to_ED(const std::vector<int> &asym_atom_list,
                   const WFN &wave,
                   std::vector<std::vector<std::complex<double>>> &sf,
                   const cell &unit_cell,
                   const hkl_list &hkl)
{
    double h2;
    hkl_list_it it;
#pragma omp parallel for private(h2, it)
    for (int s = 0; s < hkl.size(); s++)
    {
        it = next(hkl.begin(), s);
        h2 = pow(unit_cell.get_stl_of_hkl(*it), 2);
        for (int i = 0; i < asym_atom_list.size(); i++)
            sf[i][s] = std::complex<double>(constants::ED_fact * (wave.get_atom_charge(asym_atom_list[i]) - sf[i][s].real()) / h2, -constants::ED_fact * sf[i][s].imag() / h2);
    }
}

/**
 * Calculates the Thakkar scattering factors.
 *
 * This function calculates the Thakkar scattering factors based on the given options,
 * writes the results to the specified output stream, and modifies the provided WFN object.
 *
 * @param opt The options for calculating the scattering factors.
 * @param file The output stream to write the results to.
 * @param wave The WFN object to modify.
 * @return `true` if the calculation is successful, `false` otherwise.
 */
bool thakkar_sfac(
    const options &opt,
    ostream &file,
    WFN &wave)
{
    if (opt.hkl != "")
    {
        err_checkf(exists(opt.hkl), "HKL file does not exists!", file);
    }

    err_checkf(exists(opt.cif), "CIF does not exists!", file);
    file << "Number of protons: " << wave.get_nr_electrons() << endl;
    file << "Reading: " << opt.hkl;
    file.flush();

    cell unit_cell(opt.cif, file, opt.debug);

    ifstream cif_input(opt.cif.c_str(), std::ios::in);
    vector<int> atom_type_list;
    vector<int> asym_atom_to_type_list;
    vector<int> asym_atom_list;
    vector<bool> needs_grid(wave.get_ncen(), false);
    vector<string> known_atoms;

    read_atoms_from_CIF(cif_input,
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

    hkl_list hkl;
    if (!opt.read_k_pts)
    {
        if (opt.dmin != 99.0)
            if (opt.electron_diffraction)
                generate_hkl(opt.dmin / 2.0, hkl, opt.twin_law, unit_cell, file, opt.debug);
            else
                generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug);
        else
            read_hkl(opt.hkl, hkl, opt.twin_law, unit_cell, file, opt.debug);
    }

    if (opt.debug)
    {
        for (int i = 0; i < opt.Cations.size(); i++)
            file << "Cation: " << opt.Cations[i] << endl;
        for (int i = 0; i < opt.Anions.size(); i++)
            file << "Anion: " << opt.Anions[i] << endl;
    }
    vector<Thakkar> spherical_atoms;
    for (int i = 0; i < atom_type_list.size(); i++)
    {
        spherical_atoms.push_back(Thakkar(atom_type_list[i]));
    }
    // For all elements of Cations
    for (int i = 0; i < opt.Cations.size(); i++)
    {
        // look for atom that has matching label
        for (int j = 0; j < wave.atoms.size(); j++)
        {
            if (wave.atoms[j].label == opt.Cations[i])
            {
                // Look whether we already have this ion in our list
                int nr = -1;
                for (int k = 0; k < spherical_atoms.size(); k++)
                    if (spherical_atoms[k].get_atomic_number() == wave.atoms[j].charge && spherical_atoms[k].get_atomic_number() == 1)
                    {
                        nr = k;
                    }
                // If yes look for position in asym_atom_list
                if (nr > -1)
                {
                    for (int k = 0; k < asym_atom_list.size(); k++)
                    {
                        if (asym_atom_list[k] == j)
                        {
                            // and asign the nr in the type list accordingly.
                            asym_atom_to_type_list[k] = nr;
                        }
                    }
                }
                // If not append a new Cation
                else
                {
                    spherical_atoms.push_back(Thakkar_Cation(wave.atoms[j].charge));
                    atom_type_list.push_back(wave.atoms[j].charge);
                    nr = (int)spherical_atoms.size() - 1;
                    // and look for the new atom
                    for (int k = 0; k < asym_atom_list.size(); k++)
                    {
                        if (asym_atom_list[k] == j)
                        {
                            // And change the matching asym_atom_to_type entry
                            asym_atom_to_type_list[k] = nr;
                        }
                    }
                }
            }
        }
    }
    // For all elements of Anions
    for (int i = 0; i < opt.Anions.size(); i++)
    {
        // look for atom that has matching label
        for (int j = 0; j < wave.atoms.size(); j++)
        {
            if (wave.atoms[j].label == opt.Anions[i])
            {
                // Look whether we already have this ion in our list
                int nr = -1;
                for (int k = 0; k < spherical_atoms.size(); k++)
                    if (spherical_atoms[k].get_atomic_number() == wave.atoms[j].charge && spherical_atoms[k].get_atomic_number() == -1)
                    {
                        nr = k;
                    }
                // If yes look for position in asym_atom_list
                if (nr > -1)
                {
                    for (int k = 0; k < asym_atom_list.size(); k++)
                    {
                        if (asym_atom_list[k] == j)
                        {
                            // and asign the nr in the type list accordingly.
                            asym_atom_to_type_list[k] = nr;
                        }
                    }
                }
                // If not append a new Anion
                else
                {
                    spherical_atoms.push_back(Thakkar_Anion(wave.atoms[j].charge));
                    atom_type_list.push_back(wave.atoms[j].charge);
                    nr = (int)spherical_atoms.size() - 1;
                    // and look for the new atom
                    for (int k = 0; k < asym_atom_list.size(); k++)
                    {
                        if (asym_atom_list[k] == j)
                        {
                            // And change the matching asym_atom_to_type entry
                            asym_atom_to_type_list[k] = nr;
                        }
                    }
                }
            }
        }
    }
    if (opt.debug)
    {
        file << "AFTER ADDING IONS:" << endl;
        file << "There are " << atom_type_list.size() << " types of atoms" << endl;
        for (int i = 0; i < atom_type_list.size(); i++)
            file << setw(4) << atom_type_list[i];
        file << endl
             << "asym_atoms_to_type_list: " << endl;
        for (int i = 0; i < asym_atom_to_type_list.size(); i++)
            file << setw(4) << asym_atom_to_type_list[i];
        file << endl
             << "Charges of atoms:" << endl;
        for (int i = 0; i < wave.get_ncen(); i++)
            file << setw(4) << wave.get_atom_charge(i);
        file << endl
             << "Labels of atoms:" << endl;
        for (int i = 0; i < wave.get_ncen(); i++)
            file << setw(4) << wave.atoms[i].label;
        file << endl;
    }

    const int imax = (int)asym_atom_list.size();
    const int amax = (int)atom_type_list.size();

    file << "Calculating scattering factors for " << amax << " types of atoms with " << imax << " atoms in the asymmetric unit." << endl;

    vector<vector<double>> sf;
    sf.resize(asym_atom_list.size());
#pragma omp parallel for
    for (int i = 0; i < asym_atom_list.size(); i++)
        sf[i].resize(hkl.size());

    hkl_list_it it = hkl.begin();
#pragma omp parallel for private(it)
    for (int s = 0; s < hkl.size(); s++)
    {
        it = next(hkl.begin(), s);
        double k = constants::bohr2ang(constants::FOUR_PI * unit_cell.get_stl_of_hkl(*it));
        for (int i = 0; i < imax; i++)
            sf[i][s] = spherical_atoms[asym_atom_to_type_list[i]].get_form_factor(k);
    }

    if (opt.electron_diffraction)
    {
        double h2;
#pragma omp parallel for private(h2, it)
        for (int s = 0; s < hkl.size(); s++)
        {
            it = next(hkl.begin(), s);
            h2 = pow(unit_cell.get_stl_of_hkl(*it), 2);
            for (int i = 0; i < imax; i++)
                sf[i][s] = constants::ED_fact * (atom_type_list[asym_atom_to_type_list[i]] - sf[i][s]) / h2;
        }
    }

    if (opt.debug)
        file << endl
             << "SFs are made, now just write them!" << endl;
    file << endl
         << "Writing tsc file..." << endl;

    vector<string> labels;
    for (int i = 0; i < asym_atom_list.size(); i++)
        labels.push_back(wave.atoms[asym_atom_list[i]].label);

    tsc_block blocky(
        sf,
        labels,
        hkl);

    if (opt.binary_tsc)
        blocky.write_tscb_file();
    if (opt.old_tsc)
    {
        blocky.write_tsc_file(opt.cif);
    }
    file << " ... done!" << endl;
    return true;
}

/**
 * Calculates the scattering factors for multipel calcualtions happening after each other.
 *
 * @param opt The options object containing the necessary parameters.
 * @param file The output stream to write the results to.
 * @param known_atoms The vector of known atom names.
 * @param wave The vector of wave functions.
 * @param nr The number of radial points.
 * @return A tsc_block object containing the calculated scattering factors.
 */
tsc_block<int, cdouble> MTC_thakkar_sfac(
    options &opt,
    ostream &file,
    vector<string> &known_atoms,
    vector<WFN> &wave,
    const int &nr)
{
    err_checkf(exists(opt.hkl) || !(opt.dmin == 99.0), "HKL file does not exists!", file);
    err_checkf(exists(opt.cif), "CIF does not exists!", file);
    file << "Number of protons: " << wave[nr].get_nr_electrons() << endl;
    file << "Reading: " << opt.hkl;
    file.flush();

    cell unit_cell(opt.cif, file, opt.debug);
    ifstream cif_input(opt.cif.c_str(), std::ios::in);
    vector<int> atom_type_list;
    vector<int> asym_atom_to_type_list;
    vector<int> asym_atom_list;
    vector<bool> needs_grid(wave[nr].get_ncen(), false);

    read_atoms_from_CIF(cif_input,
                        opt.groups[nr],
                        unit_cell,
                        wave[nr],
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

    if (opt.debug)
        file << "made it post CIF, now make grids!" << endl;

    hkl_list hkl;
    if (opt.m_hkl_list.size() != 0)
    {
        hkl = opt.m_hkl_list;
        /*
#pragma omp parallel for
        for (int i = 0; i < known_indices[0].size(); i++) {
            ivectemp_hkl{ known_indices[0][i], known_indices[1][i], known_indices[2][i] };
#pragma omp critical
            hkl.emplace(temp_hkl);
        }*/
    }
    else if (nr == 0 && opt.read_k_pts == false)
    {
        if (opt.dmin != 99.0)
            if (opt.electron_diffraction)
                generate_hkl(opt.dmin / 2.0, hkl, opt.twin_law, unit_cell, file, opt.debug);
            else
                generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug);
        else
            read_hkl(opt.hkl, hkl, opt.twin_law, unit_cell, file, opt.debug);
        opt.m_hkl_list = hkl;
    }

    vector<Thakkar> spherical_atoms;
    for (int i = 0; i < atom_type_list.size(); i++)
        spherical_atoms.push_back(Thakkar(atom_type_list[i]));

    const int imax = (int)asym_atom_list.size();
    // const int amax = (int) atom_type_list.size();
    vector<vector<cdouble>> sf;
    sf.resize(asym_atom_list.size());
#pragma omp parallel for
    for (int i = 0; i < asym_atom_list.size(); i++)
        sf[i].resize(hkl.size());

    hkl_list_it it = hkl.begin();
#pragma omp parallel for private(it)
    for (int s = 0; s < hkl.size(); s++)
    {
        it = next(hkl.begin(), s);
        double k = constants::bohr2ang(constants::FOUR_PI * unit_cell.get_stl_of_hkl(*it));
        for (int i = 0; i < imax; i++)
            sf[i][s] = spherical_atoms[asym_atom_to_type_list[i]].get_form_factor(k);
    }

    if (opt.electron_diffraction)
    {
        double h2;
#pragma omp parallel for private(h2, it)
        for (int s = 0; s < hkl.size(); s++)
        {
            it = next(hkl.begin(), s);
            h2 = pow(unit_cell.get_stl_of_hkl(*it), 2);
            for (int i = 0; i < imax; i++)
                sf[i][s] = constants::ED_fact * ((cdouble)atom_type_list[i] - sf[i][s]) / h2;
        }
    }

    if (opt.debug)
        file << endl
             << "SFs are made, now just write them!" << endl;
    file << endl
         << "Writing tsc file..." << endl;

    vector<string> labels;
    for (int i = 0; i < asym_atom_list.size(); i++)
        labels.push_back(wave[nr].atoms[asym_atom_list[i]].label);

    tsc_block<int, cdouble> blocky(
        sf,
        labels,
        hkl);

    return blocky;
}

/**
 * Calculates the scattering factors for Hirshfeld partitioned molecules.
 *
 * @param opt The options for the calculation.
 * @param wave The WFN object representing the wavefunction.
 * @param file The output stream to write the results to.
 * @return True if the scattering factors were calculated successfully, false otherwise.
 */
bool calculate_scattering_factors_HF(
    const options &opt,
    WFN &wave,
    ostream &file)
{
#ifdef FLO_CUDA
    if (opt.pbc != 0)
    {
        file << "PERIODIC CALCULATIONS NOT IMPLEMENTED WITH CUDA YET!" << endl;
        exit(false);
    }
#endif
    err_checkf(wave.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", file);
    err_checkf(exists(opt.cif), "CIF does not exists!", file);
    file << "Number of protons: " << wave.get_nr_electrons() << endl
         << "Number of electrons: " << wave.count_nr_electrons() << endl;
    if (wave.get_has_ECPs())
        file << "Number of ECP electrons: " << wave.get_nr_ECP_electrons() << endl;
    // err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);

    time_point start = get_time();
    time_point end_becke, end_prototypes, end_spherical, end_prune, end_aspherical, end1;

    cell unit_cell(opt.cif, file, opt.debug);
    ifstream cif_input(opt.cif.c_str(), std::ios::in);
    vector<int> atom_type_list;
    vector<int> asym_atom_to_type_list;
    vector<int> asym_atom_list;
    vector<bool> needs_grid(wave.get_ncen(), false);
    vector<string> known_atoms;

    read_atoms_from_CIF(cif_input,
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

    if (opt.debug)
        file << "There are " << atom_type_list.size() << " Types of atoms and " << asym_atom_to_type_list.size() << " atoms in total" << endl;

    hkl_list hkl;
    if (!opt.read_k_pts)
    {
        if (opt.dmin != 99.0)
            if (opt.electron_diffraction)
                generate_hkl(opt.dmin / 2.0, hkl, opt.twin_law, unit_cell, file, opt.debug);
            else
                generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug);
        else
            read_hkl(opt.hkl, hkl, opt.twin_law, unit_cell, file, opt.debug);
    }

    if (opt.debug)
        file << "made it post CIF, now make grids!" << endl;
    vector<vec> d1, d2, d3, dens;

    int points = make_hirshfeld_grids(opt.pbc,
                                      opt.accuracy,
                                      unit_cell,
                                      wave,
                                      atom_type_list,
                                      asym_atom_list,
                                      needs_grid,
                                      d1, d2, d3, dens,
                                      file,
                                      start,
                                      end_becke,
                                      end_prototypes,
                                      end_spherical,
                                      end_prune,
                                      end_aspherical,
                                      opt.debug,
                                      opt.no_date);

    time_point before_kpts = get_time();

    vector<vec> k_pt;
    make_k_pts(
        opt.read_k_pts,
        opt.save_k_pts,
        points,
        unit_cell,
        hkl,
        k_pt,
        file,
        opt.debug);

    time_point after_kpts = get_time();

    vector<cvec> sf;
    calc_SF(points,
            k_pt,
            d1, d2, d3, dens,
            sf,
            file,
            start,
            end1,
            opt.debug,
            opt.no_date);

    if (wave.get_has_ECPs())
    {
        add_ECP_contribution(
            asym_atom_list,
            wave,
            sf,
            unit_cell,
            hkl,
            file,
            opt.ECP_mode,
            opt.debug);
    }

    if (opt.electron_diffraction)
    {
        convert_to_ED(asym_atom_list,
                      wave,
                      sf,
                      unit_cell,
                      hkl);
    }

    vector<string> labels;
    for (int i = 0; i < asym_atom_list.size(); i++)
        labels.push_back(wave.atoms[asym_atom_list[i]].label);

    tsc_block<int, cdouble> blocky(
        sf,
        labels,
        hkl);

    time_point end = get_time();
    if (!opt.no_date)
    {
        write_timing_to_file(file,
                             start,
                             end,
                             end_prototypes,
                             end_becke,
                             end_spherical,
                             end_prune,
                             end_aspherical,
                             before_kpts,
                             after_kpts,
                             end1);
    }

    file << "Writing tsc file... " << flush;
    blocky.write_tscb_file();
    if (opt.old_tsc)
    {
        blocky.write_tsc_file("test");
    }
    file << " ... done!" << endl;

#ifdef PEOJECT_NAME
#undef FLO_CUDA
#endif

    return true;
}

/**
 * Calculates the scattering factors for the given options, wave, and expansion coefficients.
 *
 * @param opt The options for the calculation.
 * @param wave The WFN object representing the wave.
 * @param file The output stream to write the results to.
 * @param exp_coefs The expected coefficients to use.
 * @return True if the scattering factors were calculated successfully, false otherwise.
 */
bool calculate_scattering_factors_RI(
    const options &opt,
    WFN &wave,
    ostream &file,
    const int exp_coefs)
{
    err_checkf(wave.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", file);
    err_checkf(exists(opt.cif), "CIF does not exists!", file);
    file << "Number of protons: " << wave.get_nr_electrons() << endl;
    // err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);

    time_point start = get_time();
    time_point end_becke, end_prototypes, end_spherical, end_prune, end_aspherical, end1;

    cell unit_cell(opt.cif, file, opt.debug);
    ifstream cif_input(opt.cif.c_str(), std::ios::in);
    vector<int> atom_type_list;
    vector<int> asym_atom_to_type_list;
    vector<int> asym_atom_list;
    vector<bool> needs_grid(wave.get_ncen(), false);
    vector<string> known_atoms;

    read_atoms_from_CIF(cif_input,
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

    if (opt.debug)
        file << "There are " << atom_type_list.size() << " Types of atoms and " << asym_atom_to_type_list.size() << " atoms in total" << endl;

    hkl_list hkl;
    if (!opt.read_k_pts)
    {
        if (opt.dmin != 99.0)
            if (opt.electron_diffraction)
                generate_hkl(opt.dmin / 2.0, hkl, opt.twin_law, unit_cell, file, opt.debug);
            else
                generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug);
        else
            read_hkl(opt.hkl, hkl, opt.twin_law, unit_cell, file, opt.debug);
    }

    if (opt.debug)
        file << "made it post CIF, now make grids!" << endl;
    vector<vec> d1, d2, d3, dens;

    int points = make_hirshfeld_grids_RI(
        opt.accuracy,
        wave,
        opt.coef_file,
        atom_type_list,
        asym_atom_list,
        needs_grid,
        d1, d2, d3, dens,
        exp_coefs,
        file,
        start,
        end_becke,
        end_prototypes,
        end_spherical,
        end_prune,
        end_aspherical,
        opt.debug,
        opt.no_date);

    time_point before_kpts = get_time();

    vector<vec> k_pt;
    make_k_pts(
        opt.read_k_pts,
        opt.save_k_pts,
        points,
        unit_cell,
        hkl,
        k_pt,
        file,
        opt.debug);

    time_point after_kpts = get_time();

    vector<cvec> sf;
    calc_SF(points,
            k_pt,
            d1, d2, d3, dens,
            sf,
            file,
            start,
            end1,
            opt.debug,
            opt.no_date);

    if (wave.get_has_ECPs())
    {
        add_ECP_contribution(
            asym_atom_list,
            wave,
            sf,
            unit_cell,
            hkl,
            file,
            opt.ECP_mode,
            opt.debug);
    }

    if (opt.electron_diffraction)
    {
        convert_to_ED(asym_atom_list,
                      wave,
                      sf,
                      unit_cell,
                      hkl);
    }

    vector<string> labels;
    for (int i = 0; i < asym_atom_list.size(); i++)
        labels.push_back(wave.atoms[asym_atom_list[i]].label);

    tsc_block<int, cdouble> blocky(
        sf,
        labels,
        hkl);

    time_point end = get_time();
    if (!opt.no_date)
    {
        write_timing_to_file(file,
                             start,
                             end,
                             end_prototypes,
                             end_becke,
                             end_spherical,
                             end_prune,
                             end_aspherical,
                             before_kpts,
                             after_kpts,
                             end1);
    }

    file << "Writing tsc file... " << flush;
    blocky.write_tscb_file();
    if (opt.old_tsc)
    {
        blocky.write_tsc_file("test");
    }
    file << " ... done!" << endl;

#ifdef PEOJECT_NAME
#undef FLO_CUDA
#endif

    return true;
}

bool calculate_scattering_factors_RI_No_H(
    const options &opt,
    WFN &wave,
    ostream &file,
    const int exp_coefs)
{
    err_checkf(wave.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", file);
    err_checkf(exists(opt.cif), "CIF does not exists!", file);
    file << "Number of protons: " << wave.get_nr_electrons() << endl;
    // err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);

    time_point start = get_time();
    time_point end_becke, end_prototypes, end_spherical, end_prune, end_aspherical, end1;

    cell unit_cell(opt.cif, file, opt.debug);
    ifstream cif_input(opt.cif.c_str(), std::ios::in);
    vector<int> atom_type_list;
    vector<int> asym_atom_to_type_list;
    vector<int> asym_atom_list;
    vector<bool> needs_grid(wave.get_ncen(), false);
    vector<string> known_atoms;

    read_atoms_from_CIF(cif_input,
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

    if (opt.debug)
        file << "There are " << atom_type_list.size() << " Types of atoms and " << asym_atom_to_type_list.size() << " atoms in total" << endl;

    hkl_list hkl;
    if (!opt.read_k_pts)
    {
        if (opt.dmin != 99.0)
            if (opt.electron_diffraction)
                generate_hkl(opt.dmin / 2.0, hkl, opt.twin_law, unit_cell, file, opt.debug);
            else
                generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug);
        else
            read_hkl(opt.hkl, hkl, opt.twin_law, unit_cell, file, opt.debug);
    }

    if (opt.debug)
        file << "made it post CIF, now make grids!" << endl;
    vector<vec> d1, d2, d3, dens;
    int points = 0;
    if (opt.SALTED)
    {
        file << "Making pure SALTED densities\n";
        points = make_integration_grids_SALTED(
            opt.accuracy,
            unit_cell,
            wave,
            opt.coef_file,
            atom_type_list,
            asym_atom_list,
            needs_grid,
            d1, d2, d3, dens,
            exp_coefs,
            file,
            start,
            end_becke,
            end_prototypes,
            end_prune,
            end_aspherical,
            opt.debug,
            opt.no_date);
    }
    else if (opt.SALTED_BECKE)
    {
        file << "Making SALTED densities and use Becke Integration\n";
        points = make_integration_grids(
            opt.accuracy,
            unit_cell,
            wave,
            opt.coef_file,
            atom_type_list,
            asym_atom_list,
            needs_grid,
            d1, d2, d3, dens,
            exp_coefs,
            file,
            start,
            end_becke,
            end_prototypes,
            end_spherical,
            end_prune,
            end_aspherical,
            opt.debug,
            opt.no_date);
    }
    else
        err_not_impl_f("No implementation of neither SALTED nor SALTED_BECKE", file);

    time_point before_kpts = get_time();

    vector<vec> k_pt;
    make_k_pts(
        opt.read_k_pts,
        opt.save_k_pts,
        points,
        unit_cell,
        hkl,
        k_pt,
        file,
        opt.debug);

    time_point after_kpts = get_time();

    vector<vector<complex<double>>> sf;
    calc_SF(points,
            k_pt,
            d1, d2, d3, dens,
            sf,
            file,
            start,
            end1,
            opt.debug,
            opt.no_date);

    if (wave.get_has_ECPs())
    {
        add_ECP_contribution(
            asym_atom_list,
            wave,
            sf,
            unit_cell,
            hkl,
            file,
            opt.ECP_mode,
            opt.debug);
    }

    if (opt.electron_diffraction)
    {
        convert_to_ED(asym_atom_list,
                      wave,
                      sf,
                      unit_cell,
                      hkl);
    }

    vector<string> labels;
    for (int i = 0; i < asym_atom_list.size(); i++)
        labels.push_back(wave.atoms[asym_atom_list[i]].label);

    tsc_block<int, cdouble> blocky(
        sf,
        labels,
        hkl);

    time_point end = get_time();
    if (!opt.no_date)
    {
        write_timing_to_file(file,
                             start,
                             end,
                             end_prototypes,
                             end_becke,
                             end_spherical,
                             end_prune,
                             end_aspherical,
                             before_kpts,
                             after_kpts,
                             end1);
    }

    file << "Writing tsc file... " << flush;
    blocky.write_tscb_file();
    if (opt.old_tsc)
    {
        blocky.write_tsc_file("test");
    }
    file << " ... done!" << endl;

#ifdef PEOJECT_NAME
#undef FLO_CUDA
#endif

    return true;
}

/**
 * Calculates the scattering factors for the MulTi Cif mode (MTC).
 *
 * @param opt The options for the calculation.
 * @param wave The vector of wave functions.
 * @param file The output stream to write the results to.
 * @param known_atoms The vector of known atoms.
 * @param nr The number of atoms.
 * @param kpts The pointer to the vector of k-points (optional).
 * @return The calculated scattering factors.
 */
tsc_block<int, cdouble> calculate_scattering_factors_MTC(
    options &opt,
    vector<WFN> &wave,
    ostream &file,
    vector<string> &known_atoms,
    const int &nr,
    vector<vec> *kpts)
{
#ifdef FLO_CUDA

#endif
    err_checkf(wave[nr].get_ncen() != 0, "No Atoms in the wavefunction, this will not work!!ABORTING!!", file);
    if (!opt.cif_based_combined_tsc_calc)
    {
        err_checkf(exists(opt.cif), "CIF " + opt.cif + " does not exists!", file);
    }
    else
    {
        for (int i = 0; i < opt.combined_tsc_calc_cifs.size(); i++)
        {
            err_checkf(exists(opt.combined_tsc_calc_cifs[i]), "CIF " + opt.combined_tsc_calc_cifs[i] + " does not exists!", file);
        }
    }
    err_checkf(opt.groups[nr].size() >= 1, "Not enough groups specified to work with!", file);
    file << "Number of protons: " << wave[nr].get_nr_electrons() << endl
         << "Number of electrons: " << wave[nr].count_nr_electrons() << endl;
    if (wave[nr].get_has_ECPs())
        file << "Number of ECP electrons: " << wave[nr].get_nr_ECP_electrons() << endl;
    // err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);
    if (opt.debug)
        file << "Working with: " << wave[nr].get_path() << endl;

    time_point start = get_time();
    time_point end_becke, end_prototypes, end_spherical, end_prune, end_aspherical, end1;

    string cif;
    if (opt.cif_based_combined_tsc_calc)
        cif = opt.combined_tsc_calc_cifs[nr];
    else
        cif = opt.cif;

    cell unit_cell(cif, file, opt.debug);
    ifstream cif_input(cif.c_str(), std::ios::in);
    vector<int> atom_type_list;
    vector<int> asym_atom_to_type_list;
    vector<int> asym_atom_list;
    vector<bool> needs_grid(wave[nr].get_ncen(), false);

    read_atoms_from_CIF(cif_input,
                        opt.groups[nr],
                        unit_cell,
                        wave[nr],
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

    if (opt.debug)
        file << "made it post CIF now make grids!" << endl;
    vector<vec> d1, d2, d3, dens;

    const int points = make_hirshfeld_grids(opt.pbc,
                                            opt.accuracy,
                                            unit_cell,
                                            wave[nr],
                                            atom_type_list,
                                            asym_atom_list,
                                            needs_grid,
                                            d1, d2, d3, dens,
                                            file,
                                            start,
                                            end_becke,
                                            end_prototypes,
                                            end_spherical,
                                            end_prune,
                                            end_aspherical,
                                            opt.debug,
                                            opt.no_date);

    time_point before_kpts = get_time();

    vector<vec> k_pt;
    hkl_list hkl;
    if (opt.m_hkl_list.size() != 0)
    {
        hkl = opt.m_hkl_list; /*
 #pragma omp parallel for
         for (int i = 0; i < known_indices[0].size(); i++) {
             ivectemp_hkl{ known_indices[0][i], known_indices[1][i], known_indices[2][i] };
 #pragma omp critical
             hkl.emplace(temp_hkl);
         }
         */
    }
    else if (nr == 0 && opt.read_k_pts == false)
    {
        if (opt.dmin != 99.0)
            if (opt.electron_diffraction)
                generate_hkl(opt.dmin / 2.0, hkl, opt.twin_law, unit_cell, file, opt.debug);
            else
                generate_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, file, opt.debug);
        else
            read_hkl(opt.hkl, hkl, opt.twin_law, unit_cell, file, opt.debug);
        opt.m_hkl_list = hkl;
    }
    if (kpts == NULL || kpts->size() == 0)
    {
        make_k_pts(
            nr != 0 && hkl.size() == 0,
            opt.save_k_pts,
            points,
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

    time_point after_kpts = get_time();

    vector<vector<complex<double>>> sf;
    calc_SF(points,
            k_pt,
            d1, d2, d3, dens,
            sf,
            file,
            start,
            end1,
            opt.debug,
            opt.no_date);

    if (wave[nr].get_has_ECPs())
    {
        add_ECP_contribution(
            asym_atom_list,
            wave[nr],
            sf,
            unit_cell,
            hkl,
            file,
            opt.ECP_mode,
            opt.debug);
    }

    if (opt.electron_diffraction)
    {
        convert_to_ED(asym_atom_list,
                      wave[nr],
                      sf,
                      unit_cell,
                      hkl);
    }

    vector<string> labels;
    for (int i = 0; i < asym_atom_list.size(); i++)
        labels.push_back(wave[nr].atoms[asym_atom_list[i]].label);

    tsc_block<int, cdouble> blocky(
        sf,
        labels,
        hkl);

    time_point end = get_time();
    if (!opt.no_date)
    {
        write_timing_to_file(file,
                             start,
                             end,
                             end_prototypes,
                             end_becke,
                             end_spherical,
                             end_prune,
                             end_aspherical,
                             before_kpts,
                             after_kpts,
                             end1);
    }

    return blocky;
}

/**
 * Calculates the diffuse (that is non integer hkl) scattering factors based on the given options.
 *
 * @param opt The options for calculating the diffuse scattering factors.
 * @param log_file The output stream for logging the calculation process.
 */
void calc_sfac_diffuse(const options &opt, std::ostream &log_file)
{
    using namespace std;
    std::vector<WFN> wavy;
    wavy.emplace_back(1);
    wavy[0].read_known_wavefunction_format(opt.wfn, std::cout, opt.debug);
    //set number of threads
#ifdef _OPENMP
    if (opt.threads > 0)
        omp_set_num_threads(opt.threads);
#endif
    err_checkf(wavy[0].get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
    err_checkf(exists(opt.cif), "CIF does not exists!", std::cout);
    // err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);

    time_point start = get_time();
    time_point end_becke, end_prototypes, end_spherical, end_prune, end_aspherical;

    cell unit_cell(opt.cif, std::cout, opt.debug);
    ifstream cif_input(opt.cif.c_str(), std::ios::in);
    vector<int> atom_type_list;
    vector<int> asym_atom_to_type_list;
    vector<int> asym_atom_list;
    vector<bool> needs_grid(wavy[0].get_ncen(), false);
    vector<string> known_atoms;

    read_atoms_from_CIF(cif_input,
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
    vector<vec> d1, d2, d3, dens;

    make_hirshfeld_grids(opt.pbc,
                         opt.accuracy,
                         unit_cell,
                         wavy[0],
                         atom_type_list,
                         asym_atom_list,
                         needs_grid,
                         d1, d2, d3, dens,
                         std::cout,
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
                         opt.debug,
                         opt.no_date);

    hkl_list_d hkl;
    generate_fractional_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, log_file, opt.sfac_diffuse, opt.debug);

    const long long int size = static_cast<long long int>(hkl.size());
    vector<vec> k_pt;
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
    for (const vec &hkl_ : hkl)
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
    vector<vector<complex<double>>> sf;

    const int imax = static_cast<int>(dens.size());
    const long long int smax = static_cast<long long int>(k_pt[0].size());
    long long int pmax = static_cast<long long int>(dens[0].size());
    const int step = max(static_cast<long long int>(floor(smax / 20)), 1LL);
    std::cout << "Done with making k_pt " << smax << " " << imax << " " << pmax << endl;
		sf.reserve(imax * smax);
    sf.resize(imax);
#pragma omp parallel for
    for (int i = 0; i < imax; i++)
        sf[i].resize(smax);
    double *dens_local, *d1_local, *d2_local, *d3_local;
    complex<double> *sf_local;
    const double *k1_local = k_pt[0].data();
    const double *k2_local = k_pt[1].data();
    const double *k3_local = k_pt[2].data();
    double work, rho;
    progress_bar *progress = new progress_bar{std::cout, 60u, "Calculating scattering factors"};
    for (int i = 0; i < imax; i++)
    {
        pmax = static_cast<long long int>(dens[i].size());
        dens_local = dens[i].data();
        d1_local = d1[i].data();
        d2_local = d2[i].data();
        d3_local = d3[i].data();
        sf_local = sf[i].data();
#pragma omp parallel for private(work, rho)
        for (long long int s = 0; s < smax; s++)
        {
            for (long long int p = pmax - 1; p >= 0; p--)
            {
              rho = dens_local[p];
              work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
              if (rho < 0)
              {
                rho = -rho;
                work += M_PI;
              }
#endif
#endif
              sf_local[s] += polar(rho, work);
            }
            if (i != 0 && i % step == 0)
                progress->write(i / static_cast<double>(imax));
        }
    }
    delete (progress);
    vector<string> labels;
    for (int i = 0; i < asym_atom_list.size(); i++)
        labels.push_back(wavy[0].atoms[asym_atom_list[i]].label);
    tsc_block<double, cdouble> result(sf, labels, hkl);
    result.write_tsc_file_non_integer(opt.cif);
}