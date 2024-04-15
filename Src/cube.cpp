#include "cube.h"
#include "convenience.h"

using namespace std;

cube::cube()
{
    size.resize(3, 0);
    origin.resize(3, 0.0);
    vectors.resize(3);
    for (int i = 0; i < 3; i++)
    {
        vectors[i].resize(3, 0.0);
    }
    loaded = false;
    na = 0;
    parent_wavefunction = new WFN(6);
};

cube::cube(int x, int y, int z, int g_na, bool grow_values)
{
    size.push_back(x);
    size.push_back(y);
    size.push_back(z);
    origin.resize(3, 0.0);
    vectors.resize(3);
    for (int i = 0; i < 3; i++)
        vectors[i].resize(3, 0.0);
    loaded = grow_values;
    if (grow_values)
    {
        values.resize(x);
#pragma omp parallel for
        for (int i = 0; i < x; i++)
        {
            values[i].resize(y);
            for (int j = 0; j < y; j++)
                values[i][j].resize(z, 0.0);
        }
    }
    na = g_na;
    parent_wavefunction = new WFN(6);
};

cube::cube(const string &filepath, bool read, WFN &wave, ostream &file, bool expert)
{
    err_checkf(exists(filepath), "Sorry, this file does not exist!", file);
    parent_wavefunction = &wave;
    na = parent_wavefunction->get_ncen();
    path = filepath;
    err_checkf(read_file(read, true, expert), "Sorry, something went wrong while reading!", file);
    loaded = read;
};

cube::cube(int g_na, const ivec &g_size, const vector<double> &g_origin, const vector<vec> &g_vectors, const vector<vector<vec>> &g_values)
{
    na = g_na;
    parent_wavefunction = new WFN(6);
    cout << "Assigned Nr of Atoms" << endl;
    vectors.resize(3);
    size.resize(3);
    for (int i = 0; i < 3; i++)
    {
        cout << i << ". dimension" << endl;
        size[i] = g_size[i];
        origin.push_back(g_origin[i]);
        for (int j = 0; j < 3; j++)
            vectors[i].push_back(g_vectors[i][j]);
    }
    values.resize(size[0]);
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
    {
        values[x].resize(size[1]);
        for (int y = 0; y < size[1]; y++)
        {
            values[x][y].resize(size[2]);
            for (int z = 0; z < size[2]; z++)
                values[x][y][z] = g_values[x][y][z];
        }
    }
    loaded = true;
};

cube::cube(const cube &given)
{
    na = given.get_na();
    path = given.path;
    comment1 = given.get_comment1();
    comment2 = given.get_comment2();
    vectors.resize(3);
    size.resize(3);
    for (int i = 0; i < 3; i++)
    {
        size[i] = given.get_size(i);
        origin.push_back(given.get_origin(i));
        for (int j = 0; j < 3; j++)
            vectors[i].push_back(given.get_vector(i, j));
    }
    loaded = given.get_loaded();
    if (loaded)
    {
        values.resize(size[0]);
#pragma omp parallel for
        for (int x = 0; x < size[0]; x++)
        {
            values[x].resize(size[1]);
            for (int y = 0; y < size[1]; y++)
            {
                values[x][y].resize(size[2]);
                for (int z = 0; z < size[2]; z++)
                    values[x][y][z] = given.get_value(x, y, z);
            }
        }
    }
};

bool cube::read_file(bool full, bool header, bool expert)
{
    ifstream file(path.c_str());
    string line;
    if (header)
    {
        getline(file, comment1);
        getline(file, comment2);
        getline(file, line);
        origin.resize(3);
        sscanf(line.c_str(), "%d %lf %lf %lf", &na, &origin[0], &origin[1], &origin[2]); // na=number of atoms in the cube file
        size.resize(3);
        vectors.resize(3);
        for (int i = 0; i < 3; i++)
            vectors[i].resize(3);
        getline(file, line);
        sscanf(line.c_str(), "%d%lf%lf%lf", &size[0], &vectors[0][0], &vectors[0][1], &vectors[0][2]);
        getline(file, line);
        sscanf(line.c_str(), "%d%lf%lf%lf", &size[1], &vectors[1][0], &vectors[1][1], &vectors[1][2]);
        getline(file, line);
        sscanf(line.c_str(), "%d%lf%lf%lf", &size[2], &vectors[2][0], &vectors[2][1], &vectors[2][2]);
        int atnr;
        double atp[3] = {0, 0, 0};
        bool read_atoms = true;
        if (!expert && parent_wavefunction->get_ncen() != 0)
            read_atoms = false;
        if (read_atoms)
            for (int i = 0; i < na; i++)
            { // loop of the positions of the atoms to get to the start of the values
                getline(file, line);
                double dum;
                sscanf(line.c_str(), "%d %lf %lf %lf %lf", &atnr, &dum, &atp[0], &atp[1], &atp[2]);
                parent_wavefunction->push_back_atom(atnr2letter(atnr), atp[0], atp[1], atp[2], atnr);
            }
    }
    if (full)
    {
        file.seekg(0);
        for (int i = 0; i < na + 6; i++)
            getline(file, line);
        values.resize(size[0]);
        for (int i = 0; i < size[0]; i++)
        {
            values[i].resize(size[1]);
            for (int j = 0; j < size[1]; j++)
                values[i][j].resize(size[2]);
        }
        int reads2 = 0;
        int reads1 = 0;
        int rest2 = 0;
        int run_x = 0;
        int run_y = 0;
        double tmp[6] = {0, 0, 0, 0, 0, 0};
        while (run_x < size[0] && run_y < size[1] && !file.eof())
        {
            reads2 = rest2;
            while (reads2 < size[2] && !file.eof())
            {
                getline(file, line);
                reads1 = sscanf(line.c_str(), "%lf%lf%lf%lf%lf%lf",
                                &tmp[0],
                                &tmp[1],
                                &tmp[2],
                                &tmp[3],
                                &tmp[4],
                                &tmp[5]);
                reads2 += reads1;
                rest2 = reads2 - size[2];
                if (rest2 < 6 && rest2 > 0)
                {
                    for (int j = 0; j < 6 - rest2; j++)
                        values[run_x][run_y][size[2] - (6 - rest2) + j] = tmp[j];
                    if (run_y + 1 < size[1])
                        for (int j = 0; j < rest2; j++)
                            values[run_x][run_y + 1][j] = tmp[j + (6 - rest2)];
                    else if (run_x + 1 < size[0])
                        for (int j = 0; j < rest2; j++)
                            values[run_x + 1][0][j] = tmp[j + (6 - rest2)];
                    else
                    {
                        cout << "This should not happen! Read a value outside of range!"
                             << " Run_x: " << run_x
                             << " Run_y: " << run_y
                             << " reads1: " << reads1
                             << " reads2: " << reads2 << endl;
                        return (false);
                    }
                }
                else
                    for (int i = 0; i < reads1; i++)
                        values[run_x][run_y][reads2 - reads1 + i] = tmp[i];
            }
            run_y++;
            if (run_y == size[1])
            {
                run_x++;
                if (run_x != size[0])
                    run_y = 0;
            }
        }
        if (run_x != size[0] || run_y != size[1])
        {
            cout << "This file ended before i read all expected values!" << endl;
            if (file.eof())
                cout << "ENCOUNTERED EOF!" << endl;
            cout << "x,y,reads1,reads2: " << run_x << " " << run_y << " " << reads1 << "," << reads2 << endl;
            return (false);
        }
        file.close();
        loaded = true;
    }
    return (true);
}

bool cube::write_file(bool force, bool absolute)
{
    if (exists(path))
    {
        if (force)
            remove(path.c_str());
        else
        {
            cout << "File already exists, aborting!" << endl;
            return false;
        }
    }
    /*bool end = false;
    if (!force) {
      while (exists(path) && !end) {
        cout << "File already exists, do you want to overwrite it? ";
        if (!yesno()) {
          cout << "Do you want to give a new name? ";
          if (yesno()) {
            cout << "Please give the new path where it should be saved: ";
            cin >> path;
          }
          else {
            cout << "okay, canceling!" << endl;
            return (false);
          }
        }
        else end = true;
      }
    }
    else {
      if (exists(path))
        remove(path.c_str());
    }*/
    stringstream stream;
    string temp;
    ofstream of(path.c_str(), ios::out);
    of << comment1 << endl;
    of << comment2 << endl;
    of << setw(5) << na << fixed << setw(12) << setprecision(6) << origin[0] << fixed << setw(12) << setprecision(6) << origin[1] << fixed << setw(12) << setprecision(6) << origin[2] << endl;
    for (int i = 0; i < 3; i++)
    {
        stream << setw(5) << size[i];
        for (int j = 0; j < 3; j++)
            stream << fixed << setw(12) << setprecision(6) << vectors[i][j];
        stream << endl;
    }
    of << stream.str();
    stream.str("");
    for (int i = 0; i < na; i++)
    {
        of << setw(5) << parent_wavefunction->get_atom_charge(i) << setw(5) << parent_wavefunction->get_atom_charge(i) << ".000000";
        for (int j = 0; j < 3; j++)
            of << fixed << setw(12) << setprecision(6) << parent_wavefunction->get_atom_coordinate(i, j);
        of << endl;
    }
    for (int run_x = 0; run_x < size[0]; run_x++)
    {
        for (int run_y = 0; run_y < size[1]; run_y++)
        {
            int temp_write = 0;
            while (temp_write < size[2])
            {
                if (absolute)
                    stream << uppercase << scientific << setw(13) << setprecision(5) << abs(values[run_x][run_y][temp_write]);
                else
                    stream << uppercase << scientific << setw(13) << setprecision(5) << values[run_x][run_y][temp_write];
                temp_write++;
                if (temp_write % 6 == 0)
                    stream << endl;
            }
            if (temp_write % 6 != 0)
                stream << endl;
        }
        of << stream.str();
        stream.str("");
    }
    return (true);
};

bool cube::write_file(string &given_path, bool debug)
{
    stringstream stream;
    string temp;
    ofstream of(given_path.c_str(), ios::out);
    of << comment1 << "\n";
    of << comment2 << "\n";
    of << setw(5) << na << fixed << setw(12) << setprecision(6) << origin[0] << fixed << setw(12) << setprecision(6) << origin[1] << fixed << setw(12) << setprecision(6) << origin[2] << "\n";
    for (int i = 0; i < 3; i++)
    {
        of << setw(5) << size[i];
        for (int j = 0; j < 3; j++)
            of << fixed << setw(12) << setprecision(6) << vectors[i][j];
        of << "\n";
    }
    for (int i = 0; i < na; i++)
    {
        of << setw(5) << parent_wavefunction->get_atom_charge(i) << setw(5) << parent_wavefunction->get_atom_charge(i) << ".000000";
        for (int j = 0; j < 3; j++)
            of << fixed << setw(12) << setprecision(6) << parent_wavefunction->get_atom_coordinate(i, j);
        of << "\n";
    }
    if (debug)
        cout << "Finished atoms!" << endl;
    if (get_loaded())
        for (int run_x = 0; run_x < size[0]; run_x++)
            for (int run_y = 0; run_y < size[1]; run_y++)
            {
                int temp_write = 0;
                while (temp_write < size[2])
                {
                    of << uppercase << scientific << setw(13) << setprecision(5) << values[run_x][run_y][temp_write];
                    temp_write++;
                    if (temp_write % 6 == 0)
                        of << "\n";
                }
                if (debug)
                    cout << "Write Z-line!" << endl;
                if (temp_write % 6 != 0)
                    of << "\n";
            }
    else
    {
        ifstream f(path, ios::in);
        string line_buffer;
        for (int a = 0; a < na + 6; a++)
            getline(f, line_buffer);
        while (!f.eof())
        {
            getline(f, line_buffer);
            of << line_buffer << "\n";
        }
    }
    // of << stream.str();
    of.flush();
    of.close();
    path = given_path;
    return (true);
};

bool cube::write_xdgraph(string &given_path, bool debug)
{
    stringstream stream;
    string temp;
    ofstream of(given_path.c_str(), ios::out);
    of << "2DGRDFIL  0" << endl
       << "cuQCT    FOU" << endl
       << endl;
    of << "! Gridpoints, Origin, Physical Dimensions" << endl;
    for (int i = 0; i < 3; i++)
        of << setw(14) << size[i];
    of << endl;
    of << fixed << setw(11) << setprecision(4) << origin[0] << " " << origin[1] << " " << origin[2] << endl;
    for (int j = 0; j < 3; j++)
        of << fixed << setw(12) << setprecision(6) << sqrt(pow(vectors[j][0], 2) + pow(vectors[j][1], 2) + pow(vectors[j][2], 2));
    of << endl;
    of << "! Objects" << endl;
    of << setw(10) << na << endl;
    for (int i = 0; i < na; i++)
    {
        of << parent_wavefunction->atoms[i].label;
        of << "    ";
        for (int j = 0; j < 3; j++)
            of << setw(10) << setprecision(5) << parent_wavefunction->get_atom_coordinate(i, j);
        of << " ATOM" << endl;
    }
    if (debug)
        cout << "Finished atoms!" << endl;
    of << "! Connections" << endl
       << "         0" << endl
       << "! Values" << endl;
    if (!get_loaded())
        read_file(true, false, false);
    for (int run_x = 0; run_x < size[0]; run_x++)
        for (int run_z = 0; run_z < size[2]; run_z++)
        {
            int temp_write = 0;
            while (temp_write < size[1])
            {
                of << uppercase << scientific << setw(15) << setprecision(7) << values[run_x][temp_write][run_z];
                temp_write++;
                if (temp_write % 6 == 0)
                    of << endl;
            }
            if (debug)
                cout << "Write Y-line!" << endl;
            if (temp_write % 6 != 0)
                of << endl;
        }
    of.flush();
    of.close();
    path = given_path;
    return (true);
};

bool cube::fractal_dimension(const double stepsize)
{
    double min = 100, max = -100;
    for (int i = 0; i < size[0]; i++)
        for (int j = 0; j < size[1]; j++)
            for (int k = 0; k < size[2]; k++)
            {
                if (values[i][j][k] < min)
                    min = values[i][j][k];
                if (values[i][j][k] > max)
                    max = values[i][j][k];
            }
    const double map_min = min, map_max = max;
    vec e = double_sum();
    min -= 2 * stepsize, max += 2 * stepsize;
    const int steps = int((max - min) / stepsize) + 2;
    ivec bins;
    vec df;
    vec iso;
    bins.resize(steps), df.resize(steps), iso.resize(steps);
    for (int i = 0; i < steps; i++)
        iso[i] = round((min + i * stepsize) * 100) / 100;
    const int comparisons = size[0] * size[1] * (size[2] - 1) + size[0] * size[2] * (size[1] - 1) + size[2] * size[1] * (size[0] - 1);
    double lv1, lv2;
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2] - 1; z++)
            {
                lv1 = values[x][y][z];
                lv2 = values[x][y][z + 1];
#pragma omp parallel for
                for (int i = 0; i < steps; i++)
                    if ((lv1 < iso[i] && lv2 > iso[i]) || (lv1 > iso[i] && lv2 < iso[i]))
                        bins[i]++;
            }
    for (int z = 0; z < size[2]; z++)
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1] - 1; y++)
            {
                lv1 = values[x][y][z];
                lv2 = values[x][y + 1][z];
#pragma omp parallel for
                for (int i = 0; i < steps; i++)
                    if ((lv1 < iso[i] && lv2 > iso[i]) || (lv1 > iso[i] && lv2 < iso[i]))
                        bins[i]++;
            }
    for (int y = 0; y < size[1]; y++)
        for (int z = 0; z < size[2]; z++)
            for (int x = 0; x < size[0] - 1; x++)
            {
                lv1 = values[x][y][z];
                lv2 = values[x + 1][y][z];
#pragma omp parallel for
                for (int i = 0; i < steps; i++)
                    if ((lv1 < iso[i] && lv2 > iso[i]) || (lv1 > iso[i] && lv2 < iso[i]))
                        bins[i]++;
            }
    const double third = -1.0 / 3.0;
    const double epsilon = log(1 / (pow(comparisons, third)));
    for (int i = 0; i < steps; i++)
    {
        if (bins[i] == 0)
            df[i] = 0.0;
        else
            df[i] = log(bins[i]) / epsilon;
    }
    string output(path + "_fractal_plot");
    ofstream of(output.c_str(), ios::out);
    of << setw(8) << scientific << setprecision(8) << steps
       << setw(16) << scientific << setprecision(8) << map_min
       << setw(16) << scientific << setprecision(8) << map_max
       << setw(16) << scientific << setprecision(8) << e[0] * pow(0.529177249, 3)
       << setw(16) << scientific << setprecision(8) << e[1] * pow(0.529177249, 3) << "\n";
    for (int i = 0; i < steps; i++)
        of << setw(16) << scientific << setprecision(8) << iso[i] << setw(16) << scientific << setprecision(8) << df[i] << "\n";
    of.flush();
    of.close();
    return true;
}

int cube::get_size(int direction) const
{
    if (direction < size.size() && direction >= 0)
        return (size[direction]);
    else
        return (-1);
};

double cube::get_value(int x, int y, int z) const
{
    if (x < size[0] && y < size[1] && z < size[2] && x >= 0 && y >= 0 && z >= 0)
        return (values[x][y][z]);
    else
        return (-1);
};

bool cube::set_value(int x, int y, int z, double value)
{
    if (x < size[0] && y < size[1] && z < size[2] && x >= 0 && y >= 0 && z >= 0)
        values[x][y][z] = value;
    else
        return (false);
    return (true);
};

void cube::operator=(cube &right)
{
    size.resize(3);
    for (int i = 0; i < 3; i++)
        size[i] = right.get_size(i);
    comment1 = right.get_comment1();
    comment2 = right.get_comment2();
    na = right.get_na();
    origin.resize(3);
    vectors.resize(3);
    for (int i = 0; i < 3; i++)
        vectors[i].resize(3);
    values.resize(size[0]);
#pragma omp parallel for
    for (int i = 0; i < size[0]; i++)
    {
        values[i].resize(size[1]);
        for (int j = 0; j < size[1]; j++)
            values[i][j].resize(size[2]);
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            vectors[i][j] = right.get_vector(i, j);
    }
    if (right.get_loaded())
    {
#pragma omp parallel for
        for (int x = 0; x < size[0]; x++)
        {
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    values[x][y][z] = right.get_value(x, y, z);
        }
        loaded = true;
    }
    parent_wavefunction = right.parent_wavefunction;
};

cube cube::operator+(cube &right) const
{
    cube res_cube(path, true, *parent_wavefunction, cout, false);
    res_cube.path = get_foldername_from_path(path) + get_filename_from_path(path).substr(0, get_filename_from_path(path).rfind(".cub")) + "+" + get_filename_from_path(right.path).substr(0, get_filename_from_path(right.path).rfind(".cub")) + ".cube";
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (cube());
    if (right.get_loaded())
#pragma omp parallel for
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    res_cube.set_value(x, y, z, res_cube.get_value(x, y, z) + right.get_value(x, y, z));
    else
    {
        int reads2 = 0;
        int reads1 = 0;
        int rest2 = 0;
        int run_x = 0;
        int run_y = 0;
        string line2;
        ifstream file2(right.path.c_str(), ios::in);
        double tmp[6]{0, 0, 0, 0, 0, 0};
        for (int i = 0; i < 6 + na; i++)
            if (!right.get_loaded())
                getline(file2, line2);
        while (run_x < size[0] && run_y < size[1] && !file2.eof())
        {
            reads2 = rest2;
            while (reads2 < size[2] && !file2.eof())
            {
                if (!right.get_loaded())
                {
                    getline(file2, line2);
                    reads1 = sscanf(line2.c_str(), "%lf%lf%lf%lf%lf%lf",
                                    &tmp[0],
                                    &tmp[1],
                                    &tmp[2],
                                    &tmp[3],
                                    &tmp[4],
                                    &tmp[5]);
                    reads2 += reads1;
                    rest2 = reads2 - size[2];
                    if (rest2 < 6 && rest2 > 0)
                    {
                        for (int j = 0; j < 6 - rest2; j++)
                            res_cube.set_value(run_x, run_y, size[2] - (6 - rest2) + j, res_cube.get_value(run_x, run_y, size[2] - (6 - rest2) + j) + tmp[j]);
                        if (run_y + 1 < size[1])
                            for (int j = 0; j < rest2; j++)
                                res_cube.set_value(run_x, run_y + 1, j, res_cube.get_value(run_x, run_y + 1, j) + tmp[(6 - rest2) + j]);
                        else if (run_x + 1 < size[0])
                            for (int j = 0; j < rest2; j++)
                                res_cube.set_value(run_x + 1, 0, j, res_cube.get_value(run_x + 1, 0, j) + tmp[(6 - rest2) + j]);
                        else
                        {
                            cout << "This should not happen! Read a value outside of range!"
                                 << " Run_x: " << run_x
                                 << " Run_y: " << run_y
                                 << " reads1: " << reads1
                                 << " reads2: " << reads2 << endl;
                            return (cube());
                        }
                    }
                    else
                        for (int i = 0; i < reads1; i++)
                            res_cube.set_value(run_x, run_y, reads2 - reads1 + i, res_cube.get_value(run_x, run_y, reads2 - reads1 + i) + tmp[i]);
                }
                else
                    for (int z = 0; z < size[2]; z++)
                        res_cube.set_value(run_x, run_y, z, res_cube.get_value(run_x, run_y, z) + right.get_value(run_x, run_y, z));
            }
            run_y++;
            if (run_y == size[1])
            {
                run_x++;
                if (run_x != size[0])
                    run_y = 0;
            }
        }
        if (run_x != size[0] || run_y != size[1])
        {
            cout << "This file ended before i read all expected values!" << endl;
            if (file2.eof())
                cout << "ENCOUNTERED EOF!" << endl;
            cout << "x,y,reads1,reads2: " << run_x << " " << run_y << " " << reads1 << "," << reads2 << endl;
            return (cube());
        }
        file2.close();
    }
    return (res_cube);
};

cube cube::operator-(cube &right) const
{
    cube res_cube(path, true, *parent_wavefunction, cout, false);
    res_cube.path = get_foldername_from_path(path) + get_filename_from_path(path).substr(0, get_filename_from_path(path).rfind(".cub")) + "-" + get_filename_from_path(right.path).substr(0, get_filename_from_path(right.path).rfind(".cub")) + ".cube";
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (cube());
    if (right.get_loaded())
#pragma omp parallel for
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    res_cube.set_value(x, y, z, res_cube.get_value(x, y, z) - right.get_value(x, y, z));
    else
    {
        int reads2 = 0;
        int reads1 = 0;
        int rest2 = 0;
        int run_x = 0;
        int run_y = 0;
        string line2;
        ifstream file2(right.path.c_str(), ios::in);
        double tmp[6]{0, 0, 0, 0, 0, 0};
        for (int i = 0; i < 6 + na; i++)
            if (!right.get_loaded())
                getline(file2, line2);
        while (run_x < size[0] && run_y < size[1] && !file2.eof())
        {
            reads2 = rest2;
            while (reads2 < size[2] && !file2.eof())
            {
                if (!right.get_loaded())
                {
                    getline(file2, line2);
                    reads1 = sscanf(line2.c_str(), "%lf%lf%lf%lf%lf%lf",
                                    &tmp[0],
                                    &tmp[1],
                                    &tmp[2],
                                    &tmp[3],
                                    &tmp[4],
                                    &tmp[5]);
                    reads2 += reads1;
                    rest2 = reads2 - size[2];
                    if (rest2 < 6 && rest2 > 0)
                    {
                        for (int j = 0; j < 6 - rest2; j++)
                            res_cube.set_value(run_x, run_y, size[2] - (6 - rest2) + j, res_cube.get_value(run_x, run_y, size[2] - (6 - rest2) + j) - tmp[j]);
                        if (run_y + 1 < size[1])
                            for (int j = 0; j < rest2; j++)
                                res_cube.set_value(run_x, run_y + 1, j, res_cube.get_value(run_x, run_y + 1, j) - tmp[(6 - rest2) + j]);
                        else if (run_x + 1 < size[0])
                            for (int j = 0; j < rest2; j++)
                                res_cube.set_value(run_x + 1, 0, j, res_cube.get_value(run_x + 1, 0, j) - tmp[(6 - rest2) + j]);
                        else
                        {
                            cout << "This should not happen! Read a value outside of range!"
                                 << " Run_x: " << run_x
                                 << " Run_y: " << run_y
                                 << " reads1: " << reads1
                                 << " reads2: " << reads2 << endl;
                            return cube();
                        }
                    }
                    else
                        for (int i = 0; i < reads1; i++)
                            res_cube.set_value(run_x, run_y, reads2 - reads1 + i, res_cube.get_value(run_x, run_y, reads2 - reads1 + i) - tmp[i]);
                }
                else
#pragma omp parallel for
                    for (int z = 0; z < size[2]; z++)
                        res_cube.set_value(run_x, run_y, z, res_cube.get_value(run_x, run_y, z) - right.get_value(run_x, run_y, z));
            }
            run_y++;
            if (run_y == size[1])
            {
                run_x++;
                if (run_x != size[0])
                    run_y = 0;
            }
        }
        if (run_x != size[0] || run_y != size[1])
        {
            cout << "This file ended before i read all expected values!" << endl;
            if (file2.eof())
                cout << "ENCOUNTERED EOF!" << endl;
            cout << "x,y,reads1,reads2: " << run_x << " " << run_y << " " << reads1 << "," << reads2 << endl;
            return (cube());
        }
        file2.close();
    }
    return (res_cube);
};

cube cube::operator*(cube &right) const
{
    cube res_cube(path, true, *parent_wavefunction, cout, false);
    res_cube.path = get_foldername_from_path(path) + get_filename_from_path(path).substr(0, get_filename_from_path(path).rfind(".cub")) + "*" + get_filename_from_path(right.path).substr(0, get_filename_from_path(right.path).rfind(".cub")) + ".cube";
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (cube());
    if (right.get_loaded())
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    res_cube.set_value(x, y, z, res_cube.get_value(x, y, z) + right.get_value(x, y, z));
    else
    {
        int reads2 = 0;
        int reads1 = 0;
        int rest2 = 0;
        int run_x = 0;
        int run_y = 0;
        string line2;
        ifstream file2(right.path.c_str(), ios::in);
        double tmp[6]{0, 0, 0, 0, 0, 0};
        for (int i = 0; i < 6 + na; i++)
            if (!right.get_loaded())
                getline(file2, line2);
        while (run_x < size[0] && run_y < size[1] && !file2.eof())
        {
            reads2 = rest2;
            while (reads2 < size[2] && !file2.eof())
            {
                if (!right.get_loaded())
                {
                    getline(file2, line2);
                    reads1 = sscanf(line2.c_str(), "%lf%lf%lf%lf%lf%lf",
                                    &tmp[0],
                                    &tmp[1],
                                    &tmp[2],
                                    &tmp[3],
                                    &tmp[4],
                                    &tmp[5]);
                    reads2 += reads1;
                    rest2 = reads2 - size[2];
                    if (rest2 < 6 && rest2 > 0)
                    {
                        for (int j = 0; j < 6 - rest2; j++)
                            res_cube.set_value(run_x, run_y, size[2] - (6 - rest2) + j, res_cube.get_value(run_x, run_y, size[2] - (6 - rest2) + j) * tmp[j]);
                        if (run_y + 1 < size[1])
                            for (int j = 0; j < rest2; j++)
                                res_cube.set_value(run_x, run_y + 1, j, res_cube.get_value(run_x, run_y + 1, j) * tmp[(6 - rest2) + j]);
                        else if (run_x + 1 < size[0])
                            for (int j = 0; j < rest2; j++)
                                res_cube.set_value(run_x + 1, 0, j, res_cube.get_value(run_x + 1, 0, j) * tmp[(6 - rest2) + j]);
                        else
                        {
                            cout << "This should not happen! Read a value outside of range!"
                                 << " Run_x: " << run_x
                                 << " Run_y: " << run_y
                                 << " reads1: " << reads1
                                 << " reads2: " << reads2 << endl;
                            return (cube());
                        }
                    }
                    else
                        for (int i = 0; i < reads1; i++)
                            res_cube.set_value(run_x, run_y, reads2 - reads1 + i, res_cube.get_value(run_x, run_y, reads2 - reads1 + i) * tmp[i]);
                }
                else
#pragma omp parallel for
                    for (int z = 0; z < size[2]; z++)
                        res_cube.set_value(run_x, run_y, z, res_cube.get_value(run_x, run_y, z) * right.get_value(run_x, run_y, z));
            }
            run_y++;
            if (run_y == size[1])
            {
                run_x++;
                if (run_x != size[0])
                    run_y = 0;
            }
        }
        if (run_x != size[0] || run_y != size[1])
        {
            cout << "This file ended before i read all expected values!" << endl;
            if (file2.eof())
                cout << "ENCOUNTERED EOF!" << endl;
            cout << "x,y,reads1,reads2: " << run_x << " " << run_y << " " << reads1 << "," << reads2 << endl;
            return (cube());
        }
        file2.close();
    }
    return (res_cube);
};

cube cube::operator/(cube &right) const
{
    cube res_cube(path, true, *parent_wavefunction, cout, false);
    res_cube.path = get_foldername_from_path(path) + get_filename_from_path(path).substr(0, get_filename_from_path(path).rfind(".cub")) + "_" + get_filename_from_path(right.path).substr(0, get_filename_from_path(right.path).rfind(".cub")) + ".cube";
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (cube());
    if (right.get_loaded())
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    res_cube.set_value(x, y, z, res_cube.get_value(x, y, z) + right.get_value(x, y, z));
    else
    {
        int reads2 = 0;
        int reads1 = 0;
        int rest2 = 0;
        int run_x = 0;
        int run_y = 0;
        string line2;
        ifstream file2(right.path.c_str(), ios::in);
        double tmp[6]{0, 0, 0, 0, 0, 0};
        for (int i = 0; i < 6 + na; i++)
            if (!right.get_loaded())
                getline(file2, line2);
        while (run_x < size[0] && run_y < size[1] && !file2.eof())
        {
            reads2 = rest2;
            while (reads2 < size[2] && !file2.eof())
            {
                if (!right.get_loaded())
                {
                    getline(file2, line2);
                    reads1 = sscanf(line2.c_str(), "%lf%lf%lf%lf%lf%lf",
                                    &tmp[0],
                                    &tmp[1],
                                    &tmp[2],
                                    &tmp[3],
                                    &tmp[4],
                                    &tmp[5]);
                    reads2 += reads1;
                    rest2 = reads2 - size[2];
                    if (rest2 < 6 && rest2 > 0)
                    {
                        for (int j = 0; j < 6 - rest2; j++)
                            res_cube.set_value(run_x, run_y, size[2] - (6 - rest2) + j, res_cube.get_value(run_x, run_y, size[2] - (6 - rest2) + j) / tmp[j]);
                        if (run_y + 1 < size[1])
                            for (int j = 0; j < rest2; j++)
                                res_cube.set_value(run_x, run_y + 1, j, res_cube.get_value(run_x, run_y + 1, j) / tmp[(6 - rest2) + j]);
                        else if (run_x + 1 < size[0])
                            for (int j = 0; j < rest2; j++)
                                res_cube.set_value(run_x + 1, 0, j, res_cube.get_value(run_x + 1, 0, j) / tmp[(6 - rest2) + j]);
                        else
                        {
                            cout << "This should not happen! Read a value outside of range!"
                                 << " Run_x: " << run_x
                                 << " Run_y: " << run_y
                                 << " reads1: " << reads1
                                 << " reads2: " << reads2 << endl;
                            return cube();
                        }
                    }
                    else
                        for (int i = 0; i < reads1; i++)
                            res_cube.set_value(run_x, run_y, reads2 - reads1 + i, res_cube.get_value(run_x, run_y, reads2 - reads1 + i) / tmp[i]);
                }
                else
#pragma omp parallel for
                    for (int z = 0; z < size[2]; z++)
                        res_cube.set_value(run_x, run_y, z, res_cube.get_value(run_x, run_y, z) / right.get_value(run_x, run_y, z));
            }
            run_y++;
            if (run_y == size[1])
            {
                run_x++;
                if (run_x != size[0])
                    run_y = 0;
            }
        }
        if (run_x != size[0] || run_y != size[1])
        {
            cout << "This file ended before i read all expected values!" << endl;
            if (file2.eof())
                cout << "ENCOUNTERED EOF!" << endl;
            cout << "x,y,reads1,reads2: " << run_x << " " << run_y << " " << reads1 << "," << reads2 << endl;
            return (cube());
        }
        file2.close();
    }
    return (res_cube);
};

bool cube::operator+=(cube &right)
{
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (false);
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                values[x][y][z] += right.get_value(x, y, z);
    return (true);
};

bool cube::operator-=(cube &right)
{
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (false);
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                values[x][y][z] -= right.get_value(x, y, z);
    return (true);
};

bool cube::operator*=(cube &right)
{
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (false);
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                values[x][y][z] *= right.get_value(x, y, z);
    return (true);
};

bool cube::operator/=(cube &right)
{
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (false);
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                values[x][y][z] /= right.get_value(x, y, z);
    return (true);
};

bool cube::mask(cube &right)
{
    if (size[0] != right.get_size(0) || size[1] != right.get_size(1) || size[2] != right.get_size(2))
        return (false);
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                if (right.get_value(x, y, z) == 0.0)
                    values[x][y][z] = 0.0;
    return (true);
};

bool cube::thresh(cube &right, double thresh)
{
    if (size[0] != right.get_size(0) || size[1] != right.get_size(1) || size[2] != right.get_size(2))
        return (false);
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                if (right.get_value(x, y, z) < thresh)
                    values[x][y][z] = 0.0;
    return (true);
};

bool cube::negative_mask(cube &right)
{
    if (size[0] != right.get_size(0) || size[1] != right.get_size(1) || size[2] != right.get_size(2))
        return (false);
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                if (right.get_value(x, y, z) != 0.0)
                    values[x][y][z] = 0.0;
    return (true);
};

double cube::rrs(cube &right)
{
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (-1);
    double diff_pos = 0.0;
    double diff_neg = 0.0;
#pragma omp parallel for reduction(+ : diff_pos, diff_neg)
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
            {
                diff_pos += abs((values[x][y][z] + right.get_value(x, y, z)));
                diff_neg += abs((values[x][y][z] - right.get_value(x, y, z)));
            }
    return (diff_neg / diff_pos); // RETURN Real Space R-value between this cube and the given one
};

double cube::sum()
{
    for (int i = 0; i < 3; i++)
        if (size[i] == 0)
            return (-1);
    double s = 0.0;
#pragma omp parallel for reduction(+ : s)
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                s += values[x][y][z];

    double dv = abs(vectors[0][0] * vectors[1][1] * vectors[2][2] - vectors[2][0] * vectors[1][1] * vectors[0][2] + vectors[0][1] * vectors[1][2] * vectors[2][0] - vectors[2][1] * vectors[1][2] * vectors[0][0] + vectors[0][2] * vectors[1][0] * vectors[2][1] - vectors[2][2] * vectors[1][0] * vectors[0][1]);
    return (s * dv); // RETURN Sum of values inside the cube
};

double cube::diff_sum()
{
    for (int i = 0; i < 3; i++)
        if (size[i] == 0)
            return (-1);
    double s = 0.0;
#pragma omp parallel for reduction(+ : s)
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                s += abs(values[x][y][z]) / 2;

    double dv = abs(vectors[0][0] * vectors[1][1] * vectors[2][2] - vectors[2][0] * vectors[1][1] * vectors[0][2] + vectors[0][1] * vectors[1][2] * vectors[2][0] - vectors[2][1] * vectors[1][2] * vectors[0][0] + vectors[0][2] * vectors[1][0] * vectors[2][1] - vectors[2][2] * vectors[1][0] * vectors[0][1]);
    return (s * dv); // RETURN Sum of values inside the cube
};

vec cube::double_sum()
{
    for (int i = 0; i < 3; i++)
    {
        if (size[i] == 0)
            return (vec(1, 0));
    }
    double s = 0.0;
    double s2 = 0.0;
#pragma omp parallel for reduction(+ : s, s2)
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
            {
                s += abs(values[x][y][z]) / 2;
                s2 += values[x][y][z];
            }

    double dv = abs(vectors[0][0] * vectors[1][1] * vectors[2][2] - vectors[2][0] * vectors[1][1] * vectors[0][2] + vectors[0][1] * vectors[1][2] * vectors[2][0] - vectors[2][1] * vectors[1][2] * vectors[0][0] + vectors[0][2] * vectors[1][0] * vectors[2][1] - vectors[2][2] * vectors[1][0] * vectors[0][1]);
    vec result;
    result.resize(2);
    result[0] = s * dv;
    result[1] = s2 * dv;
    return (result); // RETURN Sums of values inside the cube
};

double cube::get_vector(int i, int j) const
{
    if (i < 3 && i >= 0 && j < 3 && j >= 0)
        return (vectors[i][j]);
    else
        return (-1);
};

bool cube::set_vector(int i, int j, double value)
{
    if (i < 3 && i >= 0 && j < 3 && j >= 0)
        vectors[i][j] = value;
    else
        return (false);
    return (true);
};

double cube::get_origin(unsigned int i) const
{
    switch (i)
    {
    case 0:
        return (origin[0]);
        break;
    case 1:
        return (origin[1]);
        break;
    case 2:
        return (origin[2]);
        break;
    default:
        return (-1);
    }
};

bool cube::set_origin(unsigned int i, double value)
{
    if (i < 3 && i >= 0)
        origin[i] = value;
    else
        return (false);
    return (true);
};

string cube::super_cube()
{
    int m[3]{0, 0, 0};
    cout << "How many times in X-direction? ";
    cin >> m[0];
    while (m[0] <= 0 || m[0] > 20)
    {
        cout << "This is unreasonable, try again! (between 1-20): ";
        cin >> m[0];
    }
    cout << "Y-direction? ";
    cin >> m[1];
    while (m[1] <= 0 || m[1] > 20)
    {
        cout << "This is unreasonable, try again! (between 1-20): ";
        cin >> m[1];
    }
    cout << "Z-direction? ";
    cin >> m[2];
    while (m[2] <= 0 || m[2] > 20)
    {
        cout << "This is unreasonable, try again! (between 1-20): ";
        cin >> m[2];
    }
    string new_path(path);
    new_path += "_super";
    ofstream of(new_path.c_str(), ios::out);
    of << comment1 << " SUPER CUBE" << endl;
    of << comment2 << endl;
    of << "   " << m[0] * m[1] * m[2] * parent_wavefunction->get_ncen() << " " << origin[0] << " " << origin[1] << " " << origin[2] << endl;
    stringstream stream;
    string dummy;
    for (int i = 0; i < 3; i++)
    {
        stream << setw(6) << m[i] * size[i];
        for (int j = 0; j < 3; j++)
            stream << fixed << setw(12) << setprecision(6) << vectors[i][j];
        stream << endl;
    }
    of << stream.str();
    stream.str("");
    for (int h = 0; h < m[0]; h++)
        for (int j = 0; j < m[1]; j++)
            for (int k = 0; k < m[2]; k++)
                for (int a = 0; a < parent_wavefunction->get_ncen(); a++)
                {
                    for (int i = 0; i < 2; i++)
                        stream << setw(5) << parent_wavefunction->get_atom_charge(a);
                    stream << ".000000";
                    for (int i = 0; i < 3; i++)
                        stream << setw(12) << setprecision(6) << fixed << parent_wavefunction->get_atom_coordinate(a, i) + h * size[0] * vectors[0][i] + j * size[1] * vectors[1][i] + k * size[2] * vectors[2][i];
                    stream << endl;
                }
    of << stream.str();
    stream.str("");
    int runs = 0;
    for (int mult_x = 0; mult_x < m[0]; mult_x++)
        for (int run_x = 0; run_x < size[0]; run_x++)
            for (int mult_y = 0; mult_y < m[1]; mult_y++)
                for (int run_y = 0; run_y < size[1]; run_y++)
                {
                    runs = 0;
                    for (int mult_z = 0; mult_z < m[2]; mult_z++)
                        for (int run_z = 0; run_z < size[2]; run_z++)
                        {
                            stream << uppercase << setw(13) << setprecision(5) << scientific << values[run_x][run_y][run_z];
                            runs++;
                            if (runs % 6 == 0)
                                stream << endl;
                        }
                    if (runs % 6 != 0)
                        stream << endl;
                }
    of << stream.str();
    stream.str("");
    of.flush();
    of.close();
    return (new_path);
};

void cube::set_zero()
{
#pragma omp parallel for
    for (int i = 0; i < size[0]; i++)
        for (int j = 0; j < size[1]; j++)
            fill(values[i][j].begin(), values[i][j].end(), 0.0);
};
