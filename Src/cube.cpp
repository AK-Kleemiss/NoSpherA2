#include "pch.h"
#include "cube.h"
#include "convenience.h"
#include "constants.h"
#include "wfn_class.h"

cube::cube()
{
    loaded = false;
    na = 0;
    parent_wavefunction = new WFN(e_origin::cub);
    size = { 0, 0, 0 };
    origin = { 0.0, 0.0, 0.0 };
    dv = abs(vectors[0][0] * vectors[1][1] * vectors[2][2] - vectors[2][0] * vectors[1][1] * vectors[0][2] + vectors[0][1] * vectors[1][2] * vectors[2][0] - vectors[2][1] * vectors[1][2] * vectors[0][0] + vectors[0][2] * vectors[1][0] * vectors[2][1] - vectors[2][2] * vectors[1][0] * vectors[0][1]);
};

cube::cube(int x, int y, int z, int g_na, bool grow_values)
{
    size[0] = x;
    size[1] = y;
    size[2] = z;
    origin = { 0.0, 0.0, 0.0 };
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
    parent_wavefunction = new WFN(e_origin::cub);
    vectors = { { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } } };
    dv = abs(vectors[0][0] * vectors[1][1] * vectors[2][2] - vectors[2][0] * vectors[1][1] * vectors[0][2] + vectors[0][1] * vectors[1][2] * vectors[2][0] - vectors[2][1] * vectors[1][2] * vectors[0][0] + vectors[0][2] * vectors[1][0] * vectors[2][1] - vectors[2][2] * vectors[1][0] * vectors[0][1]);
};

cube::cube(const std::filesystem::path &filepath, bool read, WFN &wave, std::ostream &file, const bool expert, const bool header)
{
    err_checkf(exists(filepath), "Sorry, this file does not exist!", file);
    parent_wavefunction = &wave;
    na = parent_wavefunction->get_ncen();
    path = filepath;
    err_checkf(read_file(read, header, expert), "Sorry, something went wrong while reading!", file);
    loaded = read;
    dv = abs(vectors[0][0] * vectors[1][1] * vectors[2][2] - vectors[2][0] * vectors[1][1] * vectors[0][2] + vectors[0][1] * vectors[1][2] * vectors[2][0] - vectors[2][1] * vectors[1][2] * vectors[0][0] + vectors[0][2] * vectors[1][0] * vectors[2][1] - vectors[2][2] * vectors[1][0] * vectors[0][1]);
};

cube::cube(const int g_na, const ivec &g_size, const vec &g_origin, const vec2 &g_vectors, const vec3 &g_values)
{
    using namespace std;
    na = g_na;
    parent_wavefunction = new WFN(e_origin::cub);
   std::cout << "Assigned Nr of Atoms" << endl;
    for (int i = 0; i < 3; i++)
    {
       std::cout << i << ". dimension" << endl;
        size[i] = g_size[i];
        origin[i] = g_origin[i];
        for (int j = 0; j < 3; j++)
            vectors[i][j] = g_vectors[i][j];
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
    dv = abs(vectors[0][0] * vectors[1][1] * vectors[2][2] - vectors[2][0] * vectors[1][1] * vectors[0][2] + vectors[0][1] * vectors[1][2] * vectors[2][0] - vectors[2][1] * vectors[1][2] * vectors[0][0] + vectors[0][2] * vectors[1][0] * vectors[2][1] - vectors[2][2] * vectors[1][0] * vectors[0][1]);
};

cube::cube(const cube &given)
{
    na = given.get_na();
    path = given.path;
    comment1 = given.get_comment1();
    comment2 = given.get_comment2();
    dv = given.get_dv();
    for (int i = 0; i < 3; i++)
    {
        size[i] = given.get_size(i);
        origin[i] = given.get_origin(i);
        for (int j = 0; j < 3; j++)
            vectors[i][j] = given.get_vector(i, j);
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

bool cube::read_values(std::ifstream& file) {
    using namespace std;
    string line;
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
    int reads1 = 0;
    int rest2 = 0;
    int run_x = 0;
    int run_y = 0;
    int run_z = 0;
    double tmp[6] = { 0, 0, 0, 0, 0, 0 };
    while (run_x < size[0] && run_y < size[1] && !file.eof())
    {
        run_z = rest2;
        while (run_z < size[2] && !file.eof())
        {
            reads1 = 0;
            getline(file, line);
            std::istringstream iss(line);
            for (int i = 0; i < 6; ++i) {
                if (!(iss >> tmp[i]))
                    tmp[i] = std::nan(""); // Default value if there are fewer than 6 values
                else
                    reads1 = i + 1;
            }
            run_z += reads1;
            rest2 = run_z - size[2];
            if (rest2 < 6 && rest2 > 0)
            {
                for (int j = 0; j < 6 - rest2; j++) {
                    err_checkf(std::isnan(tmp[j]), "This should not happen! Read a value outside of range!", std::cout);
                    values[run_x][run_y][size[2] - (6 - rest2) + j] = tmp[j];
                }
                if (run_y + 1 < size[1])
                    for (int j = 0; j < rest2; j++) {
                        err_checkf(std::isnan(tmp[j + (6 - rest2)]), "This should not happen! Read a value outside of range!", std::cout);
                        values[run_x][run_y + 1][j] = tmp[j + (6 - rest2)];
                    }
                else if (run_x + 1 < size[0])
                    for (int j = 0; j < rest2; j++) {
                        err_checkf(std::isnan(tmp[j + (6 - rest2)]), "This should not happen! Read a value outside of range!", std::cout);
                        values[run_x + 1][0][j] = tmp[j + (6 - rest2)];
                    }
                else
                {
                   std::cout << "This should not happen! Read a value outside of range!"
                        << " Run_x: " << run_x
                        << " Run_y: " << run_y
                        << " rest2: " << rest2
                        << " run_z: " << run_z << endl;
                    return (false);
                }
            }
            else
                for (int i = 0; i < reads1; i++)
                    values[run_x][run_y][run_z - reads1 + i] = tmp[i];
        }
        run_y++;
        if (run_y == size[1])
        {
            run_x++;
            if (run_x != size[0])
                run_y = 0;
        }
    }
    if (run_x != size[0] || run_y != size[1] || run_z != size[2])
    {
       std::cout << "This file ended before i read all expected values!" << endl;
        if (file.eof())
           std::cout << "ENCOUNTERED EOF!" << endl;
       std::cout << "x,y,reads1,z: " << run_x << " " << run_y << " " << reads1 << "," << run_z << endl;
        return (false);
    }
    file.close();
    loaded = true;
    return true;
}

bool cube::read_file(bool full, bool header, bool expert)
{
    using namespace std;
    ifstream file(path);
    string line;
    if (header)
    {
        getline(file, comment1);
        getline(file, comment2);
        getline(file, line);
        std::istringstream iss(line);
        iss >> na >> origin[0] >> origin[1] >> origin[2];

        for (int i = 0; i < 3; i++) {
            getline(file, line);
            iss = std::istringstream(line);
            iss >> size[i] >> vectors[i][0] >> vectors[i][1] >> vectors[i][2];
        }

        int atnr;
        double atp[3] = { 0, 0, 0 };
        bool read_atoms = true;
        if (!expert && parent_wavefunction->get_ncen() != 0)
            read_atoms = false;
        if (read_atoms)
            for (int i = 0; i < na; i++)
            {
                getline(file, line);
                iss = std::istringstream(line);
                double dum;
                iss >> atnr >> dum >> atp[0] >> atp[1] >> atp[2];
                parent_wavefunction->push_back_atom(constants::atnr2letter(atnr), atp[0], atp[1], atp[2], atnr);
            }
    }
    if (full)
    {
        return read_values(file);
    }
    return (true);
}

bool cube::write_file(bool force, bool absolute)
{
    if (exists(path))
    {
        if (force)
            std::filesystem::remove(path);
        else
        {
            std::cout << "File already exists, aborting!" << std::endl;
            return false;
        }
    }
    using namespace std;
    std::stringstream stream;
    std::string temp;
    std::ofstream of(path, std::ios::out);
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

bool cube::write_file(const std::filesystem::path &given_path, bool debug)
{
    using namespace std;
    stringstream stream;
    string temp;
    ofstream of(given_path, ios::out);
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
       std::cout << "Finished atoms!" << endl;
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
                   std::cout << "Write Z-line!" << endl;
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

bool cube::write_xdgraph(const std::filesystem::path &given_path, bool debug)
{
    using namespace std;
    stringstream stream;
    ofstream of(given_path, ios::out);
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
        of << parent_wavefunction->get_atom_label(i);
        of << "    ";
        for (int j = 0; j < 3; j++)
            of << setw(10) << setprecision(5) << parent_wavefunction->get_atom_coordinate(i, j);
        of << " ATOM" << endl;
    }
    if (debug)
       std::cout << "Finished atoms!" << endl;
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
               std::cout << "Write Y-line!" << endl;
            if (temp_write % 6 != 0)
                of << endl;
        }
    of.flush();
    of.close();
    path = given_path;
    return (true);
};

bool cube::fractal_dimension(const double stepsize) const
{
    double min = 100, max = -100;
    for (const auto& inner_vec : values) {
      for (const auto& innerest_vec : inner_vec) {
        auto local_min_it = std::min_element(innerest_vec.begin(), innerest_vec.end());
        auto local_max_it = std::max_element(innerest_vec.begin(), innerest_vec.end());

        if (*local_min_it < min) {
          min = *local_min_it;
        }

        if (*local_max_it > max) {
          max = *local_max_it;
        }
      }
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
#pragma omp parallel private(lv1,lv2)
    {
      long long int x, y, z;
#pragma omp for
      for (long long int i = 0; i < (size_t)size[0] * size[1] * ((size_t)size[2] - 1); i++)
      {
        x = i / ((size_t)size[1] * ((size_t)size[2] - 1));
        y = (i / ((size_t)size[2] - 1)) % size[1];
        z = i % ((size_t)size[2] - 1);
        lv1 = values[x][y][z];
        lv2 = values[x][y][z + 1];
        for (int _i = 0; _i < steps; _i++)
          if ((lv1 - iso[_i]) * (lv2 - iso[_i]) < 0)
#pragma omp atomic
            bins[_i]++;
      }
#pragma omp for
      for (long long int i = 0; i < (size_t)size[0] * size[2] * ((size_t)size[1] - 1); i++)
      {
        z = i / ((size_t)size[0] * (size[1] - 1));
        x = (i / ((size_t)size[1] - 1)) % size[0];
        y = i % ((size_t)size[1] - 1);
        lv1 = values[x][y][z];
        lv2 = values[x][y + 1][z];
        for (int _i = 0; _i < steps; _i++)
          if ((lv1 - iso[_i]) * (lv2 - iso[_i]) < 0)
#pragma omp atomic
            bins[_i]++;
      }
#pragma omp for
      for (long long int i = 0; i < (size_t)size[1] * size[2] * ((size_t)size[0] - 1); i++)
      {
        y = i / (((size_t)size[0] - 1) * size[2]);
        z = (i / ((size_t)size[0] - 1)) % size[2];
        x = i % ((size_t)size[0] - 1);
        lv1 = values[x][y][z];
        lv2 = values[x + 1][y][z];
        for (int _i = 0; _i < steps; _i++)
          if ((lv1 - iso[_i]) * (lv2 - iso[_i]) < 0)
#pragma omp atomic
            bins[_i]++;
      }
    }
    const double epsilon = log(1 / (pow(comparisons, -constants::c_13)));
    for (int i = 0; i < steps; i++)
    {
        if (bins[i] == 0)
            df[i] = 0.0;
        else
            df[i] = log(bins[i]) / epsilon;
    }
    {
        using namespace std;
        std::filesystem::path out_file = path;
        out_file.replace_extension(".cube_fractal_plot");
        ofstream of(out_file, ios::out);
        of << setw(8) << scientific << setprecision(8) << steps
            << setw(16) << scientific << setprecision(8) << map_min
            << setw(16) << scientific << setprecision(8) << map_max
            << setw(16) << scientific << setprecision(8) << e[0] * pow(0.529177249, 3)
            << setw(16) << scientific << setprecision(8) << e[1] * pow(0.529177249, 3) << "\n";
        for (int i = 0; i < steps; i++)
            of << setw(16) << scientific << setprecision(8) << iso[i] << setw(16) << scientific << setprecision(8) << df[i] << "\n";
        of.flush();
        of.close();
    }
    return true;
}

int cube::get_size(int direction) const
{
    if (direction < size.size() && direction >= 0)
        return (size[direction]);
    else
        return (-1);
};

double cube::get_interpolated_value(double x, double y, double z) const
{
    if (x < origin[0] || y < origin[1] || z < origin[2] || x > origin[0] + vectors[0][0] * size[0] || y > origin[1] + vectors[1][1] * size[1] || z > origin[2] + vectors[2][2] * size[2])
        return (0.0);
    double x1 = (x - origin[0]) / vectors[0][0];
    double y1 = (y - origin[1]) / vectors[1][1];
    double z1 = (z - origin[2]) / vectors[2][2];
    int x0 = int(x1);
    int y0 = int(y1);
    int z0 = int(z1);
    x1 -= x0;
    y1 -= y0;
    z1 -= z0;
    double c00 = get_value(x0, y0, z0) * (1 - x1) + get_value(x0 + 1, y0, z0) * x1;
    double c10 = get_value(x0, y0 + 1, z0) * (1 - x1) + get_value(x0 + 1, y0 + 1, z0) * x1;
    double c01 = get_value(x0, y0, z0 + 1) * (1 - x1) + get_value(x0 + 1, y0, z0 + 1) * x1;
    double c11 = get_value(x0, y0 + 1, z0 + 1) * (1 - x1) + get_value(x0 + 1, y0 + 1, z0 + 1) * x1;
    double c0 = c00 * (1 - y1) + c10 * y1;
    double c1 = c01 * (1 - y1) + c11 * y1;
    return (c0 * (1 - z1) + c1 * z1);
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

void cube::set_dv(const double& given)
{
    dv = given;
};

void cube::calc_dv() {
    dv = abs(vectors[0][0] * vectors[1][1] * vectors[2][2] - vectors[2][0] * vectors[1][1] * vectors[0][2] + vectors[0][1] * vectors[1][2] * vectors[2][0] - vectors[2][1] * vectors[1][2] * vectors[0][0] + vectors[0][2] * vectors[1][0] * vectors[2][1] - vectors[2][2] * vectors[1][0] * vectors[0][1]);
};

// Function to compute dot product
template <typename T>
inline const double dot_(const T& a, const T& b) {
    double result = 0;
    for (int i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

bool has_converged(const double& current_value, double& previous_value, const double rel_threshold, const double rsE, std::deque<double>& history, const int window_size) {
    // Calculate relative and absolute differences
    double relative_diff = std::abs((current_value - previous_value) / current_value);
    double absolute_diff = std::abs(current_value - previous_value);

    // Update history
    history.push_back(current_value);
    if (history.size() > window_size) {
        history.pop_front();
    }
    double moving_average = 0;
    bool moving_average_converged = false;
    bool absolute_converged = false;

    // Calculate moving average
    if (history.size() != 1) {
        double sum = 0.0;
        for (double value : history) {
            sum += value;
        }
        moving_average = sum / history.size();
        moving_average_converged = std::abs((current_value - moving_average) / moving_average) < rel_threshold;
        absolute_converged = absolute_diff < rel_threshold * rsE;
    }

    // Check convergence criteria
    bool relative_converged = relative_diff < rel_threshold;

    previous_value = current_value;
    if (relative_converged || absolute_converged || moving_average_converged) {
        std::cout << "Converged rel/abs/moving " << relative_converged << "/" << absolute_converged << "/" << moving_average_converged << ": " << current_value << std::endl;
        return true;
    }
    else {
        std::cout << "Not Converged, due to rel/abs/moving " << relative_converged << "/" << absolute_converged << "/" << moving_average_converged << ": " << current_value << std::endl;
        return false;
    }
}

double cube::ewald_sum(const int kMax, const double conv) {
    calc_dv();
    const std::array<double, 3> lengths{ array_length(vectors[0])*size[0], array_length(vectors[1])*size[1], array_length(vectors[2])*size[2]};
    const double shortest_length = std::min({ lengths[0], lengths[1], lengths[2] });
    const double alpha = 2.0 * constants::sqr_pi / shortest_length;
    std::cout << "shortest length: " << shortest_length << " alpha: " << alpha << std::endl;
    double realSpaceEnergy = 0.0;
    double reciprocalSpaceEnergy = 0.0;
    double selfEnergy = 0.0;

    std::array<std::array<double, 3>, 3> cell_vectors;
    cell_vectors[0] = { vectors[0][0] * size[0], vectors[1][0] * size[0], vectors[2][0] * size[0] };
    cell_vectors[1] = { vectors[0][1] * size[1], vectors[1][1] * size[1], vectors[2][1] * size[1] };
    cell_vectors[2] = { vectors[0][2] * size[2], vectors[1][2] * size[2], vectors[2][2] * size[2] };
    std::cout << "Cell lattice:" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << std::setw(10) << cell_vectors[i][0] << std::setw(10) << cell_vectors[i][1] << std::setw(10) << cell_vectors[i][2] << std::endl;
    }

    // Compute volume of the unit cell
    const std::array<double, 3> crossProduct = cross(cell_vectors[1], cell_vectors[2]);
    const double volume = fabs(dot_(cell_vectors[0], crossProduct));
    std::cout << "Volume: " << volume << std::endl;
    const int grid_points = size[0] * size[1] * size[2];
    std::cout << "Number of grid points: " << grid_points << std::endl;
    std::cout << "dv*points: " << dv * grid_points << std::endl;

    // Compute reciprocal lattice vectors
    std::array<std::array<double, 3>, 3> reciprocalLattice = {
        cross(cell_vectors[1], cell_vectors[2]),
        cross(cell_vectors[2], cell_vectors[0]),
        cross(cell_vectors[0], cell_vectors[1])
    };
    for (auto& vec : reciprocalLattice) {
        for (double& x : vec) x *= 2 * constants::PI / volume;
    }
    std::cout << "Reciprocal cell lattice:" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << std::setw(10) << reciprocalLattice[i][0] << std::setw(10) << reciprocalLattice[i][1] << std::setw(10) << reciprocalLattice[i][2] << std::endl;
    }


    std::vector<std::array<double, 3>> ri(grid_points);
    std::vector<std::vector<std::array<double, 3>>> rij(grid_points);
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            for (int k = 0; k < size[2]; k++) {
                ri[i * size[1] * size[2] + j * size[2] + k] = get_pos(i, j, k);
                rij[i * size[1] * size[2] + j * size[2] + k].resize(grid_points);
                for (int l = 0; l < size[0]; l++) {
                    for (int m = 0; m < size[1]; m++) {
                        for (int n = 0; n < size[2]; n++) {
                            std::array<double, 3> rj = get_pos(l, m, n);
                            rij[i * size[1] * size[2] + j * size[2] + k][l * size[1] * size[2] + m * size[2] + n] =
                            { ri[i * size[1] * size[2] + j * size[2] + k][0] - rj[0],
                              ri[i * size[1] * size[2] + j * size[2] + k][1] - rj[1],
                              ri[i * size[1] * size[2] + j * size[2] + k][2] - rj[2]};
                        }
                    }
                }
            }
        }
    }

    // Real-space contribution
#pragma omp parallel for reduction(+:realSpaceEnergy)
    for (int i = 0; i < size[0]; i++) {
        double length = 0, fac = 0, v1 = 0;
        for (int j = 0; j < size[1]; j++) {
            for (int k = 0; k < size[2]; k++) {
                v1 = values[i][j][k];
#pragma ivdep
                for (int l = 0; l < size[0]; l++) {
                    for (int m = 0; m < size[1]; m++) {
                        for (int n = 0; n < size[2]; n++) {
                            length = array_length(rij[i * size[1] * size[2] + j * size[2] + k][l * size[1] * size[2] + m * size[2] + n]);
                            if (length > 6.0 || length == 0) continue;
                            fac = erfc(alpha * length) / length;
                            if (abs(fac) < 1E-10) continue;
                            realSpaceEnergy += v1 * values[l][m][n] * fac;
                        }
                    }
                }
            }
        }
    }
    realSpaceEnergy *= dv * dv * 0.5;
    std::cout << "Real-space energy: " << realSpaceEnergy << std::endl;
    double result = 0.0;
    double res_temp = 0;
    const double FOUR_alsq = 4 * alpha * alpha;
    std::vector<std::array<int, 3>> known_kVecs;
    const int window_size = 5;
    std::deque<double> history;
    for (int k_vec = 1; k_vec <= kMax; k_vec++) {
        double temp = 0;
        for (int h = -k_vec; h <= k_vec; ++h) {
            for (int k = -k_vec; k <= k_vec; ++k) {
                for (int l = -k_vec; l <= k_vec; ++l) {
                    if (h == 0 && k == 0 && l == 0) continue;
                    if (abs(h) + abs(k) + abs(l) != k_vec) continue;
                    bool known = false;
                    for (int i = 0; i < known_kVecs.size(); i++) {
                        if (((known_kVecs[i][0] == +h) && (known_kVecs[i][1] == +k) && (known_kVecs[i][2] == +l)) ||
                            ((known_kVecs[i][0] == -h) && (known_kVecs[i][1] == -k) && (known_kVecs[i][2] == -l))
                            ) {
                            known = true;
                            break;
                        }
                    }
                    if (known) continue;
                    known_kVecs.push_back({ h, k, l });

                    const std::array<double, 3> kvec = {
                        h * reciprocalLattice[0][0] + k * reciprocalLattice[1][0] + l * reciprocalLattice[2][0],
                        h * reciprocalLattice[0][1] + k * reciprocalLattice[1][1] + l * reciprocalLattice[2][1],
                        h * reciprocalLattice[0][2] + k * reciprocalLattice[1][2] + l * reciprocalLattice[2][2]
                    };
                    const double k2 = dot_(kvec, kvec);
#pragma omp parallel for reduction(+:temp)
                    for (int d1 = 0; d1 < size[0]; d1++) {
                        double v1 = 0, kDotR = 0;
#pragma ivdep
                        for (int d2 = 0; d2 < size[1]; d2++) {
                            for (int d3 = 0; d3 < size[2]; d3++) {
                                v1 = values[d1][d2][d3];
                                for (int d4 = 0; d4 < size[0]; d4++) {
                                    for (int d5 = 0; d5 < size[1]; d5++) {
                                        for (int d6 = 0; d6 < size[2]; d6++) {
                                            kDotR = dot_(kvec, rij[d1 * size[1] * size[2] + d2 * size[2] + d3][d4 * size[1] * size[2] + d5 * size[2] + d6]);
                                            temp += exp(-k2 / FOUR_alsq) / abs(k2) * v1 * values[d4][d5][d6] * cos(kDotR);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        res_temp += 4.0 * constants::PI / volume * temp * dv * dv;
        if (has_converged(res_temp, result, conv, realSpaceEnergy, history, window_size)) {
            break;
        }
    }

    reciprocalSpaceEnergy = result;

    // Self-energy term
#pragma omp parallel for reduction(+:selfEnergy)
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            for (int m = 0; m < size[2]; m++) {
                selfEnergy += values[i][j][m] * values[i][j][m];
            }
        }
    }
    //assuming the total charge is zero no need for charged system term
    selfEnergy *= -alpha / constants::sqr_pi * dv * dv;

    // Total energy
    const double totalEnergy = realSpaceEnergy + reciprocalSpaceEnergy + selfEnergy;
    std::cout << "Real-space energy: " << realSpaceEnergy << std::endl;
    std::cout << "Reciprocal-space energy: " << reciprocalSpaceEnergy << std::endl;
    std::cout << "Self-energy: " << selfEnergy << std::endl;
    std::cout << "Total energy: " << totalEnergy << std::endl;
    return totalEnergy;
}

void cube::operator=(const cube &right)
{
    for (int i = 0; i < 3; i++)
        size[i] = right.get_size(i);
    comment1 = right.get_comment1();
    comment2 = right.get_comment2();
    na = right.get_na();
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

cube cube::operator+(const cube &right) const
{
    cube res_cube(*this);
    res_cube.path = path.parent_path() / std::filesystem::path(path.stem().string() + "+" + right.get_path().stem().string() + ".cube");
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (cube());
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                res_cube.set_value(x, y, z, res_cube.get_value(x, y, z) + right.get_value(x, y, z));

    return (res_cube);
};

cube cube::operator-(const cube &right) const
{
    cube res_cube(*this);
    res_cube.path = path.parent_path() / std::filesystem::path(path.stem().string() + "-" + right.get_path().stem().string() + ".cube");
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (cube());
#pragma omp parallel for
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    res_cube.set_value(x, y, z, res_cube.get_value(x, y, z) - right.get_value(x, y, z));

    return (res_cube);
};

cube cube::operator*(const cube &right) const
{
    cube res_cube(*this);
    res_cube.path = path.parent_path() / std::filesystem::path(path.stem().string() + "*" + right.get_path().stem().string() + ".cube");
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (cube());
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                res_cube.set_value(x, y, z, res_cube.get_value(x, y, z) + right.get_value(x, y, z));

    return (res_cube);
};

cube cube::operator/(const cube &right) const
{
    cube res_cube(*this);
    res_cube.path = path.parent_path() / std::filesystem::path(path.stem().string() + "_" + right.get_path().stem().string() + ".cube");
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (cube());
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++) {
                if (right.get_value(x, y, z) == 0)
                    res_cube.set_value(x, y, z, 1E100);
                if (values[x][y][z] != 0)
                    res_cube.set_value(x, y, z, values[x][y][z] / right.get_value(x, y, z));
                else
                    res_cube.set_value(x, y, z, 0);
            }

    return (res_cube);
};

bool cube::operator+=(const cube &right)
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

bool cube::operator-=(const cube &right)
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

bool cube::operator*=(const cube &right)
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

bool cube::operator/=(const cube &right)
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

bool cube::mask(const cube &right)
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

bool cube::thresh(const double& thresh)
{
    if (size[0] == 0 || size[1] == 0 || size[2] == 0)
        return (false);
#pragma omp parallel for
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                if (values[x][y][z] < thresh)
                    values[x][y][z] = 0.0;
    return (true);
}

bool cube::thresh(const cube &right, const double& thresh)
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

bool cube::negative_mask(const cube &right)
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

double cube::rrs(const cube& right) const
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

double cube::sum() const
{
    for (int i = 0; i < 3; i++)
        if (size[i] == 0)
            return (-1);
    double _s = 0.0;
#pragma omp parallel for reduction(+ : _s)
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                _s += values[x][y][z];

    return (_s * dv); // RETURN Sum of values inside the cube
};

double cube::diff_sum() const
{
    for (int i = 0; i < 3; i++)
        if (size[i] == 0)
            return (-1);
    double _s = 0.0;
#pragma omp parallel for reduction(+ : _s)
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
                _s += abs(values[x][y][z]) / 2;

    return (_s * dv); // RETURN Sum of values inside the cube
};

vec cube::double_sum() const
{
    for (int i = 0; i < 3; i++)
    {
        if (size[i] == 0)
            return (vec(1, 0));
    }
    double _s = 0.0;
    double _s2 = 0.0;
#pragma omp parallel for reduction(+ : _s, _s2)
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++)
            {
                _s += abs(values[x][y][z]) / 2;
                _s2 += values[x][y][z];
            }

    vec result({ _s * dv, _s2 * dv });
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

std::filesystem::path cube::super_cube()
{
    using namespace std;
    int m[3]{0, 0, 0};
   std::cout << "How many times in X-direction? ";
    cin >> m[0];
    while (m[0] <= 0 || m[0] > 20)
    {
       std::cout << "This is unreasonable, try again! (between 1-20): ";
        cin >> m[0];
    }
   std::cout << "Y-direction? ";
    cin >> m[1];
    while (m[1] <= 0 || m[1] > 20)
    {
       std::cout << "This is unreasonable, try again! (between 1-20): ";
        cin >> m[1];
    }
   std::cout << "Z-direction? ";
    cin >> m[2];
    while (m[2] <= 0 || m[2] > 20)
    {
       std::cout << "This is unreasonable, try again! (between 1-20): ";
        cin >> m[2];
    }
    std::filesystem::path new_path(path);
    new_path.replace_extension(".cube_super");
    ofstream of(new_path, ios::out);
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

cube cube::super_cube(int x, int y, int z)
{
    using namespace std;
    int m[3]{ x, y, z };
    std::filesystem::path new_path(path);
    WFN super_wfn(parent_wavefunction->get_path());
    for (int h = 0; h < m[0]; h++)
        for (int j = 0; j < m[1]; j++)
            for (int k = 0; k < m[2]; k++)
                for (int a = 0; a < parent_wavefunction->get_ncen(); a++)
                {
                    super_wfn.push_back_atom(constants::atnr2letter(parent_wavefunction->get_atom_charge(a)), parent_wavefunction->get_atom_coordinate(a, 0) + h * size[0] * vectors[0][0] + j * size[1] * vectors[1][0] + k * size[2] * vectors[2][0], parent_wavefunction->get_atom_coordinate(a, 1) + h * size[0] * vectors[0][1] + j * size[1] * vectors[1][1] + k * size[2] * vectors[2][1], parent_wavefunction->get_atom_coordinate(a, 2) + h * size[0] * vectors[0][2] + j * size[1] * vectors[1][2] + k * size[2] * vectors[2][2], parent_wavefunction->get_atom_charge(a));
                }

    new_path.replace_extension(".cube_super");
    cube out = cube(m[0] * size[0], m[1] * size[1], m[2] * size[2], m[0] * m[1] * m[2] * parent_wavefunction->get_ncen(), true);
    out.set_path(new_path);
    out.parent_wavefunction = &super_wfn;
    out.set_comment1(comment1 + " SUPER CUBE");
    out.set_comment2(comment2);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            out.set_vector(i, j, vectors[i][j]);
        out.set_origin(i, origin[i]);
    }
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
                            out.set_value(run_x + mult_x * size[0], run_y + mult_y * size[1], run_z + mult_z * size[2], values[run_x][run_y][run_z]);
                        }
                }
    return (out);
};
std::vector<atom> cube::get_parent_wfn_atoms() const {
    return parent_wavefunction->get_atoms();
};

void cube::set_zero()
{
#pragma omp parallel for
    for (int i = 0; i < size[0]; i++)
        for (int j = 0; j < size[1]; j++)
            fill(values[i][j].begin(), values[i][j].end(), 0.0);
};

double cube::jaccard(const cube& right) const {
    for (int i = 0; i < 3; i++)
        if (size[i] != right.get_size(i))
            return (-1);
    double top = 0.0;
    double bot = 0.0;
    for (int x = 0; x < size[0]; x++)
        for (int y = 0; y < size[1]; y++)
            for (int z = 0; z < size[2]; z++) {
                top += std::min(values[x][y][z], right.get_value(x, y, z));
                bot += std::max(values[x][y][z], right.get_value(x, y, z));
            }
    return (top / bot); //RETURN Real Space R-value between this cube and the given one
};
