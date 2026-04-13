#pragma once

#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <type_traits>
#include "convenience.h"

class WFN;
class atom;
// Structure to hold charge and position
struct Charge {
    std::array<double, 3> r; // Fractional position (x, y, z)
    double q;                // Charge
};

class cube
{
public:
    cube();
    cube(const std::array<int, 3> xyz, int g_na = 0, bool grow_values = false);
    cube(const std::filesystem::path& filepath, bool read, WFN& wave, std::ostream& file, const bool expert = false, const bool header = true);
    cube(const int g_na, const std::vector <int>& g_size, const std::vector <double>& g_origin, const std::vector < std::vector<double> >& g_vectors, const std::vector<std::vector<std::vector<double> > >& g_values);
    cube(const cube& given);
    int get_size(int direction) const;
    i3 get_sizes() const { return size; };
    bool get_loaded() const { return loaded; };
    cube operator + (const cube& right) const;
    cube operator - (const cube& right) const;
    cube operator * (const cube& right) const;
    cube operator / (const cube& right) const;
    bool operator += (const cube& right);
    bool operator -= (const cube& right);
    bool operator *= (const cube& right);
    bool operator /= (const cube& right);
    void operator = (const cube& right);
    bool mask(const cube& right);
    bool thresh(const cube& right, const double& thresh = -1234);
    bool thresh(const double& thresh);
    bool negative_mask(const cube& right);
    double rrs(const cube& right) const;
    double sum() const;
    double diff_sum() const;
    std::vector<double> double_sum() const;
    double min_value() const {
        if (size[0] == 0 || size[1] == 0 || size[2] == 0)
            return 0.0;
        double result = values[0][0][0];
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    if (values[x][y][z] < result)
                        result = values[x][y][z];
        return result;
    }
    double max_value() const {
        if (size[0] == 0 || size[1] == 0 || size[2] == 0)
            return 0.0;
        double result = values[0][0][0];
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    if (values[x][y][z] > result)
                        result = values[x][y][z];
        return result;
    }
    template<typename T>
    std::array<double, 3> get_pos(const T& i, const T& j, const T& k) const {
        return {
            i * vectors[0][0] + j * vectors[0][1] + k * vectors[0][2] + origin[0],
            i * vectors[1][0] + j * vectors[1][1] + k * vectors[1][2] + origin[1],
            i * vectors[2][0] + j * vectors[2][1] + k * vectors[2][2] + origin[2]
        };
    };
    double get_interpolated_value(double x, double y, double z) const;
    double get_value(int x, int y, int z) const;
    bool set_value(int x, int y, int z, double value);
    bool read_file(bool full, bool header, bool expert = false);
    bool write_file(bool force = false, bool absolute = false);
    bool write_file(const std::filesystem::path& given_path, bool debug = false);
    bool write_xdgraph(const std::filesystem::path& given_path, bool debug = false);
    bool fractal_dimension(const double stepsize) const;
    double get_vector(int i, int j) const;
    std::array<d3, 3> get_vectors() const;
    bool set_vector(int i, int j, double value);
    void set_vectors(std::array<d3, 3> value);
    double get_origin(unsigned int i) const;
    double ewald_sum(const int kMax = 15, const double conv = 5E-3);
    void calc_dv();
    double get_dv() const { return dv; };
    void set_dv(const double& given);
    bool set_origin(unsigned int i, double value);
    int get_na() const { return na; };
    void set_na(int g_na) { na = g_na; };
    std::filesystem::path super_cube();
    cube super_cube(int x, int y, int z);
    void set_comment1(std::string input) { comment1 = input; };
    void set_comment2(std::string input) { comment2 = input; };
    void set_zero();
    bool evaluate_on_grid(const std::function<double(const d3&)>& func, bool wrap = false);
    bool evaluate_on_grid(const std::function<double(const d3&, const i3&, const i3&)>& func, bool wrap = false);
    void give_parent_wfn(WFN& given) { parent_wavefunction = &given; };
    std::string get_comment1() const { return comment1; };
    std::string get_comment2() const { return comment2; };
    std::filesystem::path get_path() const { return path; };
    void set_path(const std::filesystem::path& given) { path = given; };
    void resize(const i3& g_size) { size = g_size; values.resize(size[0]); for (int i = 0; i < size[0]; i++) { values[i].resize(size[1]); for (int j = 0; j < size[1]; j++) values[i][j].resize(size[2]); } }
    std::vector<atom> get_parent_wfn_atoms() const;
    bool read_values(std::ifstream& file);
    double jaccard(const cube& right) const;

    /**
     * @brief Adaptively refines the grid to reach a target integration error.
     *
     * This function iteratively increases the grid resolution (doubling it in each step)
     * until the integral of the values (calculated via sum()) converges within the
     * specified tolerance. It only evaluates the provided function at new grid points
     * where the local variation suggests it is necessary, otherwise interpolating
     * from the previous grid to save computation time (though for a regular grid structure,
     * we still store all points).
     *
     * @param func The function to evaluate at grid positions. Takes d3 and returns double.
     * @param target_error The relative error threshold for the integral convergence.
     * @param max_depth Maximum number of refinement iterations (limit total grid size).
     */
    void adaptive_refine(std::function<const double(const std::array<double, 3>)> const func, double target_error, int max_depth = 4, const int subfactor = 2);

private:
    double dv;
    int na;
    bool loaded;
    std::string comment1;
    std::string comment2;
    i3 size;
    d3 origin;
    std::array <d3, 3> vectors;
    vec3 values;
    std::filesystem::path path;
    WFN* parent_wavefunction;
    Refinepointmap refine_points;
};

template<typename T>
class cube_t {
public:
    using value_type = T;
    static_assert(std::is_same_v<T, int> || std::is_same_v<T, double>,
        "cube_t<T>: T must be int or double");

    cube_t() = delete;
    T get_t_value(int, int, int) const = delete;
    bool set_t_value(int, int, int, T) = delete;
};

template<>
class cube_t<double> : public cube {
public:
    using value_type = double;
    using cube::cube;

    double get_t_value(int x, int y, int z) const { return cube::get_value(x, y, z); }
    bool set_t_value(int x, int y, int z, double value) { return cube::set_value(x, y, z, value); }
};

template<>
class cube_t<int> {
public:
    using value_type = int;

    cube_t()
    {
        loaded = false;
        na = 0;
        parent_wavefunction = nullptr;
        size = { 0, 0, 0 };
        origin = { 0.0, 0.0, 0.0 };
        vectors = { { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } } };
        calc_dv();
    }

    cube_t(const std::array<int, 3> xyz, int g_na = 0, bool grow_values = false)
    {
        size = xyz;
        origin = { 0.0, 0.0, 0.0 };
        loaded = grow_values;
        if (grow_values)
            resize(size);
        na = g_na;
        parent_wavefunction = nullptr;
        vectors = { { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } } };
        calc_dv();
    }

    cube_t(const cube& given)
    {
        na = given.get_na();
        path = given.get_path();
        comment1 = given.get_comment1();
        comment2 = given.get_comment2();
        dv = given.get_dv();
        for (int i = 0; i < 3; i++) {
            size[i] = given.get_size(i);
            origin[i] = given.get_origin(i);
            for (int j = 0; j < 3; j++)
                vectors[i][j] = given.get_vector(i, j);
        }
        loaded = given.get_loaded();
        parent_wavefunction = nullptr;
        if (loaded) {
            resize(size);
            for (int x = 0; x < size[0]; x++)
                for (int y = 0; y < size[1]; y++)
                    for (int z = 0; z < size[2]; z++)
                        values[x][y][z] = static_cast<int>(std::llround(given.get_value(x, y, z)));
        }
    }

    int get_size(int direction) const {
        if (direction < static_cast<int>(size.size()) && direction >= 0)
            return size[direction];
        return -1;
    }

    i3 get_sizes() const { return size; }
    bool get_loaded() const { return loaded; }

    template<typename U>
    std::array<double, 3> get_pos(const U& i, const U& j, const U& k) const {
        return {
            i * vectors[0][0] + j * vectors[0][1] + k * vectors[0][2] + origin[0],
            i * vectors[1][0] + j * vectors[1][1] + k * vectors[1][2] + origin[1],
            i * vectors[2][0] + j * vectors[2][1] + k * vectors[2][2] + origin[2]
        };
    }

    int get_value(int x, int y, int z) const {
        if (x < size[0] && y < size[1] && z < size[2] && x >= 0 && y >= 0 && z >= 0)
            return values[x][y][z];
        return -1;
    }

    bool set_value(int x, int y, int z, int value) {
        if (x < size[0] && y < size[1] && z < size[2] && x >= 0 && y >= 0 && z >= 0) {
            values[x][y][z] = value;
            return true;
        }
        return false;
    }

    int get_t_value(int x, int y, int z) const { return get_value(x, y, z); }
    bool set_t_value(int x, int y, int z, int value) { return set_value(x, y, z, value); }

    cube_t<int> operator + (const cube_t<int>& right) const {
        cube_t<int> res(*this);
        for (int i = 0; i < 3; i++)
            if (size[i] != right.get_size(i))
                return cube_t<int>();
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    res.values[x][y][z] += right.get_value(x, y, z);
        return res;
    }

    cube_t<int> operator - (const cube_t<int>& right) const {
        cube_t<int> res(*this);
        for (int i = 0; i < 3; i++)
            if (size[i] != right.get_size(i))
                return cube_t<int>();
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    res.values[x][y][z] -= right.get_value(x, y, z);
        return res;
    }

    cube_t<int> operator * (const cube_t<int>& right) const {
        cube_t<int> res(*this);
        for (int i = 0; i < 3; i++)
            if (size[i] != right.get_size(i))
                return cube_t<int>();
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    res.values[x][y][z] *= right.get_value(x, y, z);
        return res;
    }

    cube_t<int> operator / (const cube_t<int>& right) const {
        cube_t<int> res(*this);
        for (int i = 0; i < 3; i++)
            if (size[i] != right.get_size(i))
                return cube_t<int>();
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++) {
                    const int rv = right.get_value(x, y, z);
                    res.values[x][y][z] = (rv == 0) ? 0 : (values[x][y][z] / rv);
                }
        return res;
    }

    bool operator += (const cube_t<int>& right) {
        for (int i = 0; i < 3; i++) if (size[i] != right.get_size(i)) return false;
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    values[x][y][z] += right.get_value(x, y, z);
        return true;
    }

    bool operator -= (const cube_t<int>& right) {
        for (int i = 0; i < 3; i++) if (size[i] != right.get_size(i)) return false;
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    values[x][y][z] -= right.get_value(x, y, z);
        return true;
    }

    bool operator *= (const cube_t<int>& right) {
        for (int i = 0; i < 3; i++) if (size[i] != right.get_size(i)) return false;
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    values[x][y][z] *= right.get_value(x, y, z);
        return true;
    }

    bool operator /= (const cube_t<int>& right) {
        for (int i = 0; i < 3; i++) if (size[i] != right.get_size(i)) return false;
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++) {
                    const int rv = right.get_value(x, y, z);
                    values[x][y][z] = (rv == 0) ? 0 : (values[x][y][z] / rv);
                }
        return true;
    }

    bool mask(const cube_t<int>& right) {
        if (size[0] != right.get_size(0) || size[1] != right.get_size(1) || size[2] != right.get_size(2))
            return false;
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    if (right.get_value(x, y, z) == 0)
                        values[x][y][z] = 0;
        return true;
    }

    bool thresh(const cube_t<int>& right, const double& thresh = -1234) {
        if (size[0] != right.get_size(0) || size[1] != right.get_size(1) || size[2] != right.get_size(2))
            return false;
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    if (right.get_value(x, y, z) < thresh)
                        values[x][y][z] = 0;
        return true;
    }

    bool negative_mask(const cube_t<int>& right) {
        if (size[0] != right.get_size(0) || size[1] != right.get_size(1) || size[2] != right.get_size(2))
            return false;
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    if (right.get_value(x, y, z) != 0)
                        values[x][y][z] = 0;
        return true;
    }

    double sum() const {
        double s = 0.0;
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    s += values[x][y][z];
        return s * dv;
    }

    int min_value() const {
        if (size[0] == 0 || size[1] == 0 || size[2] == 0)
            return 0;
        int result = values[0][0][0];
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    if (values[x][y][z] < result)
                        result = values[x][y][z];
        return result;
    }

    int max_value() const {
        if (size[0] == 0 || size[1] == 0 || size[2] == 0)
            return 0;
        int result = values[0][0][0];
        for (int x = 0; x < size[0]; x++)
            for (int y = 0; y < size[1]; y++)
                for (int z = 0; z < size[2]; z++)
                    if (values[x][y][z] > result)
                        result = values[x][y][z];
        return result;
    }

    double get_vector(int i, int j) const {
        if (i < 3 && i >= 0 && j < 3 && j >= 0)
            return vectors[i][j];
        return -1;
    }

    bool set_vector(int i, int j, double value) {
        if (i < 3 && i >= 0 && j < 3 && j >= 0) {
            vectors[i][j] = value;
            calc_dv();
            return true;
        }
        return false;
    }

    std::array<d3, 3> get_vectors() const { return vectors; }

    double get_origin(unsigned int i) const {
        if (i < 3)
            return origin[i];
        return 0.0;
    }

    bool set_origin(unsigned int i, double value) {
        if (i < 3) {
            origin[i] = value;
            return true;
        }
        return false;
    }

    int get_na() const { return na; }
    void set_na(int g_na) { na = g_na; }

    double get_dv() const { return dv; }
    void set_dv(const double& given) { dv = given; }

    std::filesystem::path get_path() const { return path; }
    void set_path(const std::filesystem::path& given) { path = given; }

    std::string get_comment1() const { return comment1; }
    std::string get_comment2() const { return comment2; }
    void set_comment1(std::string input) { comment1 = input; }
    void set_comment2(std::string input) { comment2 = input; }

    void resize(const i3& g_size) {
        size = g_size;
        values.resize(size[0]);
        for (int i = 0; i < size[0]; i++) {
            values[i].resize(size[1]);
            for (int j = 0; j < size[1]; j++)
                values[i][j].resize(size[2], 0);
        }
    }

private:
    void calc_dv() {
        dv = std::abs(vectors[0][0] * vectors[1][1] * vectors[2][2] - vectors[2][0] * vectors[1][1] * vectors[0][2] + vectors[0][1] * vectors[1][2] * vectors[2][0] - vectors[2][1] * vectors[1][2] * vectors[0][0] + vectors[0][2] * vectors[1][0] * vectors[2][1] - vectors[2][2] * vectors[1][0] * vectors[0][1]);
    }

    double dv{};
    int na{};
    bool loaded{};
    std::string comment1;
    std::string comment2;
    i3 size{};
    d3 origin{};
    std::array<d3, 3> vectors{};
    ivec3 values;
    std::filesystem::path path;
    WFN* parent_wavefunction{};
};

using cubed = cube_t<double>;
using cubei = cube_t<int>;

#include "wfn_class.h"
#include "atoms.h"
