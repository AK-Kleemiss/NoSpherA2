#pragma once

#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <filesystem>

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
  cube(int x, int y, int z, int g_na = 0, bool grow_values = false);
  cube(const std::filesystem::path& filepath, bool read, WFN& wave, std::ostream& file, const bool expert = false, const bool header = true);
  cube(const int g_na, const std::vector <int>& g_size, const std::vector <double>& g_origin, const std::vector < std::vector<double> >& g_vectors, const std::vector<std::vector<std::vector<double> > >& g_values);
  cube(const cube& given);
  int get_size(int direction) const;
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
  inline std::array<double, 3> get_pos(const int& i, const int& j, const int& k) const {
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
  bool set_vector(int i, int j, double value);
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
  void give_parent_wfn(WFN& given) { parent_wavefunction = &given; };
  std::string get_comment1() const { return comment1; };
  std::string get_comment2() const { return comment2; };
  std::filesystem::path get_path() const { return path; };
  void set_path(const std::filesystem::path& given) { path = given; };
  std::vector<atom> get_parent_wfn_atoms() const;
  bool read_values(std::ifstream& file);
  double jaccard(const cube& right) const;
private:
  double dv;
  int na;
  bool loaded;
  std::string comment1;
  std::string comment2;
  std::array <int, 3> size;
  std::array <double, 3> origin;
  std::array < std::array <double, 3>, 3> vectors;
  std::vector < std::vector < std::vector <double> > > values;
  std::filesystem::path path;
  WFN* parent_wavefunction;
};

#include "wfn_class.h"
#include "atoms.h"
