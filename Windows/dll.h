#pragma once
#include <string>
#include <vector>
#include <array>
#include <filesystem>
#include <complex>

#ifdef _WIN32
#ifdef BUILDING_DLL
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT __declspec(dllimport)
#endif
#else
#define DLL_EXPORT
#endif

// Suppose we store each face color as an (R,G,B) triple
typedef std::array<int, 3> RGB;
typedef std::vector<double> vec;
typedef std::vector<vec> vec2;
typedef std::vector<vec2> vec3;
typedef std::vector<int> ivec;
typedef std::vector<ivec> ivec2;
typedef std::vector<std::complex<double>> cvec;
typedef std::vector<cvec> cvec2;
typedef std::vector<cvec2> cvec3;
typedef std::vector<std::vector<cvec2>> cvec4;
typedef std::vector<bool> bvec;
typedef std::vector<std::string> svec;
typedef std::vector<std::filesystem::path> pathvec;

// A simple triangle struct composed of three vertices
struct Triangle {
    std::array<double, 3> v1, v2, v3;
    RGB colour;
    int colour_index;
    DLL_EXPORT double calc_area() const;
    DLL_EXPORT double calc_inner_volume() const;
    DLL_EXPORT std::array<double, 3> calc_center() const;
};

struct primitive
{
    int center, type;
    double exp, coefficient;
    double norm_const = -10;
    DLL_EXPORT void normalize();
    DLL_EXPORT void unnormalize();
    DLL_EXPORT double normalization_constant() const;
    DLL_EXPORT primitive();
    DLL_EXPORT primitive(int c, int t, double e, double coef);
    DLL_EXPORT bool operator==(const primitive& other) const;
};

struct basis_set_entry {
    double coefficient;
    double exponent;
    unsigned int type; //1=S; 2=P; 3=D; 4=F; 5=G
    unsigned int shell;
    // Default assignment operator
    DLL_EXPORT basis_set_entry& operator=(const basis_set_entry& rhs);
    DLL_EXPORT basis_set_entry();
    DLL_EXPORT basis_set_entry(double g_coefficient, double g_exponent, unsigned int g_type, unsigned int g_shell);
    // Equality operator
    DLL_EXPORT bool operator==(const basis_set_entry& other) const;
    primitive p;
};

class atom {
public:
    std::string label;
    std::string ID;
    int nr, charge, ECP_electrons;
    double x, y, z;
    vec frac_coords;
    DLL_EXPORT atom();
    DLL_EXPORT atom(const std::string& l, const std::string& id, const int& n, const double& c1, const double& c2, const double& c3, const int& ch);
    DLL_EXPORT atom(const std::string& l, const std::string& id, const int& n, const double& c1, const double& c2, const double& c3, const int& ch, const int& ECP_els);
    DLL_EXPORT atom operator=(const atom& rhs);
    DLL_EXPORT void print_values() const;
    DLL_EXPORT bool push_back_basis_set(const double& coefficient, const double& exponent, const int& type, const int& shell);
    DLL_EXPORT void print_values_long() const;
    DLL_EXPORT bool get_basis_set_loaded() const;
    DLL_EXPORT bool is_anharm() const;
    DLL_EXPORT void assign_ADPs(vec& second, vec& third, vec& fourth);
    DLL_EXPORT void assign_ADPs(vec& second);
    DLL_EXPORT void assign_ADPs(double& Uiso);
    DLL_EXPORT void assign_ID(const std::string& id);
    DLL_EXPORT void set_ID(const std::string& id);
    DLL_EXPORT std::string get_ID() const;
    std::vector<basis_set_entry> basis_set;
    int basis_set_id;
    std::vector<unsigned int> shellcount;
    //The Order is:
    //[0] = second order (U11, U22, U33, U12, U13, U23)
    //[1] = third order  (C111, C112, C113, C122, C123, C133, C222, C223, C233, C333)
    //[2] = fourth order (D1111, D1112, D1113, D1122, D1123, D1133, D1222, D1223, D1233, D1333, D2222, D2223, D2233, D2333, D3333)
    vec2 ADPs;

    DLL_EXPORT bool operator==(const atom& other) const;
};
class MO;
class cube;

class WFN {
public:
    DLL_EXPORT WFN();
    DLL_EXPORT WFN(int given_origin);
    DLL_EXPORT WFN(const std::filesystem::path& filename, const bool& debug = false);
    DLL_EXPORT WFN(const std::filesystem::path& filename, const int g_charge, const int g_mult, const bool& debug = false);
    std::vector<cube> cub;
    std::vector<atom> atoms;
    //--------------------MO handling--------------------------------------------
    DLL_EXPORT bool set_MO_coef(const int& nr_mo, const int& nr_primitive, const double& value);
    DLL_EXPORT const double get_MO_coef(const int& nr_mo, const int& nr_primtive) const;
    DLL_EXPORT const double get_MO_coef_f(const int& nr_mo, const int& nr_primtive) const;
    DLL_EXPORT const double* get_MO_coef_ptr(const int& nr_mo);
    DLL_EXPORT const int get_MO_primitive_count(const int& nr_mo) const;
    DLL_EXPORT bool push_back_MO(const int& nr, const double& occ, const double& ener);
    DLL_EXPORT bool push_back_MO(const int& nr, const double& occ, const double& ener, const int& oper);
    DLL_EXPORT bool push_back_MO(const MO& given);
    DLL_EXPORT void pop_back_MO();
    DLL_EXPORT void delete_MO(const int& nr);
    DLL_EXPORT const double get_MO_energy(const int& mo) const;
    DLL_EXPORT const int get_center(const int& nr) const;
    DLL_EXPORT const int get_type(const int& nr) const;
    DLL_EXPORT const double get_MO_occ(const int& nr) const;
    DLL_EXPORT const int get_MO_op(const int& nr) const;
    DLL_EXPORT void delete_unoccupied_MOs();
    DLL_EXPORT const MO& get_MO(const int& n) const;
    DLL_EXPORT const int get_MO_op_count(const int& op) const;
    DLL_EXPORT const void clear_MOs();

    //--------------------in and output----------------------------------------
    DLL_EXPORT void change_basis_set_name(std::string name);
    DLL_EXPORT bool add_primitive(const int& cent, const int& type, const double& e, double* values);
    DLL_EXPORT bool add_exp(const int& cent, const int& type, const double& e);
    DLL_EXPORT void read_known_wavefunction_format(const std::filesystem::path& fileName, std::ostream& file, const bool debug = false);
    DLL_EXPORT bool read_wfn(const std::filesystem::path& fileName, const bool& debug, std::ostream& file);
    DLL_EXPORT bool read_wfx(const std::filesystem::path& fileName, const bool& debug, std::ostream& file);
    DLL_EXPORT bool read_fchk(const std::filesystem::path& filename, std::ostream& log, const bool debug = false);
    DLL_EXPORT bool read_xyz(const std::filesystem::path& filename, std::ostream& file, const bool debug = false);
    DLL_EXPORT bool read_molden(const std::filesystem::path& filename, std::ostream& file, const bool debug = false);
    DLL_EXPORT bool read_gbw(const std::filesystem::path& filename, std::ostream& file, const bool debug = false, const bool has_ECPs = false);
    DLL_EXPORT bool read_ptb(const std::filesystem::path& filename, std::ostream& file, const bool debug = false);
    DLL_EXPORT bool write_wfn(const std::filesystem::path& fileName, const bool& debug, const bool occupied);
    DLL_EXPORT bool write_xyz(const std::filesystem::path& fileName, const bool& debug = false);
    DLL_EXPORT bool set_path(std::filesystem::path given_path);
    DLL_EXPORT void print_primitive(const int& nr) const;
    DLL_EXPORT void assign_charge(const int& i_charge);
    DLL_EXPORT void assign_multi(const int& i_multi);
    DLL_EXPORT const int get_charge() const;
    DLL_EXPORT const int get_multi() const;
    DLL_EXPORT void set_multi(unsigned int& in);
    DLL_EXPORT void set_charge(const int& in);
    DLL_EXPORT const int get_nex() const;
    DLL_EXPORT const int get_ncen() const;
    DLL_EXPORT const void set_ncen(const int& in);
    DLL_EXPORT const int get_nmo() const;
    DLL_EXPORT const int get_nmo(const bool& only_occ) const;
    DLL_EXPORT const int get_origin() const;
    DLL_EXPORT const int get_ECP_mode() const;
    DLL_EXPORT const std::string get_comment() const;
    DLL_EXPORT const std::string get_CIF_table(const int nr = 0) const;
    DLL_EXPORT const std::string get_basis_set_CIF(const int nr = 0) const;
    DLL_EXPORT void write_wfn_CIF(const std::filesystem::path& filename) const;
    DLL_EXPORT const double get_exponent(int nr) const;
    DLL_EXPORT const unsigned int get_nr_electrons() const;
    DLL_EXPORT const unsigned int get_nr_ECP_electrons() const;
    DLL_EXPORT double count_nr_electrons(void) const;
    DLL_EXPORT const std::string get_centers(const bool& bohr) const;
    DLL_EXPORT const std::string get_basis_set_name() const;
    DLL_EXPORT void set_basis_set_name(const std::string& input);
    DLL_EXPORT const std::filesystem::path get_path() const;
    DLL_EXPORT const std::string hdr(const bool& occupied) const;
    DLL_EXPORT void set_method(const std::string& input);
    DLL_EXPORT const std::string get_method() const;
    DLL_EXPORT bool erase_atom(const int& nr);
    DLL_EXPORT const void list_primitives() const;
    DLL_EXPORT const void list_centers() const;
    DLL_EXPORT bool remove_center(const int& nr);
    DLL_EXPORT bool remove_primitive(const int& nr);
    DLL_EXPORT void change_center(const int& nr);
    DLL_EXPORT void change_type(const int& nr);
    DLL_EXPORT void change_exponent(const int& nr);
    DLL_EXPORT void set_modified();
    DLL_EXPORT const bool get_modified() const;
    DLL_EXPORT void set_d_f_switch(const bool& in);
    DLL_EXPORT const bool get_d_f_switch() const;
    DLL_EXPORT int check_order(const bool& debug) const;
    DLL_EXPORT bool sort_wfn(const int& g_order, const bool& debug);
    DLL_EXPORT void set_dist_switch();
    DLL_EXPORT void set_dist_switch(const bool& g);
    DLL_EXPORT const bool get_dist_switch() const;
    DLL_EXPORT void set_has_ECPs(const bool& in, const bool& apply_to_aotms = true, const int& ECP_mode = 1);
    DLL_EXPORT void set_ECPs(ivec& nr, ivec& elcount);
    DLL_EXPORT const bool get_has_ECPs() const;
    DLL_EXPORT void operator=(const WFN& right);
    DLL_EXPORT int calculate_charge();
    DLL_EXPORT int calculate_charge(std::ostream& file);
    DLL_EXPORT bool guess_multiplicity(std::ostream& file);
    DLL_EXPORT const vec get_norm_const(std::ostream& file, const bool debug = false) const;
    DLL_EXPORT double get_total_energy() const;
    DLL_EXPORT double get_virial_ratio() const;
    DLL_EXPORT vec get_DensityMatrix() const;
    DLL_EXPORT vec get_SpinDensityMatrix() const;
    DLL_EXPORT int get_nfunc() const;
    DLL_EXPORT bool delete_basis_set();
    //-------------------atom handling--------------------------------------------------------------
    DLL_EXPORT const double get_atom_coordinate(const unsigned int& nr, const unsigned int& axis) const;
    DLL_EXPORT const std::string get_atom_label(const unsigned int& nr) const;
    DLL_EXPORT const int get_nr_basis_set_loaded() const;
    DLL_EXPORT const bool get_atom_basis_set_loaded(const int& nr) const;
    DLL_EXPORT const double get_atom_basis_set_exponent(const int& nr_atom, const int& nr_prim) const;
    DLL_EXPORT const double get_atom_basis_set_coefficient(const int& nr_atom, const int& nr_prim) const;
    DLL_EXPORT bool change_atom_basis_set_exponent(const int& nr_atom, const int& nr_prim, const double& value);
    DLL_EXPORT bool change_atom_basis_set_coefficient(const int& nr_atom, const int& nr_prim, const double& value);
    DLL_EXPORT const int get_atom_primitive_count(const int& nr) const;
    DLL_EXPORT const int get_atom_primitive_type(const int& nr_atom, const int& nr_prim) const;
    DLL_EXPORT bool erase_atom_primitive(const unsigned int& nr, const unsigned int& nr_prim);
    DLL_EXPORT const int get_atom_shell_count(const unsigned int& nr) const;
    DLL_EXPORT const int get_atom_shell_primitives(const unsigned int& nr_atom, const unsigned int& nr_shell) const;
    DLL_EXPORT const int get_shell_type(const unsigned int& nr_atom, const unsigned int& nr_shell) const;
    DLL_EXPORT const int get_shell_center(const unsigned int& nr_atom, const unsigned int& nr_shell) const;
    DLL_EXPORT const int get_basis_set_shell(const unsigned int& nr_atom, const unsigned int& nr_prim) const;
    DLL_EXPORT const int get_shell_start(const unsigned int& nr_atom, const unsigned int& nr_shell) const;
    DLL_EXPORT const int get_shell_start_in_primitives(const unsigned int& nr_atom, const unsigned int& nr_shell) const;
    DLL_EXPORT const int get_shell_end(const unsigned int& nr_atom, const unsigned int& nr_shell) const;
    DLL_EXPORT bool push_back_atom(const std::string& label, const double& x, const double& y, const double& z, const int& charge, const std::string& ID = "0000000000000");
    DLL_EXPORT bool push_back_atom(const atom& given);
    DLL_EXPORT bool push_back_atom_basis_set(const int& nr, const double& exponent, const double& coefficient, const int& type, const int& shell);
    DLL_EXPORT void print_atom_long(const int& nr) const;
    DLL_EXPORT const int get_atom_charge(const int& nr) const;
    DLL_EXPORT const unsigned int get_atom_integer_mass(const unsigned int& atomnr) const;
    DLL_EXPORT const double get_atom_real_mass(const int& atomnr) const;
    DLL_EXPORT atom get_atom(const int& nr) const;
    DLL_EXPORT const double compute_dens(const double& Pos1, const double& Pos2, const double& Pos3) const;
    DLL_EXPORT const double compute_dens(const double& Pos1, const double& Pos2, const double& Pos3, vec2& d, vec& phi) const;
    DLL_EXPORT const double compute_spin_dens(const double& Pos1, const double& Pos2, const double& Pos3) const;
    DLL_EXPORT const double compute_spin_dens(const double& Pos1, const double& Pos2, const double& Pos3, vec2& d, vec& phi) const;
    DLL_EXPORT const void computeValues(const std::array<double, 3>& PosGrid, double& Rho, double& normGrad, double* Hess, double& Elf, double& Eli, double& Lap) const;
    DLL_EXPORT const void computeLapELIELF(const std::array<double, 3>& PosGrid, double& Elf, double& Eli, double& Lap) const;
    DLL_EXPORT const void computeELIELF(const std::array<double, 3>& PosGrid, double& Elf, double& Eli) const;
    DLL_EXPORT const void computeLapELI(const std::array<double, 3>& PosGrid, double& Eli, double& Lap) const;
    DLL_EXPORT const double computeLap(const std::array<double, 3>& PosGrid) const;
    DLL_EXPORT const void computeELI(const std::array<double, 3>& PosGrid, double& Eli) const;
    DLL_EXPORT const void computeELF(const std::array<double, 3>& PosGrid, double& Elf) const;
    DLL_EXPORT const double computeMO(const std::array<double, 3>& PosGrid, const int& mo) const;
    DLL_EXPORT const double compute_MO_spherical(const double& Pos1, const double& Pos2, const double& Pos3, const int& MO) const;
    DLL_EXPORT const double computeESP(const std::array<double, 3>& PosGrid, const vec2& d2) const;
    DLL_EXPORT const double computeESP_noCore(const std::array<double, 3>& PosGrid, const vec2& d2) const;
    //----------DM Handling--------------------------------
    DLL_EXPORT bool build_DM(std::string basis_set_path, bool debug = false);
    DLL_EXPORT void push_back_DM(const double& value = 0.0);
    DLL_EXPORT bool set_DM(const int& nr, const double& value = 0.0);
    DLL_EXPORT const double get_DM(const int& nr) const;
    DLL_EXPORT const int get_DM_size() const;
    DLL_EXPORT void resize_DM(const int& size, const double& value = 0.0);
    DLL_EXPORT vec2 get_dm() const;
    //----------S_DM Handling--------------------------------
    DLL_EXPORT void push_back_SDM(const double& value = 0.0);
    DLL_EXPORT bool set_SDM(const int& nr, const double& value = 0.0);
    DLL_EXPORT const double get_SDM(const int& nr) const;
    DLL_EXPORT const int get_SDM_size() const;
    DLL_EXPORT void resize_SDM(const int& size, const double& value = 0.0);
    //-----------Cube handling-------------------------
    DLL_EXPORT bool push_back_cube(const std::string& filepath, const bool& full, const bool& expert = false);
    DLL_EXPORT void push_back_cube(cube given);
    DLL_EXPORT void pop_back_cube();
    //-----------Pointer to members---------------------
    DLL_EXPORT const int* get_ptr_types();
    DLL_EXPORT const int* get_ptr_centers();
    DLL_EXPORT const double* get_ptr_exponents();
    DLL_EXPORT const double* get_ptr_mo_coefficients(const int& mo);
};

class cube {
public:
    DLL_EXPORT cube();
    DLL_EXPORT cube(int x, int y, int z, int g_na = 0, bool grow_values = false);
    DLL_EXPORT cube(const std::filesystem::path& filepath, bool read, WFN& wave, std::ostream& file, const bool expert = false, const bool header = true);
    DLL_EXPORT cube(const int g_na, const std::vector<int>& g_size, const std::vector<double>& g_origin, const std::vector<std::vector<double>>& g_vectors, const std::vector<std::vector<std::vector<double>>>& g_values);
    DLL_EXPORT cube(const cube& given);
    DLL_EXPORT int get_size(int direction) const;
    DLL_EXPORT bool get_loaded() const;
    DLL_EXPORT cube operator + (cube& right) const;
    DLL_EXPORT cube operator - (cube& right) const;
    DLL_EXPORT cube operator * (cube& right) const;
    DLL_EXPORT cube operator / (cube& right) const;
    DLL_EXPORT bool operator += (cube& right);
    DLL_EXPORT bool operator -= (cube& right);
    DLL_EXPORT bool operator *= (cube& right);
    DLL_EXPORT bool operator /= (cube& right);
    DLL_EXPORT void operator = (cube& right);
    DLL_EXPORT bool mask(cube& right);
    DLL_EXPORT bool thresh(cube& right, double thresh = -1234);
    DLL_EXPORT bool negative_mask(cube& right);
    DLL_EXPORT double rrs(cube& right);
    DLL_EXPORT double sum();
    DLL_EXPORT double diff_sum();
    DLL_EXPORT std::vector<double> double_sum();
    DLL_EXPORT std::array<double, 3> get_pos(const int& i, const int& j, const int& k) const;
    DLL_EXPORT double get_interpolated_value(double x, double y, double z) const;
    DLL_EXPORT double get_value(int x, int y, int z) const;
    DLL_EXPORT bool set_value(int x, int y, int z, double value);
    DLL_EXPORT bool read_file(bool full, bool header, bool expert = false);
    DLL_EXPORT bool write_file(bool force = false, bool absolute = false);
    DLL_EXPORT bool write_file(const std::filesystem::path& given_path, bool debug = false);
    DLL_EXPORT bool write_xdgraph(const std::filesystem::path& given_path, bool debug = false);
    DLL_EXPORT bool fractal_dimension(const double stepsize);
    DLL_EXPORT double get_vector(int i, int j) const;
    DLL_EXPORT bool set_vector(int i, int j, double value);
    DLL_EXPORT double get_origin(unsigned int i) const;
    DLL_EXPORT double ewald_sum(const int kMax = 15, const double conv = 5E-3);
    DLL_EXPORT void calc_dv();
    DLL_EXPORT double get_dv() const;
    DLL_EXPORT void set_dv(const double& given);
    DLL_EXPORT bool set_origin(unsigned int i, double value);
    DLL_EXPORT int get_na() const;
    DLL_EXPORT void set_na(int g_na);
    DLL_EXPORT std::filesystem::path super_cube();
    DLL_EXPORT cube super_cube(int x, int y, int z);
    DLL_EXPORT void set_comment1(std::string input);
    DLL_EXPORT void set_comment2(std::string input);
    DLL_EXPORT void set_zero();
    DLL_EXPORT void give_parent_wfn(WFN& given);
    DLL_EXPORT std::string get_comment1() const;
    DLL_EXPORT std::string get_comment2() const;
    DLL_EXPORT std::filesystem::path get_path() const;
    DLL_EXPORT void set_path(const std::filesystem::path& given);
    DLL_EXPORT std::vector<atom> get_parent_wfn_atoms() const;
};



DLL_EXPORT std::vector<Triangle> marchingCubes(const cube& volumeData, const double isoVal, const int subdivisionLevel = 2);
DLL_EXPORT bool writeObj(const std::filesystem::path& filename, const std::vector<Triangle>& triangles);
DLL_EXPORT bool writeColourObj(const std::filesystem::path& filename, std::vector<Triangle>& triangles);
DLL_EXPORT bool writeMTL(const std::string& mtlFilename, std::vector<Triangle>& triangles);
DLL_EXPORT RGB get_colour(const Triangle& t, const cube& volumeData, std::array<std::array<int, 3>, 3> Colourcode, double low_lim, double high_lim);
DLL_EXPORT RGB get_colour(const Triangle& t, double(*func)(const std::array<double, 3>&, const WFN&), const WFN& wavy, std::array<std::array<int, 3>, 3> Colourcode, double low_lim, double high_lim);

DLL_EXPORT double calc_d_i(const std::array<double, 3>& p_t, const WFN& wavy);
DLL_EXPORT bool Calc_Spherical_Dens(cube& CubeSpher, WFN& wavy, double radius, std::ostream& file, bool wrap);
DLL_EXPORT bool Calc_Static_Def(cube& CubeDEF, cube& CubeRho, cube& CubeSpher, WFN& wavy, double radius, std::ostream& file, bool wrap);
DLL_EXPORT bool Calc_Static_Def(cube& CubeDEF, cube& CubeRho, WFN& wavy, double radius, std::ostream& file, bool wrap);
DLL_EXPORT bool Calc_Rho(cube& CubeRho, WFN& wavy, double radius, std::ostream& file, bool wrap);
DLL_EXPORT bool Calc_Prop(cube& CubeRho, cube& CubeRDG, cube& CubeElf, cube& CubeEli, cube& CubeLap, cube& CubeESP, WFN& wavy, double radius, std::ostream& file, bool test, bool wrap);
DLL_EXPORT bool Calc_ESP(cube& CubeESP, WFN& wavy, double radius, bool no_date, std::ostream& file, bool wrap);
DLL_EXPORT bool Calc_Hirshfeld_atom(cube& CubeHirsh, cube& CubeRho, cube& CubeSpher, WFN& wavy, double radius, int atom, std::ostream& file, bool wrap);
DLL_EXPORT bool Calc_Hirshfeld(cube& CubeHirsh, cube& CubeRho, cube& CubeSpher, WFN& wavy, double radius, int atom, std::ostream& file, bool wrap);
DLL_EXPORT bool Calc_MO(cube& MO, int MO_number, WFN& wavy, double radius, std::ostream& file, bool wrap);

DLL_EXPORT void readxyzMinMax_fromWFN(WFN& wavy, double* CoordMinMax, int* NbSteps, double radius, double Increments, bool no_bohr);
