#pragma once
#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
// Windows Header Files
#include <windows.h>

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
typedef std::complex<double> cdouble;
typedef std::vector<double> vec;
typedef std::vector<vec> vec2;
typedef std::vector<vec2> vec3;
typedef std::vector<int> ivec;
typedef std::vector<ivec> ivec2;
typedef std::vector<cdouble> cvec;
typedef std::vector<cvec> cvec2;
typedef std::vector<cvec2> cvec3;
typedef std::vector<std::vector<cvec2>> cvec4;
typedef std::vector<bool> bvec;
typedef std::vector<std::string> svec;
typedef std::vector<std::filesystem::path> pathvec;

// A simple triangle struct composed of three vertices
struct DLL_EXPORT Triangle
{
    double calc_area() const;
    double calc_inner_volume() const;
    std::array<double, 3> calc_center() const;
};

struct DLL_EXPORT primitive
{
    void normalize();
    void unnormalize();
    double normalization_constant() const;
    primitive();
    primitive(int c, int t, double e, double coef);
    bool operator==(const primitive &other) const;
};

struct DLL_EXPORT basis_set_entry
{
    basis_set_entry &operator=(const basis_set_entry &rhs);
};

class DLL_EXPORT atom
{
    atom();
    atom(const std::string &l, const std::string &id, const int &n, const double &c1, const double &c2, const double &c3, const int &ch);
    atom(const std::string &l, const std::string &id, const int &n, const double &c1, const double &c2, const double &c3, const int &ch, const int &ECP_els);
    atom operator=(const atom &rhs);
    void print_values() const;
    bool push_back_basis_set(const double &coefficient, const double &exponent, const int &type, const int &shell);
    void print_values_long() const;
    bool get_basis_set_loaded() const;
    bool is_anharm() const;
    void assign_ADPs(vec &second, vec &third, vec &fourth);
    void assign_ADPs(vec &second);
    void assign_ADPs(double &Uiso);
    void assign_ID(const std::string &id);
    void set_ID(const std::string &id);
    std::string get_ID() const;
    basis_set_entry get_basis_set_entry(const int &nr) const;
    double get_coordinate(const unsigned int &axis) const;
    int get_charge() const;
    int get_ECP_electrons() const;
    void set_charge(const int &ch);
    void set_ECP_electrons(const int &ECP_els);
    void set_coordinate(const unsigned int &axis, const double &value);
    void set_frac_coords(const std::array<double, 3> &frac);
    std::string get_label() const;
    int get_nr() const;
    void set_label(const std::string &l);
    void set_nr(const int &n);
    void clear_basis_set();
    int get_basis_set_size() const;
    double get_basis_set_exponent(const int &nr) const;
    double get_basis_set_coefficient(const int &nr) const;
    void set_basis_set_exponent(const int &nr, const double &value);
    void set_basis_set_coefficient(const int &nr, const double &value);
    void set_basis_set_type(const int &nr, const int &value);
    void set_basis_set_shell(const int &nr, const int &value);
    std::vector<basis_set_entry> get_basis_set() const;
    void erase_basis_set(const unsigned int &nr);
    int get_basis_set_type(const int &nr) const;
    int get_basis_set_shell(const int &nr) const;
    void set_basis_set_id(const int &id);
    int get_basis_set_id() const;
    void set_shellcount(const std::vector<unsigned int> &sc);
    std::vector<unsigned int> get_shellcount() const;
    int get_shellcount_size() const;
    void set_shellcount(const unsigned int &nr, const unsigned int &value);
    unsigned int get_shellcount(const unsigned int &nr) const;
    void clear_shellcount();
    void set_ADPs(const vec2 &adps);
    vec2 get_ADPs() const;
    bool operator==(const atom &other) const;
};

class cube;

class DLL_EXPORT WFN
{
public:
    WFN();
    WFN(int given_origin);
    WFN(const std::filesystem::path &filename, const bool &debug = false);
    WFN(const std::filesystem::path &filename, const int g_charge, const int g_mult, const bool &debug = false);

    //--------------------MO handling--------------------------------------------
    bool set_MO_coef(const int &nr_mo, const int &nr_primitive, const double &value);
    const double get_MO_coef(const int &nr_mo, const int &nr_primtive) const;
    const double get_MO_coef_f(const int &nr_mo, const int &nr_primtive) const;
    const double *get_MO_coef_ptr(const int &nr_mo);
    const int get_MO_primitive_count(const int &nr_mo) const;
    bool push_back_MO(const int &nr, const double &occ, const double &ener);
    bool push_back_MO(const int &nr, const double &occ, const double &ener, const int &oper);
    void pop_back_MO();
    void delete_MO(const int &nr);
    const double get_MO_energy(const int &mo) const;
    const int get_center(const int &nr) const;
    const int get_type(const int &nr) const;
    const double get_MO_occ(const int &nr) const;
    const int get_MO_op(const int &nr) const;
    void delete_unoccupied_MOs();
    const int get_MO_op_count(const int &op) const;
    const void clear_MOs();

    //--------------------in and output----------------------------------------
    void change_basis_set_name(std::string name);
    bool add_primitive(const int &cent, const int &type, const double &e, double *values);
    bool add_exp(const int &cent, const int &type, const double &e);
    void read_known_wavefunction_format(const std::filesystem::path &fileName, std::ostream &file, const bool debug = false);
    bool read_wfn(const std::filesystem::path &fileName, const bool &debug, std::ostream &file);
    bool read_wfx(const std::filesystem::path &fileName, const bool &debug, std::ostream &file);
    bool read_fchk(const std::filesystem::path &filename, std::ostream &log, const bool debug = false);
    bool read_xyz(const std::filesystem::path &filename, std::ostream &file, const bool debug = false);
    bool read_molden(const std::filesystem::path &filename, std::ostream &file, const bool debug = false);
    bool read_gbw(const std::filesystem::path &filename, std::ostream &file, const bool debug = false, const bool has_ECPs = false);
    bool read_ptb(const std::filesystem::path &filename, std::ostream &file, const bool debug = false);
    bool write_wfn(const std::filesystem::path &fileName, const bool &debug, const bool occupied);
    bool write_xyz(const std::filesystem::path &fileName, const bool &debug = false);
    bool set_path(std::filesystem::path given_path);
    void print_primitive(const int &nr) const;
    void assign_charge(const int &i_charge);
    void assign_multi(const int &i_multi);
    const int get_charge() const;
    const int get_multi() const;
    void set_multi(unsigned int &in);
    void set_charge(const int &in);
    const int get_nex() const;
    const int get_ncen() const;
    const void set_ncen(const int &in);
    const int get_nmo() const;
    const int get_nmo(const bool &only_occ) const;
    const int get_origin() const;
    const int get_ECP_mode() const;
    const std::string get_comment() const;
    const std::string get_CIF_table(const int nr = 0) const;
    const std::string get_basis_set_CIF(const int nr = 0) const;
    void write_wfn_CIF(const std::filesystem::path &filename) const;
    const double get_exponent(int nr) const;
    const unsigned int get_nr_electrons() const;
    const unsigned int get_nr_ECP_electrons() const;
    double count_nr_electrons(void) const;
    const std::string get_centers(const bool &bohr) const;
    const std::string get_basis_set_name() const;
    void set_basis_set_name(const std::string &input);
    const std::filesystem::path get_path() const;
    const std::string hdr(const bool &occupied) const;
    void set_method(const std::string &input);
    const std::string get_method() const;
    bool erase_atom(const int &nr);
    const void list_primitives() const;
    const void list_centers() const;
    bool remove_center(const int &nr);
    bool remove_primitive(const int &nr);
    void change_center(const int &nr);
    void change_type(const int &nr);
    void change_exponent(const int &nr);
    void set_modified();
    const bool get_modified() const;
    void set_d_f_switch(const bool &in);
    const bool get_d_f_switch() const;
    int check_order(const bool &debug) const;
    bool sort_wfn(const int &g_order, const bool &debug);
    void set_dist_switch();
    void set_dist_switch(const bool &g);
    const bool get_dist_switch() const;
    void set_has_ECPs(const bool &in, const bool &apply_to_aotms = true, const int &ECP_mode = 1);
    void set_ECPs(ivec &nr, ivec &elcount);
    const bool get_has_ECPs() const;
    void operator=(const WFN &right);
    int calculate_charge();
    int calculate_charge(std::ostream &file);
    bool guess_multiplicity(std::ostream &file);
    const vec get_norm_const(std::ostream &file, const bool debug = false) const;
    double get_total_energy() const;
    double get_virial_ratio() const;
    vec get_DensityMatrix() const;
    vec get_SpinDensityMatrix() const;
    int get_nfunc() const;
    bool delete_basis_set();
    void set_basis_set_ptr(const std::array<std::vector<primitive>, 118> &given);
    const std::array<std::vector<primitive>, 118> *get_basis_set_ptr() const;
    //-------------------atom handling--------------------------------------------------------------
    const double get_atom_coordinate(const unsigned int &nr, const unsigned int &axis) const;
    const std::string get_atom_label(const unsigned int &nr) const;
    const int get_nr_basis_set_loaded() const;
    const bool get_atom_basis_set_loaded(const int &nr) const;
    const double get_atom_basis_set_exponent(const int &nr_atom, const int &nr_prim) const;
    const double get_atom_basis_set_coefficient(const int &nr_atom, const int &nr_prim) const;
    bool change_atom_basis_set_exponent(const int &nr_atom, const int &nr_prim, const double &value);
    bool change_atom_basis_set_coefficient(const int &nr_atom, const int &nr_prim, const double &value);
    const int get_atom_primitive_count(const int &nr) const;
    const int get_atom_primitive_type(const int &nr_atom, const int &nr_prim) const;
    bool erase_atom_primitive(const unsigned int &nr, const unsigned int &nr_prim);
    const int get_atom_shell_count(const unsigned int &nr) const;
    const int get_atom_shell_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_shell_type(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_shell_center(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_basis_set_shell(const unsigned int &nr_atom, const unsigned int &nr_prim) const;
    const int get_shell_start(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_shell_start_in_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_shell_end(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    bool push_back_atom(const std::string &label, const double &x, const double &y, const double &z, const int &charge, const std::string &ID = "0000000000000");
    bool push_back_atom(const atom &given);
    const atom get_atom(const int &nr) const;
    bool push_back_atom_basis_set(const int &nr, const double &exponent, const double &coefficient, const int &type, const int &shell);
    void print_atom_long(const int &nr) const;
    const int get_atom_charge(const int &nr) const;
    const unsigned int get_atom_integer_mass(const unsigned int &atomnr) const;
    const double get_atom_real_mass(const int &atomnr) const;
    std::string get_atom_label(const int &nr) const;
    void set_atom_label(const int &nr, const std::string &label);
    std::vector<atom> get_atoms() const;
    void set_atoms(const std::vector<atom> given);
    int get_atom_ECP_electrons(const int &nr) const;
    void clear_atom_basis_set(const int &nr);
    int get_atom_basis_set_size(const int &nr) const;
    basis_set_entry get_atom_basis_set_entry(const int &nr, const int &bs) const;
    void set_atom_ADPs(const int &nr, const vec2 &adps);
    void set_atom_frac_coords(const int &nr, const std::array<double, 3> &frac);
    //----------Calcualtion of Properties-----------------
    // This second version will use phi[nmo] and d[4][ncen] as scratch instead of allocating new ones
    const double compute_dens(const double &Pos1, const double &Pos2, const double &Pos3) const;
    const double compute_dens(const double &Pos1, const double &Pos2, const double &Pos3, vec2 &d, vec &phi) const;
    const double compute_spin_dens(const double &Pos1, const double &Pos2, const double &Pos3) const;
    const double compute_spin_dens(const double &Pos1, const double &Pos2, const double &Pos3, vec2 &d, vec &phi) const;
    const void computeValues(const std::array<double, 3> &PosGrid, double &Rho, double &normGrad, double *Hess, double &Elf, double &Eli, double &Lap) const;
    const void computeLapELIELF(const std::array<double, 3> &PosGrid, double &Elf, double &Eli, double &Lap) const;
    const void computeELIELF(const std::array<double, 3> &PosGrid, double &Elf, double &Eli) const;
    const void computeLapELI(const std::array<double, 3> &PosGrid, double &Eli, double &Lap) const;
    const double computeLap(const std::array<double, 3> &PosGrid) const;
    const void computeELI(const std::array<double, 3> &PosGrid, double &Eli) const;
    const void computeELF(const std::array<double, 3> &PosGrid, double &Elf) const;
    const double computeMO(const std::array<double, 3> &PosGrid, const int &mo) const;
    const double compute_MO_spherical(const double &Pos1, const double &Pos2, const double &Pos3, const int &MO) const;
    const double computeESP(const std::array<double, 3> &PosGrid, const vec2 &d2) const;
    const double computeESP_noCore(const std::array<double, 3> &PosGrid, const vec2 &d2) const;
    //----------DM Handling--------------------------------
    bool build_DM(std::string basis_set_path, bool debug = false);
    void push_back_DM(const double &value = 0.0);
    bool set_DM(const int &nr, const double &value = 0.0);
    const double get_DM(const int &nr) const;
    const int get_DM_size() const;
    void resize_DM(const int &size, const double &value = 0.0);
    vec2 get_dm() const;
    //----------S_DM Handling--------------------------------
    void push_back_SDM(const double &value = 0.0);
    bool set_SDM(const int &nr, const double &value = 0.0);
    const double get_SDM(const int &nr) const;
    const int get_SDM_size() const;
    void resize_SDM(const int &size, const double &value = 0.0);
    //-----------Cube handling-------------------------
    bool push_back_cube(const std::string &filepath, const bool &full, const bool &expert = false);
    void push_back_cube(cube given);
    void pop_back_cube();
    const cube get_cube(const int &nr) const;
    const int get_cube_count() const;
    std::filesystem::path get_cube_path(const int &nr) const;
    //-----------Pointer to members---------------------
    const int *get_ptr_types();
    const int *get_ptr_centers();
    const double *get_ptr_exponents();
    const double *get_ptr_mo_coefficients(const int &mo);
};
class DLL_EXPORT cube
{
public:
    cube();
    cube(int x, int y, int z, int g_na = 0, bool grow_values = false);
    cube(const std::filesystem::path &filepath, bool read, WFN &wave, std::ostream &file, const bool expert = false, const bool header = true);
    cube(const int g_na, const ivec &g_size, const vec &g_origin, const vec2 &g_vectors, const vec3 &g_values);
    cube(const cube &given);
    int get_size(int direction) const;
    bool get_loaded() const;
    cube operator+(cube &right) const;
    cube operator-(cube &right) const;
    cube operator*(cube &right) const;
    cube operator/(cube &right) const;
    bool operator+=(cube &right);
    bool operator-=(cube &right);
    bool operator*=(cube &right);
    bool operator/=(cube &right);
    void operator=(cube &right);
    bool mask(cube &right);
    bool thresh(cube &right, double thresh = -1234);
    bool negative_mask(cube &right);
    double rrs(cube &right);
    double sum();
    double diff_sum();
    std::vector<double> double_sum();
    inline std::array<double, 3> get_pos(const int &i, const int &j, const int &k) const;
    double get_interpolated_value(double x, double y, double z) const;
    double get_value(int x, int y, int z) const;
    bool set_value(int x, int y, int z, double value);
    bool read_file(bool full, bool header, bool expert = false);
    bool write_file(bool force = false, bool absolute = false);
    bool write_file(const std::filesystem::path &given_path, bool debug = false);
    bool write_xdgraph(const std::filesystem::path &given_path, bool debug = false);
    bool fractal_dimension(const double stepsize);
    double get_vector(int i, int j) const;
    bool set_vector(int i, int j, double value);
    double get_origin(unsigned int i) const;
    double ewald_sum(const int kMax = 15, const double conv = 5E-3);
    void calc_dv();
    double get_dv() const;
    void set_dv(const double &given);
    bool set_origin(unsigned int i, double value);
    int get_na() const;
    void set_na(int g_na);
    std::filesystem::path super_cube();
    cube super_cube(int x, int y, int z);
    void set_comment1(std::string input);
    void set_comment2(std::string input);
    void set_zero();
    void give_parent_wfn(WFN &given);
    std::string get_comment1() const;
    std::string get_comment2() const;
    std::filesystem::path get_path() const;
    void set_path(const std::filesystem::path &given);
    std::vector<atom> get_parent_wfn_atoms() const;
};

std::vector<Triangle> DLL_EXPORT marchingCubes(const cube &volumeData, const double isoVal, const int subdivisionLevel = 2);
bool DLL_EXPORT writeObj(const std::filesystem::path &filename, const std::vector<Triangle> &triangles);
bool DLL_EXPORT writeColourObj(const std::filesystem::path &filename, std::vector<Triangle> &triangles);
bool DLL_EXPORT writeMTL(const std::string &mtlFilename, std::vector<Triangle> &triangles);
void DLL_EXPORT get_colour(Triangle &t, const cube &volumeData, std::array<std::array<int, 3>, 3> Colourcode, double low_lim, double high_lim);
void DLL_EXPORT get_colour(Triangle &t, double (*func)(const std::array<double, 3> &, const WFN &), const WFN &wavy, std::array<std::array<int, 3>, 3> Colourcode, double low_lim, double high_lim);

double DLL_EXPORT calc_d_i(const std::array<double, 3> &p_t, const WFN &wavy);
bool DLL_EXPORT Calc_Spherical_Dens(cube &CubeSpher, WFN &wavy, double radius, std::ostream &file, bool wrap);
bool DLL_EXPORT Calc_Static_Def(cube &CubeDEF, cube &CubeRho, cube &CubeSpher, WFN &wavy, double radius, std::ostream &file, bool wrap);
bool DLL_EXPORT Calc_Static_Def(cube &CubeDEF, cube &CubeRho, WFN &wavy, double radius, std::ostream &file, bool wrap);
bool DLL_EXPORT Calc_Rho(cube &CubeRho, WFN &wavy, double radius, std::ostream &file, bool wrap);
bool DLL_EXPORT Calc_Prop(cube &CubeRho, cube &CubeRDG, cube &CubeElf, cube &CubeEli, cube &CubeLap, cube &CubeESP, WFN &wavy, double radius, std::ostream &file, bool test, bool wrap);
bool DLL_EXPORT Calc_ESP(cube &CubeESP, WFN &wavy, double radius, bool no_date, std::ostream &file, bool wrap);
bool DLL_EXPORT Calc_Hirshfeld_atom(cube &CubeHirsh, cube &CubeRho, cube &CubeSpher, WFN &wavy, double radius, int atom, std::ostream &file, bool wrap);
bool DLL_EXPORT Calc_Hirshfeld(cube &CubeHirsh, cube &CubeRho, cube &CubeSpher, WFN &wavy, double radius, int atom, std::ostream &file, bool wrap);
bool DLL_EXPORT Calc_MO(cube &MO, int MO_number, WFN &wavy, double radius, std::ostream &file, bool wrap);

void DLL_EXPORT readxyzMinMax_fromWFN(WFN &wavy, double *CoordMinMax, int *NbSteps, double radius, double Increments, bool no_bohr);
