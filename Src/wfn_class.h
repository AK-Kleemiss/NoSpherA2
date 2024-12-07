#pragma once

#include "atoms.h"
#include "cube.h"

#include <vector>
#include <string>
#include <fstream>
#include <array>

class MO;

class WFN
{
private:
    int ncen;
    int nfunc;
    int nmo;
    int nex;
    int charge;
    int ECP_m;
    unsigned int multi;
    int origin; // 0=NOT_YET_DEFINED; 1=CRYSTAL; 2=WFN; 3=CUBE; 4=FFN; 5=FCHK; 6=WFX; 7=XYZ; 8=Molden; 9=gbw
    double total_energy;
    double virial_ratio;
    std::string basis_set_name;
    std::string comment;
    std::string path;
    std::string method;

    std::vector<MO> MOs;
    ivec centers;
    // For cartesian the order of types is:
    // 1                                            = S,
    // 2,3,4                                        = X,Y,Z
    // 5,6,7,8,9,10                                 = XX, YY, ZZ, XY, XZ, YZ
    // 11,12,13,14,15,16,17,118,19,20               = XXX, YYY, ZZZ, XXY, XXZ, YYZ, XYY, XZZ, YZZ, XYZ etc..
    // 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 = XXXX YYYY ZZZZ XXXY XXXZ XYYY YYYZ XZZZ YZZZ XXYY XXZZ YYZZ XXYZ XYYZ XYZZ
    //                                                0     1     2   3     4   5     6   7     8   9     10  11  12    13    14
    // for sphericals the order is:
    // 1 = S
    // 2,3,4 = 0, 1, -1
    // 5,6,7,8,9 = 0, +1, -1, +2, -2
    // 10,11,12,13,14,15,16 = 0, +1, -1, +2, -2, +3, -3 etc...
    ivec types;
    vec exponents;
    vec DensityMatrix;
    vec SpinDensityMatrix;
    const std::array<std::vector<primitive>, 118> *basis_set;

    bool erase_center(const int &g_nr);
    bool erase_type(const int &nr);
    bool erase_exponent(const int &nr);
    bool push_back_center(const int &cent);
    bool push_back_type(const int &type);
    bool push_back_exponent(const double &e);
    void assign_MO_coefs(const int &nr, vec &values);
    void push_back_MO_coef(const int& nr, const double& value);
    bool modified;
    bool d_f_switch; // true if spherical harmonics are used for the basis set
    bool distance_switch;
    bool has_ECPs;
    // precomputed factors and helper functions for ESP calc
    long long int pre[9][5][5][9];
    void fill_pre();
    long long int Afac_pre[9][5][9];
    void fill_Afac_pre();
    const double fj(int &j, int &l, int &m, double &aa, double &bb) const;
    const double Afac(int &l, int &r, int &i, double &PC, double &gamma, double &fjtmp) const;
    const double compute_dens_cartesian(const double &Pos1, const double &Pos2, const double &Pos3, vec2 &d, vec &phi) const;
    const double compute_spin_dens_cartesian(const double &Pos1, const double &Pos2, const double &Pos3, vec2 &d, vec &phi) const;
    const double compute_dens_spherical(const double &Pos1, const double &Pos2, const double &Pos3, vec2 &d, vec &phi) const;

public:
    WFN();
    WFN(int given_origin);
    std::vector<cube> cub;
    std::vector<atom> atoms;

    //--------------------MO handling--------------------------------------------
    bool set_MO_coef(const int &nr_mo, const int &nr_primitive, const double &value);
    const double get_MO_coef(const int &nr_mo, const int &nr_primtive) const;
    const double get_MO_coef_f(const int &nr_mo, const int &nr_primtive) const;
    const double *get_MO_coef_ptr(const int &nr_mo);
    const int get_MO_primitive_count(const int &nr_mo) const;
    bool push_back_MO(const int &nr, const double &occ, const double &ener);
    bool push_back_MO(const int &nr, const double &occ, const double &ener, const int &oper);
    bool push_back_MO(const MO &given);
    void pop_back_MO();
    void delete_MO(const int &nr);
    const double get_MO_energy(const int &mo) const;
    const int get_center(const int &nr) const { return centers[nr]; };
    const int get_type(const int &nr) const { return types[nr]; };
    const double get_MO_occ(const int &nr) const;
    const int get_MO_op(const int &nr) const;
    void delete_unoccupied_MOs();
    const MO &get_MO(const int &n) const;
    const int get_MO_op_count(const int &op) const;
    const void clear_MOs();

    //--------------------in and output----------------------------------------
    void change_basis_set_name(std::string name) { basis_set_name = name; };
    bool add_primitive(const int &cent, const int &type, const double &e, double *values);
    bool add_exp(const int &cent, const int &type, const double &e);
    void read_known_wavefunction_format(const std::string &fileName, std::ostream &file, const bool debug = false);
    bool read_wfn(const std::string &fileName, const bool &debug, std::ostream &file);
    bool read_wfx(const std::string &fileName, const bool &debug, std::ostream &file);
    bool read_fchk(const std::string &filename, std::ostream &log, const bool debug = false);
    bool read_xyz(const std::string &filename, std::ostream &file, const bool debug = false);
    bool read_molden(const std::string &filename, std::ostream &file, const bool debug = false);
    bool read_gbw(const std::string &filename, std::ostream &file, const bool debug = false, const bool has_ECPs = false);
    bool read_ptb(const std::string &filename, std::ostream &file, const bool debug = false);
    bool write_wfn(const std::string &fileName, const bool &debug, const bool occupied);
    bool write_xyz(const std::string& fileName, const bool& debug = false);
    bool set_path(std::string given_path)
    {
        path = given_path;
        return true;
    };
    void print_primitive(const int &nr) const;
    void assign_charge(const int &i_charge) { charge = i_charge; };
    void assign_multi(const int &i_multi) { multi = i_multi; };
    const int get_charge() const { return charge; };
    const int get_multi() const { return multi; };
    void set_multi(unsigned int &in) { multi = in; };
    void set_charge(const int &in) { charge = in; };
    const int get_nex() const { return nex; };
    const int get_ncen() const { return ncen; };
    const int get_nmo() const { return nmo; };
    const int get_nmo(const bool &only_occ) const;
    const int get_origin() const { return origin; };
    const int get_ECP_mode() const { return ECP_m; };
    const std::string get_comment() const { return comment; };
    const std::string get_CIF_table(const int nr = 0) const;
    const std::string get_basis_set_CIF(const int nr = 0) const;
    void write_wfn_CIF(const std::string &filename) const;
    const double get_exponent(int nr) const { return exponents[nr]; };
    const unsigned int get_nr_electrons() const;
    const unsigned int get_nr_ECP_electrons() const;
    double count_nr_electrons(void) const;
    const std::string get_centers(const bool &bohr) const;
    const std::string get_basis_set_name() const { return basis_set_name; };
    void set_basis_set_name(const std::string &input) { basis_set_name = input; };
    const std::string get_path() const { return path; };
    const std::string hdr(const bool &occupied) const;
    void set_method(const std::string &input) { method = input; };
    const std::string get_method() const { return method; };
    bool erase_atom(const int &nr);
    const void list_primitives() const;
    const void list_centers() const;
    bool remove_center(const int &nr);
    bool remove_primitive(const int &nr);
    void change_center(const int &nr);
    void change_type(const int &nr);
    void change_exponent(const int &nr);
    void set_modified() { modified = true; };
    const bool get_modified() const { return modified; };
    void set_d_f_switch(const bool &in) { d_f_switch = in; };
    const bool get_d_f_switch() const { return d_f_switch; };
    int check_order(const bool &debug) const;
    bool sort_wfn(const int &g_order, const bool &debug);
    void set_dist_switch() { distance_switch = true; };
    void set_dist_switch(const bool &g) { distance_switch = g; };
    const bool get_dist_switch() const { return distance_switch; };
    void set_has_ECPs(const bool &in, const bool &apply_to_aotms = true, const int &ECP_mode = 1);
    void set_ECPs(ivec &nr, ivec &elcount);
    const bool get_has_ECPs() const { return has_ECPs; };
    void operator=(const WFN &right);
    int calculate_charge();
    int calculate_charge(std::ostream &file);
    bool guess_multiplicity(std::ostream &file);
    const vec get_norm_const(std::ostream &file, const bool debug = false) const;
    double get_total_energy() const { return total_energy; };
    double get_virial_ratio() const { return virial_ratio; };
    vec get_DensityMatrix() const { return DensityMatrix; };
    vec get_SpinDensityMatrix() const { return SpinDensityMatrix; };
    int get_nfunc() const { return nfunc; };
    /**
     * Deletes the basis set information from the given WFN object.
     *
     * @param wave The WFN object to delete the basis set information from.
     * @return Returns true if the basis set information was successfully deleted, false otherwise.
     */
    bool delete_basis_set();
    void set_basis_set_ptr(const std::array<std::vector<primitive>, 118> &given) { basis_set = &given; };
    const std::array<std::vector<primitive>, 118> *get_basis_set_ptr() const { return basis_set; };
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
    const int get_atom_primitive_type(const int &nr_atom, const int &nr_prim) const
    {
        if (nr_atom < atoms.size() && nr_atom >= 0 && nr_prim < atoms[nr_atom].basis_set.size() && nr_prim >= 0)
            return atoms[nr_atom].basis_set[nr_prim].type;
        else
            return -1;
    };
    bool erase_atom_primitive(const unsigned int &nr, const unsigned int &nr_prim);
    const int get_atom_shell_count(const unsigned int &nr) const;
    const int get_atom_shell_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_shell_type(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_shell_center(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_basis_set_shell(const unsigned int &nr_atom, const unsigned int &nr_prim) const;
    const int get_shell_start(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_shell_start_in_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    const int get_shell_end(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    bool push_back_atom(const std::string &label, const double &x, const double &y, const double &z, const int &charge, const std::string& ID = "0000000000000");
    bool push_back_atom(const atom &given);
    // atom get_atom(int nr) { if(nr >= 0 && nr < ncen) return atoms[nr]; else return atom(); };
    bool push_back_atom_basis_set(const int &nr, const double &exponent, const double &coefficient, const int &type, const int &shell)
    {
        if (nr <= ncen && nr >= 0)
            return atoms[nr].push_back_basis_set(exponent, coefficient, type, shell);
        else
            return false;
    };
    void print_atom_long(const int &nr) const
    {
        if (nr <= ncen && nr >= 0)
            atoms[nr].print_values_long();
    };
    const int get_atom_charge(const int &nr) const;
    const unsigned int get_atom_integer_mass(const unsigned int &atomnr) const;
    const double get_atom_real_mass(const int &atomnr) const;
    atom get_atom(const int &nr) const;
    //----------Calcualtion of Properties-----------------
    // double compute_dens(const double* PosGrid, const int atom = -1);
    // This second version will use phi[nmo] and d[4][ncen] as scratch instead of allocating new ones
    const double compute_dens(const double &Pos1, const double &Pos2, const double &Pos3) const;
    const double compute_dens(const double &Pos1, const double &Pos2, const double &Pos3, vec2 &d, vec &phi) const;
    const double compute_spin_dens(const double &Pos1, const double &Pos2, const double &Pos3) const;
    const double compute_spin_dens(const double &Pos1, const double &Pos2, const double &Pos3, vec2 &d, vec &phi) const;
    const void computeValues(const std::array<double, 3>& PosGrid, double &Rho, double &normGrad, double *Hess, double &Elf, double &Eli, double &Lap) const;
    const void computeLapELIELF(const std::array<double, 3>& PosGrid, double &Elf, double &Eli, double &Lap) const;
    const void computeELIELF(const std::array<double, 3>& PosGrid, double &Elf, double &Eli) const;
    const void computeLapELI(const std::array<double, 3>& PosGrid, double &Eli, double &Lap) const;
    const double computeLap(const std::array<double, 3>& PosGrid) const;
    const void computeELI(const std::array<double,3>&PosGrid, double &Eli) const;
    const void computeELF(const std::array<double,3>&PosGrid, double &Elf) const;
    const double computeMO(const std::array<double, 3>& PosGrid, const int &mo) const;
    const double compute_MO_spherical(const double &Pos1, const double &Pos2, const double &Pos3, const int &MO) const;
    const double computeESP(const std::array<double, 3>& PosGrid, const vec2 &d2) const;
    const double computeESP_noCore(const std::array<double, 3>& PosGrid, const vec2 &d2) const;
    //----------DM Handling--------------------------------
    bool build_DM(std::string basis_set_path, bool debug = false);
    void push_back_DM(const double &value = 0.0);
    bool set_DM(const int &nr, const double &value = 0.0);
    const double get_DM(const int &nr) const;
    const int get_DM_size() const { return (int)DensityMatrix.size(); };
    void resize_DM(const int &size, const double &value = 0.0);
    //----------S_DM Handling--------------------------------
    void push_back_SDM(const double &value = 0.0);
    bool set_SDM(const int &nr, const double &value = 0.0);
    const double get_SDM(const int &nr) const;
    const int get_SDM_size() const { return (int)SpinDensityMatrix.size(); };
    void resize_SDM(const int &size, const double &value = 0.0);
    //-----------Cube handling-------------------------
    bool push_back_cube(const std::string &filepath, const bool &full, const bool &expert = false);
    void push_back_cube(cube given) { cub.push_back(given); };
    void pop_back_cube();
    //-----------Pointer to members---------------------
    const int *get_ptr_types() { return &types[0]; };
    const int *get_ptr_centers() { return &centers[0]; };
    const double *get_ptr_exponents() { return &exponents[0]; };
    const double *get_ptr_mo_coefficients(const int &mo);
};

#include "mo_class.h"