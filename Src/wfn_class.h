#pragma once

#include "atoms.h"
#include "cube.h"

#include <vector>
#include <string>
#include <array>
#include <filesystem>

class MO;

/**
 * @class WFN
 * @brief Container for a quantum-mechanical wavefunction including atoms, basis primitives, molecular orbitals and derived properties.
 *
 * Responsibilities:
 *  - Store atomic coordinates, charges, basis set primitives and MO coefficients
 *  - Read / write multiple file formats (wfn, wfx, fchk, molden, gbw, xyz, xtb/ptb)
 *  - Provide accessors & mutators for structural and electronic data
 *  - Construct / manipulate density and spin-density matrices
 *  - Compute real‑space properties (density, spin density, Laplacian, ELF/ELI, ESP, MOs)
 *  - Maintain ordering / normalization / contraction handling of primitives
 */
class WFN
{
private:
	//Number of centers/atoms present in the wavefunction
    int ncen;
	//Number of basis functions/primitive Gaussians present in the wavefunction
    int nfunc;
	//Number of molecular orbitals present in the wavefunction
    int nmo;
	//Number of primitive exponents present in the wavefunction
    int nex;
	//Total charge of the wavefunction
    int charge;
	//Number of the ECPs modi, as specified in the ECPs_corrections.h file
    int ECP_m;
    //Multiplicity of the wavefunction
    unsigned int multi;
    //Number of software/Filetype that was used to generate the wavefunction
    // 0=NOT_YET_DEFINED/UNKNOWN; 1=CRYSTAL; 2=WFN; 3=CUBE; 4=FFN; 5=FCHK; 6=WFX; 7=XYZ; 8=Molden; 9=gbw
    int origin;
    //Store the total energy of the wavefunction
    double total_energy;
	//Store the virial ratio of the wavefunction (if available)
    double virial_ratio;
	//Store the dipole moment of the wavefunction (if available)
    std::string basis_set_name;
	// Store some comments that might be in the wfn file
    std::string comment;
	// The path the file was read or written from/to
    std::filesystem::path path;
	// The QM method that was used to generate the wavefunction
    std::string method;
	// Vector of molecular orbitals
    std::vector<MO> MOs;
    // Vector of centeres that primitives are base on
    ivec centers;
	// Vector of types of primitives
    // For cartesian the order of types is:
    // 1                                            = S,
    // 2,3,4                                        = X,Y,Z
    // 5,6,7,8,9,10                                 = XX, YY, ZZ, XY, XZ, YZ
    // 11,12,13,14,15,16,17,18,19,20                = XXX, YYY, ZZZ, XXY, XXZ, YYZ, XYY, XZZ, YZZ, XYZ etc..
    // 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 = XXXX YYYY ZZZZ XXXY XXXZ XYYY YYYZ XZZZ YZZZ XXYY XXZZ YYZZ XXYZ XYYZ XYZZ
    //                                                0     1     2   3     4   5     6   7     8   9     10  11  12    13    14
    // for sphericals the order is:
    // 1 = S
    // 2,3,4 = 0, 1, -1
    // 5,6,7,8,9 = 0, +1, -1, +2, -2
    // 10,11,12,13,14,15,16 = 0, +1, -1, +2, -2, +3, -3 etc...
    ivec types;

	// Vector of exponents of primitives
    vec exponents;
    vec UT_DensityMatrix;
    vec UT_SpinDensityMatrix;
	// Density Matrix in mdarray
    dMatrix2 DM;
	// basis set definition (118 elements for each element in the periodic table)
    std::shared_ptr<std::array<std::vector<primitive>, 118>> basis_set;
	// Vector of cube files associated with the wavefunction (e.g. for the density or MOs)
    std::vector<cube> cub;
	// Vector of atoms/centers in the wavefunction
    std::vector<atom> atoms;

	// remove a center from the centers vector 
	// CAREFUL: also need to remove all primitives associated with that center, as well as the associated coefficients in each MO and reduce nex accordingly
	// @param g_nr the index of the center to be removed
    // @return true if successful
    bool erase_center(const int &g_nr);
    // remove a type from the type vector 
    // CAREFUL: also need to remove all primitives associated with that center, as well as the associated coefficients in each MO and reduce nex accordingly
    // @param g_nr the index of the type to be removed
    // @return true if successful
    bool erase_type(const int &nr);
    // remove an exponent from the exponents vector 
    // CAREFUL: also need to remove all primitives associated with that center, as well as the associated coefficients in each MO and reduce nex accordingly
    // @param g_nr the index of the exponent to be removed
	// @return true if successful
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
    bool isBohr = false; // True if the coordinates of the atoms are given in Bohr
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
    /** @name Constructors */
    ///@{
    /** Default constructor creates an empty wavefunction object. */
    WFN();
    /** Construct empty WFN with an explicit origin/filetype code. @param given_origin origin code */
    WFN(int given_origin);
    /** Construct by reading a file, auto-detecting filetype. @param filename path to wavefunction file @param debug enable verbose logging */
    WFN(const std::filesystem::path& filename, const bool& debug = false);
    /** Construct with forced charge / multiplicity while reading a file. */
    WFN(const std::filesystem::path& filename, const int g_charge, const int g_mult, const bool& debug = false);
    ///@}

    //--------------------MO handling--------------------------------------------
    /** Set a MO primitive coefficient. @return true if successful */
    bool set_MO_coef(const int &nr_mo, const int &nr_primitive, const double &value);
    /** Get a MO primitive coefficient (double precision / internal). */
    const double& get_MO_coef(const int &nr_mo, const int &nr_primtive) const;
    /** Get a MO primitive coefficient (fast / possibly float). */
    const double& get_MO_coef_f(const int &nr_mo, const int &nr_primtive) const;
    /** Pointer to raw coefficients of one MO (size = nex). */
    const double *get_MO_coef_ptr(const int &nr_mo);
    /** Number of primitives for a given MO (should equal nex). */
    const int get_MO_primitive_count(const int &nr_mo) const;
    /** Append a new MO (no operator / spin id). */
    bool push_back_MO(const int &nr, const double &occ, const double &ener);
    /** Append a new MO with operator id (e.g. spin channel). */
    bool push_back_MO(const int &nr, const double &occ, const double &ener, const int &oper);
    /** Append an existing MO object (copy). */
    bool push_back_MO(const MO &given);
    /** Remove last MO. */
    void pop_back_MO();
    /** Delete MO at index (reorders container). */
    void delete_MO(const int &nr);
    /** MO energy accessor. */
    const double& get_MO_energy(const int &mo) const;
    /** Primitive center index (1-based) for primitive nr. */
    const int& get_center(const int &nr) const { return centers[nr]; };
    /** Primitive type code for primitive nr. */
    const int& get_type(const int &nr) const { return types[nr]; };
    /** MO occupation number. */
    const double& get_MO_occ(const int &nr) const;
    /** Operator (spin / symmetry block) id for MO. */
    const int& get_MO_op(const int &nr) const;
    /** Remove all fully unoccupied MOs (occ==0). */
    void delete_unoccupied_MOs();
    /** Const reference to MO object. */
    const MO &get_MO(const int &n) const;
    /** Count MOs having a specified operator id. */
    const int get_MO_op_count(const int &op) const;
    /** Clear all MOs and reset counter. */
    const void clear_MOs();
	/** Get maximum absolute MO coefficient (for cutoff checks). */
    const double get_maximum_MO_coefficient(bool occu = true) const;

    //--------------------in and output----------------------------------------
    /** Change stored basis set name label. */
    void change_basis_set_name(std::string name) { basis_set_name = name; };
    /** Add fully specified primitive + coefficient values (one per MO). */
    bool add_primitive(const int &cent, const int &type, const double &e, double *values);
    /** Add primitive meta information without coefficients (e.g. building structure). */
    bool add_exp(const int &cent, const int &type, const double &e);
    /** Auto-detect file type and read wavefunction. */
    void read_known_wavefunction_format(const std::filesystem::path&fileName, std::ostream &file, const bool debug = false);
    /** Read legacy .wfn /.ffn file. */
    bool read_wfn(const std::filesystem::path &fileName, const bool &debug, std::ostream &file);
    /** Read .wfx file. */
    bool read_wfx(const std::filesystem::path&fileName, const bool &debug, std::ostream &file);
    /** Read Gaussian formatted checkpoint .fchk. */
    bool read_fchk(const std::filesystem::path&filename, std::ostream &log, const bool debug = false);
    /** Read .xyz geometry (no MOs). */
    bool read_xyz(const std::filesystem::path&filename, std::ostream &file, const bool debug = false);
    /** Read Molden format (.molden). */
    bool read_molden(const std::filesystem::path&filename, std::ostream &file, const bool debug = false);
    /** Read ORCA .gbw binary file. */
    bool read_gbw(const std::filesystem::path&filename, std::ostream &file, const bool debug = false, const bool has_ECPs = false);
    /** Read xTB / pTB binary orbital file. */
    bool read_ptb(const std::filesystem::path&filename, std::ostream &file, const bool debug = false);
    /** Write current wavefunction to .wfn file (optionally only occupied). */
    bool write_wfn(const std::filesystem::path&fileName, const bool &debug, const bool occupied) const;
    /** Write atomic geometry to .xyz file. */
    bool write_xyz(const std::filesystem::path& fileName);
    /** Set internal path field. */
    bool set_path(std::filesystem::path given_path)
    {
        path = given_path;
        return true;
    };
    /** Print primitive information to stdout. */
    void print_primitive(const int &nr) const;
    /** Assign total charge (overwrites existing). */
    void assign_charge(const int &i_charge) { charge = i_charge; };
    /** Assign multiplicity (spin (2S+1)). */
    void assign_multi(const int &i_multi) { multi = i_multi; };
    /** Get total charge. */
    const int& get_charge() const { return charge; };
    /** Get spin multiplicity. */
    const unsigned int& get_multi() const { return multi; };
    /** Set multiplicity by reference (legacy). */
    void set_multi(unsigned int &in) { multi = in; };
    /** Set charge. */
    void set_charge(const int &in) { charge = in; };
    /** Number of primitives (nex). */
    const int& get_nex() const { return nex; };
    /** Number of centers / atoms. */
    const int& get_ncen() const { return ncen; };
    /** Manually set number of centers (use with care). */
    const void set_ncen(const int& in) { ncen = in; };
    /** Number of MOs (including unoccupied). */
    const int& get_nmo() const { return nmo; };
    /** Number of (optionally only occupied) MOs. */
    const int get_nmo(const bool &only_occ) const;
    /** Origin/file type code. */
    const int& get_origin() const { return origin; };
    /** ECP mode (def2/xTB/pTB etc.). */
    const int& get_ECP_mode() const { return ECP_m; };
    /** Freeform comment header. */
    const std::string& get_comment() const { return comment; };
    //void write_wfn_CIF(const std::filesystem::path &filename) const;
    /** Get primitive exponent by index. */
    const double& get_exponent(int nr) const { return exponents[nr]; };
    /** Number of electrons (Z total - charge - ECP core electrons). */
    const unsigned int get_nr_electrons() const;
    /** Total number of core electrons represented by ECPs. */
    const unsigned int get_nr_ECP_electrons() const;
    /** Sum of MO occupations (for consistency checks). */
    double count_nr_electrons(void) const;
    /** Human-readable string listing centers and positions. */
    const std::string get_centers(const bool &bohr) const;
    /** Basis set name accessor. */
    const std::string& get_basis_set_name() const { return basis_set_name; };
    /** Set basis set name string. */
    void set_basis_set_name(const std::string &input) { basis_set_name = input; };
    /** Path to last loaded/saved file. */
    const std::filesystem::path& get_path() const { return path; };
    /** Construct file header line for .wfn writing. */
    const std::string hdr(const bool &occupied) const;
    /** Set electronic structure method label (e.g. DFT functional). */
    void set_method(const std::string &input) { method = input; };
    /** Get method label. */
    const std::string& get_method() const { return method; };
    /** Erase atom (and adjust counts). */
    bool erase_atom(const int &nr);
    /** List primitives to stdout. */
    const void list_primitives() const;
    /** List centers to stdout. */
    const void list_centers() const;
    /** Remove a center and all associated primitives. */
    bool remove_center(const int &nr);
    /** Remove primitive by index (adjust MO coefficients). */
    bool remove_primitive(const int &nr);
    /** Interactive change of primitive center (console). */
    void change_center(const int &nr);
    /** Interactive change of primitive type (console). */
    void change_type(const int &nr);
    /** Interactive change of primitive exponent (console). */
    void change_exponent(const int &nr);
    /** Mark object as modified. */
    void set_modified() { modified = true; };
    /** Modification flag. */
    const bool& get_modified() const { return modified; };
    /** Enable / disable spherical harmonic interpretation for d/f functions. */
    void set_d_f_switch(const bool &in) { d_f_switch = in; };
    /** Query spherical harmonic switch. */
    const bool& get_d_f_switch() const { return d_f_switch; };
    /** Inspect ordering of primitives vs expected canonical orders. @return composite order code */
    int check_order(const bool &debug) const;
    /** Sort primitives into canonical order based on order code. */
    bool sort_wfn(const int &g_order, const bool &debug);
    /** Force distance switch (enables distance-dependent algorithms). */
    void set_dist_switch() { distance_switch = true; };
    /** Explicitly set distance switch flag. */
    void set_dist_switch(const bool &g) { distance_switch = g; };
    /** Query distance switch. */
    const bool& get_dist_switch() const { return distance_switch; };
    /** Set ECP usage flag and optionally populate atom core electron counts. */
    void set_has_ECPs(const bool &in, const bool &apply_to_aotms = true, const int &ECP_mode = 1);
    /** Manually assign ECP core electron counts for given atomic numbers. */
    void set_ECPs(ivec &nr, ivec &elcount);
    /** Query whether ECPs are active. */
    const bool& get_has_ECPs() const { return has_ECPs; };
    /** Copy assignment (deep copy except shared basis definition pointer). */
    void operator=(const WFN &right);
    /** Compute formal charge from atom charges - electron count (ignores stored charge). */
    int calculate_charge();
    /** Same as calculate_charge() with logging. */
    int calculate_charge(std::ostream &file);
    /** Guess multiplicity based on electron parity (sets multi). */
    bool guess_multiplicity(std::ostream &file);
    /** Compute normalization constants for basis primitives / contractions. */
    const vec get_norm_const(std::ostream &file, const bool debug = false) const;
    /** Total SCF energy (if available). */
    const double& get_total_energy() const { return total_energy; };
    /** Virial ratio (if available). */
    const double& get_virial_ratio() const { return virial_ratio; };
    /** Upper-triangular density matrix (linear storage). */
    const vec& get_DensityMatrix() const { return UT_DensityMatrix; };
    /** Upper-triangular spin density matrix (linear storage). */
    const vec& get_SpinDensityMatrix() const { return UT_SpinDensityMatrix; };
    /** Number of basis functions (if tracked). */
    const int& get_nfunc() const { return nfunc; };
    /**
     * Deletes the basis set information from all atoms (shell definitions & primitives).
     * @return true if the basis set information was successfully deleted, false otherwise.
     */
    bool delete_basis_set();
    /** Replace internal shared basis set pointer. */
    void set_basis_set_ptr(std::shared_ptr<std::array<std::vector<primitive>, 118>> given) { basis_set = given; };
    /** Retrieve shared basis set pointer. */
    const std::shared_ptr<std::array<std::vector<primitive>, 118>> get_basis_set_ptr() const { return basis_set; };
    //-------------------atom handling--------------------------------------------------------------
    /** Cartesian coordinate value of atom nr along axis (0..2). */
    const double get_atom_coordinate(const unsigned int &nr, const unsigned int &axis) const;
    /** Get atom label (element name + index). */
    const std::string get_atom_label(const unsigned int &nr) const;
    /** Number of atoms with a basis set loaded. */
    const int get_nr_basis_set_loaded() const;
    /** Does atom nr have a basis set loaded? */
    const bool get_atom_basis_set_loaded(const int &nr) const;
    /** Primitive exponent for atom / primitive index. */
    const double get_atom_basis_set_exponent(const int &nr_atom, const int &nr_prim) const;
    /** Primitive coefficient for atom / primitive index. */
    const double get_atom_basis_set_coefficient(const int &nr_atom, const int &nr_prim) const;
    /** Change primitive exponent (marks modified). */
    bool change_atom_basis_set_exponent(const int &nr_atom, const int &nr_prim, const double &value);
    /** Change primitive coefficient (marks modified). */
    bool change_atom_basis_set_coefficient(const int &nr_atom, const int &nr_prim, const double &value);
    /** Number of primitives in atom basis set. */
    const int get_atom_primitive_count(const int &nr) const;
    /** Primitive type (angular momentum code) for atom primitive. */
    const int get_atom_primitive_type(const int &nr_atom, const int &nr_prim) const;
    /** Erase primitive from atom basis set. */
    bool erase_atom_primitive(const unsigned int &nr, const unsigned int &nr_prim);
    /** Number of shells on atom nr. */
    const int get_atom_shell_count(const unsigned int &nr) const;
    /** Number of primitives in given shell. */
    const int get_atom_shell_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    /** Type of shell (angular momentum) */
    const int get_shell_type(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    /** Center index for shell (1-based primitive center). */
    const int get_shell_center(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    /** Shell id for primitive (contraction index). */
    const int get_basis_set_shell(const unsigned int &nr_atom, const unsigned int &nr_prim) const;
    /** Start primitive index (relative to atom) for shell. */
    const int get_shell_start(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    /** Start primitive index (global primitive ordering) for shell. */
    const int get_shell_start_in_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    /** End primitive index (relative) for shell. */
    const int get_shell_end(const unsigned int &nr_atom, const unsigned int &nr_shell) const;
    /** Add new atom with coordinates + charge. */
    bool push_back_atom(const std::string &label, const double &x, const double &y, const double &z, const int &charge, const std::string& ID = "0000000000000");
    /** Add existing atom object (copy). */
    bool push_back_atom(const atom &given);
    /** Safe atom accessor (returns default atom if out of range). */
    const atom get_atom(const int& nr) const { if(nr >= 0 && nr < ncen) return atoms[nr]; else return atom(); };
    /** Add primitive to atom basis (uncontracted entry). */
    bool push_back_atom_basis_set(const int &nr, const double &exponent, const double &coefficient, const int &type, const int &shell)
    {
        if (nr <= ncen && nr >= 0)
            return atoms[nr].push_back_basis_set(exponent, coefficient, type, shell);
        else
            return false;
    };
    /** Verbose atom printout. */
    void print_atom_long(const int &nr) const { if (nr <= ncen && nr >= 0) atoms[nr].print_values_long(); };
    /** Atomic number (charge) accessor. */
    const int get_atom_charge(const int &nr) const;
    /** Integer isotope mass lookup. */
    const unsigned int get_atom_integer_mass(const unsigned int &atomnr) const;
    /** Real (average) isotope mass lookup. */
    const double get_atom_real_mass(const int &atomnr) const;
    /** Alternative atom label accessor (overload). */
    std::string get_atom_label(const int& nr) const;
    /** Set atom label. */
    void set_atom_label(const int& nr, const std::string& label) { atoms[nr].set_label(label); };
    /** Copy of atom vector. */
    std::vector<atom> get_atoms() const { return atoms; };
    /** Pointer to internal atom vector. */
    const std::vector<atom>* get_atoms_ptr() const { return &atoms; };
    /** Replace entire atom list (dangerous – counts not auto-rebuilt). */
    void set_atoms(const std::vector<atom> given) { atoms = given; };
    /** ECP core electrons on atom. */
    int get_atom_ECP_electrons(const int& nr) const;
    /** Erase all primitives for atom (clears basis set). */
    void clear_atom_basis_set(const int& nr) { atoms[nr].clear_basis_set(); };
    /** Basis set primitive count for atom. */
    int get_atom_basis_set_size(const int& nr) const { return atoms[nr].get_basis_set_size(); };
    /** Get full basis set entry (coefficient/exponent/type/shell). */
    basis_set_entry get_atom_basis_set_entry(const int& nr, const int& bs) const;
    /** Assign anisotropic displacement parameters (ADPs). */
    void set_atom_ADPs(const int& nr, const vec2& adps) { atoms[nr].set_ADPs(adps); };
    /** Set fractional coordinates for atom (crystallography). */
    void set_atom_frac_coords(const int& nr, const std::array<double, 3>& frac) { atoms[nr].set_frac_coords(frac); };
    int get_atom_basis_set_id(const int& nr) const { return atoms[nr].get_basis_set_id(); };
    //----------Calcualtion of Properties-----------------
    /** Density at position (helper that allocates temporaries). */
    const double compute_dens(const double &Pos1, const double &Pos2, const double &Pos3) const;
    /** Density with caller-provided scratch arrays (faster, reusable). */
    const double compute_dens(const double &Pos1, const double &Pos2, const double &Pos3, vec2 &d, vec &phi) const;
    /** Spin density at position (allocating version). */
    const double compute_spin_dens(const double &Pos1, const double &Pos2, const double &Pos3) const;
    /** Spin density with scratch arrays. */
    const double compute_spin_dens(const double &Pos1, const double &Pos2, const double &Pos3, vec2 &d, vec &phi) const;
    /** Evaluate multiple properties (density, gradient norm, Hessian, ELF/ELI/Laplacian). */
    const void computeValues(const std::array<double, 3>& PosGrid, double &Rho, double &normGrad, double *Hess, double &Elf, double &Eli, double &Lap) const;
    /** Compute Laplacian, ELI and ELF together. */
    const void computeLapELIELF(const std::array<double, 3>& PosGrid, double &Elf, double &Eli, double &Lap) const;
    /** Compute ELI and ELF only. */
    const void computeELIELF(const std::array<double, 3>& PosGrid, double &Elf, double &Eli) const;
    /** Compute Laplacian and ELI only. */
    const void computeLapELI(const std::array<double, 3>& PosGrid, double &Eli, double &Lap) const;
    /** Compute Laplacian only. */
    const double computeLap(const std::array<double, 3>& PosGrid) const;
    /** Compute ELI alone. */
    const void computeELI(const std::array<double,3>&PosGrid, double &Eli) const;
    /** Compute ELF alone. */
    const void computeELF(const std::array<double,3>&PosGrid, double &Elf) const;
    /** Molecular orbital value at position (Cartesian primitives). */
    const double computeMO(const std::array<double, 3>& PosGrid, const int &mo) const;
    /** Spherical MO value (not fully implemented). */
    const double compute_MO_spherical(const double &Pos1, const double &Pos2, const double &Pos3, const int &MO) const;
    /** Electrostatic potential including nuclear cores. */
    const double computeESP(const std::array<double, 3>& PosGrid, const vec2 &d2) const;
    /** Electrostatic potential excluding nuclear cores. */
    const double computeESP_noCore(const std::array<double, 3>& PosGrid, const vec2 &d2) const;
    //----------DM Handling--------------------------------
    /** Build density (and optionally spin density) matrix; loads basis if required. */
    bool build_DM(std::string basis_set_path, bool debug = false);
    /** Append value to upper triangular density matrix storage. */
    void push_back_DM(const double &value = 0.0);
    /** Set value in upper triangular density matrix. */
    bool set_DM(const int &nr, const double &value = 0.0);
    /** Get value from upper triangular density matrix. */
    const double get_DM(const int &nr) const;
    /** Size of upper triangular density matrix storage. */
    const int get_DM_size() const { return (int)UT_DensityMatrix.size(); };
    /** Resize density matrix storage (values initialized). */
    void resize_DM(const int &size, const double &value = 0.0);
    /** Full symmetric density matrix (mdarray form) if available. */
    dMatrix2 get_dm() const { return DM; };
    //----------S_DM Handling--------------------------------
    /** Append spin density matrix element. */
    void push_back_SDM(const double &value = 0.0);
    /** Set spin density matrix element. */
    bool set_SDM(const int &nr, const double &value = 0.0);
    /** Get spin density matrix element. */
    const double get_SDM(const int &nr) const;
    /** Size of spin density matrix container. */
    const int get_SDM_size() const { return (int)UT_SpinDensityMatrix.size(); };
    /** Resize spin density matrix storage. */
    void resize_SDM(const int &size, const double &value = 0.0);
    //-----------Cube handling-------------------------
    /** Load and register cube file (density / MO grid). */
    bool push_back_cube(const std::string &filepath, const bool &full, const bool &expert = false);
    /** Push already constructed cube object. */
    void push_back_cube(cube given) { cub.push_back(given); };
    /** Remove last cube. */
    void pop_back_cube();
    /** Return cube by value (copy). */
    const cube get_cube(const int& nr) const { return cub[nr]; };
    /** Pointer to cube (non-owning). */
    const cube* get_cube_ptr(const int& nr) { return &cub[nr]; };
    /** Number of cube objects. */
    const int get_cube_count() const { return (int)cub.size(); };
    /** Path of cube file. */
    std::filesystem::path get_cube_path(const int& nr) const;
    /** Write cube (standard format). */
    void write_cube_file(const int& nr, const std::filesystem::path& filename, const bool& debug = false);
    /** Write cube in dGrid compatible format. */
    void write_cube_dgrid(const int& nr, const std::filesystem::path& filename, const bool& debug = false);
    /** Write cube in XDGraph compatible format. */
    void write_cube_xdgraph(const int& nr, const std::filesystem::path& filename, const bool& debug = false);
    /** Check if cube data is loaded in memory. */
    bool get_cube_loaded(const int& nr) const { return cub[nr].get_loaded(); };
    /** Read cube values from disk. */
    bool read_cube(const int& nr, const bool& full, const bool& header, const bool& expert = false) { return cub[nr].read_file(full, header, expert); };
    /** Create a super-cell cube (returns new path). */
    std::filesystem::path make_super_cube(const int& nr) { return cub[nr].super_cube(); };
    /** Apply absolute threshold to cube. */
    bool apply_cube_thresh(const int& nr, const double& thresh) { return cub[nr].thresh(thresh); };
    /** Apply threshold relative to mask cube. */
    bool apply_cube_thresh(const int& nr, const cube& mask_cube, const double& thresh) { return cub[nr].thresh(mask_cube, thresh); };
    /** Multiply by mask cube (zero outside mask). */
    bool apply_cube_mask(const int& nr, const cube& mask_cube) { return cub[nr].mask(mask_cube); };
    /** Zero values where mask cube is non-zero (negative mask). */
    bool apply_cube_negative_mask(const int& nr, const cube& mask_cube) { return cub[nr].negative_mask(mask_cube); };
    /** Add cube to selected cube (in-place). */
    bool cube_add (const int& nr, const cube& right) { return cub[nr] += right; };
    /** Subtract cube from selected cube. */
    bool cube_subtract(const int& nr, const cube& right) { return cub[nr] -= right; };
    /** Multiply cube by another cube (element-wise). */
    bool cube_multiply(const int& nr, const cube& right) { return cub[nr] *= right; };
    /** Divide cube by another cube (element-wise). */
    bool cube_divide(const int& nr, const cube& right) { return cub[nr] /= right; };
    /** Writes a cubefile of the elctron density to disc */
    void write_rho_cube(const double& radius = 3., const double& increment = 0.025) const;
    /** Calculate Cubefile of electron density**/
    void calc_rho_cube(cube& cube_data) const;
    //-----------Pointer to members---------------------
    /** Raw pointer to primitive type array. */
    const int *get_ptr_types() { return &types[0]; };
    /** Raw pointer to primitive center indices. */
    const int *get_ptr_centers() { return &centers[0]; };
    /** Raw pointer to primitive exponents. */
    const double *get_ptr_exponents() { return &exponents[0]; };
};

#include "mo_class.h"
