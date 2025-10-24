#include "pch.h"
#include "convenience.h"
#include "cell.h"
#include "tsc_block.h"
#include "test_functions.h"
#include "integrator.h"
#include "properties.h"
#include "wfn_class.h"
#include "atoms.h"
#include "JKFit.h"

#ifdef _WIN32
#include <windows.h>
#endif

std::string help_message =
    ("\n----------------------------------------------------------------------------\n"
     "          These commands and arguments are known by NoSpherA2:\n"
     "----------------------------------------------------------------------------\n\n"
     ":::::::::::::::::::::: Defaults are highlighted by [] ::::::::::::::::::::::\n\n"
     "   -wfn            <FILENAME>.xxx           Read the following wavefunction file.\n"
     "                                            Supported filetypes: .wfn/wfx/ffn; .molden; .xyz; .gbw; .xtb; fch* (tested for OCC)\n"
     "   -fchk           <FILENAME>.fchk          Write a wavefunction to the given filename [requires -b and -d]\n"
     "   -b              <FILENAME>               Read this basis set\n"
     "   -d              <PATH>                   Path to basis_sets directory with basis_sets in tonto style\n"
     "   -dmin           <NUMBER>                 Minimum d-spacing to consider for scattering factors (repalaces hkl file)\n"
     "   -hkl_min_max    <6 Numbers>              Performs calculation on hkl range defined by the 6 numbers. (replaces dmin and hkl)\n"
     "   -ECP            <NUMBER>                 Defines the ECP corrections to be applied to a wavefunction. The given Number defines the ECP type:\n"
     "                                            [1]: def2-ECP\n"
     "                                            [2]: xTB\n"
     "                                            [3]: pTB\n"
     "   --help/-help/--h                         print this help\n"
     "   -cpus           <NUMBER>                 Sets the number of available threads to use. Defaults to all available CPUs"
     "   -v                                       Turn on Verbose (debug) Mode (Slow and a LOT of output!)\n"
     "   -v2                                      Even more stuff\n"
     "   -mult           <NUMBER>                 Input multiplicity of wavefunction (otherwise attempted to be read from the wfn)\n"
     "   -charge         <NUMBER>                 Input charge of wavefunction (otherwise attempted to be read from the wfn)\n"
     "   -method         <METHOD NAME>            Can be [RKS] or RHF to distinguish between DFT and HF\n"
     "   -cif            <FILENAME>.cif           CIF to get labels of atoms to use for calculation of scatteriung factors\n"
     "   -IAM                                     Make scattering factors based on Thakkar functions for atoms in CIF\n"
     "   -xyz            <FILENAME>.xyz           Read atom positions from this xyz file for IAM\n"
     "   -hkl            <FILENAME>.hkl           hkl file (ideally merged) to use for calculation of form factors. Use is discouraged!\n"
     "   -group          <LIST OF INT NUMBERS>    Disorder groups to be read from the CIF for consideration as asym unit atoms (space separated).\n"
     "   -acc            0,1,[2],3,4...           Accuracy of numerical grids used, where the number indicates a pre-defined level. 4 should be considered maximum,\n"
     "                                            anything above will most likely introduce numberical error and is just implemented for testing purposes.\n"
     "   -gbw2wfn                                 Only reads wavefucntion from .gbw specified by -wfn and prints it into .wfn format.\n"
     "   -TFVC                                    Use the Topological Fuzzy Voronoi Cells (TFVC) partitioning scheme instead of Hirshfeld for partitioning the electron density.\n"
     "   -Becke                                   Use Becke partitioning scheme instead of Hirshfeld for partitioning the electron density.\n"
     "   -tscb           <FILENAME>.tscb          Convert binary tsc file to bigger, less accurate human-readable form.\n"
     "   -twin           -1 0 0 0 -1 0 0 0 -1     3x3 floating-point-matrix in the form -1 0 0 0 -1 0 0 0 -1 which contains the twin matrix to use.\n"
     "                                            If there is more than a single twin law to be used, use the twin command multiple times.\n"
     "   -merge          <List of .tsc files>     Names/Paths to .tsc/.tscb files to be merged.\n"
     "   -merge_nocheck  <List of .tsc files>     Names/Paths to .tsc/.tscb files to be merged. They need to have identical hkl values.\n"
     "   -mtc            <List of .wfns + parts>  Performs calculation for a list of wavefunctions (=Multi-Tsc-Calc), where asymmetric unit is.\n"
     "                                            taken from given CIF. Also disorder groups are required per file as comma separated list\n"
     "                                            without spaces.\n"
     "   -salted         <Path to Model folder>   Uses a provided SALTED-ML Model to predict the electron densitie of a xyz-file\n"
     "   -ri_fit         <Aux Basis> <BETA>       Uses RI-Fitting to partition the electron density. If Aux Basis == 'auto' a optinal beta value can be given.\n"
     "   -mtc_mult       <List of multiplicity>   Matching multiplicity for -cmtc and -mtc wavefucntions\n"
     "   -mtc_charge     <List of charges>        Matching charges for -cmtc and -mtc wavefucntions\n"
     "   -mtc_ECP        <List of ECP modes>      Matching ECP modes for -cmtc and -mtc wavefucntions\n"
     "   -QCT                                     Starts the old QCT menu and options for working on wavefunctions/cubes and calcualtions\n"
     "                                            TIP: This mode can use many parameters like -radius, -b, -d, so they do not have to be mentioned later\n"
     "   -laplacian_bonds <Path to wavefunction>  Calculates the Laplacian of the electron density along the direct line between atoms that might be bonded by distance\n"
     "   -cmtc            <List of .wfns + parts> Performs calculation for a list of wavefunctions AND CIFs (=CIF-based-multi-Tsc-Calc), where asymmetric unit is defined by each CIF that matches a wfn.\n"
     "      Normal:       NoSpherA2.exe -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7\n"
     "      thakkar-tsc:  NoSpherA2.exe -cif A.cif -hkl A.hkl -xyz A.xyz -acc 1 -cpus 7 -IAM\n"
     "      Disorder:     NoSpherA2.exe -cif A.cif -hkl A.hkl -acc 1 -cpus 7 -mtc 1.wfn 0,1 2.wfn 0,2 3.wfn 0,3 -mtc_charge 0 0 0 -mtc_mult 1 1 1 -mtc_ECP 0 0 0\n"
     "      fragHAR:      NoSpherA2.exe -cif A.cif -hkl A.hkl -acc 1 -cpus 7 -cmtc 1.wfn 1.cif 0 2.wfn 2.cif 0 3_1.wfn 3_1.cif 0,1 3_2.wfn 3_2.cif 0,2\n"
     "      merging tscs: NoSpherA2.exe -merge A.tsc B.tsc C.tsc (also works for tscb files)\n"
     "      merge tsc(2): NoSpherA2.exe -merge_nocheck A.tsc B.tsc C.tsc  (MAKE SURE THEY HAVE IDENTICAL HKL INIDCES!!)\n"
     "      convert tsc:  NoSpherA2.exe -tscb A.tscb\n"
     "      convert gbw:  NoSpherA2.exe -gbw2wfn -wfn A.gbw\n"
     "      twin law:     NoSpherA2.exe -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7 -twin -1 0 0 0 -1 0 0 0 -1\n");
std::string NoSpherA2_message(bool no_date)
{
    std::string t = "    _   __     _____       __              ___   ___\n";
    t.append("   / | / /___ / ___/____  / /_  ___  _____/   | |__ \\\n");
    t.append("  /  |/ / __ \\\\__ \\/ __ \\/ __ \\/ _ \\/ ___/ /| | __/ /\n");
    t.append(" / /|  / /_/ /__/ / /_/ / / / /  __/ /  / ___ |/ __/\n");
    t.append("/_/ |_/\\____/____/ .___/_/ /_/\\___/_/  /_/  |_/____/\n");
    t.append("                /_/\n");
    t.append("This software is part of the cuQCT software suite developed by Florian Kleemiss.\n");
    t.append("Please give credit and cite corresponding pieces!\nThis Software is published with BSD-2 clause license.\n");
    if (!no_date)
    {
        t.append("List of contributors of pieces of code or functionality:\n");
        t.append("      Florian Kleemiss,\n");
        t.append("      Emmanuel Hupf,\n");
        t.append("      Alessandro Genoni,\n");
        t.append("      Lukas M. Seifert,\n");
        t.append("      Daniel Bruex,\n");
        t.append("      Marti Gimferrer,\n");
        t.append("      and many more in communications or by feedback!\n");
        t.append("NoSpherA2 uses featomic, Metatensor, and the mdspan library.\n");
        t.append("The used packages are published under BSD-3 clause License.\n");
        t.append("Please see, respectively:\n");
        t.append("   https://github.com/Luthaf/featomic\n");
        t.append("   https://github.com/lab-cosmo/metatensor\n");
        t.append("   This software utilizes Intel(c) Math Kernel Library (oneMKL), version 2025.2.0.629, for optimized mathematical computations\n");
        t.append("NoSpherA2 was published at  : Kleemiss et al. Chem. Sci., 2021, 12, 1675 - 1692.\n");
        t.append("Slater IAM was published at : Kleemiss et al. J. Appl. Cryst. 2024, 57, 161 - 174.\n");
        t.append("ECP correction functions at : Kleemiss et al. J. Appl. Cryst. 2025, 58, 374 - 382.\n");
        t.append("TFVC partitioning at        : Gimferrer et al. TBA.\n");
    }
    return t;
}

std::string build_date = ("This Executable was built on: " + std::string(__DATE__) + " " + std::string(__TIME__) + "\n");

bool is_similar_rel(const double &first, const double &second, const double &tolerance)
{
    double diff = abs(first - second);
    if (diff > abs((first + second + 0.01) * tolerance / 2))
        return false;
    else
        return true;
};

bool is_similar(const double &first, const double &second, const double &tolerance)
{
    double diff = abs(first - second);
    if (diff > pow(10, tolerance))
        return false;
    else
        return true;
};

bool is_similar_abs(const double &first, const double &second, const double &tolerance)
{
    double diff = abs(first - second);
    if (diff > abs(tolerance))
        return false;
    else
        return true;
};

/*
cosinus_annaeherung::cosinus_annaeherung() : mSize(0), mBase_values(nullptr), mStepwidth(1.0) {
    resize(100);
}

void cosinus_annaeherung::resize(size_t size)
{
    mSize = size;
    if (mBase_values) delete[] mBase_values;
    mBase_values = new double[mSize + 1];
#pragma omp parallel for
    for (auto i = 0; i < mSize + 1; i++)  // Fuer einen Werte mehr die Stueststellen speichern
    {
        double y = cos((MPI2 * i) / mSize);
        //std::cout << "resize: i="<<i<<" y=" << y << endl;
        mBase_values[i] = y;
    }
    mStepwidth = MPI2 / size;
}

double cosinus_annaeherung::calculate_error_at(double x) const
{
    return cos(x) - get(x);
}
*/
void copy_file(std::filesystem::path &from, std::filesystem::path &to)
{
    std::ifstream source(from.c_str(), std::ios::binary);
    std::ofstream dest(to.c_str(), std::ios::binary);

    dest << source.rdbuf();

    source.close();
    dest.close();
};

std::array<double, 3> cross(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]};
}

double a_dot(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};

//---------------------------Configuration files ---------------------------------------------------

std::filesystem::path get_home_path(void)
{
#ifdef _WIN32
    char* homeDrive = nullptr;
    size_t len = 0;
    std::string temp1, temp2;

    errno_t err = _dupenv_s(&homeDrive, &len, "HOMEDRIVE");
    if (err == 0 && homeDrive != nullptr) {
        temp1 = homeDrive;
        // Free the allocated memory
        free(homeDrive);
    }
    else {
        std::cerr << "Failed to retrieve the environment variable." << std::endl;
    }
    err = _dupenv_s(&homeDrive, &len, "HOMEPATH");
    if (err == 0 && homeDrive != nullptr) {
        temp2 = homeDrive;
        // Free the allocated memory
        free(homeDrive);
    }
    else {
        std::cerr << "Failed to retrieve the environment variable." << std::endl;
    }
    temp1.append(temp2);
    return temp1;
#else
    const char* home_env = getenv("HOME");
    if (home_env == nullptr) {
        std::cerr << "Warning: HOME environment variable not set." << std::endl;
        return std::filesystem::path("/tmp"); // Fallback to /tmp
    }
    
    std::string home = home_env;
    // Basic validation: check if it's a valid path and not empty
    if (home.empty() || home.find_first_of('\0') != std::string::npos) {
        std::cerr << "Warning: Invalid HOME environment variable." << std::endl;
        return std::filesystem::path("/tmp"); // Fallback to /tmp
    }
    
    return std::filesystem::path(home);
#endif
}

bool check_bohr(WFN &wave, bool debug)
{
    double min_length = 300.0;
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        double atom1[3]{0, 0, 0};
        for (int x = 0; x < 3; x++)
            atom1[x] = wave.get_atom_coordinate(i, x);
        for (int j = i + 1; j < wave.get_ncen(); j++)
        {
            double atom2[3]{0, 0, 0};
            for (int x = 0; x < 3; x++)
                atom2[x] = wave.get_atom_coordinate(j, x);
            double d[3]{0, 0, 0};
            d[0] = atom1[0] - atom2[0];
            d[1] = atom1[1] - atom2[1];
            d[2] = atom1[2] - atom2[2];
            double length = sqrt(pow(d[0], 2) + pow(d[1], 2) + pow(d[2], 2));
            if (debug)
                std::cout << "Length for: " << i << ";" << j << ": " << length << ", min_length: " << min_length << std::endl;
            if (length < min_length)
                min_length = length;
        }
    }
    if (debug)
    {
        if (min_length < 2)
            std::cout << "Decided it's written in Angstrom" << std::endl;
        else
            std::cout << "Decided it's written in Bohr" << std::endl;
    }
    return (!(min_length < 2));
};

std::string go_get_string(std::ifstream &file, std::string search, bool rewind)
{
    if (rewind)
    {
        file.clear();
        file.seekg(0, file.beg);
    }
    std::string line;
    while (line.find(search) == std::string::npos && !file.eof() && getline(file, line))
        continue;
    if (file.eof())
        return "";
    else
        return line;
}

std::string shrink_string(std::string &input)
{
    while (input.find(" ") != -1)
    {
        input.erase(input.find(" "), 1);
    }
    while (input.find("1") != -1)
    {
        input.erase(input.find("1"), 1);
    }
    while (input.find("2") != -1)
    {
        input.erase(input.find("2"), 1);
    }
    while (input.find("3") != -1)
    {
        input.erase(input.find("3"), 1);
    }
    while (input.find("4") != -1)
    {
        input.erase(input.find("4"), 1);
    }
    while (input.find("5") != -1)
    {
        input.erase(input.find("5"), 1);
    }
    while (input.find("6") != -1)
    {
        input.erase(input.find("6"), 1);
    }
    while (input.find("7") != -1)
    {
        input.erase(input.find("7"), 1);
    }
    while (input.find("8") != -1)
    {
        input.erase(input.find("8"), 1);
    }
    while (input.find("9") != -1)
    {
        input.erase(input.find("9"), 1);
    }
    while (input.find("0") != -1)
    {
        input.erase(input.find("0"), 1);
    }
    while (input.find("(") != -1)
    {
        input.erase(input.find("("), 1);
    }
    while (input.find(")") != -1)
    {
        input.erase(input.find(")"), 1);
    }
    return input;
};

std::string shrink_string_to_atom(std::string &input, const int &atom_number)
{
    while (input.find(" ") != -1)
    {
        input.erase(input.find(" "), 1);
    }
    while (input.find("1") != -1)
    {
        input.erase(input.find("1"), 1);
    }
    while (input.find("2") != -1)
    {
        input.erase(input.find("2"), 1);
    }
    while (input.find("3") != -1)
    {
        input.erase(input.find("3"), 1);
    }
    while (input.find("4") != -1)
    {
        input.erase(input.find("4"), 1);
    }
    while (input.find("5") != -1)
    {
        input.erase(input.find("5"), 1);
    }
    while (input.find("6") != -1)
    {
        input.erase(input.find("6"), 1);
    }
    while (input.find("7") != -1)
    {
        input.erase(input.find("7"), 1);
    }
    while (input.find("8") != -1)
    {
        input.erase(input.find("8"), 1);
    }
    while (input.find("9") != -1)
    {
        input.erase(input.find("9"), 1);
    }
    while (input.find("0") != -1)
    {
        input.erase(input.find("0"), 1);
    }
    while (input.find("(") != -1)
    {
        input.erase(input.find("("), 1);
    }
    while (input.find(")") != -1)
    {
        input.erase(input.find(")"), 1);
    }
    std::string temp = constants::atnr2letter(atom_number);
    err_checkf(temp != "PROBLEM", "Problem identifying atoms!", std::cout);
    if (input.find(temp) != 1)
        return temp;
    if (temp != "PROBLEM")
        while (input.size() > temp.size())
            input.pop_back();
    return input;
};

bool read_block_from_fortran_binary(std::ifstream &file, void *Target)
{
    int size_begin = 0, size_end = 0;
    file.read(reinterpret_cast<char *>(&size_begin), sizeof(int));
    file.read(reinterpret_cast<char *>(Target), size_begin);
    file.read(reinterpret_cast<char *>(&size_end), sizeof(int));
    if (size_begin != size_end)
    {
        std::cout << "Error reading block from binary file: " << size_begin << " vs. " << size_end << std::endl;
        return false;
    }
    return true;
}

primitive::primitive(int c, int t, double e, double coef) : center(c), type(t), exp(e), coefficient(coef)
{
    norm_const = pow(
        pow(2, 7 + 4 * type) * pow(exp, 3 + 2 * type) / constants::PI / pow(doublefactorial(2 * type + 1), 2),
        0.25);
};

primitive::primitive(const SimplePrimitive& other) : center(other.center), type(other.type), exp(other.exp), coefficient(other.coefficient)
{
    norm_const = pow(
        pow(2, 7 + 4 * type) * pow(exp, 3 + 2 * type) / constants::PI / pow(doublefactorial(2 * type + 1), 2),
        0.25);
};

void select_cubes(std::vector<std::vector<unsigned int>> &selection, std::vector<WFN> &wavy, unsigned int nr_of_cubes, bool wfnonly, bool debug)
{
    // asks which wfn to use, if wfnonly is set or whcih cubes up to nr of cubes to use
    // Returns values in selection[0][i] for iths selection of wavefunction and
    //  selection[1][i] for iths selection of cube
    using namespace std;
   std::cout << "Which of the following cubes to use? Need to select " << nr_of_cubes << " file";
    if (nr_of_cubes > 1)
       std::cout << "s in total." << endl;
    else
       std::cout << "." << endl;
   std::cout << endl
         << endl;
    for (int w = 0; w < wavy.size(); w++)
    {
        stringstream stream;
       std::cout << "_____________________________________________________________" << endl;
       std::cout << "WFN ";
        stream << setw(2) << w;
       std::cout << stream.str() << ") " << wavy[w].get_path().stem() << endl;
        stream.str("");
        for (int c = 0; c < wavy[w].get_cube_count(); c++)
        {
            if (c == 0)
               std::cout << "        |" << endl
                     << "Cube    |" << endl;
            else
               std::cout << "        |" << endl;
            if (!wfnonly)
            {
               std::cout << setw(2) << w;
               std::cout << ".";
               std::cout << setw(2) << c;
            }
            else
               std::cout << "     ";
           std::cout << "   |_ " << wavy[w].get_cube_path(c).stem();
            if (!exists(wavy[w].get_cube_path(c)))
               std::cout << " (MEM ONLY)";
           std::cout << endl;
        }
       std::cout << "_____________________________________________________________" << endl
             << endl
             << endl;
    }
    // bool happy = false;
    unsigned int selected_cubes = 0;
    do
    {
       std::cout << "Select " << selected_cubes + 1 << ". ";
        if (wfnonly)
           std::cout << "WFN ";
        else
           std::cout << "cube ";
       std::cout << "please: ";
        string input;
        cin >> input;
        if (!wfnonly)
        {
            if (input.find('.') == string::npos)
            {
               std::cout << "no . found in input!" << endl;
                continue;
            }
        }
        else
        {
            if (input.find('.') == string::npos)
               std::cout << "Ignoring the .!" << endl;
            unsigned int nr_wave = fromString<unsigned int>(input);
            if (nr_wave < 0 || nr_wave >= wavy.size())
            {
               std::cout << "Invalid choice!" << endl;
                continue;
            }
            selected_cubes++;
            selection[0].push_back(nr_wave);
            if (selected_cubes == nr_of_cubes)
                return;
            else
                continue;
        }
        if (debug)
        {
           std::cout << "input: " << input << endl;
           std::cout << "with . found at: " << input.find('.') << endl;
           std::cout << "substr1: " << input.substr(0, input.find('.')) << endl;
           std::cout << "substr2: " << input.substr(input.find('.') + 1) << endl;
        }
        string wave(input.substr(0, input.find('.')));
        string cube(input.substr(input.find('.') + 1));
        unsigned int nr_wave = fromString<unsigned int>(wave);
        int nr_cube = fromString<int>(cube);
        if (debug)
           std::cout << "Translated: " << nr_wave << " " << nr_cube << endl;
        if (nr_wave < 0 || nr_wave >= wavy.size() || nr_cube < 0 || nr_cube >= wavy[nr_wave].get_cube_count())
        {
           std::cout << "Invalid choice!" << endl;
            continue;
        }
        selection[0][selected_cubes] = nr_wave;
        selection[1][selected_cubes] = nr_cube;
        selected_cubes++;
        if (selected_cubes == nr_of_cubes)
        {
            if (debug)
               std::cout << "Going to return!" << endl;
            return;
        }
    } while (true);
};

bool unsaved_files(std::vector<WFN> &wavy)
{
    for (int w = 0; w < wavy.size(); w++)
        for (int c = 0; c < wavy[w].get_cube_count(); c++)
            if (!exists(wavy[w].get_cube_path(c)))
                return true;
    return false;
};

void readxyzMinMax_fromWFN(
    WFN &wavy,
    double *CoordMinMax,
    int *NbSteps,
    double Radius,
    double Increments,
    bool no_bohr)
{
    vec2 PosAtoms;
    PosAtoms.resize(3);
    for (int i = 0; i < 3; i++)
        PosAtoms[i].resize(wavy.get_ncen());
    bool bohrang = true;
    if (!no_bohr)
        bohrang = !check_bohr(wavy, false);

    for (int j = 0; j < wavy.get_ncen(); j++)
    {
        PosAtoms[0][j] = wavy.get_atom_coordinate(j, 0);
        PosAtoms[1][j] = wavy.get_atom_coordinate(j, 1);
        PosAtoms[2][j] = wavy.get_atom_coordinate(j, 2);
        if (!bohrang)
        {
            for (int i = 0; i < 3; i++)
                PosAtoms[i][j] = constants::ang2bohr(PosAtoms[i][j]);
        }
    }
    CoordMinMax[0] = *std::min_element(PosAtoms[0].begin(), PosAtoms[0].end());
    CoordMinMax[3] = *std::max_element(PosAtoms[0].begin(), PosAtoms[0].end());
    CoordMinMax[1] = *std::min_element(PosAtoms[1].begin(), PosAtoms[1].end());
    CoordMinMax[4] = *std::max_element(PosAtoms[1].begin(), PosAtoms[1].end());
    CoordMinMax[2] = *std::min_element(PosAtoms[2].begin(), PosAtoms[2].end());
    CoordMinMax[5] = *std::max_element(PosAtoms[2].begin(), PosAtoms[2].end());

    const double temp_rad = constants::ang2bohr(Radius);
    CoordMinMax[0] -= temp_rad;
    CoordMinMax[3] += temp_rad;
    CoordMinMax[1] -= temp_rad;
    CoordMinMax[4] += temp_rad;
    CoordMinMax[2] -= temp_rad;
    CoordMinMax[5] += temp_rad;

    NbSteps[0] = (int)ceil(constants::bohr2ang(CoordMinMax[3] - CoordMinMax[0]) / Increments);
    NbSteps[1] = (int)ceil(constants::bohr2ang(CoordMinMax[4] - CoordMinMax[1]) / Increments);
    NbSteps[2] = (int)ceil(constants::bohr2ang(CoordMinMax[5] - CoordMinMax[2]) / Increments);
}

void readxyzMinMax_fromCIF(
    std::filesystem::path cif,
    double *CoordMinMax,
    int *NbSteps,
    vec2 &cm,
    double Resolution)
{
    using namespace std;
    cell cell(cif);

    cm[0][0] = constants::ang2bohr(cell.get_a());
    cm[0][1] = constants::ang2bohr(cell.get_b() * cell.get_cg());
    cm[0][2] = constants::ang2bohr(cell.get_c() * cell.get_cb());
    cm[1][1] = constants::ang2bohr(cell.get_b() * cell.get_sg());
    cm[1][2] = constants::ang2bohr(cell.get_c() * (cell.get_ca() - cell.get_cb() * cell.get_cg()) / cell.get_sg());
    cm[2][2] = constants::ang2bohr(cell.get_V() / (cell.get_a() * cell.get_b() * cell.get_sg()));

    CoordMinMax[0] = 0.0;
    CoordMinMax[1] = 0.0;
    CoordMinMax[2] = 0.0;

    CoordMinMax[3] = (cell.get_a() + cell.get_b() * cell.get_cg() + cell.get_c() * cell.get_cb()) / 0.529177249;
    CoordMinMax[4] = (cell.get_b() * cell.get_sg() + cell.get_c() * (cell.get_ca() - cell.get_cb() * cell.get_cg()) / cell.get_sg()) / 0.529177249;
    CoordMinMax[5] = cm[2][2];

    NbSteps[0] = (int)ceil(cell.get_a() / Resolution);
    NbSteps[1] = (int)ceil(cell.get_b() / Resolution);
    NbSteps[2] = (int)ceil(cell.get_c() / Resolution);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cm[i][j] /= NbSteps[j];
}

bool generate_sph2cart_mat(vec2 &p, vec2 &d, vec2 &f, vec2 &g)
{
    //
    // From 3P: P0 P1 P2
    // To 3P : Z X Y (4 2 3, as in ORCA format)
    //
    p.resize(3);
#pragma omp parallel for
    for (int i = 0; i < 3; i++)
    {
        p[i].resize(3, 0.0);
    }
    p[0][2] = 1.0;
    p[1][0] = 1.0;
    p[2][1] = 1.0;

    //
    // From 5D: D 0, D + 1, D - 1, D + 2, D - 2
    // To 6D : 5  6  7  8  9 10
    // XX, YY, ZZ, XY, XZ, YZ
    //
    d.resize(6);
#pragma omp parallel for
    for (int i = 0; i < 6; i++)
    {
        d[i].resize(5, 0.0);
    }
    // XX = -0.5/SQRT(3) * D0 + 0.5 * D2
    d[0][0] = -0.5 / sqrt(3);
    d[0][3] = 0.5;
    // YY = -0.5/SQRT(3) * D0 - 0.5 * D2
    d[1][0] = -0.5 / sqrt(3);
    d[1][3] = -0.5;
    // ZZ = SQRT(1/3) * D0
    d[2][0] = sqrt(1.0 / 3.0);
    // XY = D-2
    d[3][4] = 1.0;
    // XZ = D1
    d[4][1] = 1.0;
    // YZ = D-1
    d[5][2] = 1.0;

    // From 7F: F 0, F + 1, F - 1, F + 2, F - 2, F + 3, F - 3
    // To 10F : 11   12   13   14   15   16   17   18   19  20
    // XXX, YYY, ZZZ, XXY, XXZ, YYZ, XYY, XZZ, YZZ, XYZ (AIMALL order!)
    //
    f.resize(10);
#pragma omp parallel for
    for (int i = 0; i < 10; i++)
    {
        f[i].resize(7, 0.0);
    }
    f[0][1] = -sqrt(0.025);
    f[0][5] = -sqrt(1.0 / 24.0);

    f[1][2] = -sqrt(0.025);
    f[1][6] = sqrt(1.0 / 24.0);

    f[2][0] = sqrt(1.0 / 15.0);

    f[3][2] = -sqrt(0.025);
    f[3][6] = -sqrt(0.375);

    f[4][0] = -sqrt(0.15);
    f[4][3] = 0.5;

    f[5][0] = -sqrt(0.15);
    f[5][3] = -0.5;

    f[6][1] = -sqrt(0.025);
    f[6][5] = sqrt(0.375);

    f[7][1] = sqrt(0.4);

    f[8][2] = sqrt(0.4);

    f[9][4] = 1.0;

    g.resize(15);
#pragma omp parallel for
    for (int i = 0; i < 15; i++)
    {
        g[i].resize(9, 0.0);
    }
    g[0][0] = 0.375 * sqrt(1.0 / 35.0);
    g[0][3] = -0.25 / sqrt(7);
    g[0][7] = -0.125;

    g[1][0] = g[0][0];
    g[1][3] = -g[0][3];
    g[1][7] = g[0][7];

    g[2][0] = sqrt(1.0 / 35.0);

    g[3][4] = -sqrt(1.0 / 28.0);
    g[3][8] = -0.5;

    g[4][1] = -0.5 / sqrt(98.0 / 63.0);
    g[4][5] = -1.0 / sqrt(8.0);

    g[5][4] = g[3][4];
    g[5][8] = -g[3][8];

    g[6][2] = g[4][1];
    g[6][6] = -g[4][5];

    g[7][1] = sqrt(2.0 / 7.0);

    g[8][2] = g[7][1];

    g[9][0] = 3.0 / sqrt(560);
    g[9][7] = 0.75;

    g[10][0] = -3.0 / sqrt(35);
    g[10][3] = 1.5 / sqrt(7);

    g[11][0] = g[10][0];
    g[11][3] = -g[10][3];

    g[12][2] = g[4][1];
    g[12][6] = -0.75 * sqrt(2);

    g[13][1] = g[4][1];
    g[13][5] = -g[12][6];

    g[14][4] = 3.0 / sqrt(7);
    return true;
}
bool generate_cart2sph_mat(vec2 &d, vec2 &f, vec2 &g, vec2 &h)
{
    //
    // From 5D: D 0, D + 1, D - 1, D + 2, D - 2
    // To 6D : 1  2  3  4  5  6
    // XX, YY, ZZ, XY, XZ, YZ
    //
    d.resize(6);
#pragma omp parallel for
    for (int i = 0; i < 6; i++)
    {
        d[i].resize(5, 0.0);
    }
    // D0 = -0.5 * XX - 0.5 * YY + ZZ
    d[0][0] = -0.5;
    d[1][0] = -0.5;
    d[2][0] = 1.0;
    // D + 1 = XZ
    d[4][1] = 1.0;
    // D - 1 = YZ
    d[5][2] = 1.0;
    // D + 2 = SQRT(3) / 2 * (XX - YY)
    d[0][3] = sqrt(3.0) / 2.0;
    d[1][3] = -sqrt(3.0) / 2.0;
    // D - 2 = XY
    d[3][4] = 1.0;

    // From 7F: F 0, F + 1, F - 1, F + 2, F - 2, F + 3, F - 3
    // To 10F : 1   2   3   4   5   6   7   8   9  10
    // XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ(Gaussian sequence, not identical to Multiwfn)
    //
    f.resize(10);
#pragma omp parallel for
    for (int i = 0; i < 10; i++)
    {
        f[i].resize(7, 0.0);
    }
    // F 0 = -3 / (2 * sqrt5) * (XXZ + YYZ) + ZZZ
    f[2][0] = 1.0;
    f[5][0] = -1.5 / sqrt(5.0);
    f[8][0] = -1.5 / sqrt(5.0);
    // F + 1 = -sqrt(3 / 8) * XXX - sqrt(3 / 40) * XYY + sqrt(6 / 5) * XZZ
    f[0][1] = -sqrt(3.0 / 8.0);
    f[3][1] = -sqrt(3.0 / 40.0);
    f[6][1] = sqrt(6.0 / 5.0);
    // F - 1 = -sqrt(3 / 40) * XXY - sqrt(3 / 8) * YYY + sqrt(6 / 5) * YZZ
    f[1][2] = -sqrt(3.0 / 8.0);
    f[4][2] = -sqrt(3.0 / 40.0);
    f[7][2] = sqrt(6.0 / 5.0);
    // F + 2 = sqrt3 / 2 * (XXZ - YYZ)
    f[5][3] = sqrt(3.0) / 2.0;
    f[8][3] = -sqrt(3.0) / 2.0;
    // F - 2 = XYZ
    f[9][4] = 1.0;
    // F + 3 = sqrt(5 / 8) * XXX - 3 / sqrt8 * XYY
    f[0][5] = sqrt(5.0 / 8.0);
    f[3][5] = -3.0 / sqrt(8.0);
    // F - 3 = 3 / sqrt8 * XXY - sqrt(5 / 8) * YYY
    f[1][6] = -sqrt(5.0 / 8.0);
    f[4][6] = 3.0 / sqrt(8.0);

    // From 9G: G 0, G + 1, G - 1, G + 2, G - 2, G + 3, G - 3, G + 4, G - 4
    // To 15G : 1    2    3    4    5    6    7    8
    // ZZZZ, YZZZ, YYZZ, YYYZ, YYYY, XZZZ, XYZZ, XYYZ
    // 9   10   11   12   13   14   15
    // XYYY, XXZZ, XXYZ, XXYY, XXXZ, XXXY, XXXX
    //
    g.resize(15);
#pragma omp parallel for
    for (int i = 0; i < 15; i++)
    {
        g[i].resize(9, 0.0);
    }
    // G 0 = ZZZZ + 3 / 8 * (XXXX + YYYY) - 3 * sqrt(3 / 35) * (XXZZ + YYZZ - 1 / 4 * XXYY)
    g[0][0] = 1.0;
    g[2][0] = -3.0 * sqrt(3.0 / 35.0);
    g[4][0] = 3.0 / 8.0;
    g[9][0] = -3.0 * sqrt(3.0 / 35.0);
    g[11][0] = 3.0 / 4.0 * sqrt(3.0 / 35.0);
    g[14][0] = 3.0 / 8.0;
    // G + 1 = 2 * sqrt(5 / 14) * XZZZ - 3 / 2 * sqrt(5 / 14) * XXXZ - 3 / 2 / sqrt14 * XYYZ
    g[5][1] = 2.0 * sqrt(5.0 / 14.0);
    g[7][1] = -1.5 / sqrt(14.0);
    g[12][1] = -1.5 * sqrt(5.0 / 14.0);
    // G - 1 = 2 * sqrt(5 / 14) * YZZZ - 3 / 2 * sqrt(5 / 14) * YYYZ - 3 / 2 / sqrt14 * XXYZ
    g[1][2] = 2.0 * sqrt(5.0 / 14.0);
    g[3][2] = -1.5 * sqrt(5.0 / 14.0);
    g[10][2] = -1.5 / sqrt(14.0);
    // G + 2 = 3 * sqrt(3 / 28) * (XXZZ - YYZZ) - sqrt5 / 4 * (XXXX - YYYY)
    g[2][3] = -3.0 * sqrt(3.0 / 28.0);
    g[4][3] = sqrt(5.0) / 4.0;
    g[9][3] = 3.0 * sqrt(3.0 / 28.0);
    g[14][3] = -sqrt(5.0) / 4.0;
    // G - 2 = 3 / sqrt7 * XYZZ - sqrt(5 / 28) * (XXXY + XYYY)
    g[6][4] = 3.0 / sqrt(7.0);
    g[8][4] = -sqrt(5.0 / 28.0);
    g[13][4] = -sqrt(5.0 / 28.0);
    // G + 3 = sqrt(5 / 8) * XXXZ - 3 / sqrt8 * XYYZ
    g[7][5] = -3.0 / sqrt(8.0);
    g[12][5] = sqrt(5.0 / 8.0);
    // G - 3 = -sqrt(5 / 8) * YYYZ + 3 / sqrt8 * XXYZ
    g[3][6] = -sqrt(5.0 / 8.0);
    g[10][6] = 3.0 / sqrt(8.0);
    // G + 4 = sqrt35 / 8 * (XXXX + YYYY) - 3 / 4 * sqrt3 * XXYY
    g[4][7] = sqrt(35.0) / 8.0;
    g[11][7] = -3.0 / 4.0 * sqrt(3.0);
    g[14][7] = sqrt(35.0) / 8.0;
    // G - 4 = sqrt5 / 2 * (XXXY - XYYY)
    g[8][8] = -sqrt(5.0) / 2.0;
    g[13][8] = sqrt(5.0) / 2.0;

    // From 11H: H 0, H + 1, H - 1, H + 2, H - 2, H + 3, H - 3, H + 4, H - 4, H + 5, H - 5
    // To 21H : 1     2     3     4     5     6     7     8     9    10
    // ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ
    // 11    12    13    14    15    16    17    18    19    20    21
    // XYYYY XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
    //
    h.resize(21);
#pragma omp parallel for
    for (int i = 0; i < 21; i++)
    {
        h[i].resize(11);
        std::fill(h[i].begin(), h[i].end(), 0.0);
    }
    // H 0 = ZZZZZ - 5 / sqrt21 * (XXZZZ + YYZZZ) + 5 / 8 * (XXXXZ + YYYYZ) + sqrt(15 / 7) / 4 * XXYYZ
    h[0][0] = 1.0;
    h[11][0] = -5.0 / sqrt(21.0);
    h[2][0] = -5.0 / sqrt(21.0);
    h[18][0] = 5.0 / 8.0;
    h[4][0] = 5.0 / 8.0;
    h[13][0] = sqrt(15.0 / 7.0) / 4.0;
    // H + 1 = sqrt(5 / 3) * XZZZZ - 3 * sqrt(5 / 28) * XXXZZ - 3 / sqrt28 * XYYZZ + sqrt15 / 8 * XXXXX + sqrt(5 / 3) / 8 * XYYYY + sqrt(5 / 7) / 4 * XXXYY
    h[6][1] = sqrt(5.0 / 3.0);
    h[15][1] = -3.0 * sqrt(5.0 / 28.0);
    h[8][1] = -3.0 / sqrt(28.0);
    h[20][1] = sqrt(15.0) / 8.0;
    h[10][1] = sqrt(5.0 / 3.0) / 8.0;
    h[17][1] = sqrt(5.0 / 7.0) / 4.0;
    // H - 1 = sqrt(5 / 3) * YZZZZ - 3 * sqrt(5 / 28) * YYYZZ - 3 / sqrt28 * XXYZZ + sqrt15 / 8 * YYYYY + sqrt(5 / 3) / 8 * XXXXY + sqrt(5 / 7) / 4 * XXYYY
    h[1][2] = sqrt(5.0 / 3.0);
    h[3][2] = -3.0 * sqrt(5.0 / 28.0);
    h[12][2] = -3.0 / sqrt(28.0);
    h[5][2] = sqrt(15.0) / 8.0;
    h[19][2] = sqrt(5.0 / 3.0) / 8.0;
    h[14][2] = sqrt(5.0 / 7.0) / 4.0;
    // H + 2 = sqrt5 / 2 * (XXZZZ - YYZZZ) - sqrt(35 / 3) / 4 * (XXXXZ - YYYYZ)
    h[11][3] = sqrt(5.0) / 2.0;
    h[2][3] = -sqrt(5.0) / 2.0;
    h[18][3] = -sqrt(35.0 / 3.0) / 4.0;
    h[4][3] = sqrt(35.0 / 3.0) / 4.0;
    // H - 2 = sqrt(5 / 3) * XYZZZ - sqrt(5 / 12) * (XXXYZ + XYYYZ)
    h[7][4] = sqrt(5.0 / 3.0);
    h[16][4] = -sqrt(5.0 / 12.0);
    h[9][4] = -sqrt(5.0 / 12.0);
    // H + 3 = sqrt(5 / 6) * XXXZZ - sqrt(3 / 2) * XYYZZ - sqrt(35 / 2) / 8 * (XXXXX - XYYYY) + sqrt(5 / 6) / 4 * XXXYY
    h[15][5] = sqrt(5.0 / 6.0);
    h[8][5] = -sqrt(1.5);
    h[20][5] = -sqrt(17.5) / 8.0;
    h[10][5] = sqrt(17.5) / 8.0;
    h[17][5] = sqrt(5.0 / 6.0) / 4.0;
    // H - 3 = -sqrt(5 / 6) * YYYZZ + sqrt(3 / 2) * XXYZZ - sqrt(35 / 2) / 8 * (XXXXY - YYYYY) - sqrt(5 / 6) / 4 * XXYYY
    h[3][6] = -sqrt(5.0 / 6.0);
    h[12][6] = sqrt(1.5);
    h[19][6] = -sqrt(17.5) / 8.0;
    h[5][6] = sqrt(17.5) / 8.0;
    h[14][6] = -sqrt(5.0 / 6.0) / 4.0;
    // H + 4 = sqrt35 / 8 * (XXXXZ + YYYYZ) - 3 / 4 * sqrt3 * XXYYZ
    h[18][7] = sqrt(35.0) / 8.0;
    h[4][7] = sqrt(35.0) / 8.0;
    h[13][7] = -0.75 * sqrt(3.0);
    // H - 4 = sqrt5 / 2 * (XXXYZ - XYYYZ)
    h[16][8] = sqrt(5.0) / 2.0;
    h[9][8] = -sqrt(5.0) / 2.0;
    // H + 5 = 3 / 8 * sqrt(7 / 2) * XXXXX + 5 / 8 * sqrt(7 / 2) * XYYYY - 5 / 4 * sqrt(3 / 2) * XXXYY
    h[20][9] = 3.0 / 8.0 * sqrt(3.5);
    h[10][9] = 5.0 / 8.0 * sqrt(3.5);
    h[17][9] = -1.25 * sqrt(1.5);
    // H - 5 = 3 / 8 * sqrt(7 / 2) * YYYYY + 5 / 8 * sqrt(7 / 2) * XXXXY - 5 / 4 * sqrt(3 / 2) * XXYYY
    h[5][10] = 3.0 / 8.0 * sqrt(3.5);
    h[19][10] = 5.0 / 8.0 * sqrt(3.5);
    h[14][10] = -1.25 * sqrt(1.5);
    return true;
}

bool read_fracs_ADPs_from_CIF(std::filesystem::path &cif, WFN &wavy, cell &unit_cell, std::ofstream &log3, bool debug)
{
    using namespace std;
    vec2 Uij, Cijk, Dijkl;
    ifstream asym_cif_input(cif, std::ios::in);
    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    string line;
    svec labels;
    int count_fields = 0;
    int position_field[3] = {0, 0, 0};
    int label_field = 100;
    vec2 positions;
    positions.resize(wavy.get_ncen());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        positions[i].resize(3);
    }
    bool atoms_read = false;
    while (!asym_cif_input.eof() && !atoms_read)
    {
        getline(asym_cif_input, line);
        if (line.find("loop_") != string::npos)
        {
            while (line.find("_") != string::npos)
            {
                getline(asym_cif_input, line);
                if (debug)
                    log3 << "line in loop field definition: " << line << endl;
                if (line.find("label") != string::npos)
                    label_field = count_fields;
                else if (line.find("fract_x") != string::npos)
                    position_field[0] = count_fields;
                else if (line.find("fract_y") != string::npos)
                    position_field[1] = count_fields;
                else if (line.find("fract_z") != string::npos)
                    position_field[2] = count_fields;
                else if (label_field == 100)
                {
                    if (debug)
                        log3 << "I don't think this is the atom block.. moving on!" << endl;
                    break;
                }
                count_fields++;
            }
            while (line.find("_") == string::npos && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                svec fields;
                fields.resize(count_fields);
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                if (debug)
                    log3 << "label: " << fields[label_field] << " frac_position: " << stod(fields[position_field[0]]) << " " << stod(fields[position_field[1]]) << " " << stod(fields[position_field[2]]) << endl;
                positions[labels.size()] = unit_cell.get_coords_cartesian(stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]));
                bool found_this_one = false;
                if (debug)
                    log3 << "label: " << fields[label_field] << " cartesian position: " << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
                for (int i = 0; i < wavy.get_ncen(); i++)
                {
                    if (is_similar(positions[labels.size()][0], wavy.get_atom_coordinate(i, 0), -1) && is_similar(positions[labels.size()][1], wavy.get_atom_coordinate(i, 1), -1) && is_similar(positions[labels.size()][2], wavy.get_atom_coordinate(i, 2), -1))
                    {
                        if (debug)
                            log3 << "WFN position: " << wavy.get_atom_coordinate(i, 0) << " " << wavy.get_atom_coordinate(i, 1) << " " << wavy.get_atom_coordinate(i, 2) << endl
                                 << "Found an atom: " << fields[label_field] << " Corresponding to atom charge " << wavy.get_atom_charge(i) << endl;
                        wavy.set_atom_label(i, fields[label_field]);
                        wavy.set_atom_frac_coords(i, {stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]])});
                        found_this_one = true;
                        break;
                    }
                }
                if (!found_this_one && debug)
                    log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
                labels.push_back(fields[label_field]);
                getline(asym_cif_input, line);
            }
        }
    }

    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    count_fields = 0;
    int ADP_field[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    label_field = 100;
    atoms_read = false;
    Uij.resize(wavy.get_ncen());
    while (!asym_cif_input.eof() && !atoms_read)
    {
        getline(asym_cif_input, line);
        if (line.find("loop_") != string::npos)
        {
            while (line.find("_") != string::npos)
            {
                getline(asym_cif_input, line);
                if (debug)
                    log3 << "line in loop field definition: " << line << endl;
                if (line.find("aniso_label") != string::npos)
                    label_field = count_fields;
                else if (line.find("aniso_U_11") != string::npos)
                    ADP_field[0] = count_fields;
                else if (line.find("aniso_U_22") != string::npos)
                    ADP_field[1] = count_fields;
                else if (line.find("aniso_U_33") != string::npos)
                    ADP_field[2] = count_fields;
                else if (line.find("aniso_U_12") != string::npos)
                    ADP_field[3] = count_fields;
                else if (line.find("aniso_U_13") != string::npos)
                    ADP_field[4] = count_fields;
                else if (line.find("aniso_U_23") != string::npos)
                    ADP_field[5] = count_fields;
                else if (label_field == 100)
                {
                    if (debug)
                        log3 << "I don't think this is the Uij block.. moving on!" << endl;
                    break;
                }
                count_fields++;
            }
            while (line.find("_") == string::npos && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                svec fields;
                fields.resize(count_fields);
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                if (debug)
                    log3 << "label: " << fields[label_field] << endl;
                bool found_this_one = false;
                for (int i = 0; i < wavy.get_ncen(); i++)
                {
                    if (fields[label_field] == wavy.get_atom_label(i))
                    {
                        Uij[i].resize(6);
                        for (int j = 0; j < 6; j++)
                            Uij[i][j] = stod(fields[ADP_field[j]]);
                        found_this_one = true;
                        break;
                    }
                }
                if (!found_this_one && debug)
                    log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
                getline(asym_cif_input, line);
            }
        }
    }

    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    count_fields = 0;
    label_field = 100;
    atoms_read = false;
    Cijk.resize(wavy.get_ncen());
    while (!asym_cif_input.eof() && !atoms_read)
    {
        getline(asym_cif_input, line);
        if (line.find("loop_") != string::npos)
        {
            while (line.find("_") != string::npos)
            {
                getline(asym_cif_input, line);
                if (debug)
                    log3 << "line in loop field definition: " << line << endl;
                if (line.find("C_label") != string::npos)
                    label_field = count_fields;
                else if (line.find("C_111") != string::npos)
                    ADP_field[0] = count_fields;
                else if (line.find("C_112") != string::npos)
                    ADP_field[1] = count_fields;
                else if (line.find("C_113") != string::npos)
                    ADP_field[2] = count_fields;
                else if (line.find("C_122") != string::npos)
                    ADP_field[3] = count_fields;
                else if (line.find("C_123") != string::npos)
                    ADP_field[4] = count_fields;
                else if (line.find("C_133") != string::npos)
                    ADP_field[5] = count_fields;
                else if (line.find("C_222") != string::npos)
                    ADP_field[6] = count_fields;
                else if (line.find("C_223") != string::npos)
                    ADP_field[7] = count_fields;
                else if (line.find("C_233") != string::npos)
                    ADP_field[8] = count_fields;
                else if (line.find("C_333") != string::npos)
                    ADP_field[9] = count_fields;
                else if (label_field == 100)
                {
                    if (debug)
                        log3 << "I don't think this is the Cijk block.. moving on!" << endl;
                    break;
                }
                count_fields++;
            }
            while (line.find("_") == string::npos && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                svec fields;
                fields.resize(count_fields);
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                if (debug)
                    log3 << "label: " << fields[label_field] << endl;
                bool found_this_one = false;
                for (int i = 0; i < wavy.get_ncen(); i++)
                {
                    if (fields[label_field] == wavy.get_atom_label(i))
                    {
                        Cijk[i].resize(10);
                        for (int j = 0; j < 6; j++)
                            Cijk[i][j] = stod(fields[ADP_field[j]]);
                        found_this_one = true;
                        break;
                    }
                }
                if (!found_this_one && debug)
                    log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
                getline(asym_cif_input, line);
            }
        }
    }

    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    count_fields = 0;
    label_field = 100;
    atoms_read = false;
    Dijkl.resize(wavy.get_ncen());
    while (!asym_cif_input.eof() && !atoms_read)
    {
        getline(asym_cif_input, line);
        if (line.find("loop_") != string::npos)
        {
            while (line.find("_") != string::npos)
            {
                getline(asym_cif_input, line);
                if (debug)
                    log3 << "line in loop field definition: " << line << endl;
                if (line.find("D_label") != string::npos)
                    label_field = count_fields;
                else if (line.find("D_1111") != string::npos)
                    ADP_field[0] = count_fields;
                else if (line.find("D_1112") != string::npos)
                    ADP_field[1] = count_fields;
                else if (line.find("D_1113") != string::npos)
                    ADP_field[2] = count_fields;
                else if (line.find("D_1122") != string::npos)
                    ADP_field[3] = count_fields;
                else if (line.find("D_1123") != string::npos)
                    ADP_field[4] = count_fields;
                else if (line.find("D_1133") != string::npos)
                    ADP_field[5] = count_fields;
                else if (line.find("D_1222") != string::npos)
                    ADP_field[6] = count_fields;
                else if (line.find("D_1223") != string::npos)
                    ADP_field[7] = count_fields;
                else if (line.find("D_1233") != string::npos)
                    ADP_field[8] = count_fields;
                else if (line.find("D_1333") != string::npos)
                    ADP_field[9] = count_fields;
                else if (line.find("D_2222") != string::npos)
                    ADP_field[10] = count_fields;
                else if (line.find("D_2223") != string::npos)
                    ADP_field[11] = count_fields;
                else if (line.find("D_2233") != string::npos)
                    ADP_field[12] = count_fields;
                else if (line.find("D_2333") != string::npos)
                    ADP_field[13] = count_fields;
                else if (line.find("D_3333") != string::npos)
                    ADP_field[14] = count_fields;
                else if (label_field == 100)
                {
                    if (debug)
                        log3 << "I don't think this is the Dijk block.. moving on!" << endl;
                    break;
                }
                count_fields++;
            }
            while (line.find("_") == string::npos && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                svec fields;
                fields.resize(count_fields);
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                if (debug)
                    log3 << "label: " << fields[label_field] << endl;
                bool found_this_one = false;
                for (int i = 0; i < wavy.get_ncen(); i++)
                {
                    if (fields[label_field] == wavy.get_atom_label(i))
                    {
                        Dijkl[i].resize(15);
                        for (int j = 0; j < 6; j++)
                            Dijkl[i][j] = stod(fields[ADP_field[j]]);
                        found_this_one = true;
                        break;
                    }
                }
                if (!found_this_one && debug)
                    log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
                getline(asym_cif_input, line);
            }
        }
    }

    for (int i = 0; i < wavy.get_ncen(); i++)
        wavy.set_atom_ADPs(i, {Uij[i], Cijk[i], Dijkl[i]});

    return true;
};

void swap_sort(ivec order, cvec &v)
{
    int i = 0;
    while (i < v.size() - 1)
    {
        int new_index = 0;
        for (int j = i; j < v.size(); j++)
            if (order[j] < order[i])
                new_index++;
        if (new_index > 0)
        {
            std::complex<double> temp = v[i];
            v[i] = v[i + new_index];
            v[i + new_index] = temp;
            int temp2 = order[i];
            order[i] = order[i + new_index];
            order[i + new_index] = temp2;
        }
        else
            i++;
    }
}

void swap_sort_multi(ivec order, std::vector<ivec> &v)
{
    int i = 0;
    ivec temp;
    temp.resize(v.size());
    while (i < v.size() - 1)
    {
        int new_index = 0;
#pragma omp parallel for reduction(+ : new_index)
        for (int j = i; j < v.size(); j++)
            if (order[j] < order[i])
                new_index++;
        if (new_index > 0)
        {
#pragma omp parallel for
            for (int run = 0; run < v.size(); run++)
            {
                temp[run] = v[run][i];
                v[run][i] = v[run][i + new_index];
                v[run][i + new_index] = temp[run];
            }
            int temp2 = order[i];
            order[i] = order[i + new_index];
            order[i + new_index] = temp2;
        }
        else
            i++;
    }
}

double get_lambda_1(double *a)
{
    vec bw, zw;
    // int run = 0;
    double eig1, eig2, eig3;
    const double p1 = a[1] * a[1] + a[2] * a[2] + a[5] * a[5];
    if (p1 == 0)
    {
        eig1 = a[0];
        eig2 = a[4];
        eig3 = a[8];
        if ((eig1 < eig2 && eig2 < eig3) || (eig3 < eig2 && eig2 < eig1))
            return eig2;

        else if ((eig2 < eig1 && eig1 < eig3) || (eig3 < eig1 && eig1 < eig2))
            return eig1;

        else
            return eig3;
    }
    else
    {
        const double q = (a[0] + a[4] + a[8]) / 3;
        const double p2 = pow(a[0] - q, 2) + pow(a[4] - q, 2) + pow(a[8] - q, 2) + 2 * p1;
        const double p = sqrt(p2 / 6);
        const double B[9]{
            a[0] - q,
            a[1],
            a[2],
            a[3],
            a[4] - q,
            a[5],
            a[6],
            a[7],
            a[8] - q};
        const double r = (B[0] * B[4] * B[8] + B[1] * B[5] * B[6] + B[3] * B[4] * B[7] - B[0] * B[5] * B[7] - B[1] * B[3] * B[8] - B[2] * B[4] * B[6]) / 2;
        double phi;
        if (r <= -1)
            phi = constants::PI / 3;
        else if (r >= 1)
            phi = 0;
        else
            phi = acos(r) / 3;

        eig1 = q + 2 * p * cos(phi);
        eig3 = q + 2 * p * cos(phi + 2 * constants::PI / 3);
        eig2 = 3 * q - eig1 - eig3;
        if ((eig1 < eig2 && eig2 < eig3) || (eig3 < eig2 && eig2 < eig1))
            return eig2;

        else if ((eig2 < eig1 && eig1 < eig3) || (eig3 < eig1 && eig1 < eig2))
            return eig1;

        else
            return eig3;
    }
};

const double gaussian_radial(const primitive &p, const double &r)
{
    return pow(r, p.get_type()) * std::exp(-p.get_exp() * r * r) * p.normalization_constant();
}

double get_bessel_ratio(const double nu, const double x)
{
    const double RECUR_BIG = DBL_MAX;
    const double RECUR_SMALL = DBL_MIN;
    const int maxiter = 10000;
    int n = 1;
    double Anm2 = 1.0;
    double Bnm2 = 0.0;
    double Anm1 = 0.0;
    double Bnm1 = 1.0;
    double a1 = x / (2.0 * (nu + 1.0));
    double An = Anm1 + a1 * Anm2;
    double Bn = Bnm1 + a1 * Bnm2;
    double an;
    double fn = An / Bn;
    double dn = a1;
    double s = 1.0;

    while (n < maxiter)
    {
        double old_fn;
        double del;
        n++;
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        an = -x * x / (4.0 * (nu + n - 1.0) * (nu + n));
        An = Anm1 + an * Anm2;
        Bn = Bnm1 + an * Bnm2;

        if (fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG)
        {
            An /= RECUR_BIG;
            Bn /= RECUR_BIG;
            Anm1 /= RECUR_BIG;
            Bnm1 /= RECUR_BIG;
            Anm2 /= RECUR_BIG;
        }
        else if (fabs(An) < RECUR_SMALL || fabs(Bn) < RECUR_SMALL)
        {
            An /= RECUR_SMALL;
            Bn /= RECUR_SMALL;
            Anm1 /= RECUR_SMALL;
            Bnm1 /= RECUR_SMALL;
            Anm2 /= RECUR_SMALL;
            Bnm2 /= RECUR_SMALL;
        }

        old_fn = fn;
        fn = An / Bn;
        del = old_fn / fn;

        dn = 1.0 / (2.0 * (nu + n) / x - dn);
        if (dn < 0.0)
            s = -s;

        if (fabs(del - 1.0) < 2.0 * DBL_EPSILON)
            break;
    }

    return fn;
}

double bessel_first_kind(const int l, const double x)
{
    if (l < 0 || x < 0.0)
    {
        err_not_impl_f("This is not implemented, pelase dont do this to me!", std::cout);
        return -1000;
    }
    else if (x == 0.0)
    {
        return (l > 0 ? 0.0 : 1.0);
    }
    else if (l == 0)
    {
        return sin(x) / x;
    }
    else if (l == 1)
    {
        return (sin(x) / x - cos(x)) / x;
    }
    else if (l == 2)
    {
        const double f = (3.0 / (x * x) - 1.0);
        return (f * sin(x) - 3.0 * cos(x) / x) / x;
    }
    else if (l == 3)
    {
        double x2 = x * x;
        const double f1 = (x2 - 15.0);
        const double f2 = (6. * x2 - 15.);
        return (-f2 * sin(x) + f1 * cos(x) * x) / pow(x, 4);
    }
    else if (l == 4)
    {
        double x2 = x * x;
        const double f1 = (10. * x2 - 105.0);
        const double f2 = (x2 * x2 - 45. * x2 + 105.);
        return (f2 * sin(x) + f1 * cos(x) * x) / pow(x, 5);
    }
    else if (l == 5)
    {
        double x2 = x * x;
        double x4 = x2 * x2;
        const double f1 = (-x4 + 105.0 * x2 - 945.);
        const double f2 = (15. * x4 - 420. * x2 + 945.);
        return (f2 * sin(x) + f1 * cos(x) * x) / (x4 * x2);
    }
    else if (l == 6)
    {
        double x2 = x * x;
        double x4 = x2 * x2;
        const double f1 = 21. * (x4 - 60.0 * x2 + 495.);
        const double f2 = (-x4 * x2 + 210. * x4 - 4725 * x2 + 10395.);
        return (f2 * sin(x) - f1 * cos(x) * x) / pow(x, 7);
    }
    else
    {
        double ratio = get_bessel_ratio(l + 0.5, x);
        const double smallest = DBL_MIN / DBL_EPSILON;
        double jellp1 = smallest * ratio;
        double jell = smallest;
        double jellm1;
        int ell;
        for (ell = l; ell > 0; ell--)
        {
            jellm1 = -jellp1 + (2 * ell + 1) / x * jell;
            jellp1 = jell;
            jell = jellm1;
        }

        if (fabs(jell) > fabs(jellp1))
        {
            double pre = smallest / jell;
            return bessel_first_kind(0, x) * pre;
        }
        else
        {
            double pre = smallest / jellp1;
            return bessel_first_kind(1, x) * pre;
        }
    }
}

int load_basis_into_WFN(WFN &wavy, std::shared_ptr<BasisSet> b)
{
    wavy.set_basis_set_ptr((*b).get_data());
    int nr_coefs = 0;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        int current_charge = wavy.get_atom_charge(i) - 1;
        const std::span<const SimplePrimitive> basis = (*b)[current_charge];
        int size = (int)basis.size();
        for (int e = 0; e < size; e++)
        {
            wavy.push_back_atom_basis_set(i, basis[e].exp, basis[e].coefficient, basis[e].type, basis[e].shell);
            //wavy.push_back_atom_basis_set(i, basis[e].exp, 1.0, basis[e].type, e);

            nr_coefs += 2 * basis[e].type + 1; //HIER WEITER MACHEN!!
        }
    }
    return nr_coefs;
}


double get_decimal_precision_from_CIF_number(std::string &given_string)
{
    int len = (int)given_string.length();
    int open_bracket = -1;
    int close_bracket = -1;
    int decimal_point = -1;
    // const char* gs = given_string.c_str();
    for (int i = 0; i < len; i++)
    {
        if (given_string[i] == '(' && open_bracket == -1)
        {
            open_bracket = i;
        }
        else if (given_string[i] == ')' && close_bracket == -1)
        {
            close_bracket = i;
        }
        else if (given_string[i] == '.' && decimal_point == -1)
        {
            decimal_point = i;
        }
    }
    double result = 0;
    int precision = 1;
    int size_of_precision = 1;
    if (open_bracket != -1 && close_bracket != -1)
    {
        size_of_precision = close_bracket - open_bracket - 1;
        std::string temp = given_string.substr(open_bracket + 1, size_of_precision);
        precision = std::stoi(temp);
    }
    int digits = 0;
    if (open_bracket != -1 && close_bracket != -1)
    {
        if (decimal_point != -1)
        {
            digits = open_bracket - decimal_point - 1;
        }
        else
        {
            digits = close_bracket - open_bracket - 1;
        }
        if (digits == 0)
        {
            return 0.001;
        }
        result = abs(precision * pow(10, -digits));
        return result;
    }
    else
        return 0.005;
};

void options::digest_options()
{
    using namespace std;
    // Lets print what was the command line, for debugging
    if (debug)
    {
       std::cout << " Recap of input:\nsize: " << arguments.size() << endl;
    }
    // This loop figures out command line options
    int argc = (int)arguments.size();
    for (int i = 0; i < arguments.size(); i++)
    {
        if (debug)
           std::cout << arguments[i] << endl;
        string temp = arguments[i];
        if (temp.find("-") > 0)
            continue;
        if (temp == "-acc")
            accuracy = stoi(arguments[i + 1]);
        else if (temp == "-Anion")
        {
            int n = 1;
            string store;
            if (debug)
               std::cout << "Looking for Anions!" << endl;
            while (i + n < argc && string(arguments[i + n]).find("-") != 0)
            {
                if (i + n - 1 > arguments.size())
                    break;
                store = arguments[i + n];
                svec Z = split_string<string>(store, " ");
                for (int r = 0; r < Z.size(); r++)
                {
                    if (debug)
                       std::cout << Z[r] << endl;
                    Anions.push_back(Z[r]);
                }
                n++;
            }
        }
        else if (temp == "-laplacian_bonds")
        {
            wfn = arguments[i + 1];
            bondwise_laplacian_plots(wfn);
            exit(0);
        }
        else if (temp == "-atom_dens")
        {
           std::cout << NoSpherA2_message() << endl;
            wfn = arguments[i + 1];
            err_checkf(std::filesystem::exists(wfn), "WFN doesn't exist",std::cout);
            ivec val_MOs;
            ivec val_MOs_beta;
            if (argc >= i + 3)
            {
                val_MOs = split_string<int>(arguments[i + 2], ",");
               std::cout << "Alpha MOs to keep: ";
                for (int j = 0; j < val_MOs.size(); j++)
                   std::cout << val_MOs[j] << " ";
               std::cout << endl;
                val_MOs_beta = split_string<int>(arguments[i + 3], ",");
               std::cout << "Beta MOs to keep: ";
                for (int j = 0; j < val_MOs_beta.size(); j++)
                   std::cout << val_MOs_beta[j] << " ";
               std::cout << endl;
            }
            spherically_averaged_density(*this, val_MOs, val_MOs_beta);
            exit(0);
        }
        else if (temp == "-b")
            basis_set = arguments[i + 1];
        else if (temp == "-becke" || temp == "-BECKE" || temp == "-Becke")
            partition_type = PartitionType::Becke;
        else if (temp == "-blastest")
        {
            test_openblas();
            test_solve_linear_equations();
            exit(0);
        }
        else if (temp == "-lahvatest")
        {
            //_test_lahva();
            exit(0);
        }
        else if (temp == "-Cation")
        {
            int n = 1;
            string store;
            if (debug)
               std::cout << "Looking for Cations!" << endl;
            while (i + n < argc && string(arguments[i + n]).find("-") != 0)
            {
                if (i + n - 1 > arguments.size())
                    break;
                store = arguments[i + n];
                svec Z = split_string<string>(store, " ");
                for (int r = 0; r < Z.size(); r++)
                {
                    if (debug)
                       std::cout << Z[r] << endl;
                    Cations.push_back(Z[r]);
                }
                n++;
            }
        }
        else if (temp == "-charge")
        {
            charge = stoi(arguments[i + 1]);
        }
        else if (temp == "-coef")
        {
            coef_file = arguments[i + 1];
            err_checkf(std::filesystem::exists(coef_file), "coef_file doesn't exist",std::cout);
            SALTED = true;
        }
        else if (temp == "-cif")
        {
            cif = arguments[i + 1];
            err_checkf(std::filesystem::exists(cif), "CIF doesn't exist",std::cout);
        }
        else if (temp == "-cpus")
        {
            threads = stoi(arguments[i + 1]);
            MKL_Set_Num_Threads(threads);
#ifdef _OPENMP
            omp_set_num_threads(threads);
            omp_set_dynamic(0);

#endif
        }
        else if (temp == "-cmtc")
        {
            cif_based_combined_tsc_calc = true;
            int n = 1;
            string delimiter = ",";
            groups.pop_back();
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                combined_tsc_calc_files.push_back(arguments[i + n]);
                n++;
                combined_tsc_calc_cifs.push_back(arguments[i + n]);
                n++;
                const string _temp = arguments[i + n];
                if (_temp.find(delimiter) == string::npos) {
                    if (debug)
                        std::cout << "--Delimiter not found, using ." << endl;
                    delimiter = ".";
                }
                groups.push_back(split_string<int>(_temp, delimiter));
                if (debug)
                {
                    std::cout << "--Group: " << _temp << endl << "--";
                    for (int run = 0; run < groups[groups.size() - 1].size(); run++)
                        std::cout  << groups[groups.size() - 1][run] << " ";
                    std::cout << endl;
                }
                n++;
            }
        }
        else if (temp == "-combine_mos")
        {
            combine_mo.push_back(arguments[i + 1]);
            combine_mo.push_back(arguments[i + 2]);
            do_combine_mo(*this);
            exit(0);
        }
        else if (temp == "-cmos1")
        {
            int j = 1;
            while (i + j < argc && arguments[i + j].find("-") >= 1)
            {
                cmo1.push_back(stoi(arguments[i + j]));
                j++;
            }
        }
        else if (temp == "-cmos2")
        {
            int j = 1;
            while (i + j < argc && arguments[i + j].find("-") >= 1)
            {
                cmo2.push_back(stoi(arguments[i + j]));
                j++;
            }
        }
        else if (temp == "-core_dens-corrected")
        {
            double prec = stod(arguments[i + 1]);
            ivec a = split_string<int>(arguments[i + 3], ",");
            ivec b = split_string<int>(arguments[i + 4], ",");
            test_core_dens_corrected(prec, threads, arguments[i + 2], a, b);
            exit(0);
        }
        else if (temp == "-core_tsc-corrected")
        {
            double prec = stod(arguments[i + 1]);
            ivec a = split_string<int>(arguments[i + 3], ",");
            ivec b = split_string<int>(arguments[i + 4], ",");
            test_core_sfac_corrected(prec, threads, arguments[i + 2], a, b);
            exit(0);
        }
        else if (temp == "-def" || temp == "-DEF")
            def = calc = true;
        else if (temp == "-density_difference" || temp == "-density-difference")
        {
            wfn2 = arguments[i + 1];
        }
        else if (temp == "-dmin")
            dmin = stod(arguments[i + 1]);
        else if (temp == "-d")
            basis_set_path = arguments[i + 1];
        else if (temp == "-dipole_moments")
        {
            dipole_moments(*this);
            exit(0);
        }
        // Visualize the specified orbital using spherical harmonics.
        // Call as -draw_orbits lambda,m,resolution,radius
        // Where resolution and radius are optional
        else if (temp == "-draw_orbits")
        {
            vec opts = split_string<double>(arguments[i + 1], ",");
            int l = static_cast<int>(opts[0]);
            int m = static_cast<int>(opts[1]);
            resolution = 0.025;
            radius = 3.5;
            if (opts.size() >= 3)
            {
                resolution = opts[2];
            }
            if (opts.size() == 4)
            {
                radius = opts[3];
            }

            draw_orbital(l, m, resolution, radius);
            exit(0);
        }
        else if (temp == "-e_field")
            efield = stod(arguments[i + 1]);
        else if (temp == "-ECP" || temp == "-ecp" || temp == "-Ecp")
        {
            ECP = true;
            if (argc >= i + 2 && string(arguments[i + 1]).find("-") != 0)
            {
                ECP_mode = stoi(arguments[i + 1]);
            }
        }
        else if (temp == "-ED")
            electron_diffraction = true;
        else if (temp == "-eli")
            calc = eli = true;
        else if (temp == "-elf")
            calc = elf = true;
        else if (temp == "-equi_bench") {
            exit(0);
        }
        else if (temp == "-esp")
            calc = esp = true;
        else if (temp == "-ewal_sum")
        {
            // bool read, WFN& wave, std::ostream& file,
            WFN *temp_w = new WFN(9);
            cube residual(arguments[i + 1], true, *temp_w, std::cout);
            if (argc >= i + 3)
            {
                int k_max = stoi(arguments[i + 2]);
                if (argc >= i + 4)
                    residual.ewald_sum(k_max, stod(arguments[i + 3]));
                else
                    residual.ewald_sum(k_max);
            }
            else
                residual.ewald_sum();
            delete (temp_w);
            exit(0);
        }
        else if (temp == "-fchk")
            fchk = arguments[i + 1];
        else if (temp == "-fractal")
            fract = true, fract_name = arguments[i + 1];
        else if (temp == "-gbw2wfn")
            gbw2wfn = true;
        else if (temp == "-group")
        {
            int n = 1;
            while (i + n < argc && string(arguments[i + n]).find("-") == string::npos)
            {
                int group;
                if (arguments[i + 1][0] == '+')
                    group = -stoi(arguments[i + n]);
                else
                    group = stoi(arguments[i + n]);
                groups[0].push_back(group), n++;
            }
            i += n;
        }
        else if (temp == "-HDEF")
            hdef = calc = true;
        else if (temp == "-hirsh")
            calc = hirsh = true, hirsh_number = stoi(arguments[i + 1]);
        else if (temp == "-hirshfeld_surface")
        {
            hirshfeld_surface = arguments[i + 1];
            hirshfeld_surface2 = arguments[i + 2];
        }
        else if (temp == "-hkl")
        {
            hkl = arguments[i + 1];
            err_checkf(std::filesystem::exists(hkl), "hkl doesn't exist",std::cout);
        }
        else if (temp == "-hkl_min_max")
        {
            int h_min(stoi(arguments[i + 1]));
            int h_max(stoi(arguments[i + 2]));
            int k_min(stoi(arguments[i + 3]));
            int k_max(stoi(arguments[i + 4]));
            int l_min(stoi(arguments[i + 5]));
            int l_max(stoi(arguments[i + 6]));
            hkl_min_max = {{h_min, h_max}, {k_min, k_max}, {l_min, l_max}};
        }
        else if (temp == "-IAM")
            iam_switch = true;
        else if (temp == "-lap")
            calc = lap = true;
        else if (temp == "-mem")
        {
            mem = stod(arguments[i + 1]); // In MB
            vec a;
            size_t vec_max_size = a.max_size();
            double doubel_max_size = static_cast<double>(vec_max_size * sizeof(double)) * 1e-6;
            if (mem > doubel_max_size)
            {
               std::cout << "Max memory set to " << mem << " MB, which is larger than the maximum allowed size of " << doubel_max_size << " MB. Setting max memory to " << 50000 << " MB." << endl;
                mem = 50000.0;
            }
        }
        else if (temp == "-method")
            method = arguments[i + 1];
        else if (temp == "-merge")
        {
            pathvec filenames;
            int n = 1;
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                filenames.push_back(arguments[i + n]);
                n++;
            }
            merge_tscs("combine", filenames, old_tsc);
            exit(0);
        }
        else if (temp == "-merge_nocheck")
        {
            pathvec filenames;
            int n = 1;
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                filenames.push_back(arguments[i + n]);
                n++;
            }
            merge_tscs_without_checks("combine", filenames, old_tsc);
            exit(0);
        }
        else if (temp == "-MO")
        {
            if (string(arguments[i + 1]) != "all")
                MOs.push_back(stoi(arguments[i + 1]));
            else
                all_mos = true;
            calc = true;
        }
        else if (temp == "-mtc")
        {
            combined_tsc_calc = true;
            int n = 1;
            string delimiter = ",";
            groups.pop_back();
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                combined_tsc_calc_files.push_back(arguments[i + n]);
                n++;
                const string _temp = arguments[i + n];
                if (_temp.find(delimiter) == string::npos) {
                    if (debug)
                        std::cout << "--Delimiter not found, using ." << endl;
                    delimiter = ".";
                }
                groups.push_back(split_string<int>(_temp, delimiter));
                if (debug)
                {
                    std::cout << "--Group: " << _temp << endl << "--";
                    for (int run = 0; run < groups[groups.size() - 1].size(); run++)
                        std::cout << groups[groups.size() - 1][run] << " ";
                    std::cout << endl;
                }
                n++;
            }
        }
        else if (temp == "-mtc_mult")
        {
            int n = 1;
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                combined_tsc_calc_mult.push_back(stoi(arguments[i + n]));
                n++;
            }
        }
        else if (temp == "-mtc_charge")
        {
            int n = 1;
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                if (arguments[i + n][0] == 'n')
                    combined_tsc_calc_charge.push_back(-stoi(arguments[i + n].substr(1)));
                else
                    combined_tsc_calc_charge.push_back(stoi(arguments[i + n]));
                n++;
            }
        }
        else if (temp == "-mtc_ECP")
        {
            int m = 1;
            while (i + m < argc && string(arguments[i + m]).find("-") > 0)
            {
                combined_tsc_calc_ECP.push_back(stoi(arguments[i + m]));
                m++;
            }
        }
        else if (temp == "-mult")
            mult = stoi(arguments[i + 1]);
        else if (temp == "-NNLS_TEST")
        {
            test_NNLS();
            exit(0);
        }
        else if (temp == "-no-date" || temp == "-no_date")
            no_date = true;
        else if (temp == "-pbc")
            pbc = stoi(arguments[i + 1]);
        else if (temp == "-polarizabilities")
        {
            pol_wfns = {arguments[i + 1],
                        arguments[i + 2],
                        arguments[i + 3],
                        arguments[i + 4],
                        arguments[i + 5],
                        arguments[i + 6],
                        arguments[i + 7]};
        }
        else if (temp == "-QCT" || temp == "-qct")
            qct = true;
        else if (temp == "-radius")
            radius = stod(arguments[i + 1]);
        else if (temp == "-resolution")
            resolution = stod(arguments[i + 1]);
        else if (temp == "-rdg")
            calc = rdg = true;
        else if (temp.find("-rkpts") < 1)
            read_k_pts = true;
        else if (temp == "-rho_cube")
        {
            string wfn_name = arguments[i + 1];
            std::cout << "Reading wavefunction: " << wfn_name << endl;
            WFN *wavy = new WFN(wfn_name);
            std::cout << "Assigning ECPs" << endl;
            if (ECP)
                wavy->set_has_ECPs(true);
            std::cout << "Starting cube calculation" << endl;
            calc_rho_cube(*wavy);
            delete (wavy);
            exit(0);
        }
        else if (temp == "-RI_FIT" || temp == "-ri_fit")
        {
            RI_FIT = true;
            partition_type = PartitionType::RI;
            // Check if next argument is a valid basis set name or a new argument starting with "-"
            if (i + 1 < argc && arguments[i + 1].find("-") != 0)
            {

                if (arguments[i + 1] == "auto_aux") {
                    double beta = 2.0;
                    //Check if the next argument is a valid double
                    if (i + 2 < argc && arguments[i + 2].find("-") != 0) {
                        beta = std::stod(arguments[i + 2]);
                    }
                    if (debug) cout << "Using automatic basis set generation with beta: " << beta << endl;
                    aux_basis = std::make_shared<BasisSet>(beta);
                    continue;
                }


                if (!BasisSetLibrary().check_basis_set_exists(arguments[i + 1]))
                {
                    cout << "Basis set " << arguments[i + 1] << " not found in the library. Exiting." << endl;
                    exit(1);
                }
                aux_basis = BasisSetLibrary().get_basis_set(arguments[i + 1]);
            }
            else
            {
                cout << "No basis set specified. Falling back to automatic generation using beta = 2.0!" << endl;
                aux_basis = std::make_shared<BasisSet>(2.0);
            }
        }
        else if (temp == "-RI_CUBE" || temp == "-ri_cube")
        {
            WFN wavy(wfn);
            // First name of coef_file, second name of xyz file
            // cube_from_coef_npy(arguments[i + 1], arguments[i + 2]);

            // std::string aux_basis = arguments[i + 1];
            gen_CUBE_for_RI(wavy, "def2_qzvppd_rifit", this);
            //gen_CUBE_for_RI(wavy, "def2_universal_jkfit", this);
            //gen_CUBE_for_RI(wavy, "combo-basis-fit", this);
            //gen_CUBE_for_RI(wavy, "cc-pvqz-jkfit", this);

            exit(0);
        }

        else if (temp.find("-s_rho") < 1)
            s_rho = true;
        else if (temp == "-SALTED" || temp == "-salted")
        {
            SALTED = true;
            salted_model_dir = arguments[i + 1];
        }
        else if (temp == "-SALTED_COEFS" || temp == "-salted_coefs")
        {
            salted_model_dir = arguments[i + 1];

            //Check that wfn is not empty
            if (wfn.empty())
            {
                std::cout << "No wavefunction specified! Use -wfn option BEVORE -test_RI to specify a wavefunction." << std::endl;
                exit(1);
            }

            WFN wavy(wfn);
            SALTEDPredictor SP(wavy, *this);
            string df_basis_name = SP.get_dfbasis_name();
            filesystem::path salted_model_path = SP.get_salted_filename();
            log_file << "Using " << salted_model_path << " for the prediction" << endl;
            std::shared_ptr<BasisSet> aux_basis = BasisSetLibrary().get_basis_set(df_basis_name);
            load_basis_into_WFN(SP.wavy, aux_basis);
            vec coefs = SP.gen_SALTED_densities();
            npy::npy_data<double> np_coeffs;
            np_coeffs.data = coefs;
            np_coeffs.fortran_order = false;
            np_coeffs.shape = { unsigned long(coefs.size()) };
            npy::write_npy("SALTED_COEFS.npy", np_coeffs);
            }
        else if (temp == "-test_reading_SALTED_binary") {
            test_reading_SALTED_binary_file();
            exit(0);
        }
        else if (temp == "-skpts")
            save_k_pts = true;
        else if (temp == "-sfac_scan")
        {
            d_sfac_scan = fromString<double>(arguments[i + 1]);
            cif = arguments[i + 2];
            wfn = arguments[i + 3];
            sfac_scan(*this, log_file);
            exit(0);
        }
        else if (temp == "-sfac_diffuse")
        {
            sfac_diffuse = fromString<double>(arguments[i + 1]);
            cif = arguments[i + 2];
            wfn = arguments[i + 3];
            dmin = fromString<double>(arguments[i + 4]);
            calc_sfac_diffuse(*this, std::cout);
        }
        else if (temp == "-atom_dens_diff")
        {
            filesystem::path name_wfn_1 = arguments[i + 1];
            filesystem::path name_wfn_2 = arguments[i + 2];

            subtract_dens_from_gbw(name_wfn_1, name_wfn_2, 2, 0.05);
            exit(0);
        }
        else if (temp == "-spherical_aver_fukui")
        {
            filesystem::path wfn1_name = arguments[i + 1];
            filesystem::path wfn2_name = arguments[i + 2];
            WFN *wavy1 = new WFN(wfn1_name);
            WFN *wavy2 = new WFN(wfn2_name);
            ofstream outputFile("fukui_averaged_density_wfn.dat");
            for (double r = 0.001; r < 10.0; r += 0.001)
            {
                // double dens = calc_grid_averaged_at_r_from_cube(cube_from_file, r, 360, 5800);
                double dens = calc_fukui_averaged_at_r(*wavy1, *wavy2, r, 5810, 5810);
                outputFile << r << " " << dens << "\n";
            }
            outputFile.close();
            std::cout << "Data written to output.dat" << endl;
            delete (wavy1);
            delete (wavy2);
            exit(0);
        }
        else if (temp == "-spherical_aver_hirsh")
        {
            string wfn_name = arguments[i + 1];
            std::cout << "Reading wavefunction: " << wfn_name << endl;
            WFN *wavy = new WFN(wfn_name);
            std::cout << "Assigning ECPs" << endl;
            if (ECP)
                wavy->set_has_ECPs(true);
            std::cout << "Starting spherical averaging" << endl;
            double dens;

            for (int index_atom = 0; index_atom < wavy->get_ncen(); index_atom += 1)
            {
                ofstream outputFile("hirsh_averaged_density_" + std::to_string(index_atom) + ".dat");
                for (double r = 0.001; r < 5.0; r += 0.002)
                {
                    dens = calc_hirsh_grid_averaged_at_r(*wavy, index_atom, r, 360, 5800);
                    outputFile << r << " " << dens << "\n";
                }
                outputFile.close();
            }
            std::cout << "Data written to output.dat" << endl;
            delete (wavy);
            exit(0);
        }
        else if (temp == "-spherical_harmonic")
        {
            spherical_harmonic_test();
            exit(0);
        }
        else if (temp == "-test")
            std::cout << "Running in test mode!" << endl, test = true;
        else if (temp == "-twin")
        {
            twin_law.resize(twin_law.size() + 1);
            twin_law.back().resize(9);
            for (int twl = 0; twl < 9; twl++)
                twin_law.back()[twl] = stod(arguments[i + 1 + twl]);
            if (debug)
            {
               std::cout << "twin_law: ";
                for (int twl = 0; twl < 9; twl++)
                   std::cout << setw(7) << setprecision(2) << twin_law.back()[twl];
               std::cout << endl;
            }
            i += 9;
        }
        else if (temp == "-old_tsc")
        {
            old_tsc = true;
        }
        else if (temp == "-tfvc" || temp == "-TFVC")
        {
            partition_type = PartitionType::TFVC;
        }
        else if (temp == "-tscb")
        {
            std::filesystem::path name = arguments[i + 1];
            tsc_block<int, cdouble> blocky = tsc_block<int, cdouble>(name);
            string cif_name = "test.cif";
            if (name.extension() == ".tscb")
                blocky.write_tsc_file(cif_name, name.replace_extension(".tsc"));
            else if (name.extension() == ".tsc")
                blocky.write_tscb_file(cif_name, name.replace_extension(".tscb"));
            else
                err_checkf(false, "Wrong file ending!",std::cout);
            exit(0);
        }
        else if (temp == "-test_analytical")
        {
            bool full = false;
            if ("full" == arguments[i + 1]) {
                full = true;
            }
            test_analytical_fourier(full);
            exit(0);
        }
        else if (temp == "-test_RI")
        {
            //Check that wfn is not empty
            if (wfn.empty())
            {
                std::cout << "No wavefunction specified! Use -wfn option BEVORE -test_RI to specify a wavefunction." << std::endl;
                exit(1);
            }
            if (aux_basis == nullptr)
            {
                std::cout << "No auxiliary basis set specified! Use -RI_FIT option BEVORE -test_RI to specify an auxiliary basis set." << std::endl;
                exit(1);
            }
  

            WFN wavy(wfn);

            WFN wavy_aux(0);
            wavy_aux.set_atoms(wavy.get_atoms());
            wavy_aux.set_ncen(wavy.get_ncen());
            wavy_aux.delete_basis_set();

            if ((*aux_basis).get_primitive_count() == 0) (*aux_basis).gen_auto_aux(wavy);
            load_basis_into_WFN(wavy_aux, aux_basis);
            demonstrate_enhanced_density_fitting(wavy, wavy_aux);
            //dMatrix2 dm = wavy.get_dm();
            //vec charges_sand = calculate_expected_charges(wavy, wavy_aux, "sanderson_estimate");
            //vec charges_mul = calculate_expected_charges(wavy, wavy_aux, "mulliken");

            //for (int i = 0; i < charges_mul.size(); i++)
            //    std::cout << "Charge on atom " << i << " (" << wavy.get_atom_label(i) << "): Mulliken: " << charges_mul[i] << " Sanderson: " << charges_sand[i] << std::endl;

            //aux_basis = std::make_shared<BasisSet>();
            //aux_basis->gen_auto_aux(wavy);

            //WFN wavy_aux2(0);
            //wavy_aux2.set_atoms(wavy.get_atoms());
            //wavy_aux2.set_ncen(wavy.get_ncen());
            //wavy_aux2.delete_basis_set();
            //load_basis_into_WFN(wavy_aux2, aux_basis);
            //demonstrate_enhanced_density_fitting(wavy, wavy_aux2);

            exit(0);
        }
        else if (temp == "-wfn")
        {
            wfn = arguments[i + 1];
            err_checkf(std::filesystem::exists(wfn), "Wavefunction dos not exist!",std::cout);
        }
        else if (temp == "-wfn_cif")
        {
            write_CIF = true;
        }
        else if (temp == "-xtb_test")
        {
            test_xtb_molden(*this, log_file);
            exit(0);
        }
        else if (temp == "-xyz")
        {
            xyz_file = arguments[i + 1];
        }
        else if (temp == "-partitioning_test")
        {
            calc_partition_densities();
        }
    }
};

void options::look_for_debug(int &argc, char **argv)
{
    // This loop figures out command line options
    for (int i = 0; i < argc; i++)
    {
        std::string temp = argv[i];
        arguments.push_back(temp);
        if (temp.find("-") > 0)
            continue;
        else if (temp == "-v" || temp == "-v2" || temp == "-debug")
            std::cout << "Turning on verbose mode!" << std::endl, debug = true;
        else if (temp == "--h" || temp == "-h" || temp == "-help" || temp == "--help")
        {
            std::cout << NoSpherA2_message() << help_message << build_date << std::endl;
            exit(0);
        }
    }
};

bool is_nan(double &in)
{
    return in != in;
};
bool is_nan(float &in)
{
    return in != in;
};
bool is_nan(long double &in)
{
    return in != in;
};
bool is_nan(cdouble &in)
{
    return in != in;
};

bool ends_with(const std::string &str, const std::string &suffix)
{
    if (str.length() >= suffix.length())
    {
        return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
    }
    return false;
}

cdouble hypergeometric(double a, double b, double c, cdouble x)
{
    const double TOLERANCE = 1.0e-10;
    cdouble term = a * b * x / c;
    cdouble value = 1.0 + term;
    int n = 1;

    while (std::abs(term) > TOLERANCE)
    {
        a++, b++, c++, n++;
        term *= a * b * x / c / static_cast<double>(n);
        value += term;
    }

    return value;
}

double hypergeometric(double a, double b, double c, double x)
{
    const double TOLERANCE = 1.0e-10;
    double term = a * b * x / c;
    double value = 1.0 + term;
    int n = 1;

    while (std::abs(term) > TOLERANCE)
    {
        a++, b++, c++, n++;
        term *= a * b * x / c / n;
        value += term;
    }

    return value;
}

double double_from_string_with_esd(std::string in)
{
    if (in.find('(') == std::string::npos)
        return stod(in);
    else
        return stod(in.substr(0, in.find('(')));
}

std::string trim(const std::string &s)
{
    if (s == "")
        return "";
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start))
    {
        start++;
    }

    auto end = s.end();
    do
    {
        end--;
    } while (std::distance(start, end) > 0 && std::isspace(*end));

    return std::string(start, end + 1);
}

int CountWords(const char *str)
{
    if (str == NULL)
        return -1;

    bool inSpaces = true;
    int numWords = 0;

    while (*str != '\0')
    {
        if (std::isspace(*str))
        {
            inSpaces = true;
        }
        else if (inSpaces)
        {
            numWords++;
            inSpaces = false;
        }

        ++str;
    }

    return numWords;
};

void print_duration(std::ostream &file, const std::string &description, const std::chrono::microseconds &duration, std::optional<std::chrono::microseconds> total_duration = std::nullopt)
{
    auto mins = std::chrono::duration_cast<std::chrono::minutes>(duration);
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(duration) % 60;
    auto millisecs = std::chrono::duration_cast<std::chrono::milliseconds>(duration) % 1000;

    file << std::setw(35) << std::left << std::setfill(' ') << description << ": " << std::right
         << std::setw(2) << std::setfill('0') << mins.count() << ":"
         << std::setw(2) << std::setfill('0') << secs.count() << ":"
         << std::setw(3) << std::setfill('0') << millisecs.count();
    if (total_duration.has_value())
    {
        double percentage = (double(duration.count()) / total_duration->count()) * 100.0;
        std::cout << "  (" << std::fixed << std::setprecision(2) << percentage << "%)";
    };
    file << std::endl;
    // Disable setfill 0 again
    file << std::setfill(' ');
}

void write_timing_to_file(std::ostream &file,
                          std::vector<_time_point> time_points,
                          std::vector<std::string> descriptions)
{
    using namespace std;
    // Check if either vector is empty
    if (time_points.empty() || descriptions.empty())
    {
       std::cout << "Error: Empty vector passed to write_timing_to_file" << endl;
        return;
    }

    std::chrono::microseconds total_time = std::chrono::duration_cast<std::chrono::microseconds>(time_points.back() - time_points.front());
    file << "\n\n----------------------------- Time Breakdown! -----------------------------" << endl;
    file << "                                     mm:ss:ms" << endl;
    // Time_points.size()-1 because we are comparing each time point to the next one meaning we need to stop at the second to last element
    for (int i = 0; i < time_points.size() - 1; i++)
    {
        std::chrono::microseconds dur = std::chrono::duration_cast<std::chrono::microseconds>(time_points[i + 1] - time_points[i]);
        print_duration(file, "... for " + descriptions[i], dur, total_time);
    }
    print_duration(file, "Total Time", total_time);
    file << "---------------------------------------------------------------------------" << endl;
}

void remove_empty_elements(svec &input, const std::string &empty)
{
    for (int i = (int)input.size() - 1; i >= 0; i--)
        if (input[i] == empty || input[i] == "")
            input.erase(input.begin() + i);
}

std::chrono::high_resolution_clock::time_point get_time()
{
    // gets the current time using std chrono library
    std::chrono::high_resolution_clock::time_point time = std::chrono::high_resolution_clock::now();
    return time;
}

long long int get_musec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
    // gets the time difference in microseconds
    std::chrono::microseconds musec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    return musec.count();
}

long long int get_msec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
    // gets the time difference in milliseconds
    std::chrono::milliseconds msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return msec.count();
}

long long int get_sec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
    // gets the time difference in seconds
    std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    return sec.count();
}

const int shell2function(const int &type, const int &prim)
{
    switch (type)
    {
    case (-5):
        return -32 + prim;
    case (-4):
        return -21 + prim;
    case (-3):
        return -12 + prim;
    case (-2):
        return -5 + prim;
    case (-1):
        return 1 + prim;
    case (0):
        return 1;
    case (1):
        return 2 + prim;
    case (2):
        return 5 + prim;
    case (3):
        if (prim == 0)
            return 11;
        if (prim == 1)
            return 12;
        if (prim == 2)
            return 13;
        if (prim == 3)
            return 17;
        if (prim == 4)
            return 14;
        if (prim == 5)
            return 15;
        if (prim == 6)
            return 18;
        if (prim == 7)
            return 19;
        if (prim == 8)
            return 16;
        if (prim == 9)
            return 20;
        break;
    case (4):
        return 21 + prim;
    case (5):
        return 36 + prim;
    default:
        return 0;
    }
    return 0;
}

bool open_file_dialog(std::filesystem::path& path, bool debug, std::vector <std::string> filter, const std::string& current_path) {
#ifdef _WIN32
    char filename[1024];

    OPENFILENAMEA ofn;
    ZeroMemory(&filename, sizeof(filename));
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.lpstrFilter = "Known Formats\0*.gbw;*.fchk;*.wfx;*.wfn;*.ffn;*.molden;*.molden.input;*.xtb\0gbw Files\0*.gbw\0wfn Files\0*.wfn\0wfx Files\0*.wfx\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0xtb Files\0*.xtb\0Any File\0*\0";
    ofn.lpstrFile = filename;
    ofn.nMaxFile = 1024;
    ofn.lpstrTitle = "Select a File";
    ofn.Flags = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST;

    if (GetOpenFileNameA(&ofn)) {
        if (debug) std::cout << "You chose the file \"" << filename << "\"\n";
        auto p = std::filesystem::path(filename);
        if (exists(p)) {
            path = p;
            return true;
        }
    }
    else
    {
        // All this stuff below is to tell you exactly how you messed up above.
        // Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
        switch (CommDlgExtendedError())
        {
        case CDERR_DIALOGFAILURE: std::cout << "CDERR_DIALOGFAILURE\n";   break;
        case CDERR_FINDRESFAILURE: std::cout << "CDERR_FINDRESFAILURE\n";  break;
        case CDERR_INITIALIZATION: std::cout << "CDERR_INITIALIZATION\n";  break;
        case CDERR_LOADRESFAILURE: std::cout << "CDERR_LOADRESFAILURE\n";  break;
        case CDERR_LOADSTRFAILURE: std::cout << "CDERR_LOADSTRFAILURE\n";  break;
        case CDERR_LOCKRESFAILURE: std::cout << "CDERR_LOCKRESFAILURE\n";  break;
        case CDERR_MEMALLOCFAILURE: std::cout << "CDERR_MEMALLOCFAILURE\n"; break;
        case CDERR_MEMLOCKFAILURE: std::cout << "CDERR_MEMLOCKFAILURE\n";  break;
        case CDERR_NOHINSTANCE: std::cout << "CDERR_NOHINSTANCE\n";     break;
        case CDERR_NOHOOK: std::cout << "CDERR_NOHOOK\n";          break;
        case CDERR_NOTEMPLATE: std::cout << "CDERR_NOTEMPLATE\n";      break;
        case CDERR_STRUCTSIZE: std::cout << "CDERR_STRUCTSIZE\n";      break;
        case FNERR_BUFFERTOOSMALL: std::cout << "FNERR_BUFFERTOOSMALL\n";  break;
        case FNERR_INVALIDFILENAME: std::cout << "FNERR_INVALIDFILENAME\n"; break;
        case FNERR_SUBCLASSFAILURE: std::cout << "FNERR_SUBCLASSFAILURE\n"; break;
        default: std::cout << "You cancelled.\n";
        }
    }
    return false;
#else
    std::string command;
    bool use_zenity = (system("which zenity > /dev/null 2>&1") == 0);
    bool use_kdialog = false;

    if (use_zenity) {
        command = "zenity --file-selection --title=\"Select a file to load\" --filename=\"";
        command += current_path;
        command += "/\"";
        for (const auto& f : filter) {
            command += " --file-filter='";
            command += f;
            command += "'";
        }
        command += " 2> /dev/null";
    } else {
        use_kdialog = (system("which kdialog > /dev/null 2>&1") == 0);
        if (use_kdialog) {
            command = "kdialog --getopenfilename \"";
            command += current_path;
            command += "/\" '";
            for (const auto& f : filter) {
                command += f;
                command += " ";
            }
            command += "'";
            command += " --title \"Select a file to load\" 2> /dev/null";
        } else {
            std::cout << "No suitable file dialog tool found (zenity/kdialog)." << std::endl;
            std::cout << "Please enter the full path to the file: " << std::flush;
            std::string input_path;
            std::getline(std::cin, input_path);

            // Trim leading/trailing whitespace
            input_path.erase(0, input_path.find_first_not_of(" \t\n\r"));
            input_path.erase(input_path.find_last_not_of(" \t\n\r") + 1);

            if (input_path.empty()) {
                if (debug) std::cout << "No path entered." << std::endl;
                return false;
            }

            path = input_path;
            if (std::filesystem::exists(path)) {
                if (debug) std::cout << "Selected file via manual input: " << path << std::endl;
                return true;
            } else {
                std::cerr << "Error: File not found at path: " << path << std::endl;
                return false;
            }
        }
    }

    if (debug) {
        std::cout << "Executing command: " << command << std::endl;
    }

    FILE* f = popen(command.c_str(), "r");
    if (!f) {
        std::cerr << "Error: Failed to execute file dialog command." << std::endl;
        return false;
    }

    std::string file_str;
    char buffer[1024];
    if (fgets(buffer, sizeof(buffer), f) == NULL) {
        if (debug) std::cout << "File selection cancelled." << std::endl;
        pclose(f);
        return false;
    }
    file_str = buffer;

    int pclose_status = pclose(f);
    if (pclose_status != 0) {
        if (debug) std::cout << (use_zenity ? "Zenity" : "KDialog") << " returned non-zero status: " << pclose_status << ". User might have cancelled." << std::endl;
        // This can happen on cancel, so we check if a file was actually returned.
        if (file_str.empty()) return false;
    }

    // Clean up the path string which might have a newline
    file_str.erase(file_str.find_last_not_of(" \n\r\t")+1);

    if (file_str.empty()) {
        if (debug) std::cout << "File selection cancelled or returned empty path." << std::endl;
        return false;
    }

    path = file_str;
    if (debug) {
        std::cout << "Selected file: " << path << std::endl;
    }

    return std::filesystem::exists(path);
#endif
};

bool save_file_dialog(std::filesystem::path& path, bool debug, const std::vector<std::string>& endings, const std::string& filename_given, const std::string& current_path) {
#ifdef _WIN32
    constexpr size_t MAX_FILENAME_SIZE = 4096;
    std::vector<char> filename_buf(MAX_FILENAME_SIZE);

    OPENFILENAMEA sfn;
    ZeroMemory(&filename_buf[0], filename_buf.size());
    ZeroMemory(&sfn, sizeof(sfn));
    sfn.lStructSize = sizeof(sfn);
    sfn.hwndOwner = NULL;  // If you have a window to center over, put its HANDLE here
    sfn.lpstrFilter = "Known Formats\0*.gbw;*.fchk;*.wfx;*.wfn;*.ffn;*.molden;*.molden.input;*.xtb\0gbw Files\0*.gbw\0wfn Files\0*.wfn\0wfx Files\0*.wfx\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0xtb Files\0*.xtb\0Any File\0*\0";
    sfn.lpstrFile = filename_buf.data();
    sfn.nMaxFile = static_cast<DWORD>(filename_buf.size());
    sfn.lpstrTitle = "Select a File for saving!";
    sfn.Flags = OFN_DONTADDTORECENT;
    bool end = false;
    while (!end) {
        if (GetSaveFileNameA(&sfn)) {
            std::string chosen(filename_buf.data());
            if (debug) std::cout << "You chose the file \"" << chosen << "\"\n";
            if (exists(std::filesystem::path(chosen))) {
                std::cout << chosen << " exists, do you want to overwrite it?";
                if (yesno()) {
                    path = chosen;
                    bool found = false;
                    for (int i = 0; i < endings.size(); i++) if (path.extension() == endings[i]) found = true;
                    if (found) end = true;
                }
                else return false;
            }
            else {
                path = chosen;
                bool found = false;
                for (int i = 0; i < endings.size(); i++) if (path.extension() == endings[i]) found = true;
                if (found) end = true;
            }
        }
        else
        {
            // All this stuff below is to tell you exactly how you messed up above.
            // Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
            switch (CommDlgExtendedError())
            {
            case CDERR_DIALOGFAILURE: std::cout << "CDERR_DIALOGFAILURE\n";   break;
            case CDERR_FINDRESFAILURE: std::cout << "CDERR_FINDRESFAILURE\n";  break;
            case CDERR_INITIALIZATION: std::cout << "CDERR_INITIALIZATION\n";  break;
            case CDERR_LOADRESFAILURE: std::cout << "CDERR_LOADRESFAILURE\n";  break;
            case CDERR_LOADSTRFAILURE: std::cout << "CDERR_LOADSTRFAILURE\n";  break;
            case CDERR_LOCKRESFAILURE: std::cout << "CDERR_LOCKRESFAILURE\n";  break;
            case CDERR_MEMALLOCFAILURE: std::cout << "CDERR_MEMALLOCFAILURE\n"; break;
            case CDERR_MEMLOCKFAILURE: std::cout << "CDERR_MEMLOCKFAILURE\n";  break;
            case CDERR_NOHINSTANCE: std::cout << "CDERR_NOHINSTANCE\n";     break;
            case CDERR_NOHOOK: std::cout << "CDERR_NOHOOK\n";          break;
            case CDERR_NOTEMPLATE: std::cout << "CDERR_NOTEMPLATE\n";      break;
            case CDERR_STRUCTSIZE: std::cout << "CDERR_STRUCTSIZE\n";      break;
            case FNERR_BUFFERTOOSMALL: std::cout << "FNERR_BUFFERTOOSMALL\n";  break;
            case FNERR_INVALIDFILENAME: std::cout << "FNERR_INVALIDFILENAME\n"; break;
            case FNERR_SUBCLASSFAILURE: std::cout << "FNERR_SUBCLASSFAILURE\n"; break;
            default: std::cout << "You cancelled.\n";
            }
            return false;
        }
    }
#else
    std::string command;
    command = "zenity --file-selection --title=\"Select where to save\" --filename=\"";
    command += current_path;
    command += filename_given;
    command += "/\" --save --confirm-overwrite 2> /dev/null";
    bool end = false;
    while (!end) {
        FILE* f = popen(command.c_str(), "r");
        if (!f) {
            std::cout << "ERROR" << std::endl;
            return false;
        }
        std::string file;
        char buf[256];
        while (fgets(buf, sizeof(buf), f)) {
            file += buf;
        }
        if (file.empty())
            return false;
        if (debug) 
            std::cout << "Filename: " << file << std::endl;
        path = file;
        std::stringstream ss(path);
        std::string name = path.string();
        getline(ss, name);
        if (debug) std::cout << "Path: " << path << std::endl;
        if (pclose(f) != 0) std::cout << "Zenity returned non zero, whatever that means..." << std::endl;
        bool found = false;
        for (int i = 0; i < endings.size(); i++) 
            if (path.string().find(endings[i]) != std::string::npos)
                found = true;
        if (found) 
            end = true;
    }
#endif
    return true;
};

const int sht2nbas(const int &type)
{
    const int st2bas[9]{1, 3, 6, 10, 15, 21, 28, 36};
    const int nst2bas[9]{17,15,13,11, 9, 7, 5, 4, 1};
    if (type >= 0)
        return st2bas[type];
    else
        return nst2bas[8 + type];
};

char asciitolower(char in)
{
    if (in <= 'Z' && in >= 'A')
        return in - ('Z' - 'z');
    return in;
}

int vec_sum(const bvec &in)
{
    int sum = 0;
    for (bool val : in)
    {
        sum += val;
    }
    return sum;
}

int vec_sum(const ivec &in)
{
    int sum = 0;
    for (int val : in)
    {
        sum += val;
    }
    return sum;
}

double vec_sum(const vec &in)
{
    double sum = 0.0;
    for (double val : in)
    {
        sum += val;
    }
    return sum;
}

cdouble vec_sum(const cvec &in)
{
    cdouble res = 0.0;
    for (int i = 0; i < in.size(); i++)
        res += in[i];
    return res;
}

double vec_length(const vec &in)
{
    double sum = 0.0;
    for (double val : in)
    {
        sum += val * val;
    }
    return sqrt(sum);
}

void error_check(const bool condition, const std::source_location loc, const std::string &error_message, std::ostream &log_file)
{
    if (!condition)
    {
        log_file << "Error in " << loc.function_name() << " at: " << loc.file_name() << " : " << loc.line() << " " << error_message << std::endl;
        log_file.flush();
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Error in " << loc.function_name() << " at: " << loc.file_name() << " : " << loc.line() << " " << error_message << std::endl;
        exit(-1);
    }
};
void not_implemented(const std::source_location loc, const std::string &error_message, std::ostream &log_file)
{
    log_file << loc.function_name() << " at: " << loc.file_name() << " : " << loc.line() << " " << error_message << " not yet implemented!" << std::endl;
    log_file.flush();
    std::cout.rdbuf(coutbuf); // reset to standard output again
    std::cout << "Error in " << loc.function_name() << " at: " << loc.file_name() << " : " << loc.line() << " " << error_message << " not yet implemented!" << std::endl;
    exit(-1);
};

void sha::sha256_transform(uint32_t state[8], const uint8_t block[64])
{
    uint32_t w[64];
    uint32_t a, b, c, d, e, f, g, h;

    for (int i = 0; i < 16; ++i)
    {
        w[i] = (block[i * 4] << 24) | (block[i * 4 + 1] << 16) |
               (block[i * 4 + 2] << 8) | (block[i * 4 + 3]);
    }

    for (int i = 16; i < 64; ++i)
    {
        w[i] = SIG1(w[i - 2]) + w[i - 7] + SIG0(w[i - 15]) + w[i - 16];
    }

    a = state[0];
    b = state[1];
    c = state[2];
    d = state[3];
    e = state[4];
    f = state[5];
    g = state[6];
    h = state[7];

    for (int i = 0; i < 64; ++i)
    {
        uint32_t temp1 = h + EP1(e) + CH(e, f, g) + k[i] + w[i];
        uint32_t temp2 = EP0(a) + MAJ(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + temp1;
        d = c;
        c = b;
        b = a;
        a = temp1 + temp2;
    }

    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    state[5] += f;
    state[6] += g;
    state[7] += h;
}

// SHA-256 update function
void sha::sha256_update(uint32_t state[8], uint8_t buffer[64], const uint8_t *data, size_t len, uint64_t &bitlen)
{
    for (size_t i = 0; i < len; ++i)
    {
        size_t buf_idx = (bitlen / 8) % 64;
        if (buf_idx >= 64) {
            std::cerr << "Buffer overflow detected in sha256_update!" << std::endl;
            return;
        }
        buffer[buf_idx] = data[i];
        bitlen += 8;
        if (bitlen % 512 == 0)
        {
            sha256_transform(state, buffer);
        }
    }
}

// SHA-256 padding and final hash computation
void sha::sha256_final(uint32_t state[8], uint8_t buffer[64], uint64_t bitlen, uint8_t hash[32])
{
    size_t i = bitlen / 8 % 64;

    if (i >= 64) {
        std::cerr << "Buffer overflow detected in sha256_final (initial)!" << std::endl;
        return;
    }
    buffer[i++] = 0x80;
    if (i > 56)
    {
        while (i < 64)
        {
            if (i >= 64) {
                std::cerr << "Buffer overflow detected in sha256_final (padding)!" << std::endl;
                return;
            }
            buffer[i++] = 0x00;
        }
        sha256_transform(state, buffer);
        i = 0;
    }

    while (i < 56)
    {
        if (i >= 64) {
            std::cerr << "Buffer overflow detected in sha256_final (final padding)!" << std::endl;
            return;
        }
        buffer[i++] = 0x00;
    }

    // memcpy to buffer + 56 is safe as buffer is 64 bytes
    bitlen = custom_bswap_64(bitlen);
    memcpy(buffer + 56, &bitlen, 8);
    sha256_transform(state, buffer);

    for (i = 0; i < 8; ++i)
    {
        constexpr size_t HASH_SIZE = 32;
        const size_t off = static_cast<size_t>(i) * 4;
        if (off + 4 > HASH_SIZE) {
            std::cerr << "Buffer overflow detected in sha256_final (hash write)!" << std::endl;
            return;
        }
        // convert to big-endian value locally and copy 4 bytes
        uint32_t be = custom_bswap_32(state[i]);
        std::memcpy(hash + off, &be, sizeof(be));
    }
}

// Function to calculate SHA-256 hash
std::string sha::sha256(const std::string &input)
{
    uint32_t state[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};

    uint8_t buffer[64] = {0};
    uint8_t hash[32];
    uint64_t bitlen = 0;

    sha256_update(state, buffer, reinterpret_cast<const uint8_t *>(input.c_str()), input.length(), bitlen);
    sha256_final(state, buffer, bitlen, hash);

    std::stringstream ss;
    for (int i = 0; i < 32; ++i)
    {
        ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
    }

    return ss.str();
}

/*
#ifdef _WIN32
// Function to convert a narrow string to a wide string
std::wstring s2ws(const std::string &s)
{
    int len;
    int slength = (int)s.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
    std::wstring r(len, L'\0');
    MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, &r[0], len);
    return r;
}
*/

ProgressBar::~ProgressBar()
{
    progress_ = 100.0f;
    write_progress();
    std::cout << std::endl;
#ifdef _WIN32
    if (taskbarList_)
    {
        taskbarList_->SetProgressState(GetConsoleWindow(), TBPF_NOPROGRESS);
        taskbarList_->Release();
    }
#endif
}

void ProgressBar::write_progress(std::ostream &os)
{
    // No need to write once progress is 100%
    if (progress_ > 100.0f)
        return;

    // Move cursor to the first position on the same line
    // Check if os is a file stream
    if (dynamic_cast<std::filebuf *>(std::cout.rdbuf()))
    {
        os.seekp(linestart); // Is a file stream
    }
    else
    {
        os << "\r" << std::flush; // Is not a file stream
    }

    // Start bar
    os << "[";

    const auto completed = static_cast<size_t>(progress_ * static_cast<float>(bar_width_) / 100.0);
    for (size_t i = 0; i <= completed; ++i)
    {
        os << fill_;
    }

    // End bar
    if (((progress_ < 100.0f) ? progress_ : 100.0f) == 100)
    {
        os << "] 100% " << std::flush;
#ifdef _WIN32
        if (taskbarList_)
        {
            taskbarList_->SetProgressValue(GetConsoleWindow(), 100, 100);
            taskbarList_->SetProgressState(GetConsoleWindow(), TBPF_NOPROGRESS);
        }
#endif
        return;
    }

    os << std::flush;

#ifdef _WIN32
    // Update taskbar progress
    if (taskbarList_)
    {
        taskbarList_->SetProgressValue(GetConsoleWindow(), static_cast<ULONGLONG>(progress_), 100);
    }
#endif
}

#ifdef _WIN32
void ProgressBar::initialize_taskbar_progress()
{
    if (SUCCEEDED(CoInitialize(nullptr)))
    {
        if (SUCCEEDED(CoCreateInstance(CLSID_TaskbarList, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&taskbarList_))))
        {
            taskbarList_->HrInit();
            taskbarList_->SetProgressState(GetConsoleWindow(), TBPF_NORMAL);
        }
    }
}
#endif