#include "convenience.h"
#include "cell.h"
#include "tsc_block.h"
#include "test_functions.h"
#include "properties.h"
#ifdef _WIN32
#include <windows.h>
#endif

std::string help_message()
{
    std::string t = "\n----------------------------------------------------------------------------\n";
    t.append("          These commands and arguments are known by NoSpherA2:\n");
    t.append("----------------------------------------------------------------------------\n\n");
    t.append(":::::::::::::::::::::: Defaults are highlighted by [] ::::::::::::::::::::::\n\n");
    t.append("   -wfn            <FILENAME>.xxx           Read the following wavefunction file.\n");
    t.append("                                            Supported filetypes: .wfn/wfx/ffn; .molden; .xyz; .gbw; .xtb; fch* (UNTESTED!)\n");
    t.append("   -fchk           <FILENAME>.fchk          Write a wavefunction to the given filename\n");
    t.append("   -b              <FILENAME>               Read this basis set\n");
    t.append("   -d              <PATH>                   Path to basis_sets directory with basis_sets in tonto style\n");
    t.append("   -dmin		     <NUMBER>                 Minimum d-spacing to consider for scattering factors (repalaces hkl file)\n");
    t.append("   -ECP            <NUMBER>                 Defines the ECP corrections to be applied to a wavefunction. The given Number defines the ECP type:\n");
    t.append("                                            [1]: def2-ECP\n");
    t.append("   --help/-help/--h                         print this help\n");
    t.append("   -v                                       Turn on Verbose (debug) Mode (Slow and a LOT of output!)\n");
    t.append("   -v2                                      Even more stuff\n");
    t.append("   -mult           <NUMBER>                 Input multiplicity of wavefunction (otherwise attempted to be read from the wfn)\n");
    t.append("   -method         <METHOD NAME>            Can be [RKS] or RHF to distinguish between DFT and HF\n");
    t.append("   -cif            <FILENAME>.cif           CIF to get labels of atoms to use for calculation of scatteriung factors\n");
    t.append("   -IAM                                     Make scattering factors based on Thakkar functions for atoms in CIF\n");
    t.append("   -xyz            <FILENAME>.xyz           Read atom positions from this xyz file for IAM\n");
    t.append("   -hkl            <FILENAME>.hkl           hkl file (ideally merged) to use for calculation of form factors.\n");
    t.append("   -group          <LIST OF INT NUMBERS>    Disorder groups to be read from the CIF for consideration as asym unit atoms (space separated).\n");
    t.append("   -acc            0,1,[2],3,4...           Accuracy of numerical grids used, where the number indicates a pre-defined level. 4 should be considered maximum,\n");
    t.append("                                            anything above will most likely introduce numberical error and is just implemented for testing purposes.");
    t.append("   -gbw2wfn                                 Only reads wavefucntion from .gbw specified by -wfn and prints it into .wfn format.\n");
    t.append("   -tscb           <FILENAME>.tscb          Convert binary tsc file to bigger, less accurate human-readable form.\n");
    t.append("   -twin           -1 0 0 0 -1 0 0 0 -1     3x3 floating-point-matrix in the form -1 0 0 0 -1 0 0 0 -1 which contains the twin matrix to use.\n");
    t.append("                                            If there is more than a single twin law to be used, use the twin command multiple times.\n");
    t.append("   -merge          <List of .tsc files>     Names/Paths to .tsc/.tscb files to be merged.\n");
    t.append("   -merge_nocheck  <List of .tsc files>     Names/Paths to .tsc/.tscb files to be merged. They need to have identical hkl values.\n");
    t.append("   -mtc            <List of .wfns + parts>  Performs calculation for a list of wavefunctions (=Multi-Tsc-Calc), where asymmetric unit is.\n");
    t.append("                                            taken from given CIF. Also disorder groups are required per file as comma separated list\n");
    t.append("                                            without spaces.\n   Typical use Examples:\n");
    t.append("   -SALTED         <Path to Model folder>   Uses a provided SALTED-ML Model to predict the electron densitie of a xyz-file\n");
    t.append("   -cmtc           <List of .wfns + parts>  Performs calculation for a list of wavefunctions AND CIFs (=CIF-based-multi-Tsc-Calc), where asymmetric unit is defined by each CIF that matches a wfn.\n");
    t.append("      Normal:       NoSpherA2.exe -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7\n");
    t.append("      thakkar-tsc:  NoSpherA2.exe -cif A.cif -hkl A.hkl -xyz A.xyz -acc 1 -cpus 7 -IAM\n");
    t.append("      Disorder:     NoSpherA2.exe -cif A.cif -hkl A.hkl -acc 1 -cpus 7 -mtc 1.wfn 0,1 2.wfn 0,2 3.wfn 0,3\n");
    t.append("      fragHAR:      NoSpherA2.exe -cif A.cif -hkl A.hkl -acc 1 -cpus 7 -cmtc 1.wfn 1.cif 0 2.wfn 2.cif 0 3_1.wfn 3_1.cif 0,1 3_2.wfn 3_2.cif 0,2\n");
    t.append("      merging tscs: NoSpherA2.exe -merge A.tsc B.tsc C.tsc\n");
    t.append("      merge tsc(2): NoSpherA2.exe -merge_nocheck A.tsc B.tsc C.tsc  (MAKE SURE THEY HAVE IDENTICAL HKL INIDCES!!)\n");
    t.append("      convert tsc:  NoSpherA2.exe -tscb A.tscb\n");
    t.append("      convert gbw:  NoSpherA2.exe -gbw2wfn -wfn A.gbw\n");
    t.append("      twin law:     NoSpherA2.exe -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7 -twin -1 0 0 0 -1 0 0 0 -1\n");
    return t;
}
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
        t.append("List of contributors of pieces of code or funcitonality:\n");
        t.append("      Florian Kleemiss,\n");
        t.append("      Emmanuel Hupf,\n");
        t.append("      Alessandro Genoni,\n");
        t.append("      Lukas Seifert,\n");
        t.append("      and many more in communications or by feedback!\n");
#if has_RAS
        t.append("NoSpherA2 uses Rascaline, Metatensor, OpenBLAS and the HDF5 library.\n");
        t.append("The used packages are published under BSD-3 clause License.\n");
        t.append("Please see, respectively:\n");
        t.append("   https://github.com/Luthaf/rascaline\n");
        t.append("   https://github.com/lab-cosmo/metatensor\n");
        t.append("   https://github.com/OpenMathLib/OpenBLAS\n");
        t.append("   https://github.com/HDFGroup/hdf5\n");
#endif
        t.append("NoSpherA2 was published at  : Kleemiss et al. Chem. Sci., 2021, 12, 1675 - 1692.\n");
        t.append("Slater IAM was published at : Kleemiss et al. J. Appl. Cryst 2024, 57, 161 - 174.\n");
    }
    return t;
}

std::string build_date()
{
    return ("This Executable was built on: " + std::string(__DATE__) + " " + std::string(__TIME__) + "\n");
}

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
        // cout << "resize: i="<<i<<" y=" << y << endl;
        mBase_values[i] = y;
    }
    mStepwidth = MPI2 / size;
}

double cosinus_annaeherung::calculate_error_at(double x) const
{
    return cos(x) - get(x);
}
*/
void copy_file(std::string &from, std::string &to)
{
    std::ifstream source(from.c_str(), std::ios::binary);
    std::ofstream dest(to.c_str(), std::ios::binary);

    dest << source.rdbuf();

    source.close();
    dest.close();
};

//---------------------------Configuration files ---------------------------------------------------

std::string get_home_path(void)
{
#ifdef _WIN32
    std::string temp1 = getenv("HOMEDRIVE");
    std::string temp2 = getenv("HOMEPATH");
    temp1.append(temp2);
    return temp1;
#else
    std::string home = getenv("HOME");
    return home;
#endif
}

void join_path(std::string &s1, std::string &s2)
{
#ifdef _WIN32
    s1.append("\\");
#else
    if (s1.substr(s1.length() - 1) != "/")
        s1.append("/");
#endif
    s1.append(s2);
}

void join_path(std::string &s1, std::initializer_list<std::string> s2)
{
    std::string separator;
#ifdef _WIN32
    separator = "\\";
#else
    separator = "/";
#endif
    // Ensure that the initial path has a trailing separator if needed
    if (!s1.empty() && s1.back() != '/' && s1.back() != '\\')
    {
        s1.append(separator);
    }
    // Iterate over each segment in the initializer list
    for (auto it = s2.begin(); it != s2.end(); ++it)
    {
        s1.append(*it);
        // Append the separator if it's not the last segment
        if (std::next(it) != s2.end() && it->back() != '/' && it->back() != '\\')
        {
            s1.append(separator);
        }
    }
}

std::string get_filename_from_path(const std::string &input)
{
#ifdef _WIN32
    return input.substr(input.rfind("\\") + 1);
#else
    return input.substr(input.rfind("/") + 1);
#endif
}

std::string get_foldername_from_path(const std::string &input)
{
#ifdef _WIN32
    return input.substr(0, input.rfind("\\") + 1);
#else
    return input.substr(0, input.rfind("/") + 1);
#endif
}

std::string get_basename_without_ending(const std::string &input)
{
    return input.substr(0, input.rfind("."));
}

std::string get_ending_from_filename(const std::string &input)
{
    return input.substr(input.rfind(".") + 1);
}

void write_template_confi()
{
    using namespace std;
    string line;
    string programs = get_home_path();
    string filename = ".cuQCT.conf";
    join_path(programs, filename);
    if (exists(programs))
    {
        cout << "File already exists! Aborting!" << endl;
        return;
    }
    ofstream conf(programs.c_str());
#ifdef _WIN32
    conf << "gaussian=\"D:\\g09\\g09\\\"" << endl;
    conf << "turbomole=\"D:\\turbomole\\dscf7.1\\\"" << endl;
    conf << "basis=\"D:\\tonto\\basis_sets\\\"" << endl;
#else
    conf << "gaussian=\"/usr/local/g09/g09\"" << endl;
    conf << "turbomole=\"/usr/local/bin/dscf7.1\"" << endl;
    conf << "basis=\"/basis_sets/\"" << endl;
#endif
    conf << "cpu=4" << endl;
    conf << "mem=4.0" << endl;
    conf << "rho=1" << endl;
    conf << "rdg=1" << endl;
    conf << "eli=0" << endl;
    conf << "elf=0" << endl;
    conf << "lap=0" << endl;
    conf << "esp=0" << endl;
    conf << "efv=0" << endl;
    conf << "def=0" << endl;
    conf << "hir=0" << endl;
    conf.flush();
    conf.close();
#ifdef _WIN32
    //	const wchar_t* fileLPCWSTR = programs.c_str();
    //	wstring stemp = wstring(programs.begin(), programs.end());
    //	int attr = GetFileAttributes(stemp.c_str());
    //	if ((attr & FILE_ATTRIBUTE_HIDDEN) == 0) {
    //		SetFileAttributes(stemp.c_str(), attr | FILE_ATTRIBUTE_HIDDEN);
    //	}
#endif
    return;
};

options::options(const int accuracy,
                 const int threads,
                 const int pbc,
                 const double resolution,
                 const double radius,
                 const bool becke,
                 const bool electron_diffraction,
                 const bool ECP,
                 const bool set_ECPs,
                 const int ECP_mode,
                 const bool calc,
                 const bool eli,
                 const bool esp,
                 const bool elf,
                 const bool lap,
                 const bool rdg,
                 const bool hdef,
                 const bool def,
                 const bool fract,
                 const bool hirsh,
                 const bool s_rho,
                 const bool SALTED,
                 const bool Olex2_1_3_switch,
                 const bool iam_switch,
                 const bool read_k_pts,
                 const bool save_k_pts,
                 const bool combined_tsc_calc,
                 const bool binary_tsc,
                 const bool cif_based_combined_tsc_calc,
                 const bool no_date,
                 const bool gbw2wfn,
                 const bool old_tsc,
                 const bool thakkar_d_plot,
                 const double sfac_scan,
                 const double sfac_diffuse,
                 const double dmin,
                 const int hirsh_number,
                 const ivec &MOs,
                 const ivec2 &_groups,
                 const vec2 &twin_law,
                 const ivec2 &combined_tsc_groups,
                 const bool all_mos,
                 const bool test,
                 const std::string &wfn,
                 const std::string &fchk,
                 const std::string &basis_set,
                 const std::string &hkl,
                 const std::string &cif,
                 const std::string &method,
                 const std::string &xyz_file,
                 const std::string &coef_file,
                 const std::string &fract_name,
                 const svec &combined_tsc_calc_files,
                 const svec &combined_tsc_calc_cifs,
                 const std::string &wavename,
                 const std::string &gaussian_path,
                 const std::string &turbomole_path,
                 const std::string &basis_set_path,
                 const std::string &SALTED_DIR,
                 const std::string &SALTED_DFBASIS,
                 const svec &arguments,
                 const svec &combine_mo,
                 const svec &Cations,
                 const svec &Anions,
                 const ivec &cmo1,
                 const ivec &cmo2,
                 const ivec &ECP_nrs,
                 const ivec &ECP_elcounts,
                 const double mem,
                 const unsigned int mult,
                 const bool debug,
                 const hkl_list &m_hkl_list,
                 std::ostream &log_file)
    : accuracy(accuracy), threads(threads), pbc(pbc),
      resolution(resolution), radius(radius), becke(becke),
      electron_diffraction(electron_diffraction), ECP(ECP),
      set_ECPs(set_ECPs), ECP_mode(ECP_mode), calc(calc),
      eli(eli), esp(esp), elf(elf), lap(lap), rdg(rdg),
      hdef(hdef), def(def), fract(fract), hirsh(hirsh),
      s_rho(s_rho), SALTED(SALTED), SALTED_DIR(SALTED_DIR), SALTED_DFBASIS(SALTED_DFBASIS),
      Olex2_1_3_switch(Olex2_1_3_switch), iam_switch(iam_switch),
      read_k_pts(read_k_pts), save_k_pts(save_k_pts),
      combined_tsc_calc(combined_tsc_calc), binary_tsc(binary_tsc),
      cif_based_combined_tsc_calc(cif_based_combined_tsc_calc),
      no_date(no_date), gbw2wfn(gbw2wfn), old_tsc(old_tsc),
      thakkar_d_plot(thakkar_d_plot), d_sfac_scan(sfac_scan),
      sfac_diffuse(sfac_diffuse), dmin(dmin), hirsh_number(hirsh_number),
      MOs(MOs), groups(_groups), twin_law(twin_law),
      combined_tsc_groups(combined_tsc_groups), all_mos(all_mos),
      test(test), wfn(wfn), fchk(fchk), basis_set(basis_set),
      hkl(hkl), cif(cif), method(method), xyz_file(xyz_file),
      coef_file(coef_file), fract_name(fract_name),
      combined_tsc_calc_files(combined_tsc_calc_files),
      combined_tsc_calc_cifs(combined_tsc_calc_cifs),
      wavename(wavename), gaussian_path(gaussian_path),
      turbomole_path(turbomole_path), basis_set_path(basis_set_path),
      arguments(arguments), combine_mo(combine_mo), Cations(Cations),
      Anions(Anions), cmo1(cmo1), cmo2(cmo2), ECP_nrs(ECP_nrs),
      ECP_elcounts(ECP_elcounts), mem(mem), mult(mult),
      debug(debug), m_hkl_list(m_hkl_list), log_file(log_file)
{
    groups.resize(1);
};

int program_confi(std::string &gaussian_path, std::string &turbomole_path, std::string &basis, int &ncpus, double &mem, bool debug, bool expert, unsigned int counter)
{
    using namespace std;
    counter++;
    if (counter == 3)
    {
        cout << "Too many iterations of tries to read config file, better abort..." << endl;
        return -1;
    }
    string programs = get_home_path();
    string filename = ".cuQCT.conf";
    join_path(programs, filename);
    ifstream conf(programs.c_str());
    if (debug)
        cout << programs << endl;
    string line;
    if (conf.good())
    {
        if (debug)
            cout << "File is valid, continuing..." << endl;
    }
    else
    {
        if (expert)
        {
            cout << "couldn't open or find .cuQCT.conf, in your home folder: " << programs << ", writing a template for you!" << endl;
            write_template_confi();
            if (program_confi(gaussian_path, turbomole_path, basis, ncpus, mem, debug, expert, counter) != 1)
                return -1;
            cout << "Wrote a template for you, read default values!" << endl;
            return 0;
        }
        else
        {
            write_template_confi();
            program_confi(gaussian_path, turbomole_path, basis, ncpus, mem, debug);
            return 0;
        }
    }
    conf.seekg(0);
    getline(conf, line);
    size_t length;
    char *tempchar = new char[200];
    int run = 0;
    while (!conf.eof())
    {
        switch (run)
        {
        case 0:
            length = line.copy(tempchar, line.size() - 11, 10);
            tempchar[length] = '\0';
            gaussian_path = tempchar;
            break;
        case 1:
            length = line.copy(tempchar, line.size() - 12, 11);
            tempchar[length] = '\0';
            turbomole_path = tempchar;
            break;
        case 2:
            length = line.copy(tempchar, line.size() - 8, 7);
            tempchar[length] = '\0';
            basis = tempchar;
            break;
        case 3:
            length = line.copy(tempchar, line.size() - 3, 4);
            tempchar[length] = '\0';
            ncpus = stoi(tempchar);
            break;
        case 4:
            length = line.copy(tempchar, line.size() - 3, 4);
            tempchar[length] = '\0';
            mem = stod(tempchar);
            break;
        default:
            if (debug)
                cout << "found everything i was looking for, if you miss something check the switch" << endl;
            break;
        }
        if (debug)
            cout << run << ". line: " << tempchar << endl;
        run++;
        getline(conf, line);
    }
    return 1;
};

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

int filetype_identifier(std::string &file, bool debug)
{
    /*
    List of filetypes and correpsonding values:
                    -1: unreadable keyword
    -i/o *.wfn 		2: wfn
    -i/o *.ffn 		4: ffn
    -i *.out 		1: crystal output
    -c/o *.cub(e) 	3: cube file
    -g/o *.grd      6: XDGraph grid file
    -o *.(F)fc(C)hk 5: fchk
    */
    using namespace std;
    if (debug)
    {
        std::cout << "Testing WFN:  " << file.find(".wfn") << std::endl
                  << "Testing out:  " << file.find(".out") << std::endl
                  << "Testing FFN:  " << file.find(".ffn") << std::endl
                  << "Testing CUB:  " << file.find(".cub") << std::endl
                  << "Testing CUBE: " << file.find(".cube") << std::endl
                  << "Testing Grid: " << file.find(".grd") << std::endl
                  << "Testing fchk: " << file.find(".fchk") << std::endl
                  << "Testing FChk: " << file.find(".FChk") << std::endl
                  << "Testing Fchk: " << file.find(".Fchk") << std::endl;
        std::cout << "string::npos: " << std::string::npos << std::endl;
    }
    int temp_type = 0;
    size_t found, temp;
    temp = 0;
    if (debug)
        std::cout << "Temp before any checks: " << temp << std::endl;
    svec types{".out", ".wfn", ".ffn", ".cub", ".cube", ".grd", ".fchk", ".Fchk", ".FChk"};
    if (file.find(".wfn") != std::string::npos)
    {
        if (debug)
            std::cout << "Checking for"
                      << ".wfn" << std::endl;
        temp_type = 2;
        found = file.rfind(".wfn");
        if (debug)
            std::cout << "Found: " << found << std::endl;
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != std::string::npos)
                temp = file.rfind(types[i]);
        if (debug)
            cout << "Temp: " << temp << endl;
        if (temp == found)
            return temp_type;
        else
        {
            if (debug)
                cout << "Moving on!" << endl;
        }
    }
    if (file.find(".out") != string::npos)
    {
        if (debug)
            cout << "Checking for"
                 << ".out" << endl;
        temp_type = 1;
        found = file.rfind(".out");
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
                temp = file.rfind(types[i]);
        if (temp == found)
            return temp_type;
    }
    if (file.find(".ffn") != string::npos)
    {
        if (debug)
            cout << "Checking for"
                 << ".ffn" << endl;
        temp_type = 4;
        found = file.rfind(".ffn");
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
                temp = file.rfind(types[i]);
        if (temp == found)
            return temp_type;
    }
    if (file.find(".cub") != string::npos)
    {
        if (debug)
            cout << "Checking for"
                 << ".cub" << endl;
        temp_type = 3;
        found = file.rfind(".cub");
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
                temp = file.rfind(types[i]);
        if (temp == found)
            return temp_type;
        else
        {
            if (debug)
                cout << "Moving on!" << endl;
        }
    }
    if (file.find(".cube") != string::npos)
    {
        temp_type = 3;
        found = file.rfind(".cube");
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
                temp = file.rfind(types[i]);
        if (temp == found)
            return temp_type;
    }
    if (file.find(".grd") != string::npos)
    {
        temp_type = 6;
        found = file.rfind(".grd");
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
                temp = file.rfind(types[i]);
        if (temp == found)
            return temp_type;
    }
    if (file.find(".fchk") != string::npos)
    {
        temp_type = 5;
        found = file.rfind(".fchk");
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
                temp = file.rfind(types[i]);
        if (temp == found)
            return temp_type;
    }
    if (file.find(".FChk") != string::npos)
    {
        temp_type = 5;
        found = file.rfind(".FChk");
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
                temp = file.rfind(types[i]);
        if (temp == found)
            return temp_type;
    }
    if (file.find(".Fchk") != string::npos)
    {
        temp_type = 5;
        found = file.rfind(".Fchk");
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
                temp = file.rfind(types[i]);
        if (temp == found)
            return temp_type;
    }
    return -1;
}

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

primitive::primitive(int c, int t, double e, double coef) : center(c), type(t), exp(e), coefficient(coef)
{
    norm_const = pow(
        pow(2, 7 + 4 * type) * pow(exp, 3 + 2 * type) / constants::PI / pow(doublefactorial(2 * type + 1), 2),
        0.25);
};

tonto_primitive::tonto_primitive(int c, int t, double e, double coef) : center(c), type(t), exp(e), coefficient(coef)
{
    norm_const = pow(constants::PI, -0.75) * pow(2.0, type + 0.75) * pow(exp, type * 0.5 + 0.75) / sqrt(doublefactorial(type));
};

/*bool open_file_dialog(string &path, bool debug, vector <string> filter){
    char pwd[1024];
    if(GetCurrentDir( pwd, 1024)==NULL) return false;
    string current_path(pwd);
#ifdef _WIN32
    char filename[ 1024 ];

    OPENFILENAMEA ofn;
            ZeroMemory( &filename, sizeof( filename ) );
            ZeroMemory( &ofn,      sizeof( ofn ) );
            ofn.lStructSize  = sizeof( ofn );
            ofn.hwndOwner    = NULL;  // If you have a window to center over, put its HANDLE here
            ofn.lpstrFilter  = "wfn Files\0*.wfn\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0Any File\0*.*\0";
            ofn.lpstrFile    = filename;
            ofn.nMaxFile     = 1024;
            ofn.lpstrTitle   = "Select a File";
            ofn.Flags        = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST;

    if (GetOpenFileNameA( &ofn )){
        if(debug) cout << "You chose the file \"" << filename << "\"\n";
        if(exists(filename)){
            path=filename;
            return true;
        }
    }
    else
    {
        // All this stuff below is to tell you exactly how you messed up above.
        // Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
        switch (CommDlgExtendedError())
        {
            case CDERR_DIALOGFAILURE   : cout << "CDERR_DIALOGFAILURE\n";   break;
            case CDERR_FINDRESFAILURE  : cout << "CDERR_FINDRESFAILURE\n";  break;
            case CDERR_INITIALIZATION  : cout << "CDERR_INITIALIZATION\n";  break;
            case CDERR_LOADRESFAILURE  : cout << "CDERR_LOADRESFAILURE\n";  break;
            case CDERR_LOADSTRFAILURE  : cout << "CDERR_LOADSTRFAILURE\n";  break;
            case CDERR_LOCKRESFAILURE  : cout << "CDERR_LOCKRESFAILURE\n";  break;
            case CDERR_MEMALLOCFAILURE : cout << "CDERR_MEMALLOCFAILURE\n"; break;
            case CDERR_MEMLOCKFAILURE  : cout << "CDERR_MEMLOCKFAILURE\n";  break;
            case CDERR_NOHINSTANCE     : cout << "CDERR_NOHINSTANCE\n";     break;
            case CDERR_NOHOOK          : cout << "CDERR_NOHOOK\n";          break;
            case CDERR_NOTEMPLATE      : cout << "CDERR_NOTEMPLATE\n";      break;
            case CDERR_STRUCTSIZE      : cout << "CDERR_STRUCTSIZE\n";      break;
            case FNERR_BUFFERTOOSMALL  : cout << "FNERR_BUFFERTOOSMALL\n";  break;
            case FNERR_INVALIDFILENAME : cout << "FNERR_INVALIDFILENAME\n"; break;
            case FNERR_SUBCLASSFAILURE : cout << "FNERR_SUBCLASSFAILURE\n"; break;
            default                    : cout << "You cancelled.\n";
        }
    }
    return false;
#else
    char file[1024];
    string command;
    command="zenity --file-selection --title=\"Select a file to load\" --filename=\"";
    command += current_path;
    command += "/\"";
    for(int i=0; i<filter.size(); i++){
        command += " --file-filter=\"";
        command += filter[i];
        command += "\" ";
    }
    command += " 2> /dev/null";
    FILE *f = popen(command.c_str(), "r");
    if(!f){
        cout << "ERROR" << endl;
        return false;
    }
    if(fgets(file, 1024, f)==NULL) return false;
    if (debug) cout << "Filename: " << file << endl;
    path=file;
    stringstream ss(path);
    getline(ss, path);
    if(pclose(f)!=0) cout << "Zenity returned non zero, whatever that means..." << endl;
#endif
    return true;
};

bool save_file_dialog(string &path, bool debug, const vector <string> &endings, const string &filename_given){
    char pwd[1024];
    if(GetCurrentDir( pwd, 1024)==NULL) return false;
    string current_path(pwd);
#ifdef _WIN32
    char filename[ 1024 ];

    OPENFILENAMEA sfn;
     ZeroMemory( &filename, sizeof( filename ) );
     ZeroMemory( &sfn,      sizeof( sfn ) );
     sfn.lStructSize  = sizeof( sfn );
     sfn.hwndOwner    = NULL;  // If you have a window to center over, put its HANDLE here
     sfn.lpstrFilter  = "wfn Files\0*.wfn\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0Any File\0*.*\0";
     sfn.lpstrFile    = filename;
     sfn.nMaxFile     = 1024;
     sfn.lpstrTitle   = "Select a File, yo!";
     sfn.Flags        = OFN_DONTADDTORECENT;
    bool end=false;
    while(!end){
        if (GetSaveFileNameA( &sfn )){
            if(debug) cout << "You chose the file \"" << filename << "\"\n";
            if(exists(filename)){
                cout << filename << " exists, do you want to overwrite it?";
                if(yesno()){
                    path=filename;
                    bool found=false;
                    for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
                    if(found) end=true;
                }
                else return false;
            }
            else{
                path=filename;
                bool found=false;
                for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
                if(found) end=true;
            }
        }
        else
        {
            // All this stuff below is to tell you exactly how you messed up above.
            // Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
            switch (CommDlgExtendedError())
            {
                case CDERR_DIALOGFAILURE   : cout << "CDERR_DIALOGFAILURE\n";   break;
                case CDERR_FINDRESFAILURE  : cout << "CDERR_FINDRESFAILURE\n";  break;
                case CDERR_INITIALIZATION  : cout << "CDERR_INITIALIZATION\n";  break;
                case CDERR_LOADRESFAILURE  : cout << "CDERR_LOADRESFAILURE\n";  break;
                case CDERR_LOADSTRFAILURE  : cout << "CDERR_LOADSTRFAILURE\n";  break;
                case CDERR_LOCKRESFAILURE  : cout << "CDERR_LOCKRESFAILURE\n";  break;
                case CDERR_MEMALLOCFAILURE : cout << "CDERR_MEMALLOCFAILURE\n"; break;
                case CDERR_MEMLOCKFAILURE  : cout << "CDERR_MEMLOCKFAILURE\n";  break;
                case CDERR_NOHINSTANCE     : cout << "CDERR_NOHINSTANCE\n";     break;
                case CDERR_NOHOOK          : cout << "CDERR_NOHOOK\n";          break;
                case CDERR_NOTEMPLATE      : cout << "CDERR_NOTEMPLATE\n";      break;
                case CDERR_STRUCTSIZE      : cout << "CDERR_STRUCTSIZE\n";      break;
                case FNERR_BUFFERTOOSMALL  : cout << "FNERR_BUFFERTOOSMALL\n";  break;
                case FNERR_INVALIDFILENAME : cout << "FNERR_INVALIDFILENAME\n"; break;
                case FNERR_SUBCLASSFAILURE : cout << "FNERR_SUBCLASSFAILURE\n"; break;
                default                    : cout << "You cancelled.\n";
            }
            return false;
        }
    }
    return true;
#else
    char file[1024];
    string command;
    command="zenity --file-selection --title=\"Select where to save\" --filename=\"";
    command += current_path;
    command += filename_given;
    command += "/\" --save --confirm-overwrite 2> /dev/null";
    bool end=false;
    while(!end){
        FILE *f = popen(command.c_str(), "r");
        if(!f){
            cout << "ERROR" << endl;
            return false;
        }
        if(fgets(file, 1024, f)==NULL) return false;
        if (debug) cout << "Filename: " << file << endl;
        path=file;
        stringstream ss(path);
        getline(ss, path);
        if(debug) cout << "Path: " << path << endl;
        if(pclose(f)!=0) cout << "Zenity returned non zero, whatever that means..." << endl;
        bool found=false;
        for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
        if(found) end=true;
    }
#endif
    return true;
};

bool save_file_dialog(string &path, bool debug, const vector <string> &endings){
    char pwd[1024];
    if(GetCurrentDir( pwd, 1024)==NULL) return false;
    string current_path(pwd);
#ifdef _WIN32
    char filename[ 1024 ];

    OPENFILENAMEA sfn;
     ZeroMemory( &filename, sizeof( filename ) );
     ZeroMemory( &sfn,      sizeof( sfn ) );
     sfn.lStructSize  = sizeof( sfn );
     sfn.hwndOwner    = NULL;  // If you have a window to center over, put its HANDLE here
     sfn.lpstrFilter  = "wfn Files\0*.wfn\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0Any File\0*.*\0";
     sfn.lpstrFile    = filename;
     sfn.nMaxFile     = 1024;
     sfn.lpstrTitle   = "Select a File, yo!";
     sfn.Flags        = OFN_DONTADDTORECENT;
    bool end=false;
    while(!end){
        if (GetSaveFileNameA( &sfn )){
            if(debug) cout << "You chose the file \"" << filename << "\"\n";
            if(exists(filename)){
                cout << filename << " exists, do you want to overwrite it?";
                if(yesno()){
                    path=filename;
                    bool found=false;
                    for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
                    if(found) end=true;
                }
                else return false;
            }
            else{
                path=filename;
                bool found=false;
                for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
                if(found) end=true;
            }
        }
        else
        {
            // All this stuff below is to tell you exactly how you messed up above.
            // Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
            switch (CommDlgExtendedError())
            {
                case CDERR_DIALOGFAILURE   : cout << "CDERR_DIALOGFAILURE\n";   break;
                case CDERR_FINDRESFAILURE  : cout << "CDERR_FINDRESFAILURE\n";  break;
                case CDERR_INITIALIZATION  : cout << "CDERR_INITIALIZATION\n";  break;
                case CDERR_LOADRESFAILURE  : cout << "CDERR_LOADRESFAILURE\n";  break;
                case CDERR_LOADSTRFAILURE  : cout << "CDERR_LOADSTRFAILURE\n";  break;
                case CDERR_LOCKRESFAILURE  : cout << "CDERR_LOCKRESFAILURE\n";  break;
                case CDERR_MEMALLOCFAILURE : cout << "CDERR_MEMALLOCFAILURE\n"; break;
                case CDERR_MEMLOCKFAILURE  : cout << "CDERR_MEMLOCKFAILURE\n";  break;
                case CDERR_NOHINSTANCE     : cout << "CDERR_NOHINSTANCE\n";     break;
                case CDERR_NOHOOK          : cout << "CDERR_NOHOOK\n";          break;
                case CDERR_NOTEMPLATE      : cout << "CDERR_NOTEMPLATE\n";      break;
                case CDERR_STRUCTSIZE      : cout << "CDERR_STRUCTSIZE\n";      break;
                case FNERR_BUFFERTOOSMALL  : cout << "FNERR_BUFFERTOOSMALL\n";  break;
                case FNERR_INVALIDFILENAME : cout << "FNERR_INVALIDFILENAME\n"; break;
                case FNERR_SUBCLASSFAILURE : cout << "FNERR_SUBCLASSFAILURE\n"; break;
                default                    : cout << "You cancelled.\n";
            }
            return false;
        }
    }
    return true;
#else
    char file[1024];
    string command;
    command="zenity --file-selection --title=\"Select where to save\" --filename=\"";
    command += current_path;
    command += "/\" --save --confirm-overwrite 2> /dev/null";
    bool end=false;
    while(!end){
        FILE *f = popen(command.c_str(), "r");
        if(!f){
            cout << "ERROR" << endl;
            return false;
        }
        if(fgets(file, 1024, f)==NULL) return false;
        if (debug) cout << "Filename: " << file << endl;
        path=file;
        stringstream ss(path);
        getline(ss, path);
        if(debug) cout << "Path: " << path << endl;
        if(pclose(f)!=0) cout << "Zenity returned non zero, whatever that means..." << endl;
        bool found=false;
        for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
        if(found) end=true;
    }
#endif
    return true;
};
*/

void select_cubes(std::vector<std::vector<unsigned int>> &selection, std::vector<WFN> &wavy, unsigned int nr_of_cubes, bool wfnonly, bool debug)
{
    // asks which wfn to use, if wfnonly is set or whcih cubes up to nr of cubes to use
    // Returns values in selection[0][i] for iths selection of wavefunction and
    //  selection[1][i] for iths selection of cube
    using namespace std;
    cout << "Which of the following cubes to use? Need to select " << nr_of_cubes << " file";
    if (nr_of_cubes > 1)
        cout << "s in total." << endl;
    else
        cout << "." << endl;
    cout << endl
         << endl;
    for (int w = 0; w < wavy.size(); w++)
    {
        stringstream stream;
        cout << "_____________________________________________________________" << endl;
        cout << "WFN ";
        stream << setw(2) << w;
        cout << stream.str() << ") " << get_basename_without_ending(wavy[w].get_path()) << endl;
        stream.str("");
        for (int c = 0; c < wavy[w].cub.size(); c++)
        {
            if (c == 0)
                cout << "        |" << endl
                     << "Cube    |" << endl;
            else
                cout << "        |" << endl;
            if (!wfnonly)
            {
                cout << setw(2) << w;
                cout << ".";
                cout << setw(2) << c;
            }
            else
                cout << "     ";
            cout << "   |_ " << get_basename_without_ending(wavy[w].cub[c].path);
            if (!exists(wavy[w].cub[c].path))
                cout << " (MEM ONLY)";
            cout << endl;
        }
        cout << "_____________________________________________________________" << endl
             << endl
             << endl;
    }
    // bool happy = false;
    unsigned int selected_cubes = 0;
    do
    {
        cout << "Select " << selected_cubes + 1 << ". ";
        if (wfnonly)
            cout << "WFN ";
        else
            cout << "cube ";
        cout << "please: ";
        string input;
        cin >> input;
        if (!wfnonly)
        {
            if (input.find('.') == string::npos)
            {
                cout << "no . found in input!" << endl;
                continue;
            }
        }
        else
        {
            if (input.find('.') == string::npos)
                cout << "Ignoring the .!" << endl;
            unsigned int nr_wave = fromString<unsigned int>(input);
            if (nr_wave < 0 || nr_wave >= wavy.size())
            {
                cout << "Invalid choice!" << endl;
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
            cout << "input: " << input << endl;
            cout << "with . found at: " << input.find('.') << endl;
            cout << "substr1: " << input.substr(0, input.find('.')) << endl;
            cout << "substr2: " << input.substr(input.find('.') + 1) << endl;
        }
        string wave(input.substr(0, input.find('.')));
        string cube(input.substr(input.find('.') + 1));
        unsigned int nr_wave = fromString<unsigned int>(wave);
        unsigned int nr_cube = fromString<unsigned int>(cube);
        if (debug)
            cout << "Translated: " << nr_wave << " " << nr_cube << endl;
        if (nr_wave < 0 || nr_wave >= wavy.size() || nr_cube < 0 || nr_cube >= wavy[nr_wave].cub.size())
        {
            cout << "Invalid choice!" << endl;
            continue;
        }
        selection[0][selected_cubes] = nr_wave;
        selection[1][selected_cubes] = nr_cube;
        selected_cubes++;
        if (selected_cubes == nr_of_cubes)
        {
            if (debug)
                cout << "Going to return!" << endl;
            return;
        }
    } while (true);
};

bool unsaved_files(std::vector<WFN> &wavy)
{
    for (int w = 0; w < wavy.size(); w++)
        for (int c = 0; c < wavy[w].cub.size(); c++)
            if (!exists(wavy[w].cub[c].path))
                return true;
    return false;
};

void progress_bar::write(double fraction)
{
    // clamp fraction to valid range [0,1]
    if (fraction < 0)
        fraction = 0;
    else if (fraction > 1)
        fraction = 1;

    auto width = bar_width - message.size();
    auto offset = bar_width - static_cast<unsigned>(width * round(fraction / precision) * precision);

    os << '\r' << message;
    os.write(full_bar.data() + offset, width);
    os << " [" << std::setw(3) << static_cast<int>(100 * round(fraction / precision) * precision) << "%] " << std::flush;
}

void readxyzMinMax_fromWFN(
    WFN &wavy,
    double *CoordMinMax,
    int *NbSteps,
    double Radius,
    double Increments,
    bool no_bohr)
{
    vec2 PosAtoms;
    PosAtoms.resize(wavy.get_ncen());
    for (int i = 0; i < wavy.get_ncen(); i++)
        PosAtoms[i].resize(3);
    bool bohrang = true;
    if (!no_bohr)
        bohrang = !check_bohr(wavy, false);
    for (int j = 0; j < wavy.get_ncen(); j++)
    {
        PosAtoms[j][0] = wavy.atoms[j].x;
        PosAtoms[j][1] = wavy.atoms[j].y;
        PosAtoms[j][2] = wavy.atoms[j].z;
        if (!bohrang)
        {
            for (int i = 0; i < 3; i++)
                PosAtoms[j][i] = constants::ang2bohr(PosAtoms[j][i]);
        }
        if (j == 0)
        {
            CoordMinMax[0] = PosAtoms[j][0];
            CoordMinMax[3] = PosAtoms[j][0];
            CoordMinMax[1] = PosAtoms[j][1];
            CoordMinMax[4] = PosAtoms[j][1];
            CoordMinMax[2] = PosAtoms[j][2];
            CoordMinMax[5] = PosAtoms[j][2];
        }
        else
        {
            if (CoordMinMax[0] > PosAtoms[j][0])
                CoordMinMax[0] = PosAtoms[j][0];
            if (CoordMinMax[3] < PosAtoms[j][0])
                CoordMinMax[3] = PosAtoms[j][0];

            if (CoordMinMax[1] > PosAtoms[j][1])
                CoordMinMax[1] = PosAtoms[j][1];
            if (CoordMinMax[4] < PosAtoms[j][1])
                CoordMinMax[4] = PosAtoms[j][1];

            if (CoordMinMax[2] > PosAtoms[j][2])
                CoordMinMax[2] = PosAtoms[j][2];
            if (CoordMinMax[5] < PosAtoms[j][2])
                CoordMinMax[5] = PosAtoms[j][2];
        }
    }
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
    std::string cif,
    double *CoordMinMax,
    int *NbSteps,
    vec2 &cm,
    double Resolution,
    std::ofstream &file,
    bool debug)
{
    using namespace std;
    if (debug)
        file << "starting to read cif!" << endl;
    if (!exists(cif))
    {
        file << "CIF does not exists!" << endl;
        return;
    }
    ifstream cif_input(cif.c_str(), ios::in);
    bvec found;
    string line;
    found.resize(7);
    for (int k = 0; k < 7; k++)
        found[k] = false;
    double a = 0.0, b = 0.0, c = 0.0, v = 0.0;
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    svec cell_keywords;
    cell_keywords.push_back("_cell_length_a");
    cell_keywords.push_back("_cell_length_b");
    cell_keywords.push_back("_cell_length_c");
    cell_keywords.push_back("_cell_angle_alpha");
    cell_keywords.push_back("_cell_angle_beta");
    cell_keywords.push_back("_cell_angle_gamma");
    cell_keywords.push_back("_cell_volume");
    if (debug)
        file << "Starting while !.eof()" << endl;
    while (!cif_input.eof())
    {
        if (debug)
            file << "While line! " << setw(80) << line
                 << setw(10) << a << found[0]
                 << setw(10) << b << found[1]
                 << setw(10) << c << found[2]
                 << setw(10) << alpha << found[3]
                 << setw(10) << beta << found[4]
                 << setw(10) << gamma << found[5]
                 << setw(10) << v << found[6] << endl;
        getline(cif_input, line);
        for (int k = 0; k < cell_keywords.size(); k++)
        {
            if (line.find(cell_keywords[k]) != string::npos)
            {
                switch (k)
                {
                case 0:
                    a = stod(line.substr(cell_keywords[k].length(), line.find("(")));
                    break;
                case 1:
                    b = stod(line.substr(cell_keywords[k].length(), line.find("(")));
                    break;
                case 2:
                    c = stod(line.substr(cell_keywords[k].length(), line.find("(")));
                    break;
                case 3:
                    alpha = stod(line.substr(cell_keywords[k].length(), line.find("(")));
                    break;
                case 4:
                    beta = stod(line.substr(cell_keywords[k].length(), line.find("(")));
                    break;
                case 5:
                    gamma = stod(line.substr(cell_keywords[k].length(), line.find("(")));
                    break;
                case 6:
                    v = stod(line.substr(cell_keywords[k].length(), line.find("(")));
                    break;
                default:
                    file << "This is weird... should never get here... aborting!" << endl;
                    return;
                }
                found[k] = true;
            }
        }
        if (found[0] == true && found[1] == true && found[2] == true && found[3] == true && found[4] == true && found[5] == true && found[6] == true)
            break;
    }
    double ca = cos(constants::PI_180 * alpha);
    double cb = cos(constants::PI_180 * beta);
    double cg = cos(constants::PI_180 * gamma);
    double sa = sin(constants::PI_180 * alpha);
    double sb = sin(constants::PI_180 * beta);
    double sg = sin(constants::PI_180 * gamma);
    double Vp = sqrt((1 - ca * ca - cb * cb - cg * cg) + 2 * (ca * cb * cg));
    double V = a * b * c * Vp;

    if (debug)
        file << "Making cm" << endl
             << a << " " << b << " " << c << " " << ca << " " << cb << " " << cg << " " << sa << " " << sb << " " << sg << " " << V << " vs. " << v << endl;

    cm[0][0] = a / 0.529177249;
    cm[1][0] = 0.0;
    cm[2][0] = 0.0;

    cm[0][1] = b * cg / 0.529177249;
    cm[1][1] = b * sg / 0.529177249;
    cm[2][1] = 0.0;

    cm[0][2] = c * cb / 0.529177249;
    cm[1][2] = (c * (ca - cb * cg) / sg) / 0.529177249;
    cm[2][2] = V / (a * b * sg) / 0.529177249;

    CoordMinMax[0] = 0.0;
    CoordMinMax[1] = 0.0;
    CoordMinMax[2] = 0.0;

    CoordMinMax[3] = (a + b * cg + c * cb) / 0.529177249;
    CoordMinMax[4] = (b * sg + c * (ca - cb * cg) / sg) / 0.529177249;
    CoordMinMax[5] = V / (a * b * sg) / 0.529177249;

    NbSteps[0] = (int)ceil(a / Resolution);
    NbSteps[1] = (int)ceil(b / Resolution);
    NbSteps[2] = (int)ceil(c / Resolution);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cm[i][j] /= NbSteps[j];

    cif_input.close();
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

bool read_fracs_ADPs_from_CIF(std::string cif, WFN &wavy, cell &unit_cell, std::ofstream &log3, bool debug)
{
    using namespace std;
    vec2 Uij, Cijk, Dijkl;
    ifstream asym_cif_input(cif.c_str(), std::ios::in);
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
                    if (is_similar(positions[labels.size()][0], wavy.atoms[i].x, -1) && is_similar(positions[labels.size()][1], wavy.atoms[i].y, -1) && is_similar(positions[labels.size()][2], wavy.atoms[i].z, -1))
                    {
                        if (debug)
                            log3 << "WFN position: " << wavy.atoms[i].x << " " << wavy.atoms[i].y << " " << wavy.atoms[i].z << endl
                                 << "Found an atom: " << fields[label_field] << " Corresponding to atom charge " << wavy.atoms[i].charge << endl;
                        wavy.atoms[i].label = fields[label_field];
                        wavy.atoms[i].frac_coords = {stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]])};
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
                    if (fields[label_field] == wavy.atoms[i].label)
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
                    if (fields[label_field] == wavy.atoms[i].label)
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
                    if (fields[label_field] == wavy.atoms[i].label)
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
        wavy.atoms[i].assign_ADPs(Uij[i], Cijk[i], Dijkl[i]);

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

const double gaussian_radial(primitive &p, double &r)
{
    return pow(r, p.type) * std::exp(-p.exp * r * r) * p.norm_const;
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

int load_basis_into_WFN(WFN &wavy, const std::array<std::vector<primitive>, 118> &b)
{
    wavy.set_basis_set_ptr(b);
    int nr_coefs = 0;
    for (int i = 0; i < wavy.atoms.size(); i++)
    {
        int current_charge = wavy.atoms[i].charge - 1;
        const primitive *basis = b[current_charge].data();
        int size = (int)b[current_charge].size();
        for (int e = 0; e < size; e++)
        {
            wavy.atoms[i].push_back_basis_set(basis[e].exp, 1.0, basis[e].type, e);
            nr_coefs += 2 * basis[e].type + 1;
        }
    }
    return nr_coefs;
}

int load_basis_into_WFN(WFN &wavy, BasisSet &b)
{
    wavy.set_basis_set_ptr(b.get_data());
    int nr_coefs = 0;
    for (int i = 0; i < wavy.atoms.size(); i++)
    {
        int current_charge = wavy.atoms[i].charge - 1;
        const std::vector<primitive> &basis = b[current_charge];
        int size = (int)b[current_charge].size();
        for (int e = 0; e < size; e++)
        {
            wavy.atoms[i].push_back_basis_set(basis[e].exp, 1.0, basis[e].type, e);
            nr_coefs += 2 * basis[e].type + 1;
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
        cout << " Recap of input:\nsize: " << arguments.size() << endl;
    }
    // This loop figures out command line options
    int argc = (int)arguments.size();
    for (int i = 0; i < arguments.size(); i++)
    {
        if (debug)
            cout << arguments[i] << endl;
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
                cout << "Looking for Anions!" << endl;
            while (i + n < argc && string(arguments[i + n]).find("-") != 0)
            {
                if (i + n - 1 > arguments.size())
                    break;
                store = arguments[i + n];
                svec Z = split_string<string>(store, " ");
                for (int r = 0; r < Z.size(); r++)
                {
                    if (debug)
                        cout << Z[r] << endl;
                    Anions.push_back(Z[r]);
                }
                n++;
            }
        }
        else if (temp == "-atom_dens")
        {
            wfn = arguments[i + 1];
            err_checkf(exists(wfn), "WFN doesn't exist", cout);
            ivec val_MOs;
            ivec val_MOs_beta;
            if (argc >= i + 3)
            {
                val_MOs = split_string<int>(arguments[i + 2], ",");
                cout << "Alpha MOs to keep: ";
                for (int j = 0; j < val_MOs.size(); j++)
                    cout << val_MOs[j] << " ";
                cout << endl;
                val_MOs_beta = split_string<int>(arguments[i + 3], ",");
                cout << "Beta MOs to keep: ";
                for (int j = 0; j < val_MOs_beta.size(); j++)
                    cout << val_MOs_beta[j] << " ";
                cout << endl;
            }
            spherically_averaged_density(*this, val_MOs, val_MOs_beta);
            exit(0);
        }
        else if (temp == "-b")
            basis_set = arguments[i + 1];
        else if (temp == "-blastest")
        {
            test_openblas();
        }
        else if (temp == "-Cation")
        {
            int n = 1;
            string store;
            if (debug)
                cout << "Looking for Cations!" << endl;
            while (i + n < argc && string(arguments[i + n]).find("-") != 0)
            {
                if (i + n - 1 > arguments.size())
                    break;
                store = arguments[i + n];
                svec Z = split_string<string>(store, " ");
                for (int r = 0; r < Z.size(); r++)
                {
                    if (debug)
                        cout << Z[r] << endl;
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
            err_checkf(exists(coef_file), "coef_file doesn't exist", cout);
            SALTED = true;
        }
        else if (temp == "-cif")
        {
            cif = arguments[i + 1];
            err_checkf(exists(cif), "CIF doesn't exist", cout);
        }
        else if (temp == "-cpus")
            threads = stoi(arguments[i + 1]);
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
                groups.push_back(split_string<int>(_temp, delimiter));
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
            double resolution = 0.025;
            double radius = 3.5;
            if (opts.size() >= 3)
            {
                resolution = opts[2];
            }
            if (opts.size() == 4)
            {
                radius = opts[3];
            }

            draw_orbital(opts[0], opts[1], resolution, radius);
            exit(0);
        }
        else if (temp == "-e_field")
            efield = stod(arguments[i + 1]);
        else if (temp == "-ECP")
        {
            ECP = true;
            if (argc >= i + 2 && string(arguments[i + 1]).find("-") != 0)
            {
                ECP_mode = stoi(arguments[i + 1]);
            }
        }
        else if (temp == "-set_ECPs")
        {
            set_ECPs = true;
            if (debug)
                cout << "Reading set ECPs" << endl;
            int j = 0;
            while (argc >= i + 2 * (j + 1) &&
                   string(arguments[i + j]).find("-") != 0 &&
                   string(arguments[i + j + 1]).find("-") != 0)
            {
                ECP_nrs.push_back(stoi(arguments[i + j]));
                ECP_elcounts.push_back(stoi(arguments[i + j + 1]));
                j += 2;
                if (debug)
                {
                    cout << j << " " << arguments[i + j] << " " << arguments[i + j + 1] << endl;
                    cout << ECP_nrs.size() << endl;
                    cout << ECP_elcounts.size() << endl;
                }
            }
        }
        else if (temp == "-ED")
            electron_diffraction = true;
        else if (temp == "-eli")
            calc = eli = true;
        else if (temp == "-elf")
            calc = elf = true;
        else if (temp == "-esp")
            calc = esp = true;
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
        else if (temp == "-hkl")
        {
            hkl = arguments[i + 1];
            err_checkf(exists(hkl), "hkl doesn't exist", cout);
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
        else if (temp == "-method")
            method = arguments[i + 1];
        else if (temp == "-merge")
        {
            svec filenames;
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
            svec filenames;
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
                groups.push_back(split_string<int>(_temp, delimiter));
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
                combined_tsc_calc_charge.push_back(stoi(arguments[i + n]));
                n++;
            }
        }
        else if (temp == "-mult")
            mult = stoi(arguments[i + 1]);
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
        else if (temp == "-radius")
            radius = stod(arguments[i + 1]);
        else if (temp == "-resolution")
            resolution = stod(arguments[i + 1]);
        else if (temp == "-rdg")
            calc = rdg = true;
        else if (temp.find("-rkpts") < 1)
            read_k_pts = true;
        else if (temp == "-rho_cube_test")
        {
            test_density_cubes(*this, log_file);
            exit(0);
        }
        else if (temp == "-rho_cube")
        {
            string wfn_name = arguments[i + 1];
            WFN wavy(0);
            cout << "Reading wavefunction: " << wfn_name << endl;
            wavy.read_known_wavefunction_format(wfn_name, cout, debug);
            cout << "Assigning ECPs" << endl;
            if (ECP)
                wavy.set_has_ECPs(true);
            cout << "Starting cube calculation" << endl;
            calc_rho_cube(wavy);
            exit(0);
        }
        else if (temp.find("-s_rho") < 1)
            s_rho = true;
        else if (temp == "-SALTED" || temp == "-salted")
        {
            SALTED = true;
            SALTED_DIR = arguments[i + 1];
        }
        else if (temp == "-DFBASIS" || temp == "-dfbasis")
        {
            SALTED_DFBASIS = arguments[i + 1];
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
        else if (temp == "-spherical_harmonic")
        {
            spherical_harmonic_test();
            exit(0);
        }
        else if (temp == "-test")
            cout << "Running in test mode!" << endl, test = true;
        else if (temp == "-thakkar_d_plot")
        {
            cout << "Making a table of Thakkar scattering factors and leaving!" << endl, thakkar_d_test(*this);
            exit(0);
        }
        else if (temp == "-twin")
        {
            twin_law.resize(twin_law.size() + 1);
            twin_law[twin_law.size() - 1].resize(9);
            for (int twl = 0; twl < 9; twl++)
                twin_law[twin_law.size() - 1][twl] = stod(arguments[i + 1 + twl]);
            if (debug)
            {
                cout << "twin_law: ";
                for (int twl = 0; twl < 9; twl++)
                    cout << setw(7) << setprecision(2) << twin_law[twin_law.size() - 1][twl];
                cout << endl;
            }
            i += 9;
        }
        else if (temp == "-old_tsc")
        {
            old_tsc = true;
        }
        else if (temp == "-tscb")
        {
            string name = arguments[i + 1];
            tsc_block<int, cdouble> blocky = tsc_block<int, cdouble>(name);
            string cif_name = "test.cif";
            if (get_ending_from_filename(name) == "tscb")
                blocky.write_tsc_file(cif_name, get_basename_without_ending(name) + ".tsc");
            else if (get_ending_from_filename(name) == "tsc")
                blocky.write_tscb_file(cif_name, get_basename_without_ending(name) + ".tscb");
            else
                err_checkf(false, "Wrong file ending!", cout);
            exit(0);
        }
        else if (temp == "-test_analytical")
        {
            test_analytical_fourier();
            exit(0);
        }
        else if (temp == "-wfn")
        {
            wfn = arguments[i + 1];
            err_checkf(exists(wfn), "Wavefunction dos not exist!", cout);
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
            std::cout << NoSpherA2_message() << help_message() << build_date() << std::endl;
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

bool exists(const std::string &name)
{
    if (FILE *file = fopen(name.c_str(), "r"))
    {
        fclose(file);
        return true;
    }
    else
    {
        return false;
    }
};

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
    //Disable setfill 0 again
    file << std::setfill(' ');
}

void write_timing_to_file(std::ostream &file,
                          std::vector<time_point> time_points,
                          std::vector<std::string> descriptions)
{
    using namespace std;
    // Check if either vector is empty
    if (time_points.empty() || descriptions.empty())
    {
        cout << "Error: Empty vector passed to write_timing_to_file" << endl;
        return;
    }

    std::chrono::microseconds total_time = std::chrono::duration_cast<std::chrono::microseconds>(time_points.back() - time_points.front());
    file << "\n\n------------------------------ Time Breakdown! ------------------------------" << endl;
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

const int sht2nbas(const int &type)
{
    const int st2bas[6]{1, 3, 6, 10, 15, 21};
    const int nst2bas[6]{11, 9, 7, 5, 4, 1};
    if (type >= 0)
        return st2bas[type];
    else
        return nst2bas[5 + type];
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

void error_check(const bool condition, const std::string &file, const int &line, const std::string &function, const std::string &error_mesasge, std::ostream &log_file)
{
    if (!condition)
    {
        log_file << "Error in " << function << " at: " << file << " : " << line << " " << error_mesasge << std::endl;
        log_file.flush();
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Error in " << function << " at: " << file << " : " << line << " " << error_mesasge << std::endl;
        exit(-1);
    }
};
void not_implemented(const std::string &file, const int &line, const std::string &function, const std::string &error_mesasge, std::ostream &log_file)
{
    log_file << function << " at: " << file << ":" << line << " " << error_mesasge << " not yet implemented!" << std::endl;
    log_file.flush();
    std::cout.rdbuf(coutbuf); // reset to standard output again
    std::cout << "Error in " << function << " at: " << file << " : " << line << " " << error_mesasge << " not yet implemented!" << std::endl;
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
        buffer[bitlen / 8 % 64] = data[i];
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

    buffer[i++] = 0x80;
    if (i > 56)
    {
        while (i < 64)
        {
            buffer[i++] = 0x00;
        }
        sha256_transform(state, buffer);
        i = 0;
    }

    while (i < 56)
    {
        buffer[i++] = 0x00;
    }

    bitlen = custom_bswap_64(bitlen);
    memcpy(buffer + 56, &bitlen, 8);
    sha256_transform(state, buffer);

    for (i = 0; i < 8; ++i)
    {
        state[i] = custom_bswap_32(state[i]);
        memcpy(hash + i * 4, &state[i], 4);
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

bool ExtractDLL(const std::string &dllName)
{
    // Convert the DLL name to a wide string
    std::wstring wideDllName = s2ws(dllName);

    // Check if the DLL already exists
    if (GetFileAttributes(wideDllName.c_str()) != INVALID_FILE_ATTRIBUTES)
    {
        return true; // DLL already exists
    }

    // Find the resource
    HRSRC hRes = FindResource(NULL, MAKEINTRESOURCE(IDR_OPENBLAS), L"DLL");
    if (!hRes)
        return false;

    // Load the resource
    HGLOBAL hResLoad = LoadResource(NULL, hRes);
    if (!hResLoad)
        return false;

    // Lock the resource to get a pointer to the data
    void *pResData = LockResource(hResLoad);
    if (!pResData)
        return false;

    // Get the size of the resource
    DWORD resSize = SizeofResource(NULL, hRes);
    if (resSize == 0)
        return false;

    // Write the resource data to a file
    std::ofstream outFile(dllName, std::ios::binary);
    if (!outFile)
        return false;

    outFile.write(reinterpret_cast<const char *>(pResData), resSize);
    outFile.close();

    return true;
}

bool check_OpenBLAS_DLL(const bool &debug)
{
    if (debug)
        std::cout << "Checking for OpenBLAS DLL" << std::endl;
    // Get the path of the executable
    char exePath[MAX_PATH];
    GetModuleFileNameA(NULL, exePath, MAX_PATH); // get path to NoSpherA2 executable
    if (debug)
        std::cout << "Executable path: " << exePath << std::endl;
    std::string exeDir = get_foldername_from_path(exePath);
    if (debug)
        std::cout << "Executable directory: " << exeDir << std::endl;

    // Define the DLL name
    std::string dllName = exeDir + "\\libopenblas.dll";
    if (debug)
        std::cout << "DLL name: " << dllName << std::endl;
    if (exists(dllName))
        return true; // DLL already exists
    else
    {
        if (debug)
            std::cout << "DLL does not exist, extracting it form teh executable!" << std::endl;
        // Extract the DLL if it does not exist
        if (!ExtractDLL(dllName))
        {
            std::cout << "Failed to extract DLL" << std::endl;
            return false;
        }
        if (debug)
            std::cout << "DLL extracted successfully!" << std::endl;
    }
    return true;
}
#endif