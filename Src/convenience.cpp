#include "convenience.h"
#include "cell.h"
#include "tsc_block.h"
#include "test_functions.h"

using namespace std;

string help_message()
{
    std::string t = "\n----------------------------------------------------------------------------\n";
    t.append("          These commands and arguments are known by NoSpherA2:\n");
    t.append("----------------------------------------------------------------------------\n\n");
    t.append("   -wfn            <FILENAME>.xxx           Read the following wavefunction file.\n");
    t.append("                                            Supported filetypes: .wfn/wfx/ffn; .molden; .xyz; .gbw; fch* (UNTESTED!)\n");
    t.append("   -fchk           <FILENAME>.fchk          Write a wavefunction to the given filename\n");
    t.append("   -b              <FILENAME>               Read this basis set\n");
    t.append("   -d              <PATH>                   Path to basis_sets directory with basis_sets in tonto style\n");
    t.append("   -dmin		     <NUMBER>                   Minimum d-spacing to consider for scattering factors (repalaces hkl file)\n");
    t.append("   --help/-help/--h                         print this help\n");
    t.append("   -v                                       Turn on Verbose (debug) Mode (Slow and a LOT of output!)\n");
    t.append("   -v2                                      Even more stuff\n");
    t.append("   -mult           <NUMBER>                 Input multiplicity of wavefunction (otherwise attempted to be read from the wfn)\n");
    t.append("   -method         <METHOD NAME>            Can be RKS or RHF to distinguish between DFT and HF\n");
    t.append("   -cif            <FILENAME>.cif           CIF to get labels of atoms to use for calculation of scatteriung factors\n");
    t.append("   -IAM                                     Make scattering factors based on Thakkar functions for atoms in CIF\n");
    t.append("   -xyz            <FILENAME>.xyz           Read atom positions from this xyz file for IAM\n");
    t.append("   -hkl            <FILENAME>.hkl           hkl file (ideally merged) to use for calculation of form factors.\n");
    t.append("   -group          <LIST OF INT NUMBERS>    Disorder groups to be read from the CIF for consideration as asym unit atoms (space separated).\n");
    t.append("   -acc            0,1,2,3,4...             Accuracy of numerical grids used, where the bumber indicates a pre-defined level. 4 should be considered maximum,\n");
    t.append("                                            anything above will most likely introduce numberical error and is just implemented for testing purposes.");
    t.append("   -gbw2wfn                                 Only reads wavefucntion from .gbw specified by -wfn and prints it into .wfn format.\n");
    t.append("   -tscb           <FILENAME>.tsb           Convert binary tsc file to bigger, less accurate human-readable form.\n");
    t.append("   -twin     3x3 floating-point-matrix in the form -1 0 0 0 -1 0 0 0 -1 which contains the twin matrix to use.\n");
    t.append("             If there is more than a single twin law to be used, use the twin command multiple times.\n");
    t.append("   -merge          <List of .tsc files>     Names/Paths to .tsc files to be merged.\n");
    t.append("   -merge_nocheck  <List of .tsc files>     Names/Paths to .tsc files to be merged. They need to have identical hkl values.\n");
    t.append("   -mtc            <List of .wfns + parts>  Performs calculation for a list of wavefunctions (=Multi-Tsc-Calc), where asymmetric unit is.\n");
    t.append("                                            taken from given CIF. Also disorder groups are required per file as comma separated list\n");
    t.append("                                            without spaces.\n   Typical use Examples:\n");
    t.append("      Normal:       NoSpherA2.exe -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7\n");
    t.append("      thakkar-tsc:  NoSpherA2.exe -cif A.cif -hkl A.hkl -xyz A.xyz -acc 1 -cpus 7 -IAM\n");
    t.append("      Disorder:     NoSpherA2.exe -cif A.cif -hkl A.hkl -acc 1 -cpus 7 -mtc 1.wfn 0,1 2.wfn 0,2 3.wfn 0,3\n");
    t.append("      fragHAR:      NoSpherA2.exe -cif A.cif -hkl A.hkl -acc 1 -cpus 7 -cmtc 1.wfn 1.cif 0 2.wfn 2.cif 0 3_1.wfn 3_1.cif 0,1 3_2.wfn 3_2.cif 0,2\n");
    t.append("      merging tscs: NoSpherA2.exe -merge A.tsc B.tsc C.tsc\n");
    t.append("      merge tsc(2): NoSpherA2.exe -merge_nocheck A.tsc B.tsc C.tsc  (MAKE SURE THEY HAVE IDENTICAL HKL INIDCES!!)\n");
    t.append("      convert tsc:  NoSpherA2.exe -tscb A.tsc\n");
    t.append("      convert gbw:  NoSpherA2.exe -gbw2wfn -wfn A.gbw\n");
    t.append("      twin law:     NoSpherA2.exe -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7 -twin -1 0 0 0 -1 0 0 0 -1\n");
    return t;
}
string NoSpherA2_message()
{
    string t = "    _   __     _____       __              ___   ___\n";
    t.append("   / | / /___ / ___/____  / /_  ___  _____/   | |__ \\\n");
    t.append("  /  |/ / __ \\\\__ \\/ __ \\/ __ \\/ _ \\/ ___/ /| | __/ /\n");
    t.append(" / /|  / /_/ /__/ / /_/ / / / /  __/ /  / ___ |/ __/\n");
    t.append("/_/ |_/\\____/____/ .___/_/ /_/\\___/_/  /_/  |_/____/\n");
    t.append("                /_/\n");
    t.append("This software is part of the cuQCT software suite developed by Florian Kleemiss.\n");
    t.append("Please give credit and cite corresponding pieces!\n");
    t.append("List of contributors of pieces of code or funcitonality:\n");
    t.append("      Florian Kleemiss,\n");
    t.append("      Emmanuel Hupf,\n");
    t.append("      Alessandro Genoni,\n");
    t.append("      and many more in communications or by feedback!\n");
    t.append("NoSpherA2 was published at  : Kleemiss et al. Chem.Sci., 2021, 12, 1675 - 1692.\n");
    t.append("Slater IAM was published at : Kleemiss et al. J. Appl. Cryst 2024, 57, 161 - 174.\n");
    return t;
}

string build_date()
{
    return ("This Executable was built on: " + string(__DATE__) + " " + string(__TIME__) + "\n");
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

void cls()
{
    //    cout << string( 100, '\n' );
#ifdef _WIN32
    if (system("CLS"))
        cout << "this should not happen...!" << endl;
#else
    if (system("clear"))
        cout << "this should not happen...!" << endl;
#endif
};

string atnr2letter(const int &nr)
{
    vector<string> Labels{"DM", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
    if (nr == 0)
    {
        // Exception for Q peaks in residual maps
        return "Q";
    }
    if (nr > 103 || nr < 0)
    {
        if (nr == 119)
        {
            // Exception for Q in ECPs from ORCA
            return "Q";
        }
        cout << "Only yet implemented from H-Lr, ask Florian for improvements or give a reasonable number between 1-103!" << endl;
        return ("PROBLEM");
    }
    else
        return Labels[nr];
};

int get_Z_from_label(const char *tmp)
{
    if (strcmp(tmp, "H") == 0)
        return 0;
    else if (strcmp(tmp, "D") == 0)
        return 0;
    else if (strcmp(tmp, "T") == 0)
        return 0;
    else if (strcmp(tmp, "He") == 0)
        return 1;
    else if (strcmp(tmp, "Li") == 0)
        return 2;
    else if (strcmp(tmp, "Be") == 0)
        return 3;
    else if (strcmp(tmp, "B") == 0)
        return 4;
    else if (strcmp(tmp, "C") == 0)
        return 5;
    else if (strcmp(tmp, "N") == 0)
        return 6;
    else if (strcmp(tmp, "O") == 0)
        return 7;
    else if (strcmp(tmp, "F") == 0)
        return 8;
    else if (strcmp(tmp, "Ne") == 0)
        return 9;
    else if (strcmp(tmp, "Na") == 0)
        return 10;
    else if (strcmp(tmp, "Mg") == 0)
        return 11;
    else if (strcmp(tmp, "Al") == 0)
        return 12;
    else if (strcmp(tmp, "Si") == 0)
        return 13;
    else if (strcmp(tmp, "P") == 0)
        return 14;
    else if (strcmp(tmp, "S") == 0)
        return 15;
    else if (strcmp(tmp, "Cl") == 0)
        return 16;
    else if (strcmp(tmp, "Ar") == 0)
        return 17;
    else if (strcmp(tmp, "K") == 0)
        return 18;
    else if (strcmp(tmp, "Ca") == 0)
        return 19;
    else if (strcmp(tmp, "Sc") == 0)
        return 20;
    else if (strcmp(tmp, "Ti") == 0)
        return 21;
    else if (strcmp(tmp, "V") == 0)
        return 22;
    else if (strcmp(tmp, "Cr") == 0)
        return 23;
    else if (strcmp(tmp, "Mn") == 0)
        return 24;
    else if (strcmp(tmp, "Fe") == 0)
        return 25;
    else if (strcmp(tmp, "Co") == 0)
        return 26;
    else if (strcmp(tmp, "Ni") == 0)
        return 27;
    else if (strcmp(tmp, "Cu") == 0)
        return 28;
    else if (strcmp(tmp, "Zn") == 0)
        return 29;
    else if (strcmp(tmp, "Ga") == 0)
        return 30;
    else if (strcmp(tmp, "Ge") == 0)
        return 31;
    else if (strcmp(tmp, "As") == 0)
        return 32;
    else if (strcmp(tmp, "Se") == 0)
        return 33;
    else if (strcmp(tmp, "Br") == 0)
        return 34;
    else if (strcmp(tmp, "Kr") == 0)
        return 35;
    else if (strcmp(tmp, "Rb") == 0)
        return 36;
    else if (strcmp(tmp, "Sr") == 0)
        return 37;
    else if (strcmp(tmp, "Y") == 0)
        return 38;
    else if (strcmp(tmp, "Zr") == 0)
        return 39;
    else if (strcmp(tmp, "Nb") == 0)
        return 40;
    else if (strcmp(tmp, "Mo") == 0)
        return 41;
    else if (strcmp(tmp, "Tc") == 0)
        return 42;
    else if (strcmp(tmp, "Ru") == 0)
        return 43;
    else if (strcmp(tmp, "Rh") == 0)
        return 44;
    else if (strcmp(tmp, "Pd") == 0)
        return 45;
    else if (strcmp(tmp, "Ag") == 0)
        return 46;
    else if (strcmp(tmp, "Cd") == 0)
        return 47;
    else if (strcmp(tmp, "In") == 0)
        return 48;
    else if (strcmp(tmp, "Sn") == 0)
        return 49;
    else if (strcmp(tmp, "Sb") == 0)
        return 50;
    else if (strcmp(tmp, "Te") == 0)
        return 51;
    else if (strcmp(tmp, "I") == 0)
        return 52;
    else if (strcmp(tmp, "Xe") == 0)
        return 53;
    else if (strcmp(tmp, "Cs") == 0)
        return 54;
    else if (strcmp(tmp, "Ba") == 0)
        return 55;
    else if (strcmp(tmp, "La") == 0)
        return 56;
    else if (strcmp(tmp, "Ce") == 0)
        return 57;
    else if (strcmp(tmp, "Pr") == 0)
        return 58;
    else if (strcmp(tmp, "Nd") == 0)
        return 59;
    else if (strcmp(tmp, "Pm") == 0)
        return 60;
    else if (strcmp(tmp, "Sm") == 0)
        return 61;
    else if (strcmp(tmp, "Eu") == 0)
        return 62;
    else if (strcmp(tmp, "Gd") == 0)
        return 63;
    else if (strcmp(tmp, "Tb") == 0)
        return 64;
    else if (strcmp(tmp, "Dy") == 0)
        return 65;
    else if (strcmp(tmp, "Ho") == 0)
        return 66;
    else if (strcmp(tmp, "Er") == 0)
        return 67;
    else if (strcmp(tmp, "Tm") == 0)
        return 68;
    else if (strcmp(tmp, "Yb") == 0)
        return 69;
    else if (strcmp(tmp, "Lu") == 0)
        return 70;
    else if (strcmp(tmp, "Hf") == 0)
        return 71;
    else if (strcmp(tmp, "Ta") == 0)
        return 72;
    else if (strcmp(tmp, "W") == 0)
        return 73;
    else if (strcmp(tmp, "Re") == 0)
        return 74;
    else if (strcmp(tmp, "Os") == 0)
        return 75;
    else if (strcmp(tmp, "Ir") == 0)
        return 76;
    else if (strcmp(tmp, "Pt") == 0)
        return 77;
    else if (strcmp(tmp, "Au") == 0)
        return 78;
    else if (strcmp(tmp, "Hg") == 0)
        return 79;
    else if (strcmp(tmp, "Ti") == 0)
        return 80;
    else if (strcmp(tmp, "Pb") == 0)
        return 81;
    else if (strcmp(tmp, "Bi") == 0)
        return 82;
    else if (strcmp(tmp, "Po") == 0)
        return 83;
    else if (strcmp(tmp, "At") == 0)
        return 84;
    else if (strcmp(tmp, "Rn") == 0)
        return 85;

    else if (strcmp(tmp, "Fr") == 0)
        return 86;
    else if (strcmp(tmp, "Ra") == 0)
        return 87;
    else if (strcmp(tmp, "Ac") == 0)
        return 88;
    else if (strcmp(tmp, "Th") == 0)
        return 89;
    else if (strcmp(tmp, "Pa") == 0)
        return 90;
    else if (strcmp(tmp, "U") == 0)
        return 91;
    else if (strcmp(tmp, "Np") == 0)
        return 92;
    else if (strcmp(tmp, "Pu") == 0)
        return 93;
    else if (strcmp(tmp, "Am") == 0)
        return 94;
    else if (strcmp(tmp, "Cm") == 0)
        return 95;
    else if (strcmp(tmp, "Bk") == 0)
        return 96;
    else if (strcmp(tmp, "Cf") == 0)
        return 97;
    else if (strcmp(tmp, "Es") == 0)
        return 98;
    else if (strcmp(tmp, "Fm") == 0)
        return 99;
    else if (strcmp(tmp, "Md") == 0)
        return 100;
    else if (strcmp(tmp, "No") == 0)
        return 101;
    else if (strcmp(tmp, "Lr") == 0)
        return 102;
    else if (strcmp(tmp, "Rf") == 0)
        return 103;
    else if (strcmp(tmp, "Db") == 0)
        return 104;
    else if (strcmp(tmp, "Sg") == 0)
        return 105;
    else if (strcmp(tmp, "Bh") == 0)
        return 106;
    else if (strcmp(tmp, "Hs") == 0)
        return 107;
    else if (strcmp(tmp, "Mt") == 0)
        return 108;
    else if (strcmp(tmp, "Ds") == 0)
        return 109;
    else if (strcmp(tmp, "Rg") == 0)
        return 110;
    else
        return -1;
}
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
void copy_file(string &from, string &to)
{
    ifstream source(from.c_str(), ios::binary);
    ofstream dest(to.c_str(), ios::binary);

    dest << source.rdbuf();

    source.close();
    dest.close();
};

//---------------------------Configuration files ---------------------------------------------------

string get_home_path(void)
{
#ifdef _WIN32
    string temp1 = getenv("HOMEDRIVE");
    string temp2 = getenv("HOMEPATH");
    temp1.append(temp2);
    return temp1;
#else
    string home = getenv("HOME");
    return home;
#endif
}

void join_path(string &s1, string &s2)
{
#ifdef _WIN32
    s1.append("\\");
#else
    if (s1.substr(s1.length() - 1) != "/")
        s1.append("/");
#endif
    s1.append(s2);
}

string get_filename_from_path(const string &input)
{
#ifdef _WIN32
    return input.substr(input.rfind("\\") + 1);
#else
    return input.substr(input.rfind("/") + 1);
#endif
}

string get_foldername_from_path(const string &input)
{
#ifdef _WIN32
    return input.substr(0, input.rfind("\\") + 1);
#else
    return input.substr(0, input.rfind("/") + 1);
#endif
}

string get_basename_without_ending(const string &input)
{
    return input.substr(0, input.rfind("."));
}

void write_template_confi()
{
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

int program_confi(string &gaussian_path, string &turbomole_path, string &basis, int &ncpus, double &mem, bool debug, bool expert, unsigned int counter)
{
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
                cout << "Length for: " << i << ";" << j << ": " << length << ", min_length: " << min_length << endl;
            if (length < min_length)
                min_length = length;
        }
    }
    if (debug)
    {
        if (min_length < 2)
            cout << "Decided it's written in Angstrom" << endl;
        else
            cout << "Decided it's written in Bohr" << endl;
    }
    return (!(min_length < 2));
};

int filetype_identifier(string &file, bool debug)
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
    if (debug)
    {
        cout << "Testing WFN:  " << file.find(".wfn") << endl
             << "Testing out:  " << file.find(".out") << endl
             << "Testing FFN:  " << file.find(".ffn") << endl
             << "Testing CUB:  " << file.find(".cub") << endl
             << "Testing CUBE: " << file.find(".cube") << endl
             << "Testing Grid: " << file.find(".grd") << endl
             << "Testing fchk: " << file.find(".fchk") << endl
             << "Testing FChk: " << file.find(".FChk") << endl
             << "Testing Fchk: " << file.find(".Fchk") << endl;
        cout << "string::npos: " << string::npos << endl;
    }
    int temp_type = 0;
    size_t found, temp;
    temp = 0;
    if (debug)
        cout << "Temp before any checks: " << temp << endl;
    vector<string> types{".out", ".wfn", ".ffn", ".cub", ".cube", ".grd", ".fchk", ".Fchk", ".FChk"};
    if (file.find(".wfn") != string::npos)
    {
        if (debug)
            cout << "Checking for"
                 << ".wfn" << endl;
        temp_type = 2;
        found = file.rfind(".wfn");
        if (debug)
            cout << "Found: " << found << endl;
        for (int i = 0; i < types.size(); i++)
            if (file.rfind(types[i]) >= found && file.rfind(types[i]) != string::npos)
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

string go_get_string(ifstream &file, string search, bool rewind)
{
    if (rewind)
    {
        file.clear();
        file.seekg(0, file.beg);
    }
    string line;
    while (line.find(search) == string::npos && !file.eof() && getline(file, line))
        continue;
    if (file.eof())
        return "";
    else
        return line;
}

string shrink_string(string &input)
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

string shrink_string_to_atom(string &input, const int &atom_number)
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
    string temp = atnr2letter(atom_number);
    err_checkf(temp != "PROBLEM", "Problem identifying atoms!", std::cout);
    if (input.find(temp) != 1)
        return temp;
    if (temp != "PROBLEM")
        while (input.size() > temp.size())
            input.pop_back();
    return input;
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

void select_cubes(vector<vector<unsigned int>> &selection, vector<WFN> &wavy, unsigned int nr_of_cubes, bool wfnonly, bool debug)
{
    // asks which wfn to use, if wfnonly is set or whcih cubes up to nr of cubes to use
    // Returns values in selection[0][i] for iths selection of wavefunction and
    //  selection[1][i] for iths selection of cube
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

bool unsaved_files(vector<WFN> &wavy)
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
    vector<vector<double>> PosAtoms;
    PosAtoms.resize(wavy.get_ncen());
    for (int i = 0; i < wavy.get_ncen(); i++)
        PosAtoms[i].resize(3);
    bool bohrang = !check_bohr(wavy, false);
    if (no_bohr)
        bohrang = true;
    for (int j = 0; j < wavy.get_ncen(); j++)
    {

        PosAtoms[j][0] = wavy.atoms[j].x;
        PosAtoms[j][1] = wavy.atoms[j].y;
        PosAtoms[j][2] = wavy.atoms[j].z;
        if (!bohrang)
        {
            cout << "Dividing atom positions!" << endl;
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
    double temp_rad = constants::ang2bohr(Radius);
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
    string cif,
    double *CoordMinMax,
    int *NbSteps,
    vector<vector<double>> &cm,
    double Resolution,
    ofstream &file,
    bool debug)
{
    if (debug)
        file << "starting to read cif!" << endl;
    if (!exists(cif))
    {
        file << "CIF does not exists!" << endl;
        return;
    }
    ifstream cif_input(cif.c_str(), ios::in);
    vector<bool> found;
    string line;
    found.resize(7);
    for (int k = 0; k < 7; k++)
        found[k] = false;
    double a = 0.0, b = 0.0, c = 0.0, v = 0.0;
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    vector<string> cell_keywords;
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

const double normgauss(const int &type, const double &exp)
{
    int t[3];
    if (type > 0)
    {
        type2vector(type, t);
        err_checkf(t[0] != -1, "Problem with type2vector!", std::cout);
        err_checkf(t[1] != -1, "Problem with type2vector!", std::cout);
        err_checkf(t[2] != -1, "Problem with type2vector!", std::cout);
    }
    else
        t[0] = t[1] = t[2] = 0;
    long long int temp = constants::ft[t[0]] * constants::ft[t[1]] * constants::ft[t[2]];
    long long int temp2 = constants::ft[2 * t[0]] * constants::ft[2 * t[1]] * constants::ft[2 * t[2]];
    return pow(2 * exp / constants::PI, 0.75) * sqrt(pow(8 * exp, t[0] + t[1] + t[2]) * temp / temp2);
};
bool generate_sph2cart_mat(vector<vec> &p, vector<vec> &d, vector<vec> &f, vector<vec> &g)
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
bool generate_cart2sph_mat(vector<vec> &d, vector<vec> &f, vector<vec> &g, vector<vec> &h)
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

const int type_vector[168]{
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    2, 0, 0,
    0, 2, 0,
    0, 0, 2,
    1, 1, 0,
    1, 0, 1,
    0, 1, 1,
    3, 0, 0,
    0, 3, 0,
    0, 0, 3,
    2, 1, 0,
    2, 0, 1,
    0, 2, 1,
    1, 2, 0,
    1, 0, 2,
    0, 1, 2,
    1, 1, 1,
    0, 0, 4,
    0, 1, 3,
    0, 2, 2,
    0, 3, 1,
    0, 4, 0,
    1, 0, 3,
    1, 1, 2,
    1, 2, 1,
    1, 3, 0,
    2, 0, 2,
    2, 1, 1,
    2, 2, 0,
    3, 0, 1,
    3, 1, 0,
    4, 0, 0,
    0, 0, 5,
    0, 1, 4,
    0, 2, 3,
    0, 3, 2,
    0, 4, 1,
    0, 5, 0,
    1, 0, 4,
    1, 1, 3,
    1, 2, 2,
    1, 3, 1,
    1, 4, 0,
    2, 0, 3,
    2, 1, 2,
    2, 2, 1,
    2, 3, 0,
    3, 0, 2,
    3, 1, 1,
    3, 2, 0,
    4, 0, 1,
    4, 1, 0,
    5, 0, 0};

void type2vector(
    const int &index,
    int *vector)
{
    if (index < 1 || index > 56)
    {
        vector[0] = -1;
        vector[1] = -1;
        vector[2] = -1;
        return;
    }
    const int temp = index - 1;
    vector[0] = type_vector[temp * 3];
    vector[1] = type_vector[temp * 3 + 1];
    vector[2] = type_vector[temp * 3 + 2];
}

bool read_fracs_ADPs_from_CIF(string cif, WFN &wavy, cell &unit_cell, ofstream &log3, bool debug)
{
    vector<vec> Uij, Cijk, Dijkl;
    ifstream asym_cif_input(cif.c_str(), std::ios::in);
    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    string line;
    vector<string> labels;
    int count_fields = 0;
    int position_field[3] = {0, 0, 0};
    int label_field = 100;
    vector<vector<double>> positions;
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
                vector<string> fields;
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
                vector<string> fields;
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
                vector<string> fields;
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
                vector<string> fields;
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

bool read_fchk_integer_block(ifstream &in, string heading, ivec &result, bool rewind)
{
    if (result.size() != 0)
        result.clear();
    string line = go_get_string(in, heading, rewind);
    int limit = read_fchk_integer(line);
    int run = 0;
    int temp;
    getline(in, line);
    while (run < limit)
    {
        if (in.eof())
            return false;
        temp = stoi(line.substr(12 * (run % 6), 12 * (run % 6 + 1)));
        result.push_back(temp);
        run++;
        if (run % 6 == 0)
            getline(in, line);
    }
    return true;
};
bool read_fchk_double_block(ifstream &in, string heading, vec &result, bool rewind)
{
    if (result.size() != 0)
        result.clear();
    string line = go_get_string(in, heading, rewind);
    int limit = read_fchk_integer(line);
    int run = 0;
    double temp;
    getline(in, line);
    while (run < limit)
    {
        if (in.eof())
            return false;
        temp = stod(line.substr(16 * (run % 5), 16 * (run % 5 + 1)));
        result.push_back(temp);
        run++;
        if (run % 5 == 0)
            getline(in, line);
    }
    return true;
};
int read_fchk_integer(string in)
{
    return stoi(in.substr(49, in.length() - 49));
};
double read_fchk_double(string in)
{
    return stod(in.substr(49, in.length() - 49));
};
int read_fchk_integer(std::ifstream &in, std::string search, bool rewind)
{
    string temp = go_get_string(in, search, rewind);
    return stoi(temp.substr(49, temp.length() - 49));
};
double read_fchk_double(std::ifstream &in, std::string search, bool rewind)
{
    string temp = go_get_string(in, search, rewind);
    return stod(temp.substr(49, temp.length() - 49));
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

const double calc_density_ML(double &x,
                             double &y,
                             double &z,
                             vec &coefficients,
                             std::vector<atom> &atoms,
                             const int &exp_coefs)
{
    double dens = 0, radial = 0;
    int coef_counter = 0;
    int e = 0, size = 0;
    for (int a = 0; a < atoms.size(); a++)
    {
        size = (int)atoms[a].basis_set.size();
        basis_set_entry *bf;
        double d[4]{
            x - atoms[a].x,
            y - atoms[a].y,
            z - atoms[a].z, 0.0};
        // store r in last element
        d[3] = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        // normalize distances for spherical harmonic
        for (e = 0; e < 3; e++)
            d[e] /= d[3];
        for (e = 0; e < size; e++)
        {
            bf = &atoms[a].basis_set[e];
            primitive p(a, bf->type, bf->exponent, bf->coefficient);
            radial = gaussian_radial(p, d[3]);
            for (int m = -p.type; m <= p.type; m++)
            {
                // m+p.type should yield just the running index of coefficents, since we start at -p.type
                dens += coefficients[coef_counter + m + p.type] * radial * spherical_harmonic(p.type, m, d);
            }
            coef_counter += (2 * p.type + 1);
        }
    }
    err_checkf(coef_counter == exp_coefs, "WRONG NUMBER OF COEFFICIENTS! " + std::to_string(coef_counter) + " vs. " + std::to_string(exp_coefs), std::cout);
    return dens;
}

const double calc_density_ML(double &x,
                             double &y,
                             double &z,
                             vec &coefficients,
                             std::vector<atom> &atoms,
                             const int &exp_coefs,
                             const int &atom_nr)
{
    double dens = 0, radial = 0;
    int coef_counter = 0;
    int e = 0, size = 0;
    for (int a = 0; a < atoms.size(); a++)
    {
        if (a == atom_nr)
        {
            size = (int)atoms[a].basis_set.size();
            basis_set_entry *bf;
            double d[4]{
                x - atoms[a].x,
                y - atoms[a].y,
                z - atoms[a].z, 0.0};
            // store r in last element
            d[3] = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            // normalize distances for spherical harmonic
            for (e = 0; e < 3; e++)
                d[e] /= d[3];
            for (e = 0; e < size; e++)
            {
                bf = &atoms[a].basis_set[e];
                primitive p(a, bf->type, bf->exponent, bf->coefficient);
                radial = gaussian_radial(p, d[3]);
                for (int m = -p.type; m <= p.type; m++)
                {
                    // m+p.type should yield just the running index of coefficents, since we start at -p.type
                    dens += coefficients[coef_counter + m + p.type] * radial * spherical_harmonic(p.type, m, d);
                }
                coef_counter += (2 * p.type + 1);
            }
        }
        else
        {
            size = (int)atoms[a].basis_set.size();
            for (e = 0; e < size; e++)
            {
                coef_counter += (2 * atoms[a].basis_set[e].type + 1);
            }
        }
    }
    err_checkf(coef_counter == exp_coefs, "WRONG NUMBER OF COEFFICIENTS! " + std::to_string(coef_counter) + " vs. " + std::to_string(exp_coefs), std::cout);
    return dens;
}

int load_basis_into_WFN(WFN &wavy, const std::vector<std::vector<primitive>> &b)
{
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

double get_decimal_precision_from_CIF_number(string &given_string)
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
        string temp = given_string.substr(open_bracket + 1, size_of_precision);
        precision = stoi(temp);
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
                vector<string> Z = split_string<string>(store, " ");
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
                vector<string> Z = split_string<string>(store, " ");
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
        else if (temp == "-IAM")
            iam_switch = true;
        else if (temp == "-lap")
            calc = lap = true;
        else if (temp == "-method")
            method = arguments[i + 1];
        else if (temp == "-merge")
        {
            vector<string> filenames;
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
            vector<string> filenames;
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
        else if (temp == "-ML_test")
        {
            ML_test();
            exit(0);
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
        else if (temp == "-perf_benchmark")
        {
            test_timing();
            exit(0);
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
        else if (temp.find("-s_rho") < 1)
            s_rho = true;
        else if (temp.find("-SALTED_BECKE") < 1 || temp.find("-salted_becke") < 1)
            SALTED_BECKE = true;
        else if (temp == "-SALTED" || temp == "-salted")
            SALTED = true;
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
        else if (temp == "-test-ecp")
        {
            d_sfac_scan = fromString<double>(arguments[i + 1]);
            sfac_scan_ECP(*this, log_file);
            exit(0);
        }
        else if (temp == "-test-ecp-pot")
        {
            test_esp_dens();
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
        else if (temp == "-test-core")
        {
            test_core_dens();
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
            name = "test.cif";
            blocky.write_tsc_file(name);
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
        string temp = argv[i];
        arguments.push_back(temp);
        if (temp.find("-") > 0)
            continue;
        else if (temp == "-v" || temp == "-v2")
            cout << "Turning on verbose mode!" << endl, debug = true;
        else if (temp == "--h" || temp == "-h" || temp == "-help" || temp == "--help")
        {
            cout << NoSpherA2_message() << help_message() << build_date() << endl;
            exit(0);
        }
    }
};

const double spherical_harmonic(const int &l, const int &m, const double *d)
{
    /*Here d[0] = x
                 d[1] = y
                 d[2] = z
                 d[3] = r^2 IGNORED
                 d[4] = r   IGNORED
                 */
    // Will need extension for up to l=8
    // calc spherical harmonic
    double SH = 0, x = d[0], y = d[1], z = d[2];
    switch (l)
    {
    case 0: // S
        SH = constants::c_1_4p;
        break;
    case 1:
        switch (m)
        {
        case 0: // P 0 Z
            SH = constants::c_3_4p * z;
            break;
        case 1: // P 1 X
            SH = constants::c_3_4p * x;
            break;
        case -1: // P -1 Y
            SH = constants::c_3_4p * y;
            break;
        default:
            err_not_impl_f("Wrong spherical harmonic called!", std::cout);
        }
        break;
    case 2:
        switch (m)
        {
        case 0: // D 0 Z2
            SH = constants::c_5_16p * (3 * pow(z, 2) - 1.0);
            break;
        case 1: // D 1 XZ
            SH = constants::c_15_4p * x * z;
            break;
        case -1: // D -1 YZ
            SH = constants::c_15_4p * y * z;
            break;
        case 2: // D 2 X2-Y2
            SH = constants::c_15_16p * (pow(x, 2) - pow(y, 2));
            break;
        case -2: // D -2 XY
            SH = constants::c_15_4p * y * x;
            break;
        default:
            err_not_impl_f("Wrong spherical harmonic called!", std::cout);
        }
        break;
    case 3:
        switch (m)
        {
        case 0: // F 0 Z3
            SH = constants::c_7_16p * (5 * pow(z, 3) - 3 * z);
            break;
        case 1: // F 1 XZZ
            SH = constants::c_21_32p * x * (5 * pow(z, 2) - 1.0);
            break;
        case -1: // F -1 YZZ
            SH = constants::c_21_32p * y * (5 * pow(z, 2) - 1.0);
            break;
        case 2: // F 2 Z(X2-Y2)
            SH = constants::c_105_16p * ((pow(x, 2) - pow(y, 2)) * z);
            break;
        case -2: // F -2 XYZ
            SH = constants::c_105_4p * x * y * z;
            break;
        case 3: // F 3 X(X^2-3Y^2)
            SH = constants::c_35_32p * x * (pow(x, 2) - 3 * pow(y, 2));
            break;
        case -3: // F -3 Y(3X^2-Y^2)
            SH = constants::c_35_32p * y * (3 * pow(x, 2) - pow(y, 2));
            break;
        default:
            err_not_impl_f("Wrong spherical harmonic called!", std::cout);
        }
        break;
    case 4:
        switch (m)
        {
        case 0: // G 0 Z^4
            SH = constants::c_9_256p * (35 * pow(z, 4) - 30 * pow(z, 2) + 3.0);
            break;
        case 1: // G 1 X(7Z^3-3ZR^2)
            SH = constants::c_45_32p * x * (7 * pow(z, 3) - 3 * z);
            break;
        case -1: // G -1 Y(7Z^2-3ZR^2)
            SH = constants::c_45_32p * y * (7 * pow(z, 3) - 3 * z);
            break;
        case 2: // G 2
            SH = constants::c_45_64p * (pow(x, 2) - pow(y, 2)) * (7 * pow(z, 2) - 1.0);
            break;
        case -2: // G -2
            SH = constants::c_45_16p * x * y * (7 * pow(z, 2) - 1.0);
            break;
        case 3: // G 3 XZ(X^2-3Y^2)
            SH = constants::c_315_32p * x * (pow(x, 2) - 3 * pow(y, 2)) * z;
            break;
        case -3: // G -3 XZ(3X^2-Y^2)
            SH = constants::c_315_32p * y * (3 * pow(x, 2) - pow(y, 2)) * z;
            break;
        case 4: // G 4 X^2(X^-3Y^2)-Y^2(3X^2-Y^2)
            SH = constants::c_315_256p * ((pow(x, 2) * (pow(x, 2) - 3 * pow(y, 2))) -
                                          (pow(y, 2) * (3 * pow(x, 2) - pow(y, 2))));
            break;
        case -4: // G -4 XY(X^2-Y^2)
            SH = constants::c_315_16p * x * y * (pow(x, 2) - pow(y, 2));
            break;
        default:
            err_not_impl_f("Wrong spherical harmonic called!", std::cout);
        }
        break;
    case 5:
        switch (m)
        {
        case 0: // H Z^5
            SH = constants::c_11_256p * (63 * pow(z, 5) - 70 * pow(z, 3) + 15 * z);
            break;
        case 1:
            SH = constants::c_165_256p * x * (21 * pow(z, 4) - 14 * pow(z, 2) + 1.0);
            break;
        case -1:
            SH = constants::c_165_256p * y * (21 * pow(z, 4) - 14 * pow(z, 2) + 1.0);
            break;
        case 2:
            SH = constants::c_1155_64p * (pow(x, 2) - pow(y, 2)) * (3 * pow(z, 3) - z);
            break;
        case -2:
            SH = constants::c_1155_64p * 2 * x * y * (3 * pow(z, 3) - z);
            break;
        case 3:
            SH = constants::c_385_512p * x * (pow(x, 2) - 3 * pow(y, 2)) * (9 * pow(z, 2) - 1.0);
            break;
        case -3:
            SH = constants::c_385_512p * y * (3 * pow(x, 2) - pow(y, 2)) * (9 * pow(z, 2) - 1.0);
            break;
        case 4:
            SH = constants::c_3465_256p * (pow(x, 4) - 6 * x * x * y * y + pow(y, 4)) * z;
            break;
        case -4:
            SH = -constants::c_3465_256p * (4 * x * pow(y, 3) - 4 * pow(x, 3) * y) * z;
            break;
        case 5:
            SH = constants::c_693_2048p * (2 * pow(x, 5) - 20 * pow(x, 3) * pow(y, 2) + 10 * x * pow(y, 4));
            break;
        case -5:
            SH = constants::c_693_2048p * (2 * pow(y, 5) - 20 * pow(x, 2) * pow(y, 3) + 10 * y * pow(x, 4));
            break;
        default:
            err_not_impl_f("Wrong spherical harmonic called!", std::cout);
        }
        break;
    case 6:
        switch (m)
        {
        case 0: // I Z^6
            SH = constants::c_13_1024p * (231 * pow(z, 6) - 315 * pow(z, 4) + 105 * pow(z, 2) - 5);
            break;
        case 1:
            SH = constants::c_273_256p * x * (21 * pow(z, 4) - 14 * pow(z, 2) + 1.0);
            break;
        case -1:
            SH = constants::c_165_256p * 2 * y * (21 * pow(z, 4) - 14 * pow(z, 2) + 1.0);
            break;
        case 2:
            SH = constants::c_1155_64p * (pow(x, 2) - pow(y, 2)) * (3 * pow(z, 3) - z);
            break;
        case -2:
            SH = constants::c_1155_64p * 2 * x * y * (3 * pow(z, 3) - z);
            break;
        case 3:
            SH = constants::c_385_512p * x * (pow(x, 2) - 3 * pow(y, 2)) * (9 * pow(z, 2) - 1.0);
            break;
        case -3:
            SH = constants::c_385_512p * y * (3 * pow(x, 2) - pow(y, 2)) * (9 * pow(z, 2) - 1.0);
            break;
        case 4:
            SH = constants::c_3465_256p * (pow(x, 4) - 6 * x * x * y * y + pow(y, 4)) * z;
            break;
        case -4:
            SH = -constants::c_3465_256p * (4 * x * pow(y, 3) - 4 * pow(x, 3) * y) * z;
            break;
        case 5:
            SH = constants::c_693_2048p * (2 * pow(x, 5) - 20 * pow(x, 3) * pow(y, 2) + 10 * x * pow(y, 4));
            break;
        case -5:
            SH = constants::c_693_2048p * (2 * pow(y, 5) - 20 * pow(x, 2) * pow(y, 3) + 10 * y * pow(x, 4));
            break;
        default:
            err_not_impl_f("Wrong spherical harmonic called!", std::cout);
        }
        break;
    default:
        err_not_impl_f("Higehr than l=4 not done for spherical harmonic!", std::cout);
    }
    return SH;
}