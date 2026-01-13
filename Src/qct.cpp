#include "pch.h"
#include "convenience.h"
#include "fchk.h"
#include "basis_set.h"
#include "wfn_class.h"
#include "b2c.h"
#include "bondwise_analysis.h"
#include "properties.h"

struct sel {
    int* selection;
    int* NbAtoms;
    int* MoleculeFiles;
    double* Radius;
    double* Increments;
    double* Cutplot;
    double Intermolecular;
    double Cutoffs[2];
    std::string Oname;
    bool* selection2;
    int NbFiles = 1;
    double Ligand[2];
    std::vector<int> ignore;
    int Frames;
    int Output;
    sel() {
        selection = (int*)malloc(sizeof(int) * 7);
        NbAtoms = nullptr;
        MoleculeFiles = nullptr;
        Radius = (double*)malloc(sizeof(double) * 4);
        Increments = (double*)malloc(sizeof(double) * 3);
        Cutplot = (double*)malloc(sizeof(double) * 2);
        Ligand[0] = 1;
        Ligand[1] = 4.0;
        Intermolecular = 0.9;
        Cutplot[0] = 0.05; Cutplot[1] = 0.5;
        Increments[0] = 0.1; Increments[1] = 0.1; Increments[2] = 0.1;
        Cutoffs[0] = 0.2; Cutoffs[1] = 1.0;
        Radius[0] = 0.0; Radius[1] = 0.0; Radius[2] = 0.0; Radius[3] = 3.0;
        Frames = 1;
        Oname = "out";
        Output = 3;
    }
};

sel menu(
    options& opt,
    const std::vector<WFN>& wavy,
    sel& res)
{
    using namespace std;
    bool run_o = false;
    res.NbAtoms = (int*)malloc(sizeof(int) * wavy.size());
    res.MoleculeFiles = (int*)malloc(sizeof(int) * wavy.size());
    res.NbAtoms[0] = wavy[0].get_ncen();
    res.MoleculeFiles[0] = 0;
    while (!run_o) {
        cout << "The following options are available (select preceding number to change them):" << endl;
        cout << "1) Select wavefunctions to use, the # of Files already loaded is: ";
        if (wavy.size() > 0) {
            cout << wavy.size() << endl;
            for (int i = 0; i < wavy.size(); i++) {
                cout << "   File " << i << ": " << wavy[i].get_path();
                if (res.MoleculeFiles[0] == i) cout << " *";
                if (wavy.size() > 1 && res.MoleculeFiles[i] == i) cout << " #";
                cout << endl;
            }
        }
        else cout << "No Files loaded yet!" << endl;
        if (res.NbFiles > 1) {
            cout << "2) Ligand is file # and radius (in A) around it: ";
            if (res.selection[1]) cout << res.Ligand[0] << " " << res.Ligand[1] << endl;
            else cout << "disabled" << endl;
        }
        cout << "3) Cutoff for intermolecular is: ";
        if (!res.selection[4]) cout << "disabled" << endl;
        else cout << res.Intermolecular << endl;
        cout << "4) The output name is: " << res.Oname << endl
            << "5) Level of output (1=.out only, 2= .cubes only, 3= all): " << res.Output << endl
            << "6) Point separation of the cubes is: " << res.Increments[0] << " " << res.Increments[1] << " " << res.Increments[2] << endl
            << "7) Cutplot: " << res.Cutplot[0] << " " << res.Cutplot[1] << endl
            << "8) Select radius around molecule. Current value: " << opt.properties.radius << endl
            << "9) select origin of the cube. Current value: ";
        if (!res.selection[6])
            cout << "default" << endl;
        else
            cout << "(" << res.Radius[0] << "," << res.Radius[1] << "," << res.Radius[2] << ")" << endl;
        cout << "10) Enable Calculation of ELF (Currently: ";
        if (opt.properties.elf) cout << "YES)";
        else cout << "NO)";
        cout << endl
            << "11) Enable ELI-D Calculation (Currently: ";
        if (opt.properties.eli) cout << "YES)";
        else cout << "NO)";
        cout << endl
            << "12) Enable Laplacian Calculation (Currently: ";
        if (opt.properties.lap) cout << "YES)";
        else cout << "NO)";
        cout << endl;
        if (res.NbFiles > 1) {
            cout << "13) Save Hirshfeld surface cube: ";
            if (opt.properties.hirsh) cout << "YES)";
            else cout << "NO)";
        }
        else cout << "--) Save Hirshfeld surface cube: ";
        cout << endl
            << "14) Enable ESP Calculation (Currently: ";
        if (opt.properties.esp) cout << "YES)";
        else cout << "NO)";
        /*    cout << endl;
            << "15) Enable Electric Field Calculation (Currently: ";
        if(opt.doef) cout << "YES)";
        else cout << "NO)";*/
        cout << endl
            << "16) Enable Rho Calculation (Currently: ";
        if (opt.properties.rho) cout << "YES)";
        else cout << "NO)";
        cout << endl
            << "17) Enable RDG Calculation (Currently: ";
        if (opt.properties.rdg) cout << "YES)";
        else cout << "NO)";
        cout << endl
            << "18) Enable deformation density Calculation (Currently: ";
        if (opt.properties.def) cout << "YES)";
        else cout << "NO)";
        cout << endl
            << "19) Select atoms to ignore during calculation: ";
        if (opt.ignore.size() == 0)
            cout << "None" << endl;
        else {
            for (int i = 0; i < opt.ignore.size(); i++)
                cout << std::setw(3) << opt.ignore[i];
            cout << endl;
        }
        cout << "20) Enable deformation hirshfeld density Calculation (Currently: ";
        if (opt.properties.hdef) cout << "YES)" << endl;
        else cout << "NO)" << endl;
        cout << "-1) End the program without a calculation" << endl
            << "0) Start the calculations based on the options selected above" << endl
            << "Your selection: ";
        int _sel;
        std::cin >> _sel;
        cls();
        switch (_sel) {
        case 0:
            run_o = true;
            break;
        case 1: {
            cout << "How many files do you want to include? ";
            int temp;
            cin >> temp;
            if (temp >= wavy.size()) {
                cout << "Sorry, only number of loaded files supported" << endl;
                break;
            }
            int temp2;
            for (int i = 0; i < temp; i++) {
                cout << "# of wavefunction to use as #" << i + 1 << " (from list above): ";
                cin >> temp2;
                temp2--;
                if (temp2 < 0 || temp2 >= wavy.size()) {
                    cout << "Invalid choice, try again!" << endl;
                    i--;
                    continue;
                }
                res.MoleculeFiles[i] = temp2;
            }
        }
              break;
        case 2: {
            bool run = false;
            int temp = 0;
            double temp_r = 0;
            while (!run) {
                cout << "Please select one of the following molecule(s) to be considered as the ligand:" << endl;
                for (int i = 0; i < res.NbFiles; i++)
                    cout << "   File " << res.MoleculeFiles[i] << ": " << wavy[i].get_path() << endl;
                cin >> temp;
                if (temp >= res.NbFiles || temp < 0) cout << "ERROR, invalid choice of molecule" << endl;
                else {
                    run = true;
                    res.Ligand[0] = (double)temp;
                }
            }
            run = false;
            while (!run) {
                cout << "select radius: ";
                cin >> temp_r;
                if (temp_r < 0) cout << "Negative radius is invalid, please try again" << endl;
                else if (temp_r > 20) cout << "Please stay realistic! Select smaller radius or ask the programmer to make the threshold bigger!" << endl;
                else {
                    run = true;
                    res.Ligand[1] = temp_r;
                }
            }
            res.selection[1] = true;
        }
              break;
        case 3: {
            bool run = false;
            float temp = 0;
            while (!run) {
                cout << "Please select the cutoff for intermolecular:" << endl;
                cin >> temp;
                if (temp < 0) cout << "ERROR, invalid choice of Intermolecular cutoff" << endl;
                else if (temp > 20) cout << "Please stay realistic! Select smaller radius or ask the programmer to make the threshold bigger!" << endl;
                else {
                    run = true;
                    res.Intermolecular = temp;
                }
            }
            res.selection[4] = true;
            continue;
        }
              break;
        case 4: {
            bool run = false;
            std::filesystem::path temp;
            while (!run) {
                cout << "Please select the name for the ouput files:" << endl;
                cin >> temp;
                std::filesystem::path temp2;
                temp2 = temp;
                temp2.replace_extension(".out");
                while (exists(temp2)) {
                    cout << "File " << temp2 << " already exists, do you want to overwrite it? ";
                    if (!yesno()) {
                        cout << "then enter an alternative name: ";
                        cin >> temp;
                    }
                }
                run = true;
                res.Oname = temp.string();
            }
            continue;
        }
              break;
        case 5: {
            bool run = false;
            int temp;
            while (!run) {
                cout << "select output level: ";
                cin >> temp;
                if (temp < 4 && temp>0) {
                    res.Output = temp;
                    run = true;
                }
                else cout << "invalid selection!" << endl;
            }
            continue;
        }
              break;
        case 6: {
            bool run = false;
            float temp[3];
            while (!run) {
                cout << "Please enter the separation of points in all three dimensions, starting with X: ";
                cin >> temp[0];
                cout << "Now y: ";
                cin >> temp[1];
                cout << "Now z: ";
                cin >> temp[2];
                run = true;
                for (int i = 0; i < 3; i++) {
                    if (temp[i] > 0 && temp[i] < 1) {
                        if (run)
                            res.Increments[i] = temp[i];
                        else {
                            cout << "invalid selection in dimension " << i << endl;
                            run = false;
                        }
                    }
                }
            }
            continue;
        }
              break;
        case 7: {
            bool run = false;
            while (!run) {
                cout << "Select Cutoff for plots:" << endl << "Low cutof: ";
                float temp[2];
                cin >> temp[0];
                if (temp[0] < 0 || temp[0] > 2)
                    cout << "invalid choice, cutoff must be between 0 and 2, try again!" << endl;
                else {
                    cout << "High Cutoff: ";
                    cin >> temp[1];
                    if (temp[1] < temp[0] || temp[1] > 2)
                        cout << "invalid choice, cutoff must be bigger than low cutoff and smaller than 2, try again!" << endl;
                    else {
                        res.Cutplot[0] = temp[0];
                        res.Cutplot[1] = temp[1];
                        run = true;
                    }
                }
            }
            continue;
        }
              break;
        case 10: {
            if (opt.properties.eli) opt.properties.eli = false;
            if (opt.properties.elf) opt.properties.elf = false;
            else opt.properties.elf = true;

        }
               break;
        case 11: {
            if (opt.properties.elf) opt.properties.elf = false;
            if (opt.properties.eli) opt.properties.eli = false;
            else opt.properties.eli = true;
        }
               break;
        case 12: {
            opt.properties.lap = !opt.properties.lap;
        }
               break;
        case 9: {
            cout << "Select starting point:" << endl;
            float temp[3];
            cout << "X: ";
            cin >> temp[0];
            cout << "Y: ";
            cin >> temp[1];
            cout << "Z: ";
            cin >> temp[2];
            res.Radius[0] = temp[0];
            res.Radius[1] = temp[1];
            res.Radius[2] = temp[2];
            res.selection[6] = true;
        }
              break;
        case 8: {
            bool run = false;
            while (!run) {
                cout << "Select positive radius: ";
                double temp;
                cin >> temp;
                if (temp <= 0.0 || temp > 100.0)
                    cout << "invalid choice!" << endl;
                else {
                    opt.properties.radius = temp;
                    res.selection[2] = true;
                    break;
                }
            }
        }
              break;
        case 13: {
            if (res.NbFiles > 1)
                opt.properties.hirsh = !opt.properties.hirsh;
        }
               break;
        case 14: {
            opt.properties.esp = !opt.properties.esp;
        }
               break;
               /*    case 15: {
                       opt.doef = !opt.doef;
                   }
                   break;*/
        case 16: {
            if (opt.properties.hdef)
                opt.properties.rho = true;
            opt.properties.rho = !opt.properties.rho;
        }
               break;
        case 17: {
            opt.properties.rdg = !opt.properties.rdg;
            if (opt.properties.rdg && !opt.properties.rho) opt.properties.rho = true;
        }
               break;
        case 18: {
            opt.properties.def = !opt.properties.def;
            if (opt.properties.def && !opt.properties.rho) opt.properties.rho = true;
        }
               break;
        case 19: {
            cls();
            while (true) {
                cout << "Which atoms to exlude? (marked with *)" << endl;
                for (int a = 0; a < wavy[res.MoleculeFiles[0]].get_ncen(); a++) {
                    cout << setw(4) << a << ") " << setw(4) << constants::atnr2letter(wavy[res.MoleculeFiles[0]].get_atom_charge(a))
                        << setw(14) << scientific << wavy[res.MoleculeFiles[0]].get_atom_coordinate(a, 0)
                        << setw(14) << scientific << wavy[res.MoleculeFiles[0]].get_atom_coordinate(a, 1)
                        << setw(14) << scientific << wavy[res.MoleculeFiles[0]].get_atom_coordinate(a, 2);
                    for (int i = 0; i < opt.ignore.size(); i++)
                        if (opt.ignore[i] == a)
                            cout << "*";
                    cout << endl;
                }
                cout << " -10) Accept selection and return to previous menu" << endl
                    << "Select atom number: ";
                int input;
                cin >> input;
                if (input == -10)
                    break;
                if (input >= wavy[res.MoleculeFiles[0]].get_ncen() || input < 0) {
                    cls();
                    cout << "Index out of range!" << endl;
                    continue;
                }
                if (opt.ignore.size() == 0) {
                    opt.ignore.push_back(input);
                    cls();
                    continue;
                }
                for (int i = 0; i < opt.ignore.size(); i++) {
                    if (input == opt.ignore[i]) {
                        opt.ignore.erase(opt.ignore.begin() + i);
                        cls();
                        break;
                    }
                    else if (i == opt.ignore.size() - 1) {
                        opt.ignore.push_back(input);
                        cls();
                        break;
                    }
                }
            }
        }
               break;
        case 20: {
            opt.properties.hdef = !opt.properties.hdef;
            if (opt.properties.hdef) {
                cout << "Which atom to use?" << endl;
                for (int a = 0; a < wavy[res.MoleculeFiles[0]].get_ncen(); a++)
                    cout << setw(4) << a << ") " << setw(4) << constants::atnr2letter(wavy[res.MoleculeFiles[0]].get_atom_charge(a))
                    << setw(14) << scientific << wavy[res.MoleculeFiles[0]].get_atom_coordinate(a, 0)
                    << setw(14) << scientific << wavy[res.MoleculeFiles[0]].get_atom_coordinate(a, 1)
                    << setw(14) << scientific << wavy[res.MoleculeFiles[0]].get_atom_coordinate(a, 2) << endl;
                int input;
                cin >> input;
                if (input == -10)
                    break;
                else if (input >= wavy[res.MoleculeFiles[0]].get_ncen() || input < 0) {
                    cls();
                    cout << "Index out of range!" << endl;
                    continue;
                }
                else
                    for (int a = 0; a < wavy[res.MoleculeFiles[0]].get_ncen(); a++)
                        if (a != input)
                            opt.ignore.push_back(a);
                cout << "Select positive radius: ";
                double temp;
                cin >> temp;
                if (temp <= 0.0 || temp > 100.0)
                    cout << "invalid choice!" << endl;
                else {
                    res.Radius[3] = temp;
                    res.selection[2] = true;
                }
                opt.properties.rho = true;
            }
        }
               break;
        case -1:
            return res;
        default:
            cls();
            cout << "Invalid selection!" << endl;
            continue;
        }
    }
    if (!opt.properties.rho) {
        if (opt.properties.rdg) {
            cout << "MUST calcualate Rho when using EDG, enabling RHO" << endl;
            opt.properties.rho = true;
        }
        if (opt.properties.def)
        {
            cout << "MUST calcualate Rho when using deformation density, enabling RHO" << endl;
            opt.properties.rho = true;
        }
    }
    return res;
}

int acu_nci(std::vector<WFN>& wavy, options& opt) {
    using namespace std;
    vector<std::filesystem::path> names;
    for (int i = 0; i < wavy.size(); i++)
        names.push_back(wavy[i].get_path());

    bool promolecular = false; //whether to calculate with a wavefunction or the promolecular density

    sel run;
    run.Oname = wavy[0].get_path().string();

    /*execute the interactive menu*/
    run = menu(opt, wavy, run);

    if (run.Frames != 1 || run.NbFiles > 1) promolecular = true;


    cube* persistant_cube_rho;
    cube* persistant_cube_RDG;
    cube* persistant_cube_Elf;
    cube* persistant_cube_Eli;
    cube* persistant_cube_Lap;
    cube* persistant_cube_ESP;
    //cube* persistant_cube_EF;
    //cube* persistant_cube_Hirsh;
    //cube *persistant_cube_def;

    if (run.NbFiles == 1)
        run.Intermolecular = 1;

    int counter_100 = 0;
    vector<vector<vector<int> > > sign_counter;
    vector<vector<vector<int> > > ignore_RDG;
    int negative_signs = 0;
    properties_options opts = opt.properties;
    opts.resolution = run.Increments[0];

    readxyzMinMax_fromWFN(
        wavy[run.MoleculeFiles[0]],
        opts,
        false);
    run.Oname = wavy[run.MoleculeFiles[0]].get_path().filename().string();
    cube CubeRho(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.rho || opts.eli || opts.lap),
        CubeDEF(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.def),
        CubeRDG(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.rdg),
        CubeElf(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.elf),
        CubeEli(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.eli),
        CubeLap(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.lap),
        CubeESP(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.esp),
        CubeHDEF(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.hdef);
    CubeRho.give_parent_wfn(wavy[run.MoleculeFiles[0]]);
    CubeDEF.give_parent_wfn(wavy[run.MoleculeFiles[0]]);
    CubeRDG.give_parent_wfn(wavy[run.MoleculeFiles[0]]);
    CubeElf.give_parent_wfn(wavy[run.MoleculeFiles[0]]);
    CubeEli.give_parent_wfn(wavy[run.MoleculeFiles[0]]);
    CubeLap.give_parent_wfn(wavy[run.MoleculeFiles[0]]);
    CubeESP.give_parent_wfn(wavy[run.MoleculeFiles[0]]);
    CubeHDEF.give_parent_wfn(wavy[run.MoleculeFiles[0]]);

    //CubeHirsh(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[opt.MoleculeFiles[0]].get_ncen(), opt.dohirsh);

    persistant_cube_rho = new cube(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.rho && run.Frames > 1);
    persistant_cube_RDG = new cube(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.rdg && run.Frames > 1);
    persistant_cube_Elf = new cube(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.elf && run.Frames > 1);
    persistant_cube_Eli = new cube(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.eli && run.Frames > 1);
    persistant_cube_Lap = new cube(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.lap && run.Frames > 1);
    //persistant_cube_Hirsh = new cube(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[opt.MoleculeFiles[0]].get_ncen(), opt.dohirsh);
    persistant_cube_ESP = new cube(opts.NbSteps, wavy[run.MoleculeFiles[0]].get_ncen(), opts.esp && run.Frames > 1);
    //persistant_cube_EF = new cube(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[opt.MoleculeFiles[0]].get_ncen(), opt.doef);
    //persistant_cube_def = new cube(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[opt.MoleculeFiles[0]].get_ncen(), opt.dodef);

    for (int i = 0; i < 3; i++) {
        CubeRho.set_origin(i, opts.MinMax[i]);
        CubeRDG.set_origin(i, opts.MinMax[i]);
        CubeElf.set_origin(i, opts.MinMax[i]);
        CubeEli.set_origin(i, opts.MinMax[i]);
        CubeLap.set_origin(i, opts.MinMax[i]);
        CubeESP.set_origin(i, opts.MinMax[i]);
        CubeHDEF.set_origin(i, opts.MinMax[i]);
        CubeDEF.set_origin(i, opts.MinMax[i]);
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                CubeRho.set_vector(i, j, run.Increments[i]);
                CubeRDG.set_vector(i, j, run.Increments[i]);
                CubeElf.set_vector(i, j, run.Increments[i]);
                CubeEli.set_vector(i, j, run.Increments[i]);
                CubeLap.set_vector(i, j, run.Increments[i]);
                CubeESP.set_vector(i, j, run.Increments[i]);
                CubeHDEF.set_vector(i, j, run.Increments[i]);
                CubeDEF.set_vector(i, j, run.Increments[i]);
            }
            else {
                CubeRho.set_vector(i, j, 0);
                CubeRDG.set_vector(i, j, 0);
                CubeElf.set_vector(i, j, 0);
                CubeEli.set_vector(i, j, 0);
                CubeLap.set_vector(i, j, 0);
                CubeESP.set_vector(i, j, 0);
                CubeHDEF.set_vector(i, j, 0);
                CubeDEF.set_vector(i, j, 0);
            }
        }
    }
    CubeRho.set_comment1("Calculated density using QCT");
    CubeRDG.set_comment1("Calculated reduced density gradient using QCT");
    CubeElf.set_comment1("Calculated electron localization function using QCT");
    CubeEli.set_comment1("Calculated same-spin electron localizability indicator using QCT");
    CubeLap.set_comment1("Calculated laplacian of electron density using QCT");
    CubeESP.set_comment1("Calculated electrostatic potential using QCT");
    CubeHDEF.set_comment1("Calculated Hirshfeld deformation density using QCT");
    CubeRho.set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
    CubeRDG.set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
    CubeElf.set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
    CubeLap.set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
    CubeESP.set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
    CubeHDEF.set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
    CubeRho.set_path(run.Oname + "_rho.cube");
    CubeRDG.set_path(run.Oname + "_rdg.cube");
    CubeElf.set_path(run.Oname + "_elf.cube");
    CubeEli.set_path(run.Oname + "_eli.cube");
    CubeLap.set_path(run.Oname + "_lap.cube");
    CubeESP.set_path(run.Oname + "_esp.cube");
    CubeHDEF.set_path(run.Oname + "_hdef.cube");
    CubeDEF.set_path(run.Oname + "_def.cube");

    for (int iframe = 0; iframe < run.Frames; iframe++)
    {
        std::cout << "  *                                                                       *" << std::endl;
        std::cout << "  *   Working on Frame          : " << std::setw(5) << iframe + 1 << " /" << std::setw(5) << run.Frames << "                            * " << std::endl;
        /*Read file(s) .xyz and find min and max*/

        string filenames[2];

        if (iframe == 0) { //print Header
            std::cout << "\n          ___________________________________________________________" << std::endl;
            std::cout << "  *                                                                       *" << std::endl;
            std::cout << "  *   NbFiles                   : " << std::setw(1) << run.NbFiles << "                                       *  " << std::endl;
            for (int i = 0; i < run.NbFiles; i++) {
#ifdef _WIN32
                std::wcout << L"  *   MoleculeFile[" << std::setw(2) << i << L"]          : " << std::left << std::setw(20) << wavy[run.MoleculeFiles[i]].get_path().filename().c_str() << L" / " << std::setw(5) << run.NbAtoms[i] << L" atoms      *" << std::endl;
#else
                std::cout << "  *   MoleculeFile[" << std::setw(2) << i << "]          : " << std::left << std::setw(20) << wavy[run.MoleculeFiles[i]].get_path().filename().string() << " / " << std::setw(5) << run.NbAtoms[i] << " atoms      *" << std::endl;
#endif
            }
            std::cout << "  *   OutPut filename Prefix    : " << std::setw(20) << run.Oname << "                    *" << std::endl;
            std::cout << "  *                                                                       *" << std::endl;

            std::cout << "  *   gridBox Min               : " << std::fixed << std::setprecision(6)
                << std::setw(11) << opts.MinMax[0] << " "
                << std::setw(11) << opts.MinMax[1] << " "
                << std::setw(11) << opts.MinMax[2] << "     *" << std::endl;
            std::cout << "  *   gridBox Max               : " << std::fixed << std::setprecision(6)
                << std::setw(11) << opts.MinMax[3] << " "
                << std::setw(11) << opts.MinMax[4] << " "
                << std::setw(11) << opts.MinMax[5] << "     *" << std::endl;
            std::cout << "  *   Increments(bohr)          : " << std::fixed << std::setprecision(6)
                << std::setw(11) << run.Increments[0] << " "
                << std::setw(11) << run.Increments[1] << " "
                << std::setw(11) << run.Increments[2] << "     *" << std::endl;
            std::cout << "  *   NbSteps                   : " << std::setw(11) << opts.NbSteps[0] << " "
                << std::setw(11) << opts.NbSteps[1] << " "
                << std::setw(11) << opts.NbSteps[2] << "     *" << std::endl;
            std::cout << "  *                                                                       * " << std::endl;
            if (run.NbFiles == 2)
                std::cout << "  *   Intermolecular            :" << std::fixed << std::setprecision(2) << std::setw(5) << run.Intermolecular << "                                    * " << std::endl;
            if (!promolecular) {
                std::cout << "  *   Number of primitives      :     " << std::setw(5) << wavy[run.MoleculeFiles[0]].get_nex() << "                               * " << std::endl;
                std::cout << "  *   Number of MOs             :     " << std::setw(5) << wavy[run.MoleculeFiles[0]].get_nmo() << "                               * " << std::endl;
            }
            if (run.Frames > 1) {
                ignore_RDG.resize(opts.NbSteps[0]);
                sign_counter.resize(opts.NbSteps[0]);
                for (int i = 0; i < opts.NbSteps[0]; i++) {
                    ignore_RDG[i].resize(opts.NbSteps[1]);
                    sign_counter[i].resize(opts.NbSteps[1]);
                    for (int j = 0; j < opts.NbSteps[1]; j++) {
                        ignore_RDG[i][j].resize(opts.NbSteps[2]);
                        sign_counter[i][j].resize(opts.NbSteps[2]);
                    }
                }
                for (int i = 0; i < 3; i++) {
                    persistant_cube_rho[0].set_origin(i, opts.MinMax[i]);
                    persistant_cube_RDG[0].set_origin(i, opts.MinMax[i]);
                    persistant_cube_Elf[0].set_origin(i, opts.MinMax[i]);
                    persistant_cube_Eli[0].set_origin(i, opts.MinMax[i]);
                    persistant_cube_Lap[0].set_origin(i, opts.MinMax[i]);
                    persistant_cube_ESP[0].set_origin(i, opts.MinMax[i]);
                    for (int j = 0; j < 3; j++) {
                        if (i == j) {
                            persistant_cube_rho[0].set_vector(i, j, run.Increments[i]);
                            persistant_cube_RDG[0].set_vector(i, j, run.Increments[i]);
                            persistant_cube_Elf[0].set_vector(i, j, run.Increments[i]);
                            persistant_cube_Eli[0].set_vector(i, j, run.Increments[i]);
                            persistant_cube_Lap[0].set_vector(i, j, run.Increments[i]);
                            persistant_cube_ESP[0].set_vector(i, j, run.Increments[i]);
                        }
                        else {
                            persistant_cube_rho[0].set_vector(i, j, 0);
                            persistant_cube_RDG[0].set_vector(i, j, 0);
                            persistant_cube_Elf[0].set_vector(i, j, 0);
                            persistant_cube_Eli[0].set_vector(i, j, 0);
                            persistant_cube_Lap[0].set_vector(i, j, 0);
                            persistant_cube_ESP[0].set_vector(i, j, 0);
                        }
                    }
                }
                persistant_cube_rho[0].set_comment1("Calculated density using QCT");
                persistant_cube_RDG[0].set_comment1("Calculated reduced density gradient using QCT");
                persistant_cube_Elf[0].set_comment1("Calculated electron localization function using QCT");
                persistant_cube_Eli[0].set_comment1("Calculated same-spin electron localizability indicator using QCT");
                persistant_cube_Lap[0].set_comment1("Calculated laplacian of electron density using QCT");
                persistant_cube_ESP[0].set_comment1("Calculated electrostatic potential using QCT");
                persistant_cube_rho[0].set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
                persistant_cube_RDG[0].set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
                persistant_cube_Elf[0].set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
                persistant_cube_Lap[0].set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
                persistant_cube_ESP[0].set_comment2("from " + wavy[run.MoleculeFiles[0]].get_path().string());
                persistant_cube_rho[0].set_path(run.Oname + "_rho.cube");
                persistant_cube_RDG[0].set_path(run.Oname + "_rdg.cube");
                persistant_cube_Elf[0].set_path(run.Oname + "_elf.cube");
                persistant_cube_Eli[0].set_path(run.Oname + "_eli.cube");
                persistant_cube_Lap[0].set_path(run.Oname + "_lap.cube");
                persistant_cube_ESP[0].set_path(run.Oname + "_esp.cube");
                for (int x = 0; x < opts.NbSteps[0]; x++)
                    for (int y = 0; y < opts.NbSteps[1]; y++)
                        for (int z = 0; z < opts.NbSteps[2]; z++) {
                            if (opts.rho) persistant_cube_rho[0].set_value(x, y, z, 0.0);
                            if (opts.rdg) persistant_cube_RDG[0].set_value(x, y, z, 0.0);
                            if (opts.elf || opts.eli) persistant_cube_Elf[0].set_value(x, y, z, 0.0);
                            if (opts.lap) persistant_cube_Lap[0].set_value(x, y, z, 0.0);
                            //if (opts.dohirsh) persistant_cube_Hirsh[0].set_value(x, y, z, 0.0);
                            if (opts.esp) persistant_cube_ESP[0].set_value(x, y, z, 0.0);
                            //if (opts.doef) persistant_cube_EF[0].set_value(x, y, z, 0.0);
                            //if (opts.dodef) persistant_cube_def[0].set_value(x, y, z, 0.0);
                            ignore_RDG[x][y][z] = 0;
                            sign_counter[x][y][z] = 0;
                        }
            }
        }
        if (opt.ignore.size() > 0 && opts.rho && !opts.hdef) {
            Calc_Rho(
                CubeRho,
                wavy[run.MoleculeFiles[0]],
                opts.radius,
                std::cout,
                false
            );
        }
        else if (opts.rho || opts.rdg || opts.elf || opts.eli || opts.lap || opts.def) {
            Calc_Prop(
                CubeRho,
                CubeRDG,
                CubeElf,
                CubeEli,
                CubeLap,
                CubeESP,
                wavy[run.MoleculeFiles[0]],
                opts.radius,
                std::cout,
                false,
                false
            );
        }

        if (opts.def) {
            Calc_Static_Def(
                CubeDEF,
                CubeRho,
                wavy[0],
                opts.radius,
                std::cout,
                false
            );
        }

        if (opts.hdef)
            Calc_Hirshfeld(
                CubeHDEF,
                CubeRho,
                wavy[run.MoleculeFiles[0]],
                opts.radius,
                opt.ignore[0],
                std::cout,
                false
            );

        if (opts.esp)
            Calc_ESP(
                CubeESP,
                wavy[run.MoleculeFiles[0]],
                opts.radius,
                opt.no_date,
                std::cout,
                false
            );
        int z_counter = 0;
        if (run.Frames > 1) {
            for (int x = 0; x < opts.NbSteps[0]; x++)
                for (int y = 0; y < opts.NbSteps[1]; y++)
                    for (int z = 0; z < opts.NbSteps[2]; z++) {
                        if (opts.rho) {
                            if (CubeRho.get_value(x, y, z) < 0) {
                                persistant_cube_rho[0].set_value(x, y, z, persistant_cube_rho[0].get_value(x, y, z) - CubeRho.get_value(x, y, z));
                                sign_counter[x][y][z]--;
                            }
                            else {
                                if (opt.debug && CubeRho.get_value(x, y, z) == 0.0) z_counter++;
                                persistant_cube_rho[0].set_value(x, y, z, persistant_cube_rho[0].get_value(x, y, z) + CubeRho.get_value(x, y, z));
                                sign_counter[x][y][z]++;
                            }
                        }
                        if (opts.rdg) {
                            if (CubeRDG.get_value(x, y, z) < 0)
                                ignore_RDG[x][y][z]++;
                            persistant_cube_RDG[0].set_value(x, y, z, persistant_cube_RDG[0].get_value(x, y, z) + abs(CubeRDG.get_value(x, y, z)));
                        }
                        if (opts.elf || opts.eli) persistant_cube_Elf[0].set_value(x, y, z, persistant_cube_Elf[0].get_value(x, y, z) + CubeElf.get_value(x, y, z));
                        if (opts.lap) persistant_cube_Lap[0].set_value(x, y, z, persistant_cube_Lap[0].get_value(x, y, z) + CubeLap.get_value(x, y, z));
                        if (opts.esp) persistant_cube_ESP[0].set_value(x, y, z, persistant_cube_ESP[0].get_value(x, y, z) + CubeESP.get_value(x, y, z));
                        //if (opts.dodef) persistant_cube_def[i] += CubeDEF[i];
                    }
            if (opt.debug) std::cout << "\nFinished frame " << iframe << ", count_100= " << counter_100 << ", count_0= " << z_counter << "\n";
        }
    } //END FRAMES

    if (opt.debug) std::cout << "Finished with all frames!\n";
    if (run.Frames > 1) {
        if (opt.debug) cout << "Now making averages and RDG from calculated rho and gradients\n" << endl;
        for (int x = 0; x < opts.NbSteps[0]; x++)
            for (int y = 0; y < opts.NbSteps[1]; y++)
                for (int z = 0; z < opts.NbSteps[2]; z++) {
                    if (opts.rho) {
                        if (sign_counter[x][y][z] >= 0)
                            persistant_cube_rho[0].set_value(x, y, z, persistant_cube_rho[0].get_value(x, y, z) / run.Frames);
                        else {
                            persistant_cube_rho[0].set_value(x, y, z, persistant_cube_rho[0].get_value(x, y, z) / -run.Frames);
                            negative_signs++;
                        }
                    }
                    if (opts.rdg) {
                        if (ignore_RDG[x][y][z] != run.Frames)
                            persistant_cube_RDG[0].set_value(x, y, z, persistant_cube_RDG[0].get_value(x, y, z) / (run.Frames - ignore_RDG[x][y][z]));
                        else
                            persistant_cube_RDG[0].set_value(x, y, z, 100.0);
                    }
                    if (opts.elf || opts.eli) persistant_cube_Elf[0].set_value(x, y, z, persistant_cube_Elf[0].get_value(x, y, z) / run.Frames);
                    if (opts.lap) persistant_cube_Lap[0].set_value(x, y, z, persistant_cube_Lap[0].get_value(x, y, z) / run.Frames);
                    //if (opts.dohirsh) persistant_cube_Hirsh[i] /= opt.Frames;
                    if (opts.esp)persistant_cube_ESP[0].set_value(x, y, z, persistant_cube_ESP[0].get_value(x, y, z) / run.Frames);
                    //if (opts.doef) persistant_cube_EF[i] /= opt.Frames;
                    //if (opts.dodef) persistant_cube_def[i] /= opt.Frames;
                }
        if (opt.debug) std::cout << "\n\n100_counter= " << counter_100 << ", points in vector= " << opts.n_grid_points() << ", frames= " << run.Frames << ", negative signs= " << negative_signs << "\n";
        if (opts.rdg)
            for (int x = 0; x < opts.NbSteps[0]; x++)
                for (int y = 0; y < opts.NbSteps[1]; y++)
                    for (int z = 0; z < opts.NbSteps[2]; z++)
                        if (persistant_cube_RDG[0].get_value(x, y, z) != 0.0 && persistant_cube_rho[0].get_value(x, y, z) != 0.0)
                            persistant_cube_RDG[0].set_value(x, y, z, persistant_cube_RDG[0].get_value(x, y, z) / pow(abs(persistant_cube_rho[0].get_value(x, y, z)), (double)1.3333333333333333333333));
    }
    else
        if (opts.rdg)
            for (int x = 0; x < opts.NbSteps[0]; x++)
                for (int y = 0; y < opts.NbSteps[1]; y++)
                    for (int z = 0; z < opts.NbSteps[2]; z++)
                        if (CubeRDG.get_value(x, y, z) != 0.0 && CubeRho.get_value(x, y, z) != 0.0)
                            CubeRDG.set_value(x, y, z, CubeRDG.get_value(x, y, z) / pow(abs(CubeRho.get_value(x, y, z)), (double)1.3333333333333333333333));

    /*if (opt.Output == 1 || opt.Output == 3)
    {
        if(opt.dorho && opt.dordg){
            std::cout << "  *                                                                       *\n";
            std::cout << "  *   Writing .dat file ...                                               *\n";
            if(opt.Frames>1)
                outdat(opt.Oname.c_str(),
                        opt.Cutoffs,
                        opt.NbSteps,
                        persistant_cube_rho[0],
                        persistant_cube_RDG[0]);
            else
                outdat(opt.Oname.c_str(),
                        opt.Cutoffs,
                        opt.NbSteps,
                        CubeRho,
                        CubeRDG);
        }
    }*/
    if (run.Output == 2 || run.Output == 3)
    {
        std::cout << "  *                                                                       *\n";
        std::cout << "  *   Writing .cube files ...                                             *\n";
        if (run.Frames > 1) {
            if (opts.rho && !opts.rdg) {
                persistant_cube_rho[0].set_path(run.Oname + "_rho.cube");
                persistant_cube_rho[0].write_file(true, true);
            }
            if (opts.rdg) {
                persistant_cube_rho[0].set_path(run.Oname + "_signed_rho.cube");
                persistant_cube_rho[0].write_file(true);
                persistant_cube_rho[0].set_path(run.Oname + "_rho.cube");
                persistant_cube_rho[0].write_file(true, true);
                persistant_cube_RDG[0].write_file(true);
            }
            if (opts.elf || opts.eli) persistant_cube_Elf[0].write_file(true);
            if (opts.lap) persistant_cube_Lap[0].write_file(true);
            if (opts.esp) persistant_cube_ESP[0].write_file(true);
            //if (opts.doef) outCubeEF(opts.Oname.c_str(),
            //if (opts.dodef) outCubeDEF(opts.Oname.c_str(),
            //if (opts.dohirsh) outCubeHirsh(opts.Oname.c_str(),
        }
        else {
            if (opts.rho && !opts.rdg) {
                CubeRho.set_path(run.Oname + "_rho.cube");
                CubeRho.write_file(true, true);
            }
            if (opts.rdg) {
                CubeRho.set_path(run.Oname + "_signed_rho.cube");
                CubeRho.write_file(true);
                CubeRho.set_path(run.Oname + "_rho.cube");
                CubeRho.write_file(true, true);
                CubeRDG.write_file(true);
            }
            if (opts.elf) CubeElf.write_file(true);
            if (opts.eli) CubeEli.write_file(true);
            if (opts.lap) CubeLap.write_file(true);
            if (opts.esp) CubeESP.write_file(true);
            if (opts.hdef) CubeHDEF.write_file(true);
            if (opts.def) CubeDEF.write_file(true);
            //if(opt.dohirsh) outCubeHirsh(opt.Oname.c_str(),
            //if(opt.doef) outCubeEF(opt.Oname.c_str(),
        }
        if (opts.rdg) wavy[run.MoleculeFiles[0]].push_back_cube(run.Oname + "_rdg.cube", false, false);
        if (opts.rho) {
            wavy[run.MoleculeFiles[0]].push_back_cube(run.Oname + "_rho.cube", false, false);
            if (opts.rdg)
                wavy[run.MoleculeFiles[0]].push_back_cube(run.Oname + "_signed_rho.cube", false, false);
        }
        if (opts.elf) wavy[run.MoleculeFiles[0]].push_back_cube(run.Oname + "_elf.cube", false, false);
        else if (opts.eli) wavy[run.MoleculeFiles[0]].push_back_cube(run.Oname + "_eli.cube", false, false);
        if (opts.lap) wavy[run.MoleculeFiles[0]].push_back_cube(run.Oname + "_lap.cube", false, false);
        //if(opt.dohirsh) wavy[opt.MoleculeFiles[0]].push_back_cube(opt.Oname + "-hirshfeld.cube", false, false);
        if (opts.esp) wavy[run.MoleculeFiles[0]].push_back_cube(run.Oname + "_esp.cube", false, false);
        //if(opt.doef) wavy[opt.MoleculeFiles[0]].push_back_cube(opt.Oname + "-ef_x.cube", false, false);
        //if(opt.doef) wavy[opt.MoleculeFiles[0]].push_back_cube(opt.Oname + "-ef_y.cube", false, false);
        //if(opt.doef) wavy[opt.MoleculeFiles[0]].push_back_cube(opt.Oname + "-ef_z.cube", false, false);
        //if(opt.doef) wavy[opt.MoleculeFiles[0]].push_back_cube(opt.Oname + "-ef_amp.cube", false, false);
        //if(opt.dodef) wavy[opt.MoleculeFiles[0]].push_back_cube(opt.Oname + "-def.cube", false, false);
    }

    std::cout << "  *                                                                       *" << std::endl;
    std::cout << "  *                                                                       *" << std::endl;
    std::cout << "          __________________________________________________________" << std::endl << std::endl << std::endl;
    return 0;
}

int QCT(options& opt, std::vector<WFN>& wavy)
{
    using namespace std;
    bool expert = false;
    bool end = false;
    int activewave = 0;
    string wavename;
    while (!end) {
        char sel;
        std::cout << "This Executable was built on: " + string(__DATE__) + " " + string(__TIME__) + "\n";
        std::cout << "What do you want to do?\n";
        if (wavy.size() > 0) {
            std::cout << "Active Wavefunction: " << wavy[activewave].get_path();
            string temp = wavy[activewave].get_basis_set_name();
            if (temp.length() > 2)std::cout << " (" << temp << ") ";
            else std::cout << " (no basis set loaded)";
            if (opt.debug)std::cout << temp;
            if (wavy[activewave].get_modified())std::cout << "*";
            if (expert)std::cout << " EXPERT MODE!";
            std::cout << endl;
        }
        std::cout << "R) Read in a new file" << endl;
        if (wavy.size() > 0) {
            std::cout << "M) Modify an already loaded wavefunction" << endl
                << "S) Save the active wavefunction";
            if (wavy[activewave].get_cube_count() > 0)
                std::cout << "or cube(s)" << endl;
            else std::cout << endl;
            std::cout << "O) Sort the exponents in the wavefunction" << endl
                << "N) Start NCI/Cube calculation plugin" << endl
                << "B) Read a basis set" << endl
                << "U) Check the unit of the atom positions" << endl;
            if (wavy.size() > 1)std::cout << "A) Activate another wavefunction" << endl;
            else std::cout << "-) Activate another wavefunction" << endl;
            if (wavy[activewave].get_cube_count() > 0)std::cout << "C) Work with cube files loaded" << endl;
            else std::cout << "-) Work with cube files loaded" << endl;
        }
        else {
            std::cout << "-) Modify an already loaded wavefunction" << endl
                << "-) Save the active wavefunction" << endl
                << "-) Sort the exponents in the wavefunction" << endl
                << "-) Start NCI/Cube calculation plugin" << endl
                << "-) Read a basis set" << endl
                << "-) Check the unit of the atom positions" << endl;
            if (wavy.size() > 1)std::cout << "A) Activate another wavefunction" << endl;
            else std::cout << "-) Activate another wavefunction" << endl;
            if (wavy.size() > 0 && wavy[activewave].get_cube_count() > 0)std::cout << "C) Work with cube files loaded" << endl;
            else std::cout << "-) Work with cube files loaded" << endl;
        }
        std::cout << "E) Toggle Expert mode (Disable assumptions)" << endl
            << "L) Limit the number of CPUs being used. Current value (-1 corresponds to all): " << opt.threads << endl
            << "Q) Quit the program" << endl;
        std::cin >> sel;
        cls();
        vector < vector <unsigned int> > selection;
        switch (sel) {
        case 'c':
        case 'C': {
            unsigned int nr_cubes = 0;
            selection.resize(2);
            for (int w = 0; w < wavy.size(); w++) for (int c = 0; c < wavy[w].get_cube_count(); c++) nr_cubes++;
            std::cout << "What do you want to do?" << endl
                << "1) Perform mathematic operation on one cube" << endl;
            if (nr_cubes > 1)std::cout << "2) Perform mathematic operation on two cubes" << endl;
            std::cout << "3) Analyze a cube" << endl
                << "0) Go back to main menu" << endl;
            int temp;
            std::cin >> temp;
            cls();
            switch (temp) {
            case 1: {
                if (nr_cubes <= 0) {
                    std::cout << "No cubes loaded!" << endl;
                    break;
                }
                std::cout << "Which operation do you want to perform?" << endl
                    << "1) Integrate values inside a cube" << endl
                    << "2) Periodically replicate a cube" << endl
                    << "3) Apply a threshhold to a cube" << endl
                    << "0) Get back to previous menu" << endl;
                int _sel = 0;
                std::cin >> _sel;
                selection[0].resize(1);
                selection[1].resize(1);
                if (nr_cubes >= 1 && _sel != 0)    select_cubes(selection, wavy, 1);
                if (wavy[selection[0][0]].get_cube_loaded(selection[1][0]) == false) {
                    if (opt.debug)std::cout << "Loading full file now!" << endl;
                    wavy[selection[0][0]].read_cube(selection[1][0], true, false, expert);
                }
                switch (_sel) {
                case 1: {
                    cls();
                    std::cout << "Integrated value of this cube is: ";
                    std::cout << scientific << setprecision(8) << setw(16) << wavy[selection[0][0]].get_cube_ptr(selection[1][0])->sum() << endl;
                    break;
                }
                case 2: {
                    cls();
                    std::cout << "New cube saved as: " << wavy[selection[0][0]].make_super_cube(selection[1][0]) << endl;
                    break;
                }
                case 3: {
                    double thresh = 0.0;
                    std::cout << "Please enter threshhold: ";
                    std::cin >> thresh;
                    wavy[selection[0][0]].apply_cube_thresh(selection[1][0], thresh);
                    break;
                }
                case 4: {
                    cls();
                    std::cout << "Integrated absolute_differnce value of this cube is: ";
                    std::cout << scientific << setprecision(8) << setw(16) << wavy[selection[0][0]].get_cube_ptr(selection[1][0])->diff_sum() << endl;
                    break;
                }
                case 0: {
                    end = true;
                    break;
                }
                default:
                    std::cout << "Sorry, i did not get that ..." << endl;
                    break;
                }
                break;
            }
            case 2: {
                if (nr_cubes > 1) {
                    while (!end) {
                        std::cout << "Which operation do you want to perform? ";
                        if (expert)std::cout << "(The 1. cube selected will be altered by the 2. one in cases of x=)";
                        std::cout << endl
                            << "1) +" << endl
                            << "2) -" << endl
                            << "3) *" << endl
                            << "4) /" << endl
                            << "5) RSR (Real space R-Value)" << endl
                            << "6) Mask (if(value!=0) value, else 0)" << endl
                            << "7) Negative mask (if(value==0) value, else 0" << endl
                            << "8) Mask with threshhold (if value(2) < thresh = 0, else = value(1))" << endl
                            << "9) Weighted Jaccord similarity and distance" << endl;
                        if (expert)std::cout << "11) +=" << endl
                            << "12) -=" << endl
                            << "13) *=" << endl
                            << "14) /=" << endl;
                        std::cout << "0) Go back to previous menu" << endl;
                        std::cin >> temp;
                        selection[0].resize(2);
                        selection[1].resize(2);
                        if (nr_cubes >= 2 && temp != 0) {
                            select_cubes(selection, wavy, 2);
                            if (opt.debug) {
                                std::cout << "Selection:" << endl;
                                for (int i = 0; i < 2; i++)
                                    std::cout << selection[0][i] << "." << selection[1][i] << endl;
                            }
                        }
                        cls();
                        if (temp >= 1 && temp < 5 && !wavy[selection[0][1]].get_cube_loaded(selection[1][1])) {
                            std::cout << "Loading full file now!" << endl;
                            wavy[selection[0][1]].read_cube(selection[1][1], true, false, expert);
                        }
                        switch (temp) {
                        case 1: {
                            wavy[selection[0][0]].push_back_cube(*wavy[selection[0][0]].get_cube_ptr(selection[1][0]) + *wavy[selection[0][1]].get_cube_ptr(selection[1][1]));
                            if (wavy[selection[0][0]].get_cube_ptr(wavy[selection[0][0]].get_cube_count() - 1)->get_size(0) != 0)
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                        }
                              break;
                        case 2: {
                            wavy[selection[0][0]].push_back_cube(*wavy[selection[0][0]].get_cube_ptr(selection[1][0]) - *wavy[selection[0][1]].get_cube_ptr(selection[1][1]));
                            if (wavy[selection[0][0]].get_cube_ptr(wavy[selection[0][0]].get_cube_count() - 1)->get_size(0) != 0)
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                        }
                              break;
                        case 3: {
                            wavy[selection[0][0]].push_back_cube(*wavy[selection[0][0]].get_cube_ptr(selection[1][0]) * *wavy[selection[0][1]].get_cube_ptr(selection[1][1]));
                            if (wavy[selection[0][0]].get_cube_ptr(wavy[selection[0][0]].get_cube_count() - 1)->get_size(0) != 0)
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                        }
                              break;
                        case 4: {
                            wavy[selection[0][0]].push_back_cube(*wavy[selection[0][0]].get_cube_ptr(selection[1][0]) / *wavy[selection[0][1]].get_cube_ptr(selection[1][1]));
                            if (wavy[selection[0][0]].get_cube_ptr(wavy[selection[0][0]].get_cube_count() - 1)->get_size(0) != 0)
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                        }
                              break;
                        case 5: {
                            for (int i = 0; i < 2; i++)
                                if (wavy[selection[0][i]].get_cube_loaded(selection[1][i]) == false) {
                                    if (opt.debug)std::cout << "Loading full file now!" << endl;
                                    wavy[selection[0][i]].read_cube(selection[1][i], true, false, expert);
                                }
                            double result = wavy[selection[0][0]].get_cube_ptr(selection[1][0])->rrs(*wavy[selection[0][1]].get_cube_ptr(selection[1][1]));
                            if (result != -1) {
                                std::cout << "Operation succesfull!" << endl;
                                std::cout << "RSR: " << result << endl;
                            }
                            else std::cout << "Sorry, something went wrong!" << endl;
                        }
                              break;
                        case 6: {
                            for (int i = 0; i < 2; i++)
                                if (wavy[selection[0][i]].get_cube_loaded(selection[1][i]) == false) {
                                    if (opt.debug)std::cout << "Loading full file now!" << endl;
                                    wavy[selection[0][i]].read_cube(selection[1][i], true, false, expert);
                                }
                            if (wavy[selection[0][0]].apply_cube_mask(selection[1][0], *wavy[selection[0][1]].get_cube_ptr(selection[1][1])))
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                        }
                              break;
                        case 7:
                            for (int i = 0; i < 2; i++)
                                if (wavy[selection[0][i]].get_cube_loaded(selection[1][i]) == false) {
                                    if (opt.debug)
                                        std::cout << "Loading full file now!" << endl;
                                    wavy[selection[0][i]].read_cube(selection[1][i], true, false, expert);
                                }
                            if (wavy[selection[0][0]].apply_cube_negative_mask(selection[1][0], *wavy[selection[0][1]].get_cube_ptr(selection[1][1])))
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                            break;
                        case 8: {
                            for (int i = 0; i < 2; i++)
                                if (wavy[selection[0][i]].get_cube_loaded(selection[1][i]) == false) {
                                    if (opt.debug)
                                        std::cout << "Loading full file now!" << endl;
                                    wavy[selection[0][i]].read_cube(selection[1][i], true, false, expert);
                                }
                            std::cout << "Please give threshhold to use: ";
                            double thresh = 0.0;
                            std::cin >> thresh;
                            if (wavy[selection[0][0]].apply_cube_thresh(selection[1][0], *wavy[selection[0][1]].get_cube_ptr(selection[1][1]), thresh))
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                            break;
                        }
                        case 9: {
                            for (int i = 0; i < 2; i++)
                                if (wavy[selection[0][i]].get_cube_loaded(selection[1][i]) == false) {
                                    if (opt.debug)std::cout << "Loading full file now!" << endl;
                                    wavy[selection[0][i]].read_cube(selection[1][i], true, false, expert);
                                }
                            double result = wavy[selection[0][0]].get_cube_ptr(selection[1][0])->jaccard(*wavy[selection[0][1]].get_cube_ptr(selection[1][1]));
                            if (result != -1) {
                                std::cout << "Operation succesfull!" << endl;
                                std::cout << "Jaccord similarity: " << result << " Jaccord distance: " << 1 - result << endl;
                            }
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                        }
                              break;
                        case 11:
                            for (int i = 0; i < 2; i++)
                                if (wavy[selection[0][i]].get_cube_loaded(selection[1][i]) == false) {
                                    if (opt.debug)std::cout << "Loading full file now!" << endl;
                                    wavy[selection[0][i]].read_cube(selection[1][i], true, false, expert);
                                }
                            if (wavy[selection[0][0]].cube_add(selection[1][0], *wavy[selection[0][1]].get_cube_ptr(selection[1][1])))
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                            break;
                        case 12:
                            for (int i = 0; i < 2; i++)
                                if (wavy[selection[0][i]].get_cube_loaded(selection[1][i]) == false) {
                                    if (opt.debug)std::cout << "Loading full file now!" << endl;
                                    wavy[selection[0][i]].read_cube(selection[1][i], true, false, expert);
                                }
                            if (wavy[selection[0][0]].cube_subtract(selection[1][0], *wavy[selection[0][1]].get_cube_ptr(selection[1][1])))
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                            break;
                        case 13:
                            for (int i = 0; i < 2; i++)
                                if (wavy[selection[0][i]].get_cube_loaded(selection[1][i]) == false) {
                                    if (opt.debug)std::cout << "Loading full file now!" << endl;
                                    wavy[selection[0][i]].read_cube(selection[1][i], true, false, expert);
                                }
                            if (wavy[selection[0][0]].cube_multiply(selection[1][0], *wavy[selection[0][1]].get_cube_ptr(selection[1][1])))
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                            break;
                        case 14:
                            for (int i = 0; i < 2; i++)
                                if (wavy[selection[0][i]].get_cube_loaded(selection[1][i]) == false) {
                                    if (opt.debug)std::cout << "Loading full file now!" << endl;
                                    wavy[selection[0][i]].read_cube(selection[1][i], true, false, expert);
                                }
                            if (wavy[selection[0][0]].cube_divide(selection[1][0], *wavy[selection[0][1]].get_cube_ptr(selection[1][1])))
                                std::cout << "Operation succesfull!" << endl;
                            else
                                std::cout << "Sorry, something went wrong!" << endl;
                            break;
                        case 0:
                            end = true;
                            break;
                        default:
                            std::cout << "Sorry, i did not get that..." << endl;
                            break;
                        }
                    }
                    end = false;
                }
                break;
            }
            case 3:
                std::cout << "Which analysis do you want to perform?" << endl
                    << "1) Separate into Basins according to critical point search" << endl;
                //and try new BCP implementation"
                if (expert)std::cout << "2) Separate into Basins according to critical point search and try new BCP implementation" << endl;
                std::cout << "0) Go back to selection" << endl;
                std::cin >> temp;
                selection[0].resize(1);
                selection[1].resize(1);
                if (nr_cubes >= 1 && temp != 0) select_cubes(selection, wavy, 1, false, opt.debug);
                if (opt.debug) std::cout << "selection: " << selection[0][0] << " " << selection[1][0] << endl;
                switch (temp) {
                case 1:
                    if (wavy[selection[0][0]].get_cube_loaded(selection[1][0]) == false) {
                        std::cout << "Loading full file now!" << endl;
                        if (wavy[selection[0][0]].read_cube(selection[1][0], true, false, expert) == false)
                        {
                            std::cout << "ERROR reading full file! Aborting" << endl;
                            break;
                        }
                    }
                    if (b2c(wavy[selection[0][0]].get_cube_ptr(selection[1][0]), wavy[selection[0][0]].get_atoms(), opt.debug, false) == false)
                        std::cout << "something went wrong!" << endl;
                    break;
                case 2:
                    if (expert) {
                        if (wavy[selection[0][0]].get_cube_loaded(selection[1][0]) == false) {
                            std::cout << "loading full file now!" << endl;
                            if (wavy[selection[0][0]].read_cube(selection[1][0], true, false, expert) == false) {
                                std::cout << "ERROR reading full file! Aborting!" << endl;
                                break;
                            }
                        }
                        if (b2c(wavy[selection[0][0]].get_cube_ptr(selection[1][0]), wavy[selection[0][0]].get_atoms(), opt.debug, true) == false)
                            std::cout << "something went wrong" << endl;
                        break;
                    }
                    break;
                case 0:
                    break;
                default:
                    std::cout << "sorry, i did not understand that, please try again!" << endl;
                    break;
                }
            }
            break;
        }
        case 'r':
        case 'R': {
            if (expert) {
                std::cout << "What kind of file do you want to load?" << endl
                    << "W) Wavefunction in known format (gbw/wfn/wfx/fchk/molden)" << endl
                    << "G) Grid or Cube file" << endl
                    << "E) Exit to main menu" << endl;
                std::cin >> sel;
                cls();
                switch (sel) {
                case 'G':
                case 'g': {
                    std::filesystem::path path;
                    while (!exists(path)) {
                        std::cout << "Please give Path to the cube you want to read: ";
                        std::cin >> path;
                        if (exists(path)) break;
                        else std::cout << "Sorry, couldn't find the file! Try again!" << endl;
                    }
                    if (wavy.size() == 0) {
                        activewave = 0;
                        wavy.emplace_back(path);
                    }
                    wavy[activewave].push_back_cube(cube(path, false, wavy[activewave], std::cout, expert));
                    cls();
                    break;
                }
                case 'W':
                case 'w': {
                    while (!end) {
                        string path;
                        bool new_wfn = false;
                        std::cout << "Path to the wavefunction file: ";
                        std::cin >> path;
                        int tries = 0;
                        if (path.find(".wfn") == -1 && path.find(".ffn") == -1) {
                            do {
                                std::cout << "This doesn't look like a .wfn or .ffn file, try again: " << endl;
                                tries++;
                                std::cin >> path;
                            } while (path.find(".wfn") == -1 && path.find(".ffn") == -1 && tries < 3);
                            if (tries == 3) {
                                std::cout << "Sorry, make sure you know the filename and try again!" << endl;
                                break;
                            }
                        }
                        if (activewave < wavy.size() && wavy.size()>0) {
                            std::cout << "This will delete the previously loaded wavefunction " << wavy[activewave].get_path() << "! Are you sure?";
                            if (!yesno()) {
                                std::cout << "Do you want to load it as a new wavefunction? " << endl;
                                if (yesno()) new_wfn = true;
                                end = true;
                                break;
                            }
                            wavy.erase(wavy.begin() + activewave);
                            if (path.find(".wfn") != -1) wavy.insert(wavy.begin() + activewave - 1, WFN(e_origin::wfn));
                            else wavy.insert(wavy.begin() + activewave - 1, WFN(e_origin::ffn));
                        }
                        else if (wavy.size() == 0 || activewave == wavy.size() || new_wfn) {
                            if (opt.debug) std::cout << "Making a new wavefunction!" << endl;
                            if (path.find(".wfn") != -1) wavy.emplace_back(e_origin::wfn);
                            else wavy.emplace_back(e_origin::ffn);
                            if (wavy.size() > 1) activewave++;
                            if (opt.debug) std::cout << "Size: " << wavy.size() << " active: " << activewave << endl;
                        }
                        wavy[activewave].read_known_wavefunction_format(path, std::cout, opt.debug);
                        cls();
                        if (opt.debug) std::cout << wavy[activewave].get_ncen() << endl;
                        break;
                    }
                    cls();
                    end = false;
                    break;
                }
                case 'e':
                case 'E':
                    break;
                default:
                    std::cout << "Sorry, did not understand that, going back to main menu!" << endl;
                    break;
                }
            }
            else {
                std::filesystem::path filename;
                vector <string> temp;
                temp.resize(3);
                temp[0] = "Wavefunction files (wfn,ffn,fchk,wfx) | *.wfn *.ffn *.Fchk *.fchk *.FChk *.wfx *gbw *.molden *.molden.input *.xtb";
                temp[1] = "Cube files (cub, cube, grd) | *.cube *.cub *.grd";
                temp[2] = "All filetypes | *";
                if (!open_file_dialog(filename, opt.debug, temp, opt.cwd.string())) {
                    std::cout << "Error encountered!" << endl;
                    break;
                }
                else {
                    if (wavy.size() > activewave && wavy.size() > 0) {
                        std::cout << "This will delete the previously loaded wavefunction " << wavy[activewave].get_path() << "! Are you sure?";
                        if (!yesno()) {
                            std::cout << "Do you want to load it as a new wavefunction? " << endl;
                            if (yesno()) activewave++;
                            else continue;
                        }
                        wavy.erase(wavy.begin() + activewave);
                    }
                    wavy.emplace(wavy.begin() + activewave, filename, opt.debug);
                }
            }
            break;
        }
        case 'M':
        case 'm': {
            if (wavy.size() < 1) {
                std::cout << "First you need to read a wavefunction!" << endl;
                break;
            }
            int msel = 0;
            string label = "?";
            float x, y, z = 0.0;
            int charge = 0;
            std::cout << "What do you want to do?" << endl
                << "D) Delete a center or set of values for a certain function" << endl
                << "A) Add a center" << endl
                << "C) Change a value" << endl
                << "E) Exit to Main Menu" << endl;
            std::cin >> sel;
            cls();
            switch (sel) {
            case 'D':
            case 'd':
                std::cout << "A) Atom" << endl
                    << "E) Exponent" << endl
                    << "B) Back to main menu" << endl;
                std::cin >> sel;
                end = false;
                switch (sel) {
                case 'A':
                case 'a':
                    std::cout << "The list of centers:\n";
                    wavy[activewave].list_centers();
                    std::cout << "Which center do you want to delete?\n";
                    std::cin >> msel;
                    if (msel >= 0) {
                        if (wavy[activewave].remove_center(msel)) {
                            std::cout << "Deleted center nr " << msel << " succesfully! Going back to Main Menu.\n";
                            wavy[activewave].set_modified();
                        }
                        else std::cout << "Something went wrong... Sorry, let's start again...\n";
                    }
                    else std::cout << "Wrong selection, start again!\n";
                    break;
                case 'E':
                case 'e':
                    while (!end) {
                        wavy[activewave].list_primitives();
                        std::cout << "Which exponent out of " << wavy[activewave].get_nex() << " do you want to delete?\n";
                        std::cin >> msel;
                        if (msel<0 || msel>wavy[activewave].get_nex()) {
                            std::cout << "Sorry, wrong input";
                            continue;
                        }
                        std::cout << "This is the set of information you requested:\n";
                        wavy[activewave].print_primitive(msel);
                        end = true;
                    }
                    break;
                case 'B':
                case 'b':
                    cls();
                    break;
                }
                break;
            case 'A':
            case 'a':
                std::cout << "A) Atom\nE) Exponent\nB) Back to main menu" << endl;
                std::cin >> sel;
                end = false;
                switch (sel) {
                case 'E':
                case 'e':
                    while (!end) {
                        std::cout << "Please remember that this exponent will be apended to the data structure, "
                            << "it will not be sorted in any way!\nCentre Assignement: ";
                        int temp_cen = 0;
                        std::cin >> temp_cen;
                        if (temp_cen > wavy[activewave].get_ncen() || temp_cen < 0) {
                            std::cout << "Wrong input, start over!";
                            continue;
                        }
                        std::cout << "Type: ";
                        int temp_type = 0;
                        std::cin >> temp_type;
                        std::cout << "Exponent: ";
                        double temp_exp = 0.0;
                        std::cin >> temp_exp;
                        vector<double> temp_val;
                        temp_val.resize(wavy[activewave].get_nmo());
                        for (int i = 0; i < wavy[activewave].get_nmo(); i++) {
                            std::cout << "Enter coefficient for MO " << i << ":";
                            std::cin >> temp_val[i];
                            if (temp_val[i] < -1000 || temp_val[i]>1000) {
                                std::cout << "Wrong input, please try again...\n";
                                i--;
                                continue;
                            }
                        }
                        std::cout << "Let me recapitulate: Center " << temp_cen << " type: " << temp_type << " exp: " << temp_exp
                            << " and the MO coefficients:\n";
                        for (int i = 0; i < wavy[activewave].get_nmo(); i++) {
                            std::cout << temp_val[i] << "   ";
                            if (i % 5 == 0)std::cout << endl;
                        }
                        std::cout << "is this correct?";
                        if (yesno()) end = true;
                    }
                    break;
                case 'A':
                case 'a':
                    while (!end) {
                        std::cout << "The list of centers:\n";
                        wavy[activewave].list_centers();
                        std::cout << "Which center do you want to add?\nlabel:";
                        std::cin >> label;
                        std::cout << "x: ";
                        std::cin >> x;
                        if (x < -99.999 || x > 99.999) {
                            std::cout << "Sorry, number too large\n";
                            continue;
                        }
                        std::cout << "y: ";
                        std::cin >> y;
                        if (y < -99.999 || y > 99.999) {
                            std::cout << "Sorry, number too large\n";
                            continue;
                        }
                        std::cout << "z: ";
                        std::cin >> z;
                        if (z < -99.999 || z > 99.999) {
                            std::cout << "Sorry, number too large\n";
                            continue;
                        }
                        std::cout << "charge: ";
                        std::cin >> charge;
                        if (charge <= 0 || charge > 118) {
                            std::cout << "Sorry, that atom is not yet disovered\n";
                            continue;
                        }
                        wavy[activewave].push_back_atom(label, x, y, z, charge);
                    }
                    break;
                case 'B':
                case 'b':
                    cls();
                    break;
                }
                break;
            case 'C':
            case 'c':
                if (wavy.size() < 1) continue;
                std::cout << "The center/type/exponent status until now is:\n";
                wavy[activewave].list_primitives();
                std::cout << "What do you want to change?\nC) Center Assignement\nT) Type assignement\nE) Exponent\nM) MO coefficient\nQ) Quit\n";
                std::cin >> sel;
                switch (sel) {
                case 'Q':
                case 'q':
                    cls();
                    break;
                case 'C':
                case 'c':
                    std::cout << "Which one do you want to change? (0=return to menu)";
                    std::cin >> msel;
                    if (msel > wavy[activewave].get_nex() || msel < 0) {
                        std::cout << "Wrong input, start again!\n";
                        break;
                    }
                    else if (msel == 0) break;
                    wavy[activewave].change_center(msel);
                    break;
                case 'E':
                case 'e':
                    std::cout << "What is the nr. of the exponent you want to change? (0=return to menu)" << endl;
                    std::cin >> msel;
                    if (msel > wavy[activewave].get_nex() || msel < 0) {
                        std::cout << "Wrong input, start again!\n";
                        break;
                    }
                    else if (msel == 0) break;
                    wavy[activewave].change_exponent(msel);
                    break;
                case 'T':
                case 't':
                    std::cout << "What is the nr. of the type you want to change? (0=return to menu)" << endl;
                    std::cin >> msel;
                    if (msel > wavy[activewave].get_nex() || msel < 0) {
                        std::cout << "Wrong input, start again!\n";
                        break;
                    }
                    else if (msel == 0) break;
                    wavy[activewave].change_type(msel);
                    break;
                case 'M':
                case 'm':
                    bool _end = false;
                    while (!_end) {
                        std::cout << "Which MO out of " << wavy[activewave].get_nmo() << " MOs? (0=return to menu)\n";
                        int MOsel;
                        std::cin >> MOsel;
                        if (MOsel > wavy[activewave].get_nmo() || MOsel < 0) {
                            std::cout << "This is not a valid choice...\n";
                            continue;
                        }
                        else if (MOsel == 0) break;
                        else {
                            int coef_sel = 0;
                            std::cout << "Which coefficient do you want to change? (0=return to menu)\n";
                            std::cin >> coef_sel;
                            if (coef_sel > wavy[activewave].get_nex() || coef_sel < 0) {
                                std::cout << "This is not a valid choice...\n";
                                continue;
                            }
                            else if (coef_sel == 0) break;
                            else {
                                std::cout << "Please enter the new value for the coefficient: ";
                                double new_coef;
                                std::cin >> new_coef;
                                wavy[activewave].set_MO_coef(MOsel, coef_sel, new_coef);
                            }
                            Enter();
                            _end = true;
                        }
                    }
                    break;
                }
            case 'E':
            case 'e':
                break;
            }
            break;
        }
        case 'n':
        case 'N': {
            std::cout << "Bondwise analysis (B) or NCIplot features (N)? ";
            char seln;
            std::cin >> seln;
            switch (seln) {
            case 'B':
            case 'b': {
                string inputfile;
                std::cout << "Which file to use for definition of bonds?" << endl;
                std::cin >> inputfile;
                if (autobonds(opt.debug, wavy[activewave], inputfile, false) != 1) std::cout << "Sorry, looks like something went wrong..." << endl;
                break;
            }
            case 'N':
            case 'n': {
                std::cout << "Starting acuNCI..." << endl;
                acu_nci(wavy, opt);
                break;
            }
            default:
                std::cout << "sorry, didn't understand that" << endl;
                break;
            }
            //cls();
            break;
        }
        case 'S':
        case 's': {
            if (wavy.size() < 1) {
                cls();
                continue;
            }
            vector <string> endings;
            bool w_or_c = true;
            bool convert = false;
            if (wavy[activewave].get_cube_count() > 0 && wavy[activewave].get_origin() != 3) {
                while (true) {
                    std::cout << "Do you want to save the wavefunction (W) or associated cubes (C) or convert cubes into non-cube format (N)? ";
                    char input;
                    std::cin >> input;
                    switch (input) {
                    case 'C':
                    case 'c': {
                        w_or_c = false;
                        break;
                    }
                    case 'W':
                    case 'w': {
                        w_or_c = true;
                        break;
                    }
                    case 'n':
                    case 'N': {
                        convert = true;
                        break;
                    }
                    default:
                        std::cout << "Sorry, i did not get that! Try again!" << endl;
                        continue;
                    }
                    break;
                }
            }
            else if (wavy[activewave].get_origin() == 3) {
                while (true) {
                    std::cout << "Do you want to save the cube in .cube (C) format or convert cubes into non-cube format (N)? ";
                    char input;
                    std::cin >> input;
                    switch (input) {
                    case 'C':
                    case 'c': {
                        convert = false;
                        break;
                    }
                    case 'n':
                    case 'N': {
                        convert = true;
                        break;
                    }
                    default:
                        std::cout << "Sorry, i did not get that! Try again!" << endl;
                        continue;
                    }
                    break;
                }
            }

            if (convert) {
                while (true) {
                    std::cout << "Which format do you want to convert to?\nD) DGrid\nX) XD-Graph\n";
                    char input;
                    std::cin >> input;
                    switch (input) {
                    case 'D':
                    case 'd': {
                        if (wavy[activewave].get_cube_count() >= 1) {
                            std::cout << "Which cube do you want to save?" << endl;
                            selection.resize(2);
                            selection[0].resize(1);
                            selection[1].resize(1);
                            if (!expert) select_cubes(selection, wavy, 1);
                            else {
                                int nr1 = 0;
                                for (int i = 0; i < wavy[selection[0][0]].get_cube_count(); i++) {
                                    std::cout << " " << i << ") " << wavy[activewave].get_cube_path(i);
                                    if (!exists(wavy[activewave].get_cube_path(i)))std::cout << " (MEM ONLY)";
                                    std::cout << endl;
                                }
                                std::cout << "Cube 1: "; cin >> nr1;
                                while (nr1 >= wavy[activewave].get_cube_count() || nr1 < 0) {
                                    std::cout << "Invalid choice, select again: ";
                                    std::cin >> nr1;
                                }
                                selection[0][0] = activewave;
                                selection[1][0] = nr1;
                            }
                        }
                        std::filesystem::path path = wavy[selection[0][0]].get_cube_path(selection[1][0]);
                        if (!wavy[selection[0][0]].get_cube_loaded(selection[1][0]))
                            wavy[selection[0][0]].read_cube(selection[1][0], true, false, false);
                        wavy[selection[0][0]].write_cube_dgrid(selection[1][0], path.replace_extension(".dgrid"), opt.debug);
                        break;
                    }
                    case 'X':
                    case 'x': {
                        if (wavy[activewave].get_cube_count() >= 1) {
                            std::cout << "Which cube do you want to save?" << endl;
                            selection.resize(2);
                            selection[0].resize(1);
                            selection[1].resize(1);
                            if (!expert) select_cubes(selection, wavy, 1);
                            else {
                                int nr1 = 0;
                                for (int i = 0; i < wavy[selection[0][0]].get_cube_count(); i++) {
                                    std::cout << " " << i << ") " << wavy[activewave].get_cube_path(i);
                                    if (!exists(wavy[activewave].get_cube_path(i)))std::cout << " (MEM ONLY)";
                                    std::cout << endl;
                                }
                                std::cout << "Cube 1: "; cin >> nr1;
                                while (nr1 >= wavy[activewave].get_cube_count() || nr1 < 0) {
                                    std::cout << "Invalid choice, select again: ";
                                    std::cin >> nr1;
                                }
                                selection[0][0] = activewave;
                                selection[1][0] = nr1;
                            }
                        }
                        std::filesystem::path path = wavy[selection[0][0]].get_cube_path(selection[1][0]);
                        if (!wavy[selection[0][0]].get_cube_loaded(selection[1][0]))
                            wavy[selection[0][0]].read_cube(selection[1][0], true, false, false);
                        wavy[selection[0][0]].write_cube_xdgraph(selection[1][0], path.replace_extension(".xdgrid"), opt.debug);
                        break;
                    }
                    default:
                        std::cout << "Sorry, i did not get that! Try again!" << endl;
                        continue;
                    }
                    break;
                }
            }
            else {
                if ((wavy[activewave].get_origin() == 2 || wavy[activewave].get_origin() == 4) && w_or_c) {
                    std::cout << "Which format do you want to save the wavefunction in?" << endl
                        << "W) WFN format" << endl
                        << "F) Fchk format" << endl;
                    std::cin >> sel;
                    switch (sel) {
                    case 'W':
                    case 'w': {
                        endings.push_back(".wfn");
                        endings.push_back(".ffn");
                        std::filesystem::path path;
                        if (!expert) save_file_dialog(path, opt.debug, endings, opt.cwd.string());
                        else {
                            std::cout << "Enter filename: ";
                            std::cin >> path;
                            while (exists(path)) {
                                std::cout << path << " exists, do you want to overwrite it? ";
                                if (!yesno()) {
                                    std::cout << "Then try again: ";
                                    std::cin >> path;
                                }
                            }
                        }
                        bool all = false;
                        if (path.extension() == ".wfn") {
                            std::cout << "Do you want to write all MOs?" << endl;
                            all = yesno();
                        }
                        if (!wavy[activewave].write_wfn(path, opt.debug, !all)) {
                            Enter();
                            cls();
                        }
                        else {
                            if (opt.debug) Enter();
                            cls();
                            std::cout << "Wrote Wavefunction!\n";
                        }
                        break;
                    }
                    case 'F':
                    case 'f': {
                        endings.push_back(".fchk");
                        endings.push_back(".Fchk");
                        endings.push_back(".FChk");
                        std::filesystem::path outputname = wavy[activewave].get_path();
                        if (opt.debug)std::cout << "Loaded path..." << endl;
                        if (!expert) save_file_dialog(outputname, opt.debug, endings, opt.cwd.string());
                        else {
                            std::cout << "Enter filename: ";
                            std::cin >> outputname;
                            while (exists(outputname)) {
                                std::cout << outputname << " exists, do you want to overwrite it? ";
                                if (!yesno()) {
                                    std::cout << "Then try again: ";
                                    std::cin >> outputname;
                                }
                            }
                        }
                        outputname.replace_extension(".fchk");
                        string basis_temp = wavy[activewave].get_basis_set_name();
                        if (basis_temp.length() < 3) {
                            int tries = 0;
                            bool _end = false;
                            string temp;
                            while (!_end && tries != 3) {
                                std::cout << "Please give the name of the basis set you want to use: ";
                                std::cin >> temp;
                                std::filesystem::path basis_set_file(opt.basis_set_path);
                                basis_set_file.append(temp);
                                if (opt.debug)std::cout << "looking for: " << basis_set_file << endl;
                                if (exists(basis_set_file)) _end = true;
                                else tries++;
                            }
                            if (tries == 3) {
                                std::cout << "Sorry, this takes too long... please make sure you know what you want and try again!" << endl;
                                Enter();
                                break;
                            }
                            wavy[activewave].change_basis_set_name(temp);
                        }
                        if (expert) {
                            std::cout << "What is the charge of your molecule?" << endl;
                            int temp = 0;
                            std::cin >> temp;
                            wavy[activewave].assign_charge(temp);
                        }
                        else wavy[activewave].assign_charge(wavy[activewave].calculate_charge());
                        if (wavy[activewave].get_multi() == 0) wavy[activewave].guess_multiplicity(std::cout);
                        free_fchk(std::cout, outputname, opt.basis_set_path, wavy[activewave], opt.debug);
                        break;
                    }
                    default: {
                        std::cout << "Sorry, i didn't get that, could you try it again?\n";
                        Enter();
                        cls();
                    }
                    }
                }
                else if (wavy[activewave].get_origin() == 3 || !w_or_c) {
                    std::filesystem::path path;
                    if (wavy[activewave].get_cube_count() <= 0) {
                        std::cout << "No cubes loaded!" << endl;
                        break;
                    }
                    int nr1 = 0;
                    if (wavy[activewave].get_cube_count() >= 1) {
                        std::cout << "Which cube do you want to save?" << endl;
                        selection.resize(2);
                        selection[0].resize(1);
                        selection[1].resize(1);
                        if (!expert) select_cubes(selection, wavy, 1);
                        else {
                            for (int i = 0; i < wavy[selection[0][0]].get_cube_count(); i++) {
                                std::cout << " " << i << ") " << wavy[activewave].get_cube_path(i);
                                if (!exists(wavy[activewave].get_cube_path(i)))std::cout << " (MEM ONLY)";
                                std::cout << endl;
                            }
                            std::cout << "Cube 1: "; cin >> nr1;
                            while (nr1 >= wavy[activewave].get_cube_count() || nr1 < 0) {
                                std::cout << "Invalid choice, select again: ";
                                std::cin >> nr1;
                            }
                            selection[0][0] = activewave;
                            selection[1][0] = nr1;
                        }
                    }

                    endings.push_back(".cube");
                    endings.push_back(".cub");
                    if (!expert)
                        save_file_dialog(path, opt.debug, endings, opt.cwd.string());
                    else {
                        std::cout << "Give filepath please: ";
                        std::cin >> path;
                    }
                    wavy[selection[0][0]].write_cube_file(selection[1][0], path, opt.debug);
                }
            }
            endings.resize(0);
            break;
        }
        case 'O':
        case 'o': {
            if (wavy.size() < 1) {
                std::cout << "First you need to read a wavefunction!" << endl;
                break;
            }
            if (wavy[activewave].get_origin() == 2 || wavy[activewave].get_origin() == 4) {
                std::cout << "Sorting wavefunction!" << endl;
                wavy[activewave].sort_wfn(wavy[activewave].check_order(opt.debug), opt.debug);
                if (opt.debug) Enter();
                cls();
            }
            else {
                std::cout << "I can only sort .wfn/ffn files!" << endl;
                Enter();
                std::cout << endl;
            }
            break;
        }
        case 'B':
        case 'b': {
            if (wavy.size() < 1) {
                std::cout << "First you need to read a wavefunction!" << endl;
                break;
            }
            if (!read_basis_set_vanilla(opt.basis_set_path, wavy[activewave], opt.debug))std::cout << "Problem during reading of the basis set!" << endl;
            if (opt.debug) Enter();
            cls();
            break;
        case 'a':
        case 'A':
            if (wavy.size() < 1) {
                cls();
                break;
            }
            vector < vector < unsigned int > > _selection;
            _selection.resize(1);
            select_cubes(_selection, wavy, 1, true);
            activewave = _selection[0][0];
            cls();
            break;
        }
        case 'e':
        case 'E':
            if (expert) expert = false;
            else expert = true;
            break;
        case 'q':
        case 'Q':
            if (!unsaved_files(wavy))
                end = true;
            else {
                std::cout << "There are unsaved files! Are you sure you want to exit? ";
                if (yesno())
                    end = true;
            }
            break;
        case 'l':
        case 'L':
            cout << "Please give the number of cpus to use (-1 = all): ";
            int threads;
            cin >> threads;
            if (threads < -1 || threads == 0) {
                cout << "Sorry, that is not a valid number of threads!" << endl;
                break;
            }
            opt.threads = threads;
            if (opt.threads != -1)
            {
#ifdef _OPENMP
                omp_set_num_threads(opt.threads);
#endif
            }
            std::cout << "Number of threads set to " << opt.threads << endl;
            break;
        case 'd':
        case 'D':
            opt.debug = true;
            cls();
            break;
        case 'u':
        case 'U':
            if (check_bohr(wavy[activewave], opt.debug)) {
                cls();
                std::cout << "Appears to be in bohr!" << endl;
            }
            else {
                cls();
                std::cout << "Appears to be in Angström!" << endl;
            }
            break;
        default:
            if (opt.debug) {
                std::cout << "This command is unknown!" << endl;
                Enter();
            }
            cls();
            std::cout << "Sorry, i didn't get that, could you try it again?\n";
            break;
        }
    }
    cls();
    return 0;
}
