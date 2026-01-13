#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "bondwise_analysis.h"
#include "properties.h"
#include "libCintMain.h"
#include "nos_math.h"
#include "integration_params.h"

#ifdef NSA2DEBUG
void print_dmatrix2(const dMatrix2& EVC2, const std::string name) {
    std::cout << std::endl << name << ":\n";
    for (int i = 0; i < EVC2.extent(0); i++) {
        for (int j = 0; j < EVC2.extent(1); j++)
            std::cout << std::setw(14) << std::setprecision(8) << std::fixed << EVC2(i, j) << " ";
        std::cout << std::endl;
    }
}
#endif


int compute_dens(WFN& wavy, bool debug, int* np, double* origin, double* gvector, double* incr, std::string& outname, bool rho, bool rdg, bool eli, bool lap) {
    options opt;

    cube CubeRho({ np[0], np[1], np[2] }, wavy.get_ncen(), rho),
        CubeRDG({ np[0], np[1], np[2] }, wavy.get_ncen(), rdg),
        CubeEli({ np[0], np[1], np[2] }, wavy.get_ncen(), eli),
        CubeLap({ np[0], np[1], np[2] }, wavy.get_ncen(), lap);

    std::string Oname = outname;
    std::vector<int> ntyp;
    for (int i = 0; i < 3; i++) {
        opt.properties.NbSteps[i] = np[i];
        if (debug) {
            std::cout << "gvector before: ";
            for (int j = 0; j < 3; j++) std::cout << gvector[j + 3 * i] << " ";
            std::cout << std::endl;
        }
        for (int j = 0; j < 3; j++)
            gvector[j + 3 * i] = incr[i] * gvector[j + 3 * i];
        if (debug) {
            for (int j = 0; j < 3; j++) std::cout << gvector[j + 3 * i] << " ";
            std::cout << std::endl;
        }
    }
    for (int i = 0; i < 3; i++) {
        CubeRho.set_origin(i, origin[i]);
        CubeRDG.set_origin(i, origin[i]);
        CubeEli.set_origin(i, origin[i]);
        CubeLap.set_origin(i, origin[i]);
        for (int j = 0; j < 3; j++) {
            CubeRho.set_vector(i, j, gvector[j + 3 * i]);
            CubeRDG.set_vector(i, j, gvector[j + 3 * i]);
            CubeEli.set_vector(i, j, gvector[j + 3 * i]);
            CubeLap.set_vector(i, j, gvector[j + 3 * i]);
        }
    }
    CubeRho.set_path(Oname + "_rho.cube");
    CubeRDG.set_path(Oname + "_rdg.cube");
    CubeEli.set_path(Oname + "_eli.cube");
    CubeLap.set_path(Oname + "_lap.cube");

    //opt.NbAtoms[0]=wavy.get_ncen();
    std::cout << "\n   .      ___________________________________________________________      .\n";
    std::cout << "  *.                                                                      *.\n";
    std::cout << "  *.  Wavefunction              : " << std::setw(20) << wavy.get_path().filename() << " / " << std::setw(5) << wavy.get_ncen() << " atoms      *.\n";
    std::cout << "  *.  OutPut filename Prefix    : " << std::setw(40) << Oname << "*.\n";
    std::cout << "  *.                                                                      *.\n";
    std::cout << "  *.  gridBox Min               : " << std::setw(11) << std::setprecision(6) << origin[0] << std::setw(12) << origin[1] << origin[2] << "     *.\n";
    std::cout << "  *.  gridBox Max               : " << std::setw(11) << std::setprecision(6) << (origin[0] + gvector[0] * np[0] + gvector[3] * np[1] + gvector[6] * np[2]) << std::setw(12) << (origin[1] + gvector[1] * np[0] + gvector[4] * np[1] + gvector[7] * np[2]) << (origin[2] + gvector[2] * np[0] + gvector[5] * np[1] + gvector[8] * np[2]) << "     *.\n";
    std::cout << "  *.  Increments(bohr)          : " << std::setw(11) << std::setprecision(6) << sqrt(pow(gvector[0], 2) + pow(gvector[1], 2) + pow(gvector[2], 2)) << std::setw(12) << sqrt(pow(gvector[3], 2) + pow(gvector[4], 2) + pow(gvector[5], 2)) << std::setw(12) << sqrt(pow(gvector[6], 2) + pow(gvector[7], 2) + pow(gvector[8], 2)) << "     *.\n";
    std::cout << "  *.  NbSteps                   : " << std::setw(11) << opt.properties.NbSteps[0] << std::setw(12) << opt.properties.NbSteps[1] << std::setw(12) << opt.properties.NbSteps[2] << "     *.\n";
    std::cout << "  *.                                                                      *.\n";
    std::cout << "  *.  Number of primitives      :     " << std::setw(5) << wavy.get_nex() << "                               *.\n";
    std::cout << "  *.  Number of MOs             :       " << std::setw(3) << wavy.get_nmo() << "                               *.\n";

    cube dummy({ 0, 0, 0 });
    Calc_Prop(
        CubeRho,
        CubeRDG,
        dummy,
        CubeEli,
        CubeLap,
        dummy,
        wavy,
        20.0,
        std::cout,
        false,
        false
    );

    std::cout << "  *.                                                                      *.\n";
    std::cout << "  *.  Writing .cube files ...                                             *.\n";
    if (rho && !rdg) {
        CubeRho.set_path(Oname + "_rho.cube");
        CubeRho.write_file(true, true);
    }
    if (rdg) {
        CubeRho.set_path(Oname + "_signed_rho.cube");
        CubeRho.write_file(true);
        CubeRho.set_path(Oname + "_rho.cube");
        CubeRho.write_file(true, true);
        CubeRDG.write_file(true);
    }
    if (eli) CubeEli.write_file(true);
    if (lap) CubeLap.write_file(true);

    std::cout << "  *                                                                       *\n";
    std::cout << "          ___________________________________________________________\n";
    return 0;
};

bond do_bonds(WFN& wavy,
    int mode_sel, bool mode_leng, bool mode_res,
    double res[], bool cub, double boxsize[],
    int atom1, int atom2, int atom3,
    const bool& debug, const bool& bohr,
    int runnumber,
    bool rho, bool rdg, bool eli, bool lap) {
    bond results{ "","","","",false,false,false };
    std::string line("");
    std::vector <std::string> label;
    label.resize(3);
    double coords1[3], coords2[3], coords3[3];
    int na = 0;
    double z[3], x[3], y[3], help[3], size[3], s2[3];
    int np[3];
    double ang2bohr;
    if (bohr)ang2bohr = 0.52917720859;
    else ang2bohr = 1.0;
    na = wavy.get_ncen();
    if (atom1 <= 0 || atom2 <= 0 || atom3 <= 0 || mode_sel <= 0 || mode_sel > 4 || atom1 > na || atom2 > na || atom3 > na || atom1 == atom2 || atom2 == atom3 || atom1 == atom3)
    {
        std::cout << "Invalid selections of atoms or mode_sel, please try again!\n";
        return(results);
    }
    if (debug) std::cout << "No. of Atoms selected: " << atom1 << atom2 << atom2 << "\n";
    for (int i = 0; i < 3; i++)
        coords1[i] = wavy.get_atom_coordinate(atom1 - 1, i);
    for (int i = 0; i < 3; i++)
        coords2[i] = wavy.get_atom_coordinate(atom2 - 1, i);
    for (int i = 0; i < 3; i++)
        coords3[i] = wavy.get_atom_coordinate(atom3 - 1, i);
    label[0] = wavy.get_atom_label(atom1 - 1);
    label[1] = wavy.get_atom_label(atom2 - 1);
    label[2] = wavy.get_atom_label(atom3 - 1);
    if (debug)
    {
        std::cout << "The Atoms found corresponding to your selection are:\n";
        std::cout << label[0] << " " << coords1[0] << " " << coords1[1] << " " << coords1[2] << "\n";
        std::cout << label[1] << " " << coords2[0] << " " << coords2[1] << " " << coords2[2] << "\n";
        std::cout << label[2] << " " << coords3[0] << " " << coords3[1] << " " << coords3[2] << "\n";
    }
    double znorm = sqrt((coords2[0] - coords1[0]) * (coords2[0] - coords1[0]) + (coords2[1] - coords1[1]) * (coords2[1] - coords1[1]) + (coords2[2] - coords1[2]) * (coords2[2] - coords1[2]));
    z[0] = (coords2[0] - coords1[0]) / znorm;
    z[1] = (coords2[1] - coords1[1]) / znorm;
    z[2] = (coords2[2] - coords1[2]) / znorm;
    help[0] = coords3[0] - coords2[0];
    help[1] = coords3[1] - coords2[1];
    help[2] = coords3[2] - coords2[2];
    y[0] = z[1] * help[2] - z[2] * help[1];
    y[1] = z[2] * help[0] - z[0] * help[2];
    y[2] = z[0] * help[1] - z[1] * help[0];
    double ynorm = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    y[0] = y[0] / ynorm;
    y[1] = y[1] / ynorm;
    y[2] = y[2] / ynorm;
    x[0] = (z[1] * y[2] - z[2] * y[1]);
    x[1] = (z[2] * y[0] - z[0] * y[2]);
    x[2] = (z[0] * y[1] - z[1] * y[0]);
    z[0] = -z[0];
    z[1] = -z[1];
    z[2] = -z[2];

    double hnorm = 0.0;
    double ringhelp[3] = {};
    if (mode_sel == 2 || mode_sel == 1) hnorm = znorm;
    else if (mode_sel == 3) hnorm = sqrt((coords3[0] - coords1[0]) * (coords3[0] - coords1[0]) + (coords3[1] - coords1[1]) * (coords3[1] - coords1[1]) + (coords3[2] - coords1[2]) * (coords3[2] - coords1[2]));
    else if (mode_sel == 4)
    {
        for (int r = 0; r < 3; r++) ringhelp[r] = (coords3[r] + coords1[r] + coords2[r]) / 3;
        hnorm = (sqrt((coords3[0] - ringhelp[0]) * (coords3[0] - ringhelp[0]) + (coords3[1] - ringhelp[1]) * (coords3[1] - ringhelp[1]) + (coords3[2] - ringhelp[2]) * (coords3[2] - ringhelp[2]))
            + sqrt((coords1[0] - ringhelp[0]) * (coords1[0] - ringhelp[0]) + (coords1[1] - ringhelp[1]) * (coords1[1] - ringhelp[1]) + (coords1[2] - ringhelp[2]) * (coords1[2] - ringhelp[2]))
            + sqrt((coords2[0] - ringhelp[0]) * (coords2[0] - ringhelp[0]) + (coords2[1] - ringhelp[1]) * (coords2[1] - ringhelp[1]) + (coords2[2] - ringhelp[2]) * (coords2[2] - ringhelp[2]))) / 3;
    }
    if (debug)
    {
        std::cout << "Your three vectors are:\n";
        std::cout << "X= " << x[0] << " " << x[1] << " " << x[2] << "\n";
        std::cout << "Y= " << y[0] << " " << y[1] << " " << y[2] << "\n";
        std::cout << "Z= " << z[0] << " " << z[1] << " " << z[2] << "\n";
        std::cout << "From here on mode_leng and cub are used\n";
    }
    for (int r = 0; r < 3; r++)
    {
        if (res[r] < 0)
        {
            std::cout << "Wrong input in res! Try again!\n";
            return(results);
        }
        if (boxsize[r] < 0)
        {
            std::cout << "Wrong input for box scaling! Try again!\n";
            return(results);
        }
        else if (boxsize[r] >= 50)
        {
            std::cout << "Come on, be realistic!\nTry Again!\n";
            return(results);
        }
    }
    if (!mode_res)
    {
        if (debug) std::cout << "Determining resolution\n";
        if (mode_leng)
        {
            if (debug) std::cout << "mres=true; mleng=true\n";
            switch (mode_sel)
            {
            case 1:
                size[0] = ceil(9 * znorm / ang2bohr) / 10;
                break;
            case 2:
                size[0] = ceil(15 * znorm / ang2bohr) / 10;
                break;
            case 3:
                size[0] = ceil(15 * hnorm / ang2bohr) / 10;
                break;
            case 4:
                size[0] = ceil(30 * hnorm / ang2bohr) / 10;
                break;
            }
            for (int r = 0; r < 3; r++)
            {
                if (boxsize[r] != 0) size[r] = boxsize[r] / ang2bohr;
                else if (r > 0 || cub == true) size[r] = size[0];
                s2[r] = size[r] / 2;
            }
        }
        else
        {
            if (debug) std::cout << "mres=true; mleng=false\n";
            for (int r = 0; r < 3; r++)
            {
                if (cub && r > 0)
                {
                    if (debug) std::cout << "I'm making things cubic!\n";
                    s2[r] = s2[0];
                    continue;
                }
                switch (mode_sel)
                {
                case 1:
                case 2:
                    size[r] = znorm / ang2bohr;
                    break;
                case 3:
                case 4:
                    size[r] = hnorm / ang2bohr;
                    break;
                }
                if (boxsize[r] != 0) size[r] = size[r] * boxsize[r];
                else {
                    switch (mode_sel)
                    {
                    case 1:
                        size[r] = ceil(9 * znorm / ang2bohr) / 10;
                        break;
                    case 2:
                        size[r] = ceil(15 * znorm / ang2bohr) / 10;
                        break;
                    case 3:
                        size[r] = ceil(15 * hnorm / ang2bohr) / 10;
                        break;
                    case 4:
                        size[r] = ceil(30 * hnorm / ang2bohr) / 10;
                        break;
                    }
                }
                s2[r] = size[r] / 2;
            }
        }
        for (int r = 0; r < 3; r++) np[r] = (int)round((2 * s2[r]) / res[r]) + 1;
    }
    else
    {
        if (debug) std::cout << "Boxsize is used\n";
        for (int r = 0; r < 3; r++)
        {
            if (cub) np[r] = (int)res[0];
            else np[r] = (int)res[r];
        }
        if (mode_leng)
        {
            if (debug) std::cout << "mode_leng=true; using boxsize\n";
            switch (mode_sel)
            {
            case 1:
                size[0] = ceil(9 * znorm / ang2bohr) / 10;
                break;
            case 2:
                size[0] = ceil(15 * znorm / ang2bohr) / 10;
                break;
            case 3:
                size[0] = ceil(15 * hnorm / ang2bohr) / 10;
                break;
            case 4:
                size[0] = ceil(30 * hnorm / ang2bohr) / 10;
                break;
            }
            for (int r = 0; r < 3; r++)
            {
                if (debug) std::cout << r + 1 << ". Axis:";
                if (boxsize[r] != 0) size[r] = boxsize[r];
                else if (r > 0 || cub == true) size[r] = size[0];
                s2[r] = size[r] / 2;
            }
        }
        else
        {
            if (debug) std::cout << "Enter multiplicator for length (z,y,x)\n";
            for (int r = 0; r < 3; r++)
            {
                if (cub == 1 && r > 0)
                {
                    if (debug) std::cout << "I'm making things cubic!\n";
                    s2[r] = s2[0];
                    size[r] = size[0];
                    continue;
                }
                if (debug) std::cout << r + 1 << ". Axis:";
                switch (mode_sel)
                {
                case 1:
                case 2:
                    size[r] = znorm / ang2bohr;
                    break;
                case 3:
                case 4:
                    size[r] = hnorm / ang2bohr;
                    break;
                }
                if (boxsize[r] != 0) size[r] = size[r] * boxsize[r];
                else {
                    switch (mode_sel)
                    {
                    case 1:
                        size[r] = ceil(9 * znorm / ang2bohr) / 10;
                        break;
                    case 2:
                        size[r] = ceil(15 * znorm / ang2bohr) / 10;
                        break;
                    case 3:
                        size[r] = ceil(15 * hnorm / ang2bohr) / 10;
                        break;
                    case 4:
                        size[r] = ceil(30 * hnorm / ang2bohr) / 10;
                        break;
                    }
                }
                s2[r] = size[r] / 2;
            }
        }
    }
    double o[3] = {};
    if (debug) std::cout << "This is the origin:\n";
    for (int d = 0; d < 3; d++)
    {
        switch (mode_sel)
        {
        case 1:
            o[d] = coords2[d] - (s2[0] * x[d] + s2[1] * y[d] + s2[2] * z[d]) / ang2bohr;
            break;
        case 2:
            o[d] = (coords1[d] - (coords1[d] - coords2[d]) / 2) - (s2[0] * x[d] + s2[1] * y[d] + s2[2] * z[d]) / ang2bohr;
            break;
        case 3:
            o[d] = (coords3[d] - (coords3[d] - coords1[d]) / 2) - (s2[0] * x[d] + s2[1] * y[d] + s2[2] * z[d]) / ang2bohr;
            break;
        case 4:
            o[d] = ringhelp[d] - (s2[0] * x[d] + s2[1] * y[d] + s2[2] * z[d]) / ang2bohr;
            break;
        }
        if (debug) std::cout << o[d] << "\n";
    }
    std::string outname = { "" };
    outname += wavy.get_path().generic_string();
    outname += "_";
    outname += label[0];
    outname += std::to_string(atom1);
    outname += "_";
    outname += label[1];
    outname += std::to_string(atom2);
    outname += "_";
    outname += label[2];
    outname += std::to_string(atom3);
    outname += "_";
    outname += std::to_string(runnumber);
    double v[9];
    double incr[3];
    for (int i = 0; i < 3; i++) {
        v[i] = x[i];
        v[i + 3] = y[i];
        v[i + 6] = z[i];
        incr[i] = size[i] / np[i];
        //incr[i]=res[i];
    }
    if (compute_dens(wavy, debug, np, o, v, incr, outname, rho, rdg, eli, lap) == 0) {
        results.success = true;
        results.filename = outname;
    }
    return(results);
}

int autobonds(bool debug, WFN& wavy, const std::filesystem::path& inputfile, const bool& bohr) {
    char inputFile[2048] = "";
    if (!exists(inputfile))
    {
        std::cout << "No input file specified! I am going to make an example for you in input.example!\n";
        std::ofstream example("input.example");
        example << "!COMMENT LINES START WITH ! AND CAN ONLY BE ABOVE THE FIRST SWITCHES AND NUMBERS!\n!First row contains the following switches: rho (dens), Reduced Density Gradient (RDG), ELI-d (eli) and Laplacian of rho (lap)\n!Following rows each contain a bond you want to investigate.\n!The key to read these numbers is:\n!orientation_selection(1-4)\n!      1=atom1 centered\n!      2=bond atom1 atom2\n!      3=bond atom1 and atom3\n!      4=ring centroid of the three atoms\n!length selection(0/1)\n!      0=box will contain multiplicator of bondlength between atom1 and atom2\n!      1=box will contain length in angstrom (?)\n!resolution selection(0/1)\n!      0=res will contain number of gridpoints\n!      1=res will contain distance between gridpoints\n!res1 res2 res3\n!      based on selection above resolution in x,y,z direction (either gridpoints or distance between points)\n!cube selection (0/1)\n!      0= all selections are considered\n!      1= all selections in x-direction will be applied in the y and z direction, as well, making it a cube\n!box1 box2 box3\n!      based on selection above size of the box/cube (either multiplicator of bondlength or absolute length)\n!atom1 atom2 atom3\n!      the atoms in your wavefunction file (counting from 1) to be used as references\n! BELOW THE INPUT SECTION STARTS\n!rho, rdg, eli, lap\n!mode_sel(int) mode_leng(bool) mode_res(bool) res1(double) res2(double) res3(double) cube(bool) boxsize1(float) boxsize2(float) boxsize3(float) atom1(int) atom2(int) atom3(int)\n 1    1    1    1\n            2             1               1              20.0         20.0         20.0         0          5.0             5.1             5.2             1          2          3\n";
        example.close();
        return 0;
    }
    std::ifstream input(inputfile.c_str());
    if (!input.good())
    {
        std::cout << inputFile << " does not exist or is not readable!\n";
        return 0;
    }
    input.seekg(0);
    std::string line("");
    getline(input, line);
    std::string comment("!");
    while (line.compare(0, 1, comment) == 0) getline(input, line);
    int rho, rdg, eli, lap;
    if (line.length() < 10) return 0;
    else
    {
        std::istringstream iss(line);
        if (!(iss >> rho >> rdg >> eli >> lap)) {
            std::cerr << "Error parsing line for rho, rdg, eli, and lap values." << std::endl;
            return 0; // Handle the error appropriately
        }
    }
    int errorcount = 0, runnumber = 0;
    do {
        getline(input, line);
        if (line.length() < 10) continue;
        runnumber++;
        int sel = 0, leng = 0, cube = 0, mres = 0, a1 = 0, a2 = 0, a3 = 0;
        double res[3], box[3];
        bool bleng = false, bres = false, bcube = false;
        if (line.length() > 1)
        {
            std::istringstream iss(line);
            if (!(iss >> sel >> leng >> mres >> res[0] >> res[1] >> res[2] >> cube >> box[0] >> box[1] >> box[2] >> a1 >> a2 >> a3)) {
                std::cerr << "Error parsing line for bond parameters." << std::endl;
                continue; // Skip this line and move to the next
            }
        }
        if (leng == 1) { bleng = true; if (debug) std::cout << "leng=true\n"; }
        else { bleng = false; if (debug) std::cout << "leng=false\n"; }
        if (cube == 1) { bcube = true; if (debug) std::cout << "cube=true\n"; }
        else { bcube = false; if (debug) std::cout << "cube=false\n"; }
        if (mres == 1) { bres = true; if (debug) std::cout << "mres=true\n"; }
        else { bres = false; if (debug) std::cout << "mres=false\n"; }
        if (debug) std::cout << "running calculations for line " << runnumber << " of the input file:\n\n";
        bond work{ "","","","",false,false,false };
        work = do_bonds(wavy, sel, bleng, bres, res, bcube, box, a1, a2, a3, debug, bohr, runnumber, rho == 1, rdg == 1, eli == 1, lap == 1);
        if (!work.success)
        {
            std::cout << "!!!!!!!!!!!!problem somewhere during the calculations, see messages above!!!!!\n";
            errorcount++;
        }
        else {
            if (rho == 1) wavy.push_back_cube(work.filename + "_rho.cube", false, false);
            if (rdg == 1) wavy.push_back_cube(work.filename + "_rdg.cube", false, false);
            if (eli == 1) wavy.push_back_cube(work.filename + "_eli.cube", false, false);
            if (lap == 1) wavy.push_back_cube(work.filename + "_lap.cube", false, false);
            if (rho == 1 && rdg == 1) wavy.push_back_cube(work.filename + "_signed_rho.cube", false, false);
        }
    } while (!input.eof());
    std::cout << "\n  *   Finished all calculations! " << runnumber - errorcount << " out of " << runnumber << " were successful!            *\n";
    return 1;
}

std::vector<std::pair<int, int>> get_bonded_atom_pairs(const WFN& wavy) {
    std::vector<std::pair<int, int>> bonds;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        for (int j = i + 1; j < wavy.get_ncen(); j++)
        {
            const double distance = std::hypot(wavy.get_atom_coordinate(i, 0) - wavy.get_atom_coordinate(j, 0),
                wavy.get_atom_coordinate(i, 1) - wavy.get_atom_coordinate(j, 1),
                wavy.get_atom_coordinate(i, 2) - wavy.get_atom_coordinate(j, 2));
            const double svdW = constants::ang2bohr(constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)]);
            if (distance < 1.25 * svdW)
            {
#ifdef NSA2DEBUG
                std::cout << "Bond between " << i << " (" << wavy.get_atom_charge(i) << ") and " << j << " (" << wavy.get_atom_charge(j) << ") with distance " << distance << " and svdW " << svdW << std::endl;
#endif // NSA2DEBUG
                bonds.push_back(std::make_pair(i, j));
            }
        }
    }
    return bonds;
}

//Assuming square matrices
vec change_basis_sq(const vec& in, const vec& transformation, int size) {

    // new = transform^T * in * tranform
    vec result(size * size);
    vec temp(size * size);
    // first we do temp = t^T * i
    cblas_dgemm(CblasRowMajor,
        CblasTrans, CblasNoTrans,
        size, size, size,
        1.0,
        transformation.data(), size,
        in.data(), size,
        0.0,
        temp.data(), size);

    //Then we do res = temp * t
    cblas_dgemm(CblasRowMajor,
        CblasNoTrans, CblasNoTrans,
        size, size, size,
        1.0,
        temp.data(), size,
        transformation.data(), size,
        0.0,
        result.data(), size);
    return result;
}

//For non-square transformation matrices
//performs res = trans * in * trans^T in case of forward
//and res = trans^T * int * trans in case of !forward
dMatrix2 change_basis_general(const dMatrix2& in, const dMatrix2& transformation, bool forward = true) {
    //Checks should be handled by dot itself
    //int a, b, c, d;
    //c = transformation.extent(1); //cols of t
    //d = transformation.extent(0); //rows of t
    //a = in.extent(1); //cols of in
    //b = in.extent(0); //rows of in
    //err_checkf(a == b, "Input matrix must be square.", std::cout);
    //err_checkf(c == b, "Incompatible matrix dimensions for basis change.", std::cout);
    if (!forward) {
        dMatrix2 temp = dot<dMatrix2>(transformation, in, false, false); // temp = t * in
        dMatrix2 result = dot<dMatrix2>(temp, transformation, false, true); // res = temp * t^T
        return result;
    }
    else {
        dMatrix2 temp = dot<dMatrix2>(transformation, in, true, false); // temp = t^T * in
        dMatrix2 result = dot<dMatrix2>(temp, transformation, false, false); // res = temp * t
        return result;
    }
}

/**
 * Calculates Atomic Natural Orbitals for a specific atom.
 *
 * @param D_full      Pointer to the full Density Matrix (N_basis x N_basis)
 * @param S_full      Pointer to the full Overlap Matrix (N_basis x N_basis)
 * @param full_stride The leading dimension of the full matrices (usually N_basis)
 * @param atom_indices A vector containing the indices (0-based) of the basis functions for this atom
 * @return NAOResult containing sorted occupancies and coefficients
 */
Roby_information::NAOResult Roby_information::calculateAtomicNAO(const dMatrix2& D_full,
    const dMatrix2& S_full,
    const std::vector<int>& atom_indices) {

    err_checkf(D_full.extent(0) == D_full.extent(1), "Density matrix D must be square.", std::cout);
    err_checkf(S_full.extent(0) == S_full.extent(1), "Overlap matrix S must be square.", std::cout);
    err_checkf(D_full.extent(0) == S_full.extent(0), "Density and Overlap matrices must be of the same size.", std::cout);

    const int n = static_cast<int>(atom_indices.size());

    // 1. Memory Allocation for Submatrices
    // Using flat std::vectors to ensure contiguous memory for MKL
    vec D_sub(n * n, 0.0); // called P in tonto
    vec S_sub(n * n, 0.0); // called S in tonto
    get_submatrices(D_full, S_full, D_sub, S_sub, atom_indices);
    vec Rho(n * n);        // To store target density

    vec V = S_sub;
    vec W(n);
    // make V = Sqrt(S)
    const vec Temp = mat_sqrt(V, W);

#ifdef NSA2DEBUG
    print_dmatrix2(reshape<dMatrix2>(Temp, Shape2D(n, n)), "projection matrix V");
#endif

    const vec X = change_basis_sq(D_sub, Temp, n);
#ifdef NSA2DEBUG
    print_dmatrix2(reshape<dMatrix2>(X, Shape2D(n, n)), "projected density X");
#endif

    vec occu(n, 0);
    vec P = X;
#ifdef NSA2DEBUG
    err_checkf(isSymmetricViaEigenvalues<vec>(P, n), "Transformed matrix not symmetric!", std::cout);
#endif
    make_Eigenvalues(P, occu);
#ifdef NSA2DEBUG
    std::cout << "Eigenvalues of projected density P:\n";
    for (int i = 0; i < n; ++i) {
        std::cout << std::setw(14) << std::setprecision(8) << std::fixed << W[i] << " ";
    }
    print_dmatrix2(reshape<dMatrix2>(P, Shape2D(n, n)), "Projected density P");
#endif

    for (int i = 0; i < n; ++i)
        W[i] = abs(W[i]) < 1E-10 ? 0.0 : 1.0 / W[i];

    vec Temp2(n * n, 0.0);
    double* T;
    int in, jn;
    for (int i = 0; i < n; i++) {
        in = i * n;
        for (int j = 0; j < n; j++) {
            jn = j * n;
            T = &Temp2[in + j];
#pragma simd
            for (int k = 0; k < n; k++)
                *T += V[in + k] * W[k] * V[jn + k];
        }
    }

    V.clear();
    W.clear();

#ifdef NSA2DEBUG
    print_dmatrix2(reshape<dMatrix2>(Temp2, Shape2D(n, n)), "Back projection S^-0.5");
#endif

    cblas_dgemm(CblasRowMajor,
        CblasNoTrans, CblasNoTrans,
        n, n, n,
        1.0,
        Temp2.data(), n,
        P.data(), n,
        0.0,
        Rho.data(), n);

#ifdef NSA2DEBUG
    print_dmatrix2(reshape<dMatrix2>(Rho, Shape2D(n, n)), "resulting NAO");
#endif

    NAOResult result;
    result.eigenvalues = occu;
    result.eigenvectors = Rho; // Rho now contains eigenvectors
    result.matrix_elements = atom_indices;
    result.sub_DM = D_sub;
    result.sub_OM = S_sub;

    // Create an index vector to sort descending
    ivec idx(n);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), [&](int i1, int i2) {
        return result.eigenvalues[i1] > result.eigenvalues[i2];
        });

    // Reorder based on sorted indices
    vec sorted_evals;
    vec sorted_evecs;
    sorted_evecs.reserve(n * n);

    for (int i = 0; i < n; ++i) {
        int original_idx = idx[i];
        if (result.eigenvalues[original_idx] < 0.075)
            continue;
        sorted_evals.emplace_back(result.eigenvalues[original_idx]);

        for (int row = 0; row < n; ++row) {
            sorted_evecs.emplace_back(result.eigenvectors[row * n + original_idx]);
        }
    }

    result.eigenvalues = sorted_evals;
    result.eigenvectors = sorted_evecs;

    return result;
}

/**
 * Computes Final NAOs by symmetrically orthogonalizing the Pre-NAOs.
 *
 * @param C_PNAO    Global matrix of Pre-NAOs (N x N). Columns are orbitals.
 * @param S_AO      Original Overlap Matrix in AO basis (N x N).
 * @param n         Number of basis functions (N).
 * @return          Matrix of final NAOs (N x N), Columns are orbitals.
 */

 /* CURRENTLY NOT IN USE
 std::vector<double> orthogonalizePNAOs(const vec& C_PNAO,
     const vec& S_AO,
     int n) {

     vec S_PNAO(n * n);
     vec Temp(n * n); // Intermediate buffer
     vec C_NAO(n * n); // Result

     // 1. Compute Overlap in PNAO basis: S_PNAO = C_PNAO^T * S_AO * C_PNAO

     // Step A: Temp = S_AO * C_PNAO
     // S_AO is symmetric.
     cblas_dsymm(CblasRowMajor,
         CblasLeft, CblasUpper,
         n, n,
         1.0, S_AO.data(), n,
         C_PNAO.data(), n,
         0.0, Temp.data(), n);

     // Step B: S_PNAO = C_PNAO^T * Temp
     // C_PNAO is not symmetric, so we use dgemm.
     // Transpose the first matrix (C_PNAO^T).
     cblas_dgemm(CblasRowMajor,
         CblasTrans, CblasNoTrans,
         n, n, n,
         1.0, C_PNAO.data(), n,
         Temp.data(), n,
         0.0, S_PNAO.data(), n);

     // 2. Compute S_PNAO^(-1/2) using Eigendecomposition
     // S_PNAO is symmetric (and positive definite).

     vec W(n); // Eigenvalues
     // We can overwrite S_PNAO with eigenvectors to save memory,
     // but let's keep it clear. Copy S_PNAO to 'U' (Eigenvectors).
     vec U = S_PNAO;

     // LAPACKE_dsyevd: Computes all eigenvalues and eigenvectors
     err_checkf(LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, U.data(), n, W.data()) == 0, "Eigenvalue computation failed.", std::cout);

     // Construct S^-1/2 = U * Lambda^(-1/2) * U^T
     std::fill(S_PNAO.begin(), S_PNAO.end(), 0.0); // Reuse S_PNAO to store S^-1/2

     for (int k = 0; k < n; ++k) {
         double scale = 1.0 / std::sqrt(W[k]);
         // Add contribution of k-th eigenvector: scale * (v_k * v_k^T)
         // v_k is the k-th ROW of U.
         const double* v_k = &U[k * n];

         // This is a rank-1 update (dger), but we can just sum manually or loop
         // Since we need the full matrix for the next multiplication
         for (int i = 0; i < n; ++i) {
             for (int j = 0; j < n; ++j) {
                 S_PNAO[i * n + j] += scale * v_k[i] * v_k[j];
             }
         }
     }

     // 3. Transform PNAOs to NAOs
     // C_NAO = C_PNAO * S^(-1/2)
     // C_PNAO (Columns are PNAOs) * S_inv_sqrt (transformation matrix)

     cblas_dgemm(CblasRowMajor,
         CblasNoTrans, CblasNoTrans,
         n, n, n,
         1.0, C_PNAO.data(), n,
         S_PNAO.data(), n, // This is now S^-1/2
         0.0,
         C_NAO.data(), n);

     return C_NAO;
 }
 */

double Roby_information::projection_matrix_and_expectation(const ivec& indices, const ivec& eigvals, const ivec& eigvecs, dMatrix2* given_NAO) {
    const int n = indices.size();
    //vec D_Sub(n * n, 0.0);
    vec S_Sub(n * n, 0.0);
    get_submatrix(overlap_matrix, S_Sub, indices);
    dMatrix2 S = reshape<dMatrix2>(S_Sub, Shape2D(n, n));
    int atom = -1;
    //dMatrix2 D = reshape<dMatrix2>(D_Sub, Shape2D(n, n));

    dMatrix2 NAOs;
    if (n == total_NAOs.extent(1))
        NAOs = total_NAOs;
    else {
        //TODO: assign subspace NAOs from NAOResults for a given atom
        for (auto NAO : this->NAOs) {
            // if matrix_elemnts of this NAO are identical to indices, resahpe NAO.eigenvectors to the correct shape
            if (NAO.matrix_elements == indices)
                NAOs = reshape<dMatrix2>(NAO.eigenvectors, Shape2D(NAO.eigenvalues.size(), n));
            atom = NAO.atom_index;
        }
    }
    //For atom groups we can use submatrices of total_NAOs
    if (NAOs.size() == 0) {
        const int n1 = eigvals.size();
        const int n2 = eigvecs.size();
        vec NAO_sub(n1 * n2);
        if (given_NAO == nullptr)
            given_NAO = &total_NAOs;
        get_submatrix(*given_NAO, NAO_sub, eigvals, eigvecs);
        NAOs = reshape<dMatrix2>(NAO_sub, Shape2D(n1, n2));
    }
#ifdef NSA2DEBUG
    print_dmatrix2(NAOs, "W in projection matrix making");
#endif

    auto X = change_basis_general(S, transpose(NAOs));
#ifdef NSA2DEBUG
    print_dmatrix2(X, "new Basis");
#endif

    auto Y = LAPACKE_invert(X);
#ifdef NSA2DEBUG
    print_dmatrix2(Y, "Pseudo inverse of Y");
#endif

    //making the projection matrix
    X = change_basis_general(Y, transpose(NAOs), false);
#ifdef NSA2DEBUG
    print_dmatrix2(X, "Backtransformed X");
#endif

    if (atom >= 0) {
        projection_matrices.push_back(X);
        overlap_matrices.push_back(S);
    }

    dMatrix2 S_rect = get_rectangle(overlap_matrix, indices);
#ifdef NSA2DEBUG
    print_dmatrix2(S_rect, "S_rect:");
#endif

    //overlap projection
    auto W = change_basis_general(X, S_rect);
#ifdef NSA2DEBUG
    print_dmatrix2(W, "Overlap transformed W");
#endif

    const double expect = trace_product<double>(W, density_matrix);
    return expect;
}

double Roby_information::Roby_population_analysis(const ivec atoms) {
    ivec bf_indices;
    if (atoms.size() == 0) {
        for (NAOResult& NAO : NAOs) {
            for (auto index : NAO.matrix_elements)
                bf_indices.push_back(index);
        }
    }
    else {
        bf_indices = atoms;
    }
    double P = projection_matrix_and_expectation(bf_indices);
    return P;
}

void Roby_information::computeAllAtomicNAOs(WFN& wavy) {
    const int N_atoms = wavy.get_ncen();
    const std::vector<atom> ats = wavy.get_atoms();
    NAOs.reserve(N_atoms);

    density_matrix = wavy.get_dm();

    Int_Params basis(wavy);
    vec S_full;
    if (wavy.get_d_f_switch())
        compute2c_Overlap_Cart(basis, S_full);
    else
        compute2C<Overlap2C>(basis, S_full);

    overlap_matrix = reshape<dMatrix2>(S_full, Shape2D(density_matrix.extent(0), density_matrix.extent(1)));

#ifdef NSA2DEBUG
    print_dmatrix2(overlap_matrix, "Overlap matrix");
#endif

    //err_checkf()

    int last_index = 0;
    ivec2 indices(wavy.get_ncen());
    for (auto& a : ats) {
        indices[a.get_nr() - 1].reserve(density_matrix.extent(0) / N_atoms); // Rough estimate
        int current_shell = -1;
        int nr_indices = 0;
        std::vector<basis_set_entry> basis_set = a.get_basis_set();
        for (auto& bf : basis_set) {
            if (bf.get_shell() != current_shell) {
                current_shell++;
                if (wavy.get_origin() != e_origin::tonto)
                    bf.get_type() == 1 ? nr_indices = 1 : (bf.get_type() == 2 ? nr_indices = 3 : (bf.get_type() == 3 ? nr_indices = 5 : nr_indices = 7));
                else {
                    bf.get_type() == 1 ? nr_indices = 1 : (bf.get_type() == 2 ? nr_indices = 3 : (bf.get_type() == 3 ? nr_indices = 6 : nr_indices = 10));
                    if (bf.get_type() == 3) {
                        //2 - 4; 3 - 6; 5 - 6
                        swap_rows_cols_symm(overlap_matrix, last_index + 1, last_index + 3);
                        swap_rows_cols_symm(overlap_matrix, last_index + 2, last_index + 5);
                        swap_rows_cols_symm(overlap_matrix, last_index + 4, last_index + 5);
                    }
                }
                for (int i = 0; i < nr_indices; i++) {
                    indices[a.get_nr() - 1].push_back(last_index);
                    last_index++;
                }
            }
        }

        //pNAO.eigenvectors = orthogonalizePNAOs(pNAO.eigenvectors, S_full, static_cast<int>(indices[a.get_nr()].size()));
        NAOs.emplace_back(calculateAtomicNAO(density_matrix, overlap_matrix, indices[a.get_nr() - 1]));
        NAOs.back().atom_index = a.get_nr() - 1;
    }
#ifdef NSA2DEBUG
    print_dmatrix2(overlap_matrix, "Overlap matrix repaired");
#endif
}

ivec Roby_information::find_eigenvalue_pairs(const vec& eigvals, const double tolerance) {
    const int n = eigvals.size();
    ivec pairs(n, -1);
    for (int i = 0; i < n; i++) {
        if (pairs[i] >= 0)
            continue; // already paired
        for (int j = 0; j < n; j++) {
            if (pairs[j] >= 0)
                continue; // already paired
            if (abs(abs(eigvals[j]) - 1.0) < tolerance)
                continue; // skip values close ot +/- one, as they reside only on one atom
            if (abs(eigvals[i] + eigvals[j]) < tolerance) {
                pairs[i] = j;
                pairs[j] = i;
                break;
            }
        }
        if (pairs[i] == -1)
            pairs[i] = i;
    }
    return pairs;
}

void Roby_information::transform_Ionic_eigenvectors_to_Ionic_orbitals(
    dMatrix2& EVC,
    const vec& eigvals,
    const ivec& pairs,
    const int index_a,
    const int index_b,
    const ivec& pair_matrix_indices)
{
    double fp, fm, fa, fb, s, c, s2;
    const int n_ab = EVC.extent(0);
    const int n_eigvals = eigvals.size();
    const int n_a = projection_matrices[index_a].extent(0);
    const int n_b = projection_matrices[index_b].extent(0);
    err_checkf(n_ab == n_a + n_b, "Inconsitent size in projection matrices?!", std::cout);
    err_checkf(n_eigvals == EVC.extent(1), "Inconsitency between EVC and eigvals", std::cout);

    dMatrix1 A(n_a), B(n_b);
    dMatrix1 EVC_column(n_ab);
    dMatrix2 PAS(n_a, n_ab), PBS(n_b, n_ab);

    vec Sub_overlap(n_a * n_ab);
    get_submatrix(overlap_matrix, Sub_overlap, NAOs[index_a].matrix_elements, pair_matrix_indices);
    dMatrix2 Sa = reshape<dMatrix2>(Sub_overlap, Shape2D(n_a, n_ab));
    PAS = dot<dMatrix2>(projection_matrices[index_a], Sa, false, false);
    Sub_overlap.clear(); Sub_overlap.resize(n_b * n_ab);
    get_submatrix(overlap_matrix, Sub_overlap, NAOs[index_b].matrix_elements, pair_matrix_indices);
    dMatrix2 Sb = reshape<dMatrix2>(Sub_overlap, Shape2D(n_b, n_ab));
    PBS = dot<dMatrix2>(projection_matrices[index_b], Sb, false, false);
#ifdef NSA2DEBUG
    print_dmatrix2(PAS, "PAS");
    print_dmatrix2(PBS, "PBS");
#endif // NSA2DEBUG


    for (int i = 0; i < n_eigvals; i++) {
        if (pairs[i] < 0) continue;
        if (pairs[i] == i) continue;
        if (eigvals[i] < eigvals[pairs[i]]) continue;
#ifdef NSA2DEBUG
        std::cout << "Doing i=" << i << std::endl;
#endif
        for (int a = 0; a < n_ab; a++)
            EVC_column(a) = EVC(a, i);

#ifdef NSA2DEBUG
        print_dmatrix2(reshape<dMatrix2>(EVC_column, Shape2D(n_ab, 1)), "slice used");
#endif

        s = eigvals[i];
        s2 = s * s;
        if (abs(s2 - 1.0) < 1E-8) c = 0.0;
        else c = sqrt(1.0 - s2);

        fm = sqrt(1 - c) / s;
        fp = sqrt(1 + c) / s;
        fa = 0.5 * ((fm + fp) + c * (fm - fp));
        fb = 0.5 * (c * (fm + fp) + (fm - fp));

        if (abs(fa - 1.0) > 1E-8) {
            A = dot_BLAS<dMatrix1, dMatrix2>(PAS, EVC_column, false);
            for (int a = 0; a < n_a; a++)
                A(a) /= fa;
        }
        if (abs(fa - 1.0) > 1E-8) {
            B = dot_BLAS<dMatrix1, dMatrix2>(PBS, EVC_column, false);
            for (int b = 0; b < n_b; b++)
                B(b) /= fb;
        }

#ifdef NSA2DEBUG
        std::cout << "fa: " << fa << std::endl << "fb: " << fb << std::endl;
        print_dmatrix2(reshape<dMatrix2>(A, Shape2D(n_a, 1)), "A");
        print_dmatrix2(reshape<dMatrix2>(B, Shape2D(n_b, 1)), "B");
#endif

        //build antibonding state in pairs[i]
        fa = 0.5 * (fm - fp);
        fb = 0.5 * (fm + fp);

        for (int a = 0; a < n_a; a++) {
            EVC(a, pairs[i]) = fa * A(a);
        }
        for (int b = n_a; b < n_ab; b++) {
            EVC(b, pairs[i]) = fb * B(b - n_a);
        }
    }
}

std::map<char, dMatrix2> Roby_information::make_covalent_from_ionic(
    const dMatrix2& theta_I,
    const vec& eigvals,
    const ivec& pairs) {
    std::map<char, dMatrix2> res; // A = angle, V = eigen_value, T = Theta_vector
    const int size = eigvals.size();
    res.emplace('A', dMatrix2(size, 1));
    res.emplace('V', dMatrix2(size, 1));
    res.emplace('T', dMatrix2(theta_I.extent(0), theta_I.extent(1)));

    for (int i = 0; i < size; i++) {
        res['A'](i, 0) = 90.0;
        res['V'](i, 0) = 0.0;
    }

    for (int val = 0; val < size; val++) {
        if (pairs[val] < 0) continue;
        if (pairs[val] == val) continue;
        const double s = eigvals[val];
        if (s < eigvals[pairs[val]]) continue;

        for (int i = 0; i < theta_I.extent(0); i++) {
            res['T'](i, val) = constants::INV_SQRT2 * (theta_I(i, val) + theta_I(i, pairs[val]));
            res['T'](i, pairs[val]) = constants::INV_SQRT2 * (theta_I(i, val) - theta_I(i, pairs[val]));
        }

        const double c = sqrt(1.0 - s * s);
        res['V'](val, 0) = c;
        res['V'](pairs[val], 0) = -c;

        res['A'](val, 0) = atan2(s, c) * constants::INV_PI_180;
    }

    return res;
}

Roby_information::Roby_information(WFN& wavy) {
    auto bonds = get_bonded_atom_pairs(wavy);
    std::cout << "Calculating NAOs for all atoms...                 " << std::flush;
    computeAllAtomicNAOs(wavy);
    std::cout << " ...done!" << std::endl;
    Shape2D NAOs_size;
    for (size_t atom_idx = 0; atom_idx < NAOs.size(); atom_idx++) {
#ifdef NSA2DEBUG
        std::cout << "Atom " << atom_idx + 1 << " NAO Occupancies:\n";
#endif
        for (size_t i = 0; i < NAOs[atom_idx].eigenvalues.size(); i++) {
#ifdef NSA2DEBUG
            std::cout << "  NAO " << i + 1 << ": " << std::setprecision(6) << NAOs[atom_idx].eigenvalues[i] << "\n";
#endif
            NAOs_size.rows++;
        }
        const int cols = NAOs[atom_idx].eigenvectors.size() / NAOs[atom_idx].eigenvalues.size();
        NAOs_size.cols += cols;
#ifdef NSA2DEBUG
        std::cout << std::endl;
#endif
    }
    //Create a matrix of size n_NAOs x n_bf and fill with subblocks of the evecs
    vec NAO_matrix(NAOs_size.cols * NAOs_size.rows, 0.0);
    total_NAOs = reshape<dMatrix2>(NAO_matrix, Shape2D(NAOs_size.rows, NAOs_size.cols));
    Shape2D temp = { 0,0 };
    for (auto NAO : NAOs) {
        const int ME = NAO.matrix_elements.size();
        const int NEV = NAO.eigenvalues.size();
        for (int i = 0; i < NEV; i++) {
            const int row = temp.rows + i;
            for (int j = 0; j < ME; j++) {
                const int index_eigenvector = i * ME + j;
                const int col = temp.cols + j;
                total_NAOs(row, col) = NAO.eigenvectors[index_eigenvector];
            }
        }
        temp.rows += NAO.eigenvalues.size();
        temp.cols += NAO.matrix_elements.size();
    }
#ifdef NSA2DEBUG
    print_dmatrix2(total_NAOs, "Global NAO Matrix");
#endif
    // total_NAOs has n_NAOs rows (each row is an individual NAO in the basis set)
    // and n_bfs columns, where each value corresponds to a basis function from here on
    std::cout << "Calculating Population for all atoms...           " << std::flush;
    const double all_atom_population = Roby_population_analysis({});
    std::cout << " ...done." << std::endl;
    std::cout << std::endl << "Total Population: " << all_atom_population << "\n\n";
    vec atom_pops(NAOs.size(), 0.0);
#ifndef NSA2DEBUG
    ProgressBar* pb = new ProgressBar(NAOs.size(), 40, "-", " ", "Calculating Atomic Populations");
#endif
    for (auto NAO : NAOs) {
        atom_pops[NAO.atom_index] = Roby_population_analysis(NAO.matrix_elements);
#ifndef NSA2DEBUG
        pb->update();
#endif
    }
#ifndef NSA2DEBUG
    delete pb;
#endif
    for (int i = 0; i < atom_pops.size(); i++) {
        std::cout << "Population of atom " << i << ": " << atom_pops[i] << std::endl;
    }

#ifndef NSA2DEBUG
    pb = new ProgressBar(bonds.size(), 40, "-", " ", "Calculating Bond Populations");
#endif
    //now perform bond analysis for all bonded atoms
    for (auto bond : bonds) {
        //std::cout << std::endl << "---------------------------- Atom Pair: " << bond.first << " " << bond.second << " ----------------------\n";
        ivec bond_indices, bond_eigenvecs, bond_eigenvals;
        //gather basis function indices for both atoms
        for (auto idx : NAOs[bond.first].matrix_elements)
            bond_indices.push_back(idx);
        for (auto idx : NAOs[bond.second].matrix_elements)
            bond_indices.push_back(idx);
        //just in case: sort bond indices, easy since each atom's indices are already assumed sorted and each basis function belongs to only one atom once
        std::sort(bond_indices.begin(), bond_indices.end());

        //now determine the bond atom NAO indices
        int start_vec = 0, start_val = 0;
        for (auto NAO : NAOs) {
            if (NAO.atom_index == bond.first || NAO.atom_index == bond.second) {
                const int n1 = NAO.eigenvalues.size();
                const int n2 = NAO.eigenvectors.size() / n1;
                for (int i = 0; i < n1; i++) {
                    bond_eigenvals.push_back(start_val + i);
                }
                for (int i = 0; i < n2; i++) {
                    bond_eigenvecs.push_back(start_vec + i);
                }
            }
            start_vec += NAO.matrix_elements.size();
            start_val += NAO.eigenvalues.size();
        }

        //calcualte population using data from both atoms
        const double bond_population = projection_matrix_and_expectation(bond_indices, bond_eigenvals, bond_eigenvecs);
        //atom_pair_populations(bond.first, bond.second) = bond_population;
        //atom_pair_populations(bond.second, bond.first) = bond_population;
        std::cout << "Bond population between atom " << bond.first + 1 << " and atom " << bond.second + 1 << ": " << bond_population << std::endl;
        const int size_ion1 = projection_matrices[bond.first].extent(0) + projection_matrices[bond.second].extent(0),
            size_ion2 = projection_matrices[bond.first].extent(1) + projection_matrices[bond.second].extent(1),
            size1_1 = projection_matrices[bond.first].extent(0),
            size1_2 = projection_matrices[bond.first].extent(1),
            size2_1 = projection_matrices[bond.second].extent(0),
            size2_2 = projection_matrices[bond.second].extent(1);
        dMatrix2 Ionic_Operator(size_ion1, size_ion2);
        for (int i = 0; i < size1_1; i++) {
            for (int j = 0; j < size1_2; j++) {
                Ionic_Operator(i, j) = projection_matrices[bond.first](i, j);
            }
        }
        for (int i = 0; i < size2_1; i++) {
            for (int j = 0; j < size2_2; j++) {
                Ionic_Operator(i + size1_1, j + size1_2) = -projection_matrices[bond.second](i, j);
            }
        }
#ifdef NSA2DEBUG
        print_dmatrix2(Ionic_Operator, "Ionic Operator");
#endif
        // get matching suboverlap matrix
        const int n = bond_indices.size();
        vec S_Sub(n * n, 0.0);
        get_submatrix(overlap_matrix, S_Sub, bond_indices);
        //dMatrix2 S = reshape<dMatrix2>(S_Sub, Shape2D(n, n));

        vec V = S_Sub;
        vec W(n);
        // make V = Sqrt(S)
        const vec Temp = mat_sqrt(V, W);

        dMatrix2 A = reshape<dMatrix2>(Temp, Shape2D(n, n));
#ifdef NSA2DEBUG
        print_dmatrix2(A, "Overlap Sqrt SH");
#endif

        dMatrix2 SI = LAPACKE_invert(A);

#ifdef NSA2DEBUG
        print_dmatrix2(SI, "Overlap Pseudo Inverse");
#endif

        auto X = change_basis_general(Ionic_Operator, transpose(A), true);
#ifdef NSA2DEBUG
        print_dmatrix2(X, "Overlap Eigenproblem");
#endif

        // solve symmetric eigenproblem of X
        vec ionic_eigenvals(X.extent(0));
        make_Eigenvalues(X.container(), ionic_eigenvals);

#ifdef NSA2DEBUG
        std::cout << "Ionic eigenvalues between atom " << bond.first + 1 << " and atom " << bond.second + 1 << ":\n";
        for (size_t i = 0; i < ionic_eigenvals.size(); i++) {
            std::cout << std::setw(5) << i + 1 << ": " << std::setw(10) << std::setprecision(6) << ionic_eigenvals[i] << "\n";
        }
        std::cout << std::endl;
#endif

        auto EVC = dot<dMatrix2>(SI, X);
#ifdef NSA2DEBUG
        print_dmatrix2(EVC, "theta_I");
#endif

        ivec non_zero_indices;
        for (int i = 0; i < ionic_eigenvals.size(); i++) {
            if (abs(ionic_eigenvals[i]) > 1E-5)
                non_zero_indices.push_back(i);
        }
        const int n0 = non_zero_indices.size();
#ifdef NSA2DEBUG
        std::cout << "Non zero eigenvalues:\n";
        for (size_t i = 0; i < n0; i++) {
            std::cout << std::setw(5) << i + 1 << ": " << std::setw(10) << std::setprecision(6) << non_zero_indices[i] << "\n";
        }
        std::cout << std::endl;
#endif

        vec pruned_eigvals;
        for (int nzv = 0; nzv < n0; nzv++)
            pruned_eigvals.push_back(ionic_eigenvals[non_zero_indices[nzv]]);
        EVC = transpose(EVC);
        auto EVC2 = transpose(get_rectangle(EVC, non_zero_indices));
#ifdef NSA2DEBUG
        print_dmatrix2(EVC2, "theta_I after pruning:");
#endif
        EVC.container().clear();

        auto pairs = find_eigenvalue_pairs(pruned_eigvals);
#ifdef NSA2DEBUG
        std::cout << "Pairs:\n";
        for (size_t i = 0; i < n0; i++) {
            std::cout << std::setw(3) << i << ": " << std::setw(5) << std::setprecision(6) << pairs[i] << "\n";
        }
        std::cout << std::endl;
#endif

        transform_Ionic_eigenvectors_to_Ionic_orbitals(EVC2, pruned_eigvals, pairs, bond.first, bond.second, bond_indices);
#ifdef NSA2DEBUG
        print_dmatrix2(EVC2, "theta_I after Ionic");
#endif
        //EVC2 is theta_Ionic, now make theta_Covalent
        auto covalent_info = make_covalent_from_ionic(EVC2, pruned_eigvals, pairs);
#ifdef NSA2DEBUG
        print_dmatrix2(covalent_info['A'], "theta_Angles");
        print_dmatrix2(covalent_info['V'], "eigen_C");
        print_dmatrix2(covalent_info['T'], "theta_C");
#endif
        vec covalent_popul(n0), ionic_popul(n0);
        ivec vals;
        for (int i = 0; i < EVC2.extent(0); i++)
            vals.emplace_back(i);
        for (int i = 0; i < n0; i++) {
            //make covalent populations
            auto temp = transpose(covalent_info['T']);
            covalent_popul[i] = projection_matrix_and_expectation(bond_indices, { i }, vals, &(temp));
            temp = transpose(EVC2);
            ionic_popul[i] = projection_matrix_and_expectation(bond_indices, { i }, vals, &(temp));
        }
#ifdef NSA2DEBUG
        std::cout << "Covalent populations:\n";
        for (int i = 0; i < n0; i++)
            std::cout << "\t" << i << ":\t" << covalent_popul[i] << std::endl;
        std::cout << "Ionic populations:\n";
        for (int i = 0; i < n0; i++)
            std::cout << "\t" << i << ":\t" << ionic_popul[i] << std::endl;
#endif

        const double zero_angle_cutoff = 1E-2 * constants::INV_PI_180;

        vec cov_index(n0, 0.0), ion_index(n0, 0.0);
        double b_ab = 0;
        bond_index_result results(0, 0, 0, 0, 0, 0, 0, 0, {});
        for (int i = 0; i < n0; i++) {
            if (covalent_info['A'](i, 0) < zero_angle_cutoff || covalent_info['A'](i, 0) > 90 - zero_angle_cutoff)
                continue; // skip lone pairs

            if (pruned_eigvals[i] < pruned_eigvals[pairs[i]])
                continue; // skip antibonding

            if (pairs[i] != i) {
                cov_index[i] = 0.5 * (covalent_popul[i] - covalent_popul[pairs[i]]);
                ion_index[i] = 0.5 * (ionic_popul[i] - ionic_popul[pairs[i]]);
                results.covalent += cov_index[i];
                results.ionic += ion_index[i];
            }
            else if (pruned_eigvals[i] > 0.0)
                ion_index[i] = 0.5 * ionic_popul[i];
            else if (pruned_eigvals[i] < 0.0)
                ion_index[i] = -0.5 * ionic_popul[i];
        }

        b_ab = results.covalent * results.covalent + results.ionic * results.ionic;
        results.percent_covalent_Pyth = 100 * (results.covalent * results.covalent / b_ab);
        b_ab = sqrt(b_ab);
        results.percent_covalent_Arakai = 200 * abs(asin(results.covalent / b_ab)) / constants::PI;

        int el_a = wavy.get_atom_charge(bond.first);
        int el_b = wavy.get_atom_charge(bond.second);
        if (el_a > el_b) {
            results.atom_indices = bond;
            results.atom_element_nr = std::make_pair(el_a, el_b);
            results.total = b_ab;
            results.pair_population = bond_population;
            results.population_first = atom_pops[bond.first];
            results.population_second = atom_pops[bond.second];
        }
        else {
            results.atom_indices = std::make_pair(bond.second, bond.first);
            results.atom_element_nr = std::make_pair(el_b, el_a);
            results.total = b_ab;
            results.pair_population = bond_population;
            results.population_first = atom_pops[bond.second];
            results.population_second = atom_pops[bond.first];
            results.ionic = -results.ionic;
        }
        RGBI.push_back(results);
#ifndef NSA2DEBUG
        pb->update();
#endif
    }
#ifndef NSA2DEBUG
    delete pb;
#endif
    const double number_of_electrons = wavy.get_nr_electrons();

    std::cout << "\nRoby-Gould Bond Indices (RGBI) Analysis\n----------------------------------------------\n";
    std::cout << "Number of electrons in system:         " << number_of_electrons << "\n";
    std::cout << "Number of electrons in Roby Analysis:  " << all_atom_population << "\n";
    std::cout << "Percentage of electrons accounted for: " << std::setprecision(4) << (100.0 * all_atom_population / number_of_electrons) << " %\n";
    std::cout << "----------------------------------------------\n";

    std::cout << "\n\nAtom Nr        Els  "
        << std::setw(8) << "n_A"
        << std::setw(8) << "n_B"
        << std::setw(8) << "n_AB"
        << std::setw(8) << "s_AB"
        << std::setw(8) << "Cov."
        << std::setw(8) << "Ion."
        << std::setw(8) << "Tot."
        << std::setw(8) << "Pyth."
        << std::setw(8) << "Arak."
        << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------\n";

    //Sort results by Element of the heavier atom of the bonds and then within each group by the heavier of the second atom
    std::sort(RGBI.begin(), RGBI.end(), [](const bond_index_result& a, const bond_index_result& b) {
        if (a.atom_element_nr.first != b.atom_element_nr.first) {
            return a.atom_element_nr.first > b.atom_element_nr.first;
        }
        return a.atom_element_nr.second > b.atom_element_nr.second;
        });

    for (auto res : RGBI) {
        std::cout << std::setw(4) << res.atom_indices.first << " -" << std::setw(4) << res.atom_indices.second << "  "
            << std::setw(3) << constants::atnr2letter(res.atom_element_nr.first) << " -" << std::setw(3) << constants::atnr2letter(res.atom_element_nr.second)
            << std::fixed << std::setprecision(3) << std::setw(8) << res.population_first
            << std::fixed << std::setprecision(3) << std::setw(8) << res.population_second
            << std::fixed << std::setprecision(3) << std::setw(8) << res.pair_population
            << std::fixed << std::setprecision(3) << std::setw(8) << res.population_first + res.population_second - res.pair_population
            << std::fixed << std::setprecision(3) << std::setw(8) << res.covalent
            << std::fixed << std::setprecision(3) << std::setw(8) << res.ionic
            << std::fixed << std::setprecision(3) << std::setw(8) << res.total
            << std::fixed << std::setprecision(3) << std::setw(8) << res.percent_covalent_Pyth
            << std::fixed << std::setprecision(3) << std::setw(8) << res.percent_covalent_Arakai << std::endl;
    }
    std::cout << "--------------------------------------------------------------------------------------------\n";

    double null = 0;
}
