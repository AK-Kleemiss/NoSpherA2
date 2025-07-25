#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "bondwise_analysis.h"
#include "properties.h"

int compute_dens(WFN& wavy, bool debug, int* np, double* origin, double* gvector, double* incr, std::string& outname, bool rho, bool rdg, bool eli, bool lap) {
    options opt;

    cube CubeRho(np[0], np[1], np[2], wavy.get_ncen(), rho),
        CubeRDG(np[0], np[1], np[2], wavy.get_ncen(), rdg),
        CubeEli(np[0], np[1], np[2], wavy.get_ncen(), eli),
        CubeLap(np[0], np[1], np[2], wavy.get_ncen(), lap);

    std::string Oname = outname;
    std::vector<int> ntyp;
    for (int i = 0; i < 3; i++) {
        opt.NbSteps[i] = np[i];
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
    std::cout << "  *.  NbSteps                   : " << std::setw(11) << opt.NbSteps[0] << std::setw(12) << opt.NbSteps[1] << std::setw(12) << opt.NbSteps[2] << "     *.\n";
    std::cout << "  *.                                                                      *.\n";
    std::cout << "  *.  Number of primitives      :     " << std::setw(5) << wavy.get_nex() << "                               *.\n";
    std::cout << "  *.  Number of MOs             :       " << std::setw(3) << wavy.get_nmo() << "                               *.\n";

    cube dummy(0,0,0);
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

bond do_bonds(WFN &wavy, 
    int mode_sel, 
    bool mode_leng, 
    bool mode_res, 
    double res[], 
    bool cub, 
    double boxsize[], 
    int atom1, 
    int atom2, 
    int atom3, 
    const bool& debug, 
    const bool& bohr, 
    int runnumber, 
    bool rho, 
    bool rdg, 
    bool eli,
    bool lap){
    bond results{"","","","",false,false,false};
    std::string line("");
    std::vector <std::string> label;
    label.resize(3);
    double coords1[3],coords2[3],coords3[3];
    int na=0;
    double z[3],x[3],y[3],help[3],size[3],s2[3];
    int np[3];
    double ang2bohr;
    if(bohr)ang2bohr=0.52917720859;
    else ang2bohr=1.0;
    na = wavy.get_ncen();
    if(atom1<=0 || atom2<=0 || atom3<=0 || mode_sel<=0 || mode_sel> 4 || atom1>na || atom2>na || atom3>na || atom1==atom2 || atom2==atom3 || atom1==atom3)
    {
        std::cout << "Invalid selections of atoms or mode_sel, please try again!\n";
        return(results);
    }
    if(debug) std::cout << "No. of Atoms selected: " << atom1 << atom2 << atom2 << "\n";
    for (int i = 0; i < 3; i++)
        coords1[i] = wavy.get_atom_coordinate(atom1 - 1, i);
    for (int i = 0; i < 3; i++)
        coords2[i] = wavy.get_atom_coordinate(atom2 - 1, i);
    for (int i = 0; i < 3; i++)
        coords3[i] = wavy.get_atom_coordinate(atom3 - 1, i);
    label[0]=wavy.get_atom_label(atom1-1);
    label[1]=wavy.get_atom_label(atom2-1);
    label[2]=wavy.get_atom_label(atom3-1);
    if(debug)
    {
        std::cout << "The Atoms found corresponding to your selection are:\n";
        std::cout << label[0]<< " " <<coords1[0]<< " " <<coords1[1] << " " << coords1[2] << "\n";
        std::cout << label[1]<< " " <<coords2[0]<< " " <<coords2[1] << " " << coords2[2] << "\n";
        std::cout << label[2]<< " " <<coords3[0]<< " " <<coords3[1] << " " << coords3[2] << "\n";
    }
    double znorm=sqrt((coords2[0]-coords1[0])*(coords2[0]-coords1[0])+(coords2[1]-coords1[1])*(coords2[1]-coords1[1])+(coords2[2]-coords1[2])*(coords2[2]-coords1[2]));
    z[0]=(coords2[0]-coords1[0])/znorm;
    z[1]=(coords2[1]-coords1[1])/znorm;
    z[2]=(coords2[2]-coords1[2])/znorm;
    help[0]=coords3[0]-coords2[0];
    help[1]=coords3[1]-coords2[1];
    help[2]=coords3[2]-coords2[2];
    y[0]=z[1]*help[2]-z[2]*help[1];
    y[1]=z[2]*help[0]-z[0]*help[2];
    y[2]=z[0]*help[1]-z[1]*help[0];
    double ynorm=sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
    y[0]=y[0]/ynorm;
    y[1]=y[1]/ynorm;
    y[2]=y[2]/ynorm;
    x[0]=(z[1]*y[2]-z[2]*y[1]);
    x[1]=(z[2]*y[0]-z[0]*y[2]);
    x[2]=(z[0]*y[1]-z[1]*y[0]);
    z[0]=-z[0];
    z[1]=-z[1];
    z[2]=-z[2];

    double hnorm=0.0;
    double ringhelp[3]={};
    if(mode_sel==2 || mode_sel==1) hnorm=znorm;
    else if(mode_sel==3) hnorm=sqrt((coords3[0]-coords1[0])*(coords3[0]-coords1[0])+(coords3[1]-coords1[1])*(coords3[1]-coords1[1])+(coords3[2]-coords1[2])*(coords3[2]-coords1[2]));
    else if(mode_sel==4)
    {
        for(int r=0; r<3;r++) ringhelp[r]=(coords3[r]+coords1[r]+coords2[r])/3;
        hnorm=(sqrt( (coords3[0]-ringhelp[0])*(coords3[0]-ringhelp[0])+(coords3[1]-ringhelp[1])*(coords3[1]-ringhelp[1])+(coords3[2]-ringhelp[2])*(coords3[2]-ringhelp[2]))
                +sqrt( (coords1[0]-ringhelp[0])*(coords1[0]-ringhelp[0])+(coords1[1]-ringhelp[1])*(coords1[1]-ringhelp[1])+(coords1[2]-ringhelp[2])*(coords1[2]-ringhelp[2]))
                +sqrt( (coords2[0]-ringhelp[0])*(coords2[0]-ringhelp[0])+(coords2[1]-ringhelp[1])*(coords2[1]-ringhelp[1])+(coords2[2]-ringhelp[2])*(coords2[2]-ringhelp[2])))/3;
    }
    if(debug)
    {
        std::cout << "Your three vectors are:\n";
        std::cout << "X= " << x[0] << " " << x[1] << " " << x[2] << "\n";
        std::cout << "Y= " << y[0] << " " << y[1] << " " << y[2] << "\n";
        std::cout << "Z= " << z[0] << " " << z[1] << " " << z[2] << "\n";
        std::cout << "From here on mode_leng and cub are used\n";
    }
    for(int r=0;r<3;r++)
    {
        if(res[r]<0)
        {
            std::cout << "Wrong input in res! Try again!\n";
                return(results);
            }
        if(boxsize[r]<0)
        {
            std::cout << "Wrong input for box scaling! Try again!\n";
            return(results);
        }
        else if(boxsize[r]>=50)
        {
            std::cout << "Come on, be realistic!\nTry Again!\n";
            return(results);
        }
    }
    if(!mode_res)
    {
        if(debug) std::cout << "Determining resolution\n";
        if(mode_leng)
        {
            if(debug) std::cout << "mres=true; mleng=true\n";
            switch(mode_sel)
            {
                case 1:
                    size[0]=ceil(9*znorm/ang2bohr)/10;
                    break;
                case 2:
                    size[0]=ceil(15*znorm/ang2bohr)/10;
                    break;
                case 3:
                    size[0]=ceil(15*hnorm/ang2bohr)/10;
                    break;
                case 4:
                    size[0]=ceil(30*hnorm/ang2bohr)/10;
                    break;
            }
            for(int r=0;r<3;r++)
            {
                if(boxsize[r]!=0) size[r]=boxsize[r]/ang2bohr;
                else if(r>0||cub==true) size[r]=size[0];
                s2[r]=size[r]/2;
            }
        }
        else
        {
            if(debug) std::cout << "mres=true; mleng=false\n";
            for(int r=0;r<3;r++)
            {
                if(cub&&r>0)
                {
                    if(debug) std::cout << "I'm making things cubic!\n";
                    s2[r]=s2[0];
                    continue;
                }
                switch(mode_sel)
                {
                    case 1:
                    case 2:
                        size[r]=znorm/ang2bohr;
                        break;
                    case 3:
                    case 4:
                        size[r]=hnorm/ang2bohr;
                        break;
                }
                if(boxsize[r]!=0) size[r]=size[r]*boxsize[r];
                else{
                    switch(mode_sel)
                    {
                        case 1:
                            size[r]=ceil(9*znorm/ang2bohr)/10;
                            break;
                        case 2:
                            size[r]=ceil(15*znorm/ang2bohr)/10;
                            break;
                        case 3:
                            size[r]=ceil(15*hnorm/ang2bohr)/10;
                            break;
                        case 4:
                            size[r]=ceil(30*hnorm/ang2bohr)/10;
                            break;
                    }
                }
                s2[r]=size[r]/2;
            }
        }
        for(int r=0;r<3;r++) np[r]=(int) round((2*s2[r])/res[r])+1;
    }
    else
    {
        if(debug) std::cout << "Boxsize is used\n";
        for (int r=0;r<3;r++)
        {
            if(cub) np[r]=(int) res[0];
            else np[r]=(int) res[r];
        }
        if(mode_leng)
        {
            if(debug) std::cout << "mode_leng=true; using boxsize\n";
            switch(mode_sel)
            {
                case 1:
                    size[0]=ceil(9*znorm/ang2bohr)/10;
                    break;
                case 2:
                    size[0]=ceil(15*znorm/ang2bohr)/10;
                    break;
                case 3:
                    size[0]=ceil(15*hnorm/ang2bohr)/10;
                    break;
                case 4:
                    size[0]=ceil(30*hnorm/ang2bohr)/10;
                    break;
            }
            for(int r=0;r<3;r++)
            {
                if(debug) std::cout << r+1 << ". Axis:";
                if(boxsize[r]!=0) size[r]=boxsize[r];
                else if(r>0||cub==true) size[r]=size[0];
                s2[r]=size[r]/2;
            }
        }
        else
        {
            if(debug) std::cout << "Enter multiplicator for length (z,y,x)\n";
            for(int r=0;r<3;r++)
            {
                if(cub==1&&r>0)
                {
                    if(debug) std::cout << "I'm making things cubic!\n";
                    s2[r]=s2[0];
                    size[r]=size[0];
                    continue;
                }
                if(debug) std::cout << r+1 << ". Axis:";
                switch(mode_sel)
                {
                    case 1:
                    case 2:
                        size[r]=znorm/ang2bohr;
                        break;
                    case 3:
                    case 4:
                        size[r]=hnorm/ang2bohr;
                        break;
                }
                if(boxsize[r]!=0) size[r]=size[r]*boxsize[r];
                else{
                    switch(mode_sel)
                    {
                        case 1:
                            size[r]=ceil(9*znorm/ang2bohr)/10;
                            break;
                        case 2:
                            size[r]=ceil(15*znorm/ang2bohr)/10;
                            break;
                        case 3:
                            size[r]=ceil(15*hnorm/ang2bohr)/10;
                            break;
                        case 4:
                            size[r]=ceil(30*hnorm/ang2bohr)/10;
                            break;
                    }
                }
                s2[r]=size[r]/2;
            }
        }
    }
    double o[3]={};
    if(debug) std::cout << "This is the origin:\n";
    for(int d=0;d<3;d++)
    {
        switch(mode_sel)
        {
            case 1:
                o[d]=coords2[d]-(s2[0]*x[d]+s2[1]*y[d]+s2[2]*z[d])/ang2bohr;
                break;
            case 2:
                o[d]=(coords1[d]-(coords1[d]-coords2[d])/2)-(s2[0]*x[d]+s2[1]*y[d]+s2[2]*z[d])/ang2bohr;
                break;
            case 3:
                o[d]=(coords3[d]-(coords3[d]-coords1[d])/2)-(s2[0]*x[d]+s2[1]*y[d]+s2[2]*z[d])/ang2bohr;
                break;
            case 4:
                o[d]=ringhelp[d]-(s2[0]*x[d]+s2[1]*y[d]+s2[2]*z[d])/ang2bohr;
                break;
        }
        if(debug) std::cout << o[d] << "\n";
    }
    std::string outname={""};
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
    for (int i=0; i<3; i++){
        v[i]=x[i];
        v[i+3]=y[i];
        v[i+6]=z[i];
        incr[i]=size[i]/np[i];
        //incr[i]=res[i];
    }
    if(compute_dens(wavy, debug, np, o, v, incr, outname,rho,rdg,eli,lap)==0){
        results.success=true;
        results.filename=outname;
    }
    return(results);
}

int autobonds(bool debug, WFN &wavy, const std::filesystem::path& inputfile, const bool& bohr){
    char inputFile[2048]="";
    if (!exists(inputfile))
    {
        std::cout << "No input file specified! I am going to make an example for you in input.example!\n";
        std::ofstream example("input.example");
        example << "!COMMENT LINES START WITH ! AND CAN ONLY BE ABOVE THE FIRST SWITCHES AND NUMBERS!\n!First row contains the following switches: rho (dens), Reduced Density Gradient (RDG), ELI-d (eli) and Laplacian of rho (lap)\n!Following rows each contain a bond you want to investigate.\n!The key to read these numbers is:\n!orientation_selection(1-4)\n!      1=atom1 centered\n!      2=bond atom1 atom2\n!      3=bond atom1 and atom3\n!      4=ring centroid of the three atoms\n!length selection(0/1)\n!      0=box will contain multiplicator of bondlength between atom1 and atom2\n!      1=box will contain length in angstrom (?)\n!resolution selection(0/1)\n!      0=res will contain number of gridpoints\n!      1=res will contain distance between gridpoints\n!res1 res2 res3\n!      based on selection above resolution in x,y,z direction (either gridpoints or distance between points)\n!cube selection (0/1)\n!      0= all selections are considered\n!      1= all selections in x-direction will be applied in the y and z direction, as well, making it a cube\n!box1 box2 box3\n!      based on selection above size of the box/cube (either multiplicator of bondlength or absolute length)\n!atom1 atom2 atom3\n!      the atoms in your wavefunction file (counting from 1) to be used as references\n! BELOW THE INPUT SECTION STARTS\n!rho, rdg, eli, lap\n!mode_sel(int) mode_leng(bool) mode_res(bool) res1(double) res2(double) res3(double) cube(bool) boxsize1(float) boxsize2(float) boxsize3(float) atom1(int) atom2(int) atom3(int)\n 1    1    1    1\n            2             1               1              20.0         20.0         20.0         0          5.0             5.1             5.2             1          2          3\n";
        example.close();
        return 0;
    }
    std::ifstream input(inputfile.c_str());
    if(!input.good())
    {
        std::cout << inputFile << " does not exist or is not readable!\n";
        return 0;
    }
    input.seekg(0);
    std::string line("");
    getline(input,line);
    std::string comment("!");
    while(line.compare(0,1,comment)==0) getline(input,line);
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
    int errorcount=0,runnumber=0;
    do{
        getline(input,line);
        if(line.length()<10) continue;
        runnumber++;
        int sel=0, leng=0, cube=0, mres=0, a1=0, a2=0, a3=0;
        double res[3],box[3];
        bool bleng=false,bres=false,bcube=false;
        if(line.length() > 1)
        {
            std::istringstream iss(line);
            if (!(iss >> sel >> leng >> mres >> res[0] >> res[1] >> res[2] >> cube >> box[0] >> box[1] >> box[2] >> a1 >> a2 >> a3)) {
                std::cerr << "Error parsing line for bond parameters." << std::endl;
                continue; // Skip this line and move to the next
            }
        }
        if(leng==1){ bleng=true; if(debug) std::cout << "leng=true\n";}
        else { bleng=false; if(debug) std::cout << "leng=false\n";}
        if(cube==1){ bcube=true; if(debug) std::cout << "cube=true\n";}
        else { bcube=false; if(debug) std::cout << "cube=false\n";}
        if(mres==1){ bres=true; if(debug) std::cout << "mres=true\n";}
        else { bres=false; if(debug) std::cout << "mres=false\n";}
        if (debug) std::cout << "running calculations for line " << runnumber<< " of the input file:\n\n";
        bond work{"","","","",false,false,false};
        work=do_bonds(wavy,sel,bleng,bres,res,bcube,box,a1,a2,a3,debug,bohr,runnumber,rho==1,rdg==1,eli==1,lap==1);
        if(!work.success)
        {
            std::cout << "!!!!!!!!!!!!problem somewhere during the calculations, see messages above!!!!!\n";
            errorcount++;
        }
        else {
            if(rho==1) wavy.push_back_cube(work.filename + "_rho.cube", false, false);
            if(rdg==1) wavy.push_back_cube(work.filename + "_rdg.cube", false, false);
            if(eli==1) wavy.push_back_cube(work.filename + "_eli.cube", false, false);
            if(lap==1) wavy.push_back_cube(work.filename + "_lap.cube", false, false);
            if(rho==1&&rdg==1) wavy.push_back_cube(work.filename + "_signed_rho.cube", false, false);
        }
    }while (!input.eof());
    std::cout << "\n  *   Finished all calculations! " << runnumber-errorcount << " out of " << runnumber << " were successful!            *\n";
    return 1;
}
