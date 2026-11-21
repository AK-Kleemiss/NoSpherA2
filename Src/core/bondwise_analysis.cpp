#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "bondwise_analysis.h"
#include "properties.h"
#include "topology.h"
#include <filesystem>
#include <fstream>

int compute_dens(WFN &wavy, bool debug, int *np, double *origin, double *gvector, double *incr, std::string &outname, bool rho, bool rdg, bool eli, bool lap) {
	properties_options opts;
	opts.lap = lap;
	opts.eli = eli;
	opts.rdg = rdg;
	opts.rho = rho;
	opts.resolution = *incr;
	opts.NbSteps = { np[0], np[1], np[2] };

	//indexed by cube_type: Rho, RDG, Elf, Eli, Lap, then the eleven properties this run never asks for
	std::vector<cube> cubes = { {opts.NbSteps, wavy.get_ncen(), opts.rho || opts.rdg || opts.eli || opts.lap},
		{opts.NbSteps, wavy.get_ncen(), opts.rdg},
		{},
		{opts.NbSteps, wavy.get_ncen(), opts.eli},
		{opts.NbSteps, wavy.get_ncen(), opts.lap},
		{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {} };
	cubes[cube_type::Rho].give_parent_wfn(wavy);
	cubes[cube_type::RDG].give_parent_wfn(wavy);
	cubes[cube_type::Eli].give_parent_wfn(wavy);
	cubes[cube_type::Lap].give_parent_wfn(wavy);

	std::string Oname = outname;
	ivec ntyp;
	for (int i = 0; i < 3; i++) {
		if (debug) {
			std::cout << "gvector before: ";
			for (int j = 0; j < 3; j++) std::cout << gvector[j + 3 * i] << " ";
			std::cout << "\n";
		}
		for (int j = 0; j < 3; j++)
			gvector[j + 3 * i] = incr[i] * gvector[j + 3 * i];
		if (debug) {
			for (int j = 0; j < 3; j++) std::cout << gvector[j + 3 * i] << " ";
			std::cout << "\n";
		}
	}
	for (int i = 0; i < 3; i++) {
		cubes[cube_type::Rho].set_origin(i, origin[i]);
		cubes[cube_type::RDG].set_origin(i, origin[i]);
		cubes[cube_type::Eli].set_origin(i, origin[i]);
		cubes[cube_type::Lap].set_origin(i, origin[i]);
		for (int j = 0; j < 3; j++) {
			cubes[cube_type::Rho].set_vector(i, j, gvector[j + 3 * i]);
			cubes[cube_type::RDG].set_vector(i, j, gvector[j + 3 * i]);
			cubes[cube_type::Eli].set_vector(i, j, gvector[j + 3 * i]);
			cubes[cube_type::Lap].set_vector(i, j, gvector[j + 3 * i]);
		}
	}
	cubes[cube_type::Rho].set_path(Oname + "_rho.cube");
	cubes[cube_type::RDG].set_path(Oname + "_rdg.cube");
	cubes[cube_type::Eli].set_path(Oname + "_eli.cube");
	cubes[cube_type::Lap].set_path(Oname + "_lap.cube");

	//opt.NbAtoms[0]=wavy.get_ncen();
	std::cout << "\n   .      ___________________________________________________________      .\n";
	std::cout << "  *.                                                                      *.\n";
	std::cout << "  *.  Wavefunction              : " << std::setw(20) << wavy.get_path().filename() << " / " << std::setw(5) << wavy.get_ncen() << " atoms      *.\n";
	std::cout << "  *.  OutPut filename Prefix    : " << std::setw(40) << Oname << "*.\n";
	std::cout << "  *.                                                                      *.\n";
	std::cout << "  *.  gridBox Min               : " << std::setw(11) << std::setprecision(6) << origin[0] << std::setw(12) << origin[1] << origin[2] << "     *.\n";
	std::cout << "  *.  gridBox Max               : " << std::setw(11) << std::setprecision(6) << (origin[0] + gvector[0] * np[0] + gvector[3] * np[1] + gvector[6] * np[2]) << std::setw(12) << (origin[1] + gvector[1] * np[0] + gvector[4] * np[1] + gvector[7] * np[2]) << (origin[2] + gvector[2] * np[0] + gvector[5] * np[1] + gvector[8] * np[2]) << "     *.\n";
	std::cout << "  *.  Increments(bohr)          : " << std::setw(11) << std::setprecision(6) << array_length(d3{ gvector[0], gvector[1], gvector[2] }) << std::setw(12) << array_length(d3{ gvector[3], gvector[4], gvector[5] }) << std::setw(12) << array_length(d3{ gvector[6], gvector[7], gvector[8] }) << "     *.\n";
	std::cout << "  *.  NbSteps                   : " << std::setw(11) << opts.NbSteps[0] << std::setw(12) << opts.NbSteps[1] << std::setw(12) << opts.NbSteps[2] << "     *.\n";
	std::cout << "  *.                                                                      *.\n";
	std::cout << "  *.  Number of primitives      :     " << std::setw(5) << wavy.get_nex() << "                               *.\n";
	std::cout << "  *.  Number of MOs             :       " << std::setw(3) << wavy.get_nmo() << "                               *.\n";

	//Calc_Prop leaves the Rho cube untouched unless RDG is requested (it hands back sign(lambda2)*rho then)
	if (rho && !rdg) Calc_Rho(cubes[cube_type::Rho], wavy, 20.0, std::cout, false);
	if (rdg || eli || lap) Calc_Prop(cubes, wavy, 20.0, std::cout, false, false);

	std::cout << "  *.                                                                      *.\n";
	std::cout << "  *.  Writing .cube files ...                                             *.\n";
	if (rho && !rdg) {
		cubes[cube_type::Rho].set_path(Oname + "_rho.cube");
		cubes[cube_type::Rho].write_file(true, true);
	}
	if (rdg) {
		cubes[cube_type::Rho].set_path(Oname + "_signed_rho.cube");
		cubes[cube_type::Rho].write_file(true);
		cubes[cube_type::Rho].set_path(Oname + "_rho.cube");
		cubes[cube_type::Rho].write_file(true, true);
		cubes[cube_type::RDG].write_file(true);
	}
	if (eli) cubes[cube_type::Eli].write_file(true);
	if (lap) cubes[cube_type::Lap].write_file(true);

	std::cout << "  *                                                                       *\n";
	std::cout << "          ___________________________________________________________\n";
	return 0;
};

bond do_bonds(WFN &wavy,
	int mode_sel, bool mode_leng, bool mode_res,
	double res[], bool cub, double boxsize[],
	int atom1, int atom2, int atom3,
	const bool &debug, const bool &bohr,
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
	//bohr: the user's lengths (box lengths when mode_leng, spacings when !mode_res) are Angstrom; everything below is bohr
	double res_b[3], boxsize_b[3];
	for (int r = 0; r < 3; r++)
	{
		res_b[r] = (bohr && !mode_res) ? constants::ang2bohr(res[r]) : res[r];
		boxsize_b[r] = (bohr && mode_leng) ? constants::ang2bohr(boxsize[r]) : boxsize[r];
	}
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
				size[0] = ceil(9 * znorm) / 10;
				break;
			case 2:
				size[0] = ceil(15 * znorm) / 10;
				break;
			case 3:
				size[0] = ceil(15 * hnorm) / 10;
				break;
			case 4:
				size[0] = ceil(30 * hnorm) / 10;
				break;
			}
			for (int r = 0; r < 3; r++)
			{
				if (boxsize_b[r] != 0) size[r] = boxsize_b[r];
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
					size[r] = size[0];
					continue;
				}
				switch (mode_sel)
				{
				case 1:
				case 2:
					size[r] = znorm;
					break;
				case 3:
				case 4:
					size[r] = hnorm;
					break;
				}
				if (boxsize[r] != 0) size[r] = size[r] * boxsize[r];
				else {
					switch (mode_sel)
					{
					case 1:
						size[r] = ceil(9 * znorm) / 10;
						break;
					case 2:
						size[r] = ceil(15 * znorm) / 10;
						break;
					case 3:
						size[r] = ceil(15 * hnorm) / 10;
						break;
					case 4:
						size[r] = ceil(30 * hnorm) / 10;
						break;
					}
				}
				s2[r] = size[r] / 2;
			}
		}
		for (int r = 0; r < 3; r++) np[r] = (int)round((2 * s2[r]) / res_b[cub ? 0 : r]) + 1;
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
				size[0] = ceil(9 * znorm) / 10;
				break;
			case 2:
				size[0] = ceil(15 * znorm) / 10;
				break;
			case 3:
				size[0] = ceil(15 * hnorm) / 10;
				break;
			case 4:
				size[0] = ceil(30 * hnorm) / 10;
				break;
			}
			for (int r = 0; r < 3; r++)
			{
				if (debug) std::cout << r + 1 << ". Axis:";
				if (boxsize_b[r] != 0) size[r] = boxsize_b[r];
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
					size[r] = znorm;
					break;
				case 3:
				case 4:
					size[r] = hnorm;
					break;
				}
				if (boxsize[r] != 0) size[r] = size[r] * boxsize[r];
				else {
					switch (mode_sel)
					{
					case 1:
						size[r] = ceil(9 * znorm) / 10;
						break;
					case 2:
						size[r] = ceil(15 * znorm) / 10;
						break;
					case 3:
						size[r] = ceil(15 * hnorm) / 10;
						break;
					case 4:
						size[r] = ceil(30 * hnorm) / 10;
						break;
					}
				}
				s2[r] = size[r] / 2;
			}
		}
	}
	//grid spacing: the requested one when res is a spacing, box length / points otherwise
	double incr[3];
	for (int i = 0; i < 3; i++)
		incr[i] = mode_res ? size[i] / np[i] : res_b[cub ? 0 : i];
	//origin = centre - half the grid extent (np - 1 steps per axis), so the grid points are symmetric about the centre
	double o[3] = {};
	if (debug) std::cout << "This is the origin:\n";
	for (int d = 0; d < 3; d++)
	{
		switch (mode_sel)
		{
		case 1:
			o[d] = coords1[d];
			break;
		case 2:
			o[d] = (coords1[d] + coords2[d]) / 2;
			break;
		case 3:
			o[d] = (coords3[d] + coords1[d]) / 2;
			break;
		case 4:
			o[d] = ringhelp[d];
			break;
		}
		o[d] -= 0.5 * ((np[0] - 1) * incr[0] * x[d] + (np[1] - 1) * incr[1] * y[d] + (np[2] - 1) * incr[2] * z[d]);
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
	for (int i = 0; i < 3; i++) {
		v[i] = x[i];
		v[i + 3] = y[i];
		v[i + 6] = z[i];
	}
	if (compute_dens(wavy, debug, np, o, v, incr, outname, rho, rdg, eli, lap) == 0) {
		results.success = true;
		results.filename = outname;
	}
	return(results);
}

int autobonds(bool debug, WFN &wavy, const std::filesystem::path &inputfile, const bool &bohr) {
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
	getline_universal(input, line);
	std::string comment("!");
	while (line.compare(0, 1, comment) == 0) getline_universal(input, line);
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
		getline_universal(input, line);
		if (line.length() < 10) continue;
		runnumber++;
		int sel = 0, leng = 0, cube = 0, mres = 0, a1 = 0, a2 = 0, a3 = 0;
		double res[3], box[3];
		bool bleng = false, bres = false, bcube = false;
		if (line.length() > 1)
		{
			std::istringstream iss(line);
			if (!(iss >> sel >> leng >> mres >> res[0] >> res[1] >> res[2] >> cube >> box[0] >> box[1] >> box[2] >> a1 >> a2 >> a3)) {
				std::cerr << "Error parsing line for bond parameters.\n";
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


void write_bond_path_profiles(const WFN &wavy, const std::vector<topology::bond_path> &paths, const std::string &stem, std::ostream *log)
{
	std::map<std::pair<int, int>, int> seen; //a second path between one pair must not overwrite the first
	for (const topology::bond_path &bp : paths) {
		if (bp.a < 0) continue;
		std::string name = stem + "_bond_path_" + std::to_string(bp.a) + "_" + std::to_string(bp.b);
		if (const int n = seen[{ bp.a, bp.b }]++) name += "_" + std::to_string(n + 1);
		if (log) {
			*log << "Bond path " << bp.a << " (" << wavy.get_atom_charge(bp.a) << ") - " << bp.b << " (" << wavy.get_atom_charge(bp.b)
				<< "): BPL " << std::fixed << std::setprecision(4) << bp.bpl << " GPL " << bp.gpl << " BPL-GPL " << bp.bpl - bp.gpl << " bohr";
			if (bp.gap_a > 0.0 || bp.gap_b > 0.0)
				*log << ", end maxima " << bp.gap_a << " / " << bp.gap_b << " bohr off the nuclei";
			*log << "\n";
		}
		std::ofstream result(std::filesystem::current_path() / std::filesystem::path(name + ".dat"), std::ios::out);
		result << "# s x y z rho laplacian, bohr and a.u.; the bond critical point is point " << bp.bcp_index + 1 << "\n";
		double s = 0.0;
		for (size_t k = 0; k < bp.points.size(); k++) {
			if (k) s += array_length(bp.points[k], bp.points[k - 1]);
			d3 grad{ 0.0, 0.0, 0.0 };
			double lap = 0.0;
			const double rho = calculate_density(wavy, bp.points[k], grad, lap);
			result << std::scientific << std::setprecision(6) << std::setw(14) << s;
			for (int c = 0; c < 3; c++) result << std::setw(14) << bp.points[k][c];
			result << std::setw(14) << rho << std::setw(14) << lap << "\n";
		}
	}
}

void bondwise_laplacian_plots(std::filesystem::path &wfn_name, const bool along_path)
{
	WFN wavy(wfn_name);
	wavy.delete_unoccupied_MOs();

	err_checkf(wavy.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);

	if (along_path) {
		topology::options topo_opt;
		topo_opt.bond_scale = 1.35; //the pairs the straight-line mode plots
		const std::vector<topology::nucleus> nuclei = topology::nuclei_of(wavy);
		const topology::result top = topology::analyze_topology(wavy, nuclei, topo_opt);
		std::vector<d3> bcps;
		for (const topology::cp &p : top.points)
			if (p.kind == topology::cp_kind::bond) bcps.push_back(p.position);
		std::vector<topology::bond_path> paths(bcps.size());
#pragma omp parallel for schedule(dynamic)
		for (int k = 0; k < (int)bcps.size(); k++)
			paths[k] = topology::trace_bond_path(wavy, bcps[k], nuclei);
		write_bond_path_profiles(wavy, paths, wfn_name.string(), &std::cout);
		return;
	}

	int points = 1001;

	for (int i = 0; i < wavy.get_ncen(); i++)
	{
		for (int j = i + 1; j < wavy.get_ncen(); j++)
		{
			std::filesystem::path path = std::filesystem::current_path();
			double distance = sqrt(pow(wavy.get_atom_coordinate(i, 0) - wavy.get_atom_coordinate(j, 0), 2) + pow(wavy.get_atom_coordinate(i, 1) - wavy.get_atom_coordinate(j, 1), 2) + pow(wavy.get_atom_coordinate(i, 2) - wavy.get_atom_coordinate(j, 2), 2));
			double svdW = constants::ang2bohr(constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)]);
			if (distance < 1.35 * svdW)
			{
				std::cout << "Bond between " << i << " (" << wavy.get_atom_charge(i) << ") and " << j << " (" << wavy.get_atom_charge(j) << ") with distance " << distance << " and svdW " << svdW << "\n";
				const vec bond_vec = { (wavy.get_atom_coordinate(j, 0) - wavy.get_atom_coordinate(i, 0)) / points, (wavy.get_atom_coordinate(j, 1) - wavy.get_atom_coordinate(i, 1)) / points, (wavy.get_atom_coordinate(j, 2) - wavy.get_atom_coordinate(i, 2)) / points };
				const double dr = distance / points;
				vec lapl(points, 0.0);
				const vec pos = { wavy.get_atom_coordinate(i, 0), wavy.get_atom_coordinate(i, 1), wavy.get_atom_coordinate(i, 2) };
#pragma omp parallel for schedule(dynamic)
				for (int k = 0; k < points; k++)
				{
					d3 t_pos = { pos[0], pos[1], pos[2] };
					t_pos[0] += k * bond_vec[0];
					t_pos[1] += k * bond_vec[1];
					t_pos[2] += k * bond_vec[2];
					lapl[k] = wavy.computeLap(t_pos);
				}
				std::filesystem::path outname(wfn_name.string() + "_bondwise_laplacian_" + std::to_string(i) + "_" + std::to_string(j) + ".dat");
				path = path / std::filesystem::path(outname);
				std::ofstream result(path, std::ios::out);
				for (int k = 0; k < points; k++)
				{
					result << std::setw(10) << std::scientific << std::setprecision(6) << dr * k << " " << std::setw(10) << std::scientific << std::setprecision(6) << lapl[k] << "\n";
				}
				result.close();
			}
			else
			{
				std::cout << "No bond between " << i << " and " << j << " with distance " << distance << " and svdW " << svdW << "\n";
			}
		}
	}
}
