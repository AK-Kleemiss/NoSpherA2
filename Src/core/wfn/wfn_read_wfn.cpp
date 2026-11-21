#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "mo_class.h"
#include "cube.h"
#include "constants.h"
#include "fchk.h"
#include "basis_set.h"
#include "nos_math.h"
#include "libCintMain.h"
#include "integrator.h"
#include "cell.h"

bool WFN::read_wfn(const std::filesystem::path &fileName, const bool &debug, std::ostream &file)
{
	using namespace std;
	if (ncen > 0)
	{
		file << "There is already a wavefunction loaded, aborting!" << endl;
		return false;
	}
	origin = e_origin::wfn;
	d_f_switch = true;
	err_checkf(std::filesystem::exists(fileName), "Couldn't open or find " + fileName.string() + ", leaving", file);
	ifstream rf(fileName);
	if (rf.good())
		path = fileName;
	string line;
	read_line_or_fail(rf, line, "wfn title line", file);
	comment = line;
	read_line_or_fail(rf, line, "wfn header line", file);
	stringstream stream(line);
	string header_tmp;
	int e_nmo, e_nex, e_nuc = 0; // number of expected MOs, Exponents and nuclei
	stream >> header_tmp >> e_nmo >> header_tmp >> header_tmp >> e_nex >> header_tmp >> e_nuc;
	err_checkf(!stream.fail() && e_nmo > 0 && e_nex > 0 && e_nuc > 0, "Bad wfn header line: '" + line + "'", file);
	if (debug)
		file << "e_nmo: " << e_nmo << ", e_nex: " << e_nex << ", e_nuc : " << e_nuc << endl;
	//The wfn blocks are fixed-width Fortran columns (20I3 assignments, 5E14 exponents, 5E16 MO
	//coefficients); fields can touch, so they are cut by width, not by whitespace. Fortran 'D'
	//exponents are accepted. Reads lines starting with tag until n values are collected
	auto read_block = [&](const string &tag, const size_t start, const size_t width, const int n, auto &out)
	{
		const string what = tag.empty() ? "MO coefficients" : tag;
		out.clear();
		while (static_cast<int>(out.size()) < n)
		{
			read_line_or_fail(rf, line, what, file);
			err_checkf(line.compare(0, tag.size(), tag) == 0 && line.size() > start, "Expected a " + what + " line but found: '" + line + "'", file);
			for (size_t pos = start; pos < line.size() && static_cast<int>(out.size()) < n; pos += width)
			{
				string field = line.substr(pos, width);
				replace(field.begin(), field.end(), 'D', 'E');
				append_numbers(field, out, what, file);
			}
		}
	};
	//A fixed-width field of line as a number
	auto field = [&](const size_t pos, const size_t len)
	{
		err_checkf(line.size() > pos, "wfn line too short: '" + line + "'", file);
		return stod(line.substr(pos, len));
	};
	//----------------------------- Read Atoms ------------------------------------------------------------
	for (int i = 0; i < e_nuc; i++)
	{
		read_line_or_fail(rf, line, "atom " + to_string(i + 1) + " of " + to_string(e_nuc), file);
		err_checkf(line.size() >= 73, "wfn atom line " + to_string(i + 1) + " is too short: '" + line + "'", file);
		//Columns are fixed only relative to the label, whose width varies between writers, and the 12-character coordinate
		//fields touch when one is negative, so they are cut by width, anchored two characters past the ')' of "(CENTRE n)".
		//The charge follows the last '=', the label is the first word.
		const size_t centre_end = line.find(')');
		err_checkf(centre_end != string::npos && centre_end + 2 + 36 <= line.size(), "wfn atom line " + to_string(i + 1) + " carries no '(CENTRE n)' with 3 coordinates after it: '" + line + "'", file);
		const size_t chg = line.find("CHARGE");
		err_checkf(chg != string::npos && chg >= centre_end + 2 + 36, "wfn atom line " + to_string(i + 1) + " carries no CHARGE field after its coordinates: '" + line + "'", file);
		const size_t eq = line.find('=', chg);
		err_checkf(eq != string::npos && eq + 1 < line.size(), "wfn atom line " + to_string(i + 1) + " carries no '= <charge>': '" + line + "'", file);
		const int charge = static_cast<int>(stod(line.substr(eq + 1)));
		string label;
		{
			std::istringstream first(line);
			first >> label;
		}
		err_checkf(!label.empty(), "wfn atom line " + to_string(i + 1) + " carries no label: '" + line + "'", file);
		err_checkf(push_back_atom(shrink_string_to_atom(label, charge), field(centre_end + 2, 12), field(centre_end + 14, 12), field(centre_end + 26, 12), charge), "Error while making atoms!!\n", file);
	}
	//------------------------------ Read basis: centres, types, exponents ------------------------------
	ivec centre, type;
	vec exponent;
	read_block("CENTRE ASSIGNMENTS", 20, 3, e_nex, centre);
	read_block("TYPE ASSIGNMENTS", 20, 3, e_nex, type);
	read_block("EXPONENTS", 10, 14, e_nex, exponent);
	for (int j = 0; j < e_nex; j++)
	{
		err_checkf(centre[j] >= 1 && centre[j] <= e_nuc, "Primitive " + to_string(j + 1) + " sits on centre " + to_string(centre[j]) + " of " + to_string(e_nuc), file);
		err_checkf(add_exp(centre[j], type[j], exponent[j]), "Error while adding primitive " + to_string(j + 1), file);
	}
	isBohr = true;
	//-------------------------------- Read MOs --------------------------------------
	//A .wfn carries no spin labels, so a second spin set shows only as energies starting over. The energy is read after
	//the last '=' (a fixed field truncates digits), only a strict decrease starts a new set because degenerate orbitals
	//print equal energies, and an occupation above 1 rules out a beta set.
	int oper = 0;
	double last_ener = -DBL_MAX;
	bool spin_orbitals = true;
	vec coef;
	for (int monum = 0; monum < e_nmo; monum++)
	{
		read_line_or_fail(rf, line, "MO " + to_string(monum + 1) + " of " + to_string(e_nmo), file);
		err_checkf(line.compare(0, 2, "MO") == 0, "Expected MO " + to_string(monum + 1) + " but found: '" + line + "'", file);
		const int nr = static_cast<int>(field(2, 6));
		const double occ = field(36, 12);
		if (occ > 1.0 + 1e-6)
			spin_orbitals = false;
		//whatever follows the last '=', which also covers the writers that glue it to the number
		const size_t eq = line.rfind('=');
		err_checkf(eq != string::npos && eq + 1 < line.size(), "MO " + to_string(monum + 1) + " carries no '= <energy>': '" + line + "'", file);
		const double ener = stod(line.substr(eq + 1));
		if (spin_orbitals && ener < last_ener)
		{
			oper++;
			is_unrestricted = true;
		}
		last_ener = ener;
		push_back_MO(nr, occ, ener, oper);
		read_block("", 0, 16, e_nex, coef);
		for (const double c : coef)
			MOs.back().push_back_coef(c);
	}
	//ponytail: what follows (END DATA, energy, virial) is not read; PySCF unrestricted files list a
	//second spin set past the header count and the old reader ignored it too
	set_exp_cutoff();
	return true;
};

bool WFN::read_xyz(const std::filesystem::path &filename, std::ostream &file, const bool debug)
{
	using namespace std;
	err_checkf(filesystem::exists(filename), "Couldn't open or find " + filename.string() + ", leaving", file);
	origin = e_origin::xyz;
	ifstream rf(filename.c_str());
	if (rf.good())
		path = filename;
	string line;
	rf.seekg(0);

	// Hand-edited xyz files arrive with the count line missing, an atom too many or
	// too few, or a label that is not an element; every one of those used to end in
	// a bare "invalid stod argument" or a charge-0 atom. Name the file and the line.
	const string where = filename.string() + ": ";
	auto strict_double = [&](const string &token, const string &what) {
		size_t used = 0;
		double value = 0.0;
		try { value = stod(token, &used); } catch (...) { used = 0; }
		err_checkf(used == token.size(), where + what + " '" + token + "' is not a number", file);
		return value;
	};
	auto tokens_of = [](const string &l) {
		svec t = split_string<string>(l, " ");
		remove_empty_elements(t);
		return t;
	};

	getline_universal(rf, line);
	svec head = tokens_of(line);
	err_checkf(head.size() == 1, where + "first line must be the atom count alone, found '" + line + "'", file);
	const double count = strict_double(head[0], "atom count");
	err_checkf(count > 0 && count == floor(count), where + "atom count '" + head[0] + "' must be a positive integer", file);
	int e_nuc = static_cast<int>(count); // number of expected nuclei
	if (debug)
		file << "e_nuc: " << e_nuc << endl;
	getline_universal(rf, line);
	comment = line;
	//----------------------------- Read Atoms ------------------------------------------------------------
	ivec dum_nr, dum_ch;
	dum_nr.resize(e_nuc);
	dum_ch.resize(e_nuc);
	svec dum_label;
	vec dum_x, dum_y, dum_z;
	dum_x.resize(e_nuc);
	dum_y.resize(e_nuc);
	dum_z.resize(e_nuc);
	dum_label.resize(e_nuc);
	for (int i = 0; i < e_nuc; i++)
	{
		const string atom_i = "atom " + to_string(i + 1) + " of " + to_string(e_nuc);
		err_checkf(static_cast<bool>(getline_universal(rf, line)), where + "file ends before " + atom_i + " (atom count too large?)", file);
		if (debug)
			file << i << ".run, line:" << line << "\n";
		dum_nr[i] = i;
		svec temp = tokens_of(line);
		err_checkf(temp.size() >= 4, where + atom_i + " needs 'El x y z', found '" + line + "'", file);
		dum_label[i] = temp[0];
		dum_x[i] = constants::ang2bohr(strict_double(temp[1], atom_i + " x"));
		dum_y[i] = constants::ang2bohr(strict_double(temp[2], atom_i + " y"));
		dum_z[i] = constants::ang2bohr(strict_double(temp[3], atom_i + " z"));
		const int Z = constants::get_Z_from_label(dum_label[i].c_str());
		err_checkf(Z >= 0, where + atom_i + " has unknown element label '" + dum_label[i] + "'", file);
		dum_ch[i] = Z + 1;
		if (debug)
		{
			file << "label:" << dum_label[i]
				<< " nr: " << dum_nr[i]
				<< " x: " << dum_x[i]
				<< " y: " << dum_y[i]
				<< " z: " << dum_z[i]
				<< " charge: " << dum_ch[i] << "\n";
		}
	}
	// Anything left must be a further frame (its own count line) or blank; a
	// trailing atom line means the count is too small.
	while (getline_universal(rf, line))
	{
		svec rest = tokens_of(line);
		if (rest.empty())
			continue;
		err_checkf(rest.size() == 1, where + "more atom lines than the atom count " + to_string(e_nuc) + " announces, first extra line '" + line + "'", file);
		break;
	}
	isBohr = true;
	//---------------------Start writing everything from the temp arrays into wave ---------------------
	if (debug)
		file << "finished with reading the file, now i'm going to make everything permantent in the wavefunction...\n";

	for (int i = 0; i < e_nuc; i++)
		err_checkf(push_back_atom(dum_label[i], dum_x[i], dum_y[i], dum_z[i], dum_ch[i]), "Error while making atoms!!", file);
	return true;
};

bool WFN::read_wfx(const std::filesystem::path &fileName, const bool &debug, std::ostream &file)
{
	origin = e_origin::wfx;
	d_f_switch = true;
	using namespace std;
	err_checkf(std::filesystem::exists(fileName), "Couldn't open or find " + fileName.string() + ", leaving", file);
	ifstream rf(fileName.c_str());
	path = fileName;
	string line;
	//Every block is looked up from the top, so the order of the blocks in the file does not matter
	auto seek = [&](const string &tag) {
		rf.clear();
		rf.seekg(0);
		line.clear();
		seek_line(rf, line, "<" + tag + ">", file);
	};
	auto read_block = [&](const string &tag, auto &out) {
		seek(tag);
		const string end = "</" + tag + ">";
		size_t pos;
		do
		{
			read_line_or_fail(rf, line, end, file);
			pos = line.find(end);
			append_numbers(line.substr(0, pos), out, tag, file);
		} while (pos == string::npos);
	};
	//Read as double: some writers put the net charge as "0.0"
	auto read_int = [&](const string &tag) {
		vec v;
		read_block(tag, v);
		err_checkf(v.size() == 1, "Expected one number in <" + tag + ">, found " + to_string(v.size()), file);
		return static_cast<int>(v[0]);
	};
	seek("Title");
	read_line_or_fail(rf, line, "the title", file);
	comment = line;
	if (debug)
		file << "comment line " << line << endl;
	const int temp_ncen = read_int("Number of Nuclei");
	const int temp_nex = read_int("Number of Primitives");
	const int temp_nmo = read_int("Number of Occupied Molecular Orbitals");
	err_checkf(temp_ncen > 0 && temp_nex > 0 && temp_nmo > 0, "wfx announces " + to_string(temp_ncen) + " nuclei, " + to_string(temp_nex) + " primitives and " + to_string(temp_nmo) + " MOs", file);
	ivec nrs;
	read_block("Atomic Numbers", nrs);
	err_checkf(nrs.size() == temp_ncen, "Mismatch in atom number size", file);
	vec pos;
	read_block("Nuclear Cartesian Coordinates", pos);
	err_checkf(pos.size() == 3 * temp_ncen, "Mismatch in atom position size", file);
	for (int i = 0; i < temp_ncen; i++)
	{
		const string label = constants::atnr2letter(nrs[i]);
		err_checkf(label != "PROBLEM", "Unknown atomic number " + to_string(nrs[i]) + " for atom " + to_string(i + 1), file);
		push_back_atom(label + to_string(i + 1), pos[3 * i], pos[3 * i + 1], pos[3 * i + 2], nrs[i]);
	}
	err_checkf(ncen == temp_ncen, "Mismatch in atom position size", file);
	isBohr = true;
	charge = read_int("Net Charge");
	multi = read_int("Electronic Spin Multiplicity");
	read_block("Primitive Centers", centers);
	read_block("Primitive Types", types);
	read_block("Primitive Exponents", exponents);
	err_checkf(exponents.size() == temp_nex && centers.size() == temp_nex && types.size() == temp_nex, "Mismatch in numbers! aborting!", file);
	for (int i = 0; i < temp_nex; i++)
		err_checkf(centers[i] >= 1 && centers[i] <= temp_ncen, "Primitive " + to_string(i + 1) + " sits on centre " + to_string(centers[i]) + " of " + to_string(temp_ncen), file);
	nex = temp_nex;
	vec occ, ener;
	read_block("Molecular Orbital Occupation Numbers", occ);
	read_block("Molecular Orbital Energies", ener);
	err_checkf(occ.size() == temp_nmo && ener.size() == temp_nmo, "Found " + to_string(occ.size()) + " occupations and " + to_string(ener.size()) + " energies for " + to_string(temp_nmo) + " MOs", file);
	//A wfx labels every orbital's spin; the energy heuristic below is only for files without that block.
	std::vector<std::string> spins;
	{
		rf.clear();
		rf.seekg(0);
		string l;
		bool inside = false;
		while (getline(rf, l))
		{
			if (!inside)
			{
				inside = l.find("<Molecular Orbital Spin Types>") != string::npos;
				continue;
			}
			if (l.find("</Molecular Orbital Spin Types>") != string::npos)
				break;
			const size_t a = l.find_first_not_of(" \t\r\n");
			if (a == string::npos)
				continue;
			spins.push_back(l.substr(a, l.find_last_not_of(" \t\r\n") - a + 1));
		}
	}
	const bool labelled = spins.size() == static_cast<size_t>(temp_nmo);
	err_checkf(spins.empty() || labelled, "<Molecular Orbital Spin Types> lists " + to_string(spins.size()) + " labels for " + to_string(temp_nmo) + " orbitals", file);
	double last_ener = -DBL_MAX;
	int oper = 0;
	bool spin_orbitals = true;
	for (int i = 0; i < temp_nmo; i++)
	{
		if (labelled)
		{
			//"Alpha and Beta" is one spin-averaged orbital, not two, and belongs to the first operator
			err_checkf(spins[i] == "Alpha" || spins[i] == "Beta" || spins[i] == "Alpha and Beta",
					   "Orbital " + to_string(i + 1) + " carries the spin label '" + spins[i] + "', which is none of Alpha, Beta or 'Alpha and Beta'", file);
			oper = spins[i] == "Beta" ? 1 : 0;
			if (oper == 1)
				is_unrestricted = true;
		}
		else
		{
			if (occ[i] > 1.0 + 1e-6)
				spin_orbitals = false;
			if (spin_orbitals && ener[i] < last_ener)
			{
				oper++;
				is_unrestricted = true;
			}
			last_ener = ener[i];
		}
		err_checkf(push_back_MO(i + 1, occ[i], ener[i], oper), "Error poshing back MO! MO: " + to_string(i), file);
	}
	seek("Molecular Orbital Primitive Coefficients");
	const string end = "</Molecular Orbital Primitive Coefficients>";
	vec coef;
	int read_MOs = 0;
	while (read_line_or_fail(rf, line, end, file), line.find(end) == string::npos)
	{
		if (line.find("<MO Number>") == string::npos)
			continue;
		ivec nr;
		read_line_or_fail(rf, line, "MO Number", file);
		append_numbers(line, nr, "MO Number", file);
		err_checkf(nr.size() == 1 && nr[0] >= 1 && nr[0] <= temp_nmo, "Bad <MO Number> line: '" + line + "'", file);
		seek_line(rf, line, "</MO Number>", file);
		const string what = "coefficients of MO " + to_string(nr[0]);
		coef.clear();
		while (coef.size() < nex)
		{
			read_line_or_fail(rf, line, what, file);
			append_numbers(line, coef, what, file);
		}
		err_checkf(coef.size() == nex, "More coefficients than primitives in " + what, file);
		for (int i = 0; i < nex; i++)
			MOs[nr[0] - 1].push_back_coef(coef[i]);
		read_MOs++;
	}
	err_checkf(read_MOs == temp_nmo, "Read coefficients for " + to_string(read_MOs) + " of " + to_string(temp_nmo) + " MOs", file);
	seek_line(rf, line, "<Energy =", file);
	read_line_or_fail(rf, line, "the energy", file);
	total_energy = stod(line);
	seek_line(rf, line, "<Virial Ratio", file);
	read_line_or_fail(rf, line, "the virial ratio", file);
	virial_ratio = stod(line);
	rf.close();
	set_exp_cutoff();
	return true;
};
