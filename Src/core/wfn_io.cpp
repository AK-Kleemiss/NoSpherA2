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


const std::string WFN::hdr(const bool &occupied) const
{
	std::stringstream temp;
	temp << "GTO";
	if (!occupied)
	{
		temp << std::setw(20) << nmo;
	}
	else
	{
		const int occupied_mos = get_nmo(true);
		temp << std::setw(20) << occupied_mos;
	}
	temp << " MOL ORBITALS";
	temp << std::setw(7) << nex;
	temp << " PRIMITIVES";
	temp << std::setw(9) << ncen;
	temp << " NUCLEI\n";
	return temp.str();
};

void WFN::read_known_wavefunction_format(const std::filesystem::path &fileName, std::ostream &file, const bool debug)
{
	if (fileName.extension() == ".wfn")
		err_checkf(read_wfn(fileName, debug, file), "Problem reading wfn", file);
	else if (fileName.extension() == ".ffn")
		err_checkf(read_wfn(fileName, debug, file), "Problem reading ffn", file);
	else if (fileName.extension() == ".wfx")
		err_checkf(read_wfx(fileName, debug, file), "Problem reading wfx", file);
	else if (fileName.extension() == ".fch" || fileName.extension() == ".fchk" || fileName.extension() == ".FCh" || fileName.extension() == ".FChK" || fileName.extension() == ".FChk")
		err_checkf(read_fchk(fileName, file, debug), "Problem reading fchk", file);
	else if (fileName.extension() == ".xyz")
		err_checkf(read_xyz(fileName, file, debug), "Problem reading xyz", file);
	else if (fileName.extension() == ".molden")
		err_checkf(read_molden(fileName, file, debug), "Problem reading molden file", file);
	else if (fileName.extension() == ".gbw")
		err_checkf(read_gbw(fileName, file, debug), "Problem reading gbw file", file);
	else if (fileName.extension() == ".xtb")
		err_checkf(read_ptb(fileName, file, debug), "Problem reading xtb file", file);
	else if (fileName.extension() == ".stda")
		err_checkf(read_ptb(fileName, file, debug), "Problem reading ptb file", file);
	else if (fileName.extension() == ".orbital_energies,restricted" || fileName.extension() == ".MO_energies,r"
		|| fileName.extension() == ".molecular_orbitals,restricted" || fileName.extension() == ".MOs,r"
		|| fileName.string().find("stdout") != std::string::npos)
		err_checkf(read_tonto(fileName, file, debug), "Problem reading tonto file", file);
	else
		err_checkf(false, "Unknown filetype!", file);
};

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
		string label = line.substr(0, 4);
		const int charge = static_cast<int>(field(70, 3));
		err_checkf(push_back_atom(shrink_string_to_atom(label, charge), field(24, 12), field(36, 12), field(48, 12), charge), "Error while making atoms!!\n", file);
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
	int oper = 0;
	double last_ener = -DBL_MAX;
	vec coef;
	for (int monum = 0; monum < e_nmo; monum++)
	{
		read_line_or_fail(rf, line, "MO " + to_string(monum + 1) + " of " + to_string(e_nmo), file);
		err_checkf(line.compare(0, 2, "MO") == 0, "Expected MO " + to_string(monum + 1) + " but found: '" + line + "'", file);
		const int nr = static_cast<int>(field(2, 6));
		const double occ = field(36, 12);
		// some writers glue the '=' to the energy
		if (line.size() > 61 && line[61] == '=') line[61] = ' ';
		const double ener = field(61, 12);
		//the energies restart from the bottom for the second spin
		if (ener > last_ener)
			last_ener = ener;
		else
		{
			last_ener = -DBL_MAX;
			oper++;
			is_unrestricted = true;
		}
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
	double last_ener = -DBL_MAX;
	int oper = 0;
	for (int i = 0; i < temp_nmo; i++)
	{
		if (ener[i] > last_ener)
		{
			last_ener = ener[i];
		}
		else
		{
			last_ener = -DBL_MAX;
			oper++;
			is_unrestricted = true;
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

//shell[m][s] holds the m-th spherical coefficient times the s-th contraction coefficient of one shell starting at prims[start]
static double odd_ft(const int n) { return n < 1 ? 1.0 : (double)constants::double_ft[n]; }
//norm of x^a y^b z^c relative to x^l: sqrt((2l-1)!! / ((2a-1)!!(2b-1)!!(2c-1)!!)), Gaussian, molden and fchk give only the x^l norm
static vec cart_norm(const int l)
{
	vec n(constants::n_cart(l));
	int v[3];
	for (int cart = 0; cart < (int)n.size(); cart++)
	{
		constants::type2vector(constants::first_type[l] + cart, v);
		n[cart] = sqrt(odd_ft(2 * l - 1) / (odd_ft(2 * v[0] - 1) * odd_ft(2 * v[1] - 1) * odd_ft(2 * v[2] - 1)));
	}
	return n;
}
//source column of each WFN Cartesian type; Gaussian/molden f: xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz, tonto f swaps types 16/17
static const int gaussian_f_order[10] = { 0, 1, 2, 4, 5, 8, 3, 6, 7, 9 };
static const int molden_g_order[15] = { 2, 8, 11, 6, 1, 7, 14, 13, 5, 10, 12, 9, 4, 3, 0 };
static const int tonto_f_order[10] = { 0, 1, 2, 3, 4, 6, 5, 7, 8, 9 };
static const int* const molden_order[5] = { nullptr, nullptr, nullptr, gaussian_f_order, molden_g_order };
static const double tonto_scale[5] = { 1.0, 1.0, 1.0 / sqrt(1.5), 1.0 / sqrt(5.0), 1.0 / sqrt(13.125) };
//shell[cart][s] in the source order; order[cart] picks the source column of WFN type first_type[l] + cart, scale[cart] its normalisation
void WFN::push_back_cartesian_shell(const int mo, const int l, const vec2& shell, const std::vector<primitive>& prims, const int start, const int size, const int* order, const double* scale)
{
	for (int s = 0; s < size; s++)
		for (int cart = 0; cart < constants::n_cart(l); cart++)
		{
			const double t = shell[order ? order[cart] : cart][s] * (scale ? scale[cart] : 1.0);
			push_back_MO_coef(mo, abs(t) < 1E-10 ? 0 : t);
			if (mo == 0)
			{
				push_back_exponent(prims[start + s].get_exp());
				push_back_center(prims[start].get_center());
				push_back_type(constants::first_type[l] + cart);
				nex++;
			}
		}
}
void WFN::push_back_spherical_shell(const int mo, const int l, const vec2& shell, const std::vector<primitive>& prims, const int start, const int size)
{
	const int nsph = constants::n_spher(l);
	vec2 c(constants::n_cart(l), vec(size, 0.0));
	for (int cart = 0; cart < (int)c.size(); cart++)
		for (int m = 0; m < nsph; m++)
			for (int s = 0; s < size; s++)
				c[cart][s] += constants::sph2cart(l)[cart * nsph + m] * shell[m][s];
	push_back_cartesian_shell(mo, l, c, prims, start, size);
}
//the neglected tail of c x^a y^b z^c exp(-ar^2) is bounded by c (u/a)^(l/2) exp(-u) for u = a r^2 beyond l/2,
//so the cutoff on -a r^2 alone loses 1e-3 electrons per l = 10 orbital; three fixed-point steps per primitive, the minimum wins
void WFN::set_exp_cutoff() const {
	const double cut0 = std::log(constants::density_accuracy / get_maximum_MO_coefficient());
	double cut = cut0;
	for (int i = 0; i < nex; i++) {
		int v[3];
		constants::type2vector(get_type(i), v);
		const int l = v[0] + v[1] + v[2];
		double ci = cut0;
		for (int it = 0; it < 3 && l > 0; it++)
			ci = cut0 - 0.5 * l * std::log(std::max(-ci / get_exponent(i), 0.5 * l));
		cut = std::min(cut, ci);
	}
	constants::exp_cutoff = cut;
}

bool WFN::read_molden(const std::filesystem::path &filename, std::ostream &file, const bool debug)
{
	using namespace std;
	err_checkf(std::filesystem::exists(filename), "couldn't open or find " + filename.string() + ", leaving", file);
	if (debug)
		file << "File is valid, continuing...\n"
		<< GetCurrentDir << endl;
	origin = e_origin::molden;
	isBohr = true;
	ifstream rf(filename.c_str());
	if (rf.good())
		path = filename;
	string line;
	read_line_or_fail(rf, line, "the first line", file);
	err_checkf(line.find("Molden Format") != string::npos, "Does not look like proper molden format file!", file);
	//The whitespace-separated fields of the current line
	auto fields = [&]()
	{
		svec f = split_string<string>(line, " ");
		remove_empty_elements(f);
		return f;
	};
	//Everything up to [Atoms] that is not a section header is the title
	while (read_line_or_fail(rf, line, "[Atoms]", file), line.find("[Atoms]") == string::npos)
		if (line.find("[") == string::npos)
			comment += trim(line);
	const bool angstrom = line.find("ngs") != string::npos;
	//----------------------------- Atoms: label index charge x y z ------------------------------
	while (read_line_or_fail(rf, line, "the atoms", file), line.find("[") == string::npos)
	{
		const svec f = fields();
		if (f.empty())
			continue;
		err_checkf(f.size() >= 6, "Molden atom line with fewer than 6 fields: '" + line + "'", file);
		const double scale = angstrom ? constants::ang2bohr(1.0) : 1.0;
		err_checkf(push_back_atom(f[0], scale * stod(f[3]), scale * stod(f[4]), scale * stod(f[5]), stoi(f[2])), "Error pushing back atom", file);
	}
	err_checkf(line.find("[STO]") == string::npos, "ERROR: STOs are not yet suupported!", file);
	err_checkf(ncen > 0, "No atoms in molden file", file);
	err_checkf(line.find("[GTO]") != string::npos, "Expected [GTO] after the atoms but found: '" + line + "'", file);
	//----------------------------- Basis: per atom "index 0", shells "type nprim 1.0", primitives, blank line ------------------------------
	//The contraction coefficients are taken as multiplying bare x^l exp(-a r^2), with the
	//contracted shell already normalised - what ORCA and orca_2mkl write. The format's own
	//documentation describes the other convention (coefficients multiply individually normalised
	//primitives, the contracted shell renormalised afterwards) and nothing in a molden file says
	//which one it is in: no file we have carries a "program=" keyword, and the [Title] line is
	//not evidence either - F2.molden has an empty title and is ORCA-convention.
	//
	//Do not add a norm-based detector without reading this first: the contracted self-overlap is
	//not an l-independent discriminator. Co2.molden's contracted d shell has bare self-overlap
	//1.000 for the xx component while Ce_full.molden's uncontracted d shells have 3.000 - both
	//ORCA files, differing only in which cartesian component the writer normalised. What does
	//discriminate is per-MO: the file's own MO vectors are orthonormal in whatever convention it
	//was written in. tests/src/MoldenConventionTests.cpp records those numbers for F2.molden
	//(1.00000 per MO and 14.000000 electrons under this reading, 0.889..1.344 and 13.727 under
	//the other) and pins the density against an evaluator independent of this code. Applying that
	//test inside the reader means building the contracted AO overlap, and not one file in the
	//corpus or the test set is in the other convention, so it is not built. Symptom if one ever
	//turns up: a density wrong by a factor of thousands, not by a little.
	int atoms_with_basis = 0;
	while (atoms_with_basis < ncen && (read_line_or_fail(rf, line, "the basis set", file), line.find("[") == string::npos))
	{
		svec f = fields();
		if (f.empty())
			continue;
		const int atom_based = stoi(f[0]) - 1;
		err_checkf(atom_based >= 0 && atom_based < ncen, "Basis set for atom " + to_string(atom_based + 1) + " of " + to_string(ncen), file);
		const string where = "the basis set of atom " + to_string(atom_based + 1);
		int shell = 0;
		while (read_line_or_fail(rf, line, where, file), !trim(line).empty())
		{
			f = fields();
			err_checkf(f.size() >= 2, "Bad shell line in " + where + ": '" + line + "'", file);
			//spectroscopic letters without j, as ORCA and pyscf write h and i shells
			const size_t shell_type = string("spdfghiklmn").find(static_cast<char>(tolower(f[0][0]))) + 1;
			err_checkf(f[0].size() == 1 && shell_type > 0, "Unknown shell type in " + where + ": '" + line + "'", file);
			const int number_of_functions = stoi(f[1]);
			err_checkf(number_of_functions > 0, "Shell without primitives in " + where + ": '" + line + "'", file);
			for (int i = 0; i < number_of_functions; i++)
			{
				read_line_or_fail(rf, line, where, file);
				vec v;
				append_numbers(line, v, where, file);
				err_checkf(v.size() == 2, "Expected 'exponent coefficient' in " + where + " but found: '" + line + "'", file);
				err_checkf(atoms[atom_based].push_back_basis_set(v[0], v[1], static_cast<int>(shell_type), shell), "Error pushing back basis", file);
			}
			shell++;
		}
		err_checkf(shell > 0, "No shells given for atom " + to_string(atom_based + 1), file);
		atoms_with_basis++;
	}
	err_checkf(atoms_with_basis == ncen, "Basis set given for " + to_string(atoms_with_basis) + " of " + to_string(ncen) + " atoms", file);
	//----------------------------- Flags up to [MO] ------------------------------
	bool d5 = false, f7 = false, g9 = false;
	while (line.find("[MO]") == string::npos)
	{
		string flag = line;
		transform(flag.begin(), flag.end(), flag.begin(), [](unsigned char c) { return static_cast<char>(toupper(c)); });
		d5 |= flag.find("5D") != string::npos;
		f7 |= flag.find("7F") != string::npos;
		g9 |= flag.find("9G") != string::npos;
		read_line_or_fail(rf, line, "[MO]", file);
	}
	const bool spherical = d5 && f7 && g9;
	err_checkf(spherical || (!d5 && !f7 && !g9), "Mixed cartesian and spherical shells in molden file are not supported", file);
	d_f_switch = !spherical;
	auto nfunc = [spherical](const int l) { return spherical ? constants::n_spher(l) : constants::n_cart(l); };
	//The primitives in file order, the size of the shell each belongs to and the coefficients per MO
	vector<primitive> prims;
	ivec shellsizes;
	int expected_coefs = 0;
	for (int a = 0; a < ncen; a++)
	{
		int current_shell = -1;
		for (unsigned int s = 0; s < atoms[a].get_basis_set_size(); s++)
		{
			if ((int)atoms[a].get_basis_set_shell(s) != current_shell)
			{
				const int l = atoms[a].get_basis_set_type(s) - 1;
				//the format defines the cartesian order only up to g; pure shells follow the gbw tables to l = 10
				err_checkf(l <= 10 && (spherical || l <= 4), "Molden shells beyond g are only supported as spherical harmonics", file);
				expected_coefs += nfunc(l);
				current_shell++;
			}
			shellsizes.push_back(atoms[a].get_shellcount(current_shell));
			prims.push_back(primitive(a + 1, atoms[a].get_basis_set_type(s), atoms[a].get_basis_set_exponent(s), atoms[a].get_basis_set_coefficient(s)));
		}
	}
	//----------------------------- MOs: "Key= value" header lines, then "index coefficient" lines ------------------------------
	//One row per MO of either spin: the density matrix is sum_i occ_i c_i c_i^T whatever the spin
	vec2 coefficients;
	vec occ, occ_beta;
	int nmo = 0;
	while (getline_universal(rf, line) && line.find("[") == string::npos)
	{
		if (trim(line).empty())
			continue;
		const string mo = "MO " + to_string(nmo + 1);
		double ene = 0.0, occup = 0.0;
		bool spin = false; // alpha = false, beta = true
		for (size_t eq = line.find('='); eq != string::npos; eq = line.find('='))
		{
			const string key = trim(line.substr(0, eq)), value = trim(line.substr(eq + 1));
			if (key == "Ene")
				ene = stod(value);
			else if (key == "Occup")
				occup = stod(value);
			else if (key == "Spin")
				spin = value != "Alpha" && value != "alpha";
			else
				err_checkf(key == "Sym", "Unknown key in the header of " + mo + ": '" + line + "'", file);
			read_line_or_fail(rf, line, "the header of " + mo, file);
		}
		if (spin)
			is_unrestricted = true;
		push_back_MO(nmo + 1, occup, ene, spin);
		occ.push_back(occup);
		occ_beta.push_back(spin ? occup : 0.0);
		coefficients.push_back(vec());
		int run = 0, basis_run = 0;
		vec2 shell;
		for (int i = 0; i < expected_coefs; i++)
		{
			if (i > 0)
				read_line_or_fail(rf, line, "coefficient " + to_string(i + 1) + " of " + mo, file);
			vec v;
			append_numbers(line, v, "coefficient " + to_string(i + 1) + " of " + mo, file);
			err_checkf(v.size() == 2 && static_cast<int>(v[0]) == i + 1, "Expected coefficient " + to_string(i + 1) + " of " + mo + " but found: '" + line + "'", file);
			coefficients.back().push_back(v[1]);
			const int l = prims[basis_run].get_type() - 1, size = shellsizes[basis_run], n = nfunc(l);
			if (run == 0)
				shell.assign(n, vec(size));
			for (int s = 0; s < size; s++)
				shell[run][s] = v[1] * prims[basis_run + s].get_coef();
			if (++run < n)
				continue;
			//p shells are x, y, z in the file whatever the [5D] flags (the gbw tables would read them as z, x, y)
			if (spherical && l > 1)
				push_back_spherical_shell(nmo, l, shell, prims, basis_run, size);
			else
				push_back_cartesian_shell(nmo, l, shell, prims, basis_run, size, molden_order[l], cart_norm(l).data());
			run = 0;
			basis_run += size;
		}
		nmo++;
	}
	//DM = C^T occ C over the nmo x nbf coefficient matrix, nbf x nbf whatever the MO count
	err_checkf(nmo > 0, "No MOs in molden file", file);
	vec _coefficients = flatten<double>(coefficients);
	dMatrix2 m_coefs = reshape<dMatrix2>(_coefficients, Shape2D(nmo, expected_coefs));
	dMatrix2 temp_co = diag_dot(m_coefs, occ, true);
	DM = dot(temp_co, m_coefs);
	if (is_unrestricted) {
		dMatrix2 temp_b = diag_dot(m_coefs, occ_beta, true);
		DM_beta = dot(temp_b, m_coefs);
	}
	set_exp_cutoff();
	return true;
};

bool WFN::read_tonto(const std::filesystem::path &filename, std::ostream &file, const bool debug, const std::filesystem::path &energies_filename, const std::filesystem::path &orbitals_filename)
{
	using namespace std;
	d_f_switch = true;
	err_checkf(std::filesystem::exists(filename), "couldn't open or find " + filename.string() + ", leaving", file);
	std::filesystem::path energies_file, orbitals_file, stdout_file;
	ifstream rf;
	string line;
	string scf_kind;
	bool restricted_search = true;
	if (energies_filename == "" && orbitals_filename == "") {
		if (filename.string().find("stdout") != std::string::npos) {
			string jobname;
			stdout_file = filename;
			rf.open(stdout_file.string().c_str(), ios::in);
			seek_line(rf, line, "Name ...", file);
			jobname = split_string<string>(line, " ")[2];
			seek_line(rf, line, "SCF kind ....", file);
			scf_kind = split_string<string>(line, " ")[3];
			if (scf_kind == "rhf" || scf_kind == "rks" || scf_kind == "xray_rhf" || scf_kind == "xray_rks")
				restricted_search = true;
			else //there could be other kinds like ghf, but at the moment let me expect either uhf/uks or rhf/rks
				restricted_search = false;

			method = scf_kind;

			energies_file = filename.parent_path() / (jobname + ".orbital_energies,restricted");
			if (!std::filesystem::exists(energies_file))
				energies_file = filename.parent_path() / (jobname + ".MO_energies,r");
			if (!std::filesystem::exists(energies_file)) {
				energies_file = filename.parent_path() / (jobname + ".orbital_energies,alpha");
				restricted_search = false;
			}
			if (!std::filesystem::exists(energies_file)) {
				energies_file = filename.parent_path() / (jobname + ".MO_energies,a");
				restricted_search = false;
			}

			orbitals_file = filename.parent_path() / (jobname + ".molecular_orbitals,restricted");
			if (!std::filesystem::exists(orbitals_file))
				orbitals_file = filename.parent_path() / (jobname + ".MOs,r");
			if (!std::filesystem::exists(orbitals_file)) {
				orbitals_file = filename.parent_path() / (jobname + ".molecular_orbitals,alpha");
				restricted_search = false;
			}
			if (!std::filesystem::exists(orbitals_file)) {
				orbitals_file = filename.parent_path() / (jobname + ".MOs,a");
				restricted_search = false;
			}
			err_checkf(std::filesystem::exists(orbitals_file), "couldn't open or find " + orbitals_file.string() + ", leaving", file);
			err_checkf(std::filesystem::exists(energies_file), "couldn't open or find " + energies_file.string() + ", leaving", file);
			//check if there is a stdout file in the same folder
			err_checkf(std::filesystem::exists(stdout_file), "couldn't open or find " + stdout_file.string() + ", leaving", file);
		}
		else if (filename.extension() == ".orbital_energies,restricted" || filename.extension() == ".MO_energies,r") {
			energies_file = filename;
			orbitals_file = filename;
			if (filename.extension() == ".orbital_energies,restricted") {
				orbitals_file.replace_extension(".molecular_orbitals,restricted");
			}
			else if (filename.extension() == ".MO_energies,r") {
				orbitals_file.replace_extension(".MOs,r");
			}
			stdout_file = filename.parent_path() / "stdout";
			err_checkf(std::filesystem::exists(orbitals_file), "couldn't open or find " + orbitals_file.string() + ", leaving", file);
			err_checkf(std::filesystem::exists(energies_file), "couldn't open or find " + energies_file.string() + ", leaving", file);
			//check if there is a stdout file in the same folder
			err_checkf(std::filesystem::exists(stdout_file), "couldn't open or find " + stdout_file.string() + ", leaving", file);
			rf.open(stdout_file.string().c_str(), ios::in);
		}
		else if (filename.extension() == ".molecular_orbitals,restricted" || filename.extension() == ".MOs,r") {
			orbitals_file = filename;
			energies_file = filename;
			if (filename.extension() == ".molecular_orbitals,restricted") {
				energies_file.replace_extension(".orbital_energies,restricted");
			}
			else if (filename.extension() == ".MOs,r") {
				energies_file.replace_extension(".MO_energies,r");
			}
			stdout_file = filename.parent_path() / "stdout";
			err_checkf(std::filesystem::exists(orbitals_file), "couldn't open or find " + orbitals_file.string() + ", leaving", file);
			err_checkf(std::filesystem::exists(energies_file), "couldn't open or find " + energies_file.string() + ", leaving", file);
			//check if there is a stdout file in the same folder
			err_checkf(std::filesystem::exists(stdout_file), "couldn't open or find " + stdout_file.string() + ", leaving", file);
			rf.open(stdout_file.string().c_str(), ios::in);
		}
		else {
			err_checkf(false, "Filename extension not recognized for tonto files! Please provide either .orbital_energies,restricted or .molecular_orbitals,restricted (or their short forms).", file);
		}
	}
	else {
		energies_file = energies_filename;
		orbitals_file = orbitals_filename;
		if (energies_file.string().find("alpha") != string::npos || energies_file.string().find(",a") != string::npos)
			restricted_search = false;
		stdout_file = filename;
		err_checkf(std::filesystem::exists(orbitals_file), "couldn't open or find " + orbitals_file.string() + ", leaving", file);
		err_checkf(std::filesystem::exists(energies_file), "couldn't open or find " + energies_file.string() + ", leaving", file);
		err_checkf(std::filesystem::exists(stdout_file), "couldn't open or find " + stdout_file.string() + ", leaving", file);
		rf.open(stdout_file.string().c_str(), ios::in);
	}

	const bool restricted = restricted_search;


	if (debug)
		file << "File is valid, continuing...\n" << GetCurrentDir << endl;
	origin = e_origin::tonto;
	isBohr = true;
	//open the files as read-only binary files
	ifstream rf_e(energies_file.c_str(), ios::binary);
	ifstream rf_o(orbitals_file.c_str(), ios::binary);
	err_checkf(rf_e.good(), "couldn't open " + energies_file.string() + ", leaving", file);
	err_checkf(rf_o.good(), "couldn't open " + orbitals_file.string() + ", leaving", file);
	if (rf_e.good())
		path = filename;
	rf_e.seekg(0);
	rf_o.seekg(0);

	//Read the energies and orbital coefficients from binary files
	vec energies, orbitals;
	read_block_from_fortran_binary(rf_e, energies);
	read_block_from_fortran_binary(rf_o, orbitals);

	vec energies_beta, orbitals_beta;
	if (!restricted) {
		rf_e.close();
		rf_o.close();
		std::string s = energies_file.string();
		if (auto pos = s.rfind("alpha"); pos != std::string::npos) {
			s.replace(pos, std::strlen("alpha"), "beta");
		}
		energies_file = std::filesystem::path(std::move(s));
		s = orbitals_file.string();
		if (auto pos = s.rfind("alpha"); pos != std::string::npos) {
			s.replace(pos, std::strlen("alpha"), "beta");
		}
		orbitals_file = std::filesystem::path(std::move(s));
		if (!std::filesystem::exists(energies_file)) {
			s = energies_file.string();
			if (auto pos = s.rfind(",a"); pos != std::string::npos) {
				s.replace(pos, std::strlen(",a"), ",b");
			}
			energies_file = std::filesystem::path(std::move(s));
			restricted_search = false;
		}
		if (!std::filesystem::exists(orbitals_file)) {
			s = orbitals_file.string();
			if (auto pos = s.rfind(",a"); pos != std::string::npos) {
				s.replace(pos, std::strlen(",a"), ",b");
			}
			orbitals_file = std::filesystem::path(std::move(s));;
			restricted_search = false;
		}
		rf_e.open(energies_file.string().c_str(), ios::binary);
		rf_o.open(orbitals_file.string().c_str(), ios::binary);
		err_checkf(rf_e.good(), "couldn't open " + energies_file.string() + ", leaving", file);
		err_checkf(rf_o.good(), "couldn't open " + orbitals_file.string() + ", leaving", file);
		read_block_from_fortran_binary(rf_e, energies_beta);
		read_block_from_fortran_binary(rf_o, orbitals_beta);
		is_unrestricted = true;
	}

	rf_e.close();
	rf_o.close();

	//Now read the stdout file to get the atomic positions and basis set
	rf.seekg(0);
	err_checkf(rf.good(), "couldn't open " + stdout_file.string() + ", leaving", file);
	//Tonto prints "Label .... value": skip lines, then the last field of the line is the value
	auto last_field = [&](const int skip, const string &what)
	{
		for (int i = 0; i < skip; i++)
			read_line_or_fail(rf, line, what, file);
		svec f = split_string<string>(line, " ");
		remove_empty_elements(f);
		err_checkf(!f.empty(), "Empty line where " + what + " was expected", file);
		return f.back();
	};
	seek_line(rf, line, "Molecule information", file);
	charge = stoi(last_field(8, "the charge"));
	multi = stoi(last_field(1, "the multiplicity"));
	const int expected_atoms = stoi(last_field(2, "the atom count"));
	const int expected_electrons = stoi(last_field(1, "the electron count"));
	seek_line(rf, line, "Atom coordinates", file);
	svec line_digest;
	//skip 11 lines to get to the atom list
	for (int i = 0; i < 11; i++) {
		getline_universal(rf, line);
		if (line.find("This molecule has non trivial group") != string::npos)
			i -= 3;
	}
	line_digest = split_string<string>(line, " ");
	remove_empty_elements(line_digest);
	err_checkf(line_digest[0] == "1", "Atom list does not start with one? let me stop...", std::cout);

	while (!line.empty() && line.front() != '_')
	{
		line_digest = split_string<string>(line, " ");
		remove_empty_elements(line_digest);
		string label = line_digest[1];
		int atomic_number = static_cast<int>(stod(line_digest[2]));
		double x, y, z;
		if (line_digest.size() == 6) {
			x = constants::ang2bohr(stod(line_digest[3]));
			y = constants::ang2bohr(stod(line_digest[4]));
			z = constants::ang2bohr(stod(line_digest[5]));

		}
		else if (line_digest.size() == 7) {
			x = constants::ang2bohr(stod(line_digest[4]));
			y = constants::ang2bohr(stod(line_digest[5]));
			z = constants::ang2bohr(stod(line_digest[6]));
		}
		err_checkf(push_back_atom(label, x, y, z, atomic_number), "Error pushing back an atom!", std::cout);
		getline_universal(rf, line);
	}
	err_checkf(ncen == expected_atoms, "Did not read expected number of atoms!", std::cout);

	int alpha_els = 0, beta_els = 0, temp_els = get_nr_electrons();
	while (temp_els > 1)
	{
		alpha_els++;
		beta_els++;
		temp_els -= 2;
		if (debug)
			file << temp_els << "\n";
		err_checkf(alpha_els >= 0 && beta_els >= 0, "Error setting alpha and beta electrons! a or b are negative!", file);
		err_checkf(alpha_els + beta_els <= get_nr_electrons(), "Error setting alpha and beta electrons! Sum a + b > elcount!", file);
		err_checkf(temp_els > -int(get_nr_electrons()), "Error setting alpha and beta electrons! Ran below -elcount!", file);
	}
	alpha_els += temp_els;
	if (debug)
		file << "al/be els:" << alpha_els << " " << beta_els << std::endl;
	const int mult = get_multi();
	int diff = 0;
	if (mult != 0)
		diff = get_multi() - 1;
	if (debug)
		file << "diff: " << diff << std::endl;
	while (alpha_els - beta_els != diff)
	{
		alpha_els++;
		beta_els--;
		err_checkf(alpha_els >= 0 && beta_els >= 0, "Error setting alpha and beta electrons: " + std::to_string(alpha_els) + "/" + std::to_string(beta_els), file);
	}

	seek_line(rf, line, "Gaussian basis sets", file);
	basis_set_name = last_field(3, "the basis set name");
	const int no_basis_sets = stoi(last_field(2, "the number of basis sets"));
	const int no_shells = stoi(last_field(1, "the number of shells"));
	const int no_bf = stoi(last_field(2, "the number of basis functions")); //the shell-pair line before it is skipped
	const int no_prim = stoi(last_field(1, "the number of primitives"));
	std::vector<std::pair<std::string, atom>> basis_set_data;
	std::map<char, int> l_map = { {'S',0}, {'P',1}, {'D',2}, {'F',3}, {'G',4}, {'H',5}, {'I',6}, {'s',0}, {'p',1}, {'d',2}, {'f',3}, {'g',4}, {'h',5}, {'i',6} };
	for (int nbs = 0; nbs < no_basis_sets; nbs++)
	{
		//two empty lines
		getline_universal(rf, line);
		getline_universal(rf, line);
		getline_universal(rf, line);//looks like "Basis set H:3-21G"
		line_digest = split_string<string>(line, " ");
		const string atom_type = split_string<string>(line_digest[2], ":")[0];
		getline_universal(rf, line); // empty line
		getline_universal(rf, line);//looks like "No. of shells .... N"
		line_digest = split_string<string>(line, " ");
		const int shells_local = stoi(line_digest[4]);
		getline_universal(rf, line);//looks like "No. of basis functions .... N"
		line_digest = split_string<string>(line, " ");
		const int bfs_local = stoi(line_digest[5]);
		getline_universal(rf, line);//looks like "No. of primitives .... N"
		line_digest = split_string<string>(line, " ");
		const int prims_local = stoi(line_digest[4]);
		/*
__________________________________

 -L-   Fn    Exponent  Contraction
		#         /au          /au
__________________________________

*/
		for (int i = 0; i < 7; i++) getline_universal(rf, line); //skip 6 lines to get to the shells
		atom temp_at(atom_type, {}, 0, 0, 0, 0, constants::get_Z_from_label(atom_type.c_str()));
		for (int s = 0; s < shells_local; s++)
		{
			//getline_universal(rf, line); //get shell line
			line_digest = split_string<string>(line, " ");
			remove_empty_elements(line_digest);
			err_checkf(l_map.contains(line_digest[0][0]), "Angular momentum not found: " + line_digest[0], std::cout);
			const int angul = l_map.at(line_digest[0][0]); // safe because contains returned true
			const int n_prim = stoi(line_digest[1]);
			do {
				const double exponent = stod(line_digest[line_digest.size() == 2 ? 0 : 2]);
				const double coefficient = stod(line_digest[line_digest.size() == 2 ? 1 : 3]);
				const double norm_fac = pow(pow(2, 3 + 4 * angul) * pow(exponent, 2 * angul + 3) / constants::PI3 / pow(constants::double_ft[angul], 2), 0.25);
				temp_at.push_back_basis_set(exponent, coefficient * norm_fac, angul + 1, s);
				getline_universal(rf, line);
				line_digest = split_string<string>(line, " ");
				remove_empty_elements(line_digest);
			} while (line_digest.size() == 2);
		}
		err_checkf(line[0] == '_', "Expected a line of underscores after basis set for atom " + atom_type, file);
		basis_set_data.push_back(std::make_pair(atom_type, temp_at));
	}
	//Now assign the basis set data to the atoms
	for (int a = 0; a < ncen; a++)
	{
		for (const auto &[atom_type, atom_template] : basis_set_data)
		{
			if (atom_type == constants::atnr2letter(atoms[a].get_charge()))
			{
				//copy atom basis set information
				atoms[a].set_basis_set(atom_template.get_basis_set());
				atoms[a].set_shellcount(atom_template.get_shellcount());
				break;
			}
		}
	}

	int expected_coefs = 0;
	vector<primitive> prims;
	ivec temp_shellsizes;
	for (int a = 0; a < ncen; a++)
	{
		int current_shell = -1;
		for (unsigned int s = 0; s < atoms[a].get_basis_set_size(); s++)
		{
			if ((int)atoms[a].get_basis_set_shell(s) != current_shell)
			{
				if (atoms[a].get_basis_set_type(s) == 1)
				{
					expected_coefs++;
				}
				else if (atoms[a].get_basis_set_type(s) == 2)
				{
					expected_coefs += 3;
				}
				else if (atoms[a].get_basis_set_type(s) == 3)
				{
					expected_coefs += 6;
				}
				else if (atoms[a].get_basis_set_type(s) == 4)
				{
					expected_coefs += 10;
				}
				else if (atoms[a].get_basis_set_type(s) == 5)
				{
					expected_coefs += 15;
				}
				current_shell++;
			}
			temp_shellsizes.push_back(atoms[a].get_shellcount(current_shell));
			prims.push_back(primitive(a + 1,
				atoms[a].get_basis_set_type(s),
				atoms[a].get_basis_set_exponent(s),
				atoms[a].get_basis_set_coefficient(s)));
		}
	}
	err_checkf(expected_coefs == no_bf, "Expected number of basis functions (" + to_string(expected_coefs) + ") does not match number in file (" + to_string(no_bf) + ")!", file);

	const unsigned int nr_operators = restricted ? 1 : 2;

	for (int op = 0; op < nr_operators; op++) {
		if (debug)
			file << "Reading " << (op == 0 ? "alpha/restricted" : "beta") << " orbitals..." << std::endl;

		dMatrix2 coefficients;
		vec occ(expected_coefs, 0.0);
		if (op == 0)
			coefficients = reshape<dMatrix2>(orbitals, Shape2D(expected_coefs, expected_coefs));
		else
			coefficients = reshape<dMatrix2>(orbitals_beta, Shape2D(expected_coefs, expected_coefs));

		for (int MO_run = 0; MO_run < expected_coefs; MO_run++)
		{
			if (op == 0) {
				if (restricted) {
					if (MO_run < alpha_els && multi == 1) { //RHF all paired
						push_back_MO(MO_run, 2.0, energies[MO_run], 0);
						occ[MO_run] = 2.0;
					}
					else if (multi != 1 && MO_run < beta_els) { // RHF only 2 until beta electrons
						push_back_MO(MO_run, 2.0, energies[MO_run], 0);
						occ[MO_run] = 2.0;
					}
					else if (multi != 1 && MO_run >= beta_els && MO_run < alpha_els) { // RHF only 2 until beta electrons
						push_back_MO(MO_run, 1.0, energies[MO_run], 0);
						occ[MO_run] = 1.0;
					}
					else
						push_back_MO(MO_run, 0.0, energies[MO_run], 0);
				}
				else {
					if (MO_run < alpha_els) { // UHF
						push_back_MO(MO_run, 1.0, energies[MO_run], 0);
						occ[MO_run] = 1.0;
					}
					else
						push_back_MO(MO_run, 0.0, energies[MO_run], 0);
				}

			}
			else {
				if (MO_run < beta_els) {
					push_back_MO(expected_coefs + MO_run, 1.0, energies_beta[MO_run], 1);
					occ[MO_run] = 1.0;
				}
				else
					push_back_MO(expected_coefs + MO_run, 0.0, energies_beta[MO_run], 1);
			}
			const int mo = op == 0 ? MO_run : MO_run + expected_coefs;
			int basis_run = 0, run = 0;
			vec2 shell;
			for (int i = 0; i < expected_coefs; i++)
			{
				const int l = prims[basis_run].get_type() - 1, size = temp_shellsizes[basis_run];
				err_checkf(l <= 4, "tonto shells beyond g are not supported", file);
				if (run == 0)
					shell.assign(constants::n_cart(l), vec(size));
				for (int s = 0; s < size; s++)
					shell[run][s] = coefficients(MO_run, i) * prims[basis_run + s].get_coef() * tonto_scale[l];
				if (++run < constants::n_cart(l))
					continue;
				push_back_cartesian_shell(mo, l, shell, prims, basis_run, size, l == 3 ? tonto_f_order : nullptr);
				run = 0;
				basis_run += size;
			}
			err_checkf(run == 0, "There should not be any unfinished shells! Aborting reading tonto file after MO " + to_string(MO_run) + "!", file);
		}
		dMatrix2 temp_co = diag_dot(coefficients, occ, true);
		if (op == 0)
			DM = dot(temp_co, coefficients);
		else {
			dMatrix2 DM_beta = dot(temp_co, coefficients);
			for (int i = 0; i < expected_coefs; i++)
				for (int j = 0; j < expected_coefs; j++)
					DM(i, j) += DM_beta(i, j);
		}
	}
	set_exp_cutoff();
	return true;
};

template<typename T>
void move_columns(T &matrix, int dimension, int from_col, int num_cols, int to_col) {
	if (from_col == to_col || num_cols == 0) return;

	int block_size = num_cols * dimension;
	T tmp_block(block_size);

	// Copy the block to be moved
	std::copy(matrix.data() + from_col * dimension,
		matrix.data() + (from_col + num_cols) * dimension,
		tmp_block.data());

	if (from_col < to_col) {
		// Shift the intermediate block leftward
		std::copy(matrix.data() + (from_col + num_cols) * dimension,
			matrix.data() + to_col * dimension,
			matrix.data() + from_col * dimension);

		// Insert the moved block
		std::copy(tmp_block.begin(),
			tmp_block.end(),
			matrix.data() + (to_col - num_cols) * dimension);
	}
	else {
		// Shift the intermediate block rightward
		std::copy_backward(matrix.data() + to_col * dimension,
			matrix.data() + from_col * dimension,
			matrix.data() + (from_col + num_cols) * dimension);
		// Insert the moved block
		std::copy(tmp_block.begin(),
			tmp_block.end(),
			matrix.data() + to_col * dimension);
	}
}


bool WFN::read_gbw(const std::filesystem::path &filename, std::ostream &file, const bool debug, const bool _has_ECPs)
{
	using namespace std;
	// Details form https://orcaforum.kofo.mpg.de/viewtopic.php?f=8&t=3299&start=20
	err_checkf(std::filesystem::exists(filename), "couldn't open or find " + filename.string() + ", leaving", file);
	if (debug)
		file << "File is valid, continuing...\n"
		<< GetCurrentDir << endl;
	origin = e_origin::gbw;
	isBohr = true;
	ifstream rf(filename.c_str(), ios::binary);
	if (rf.good())
		path = filename;
	const int64_t file_size = static_cast<int64_t>(std::filesystem::file_size(filename));
	//A section pointer or count from a damaged file used to be followed blindly
	auto check_offset = [&](const int64_t offset, const string &what) {
		err_checkf(offset > 0 && offset < file_size, what + " points to byte " + to_string(offset) + " of a " + to_string(file_size) + " byte file", file);
	};
	auto check_count = [&](const int64_t count, const string &what) {
		err_checkf(count > 0 && count < file_size, what + " of " + to_string(count) + " in a " + to_string(file_size) + " byte file", file);
	};
	//One binary field; every read is checked so a short file fails here and not on the garbage it would leave behind
	auto rd = [&](void *dst, const std::streamsize bytes, const string &what) {
		rf.read(static_cast<char *>(dst), bytes);
		err_checkf(rf.good(), "Error reading " + what, file);
	};
	string line;
	int geo_start_bit = 8;
	int basis_start_bit = 16;
	int MO_start_bit = 24;
	int ECP_start_bit = 32;
	int soi = constants::soi;
	int geo_int_lim = 5;

	try
	{
		rf.seekg(0, ios::beg);
		int64_t magic = 0;
		rd(&magic, sizeof(magic), "magic number");
		if (magic == -1) {
			geo_start_bit += 24;
			basis_start_bit += 24;
			MO_start_bit += 24;
			ECP_start_bit += 24;
			soi = 8;
			geo_int_lim = 1;
		}
		// Reading geometry
		rf.seekg(geo_start_bit, ios::beg);
		int64_t geo_start = 0;
		rd(&geo_start, sizeof(geo_start), "geo_start");
		check_offset(geo_start, "geometry section");
		if (debug)
			file << "I read the pointer of geometry successfully" << endl;
		rf.seekg(geo_start, ios::beg);
		int at = 0;
		rd(&at, constants::soi, "atom count");
		check_count(at, "atom count");
		double geo_vals[6]{ 0, 0, 0, 0, 0, 0 }; // x,y,z, ch, exp_fin_nuc, mass
		// Use int64_t to safely hold soi-sized reads (soi may be 4 or 8 bytes)
		int64_t geo_ints[5]{ 0, 0, 0, 0, 0 };
		for (int a = 0; a < at; a++)
		{
			for (int i = 0; i < 6; i++)
			{
				rd(&geo_vals[i], constants::sod, "geo_val");
			}
			for (int i = 0; i < geo_int_lim; i++)
			{
				geo_ints[i] = 0;
				rd(&geo_ints[i], soi, "geo_int");
			}
			string temp = constants::atnr2letter(static_cast<int>(geo_ints[0]));
			err_checkf(temp != "PROBLEM", "Problem identifying atoms!", std::cout);
			err_checkf(push_back_atom(temp,
				geo_vals[0],
				geo_vals[1],
				geo_vals[2],
				static_cast<int>(geo_ints[0])),
				"Error pushing back atom", file);
		}
		if (debug)
			file << "I read the geometry of " << at << " atoms successfully" << endl;

		rf.seekg(basis_start_bit, ios::beg);
		int64_t basis_start = 0;
		rd(&basis_start, constants::soli, "basis_start");
		check_offset(basis_start, "basis set section");
		if (debug)
			file << "I read the pointer of basis set successfully" << endl;
		rf.seekg(basis_start, ios::beg);
		int atoms2 = 0, temp = 0;
		rd(&temp, constants::soi, "basis header");
		rd(&atoms2, constants::soi, "atoms2");
		err_checkf(atoms2 == at, "Basis set for " + to_string(atoms2) + " atoms but geometry for " + to_string(at), file);
		// long unsigned int atoms_with_basis = 0;
		vec exp(37, 0);
		vec con(37, 0);
		for (int a = 0; a < atoms2; a++)
		{
			int atom_based = 0, nr_shells = 0;
			rd(&atom_based, constants::soi, "atom_based");
			err_checkf(atom_based >= 0 && atom_based < ncen, "Basis set for atom " + to_string(atom_based + 1) + " of " + to_string(ncen), file);
			rd(&nr_shells, constants::soi, "nr_shells");
			check_count(nr_shells, "shell count");
			int shell = 0;
			for (int p = 0; p < nr_shells; p++)
			{
				int ang_mom = 0, coeff_ind = 0, nr_funct = 0, center = 0;
				rd(&ang_mom, constants::soi, "ang_mom");
				err_checkf(ang_mom <= 10, "Higher angular momentum basis functions than l = 10", file);
				rd(&coeff_ind, constants::soi, "coeff_ind");
				rd(&nr_funct, constants::soi, "nr_func");
				rd(&center, constants::soi, "center");
				for (int b = 0; b < 37; b++)
				{
					rd(&exp[b], constants::sod, "exp");
				}
				for (int b = 0; b < 37; b++)
				{
					rd(&con[b], constants::sod, "con");
					if (exp[b] != 0 && con[b] != 0)
					{
						err_checkf(atoms[atom_based].push_back_basis_set(exp[b], con[b], ang_mom + 1, shell), "Error pushing back basis", file);
					}
				}
				shell++;
			}
			// atoms_with_basis++;
		}
		int expected_coefs = 0;
		vector<primitive> prims;
		ivec temp_shellsizes;
		for (int a = 0; a < ncen; a++)
		{
			int current_shell = -1;
			for (unsigned int s = 0; s < atoms[a].get_basis_set_size(); s++)
			{
				if ((int)atoms[a].get_basis_set_shell(s) != current_shell)
				{
					expected_coefs += 2 * atoms[a].get_basis_set_type(s) - 1;
					current_shell++;
				}
				temp_shellsizes.push_back(atoms[a].get_shellcount(current_shell));
				prims.push_back(primitive(a + 1,
					atoms[a].get_basis_set_type(s),
					atoms[a].get_basis_set_exponent(s),
					atoms[a].get_basis_set_coefficient(s)));
			}
		}
		basis_set_name = "GBW read basis set";
		// int norm_const_run = 0;
		int MO_run = 0;
		if (debug)
			file << "I read the basis of " << atoms2 << " atoms successfully" << endl;

		rf.seekg(MO_start_bit, ios::beg);
		int64_t MOs_start = 0;
		rd(&MOs_start, constants::soli, "MO_start");
		check_offset(MOs_start, "MO section");
		if (debug)
			file << "I read the pointer of MOs successfully" << endl;
		rf.seekg(MOs_start, ios::beg);
		int operators = 0;
		int64_t dimension_i64 = 0;
		rd(&operators, constants::soi, "operators");
		rd(&dimension_i64, soi, "MO dimension");
		int dimension = static_cast<int>(dimension_i64);
		err_checkf(operators == 1 || operators == 2, "gbw with " + to_string(operators) + " operators", file);
		err_checkf(dimension == expected_coefs, "MO matrix dimension " + to_string(dimension) + " but the basis set has " + to_string(expected_coefs) + " functions", file);
		size_t coef_nr = size_t(dimension) * size_t(dimension);
		vec2 coefficients(operators);
		vec2 occupations(operators);
		vec2 energies(operators);
		ivec2 irreps(operators);
		ivec2 cores(operators);
		if (operators > 1)
			is_unrestricted = true;
		for (int i = 0; i < operators; i++)
		{
			coefficients[i].resize(coef_nr, 0);
			occupations[i].resize(dimension, 0);
			energies[i].resize(dimension, 0);
			irreps[i].resize(dimension, 0);
			cores[i].resize(dimension, 0);
			if (debug)
				file << "operators: " << operators << " coef_nr: " << coef_nr << " dimension: " << dimension << "\n";
			rd(coefficients[i].data(), constants::sod * coef_nr, "coefficients");
			if (debug)
				file << "I read the coefficients successfully\n";
			rd(occupations[i].data(), constants::sod * dimension, "occupations");
			if (debug)
				file << "I read the occupations successfully\n";
			rd(energies[i].data(), constants::sod * dimension, "energies");
			if (debug)
				file << "I read the energies successfully\n";
			rd(irreps[i].data(), constants::soi * dimension, "irreps");
			if (debug)
				file << "I read the irreps successfully\n";
			rd(cores[i].data(), constants::soi * dimension, "cores");
			if (debug)
			{
				file << "I read the cores successfully\nI am expecting " << expected_coefs << " coefficients per MO\n";
			}
			//Without this, expected_coefs is set to 0 after push_back_MO... no idea why
			const double constant_expected_coefs = expected_coefs;
			for (int j = 0; j < dimension; j++)
			{
				push_back_MO(i * dimension + j + 1, occupations[i][j], energies[i][j], i);
				// int run_coef = 0;
				int run = 0, basis_run = 0;
				vec2 shell;
				for (int p = 0; p < constant_expected_coefs; p++)
				{
					const int l = prims[basis_run].get_type() - 1, nsph = constants::n_spher(l), size = temp_shellsizes[basis_run];
					err_checkf(l <= 10, "Types higher than l = 10 in gbws", file);
					if (run == 0) shell.assign(nsph, vec(size));
					for (int s = 0; s < size; s++)
						shell[run][s] = coefficients[i][j + p * dimension] * prims[basis_run + s].get_coef();
					if (++run < nsph) continue;
					push_back_spherical_shell(MO_run, l, shell, prims, basis_run, size);
					run = 0;
					basis_run += size;
				}
				err_checkf(run == 0, "There should not be any unfinished shells! Aborting reading gbw file after MO " + to_string(MO_run) + "!", file);
				MO_run++;
			}
		}


		dMatrix2 reorderd_coefs_s1(dimension, dimension), reorderd_coefs_s2;
		if (operators == 2) reorderd_coefs_s2 = dMatrix2(dimension, dimension);

		dMatrixRef2 coefs_2D_s1_span(coefficients[0].data(), dimension, dimension);
		// coefficients[1] only exists when operators==2; use a harmless fallback otherwise
		// to avoid a Debug-mode vector bounds assertion (the span is never accessed for op==1).
		dMatrixRef2 coefs_2D_s2_span(
			operators == 2 ? coefficients[1].data() : coefficients[0].data(),
			dimension, dimension);
		int index = 0;
		for (const atom &_atom : atoms) {
			std::vector<basis_set_entry> basis = _atom.get_basis_set();
			int temp_bas_idx = 0;
			for (unsigned int shell = 0; shell < _atom.get_shellcount_size(); shell++) {
				int type = basis[temp_bas_idx].get_type() - 1;
				temp_bas_idx += _atom.get_shellcount(shell);
				for (int m_idx = 0; m_idx < 2 * type + 1; m_idx++) {
					int offset = constants::orca_2_pySCF(type, m_idx).value();

					auto coefs_2D_s1_slice = Kokkos::submdspan(coefs_2D_s1_span, index + m_idx, Kokkos::full_extent);
					auto reord_coefs_slice = Kokkos::submdspan(reorderd_coefs_s1.to_mdspan(), index + offset, Kokkos::full_extent);
					std::copy(coefs_2D_s1_slice.data_handle(), coefs_2D_s1_slice.data_handle() + dimension, reord_coefs_slice.data_handle());
					if (operators == 2) {
						auto coefs_2D_s2_slice = Kokkos::submdspan(coefs_2D_s2_span, index + m_idx, Kokkos::full_extent);
						reord_coefs_slice = Kokkos::submdspan(reorderd_coefs_s2.to_mdspan(), index + offset, Kokkos::full_extent);
						std::copy(coefs_2D_s2_slice.data_handle(), coefs_2D_s2_slice.data_handle() + dimension, reord_coefs_slice.data_handle());
					}
				}
				index += 2 * type + 1;
			}
		}

		//Map to collect the end index of every type
		index = 0;
		int atom_index = 0;
		for (atom &_atom : atoms) {
			std::map<int, int> type_end;
			std::vector<basis_set_entry> basis = _atom.get_basis_set();
			int temp_bas_idx = 0;
			for (unsigned int shell = 0; shell < _atom.get_shellcount_size(); shell++) {
				int type = basis[temp_bas_idx].get_type() - 1;
				temp_bas_idx += _atom.get_shellcount(shell);
				if (type_end.find(type + 1) == type_end.end()) {
					type_end[type] = index + 2 * type + 1;
				}
				else {
					if (debug) {
						file << "GBW reader: moving out-of-order angular momentum block l=" << type
							 << " for atom " << atom_index + 1 << ", shell " << shell + 1
							 << " into the internal coefficient order.\n";
					}
					move_columns(reorderd_coefs_s1.container(), dimension, index, 2 * type + 1, type_end[type]);

					if (operators == 2) {
						move_columns(reorderd_coefs_s2.container(), dimension, index, 2 * type + 1, type_end[type]);
					}

					//Add one to each value in the map
					for (int i = type; i <= type + 1; i++) {
						type_end[i] += 2 * type + 1;
					}
				}
				index += 2 * type + 1;
			}
			atom_index++;
		}



		int n_occ = 0;
		for (int i = 0; i < occupations[0].size(); i++) { if (occupations[0][i] > 0.0) n_occ++; }

		dMatrix2 coeff_mo_s1(dimension, dimension), coeff_small_s1(dimension, dimension);
		dMatrix2 coeff_mo_s2, coeff_small_s2;
		if (operators == 2)  coeff_mo_s2 = dMatrix2(dimension, dimension); coeff_small_s2 = dMatrix2(dimension, dimension);

		for (int i = 0; i < dimension; i++) {
			for (int oc = 0; oc < occupations[0].size(); oc++) {
				if (occupations[0][oc] <= 0.0) continue;
				coeff_mo_s1(i, oc) = reorderd_coefs_s1(i, oc) * occupations[0][oc];
				coeff_small_s1(i, oc) = reorderd_coefs_s1(i, oc);

				if (operators == 2) coeff_mo_s2(i, oc) = reorderd_coefs_s2(i, oc) * occupations[1][oc];
				if (operators == 2) coeff_small_s2(i, oc) = reorderd_coefs_s2(i, oc);
			}
		}

		if (operators == 1) {
			DM = dot(coeff_mo_s1, coeff_small_s1, false, true);
		}
		else {
			dMatrix2 DM_s1 = dot(coeff_mo_s1, coeff_small_s1, false, true);
			dMatrix2 DM_s2 = dot(coeff_mo_s2, coeff_small_s2, false, true);

			std::transform(DM_s1.container().begin(), DM_s1.container().end(), DM_s2.data(), DM_s1.data(), std::plus<double>());

			DM = DM_s1;
			DM_beta = DM_s2;
		}

		if (debug)
		{
			file << "\nI read " << MO_run << "/" << dimension << " MOs of " << operators << " operators successfully" << endl;
			file << "There are " << nex << " primitives after conversion" << endl;
		}
		if (_has_ECPs)
		{
			has_ECPs = true;
			vector<ECP_primitive> ECP_prims;
			// Reading ECPs? The pointer sits after the MO pointer in both layouts (byte 32, or 56 with magic -1)
			rf.seekg(ECP_start_bit, ios::beg);
			int64_t ECP_start = 0;
			rd(&ECP_start, sizeof(ECP_start), "ECP pointer");
			err_checkf(ECP_start != 0, "Could not read ECP information location from GBW file!", file);
			check_offset(ECP_start, "ECP section");
			if (debug)
				file << "I read the pointer of ECP successfully" << endl;
			rf.seekg(ECP_start, ios::beg);
			int64_t i1 = 0; //8 bytes in the file, a long is 4 on MSVC
			int i2 = 0;
			const int soi = 4;
			const int sod = 8;
			rd(&i1, 8, "ECP count");
			check_count(i1, "ECP count");
			err_checkf(i1 <= (int64_t)atoms.size(), "ECP block lists more atoms than the geometry", file);
			file << "First line: " << i1 << endl;
			for (int i = 0; i < i1; i++)
			{
				rd(&i2, 1, "ECP flag");
				//an atom without ECP is the flag byte alone, the record fields only follow flag 1
				if (i2 == 0)
					continue;
				int Z = 0;
				int nr_core = 0;
				int temp_0 = 0;
				int max_contract = 0;
				int max_angular = 0;
				int exps = 0;
				double n = 0;
				int center = 0;
				int type = 0;
				double e = 0;
				double c = 0;
				rd(&Z, soi, "Z");
				err_checkf(Z > 0, "Error reading Z in ECPs", file);
				rd(&temp_0, soi, "temp_0");
				err_checkf(temp_0 > 0, "Error reading temp_0 in ECPs", file);
				char *temp_c = new char[temp_0];
				rd(temp_c, temp_0, "temp_c in ECPs");
				rd(&nr_core, soi, "nr_core");
				err_checkf(nr_core >= 0, "Error reading nr_core in ECPs", file);
				atoms[i].set_ECP_electrons(nr_core);
				rd(&max_contract, soi, "max_contract");
				err_checkf(max_contract > 0, "Error reading max_contract in ECPs", file);
				rd(&max_angular, soi, "max_angular");
				err_checkf(max_angular > 0, "Error reading max_angular in ECPs", file);
				rd(&center, soi, "center");
				err_checkf(center >= 0 && center < (int)atoms.size(), "Error reading center in ECPs", file); //0-based atom index
				file << "I read " << Z << " " << temp_0 << " " << nr_core << " " << max_contract << " " << max_angular << "\n";
				for (int l = 0; l < max_angular; l++)
				{
					rd(&exps, soi, "exps");
					err_checkf(exps > 0, "Error reading exps in ECPs", file);
					rd(&type, soi, "type");
					err_checkf(type >= 0, "Error reading type in ECPs", file);
					err_checkf(type < 200, "This type will give me a headache...", file);
					file << "There are " << exps << " exponents of type " << type << " for angular momentum " << l << "\n";
					for (int fun = 0; fun < exps; fun++)
					{

						//Stuttgart order: r power, exponent, coefficient (Hg def2 f term: 2, 3.8857911, -30.3649964);
						//coefficients reach 275 and are often negative, so only garbage is rejected
						rd(&n, sod, "n");
						err_checkf(std::isfinite(n) && n >= 0 && n < 10, "Unreasonable r power in ECPs: " + to_string(n), file);
						rd(&e, sod, "e");
						err_checkf(std::isfinite(e) && e > 0, "Unreasonable exponent in ECPs: " + to_string(e), file);
						rd(&c, sod, "c");
						err_checkf(std::isfinite(c), "Unreasonable coefficient in ECPs", file);
						file << fun << " " << c << " " << e << " " << n << "\n";
						ECP_prims.push_back(ECP_primitive(center, type, e, c, static_cast<int>(n)));
					}
				}
				for (int _i = 0; _i < ncen; _i++)
					if (atoms[_i].get_charge() == Z)
						atoms[_i].set_ECP_electrons(nr_core);
			}
			if (debug)
			{
				file << "Ended reading" << endl;
			}
		}
	}
	catch (const exception &e)
	{
		err_checkf(false, "Error during reading of the gbw file! " + string(e.what()), file);
	}
	set_exp_cutoff();
	return true;
};

const vec WFN::get_norm_const(std::ostream &file, bool debug) const
{
	err_checkf(get_nr_basis_set_loaded() != 0, "No basis set loaded!", file);
	err_checkf(get_nr_basis_set_loaded() == get_ncen(), "Not all atoms have a basis set loaded!", file);
	vec norm_const;
	//-------------------normalize the basis set shell wise into a copy vector---------
	vec2 basis_coefficients;
	basis_coefficients.resize(ncen);
	for (int a = 0; a < ncen; a++)
		for (int p = 0; p < get_atom_primitive_count(a); p++)
		{
			double temp = get_atom_basis_set_exponent(a, p);
			switch (get_atom_primitive_type(a, p))
			{
			case 1:
				temp = 2 * temp / constants::PI;
				temp = pow(temp, 0.75);
				temp = temp * get_atom_basis_set_coefficient(a, p);
				basis_coefficients[a].push_back(temp);
				break;
			case 2:
				temp = 128 * pow(temp, 5);
				temp = temp / constants::PI3;
				temp = pow(temp, 0.25);
				temp = get_atom_basis_set_coefficient(a, p) * temp;
				basis_coefficients[a].push_back(temp);
				break;
			case 3:
				temp = 2048 * pow(temp, 7);
				temp = temp / (9 * constants::PI3);
				temp = pow(temp, 0.25);
				temp = get_atom_basis_set_coefficient(a, p) * temp;
				basis_coefficients[a].push_back(temp);
				break;
			case 4:
				temp = 32768 * pow(temp, 9);
				temp = temp / (225 * constants::PI3);
				temp = pow(temp, 0.25);
				temp = get_atom_basis_set_coefficient(a, p) * temp;
				basis_coefficients[a].push_back(temp);
				break;
			case -1:
				std::cout << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy...\n";
				break;
			}
		}
	for (int a = 0; a < ncen; a++)
	{
		double factor = 0.0;
		for (int s = 0; s < get_atom_shell_count(a); s++)
		{
			int type_temp = get_shell_type(a, s);
			if (type_temp == -1)
			{
				std::cout << "ERROR in type assignement!!\n";
			}
			if (debug)
			{
				std::cout << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << "\n"
					<< "start: " << get_shell_start(a, s)
					<< " stop: " << get_shell_end(a, s) << "\n"
					<< "factor: ";
			}
			switch (type_temp)
			{
			case 1:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						double aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						double term = (constants::PI / aiaj);
						term = pow(term, 1.5);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
					}
				}
				err_checkf(factor != 0, "Factor of 0 is unphysical!", file);
				factor = pow(factor, -0.5);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					if (debug)
					{
						std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << "\n"
							<< "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << "\n";
					}
					basis_coefficients[a][i] *= factor;
					norm_const.push_back(basis_coefficients[a][i]);
				}
				break;
			case 2:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						double aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						double term = 4 * pow(aiaj, 5);
						term = constants::PI3 / term;
						term = pow(term, 0.5);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
					}
				}
				err_checkf(factor != 0, "Factor of 0 is unphysical!", file);
				factor = pow(factor, -0.5);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					if (debug)
					{
						std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << "\n"
							<< "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << "\n";
					}
					basis_coefficients[a][i] *= factor;
					for (int k = 0; k < 3; k++)
						norm_const.push_back(basis_coefficients[a][i]);
				}
				break;
			case 3:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						double aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						double term = 16 * pow(aiaj, 7);
						term = constants::PI3 / term;
						term = pow(term, 0.5);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
					}
				}
				err_checkf(factor != 0, "Factor of 0 is unphysical!", file);
				factor = (pow(factor, -0.5)) / sqrt(3);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					if (debug)
					{
						std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << "\n"
							<< "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << "\n";
					}
					basis_coefficients[a][i] *= factor;
					for (int k = 0; k < 3; k++)
						norm_const.push_back(basis_coefficients[a][i]);
					for (int k = 0; k < 3; k++)
						norm_const.push_back(sqrt(3) * basis_coefficients[a][i]);
				}
				break;
			case 4:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						double aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						double term = 64 * pow((aiaj), 9);
						term = constants::PI3 / term;
						term = pow(term, 0.5);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
					}
				}
				err_checkf(factor != 0, "Factor of 0 is unphysical!", file);
				factor = pow(factor, -0.5) / sqrt(15);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					if (debug)
					{
						std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << "\n"
							<< "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << "\n";
					}
					basis_coefficients[a][i] *= factor;
					for (int l = 0; l < 3; l++)
						norm_const.push_back(basis_coefficients[a][i]);
					for (int l = 0; l < 6; l++)
						norm_const.push_back(sqrt(5) * basis_coefficients[a][i]);
					norm_const.push_back(sqrt(15) * basis_coefficients[a][i]);
				}
				break;
			}
			if (debug)
				std::cout << "This shell has: " << get_shell_end(a, s) - get_shell_start(a, s) + 1 << " primitives\n";
		}
	}
	return norm_const;
}

bool WFN::write_wfn(const std::filesystem::path &fileName, const bool &debug, const bool occupied) const
{
	using namespace std;
	ofstream rf(fileName, ios::out);
	if (!rf.is_open())
	{
		std::cout << "Sorry, can't open the file...\n";
		return false;
	}
	if (debug)
		std::cout << "Writing " << fileName << ": ncen " << ncen << " nex " << nex << " nmo " << nmo << endl;
	rf << comment << '\n' << hdr(occupied);
	for (int i = 0; i < ncen; i++)
	{
		rf << setw(3) << atoms[i].get_label() << ' ';
		rf << left << setw(8) << i + 1 << right << "(CENTRE";
		rf << setw(3) << i + 1 << ") ";
		rf << fixed << showpoint << setprecision(8);
		rf << setw(12) << get_atom_coordinate(i, 0);
		rf << setw(12) << get_atom_coordinate(i, 1);
		rf << setw(12) << get_atom_coordinate(i, 2);
		rf << "  CHARGE = " << setw(2) << get_atom_charge(i) << ".0\n";
	}
	//Fortran fixed-width blocks: every line starts with tag and holds per values, value(i) prints the i-th
	char buf[32];
	auto write_block = [&](const string &tag, const int per, auto value)
	{
		for (int i = 0; i < nex; i++)
		{
			if (i % per == 0)
				rf << (i ? "\n" : "") << tag;
			value(i);
		}
		if (nex > 0)
			rf << '\n';
	};
	write_block("CENTRE ASSIGNMENTS  ", 20, [&](const int i) { rf << setw(3) << centers[i]; });
	write_block("TYPE ASSIGNMENTS    ", 20, [&](const int i) { rf << setw(3) << types[i]; });
	write_block("EXPONENTS ", 5, [&](const int i) { snprintf(buf, sizeof(buf), "%14.7E", exponents[i]); rf << buf; });
	int mo_run = 1;
	for (int m = 0; m < nmo; m++)
	{
		if (occupied && MOs[m].get_occ() == 0)
			continue;
		rf << "MO" << setw(3) << mo_run++ << setw(29) << "OCC NO =" << setw(13) << fixed << setprecision(8) << MOs[m].get_occ()
		   << setw(14) << "ORB. ENERGY =" << setw(13) << fixed << setprecision(8) << MOs[m].get_energy() << '\n';
		write_block("", 5, [&](const int i) { snprintf(buf, sizeof(buf), "%16.8E", MOs[m].get_coefficient(i)); rf << buf; });
	}
	rf << "END DATA\n";
	rf << " THE SCF ENERGY =" << setw(20) << fixed << setprecision(12) << total_energy << " THE VIRIAL(-V/T)=   0.00000000" << endl;
	return rf.good();
};

//AIM wfx: every quantity in its own <Tag> ... </Tag> block, the ORCA/AIMAll layout that read_wfx expects
bool WFN::write_wfx(const std::filesystem::path &fileName, const bool occupied) const
{
	using namespace std;
	ofstream rf(fileName, ios::out);
	if (!rf.is_open())
	{
		std::cout << "Sorry, can't open the file...\n";
		return false;
	}
	char buf[32];
	auto block = [&](const string &tag, auto body) { rf << '<' << tag << ">\n"; body(); rf << "</" << tag << ">\n"; };
	//per values per line; value(i) prints the i-th of n
	auto numbers = [&](const string &tag, const int n, const int per, auto value)
	{
		block(tag, [&]() { for (int i = 0; i < n; i++) { value(i); rf << ((i + 1) % per == 0 || i + 1 == n ? "\n" : " "); } });
	};
	auto sci = [&](const double x) { snprintf(buf, sizeof(buf), "%16.8E", x); rf << buf; };
	ivec sel;
	double nel = 0, nalpha = 0;
	for (int m = 0; m < nmo; m++)
	{
		if (occupied && MOs[m].get_occ() == 0)
			continue;
		sel.push_back(m);
		nel += MOs[m].get_occ();
		nalpha += is_unrestricted ? (MOs[m].get_op() == 0 ? MOs[m].get_occ() : 0) : MOs[m].get_occ() / 2;
	}
	const int n = static_cast<int>(sel.size()), i_nel = static_cast<int>(round(nel)), i_nalpha = static_cast<int>(round(nalpha));
	block("Title", [&]() { rf << comment << '\n'; });
	block("Keywords", [&]() { rf << "GTO\n"; });
	block("Number of Nuclei", [&]() { rf << ncen << '\n'; });
	block("Number of Primitives", [&]() { rf << nex << '\n'; });
	block("Number of Occupied Molecular Orbitals", [&]() { rf << n << '\n'; });
	block("Number of Perturbations", [&]() { rf << "0\n"; });
	numbers("Nuclear Names", ncen, 1, [&](const int i) { rf << (atoms[i].get_label().empty() ? constants::atnr2letter(get_atom_charge(i)) + to_string(i + 1) : atoms[i].get_label()); });
	numbers("Atomic Numbers", ncen, 1, [&](const int i) { rf << get_atom_charge(i); });
	numbers("Nuclear Charges", ncen, 1, [&](const int i) { sci(get_atom_charge(i)); });
	numbers("Nuclear Cartesian Coordinates", 3 * ncen, 3, [&](const int i) { sci(get_atom_coordinate(i / 3, i % 3)); });
	block("Net Charge", [&]() { rf << charge << '\n'; });
	block("Number of Electrons", [&]() { rf << i_nel << '\n'; });
	block("Number of Alpha Electrons", [&]() { rf << i_nalpha << '\n'; });
	block("Number of Beta Electrons", [&]() { rf << i_nel - i_nalpha << '\n'; });
	block("Electronic Spin Multiplicity", [&]() { rf << (multi > 0 ? static_cast<int>(multi) : 2 * i_nalpha - i_nel + 1) << '\n'; });
	numbers("Primitive Centers", nex, 20, [&](const int i) { rf << centers[i]; });
	numbers("Primitive Types", nex, 20, [&](const int i) { rf << types[i]; });
	numbers("Primitive Exponents", nex, 5, [&](const int i) { sci(exponents[i]); });
	numbers("Molecular Orbital Occupation Numbers", n, 1, [&](const int i) { sci(MOs[sel[i]].get_occ()); });
	numbers("Molecular Orbital Energies", n, 1, [&](const int i) { sci(MOs[sel[i]].get_energy()); });
	numbers("Molecular Orbital Spin Types", n, 1, [&](const int i) { rf << (!is_unrestricted ? "Alpha and Beta" : MOs[sel[i]].get_op() == 0 ? "Alpha" : "Beta"); });
	block("Molecular Orbital Primitive Coefficients", [&]()
	{
		for (int i = 0; i < n; i++)
		{
			block("MO Number", [&]() { rf << i + 1 << '\n'; });
			for (int j = 0; j < nex; j++)
			{
				sci(MOs[sel[i]].get_coefficient(j));
				rf << ((j + 1) % 5 == 0 || j + 1 == nex ? "\n" : " ");
			}
		}
	});
	block("Energy = T + Vne + Vee + Vnn", [&]() { snprintf(buf, sizeof(buf), "%22.14E", total_energy); rf << buf << '\n'; });
	block("Virial Ratio (-V/T)", [&]() { sci(virial_ratio); rf << '\n'; });
	return rf.good();
};

bool WFN::write_nbo(const std::filesystem::path &fileName, const bool &debug, std::ostream* progress_log, const std::string& nbo_keywords)
{
	using namespace std;

	//A FILE47 needs the contracted shell structure ($BASIS/$CONTRACT) and an AO overlap
	//computed over it. A .wfn/.wfx carries primitives only - the shells, their contraction
	//coefficients and the primitive-to-shell ordering are all gone - so no archive can be
	//built from one without guessing, and a guessed archive produces plausible-looking but
	//wrong NBO output. Convert through .molden/.fchk/.gbw instead.
	err_checkf(get_nr_basis_set_loaded() == ncen,
		"Can only write a .47 file when a contracted basis set is present. A primitive-only source"
		" (.wfn/.wfx) does not carry one - use the .gbw, .fchk or .molden of the same calculation.",
		std::cout);
	const auto nbo_start_time = std::chrono::high_resolution_clock::now();
	auto progress_elapsed_seconds = [&]() {
		return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - nbo_start_time).count();
	};
	auto progress = [&](const std::string& message) {
		if (progress_log != nullptr) {
			*progress_log << "[FILE47] " << message << " (" << progress_elapsed_seconds() << " s)" << std::endl;
			progress_log->flush();
		}
	};
	auto progress_percent = [&](const std::string& label, const int completed, const int total, int& next_percent) {
		if (progress_log == nullptr || total <= 0)
			return;
		const int percent = static_cast<int>((100LL * completed) / total);
		if (percent >= next_percent || completed == total) {
			*progress_log << "[FILE47] " << label << " " << percent << "% (" << completed << "/" << total << ", "
						  << progress_elapsed_seconds() << " s)" << std::endl;
			progress_log->flush();
			next_percent = percent + 10;
		}
	};

	progress("Starting .47 conversion to " + fileName.string());

	struct NboShell {
		int atom = 0;
		int shell = 0;
		int type = 0;
		int nprim = 0;
		int atom_prim_start = 0;
		int internal_start = 0;
		int nbo_start = 0;
		int cart_components = 0;
		int nbo_components = 0;
	};


	auto nbo_labels = [](const int type) -> ivec {
		switch (type) {
		case 1: return { 1 };
		case 2: return { 103, 101, 102 };
		case 3: return { 255, 252, 253, 254, 251 };
		case 4: return { 351, 352, 353, 354, 355, 356, 357 };
		case 5: return { 451, 452, 453, 454, 455, 456, 457, 458, 459 };
		default: return {};
		}
	};

	auto sph2cart_coefficient = [](const int type, const int cart, const int spher) {
		return constants::sph2cart(type - 1)[cart * constants::n_spher(type - 1) + spher];
	};

	auto project_to_nbo = [&](const int type, const vec& cart_values) {
		const int cart_count = constants::n_cart(type - 1);
		const int nbo_count = constants::n_spher(type - 1);
		vec2 normal(nbo_count, vec(nbo_count, 0.0));
		vec rhs(nbo_count, 0.0);
		for (int c = 0; c < cart_count; c++) {
			for (int i = 0; i < nbo_count; i++) {
				const double transform = sph2cart_coefficient(type, c, i);
				rhs[i] += transform * cart_values[c];
				for (int j = 0; j < nbo_count; j++)
					normal[i][j] += transform * sph2cart_coefficient(type, c, j);
			}
		}
		solve_linear_system(normal, rhs);
		return rhs;
	};

	vector<NboShell> shells;
	int nbo_nao = 0;
	int nbo_nexp = 0;
	int highest_angular = -1;
	for (int a = 0; a < get_ncen(); a++) {
		for (int s = 0; s < get_atom_shell_count(a); s++) {
			const int type = get_shell_type(a, s);
			const int cart_count = constants::n_cart(type - 1);
			const int nbo_count = constants::n_spher(type - 1);
			//FILE47 itself goes further than g: it has label codes for h (Cartesian 501-521,
			//spherical 551-563) and i (601-628 / 651-665) and $CONTRACT arrays CH and CI. The
			//ceiling here is NoSpherA2's, not the archive's - constants::n_cart / n_spher and
			//constants::sph2cart stop at g, and no basis used with NoSpherA2 (def2, cc-pVnZ up
			//to quadruple zeta, jorge, x2c) carries an h shell. Add CH/CI here if one ever does.
			err_checkf(type <= 5, "Unsupported basis shell in .47 writer: shells beyond g need"
				" constants::sph2cart extended first (FILE47 itself supports h and i)", std::cout);
			NboShell shell;
			shell.atom = a;
			shell.shell = s;
			shell.type = type;
			shell.nprim = get_atom_shell_primitives(a, s);
			shell.atom_prim_start = get_shell_start(a, s);
			shell.internal_start = 0;
			shell.nbo_start = nbo_nao;
			shell.cart_components = cart_count;
			shell.nbo_components = nbo_count;
			shells.push_back(shell);
			nbo_nao += nbo_count;
			nbo_nexp += shell.nprim;
			highest_angular = std::max(highest_angular, type - 1);
		}
	}
	//Int_Params, the gbw reader and OCC group an atom's shells by angular momentum, the loop
	//above follows the wavefunction's shell order. internal_start indexes the spherical
	//overlap, density and coefficient matrices in the grouped order, nbo_start the .47 order.
	{
		int internal_run = 0;
		for (int a = 0; a < get_ncen(); a++) {
			int max_l = 0;
			for (auto& sh : shells)
				if (sh.atom == a) max_l = std::max(max_l, sh.type - 1);
			for (int l = 0; l <= max_l; l++)
				for (auto& sh : shells)
					if (sh.atom == a && sh.type - 1 == l) {
						sh.internal_start = internal_run;
						internal_run += sh.nbo_components;
					}
		}
		err_checkf(internal_run == nbo_nao, "Angular-momentum regrouping lost basis functions in the .47 writer", std::cout);
	}
	progress("Built NBO shell model: atoms=" + std::to_string(get_ncen()) +
		", shells=" + std::to_string(shells.size()) +
		", NBO basis functions=" + std::to_string(nbo_nao) +
		", primitives=" + std::to_string(nbo_nexp));

	progress("Normalizing basis coefficients");
	vec2 basis_coefficients(get_ncen());
#pragma omp parallel for
	for (int a = 0; a < get_ncen(); a++)
	{
		for (int p = 0; p < get_atom_primitive_count(a); p++)
		{
			double temp_c = get_atom_basis_set_exponent(a, p);
			switch (get_atom_primitive_type(a, p))
			{
			case 1:
				temp_c = 8 * pow(temp_c, 3) / constants::PI3;
				break;
			case 2:
				temp_c = 128 * pow(temp_c, 5) / constants::PI3;
				break;
			case 3:
				temp_c = 2048 * pow(temp_c, 7) / (9 * constants::PI3);
				break;
			case 4:
				temp_c = 32768 * pow(temp_c, 9) / (225 * constants::PI3);
				break;
			case -1:
				std::cout << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy...\n";
				break;
			}
			temp_c = pow(temp_c, 0.25) * get_atom_basis_set_coefficient(a, p);
			basis_coefficients[a].push_back(temp_c);
		}
	}
	for (int a = 0; a < get_ncen(); a++)
	{
		double aiaj = 0.0;
		double factor = 0.0;
		for (int s = 0; s < get_atom_shell_count(a); s++)
		{
			int type_temp = get_shell_type(a, s);
			err_checkf(type_temp != -1, "ERROR in type assignement!!", std::cout);
			if (debug)
			{
				std::cout << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << "\n"
					<< "start: " << get_shell_start(a, s)
					<< " stop: " << get_shell_end(a, s) << "\n"
					<< "factor: ";
			}
			switch (type_temp)
			{
			case 1:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * pow(constants::PI3 / pow(aiaj, 3), 0.5);
					}
				}
				if (factor == 0)
					return false;
				factor = pow(factor, -0.5);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					basis_coefficients[a][i] *= factor;
				}
				break;
			case 2:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * pow(constants::PI3 / (4 * pow(aiaj, 5)), 0.5);
					}
				}
				if (factor == 0)
					return false;
				factor = pow(factor, -0.5);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					basis_coefficients[a][i] *= factor;
				}
				break;
			case 3:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * pow(constants::PI3 / (16 * pow(aiaj, 7)), 0.5);
					}
				}
				if (factor == 0)
					return false;
				factor = (pow(factor, -0.5)) / sqrt(3);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					basis_coefficients[a][i] *= factor;
				}
				break;
			case 4:
				factor = 0;
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					for (int j = get_shell_start(a, s); j <= get_shell_end(a, s); j++)
					{
						aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
						factor += basis_coefficients[a][i] * basis_coefficients[a][j] * pow(constants::PI3 / (64 * pow((aiaj), 9)), 0.5);
					}
				}
				if (factor == 0)
					return false;
				factor = pow(factor, -0.5) / sqrt(15);
				if (debug)
					std::cout << factor << "\n";
				for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
				{
					basis_coefficients[a][i] *= factor;
				}
				break;
			}
			if (debug)
				std::cout << "This shell has: " << get_shell_end(a, s) - get_shell_start(a, s) + 1 << " primitives\n";
		}
	}
	//NBO slot c of a shell (m = 0, +1, -1, +2, -2, ... as ORCA orders it) holds pure function
	//orca_2_pySCF(l, c) of the m = -l..l order libcint, OCC and the gbw reader keep their
	//matrices in; sph_index maps a slot to its row in those matrices.
	ivec sph_index(nbo_nao, 0);
	for (const auto& shell : shells)
		for (int c = 0; c < shell.nbo_components; c++) {
			const auto offset = constants::orca_2_pySCF(shell.type - 1, c);
			err_checkf(offset.has_value(), "Unsupported NBO AO order in .47 writer", std::cout);
			sph_index[shell.nbo_start + c] = shell.internal_start + static_cast<int>(offset.value());
		}

	const int alpha_mos = get_MO_op_count(0);
	const int beta_mos = get_MO_op_count(1);
	vec2 CMO(alpha_mos, vec(nbo_nao, 0.0));
	vec2 CMO_beta(beta_mos, vec(nbo_nao, 0.0));
	//An XCW_fit wavefunction still holds the spherical coefficients OCC converged; rebuilding
	//them from the primitives below needs OCC's normalization convention, which the WFN does
	//not carry. Other origins that cache their coefficients could take this path as well.
	const int spin_blocks = get_is_unrestricted() ? 2 : 1;
	const bool cached_mos = get_origin() == e_origin::XCW_fit
		&& static_cast<int>(MO_sph.extent(0)) == spin_blocks * nbo_nao
		&& static_cast<int>(MO_sph.extent(1)) >= std::max(alpha_mos, beta_mos);
	if (cached_mos) {
		progress("Taking MO coefficients from the cached spherical coefficient matrix");
		for (int m = 0; m < alpha_mos; m++)
			for (int i = 0; i < nbo_nao; i++)
				CMO[m][i] = MO_sph(sph_index[i], m);
		for (int m = 0; m < beta_mos; m++)
			for (int i = 0; i < nbo_nao; i++)
				CMO_beta[m][i] = MO_sph(nbo_nao + sph_index[i], m);
	}
	else
		progress("Reconstructing MO coefficients in NBO AO order: alpha=" + std::to_string(alpha_mos) +
			", beta=" + std::to_string(beta_mos));
	int alpha_run = 0;
	int beta_run = 0;
	int mo_progress_next = 10;
	int mo_seen = 0;
	for (int m = 0; m < get_nmo() && !cached_mos; m++)
	{
		if (MOs[m].get_op() != 0 && MOs[m].get_op() != 1)
			continue;
		mo_seen++;
		vec* target = nullptr;
		if (MOs[m].get_op() == 0) {
			target = &CMO[alpha_run++];
		}
		else {
			target = &CMO_beta[beta_run++];
		}
		for (const auto& shell : shells) {
			vec cart_values(shell.cart_components, 0.0);
			const int primitive_start = get_shell_start_in_primitives(shell.atom, shell.shell);
			//Every primitive of a contracted shell is AO coefficient times its stored contraction
			//coefficient, so the largest one is the safest divisor. The gbw reader stores a
			//shell primitive-major (x, y, z of primitive 0, then of primitive 1), the OCC
			//constructor component-major (all primitives of x, then of y); equal consecutive
			//types tell the two apart.
			const int shell_start = get_shell_start(shell.atom, shell.shell);
			int rep = 0;
			for (int p = 1; p < shell.nprim; p++)
				if (std::abs(get_atom_basis_set_coefficient(shell.atom, shell_start + p)) >
					std::abs(get_atom_basis_set_coefficient(shell.atom, shell_start + rep)))
					rep = p;
			const double contraction = get_atom_basis_set_coefficient(shell.atom, shell_start + rep);
			err_checkf(std::abs(contraction) > 1E-14, "Cannot reconstruct pure MO coefficients for .47 output", std::cout);
			const bool component_major = shell.nprim > 1 && shell.cart_components > 1
				&& get_type(primitive_start) == get_type(primitive_start + 1);
			const int stride = component_major ? shell.nprim : 1;
			const int rep_offset = component_major ? rep : rep * shell.cart_components;
			//the components may be permuted as well (gbw stores p as z, x, y): the type says which row
			for (int c = 0; c < shell.cart_components; c++) {
				const int prim = primitive_start + c * stride + rep_offset;
				cart_values[get_type(prim) - constants::first_type[shell.type - 1]] = get_MO_coef(m, prim) / contraction;
			}
			const vec nbo_values = project_to_nbo(shell.type, cart_values);
			for (int c = 0; c < shell.nbo_components; c++)
				(*target)[shell.nbo_start + c] = nbo_values[c];
		}
		progress_percent("MO coefficient reconstruction", mo_seen, alpha_mos + beta_mos, mo_progress_next);
	}
	progress("Building density matrix");
	int naotr = nbo_nao * (nbo_nao + 1) / 2;
	vec CDM(naotr, 0.0);
	//Occupations and energies per spin, in the row order of CMO / CMO_beta. An open-shell
	//FILE47 carries $DENSITY, $FOCK and $LCAOMO twice - alpha block then beta block - while
	//$OVERLAP stays single; a closed-shell one carries the spin sum once.
	vec occ_spin[2], energy_spin[2];
	for (int m = 0; m < get_nmo(); m++) {
		const int op = MOs[m].get_op();
		if (op != 0 && op != 1)
			continue;
		occ_spin[op].push_back(get_MO_occ(m));
		energy_spin[op].push_back(get_MO_energy(m));
	}
	const bool open_shell = get_is_unrestricted() && beta_mos > 0;
	auto build_mo_density = [&](const int op) {
		const vec2& C = op == 0 ? CMO : CMO_beta;
		const vec& occs = occ_spin[op];
		const int nmo_spin = std::min(static_cast<int>(C.size()), static_cast<int>(occs.size()));
		vec density(naotr, 0.0);
		int density_progress_next = 10;
#pragma omp parallel for schedule(dynamic)
		for (int iu = 0; iu < nbo_nao; iu++) {
			for (int iv = 0; iv <= iu; iv++) {
				const int iuv = (iu * (iu + 1) / 2) + iv;
				for (int m = 0; m < nmo_spin; m++)
					if (occs[m] != 0.0)
						density[iuv] += occs[m] * C[m][iu] * C[m][iv];
			}
#pragma omp critical(nbo_progress)
			progress_percent("Density build", iu + 1, nbo_nao, density_progress_next);
		}
		return density;
	};
	auto build_total_density = [&]() {
		vec density = build_mo_density(0);
		if (beta_mos > 0) {
			const vec beta_density = build_mo_density(1);
			for (int i = 0; i < naotr; i++)
				density[i] += beta_density[i];
		}
		return density;
	};
	vec CDM_alpha, CDM_beta;
	bool density_from_cached_dm = false;
	if (open_shell) {
		progress("Building spin-resolved density matrices for the open-shell FILE47");
		CDM_alpha = build_mo_density(0);
		CDM_beta = build_mo_density(1);
		for (int i = 0; i < naotr; i++)
			CDM[i] = CDM_alpha[i] + CDM_beta[i];
	}
	//A cached DM is the spin sum, so it can only serve the closed-shell layout.
	else if (static_cast<int>(DM.extent(0)) == nbo_nao && static_cast<int>(DM.extent(1)) == nbo_nao) {
		density_from_cached_dm = true;
		for (int iu = 0; iu < nbo_nao; iu++) {
			for (int iv = 0; iv <= iu; iv++) {
				const int iuv = (iu * (iu + 1) / 2) + iv;
				CDM[iuv] = DM(sph_index[iu], sph_index[iv]);
			}
		}
	}
	else {
		CDM = build_total_density();
	}

	vec OVLP_matrix = {};
	//Int_Params reads each shell's angular momentum by origin (OCC and NOT_YET_DEFINED
	//store l, everything else l + 1) and normalises by origin; an XCW_fit wavefunction
	//says which convention it has, so it is handed over as it is.
	Int_Params int_params(*this);
	progress("Computing spherical AO overlap integrals");
	compute2C<Overlap2C_SPH>(int_params, OVLP_matrix);
	err_checkf(static_cast<int>(OVLP_matrix.size()) == nbo_nao * nbo_nao, "Spherical overlap has the wrong size in the .47 writer", std::cout);
	dMatrixRef2 OVLP_sph(OVLP_matrix.data(), nbo_nao, nbo_nao);
	//libcint and OCC use the standard phases; ORCA's f(+-3), g(+-3), g(+-4) have the opposite
	//sign and the gbw reader keeps them, so the overlap takes ORCA's sign there. Other origins
	//with ORCA-like conventions may need the same and are not checked.
	vec phase(nbo_nao, 1.0);
	if (get_origin() == e_origin::gbw)
		for (const auto& shell : shells)
			for (int c = 5; c < shell.nbo_components; c++)
				phase[shell.nbo_start + c] = -1.0;
	vec2 OVLP_nbo(nbo_nao, vec(nbo_nao, 0.0));
	for (int i = 0; i < nbo_nao; i++)
		for (int j = 0; j < nbo_nao; j++)
			OVLP_nbo[i][j] = phase[i] * phase[j] * OVLP_sph(sph_index[i], sph_index[j]);
	//A diagonal off 1 means Int_Params normalised this origin with the wrong convention; the
	//density beside it is in a normalised basis, so rescaling S alone would only hide it.
	int unnormalised = 0;
	for (int i = 0; i < nbo_nao; i++)
		if (std::abs(OVLP_nbo[i][i] - 1.0) > 1e-8) unnormalised++;
	if (unnormalised > 0)
		progress(std::to_string(unnormalised) + " of " + std::to_string(nbo_nao)
			+ " AOs are not normalised in the .47 overlap; check the origin's normalisation in Int_Params");
	auto packed_trace_product = [&](const vec& density) {
		double trace = 0.0;
		for (int i = 0; i < nbo_nao; i++) {
			for (int j = 0; j <= i; j++) {
				const int ij = (i * (i + 1) / 2) + j;
				trace += density[ij] * OVLP_nbo[i][j] * (i == j ? 1.0 : 2.0);
			}
		}
		return trace;
	};
	double expected_electrons = 0.0;
	for (int m = 0; m < get_nmo(); m++)
		expected_electrons += get_MO_occ(m);
	double density_electrons = packed_trace_product(CDM);
	if (density_from_cached_dm && std::abs(density_electrons - expected_electrons) > 1.0E-4) {
		progress("Cached GBW density is inconsistent with FILE47 overlap: Tr(P*S)=" +
			std::to_string(density_electrons) + ", expected=" + std::to_string(expected_electrons) +
			". Rebuilding density from NBO-ordered MO coefficients");
		CDM = build_total_density();
		density_from_cached_dm = false;
		density_electrons = packed_trace_product(CDM);
	}
	if (progress_log != nullptr) {
		std::ostringstream density_check;
		density_check << std::fixed << std::setprecision(6) << density_electrons
					  << ", expected=" << expected_electrons;
		*progress_log << "[FILE47] Density electron check Tr(P*S)=" << density_check.str()
					  << " (" << progress_elapsed_seconds() << " s)" << std::endl;
		progress_log->flush();
	}
	auto build_fock = [&](const vec2& C, const vec& energies, const std::string& label) {
		vec2 result;
		if (static_cast<int>(C.size()) != nbo_nao || static_cast<int>(energies.size()) < nbo_nao) {
			progress("Skipping " + label + " Fock matrix: MO count does not match NBO basis size");
			return result;
		}
		progress("Building " + label + " Fock matrix from MO energies with BLAS");
		result = vec2(nbo_nao, vec(nbo_nao, 0.0));
		dMatrix2 overlap(nbo_nao, nbo_nao);
		dMatrix2 cmo(nbo_nao, nbo_nao);
		for (int i = 0; i < nbo_nao; i++) {
			for (int j = 0; j < nbo_nao; j++)
				overlap(i, j) = OVLP_nbo[i][j];
			for (int m = 0; m < nbo_nao; m++)
				cmo(m, i) = C[m][i];
		}
		dMatrix2 eps_cmo(nbo_nao, nbo_nao);
#pragma omp parallel for schedule(dynamic)
		for (int m = 0; m < nbo_nao; m++)
			for (int i = 0; i < nbo_nao; i++)
				eps_cmo(m, i) = energies[m] * cmo(m, i);
		progress("Fock build: C^T * eps * C");
		dMatrix2 ctc = dot<dMatrix2>(cmo, eps_cmo, true, false);
		progress("Fock build: S * (C^T * eps * C)");
		dMatrix2 left = dot<dMatrix2>(overlap, ctc, false, false);
		progress("Fock build: S * (C^T * eps * C) * S");
		dMatrix2 fock = dot<dMatrix2>(left, overlap, false, false);
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nbo_nao; i++)
			for (int j = 0; j < nbo_nao; j++)
				result[i][j] = fock(i, j);
		return result;
	};
	vec2 FOCK_nbo = build_fock(CMO, energy_spin[0], open_shell ? "alpha" : "total");
	vec2 FOCK_beta;
	if (open_shell) {
		FOCK_beta = build_fock(CMO_beta, energy_spin[1], "beta");
		//An open-shell $FOCK is read as two blocks; one alone would be parsed as the alpha
		//block and leave NBO reading the next section as beta, so it is both or neither.
		if (FOCK_nbo.empty() || FOCK_beta.empty()) {
			FOCK_nbo.clear();
			FOCK_beta.clear();
			progress("Skipping $FOCK entirely: an open-shell FILE47 needs both spin blocks");
		}
	}

	ofstream rf(fileName, ios::out);
	if (!rf.is_open())
	{
		std::cout << "Sorry, can't open the file...\n";
		return false;
	}
	progress("Writing FILE47 sections");

	auto write_int_array = [&](const string& name, const ivec& values, const int per_line, const int continuation_indent) {
		rf << name;
		for (int i = 0; i < values.size(); i++) {
			if (i > 0 && i % per_line == 0)
				rf << "\n" << string(continuation_indent, ' ');
			rf << setw(5) << values[i];
		}
		rf << endl;
	};

	auto format_precise_value = [](const double value) {
		if (value == 0.0)
			return string("  0.000000000000E+00");

		const double abs_value = std::abs(value);
		const int exponent = static_cast<int>(std::floor(std::log10(abs_value))) + 1;
		double mantissa = value / std::pow(10.0, exponent);
		int adjusted_exponent = exponent;
		if (std::abs(mantissa) >= 1.0) {
			mantissa /= 10.0;
			adjusted_exponent++;
		}

		stringstream local;
		local << fixed << setprecision(12) << mantissa;
		if (local.str().rfind("1.", 0) == 0 || local.str().rfind("-1.", 0) == 0) {
			mantissa /= 10.0;
			adjusted_exponent++;
			local.str("");
			local.clear();
			local << fixed << setprecision(12) << mantissa;
		}
		string result = local.str() + "E" + (adjusted_exponent >= 0 ? "+" : "-");
		result += (std::abs(adjusted_exponent) < 10 ? "0" : "") + std::to_string(std::abs(adjusted_exponent));
		return string(std::max(0, 20 - static_cast<int>(result.size())), ' ') + result;
	};

	auto write_real_array = [&](const string& name, const vec& values) {
		rf << name;
		for (int i = 0; i < values.size(); i++) {
			if (i > 0 && i % 3 == 0)
				rf << "\n" << string(10, ' ');
			rf << format_precise_value(values[i]);
		}
		rf << endl;
	};

	rf << " $GENNBO NATOMS=" << ncen << " NBAS=" << nbo_nao << (open_shell ? " OPEN" : "") << " UPPER BODM FORMAT=PRECISE $END" << endl;
	rf << " $NBO" << (nbo_keywords.empty() ? "" : " " + nbo_keywords) << " $END" << endl;
	rf << " $COORD" << endl;
	rf << " .47 file generated by NoSpherA2 based on " << path << endl;
	for (int i = 0; i < ncen; i++)
	{
		const int atomic_number = get_atom_charge(i);
		const int nuclear_charge = atomic_number - get_atom_ECP_electrons(i);
		rf << " " << setw(5) << atomic_number << setw(5) << nuclear_charge
			<< fixed << setprecision(6) << setw(15) << constants::bohr2ang(get_atom_coordinate(i, 0))
			<< fixed << setprecision(6) << setw(15) << constants::bohr2ang(get_atom_coordinate(i, 1))
			<< fixed << setprecision(6) << setw(15) << constants::bohr2ang(get_atom_coordinate(i, 2)) << "\n";
	}
	rf << " $END" << endl;

	ivec nbo_centers;
	ivec labels;
	for (const auto& shell : shells) {
		const ivec shell_labels = nbo_labels(shell.type);
		for (int label : shell_labels) {
			nbo_centers.push_back(shell.atom + 1);
			labels.push_back(label);
		}
	}
	rf << " $BASIS" << endl;
	write_int_array("  CENTER =", nbo_centers, 11, 10);
	write_int_array("   LABEL =", labels, 11, 10);
	rf << " $END" << endl;

	ivec ncomp;
	ivec nprim;
	ivec nptr;
	vec exp_values;
	vec2 contract_values(5, vec(nbo_nexp, 0.0));
	int exp_index = 0;
	for (const auto& shell : shells) {
		ncomp.push_back(shell.nbo_components);
		nprim.push_back(shell.nprim);
		nptr.push_back(exp_index + 1);
		for (int p = 0; p < shell.nprim; p++) {
			exp_values.push_back(get_atom_basis_set_exponent(shell.atom, shell.atom_prim_start + p));
			contract_values[shell.type - 1][exp_index] = get_atom_basis_set_coefficient(shell.atom, shell.atom_prim_start + p);
			exp_index++;
		}
	}

	rf << " $CONTRACT" << endl;
	rf << "  NSHELL =" << setw(6) << shells.size() << endl;
	rf << "    NEXP =" << setw(6) << nbo_nexp << endl;
	write_int_array("   NCOMP =", ncomp, 11, 10);
	write_int_array("   NPRIM =", nprim, 11, 10);
	write_int_array("    NPTR =", nptr, 11, 10);
	write_real_array("     EXP =", exp_values);
	if (highest_angular >= 0) write_real_array("      CS =", contract_values[0]);
	if (highest_angular >= 1) write_real_array("      CP =", contract_values[1]);
	if (highest_angular >= 2) write_real_array("      CD =", contract_values[2]);
	if (highest_angular >= 3) write_real_array("      CF =", contract_values[3]);
	if (highest_angular >= 4) write_real_array("      CG =", contract_values[4]);
	rf << " $END" << endl;
	auto write_precise_value = [&](const double value, int& count) {
		rf << format_precise_value(value);
		count++;
		if (count % 4 == 0)
			rf << "\n";
	};
	//NBO reads each spin block with its own Fortran READ, so the beta block has to start on a
	//fresh record. Streaming both blocks as one run of values makes TINP report
	//"error reading $LCAOMO" whenever a block length is not a multiple of four.
	auto end_block = [&](int& count) {
		if (count % 4 != 0)
			rf << "\n";
		count = 0;
	};

	rf << " $OVERLAP" << endl;
	int runner = 0;
	for (int i = 0; i < nbo_nao; i++)
		for (int j = 0; j <= i; j++)
			write_precise_value(OVLP_nbo[i][j], runner);
	end_block(runner);
	rf << " $END" << endl;
	rf << " $DENSITY" << endl;
	for (int i = 0; i < naotr; i++)
		write_precise_value(open_shell ? CDM_alpha[i] : CDM[i], runner);
	end_block(runner);
	if (open_shell) {
		for (int i = 0; i < naotr; i++)
			write_precise_value(CDM_beta[i], runner);
		end_block(runner);
	}
	rf << " $END" << endl;
	if (!FOCK_nbo.empty()) {
		auto write_fock_block = [&](const vec2& fock) {
			for (int i = 0; i < nbo_nao; i++)
				for (int j = 0; j <= i; j++)
					write_precise_value(fock[i][j], runner);
			end_block(runner);
		};
		rf << " $FOCK" << endl;
		write_fock_block(FOCK_nbo);
		if (open_shell)
			write_fock_block(FOCK_beta);
		rf << " $END" << endl;
	}
	//$LCAOMO is nbas x nbas per spin block whatever the MO count, so a wavefunction that
	//carries fewer MOs than basis functions is zero-padded - a short block would otherwise
	//shift every value after it (and, open shell, the whole beta block).
	auto write_lcaomo_block = [&](const vec2& C) {
		for (int mo_counter = 0; mo_counter < nbo_nao; mo_counter++)
		{
			if (debug)
				std::cout << "Writing MO #" << mo_counter + 1 << "...\n";
			for (int i = 0; i < nbo_nao; i++)
				write_precise_value(mo_counter < static_cast<int>(C.size()) ? C[mo_counter][i] : 0.0, runner);
		}
		end_block(runner);
	};
	rf << " $LCAOMO" << endl;
	write_lcaomo_block(CMO);
	if (open_shell)
		write_lcaomo_block(CMO_beta);
	rf << " $END" << endl;
	rf.close();
	progress("Finished .47 conversion");
	return true;
};

bool WFN::write_xyz(const std::filesystem::path &fileName)
{
	using namespace std;
	ofstream f(fileName, ios::out);
	if (!f.is_open())
	{
		err("Error writing the xyz file! Aborting!", std::cout);
		return false;
	}
	f << ncen << '\n' << "XYZ File written by NoSpherA2 based on " << path << '\n';
	for (int i = 0; i < ncen; i++)
	{
		f << (atoms[i].get_label().empty() ? constants::atnr2letter(get_atom_charge(i)) : atoms[i].get_label());
		for (int c = 0; c < 3; c++)
			f << setw(14) << setprecision(8) << (isBohr ? constants::bohr2ang(get_atom_coordinate(i, c)) : get_atom_coordinate(i, c));
		f << '\n';
	}
	return f.good();
};

bool WFN::read_fchk(const std::filesystem::path &filename, std::ostream &log, const bool debug)
{
	int r_u_ro_switch = 0;
	std::ifstream fchk(filename, std::ios::in);
	if (!fchk.is_open())
	{
		log << "ERROR while opening .fchk file!" << std::endl;
		return false;
	}
	origin = e_origin::fchk;
	isBohr = true;
	// Every other reader (read_wfn, read_wfx, read_xyz, read_molden, ...) records
	// the source path here; read_fchk did not. The path is what the property code
	// builds cube filenames from, so without it an fchk-driven run wrote
	// "_rho.cube", "_lap.cube", "_fukui_plus.cube" etc. with an empty stem - which
	// silently collide when more than one structure is processed in one directory.
	path = filename;
	std::string line;
	getline_universal(fchk, line);
	std::string title = line;
	getline_universal(fchk, line);
	std::string calculation_level = line;
	if (line[10] == 'R')
	{
		if (line[11] == 'O')
		{
			if (line[12] == '3')
				r_u_ro_switch = 0;
			else
				r_u_ro_switch = 2;
		}
	}
	else if (line[10] == 'U') // Unrestricted
		r_u_ro_switch = 1;
	charge = read_fchk_integer(fchk, "Charge", false);
	multi = read_fchk_integer(fchk, "Multiplicity", false);
	const int el = read_fchk_integer(fchk, "Number of electrons", false);
	getline_universal(fchk, line);
	const int ael = read_fchk_integer(line);
	getline_universal(fchk, line);
	const int bel = read_fchk_integer(line);
	err_checkf(el == ael + bel, "Error in number of electrons!", log);
	if (ael != bel && r_u_ro_switch == 0)
		r_u_ro_switch = 1; // If U was not correctly recognized
	if (calculation_level.find("CASSCF") != std::string::npos && ael != bel)
		r_u_ro_switch = 2; // CASSCF requires open shell treatment
	const int nbas = read_fchk_integer(fchk, "Number of basis functions");
	line = go_get_string(fchk, "Virial Ratio");
	if (line != "")
		virial_ratio = read_fchk_double(line);
	line = go_get_string(fchk, "Total Energy");
	if (line != "")
		total_energy = read_fchk_double(line);
	ivec atnbrs;
	err_checkf(read_fchk_integer_block(fchk, "Atomic numbers", atnbrs), "Error reading atnbrs", log);
	ncen = static_cast<int>(atnbrs.size());
	atoms.resize(ncen);
	for (int i = 0; i < ncen; i++)
		atoms[i].set_label(constants::atnr2letter(atnbrs[i]));
	vec charges;
	err_checkf(read_fchk_double_block(fchk, "Nuclear charges", charges), "Error reading charges", log);
	for (int i = 0; i < charges.size(); i++)
		atoms[i].set_charge(static_cast<int>(charges[i]));
	vec coords;
	err_checkf(read_fchk_double_block(fchk, "Current cartesian coordinates", coords), "Error reading coordinates", log);
	if (coords.size() != ncen * 3)
	{
		log << "Inconsistant number of atoms and coordinates" << std::endl;
		return false;
	}
	for (int i = 0; i < ncen; i++)
	{
		atoms[i].set_coordinate(0, coords[3 * i]);
		atoms[i].set_coordinate(1, coords[3 * i + 1]);
		atoms[i].set_coordinate(2, coords[3 * i + 2]);
	}
	ivec shell_types;
	err_checkf(read_fchk_integer_block(fchk, "Shell types", shell_types, false), "Error reading shell types", log);
	bool is_spherical = false;
	for (int i = 0; i < shell_types.size(); i++)
		if (shell_types[i] < -1)
			is_spherical = true;
	if (debug)
		log << "This fchk contains spherical harmonics, which will be transformed into cartesian functions!" << std::endl
		<< "Loading basis set information..." << std::endl;
	ivec nr_prims_shell;
	err_checkf(read_fchk_integer_block(fchk, "Number of primitives per shell", nr_prims_shell), "Error reading primitives per shell", log);
	ivec shell2atom;
	err_checkf(read_fchk_integer_block(fchk, "Shell to atom map", shell2atom), "Error reading shell2atom", log);
	vec exp;
	err_checkf(read_fchk_double_block(fchk, "Primitive exponents", exp), "Error reading Primitive exponents", log);
	vec con;
	err_checkf(read_fchk_double_block(fchk, "Contraction coefficients", con), "Error reading Contraction coefficients", log);
	vec2 coef(2);
	vec2 MOocc(2), MOene(2);
	if (r_u_ro_switch == 0 || r_u_ro_switch == 2)
	{ // Restricted or Restricted-Open-Shell
		err_checkf(read_fchk_double_block(fchk, "Alpha Orbital Energies", MOene[0]), "Error during reading of Alpha Energies", log);
		err_checkf(read_fchk_double_block(fchk, "MO coefficients", coef[0]), "Error during reading of Alpha MOs", log);
		MOocc[0].resize(MOene[0].size());
		if (r_u_ro_switch == 0)
		{
#pragma omp parallel for
			for (int i = 0; i < MOocc[0].size(); i++)
			{
				if (i < ael)
					MOocc[0][i] = 2.0;
				else
					MOocc[0][i] = 0.0;
			}
		}
		else
		{
#pragma omp parallel for
			for (int i = 0; i < MOocc[0].size(); i++)
			{
				if (i < bel)
					MOocc[0][i] = 2.0;
				else if (i < ael)
					MOocc[0][i] = 1.0;
				else
					MOocc[0][i] = 0.0;
			}
		}
	}
	else
	{ // Unrestricted
		is_unrestricted = true;
		err_checkf(read_fchk_double_block(fchk, "Alpha Orbital Energies", MOene[0]), "Error during reading of Alpha Energies", log);
		err_checkf(read_fchk_double_block(fchk, "Beta Orbital Energies", MOene[1]), "Error during reading of Beta Energies", log);
		err_checkf(read_fchk_double_block(fchk, "Alpha MO coefficients", coef[0]), "Error during reading of Alpha MOs", log);
		err_checkf(read_fchk_double_block(fchk, "Beta MO coefficients", coef[1]), "Error during reading of Beta MOs", log);
		MOocc[0].resize(MOene[0].size());
		MOocc[1].resize(MOene[1].size());
#pragma omp parallel for
		for (int i = 0; i < static_cast<int>(MOene[0].size()); i++)
		{
			if (i < ael)
				MOocc[0][i] = 1.0;
			else
				MOocc[0][i] = 0.0;
		}
#pragma omp parallel for
		for (int i = 0; i < static_cast<int>(MOene[1].size()); i++)
		{
			if (i < bel)
				MOocc[1][i] = 1.0;
			else
				MOocc[1][i] = 0.0;
		}
	}
	if (debug)
		log << "Finished reading the file! Transferring to WFN object!" << std::endl;

	std::vector<primitive> prims;
	ivec shells_of_atom(ncen, 0); //the contracted basis goes onto the atoms too, so free_fchk can write this wavefunction again
	for (int a = 0, e = 0; a < shell_types.size(); a++)
	{
		const int l = abs(shell_types[a]), n = nr_prims_shell[a];
		err_checkf(shell_types[a] != -1 && l <= 10, "SP shells and l > 10 are not supported in fchk", log);
		//the MO coefficients refer to normalised contracted functions; NoSpherA2's own fchk writer leaves the contraction unnormalised
		double norm = 0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				norm += con[e + i] * con[e + j] * pow(2 * sqrt(exp[e + i] * exp[e + j]) / (exp[e + i] + exp[e + j]), l + 1.5);
		//primitive norm of x^l for Cartesian shells, the sph2cart tables expect ORCA's pure scaling
		norm *= shell_types[a] < 0 ? constants::sph2cart_norm2[l] : odd_ft(2 * l - 1);
		for (int i = 0; i < n; i++, e++)
		{
			const double c = con[e] / sqrt(norm) * pow(pow(2, 4 * l + 3) * pow(exp[e], 2 * l + 3) / constants::PI3, 0.25);
			prims.emplace_back(shell2atom[a], constants::first_type[l], exp[e], c);
			err_checkf(atoms[shell2atom[a] - 1].push_back_basis_set(exp[e], c, l + 1, shells_of_atom[shell2atom[a] - 1]), "Error pushing back basis", log);
		}
		shells_of_atom[shell2atom[a] - 1]++;
	}
	if (debug)
		log << "I read the basis of " << ncen << " atoms successfully" << std::endl;
	nex = 0;
	for (int i = 0; i < 2; i++) {
		if (MOocc[i].size() == 0)
			break;
		for (int j = 0; j < nbas; j++) {
			push_back_MO(i * nbas + j + 1, MOocc[i][j], MOene[i][j], 0);
			int coef_run = 0, basis_run = 0;
			for (int p = 0; p < nr_prims_shell.size(); p++)
			{
				const int l = abs(shell_types[p]), size = nr_prims_shell[p], n = shell_types[p] < 0 ? constants::n_spher(l) : constants::n_cart(l);
				vec2 shell(n, vec(size));
				for (int m = 0; m < n; m++)
				{
					//pure fchk functions m = 0, +1, -1, ... carry the Gaussian/libcint phase, the sph2cart tables ORCA's: |m| = 3, 4, 7, 8 change sign
					const double phase = shell_types[p] < 0 && ((m + 1) / 2 % 4 == 3 || (m + 1) / 2 % 4 == 0) && m > 0 ? -1.0 : 1.0;
					for (int s = 0; s < size; s++)
						shell[m][s] = phase * coef[i][j * nbas + coef_run + m] * prims[basis_run + s].get_coef();
				}
				if (shell_types[p] < 0)
					push_back_spherical_shell(i * nbas + j, l, shell, prims, basis_run, size);
				else
					push_back_cartesian_shell(i * nbas + j, l, shell, prims, basis_run, size, l == 3 ? gaussian_f_order : nullptr, cart_norm(l).data());
				coef_run += n;
				basis_run += size;
			}
		}
	}
	set_exp_cutoff();
	return true;
};

bool WFN::read_ptb(const std::filesystem::path &filename, std::ostream &file, const bool debug)
{
	origin = e_origin::ptb;
	isBohr = true;
	path = filename;
	if (debug)
		file << "Reading pTB file: " << filename << std::endl;
	std::ifstream inFile(filename, std::ios::binary | std::ios::in);
	if (!inFile)
	{
		std::cerr << "File could not be opened!\n";
		return false;
	}
	inFile.seekg(0, std::ios::beg);
	int one = 2;
	err_checkf(read_block_from_fortran_binary(inFile, &one, sizeof(one)), "Error reading initial number", std::cout);
	err_checkf(one != 2, "Error reading first number in the xtb file!", std::cout);

	int infos[4] = { 0, 0, 0, 0 }; //ncent nbf nmomax nprims
	err_checkf(read_block_from_fortran_binary(inFile, infos, sizeof(infos)), "Error reading sizes of data", std::cout);
	int ncent = infos[0];
	int nbf = infos[1];
	int nmomax = infos[2];
	int nprims = infos[3];
	err_checkf(ncent > 0 && nbf > 0 && nmomax > 0 && nprims > 0, "xtb file announces " + std::to_string(ncent) + " atoms, " + std::to_string(nbf) + " basis functions, " + std::to_string(nmomax) + " MOs and " + std::to_string(nprims) + " primitives", std::cout);

	svec atyp(ncent);
	char temp[3]{ 0, 0, '\0' };
	for (int i = 0; i < ncent; ++i)
	{
		err_checkf(read_block_from_fortran_binary(inFile, temp, sizeof(temp) - 1), "Error reading atom label " + std::to_string(i), std::cout);
		atyp[i] = temp;
		atyp[i].erase(remove(atyp[i].begin(), atyp[i].end(), ' '), atyp[i].end());
	}

	vec x(ncent), y(ncent), z(ncent);
	ivec _charge(ncent);
	for (int i = 0; i < ncent; ++i)
	{
		err_checkf(read_block_from_fortran_binary(inFile, &x[i], sizeof(double)), "Error reading atom data for atom " + std::to_string(i), std::cout);
		err_checkf(read_block_from_fortran_binary(inFile, &y[i], sizeof(double)), "Error reading atom data for atom " + std::to_string(i), std::cout);
		err_checkf(read_block_from_fortran_binary(inFile, &z[i], sizeof(double)), "Error reading atom data for atom " + std::to_string(i), std::cout);
		err_checkf(read_block_from_fortran_binary(inFile, &_charge[i], sizeof(int)), "Error reading atom data for atom " + std::to_string(i), std::cout);
	}

	// making it into the wavefunction data
	for (int i = 0; i < ncent; i++)
	{
		err_checkf(push_back_atom(atom(atyp[i], {}, i, x[i], y[i], z[i], _charge[i])), "Error adding atom to WFN!", file);
	}
	err_checkf(ncen == ncent, "Error adding atoms to WFN!", file);

	ivec lao(nprims), aoatcart(nprims), ipao(nprims);
	for (int i = 0; i < nprims; ++i) err_checkf(read_block_from_fortran_binary(inFile, &lao[i], sizeof(int)), "Error reading basis set information lao of primitive " + std::to_string(i), std::cout);
	for (int i = 0; i < nprims; ++i) err_checkf(read_block_from_fortran_binary(inFile, &aoatcart[i], sizeof(int)), "Error reading basis set information aotcart of primitive " + std::to_string(i), std::cout);
	for (int i = 0; i < nprims; ++i) err_checkf(read_block_from_fortran_binary(inFile, &ipao[i], sizeof(int)), "Error reading basis set information ipao of primitive " + std::to_string(i), std::cout);

	vec exps(nprims), contr(nprims);
	err_checkf(read_block_from_fortran_binary(inFile, exps.data(), exps.size() * sizeof(double)), "Error reading exponents!", std::cout);
	err_checkf(read_block_from_fortran_binary(inFile, contr.data(), contr.size() * sizeof(double)), "Error reading contraction coefs!", std::cout);
	vec occ(nmomax), eval(nmomax);
	err_checkf(read_block_from_fortran_binary(inFile, occ.data(), occ.size() * sizeof(double)), "Error reading occupancies!", std::cout);
	err_checkf(read_block_from_fortran_binary(inFile, eval.data(), eval.size() * sizeof(double)), "Error reading energies!", std::cout);

	vec tempvec((size_t)nbf * (size_t)nmomax);
	err_checkf(read_block_from_fortran_binary(inFile, tempvec.data(), tempvec.size() * sizeof(double)), "Error reading MO coefficients!", std::cout);
	dMatrix2 momat = reshape<dMatrix2>(tempvec, Shape2D(nmomax, nbf));

	//vec tempvec2((size_t)nmomax * (size_t)nmomax);
	//err_checkf(read_block_from_fortran_binary(inFile, tempvec2.data()), "Error reading spherical MO coefficients!", std::cout);

	//Not every xtb version writes the density matrix record; without it DM stays empty and the
	//density is built from the MOs (the old reader silently accepted an all-zero DM here)
	vec Pmat;
	if (inFile.peek() != std::ifstream::traits_type::eof())
	{
		Pmat.resize((size_t)nmomax * (size_t)(nmomax + 1) / 2);
		err_checkf(read_block_from_fortran_binary(inFile, Pmat.data(), Pmat.size() * sizeof(double)), "Error reading density matrix!", std::cout);
	}

	//  Add Basis set information to atoms
	//  This is a cartesian basis
	//  Exp and Contr are given in therms of the correponding functions, thus for p-type basis functions, we get 3-times the same block
	//  lao tells us what type of function we got (s,p,d...)  (s = 1, p = 2-4, d = 5-10, f = 11-20, g = 21-35)
	//  aoatcart is the corresponding atom which we have to assign everything to
	//  ipao tells us the shells for every basis function
		//int shell = 0;
		//for (int prim = 0; prim < nprims;) {
		//    int function_type = 1;
		//    if (lao[prim] == 1) function_type = 1;
		//    else if (lao[prim] >= 2 && lao[prim] <= 4) function_type = 2;
		//    else if (lao[prim] >= 5 && lao[prim] <= 10) function_type = 3;
		//    else if (lao[prim] >= 11 && lao[prim] <= 20) function_type = 4;
		//    else if (lao[prim] >= 21 && lao[prim] <= 35) function_type = 5;
		//    else err_checkf(true, "Error interpreting basis function type in ptb file!", file);
		//    const int n_prim_type = (function_type * (function_type + 1)) / 2;  //Number of cartesian functions per type
		//    int prims_in_this_shell = 1;
		//    while (ipao[prim] == ipao[prim + 1]) {
		//        atoms[aoatcart[prim] - 1].push_back_basis_set(exps[prim], contr[prim], function_type, shell);
		//        prim++;
		//        prims_in_this_shell++;
		//    }
		//    atoms[aoatcart[prim] - 1].push_back_basis_set(exps[prim], contr[prim], function_type, shell); // One extra time to catch the last one
		//    prim++;
		//    shell++;
		//    if (function_type != 1) { //Skip all the repetition
		//        prim += prims_in_this_shell * (n_prim_type-1);
		//    }
		//    if (aoatcart[prim - 1] != aoatcart[prim]) { // Reste the shellcounter, if at the
		//        shell = 0;
		//    }
		//
		//}

	std::shared_ptr<BasisSet> aux_basis = BasisSetLibrary::get_basis_set("ptb-vdzp");
	set_basis_set_ptr(aux_basis->get_data());
	int nr_coefs = 0;
	for (int i = 0; i < get_ncen(); i++)
	{
		int current_charge = get_atom_charge(i) - 1;
		const std::span<const SimplePrimitive> basis = (*aux_basis)[current_charge];
		int size = (int)basis.size();

		//Different loop to keep the original contraction coefficients
		for (int e = 0; e < size; e++)
		{
			push_back_atom_basis_set(i, basis[e].exp, basis[e].coefficient, basis[e].type + 1, basis[e].shell);
		}
		//for (int e = 0; e < 5; e++)
		//{
		//    push_back_atom_basis_set(i, basis[e].exp, basis[e].coefficient, 2, basis[e].shell);
		//}
	}

	int elcount = -get_charge();
	if (debug)
		file << "elcount: " << elcount << std::endl;
	for (int i = 0; i < ncen; i++)
	{
		elcount += get_atom_charge(i);
		elcount -= constants::ECP_electrons_pTB[get_atom_charge(i)];
		atoms[i].set_ECP_electrons(constants::ECP_electrons_pTB[get_atom_charge(i)]);
	}
	if (debug)
		file << "elcount after: " << elcount << std::endl;
	if (multi == 0)
		multi = elcount % 2 + 1;
	err_checkf((elcount % 2 == 0 && multi % 2 == 1) || elcount % 2 == 1 && multi % 2 == 0, "Impossible combination of number of electrons and multiplicity! " + std::to_string(elcount) + " " + std::to_string(multi), std::cout);

	int alpha_els = 0, beta_els = 0, temp_els = elcount;
	while (temp_els > 1)
	{
		alpha_els++;
		beta_els++;
		temp_els -= 2;
		if (debug)
			file << temp_els << "\n";
		err_checkf(alpha_els >= 0 && beta_els >= 0, "Error setting alpha and beta electrons! a or b are negative!", file);
		err_checkf(alpha_els + beta_els <= elcount, "Error setting alpha and beta electrons! Sum a + b > elcount!", file);
		err_checkf(temp_els > -elcount, "Error setting alpha and beta electrons! Ran below -elcount!", file);
	}
	alpha_els += temp_els;
	if (debug)
		file << "al/be els:" << alpha_els << " " << beta_els << std::endl;
	const int mult = get_multi();
	int diff = 0;
	if (mult != 0)
		diff = get_multi() - 1;
	if (debug)
		file << "diff: " << diff << std::endl;
	while (alpha_els - beta_els != diff)
	{
		alpha_els++;
		beta_els--;
		err_checkf(alpha_els >= 0 && beta_els >= 0, "Error setting alpha and beta electrons: " + std::to_string(alpha_els) + "/" + std::to_string(beta_els), file);
	}

	if (debug)
	{
		file << "al/be els after:" << alpha_els << " " << beta_els << std::endl;
		file << "occs: ";
		for (int i = 0; i < nmomax; i++)
		{
			file << occ[i] << " ";
		}
		file << std::endl;
	}

	for (int i = 0; i < nmomax; i++)
	{
		if (i < alpha_els)
			err_checkf(push_back_MO(MO(i, occ[i], eval[i], 0)), "Error adding MO to WFN!", file);
		else {
			err_checkf(push_back_MO(MO(i, occ[i], eval[i], 1)), "Error adding MO to WFN!", file);
			is_unrestricted = true;
		}
	}
	err_checkf(nmo == nmomax, "Error adding MOs to WFN!", file);

	// we need to generate the primitive coefficients from the contr and exp from the momat, then we cann add them MO-wise
	// pTB orders cartesian f as xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz, so its 16 and 17
	// are xyy and yyz where type_vector has them the other way round
	for (int i = 0; i < nprims; i++)
	{
		vec values;
		for (int j = 0; j < nmomax; j++)
		{
			values.push_back(momat(j, ipao[i] - 1) * contr[i]);
		}
		int type = lao[i];
		if (type == 16) type = 17;
		else if (type == 17) type = 16;
		add_primitive(aoatcart[i], type, exps[i], values.data());
	}

	//Now turn Pmat into a full matrix
	if (!Pmat.empty())
	{
		DM = dMatrix2(nmomax, nmomax);
		double *pmat_ptr = Pmat.data();
		for (int j = 0; j < nmomax; j++) {
			for (int i = 0; i < j; i++) {
				const double v = *pmat_ptr++;
				DM(i, j) = v;
				DM(j, i) = v;
			}
			DM(j, j) = *pmat_ptr++;
		}
	}

	////If i ever need it again, we can reorder the orbitals
	////""" Reorder L=1 components from +1,-1,0 to -1,0,+1 in the overlap matrix"""
	//auto get_new_index = [](const int l, const int m_idx) {
	//    switch (l) {
	//    case 1: { constexpr std::array<int, 3>  map = { 2,0,1 }; return map[m_idx]; }
	//    default: return m_idx;
	//    }
	//    };
	//ivec permutations(nmomax);
	//size_t ao = 0;
	//for (const atom& at : atoms) {
	//    int prim = 0;
	//    for (unsigned int shell = 0; shell < at.get_shellcount_size(); ++shell) {
	//        const int l = at.get_basis_set_type(prim) - 1;
	//        const size_t shell_start = ao;
	//        const int l21 = 2 * l + 1;
	//        for (int m_idx = 0; m_idx < l21; m_idx++) {
	//            const size_t old_idx = shell_start + size_t(m_idx);
	//            const size_t new_idx = shell_start + size_t(get_new_index(l, m_idx));

	//            // perm[old] = new
	//            permutations[old_idx] = int(new_idx);
	//        }

	//        ao += size_t(l21);
	//        prim += at.get_shellcount(shell);
	//    }
	//}
	//for (int j = 0; j < nmomax; j++) {
	//    const int pj = permutations[j];
	//    // Handle diagonal separately (no redundant assignment)
	//    for (int i = 0; i < j; i++) {
	//        const int pi = permutations[i];
	//        const double v = *pmat_ptr++;
	//        DM(pi, pj) = v;
	//        DM(pj, pi) = v;
	//    }
	//    // Diagonal element
	//    DM(pj, pj) = *pmat_ptr++;
	//}



	err_checkf(nprims == nex, "Error adding primitives to WFN!", file);
	inFile.close();
	if (debug)
		this->write_wfn("test_convert_from_xtb.wfn", false, false);
	set_exp_cutoff();
	return true;
}
