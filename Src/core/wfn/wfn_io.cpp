#include "pch.h"
#include "tuning.h"
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
#include "citations.h"
#include <Eigen/Dense>
#include <nlohmann/json.hpp>


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

//base.property.txt names the method and level of the density ORCA put into base.mdci.nat or base.mp2nat
static void report_orca_density(const std::filesystem::path &nat, std::ostream &file)
{
	const std::string name = nat.filename().string(), mdci = ".mdci.nat";
	const bool is_mdci = name.size() > mdci.size() && name.compare(name.size() - mdci.size(), mdci.size(), mdci) == 0;
	std::filesystem::path prop = nat.parent_path() / (is_mdci ? name.substr(0, name.size() - mdci.size()) : nat.stem().string());
	prop += ".property.txt";
	std::ifstream in(prop);
	if (!in.good())
		return;
	const auto quoted = [](const std::string &l) {
		const size_t e = l.find_last_of('"');
		const size_t b = e == std::string::npos || e == 0 ? std::string::npos : l.find_last_of('"', e - 1);
		return b == std::string::npos ? std::string() : l.substr(b + 1, e - b - 1);
	};
	std::string line, method, density_method, level;
	bool triples = false;
	while (std::getline(in, line))
	{
		if (line.find("&Method") != std::string::npos)
			method = quoted(line);
		else if (line.find("&Level") != std::string::npos && !method.empty() && method != "SCF")
		{
			density_method = method;
			level = quoted(line);
		}
		triples = triples || method == "MDCI(SD(T))";
	}
	if (density_method.empty())
		return;
	file << "\nORCA density in " << nat.filename().string() << ": " << density_method << ", " << level << "\n";
	if (triples && density_method != "MDCI(SD(T))")
		file << "WARNING: ORCA's MDCI CCSD(T) stores the CCSD density, without (T); -occ with method = \"ccsd(t)\" gives the CCSD(T) density\n";
	if (density_method.rfind("MDCI", 0) == 0 && level != "Unrelaxed density")
		file << "WARNING: this is not the unrelaxed CC density; %mdci Density Unrelaxed gives it\n";
	if (density_method == "MP2" && level != "Relaxed density")
		file << "WARNING: the unrelaxed MP2 density lacks the orbital response; %mp2 Density Relaxed gives it\n";
}

//Eigenvectors U (U(k, j) at [k * n + j]) and eigenvalues of an MO-basis 1-RDM, largest occupation first
static std::pair<vec, vec> natural_orbitals_of(const Eigen::MatrixXd &G)
{
	const int n = static_cast<int>(G.rows());
	const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G);
	std::pair<vec, vec> natural{vec(size_t(n) * n), vec(n)};
	for (int j = 0; j < n; j++)
	{
		const int s = n - 1 - j;
		natural.second[j] = es.eigenvalues()(s);
		for (int k = 0; k < n; k++) natural.first[size_t(k) * n + j] = es.eigenvectors()(k, s);
	}
	return natural;
}

//An ORCA full CI (every orbital active, "Skipping orbital update") writes no natural orbitals: the gbw keeps the guess
//orbitals with diag(gamma) as occupations, ignoring ActOrbs NatOrbs. The CI 1-RDM is base.scfp in base.densities, in
//the gbw's own AO basis. True when it belongs to these orbitals (diagonal = occupations) but is not diagonal in them.
static bool orca_densities_natural(const std::filesystem::path &gbw, const vec &coefs, const vec &occ, const int n, std::pair<vec, vec> &natural, std::ostream &file)
{
	std::filesystem::path dens = gbw, info = gbw;
	dens.replace_extension(".densities");
	info.replace_extension(".densitiesinfo");
	if (!std::filesystem::exists(dens) || !std::filesystem::exists(info))
	{
		file << "NOTE: fractional occupations in " << gbw.filename().string() << " but no " << dens.filename().string()
			<< " to check them against; an ORCA full CI (every orbital active) leaves its guess orbitals, not natural orbitals, in the gbw\n";
		return false;
	}
	//densitiesinfo names the stored matrices "base.<name>" in file order
	std::ifstream in(info, std::ios::binary);
	const std::string raw((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	const std::string stem = gbw.stem().string() + ".";
	int count = 0, scfp = -1;
	for (size_t at = raw.find(stem); at != std::string::npos; at = raw.find(stem, at + 1), count++)
		if (raw.compare(at + stem.size(), 4, "scfp") == 0 && (at + stem.size() + 4 == raw.size() || !std::isgraph(static_cast<unsigned char>(raw[at + stem.size() + 4]))))
			scfp = count;
	const uintmax_t nn = uintmax_t(n) * n;
	if (scfp < 0 || std::filesystem::file_size(dens) != count * nn * sizeof(double))
		return false;
	Eigen::MatrixXd P(n, n), C(n, n);
	std::ifstream d(dens, std::ios::binary);
	d.seekg(std::streamoff(scfp * nn * sizeof(double)));
	d.read(reinterpret_cast<char *>(P.data()), std::streamsize(nn * sizeof(double)));
	if (!d)
		return false;
	for (int p = 0; p < n; p++)
		for (int k = 0; k < n; k++) C(p, k) = coefs[size_t(p) * n + k];
	const Eigen::MatrixXd X = C.partialPivLu().inverse(), G = X * P * X.transpose();
	double diag = 0, off = 0;
	for (int j = 0; j < n; j++)
		for (int k = 0; k < n; k++)
			(j == k ? diag : off) = std::max(j == k ? diag : off, std::abs(G(j, k) - (j == k ? occ[j] : 0.0)));
	if (diag > 1e-6 || off < 1e-6)
		return false;
	natural = natural_orbitals_of(G);
	file << "WARNING: the orbitals in " << gbw.filename().string() << " do not diagonalise the density " << stem << "scfp in "
		<< dens.filename().string() << " (largest off-diagonal " << std::scientific << std::setprecision(2) << off << std::defaultfloat
		<< "; ORCA writes no natural orbitals for a full CI). Using the natural orbitals of that density.\n";
	return true;
}

bool WFN::read_orca_json_density(const std::filesystem::path &filename, std::ostream &file, const bool debug)
{
	std::ifstream in(filename);
	err_checkf(in.good(), "Could not open " + filename.string(), file);
	const nlohmann::json j = nlohmann::json::parse(in);
	err_checkf(j.contains("ORCA Header") && j.contains("Molecule") && j["Molecule"].contains("Densities") && j["Molecule"].contains("MolecularOrbitals"),
		filename.string() + " is no orca_2json export with densities and MOs (base.json.conf: \"Densities\": [\"all\"], \"MOCoefficients\": true)", file);
	const nlohmann::json &mol = j["Molecule"];
	err_checkf(mol.value("HFTyp", std::string()) == "RHF", "ORCA json densities are read for RHF references only", file);
	const char *names[] = { "autocipur", "autocipre", "mdcip" };
	std::string name;
	for (int k = 0; k < 3 && name.empty(); k++)
		if (mol["Densities"].contains(names[k])) name = names[k];
	err_checkf(!name.empty(), "No correlated density (autocipur, autocipre, mdcip) in " + filename.string(), file);
	const nlohmann::json &dens = mol["Densities"][name], &mos = mol["MolecularOrbitals"]["MOs"];
	const int nbf = static_cast<int>(dens.size()), nmo = static_cast<int>(mos.size());
	Eigen::MatrixXd P(nbf, nbf), C(nbf, nmo);
	for (int p = 0; p < nbf; p++)
		for (int q = 0; q < nbf; q++) P(p, q) = dens[p][q].get<double>();
	for (int k = 0; k < nmo; k++)
		for (int p = 0; p < nbf; p++) C(p, k) = mos[k]["MOCoefficients"][p].get<double>();
	//P in the MO basis, X P X^T with X = (C^T C)^-1 C^T: independent of the export's AO order and phases
	const Eigen::MatrixXd X = (C.transpose() * C).ldlt().solve(C.transpose());
	const std::pair<vec, vec> natural = natural_orbitals_of(X * P * X.transpose());
	std::filesystem::path gbw = filename;
	gbw.replace_extension(".gbw");
	file << "\nORCA density " << name << " from " << filename.filename().string() << " as natural orbitals of " << gbw.filename().string() << "\n";
	return read_gbw(gbw, file, debug, false, &natural);
}

bool WFN::known_wavefunction_format(const std::filesystem::path &fileName)
{
	static const svec ext = {".wfn", ".ffn", ".wfx", ".fch", ".fchk", ".FCh", ".FChK", ".FChk", ".xyz", ".molden", ".gbw", ".nat", ".mp2nat",
		".json", ".xtb", ".stda", ".orbital_energies,restricted", ".MO_energies,r", ".molecular_orbitals,restricted", ".MOs,r"};
	return std::find(ext.begin(), ext.end(), fileName.extension().string()) != ext.end() || fileName.string().find("stdout") != std::string::npos;
}

void WFN::read_known_wavefunction_format(const std::filesystem::path &fileName, std::ostream &file, const bool debug)
{
	err_checkf(known_wavefunction_format(fileName), "Unknown filetype!", file);
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
	else if (fileName.extension() == ".gbw" || fileName.extension() == ".nat" || fileName.extension() == ".mp2nat")  //ORCA natural orbitals (base.mdci.nat, base.mp2nat) are gbw format
	{
		err_checkf(read_gbw(fileName, file, debug), "Problem reading gbw file", file);
		if (fileName.extension() != ".gbw")
			report_orca_density(fileName, file);
	}
	else if (fileName.extension() == ".json")
		err_checkf(read_orca_json_density(fileName, file, debug), "Problem reading ORCA json density", file);
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
	declare_ECPs_if_core_electrons_are_missing(file);
};

//An ECP wavefunction describes fewer electrons than its nuclei carry and no format above says so, so every analysis
//filling orbitals from Z would use electrons the basis lacks. The core is declared only when the shortfall matches the
//def2 core counts exactly (the -ECP table); any other shortfall leaves the wavefunction alone.
void WFN::declare_ECPs_if_core_electrons_are_missing(std::ostream &file)
{
	if (has_ECPs || nmo <= 0 || ncen <= 0)
		return;
	const double occupied = count_nr_electrons();
	if (occupied <= 0.0)
		return; //a geometry without orbitals says nothing about electrons
	int table_core = 0;
	for (int i = 0; i < ncen; i++)
		table_core += constants::ECP_core_electrons(constants::ECP_electrons, get_atom_charge(i));
	if (table_core == 0)
		return;
	const long long missing = static_cast<long long>(get_nr_electrons()) - std::llround(occupied);
	if (missing != table_core)
		return;
	//Newline first: this runs between the caller's "Reading: x" and its " done!"
	file << "\nThe orbitals hold " << std::llround(occupied) << " electrons, " << missing
		<< " fewer than the nuclei carry, and that is exactly the def2 ECP core of these atoms: "
		<< "treating them as ECP atoms, as -ECP would.\n";
	set_has_ECPs(true, true, 1);
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
	//The accuracy the primitive screening holds. ELI-D amplifies any truncation, since g = rho tau - |grad rho|^2/4 is a
	//difference of large terms, so a looser default moves basin populations and invents spurious maxima.
	double acc = constants::density_accuracy;
	if (const char *e = tuning("NOS_DENSITY_ACCURACY")) {
		try {
			const double v = std::stod(e);
			if (std::isfinite(v) && v > 0.0 && v < 1.0) acc = v;
		}
		catch (...) {}
	}
	const double cut0 = std::log(acc / get_maximum_MO_coefficient());
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
	//Queued, not printed: the caller is mid-line in "Reading: <file> ... done!".
	citations::queue(citations::Method::Molden);
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
	//Contraction coefficients multiply bare x^l exp(-a r^2) with the contracted shell normalised, as ORCA and orca_2mkl
	//write; the format's documentation describes the other convention and nothing in a file says which. A norm-based
	//detector does not work, the contracted self-overlap depends on which cartesian component the writer normalised; only
	//MO orthonormality discriminates (MoldenConventionTests). A file in the other convention gives a density off by orders of magnitude.
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
	//Every consumer of DM pairs it with an Int_Params overlap, whose shells are sorted by l with libcint component order;
	//the molden lists shells in file order and a pure shell in ORCA's m order, so the coefficients are permuted as in the gbw reader.
	if (spherical) {
		ivec perm(expected_coefs, -1), sph_row(expected_coefs, -1);
		int file_idx = 0, internal = 0;
		for (int a = 0; a < ncen; a++) {
			ivec shell_l;
			int current_shell = -1;
			for (unsigned int s = 0; s < atoms[a].get_basis_set_size(); s++)
				if ((int)atoms[a].get_basis_set_shell(s) != current_shell) {
					shell_l.push_back(atoms[a].get_basis_set_type(s) - 1);
					current_shell++;
				}
			//the internal index of each of this atom's shells: all s shells, then all p, then all d
			ivec internal_off(shell_l.size(), 0);
			const int max_l = shell_l.empty() ? 0 : *std::max_element(shell_l.begin(), shell_l.end());
			for (int l = 0; l <= max_l; l++)
				for (size_t sh = 0; sh < shell_l.size(); sh++)
					if (shell_l[sh] == l) {
						internal_off[sh] = internal;
						internal += 2 * l + 1;
					}
			for (size_t sh = 0; sh < shell_l.size(); sh++) {
				const int l = shell_l[sh], shell_start = file_idx;
				for (int idx = 0; idx <= 2 * l; idx++) {
					//l = 1 is the one shell the file writes cartesian whatever the [5D] flags
					const std::optional<std::size_t> to = (l == 1) ? std::optional<std::size_t>(std::array<std::size_t, 3>{2, 0, 1}[idx])
						: constants::orca_2_pySCF(l, idx);
					err_checkf(to.has_value(), "No component order known for l = " + to_string(l) + " in " + filename.string(), file);
					sph_row[file_idx] = shell_start + static_cast<int>(to.value());
					perm[file_idx++] = internal_off[sh] + static_cast<int>(to.value());
				}
			}
		}
		err_checkf(file_idx == expected_coefs && internal == expected_coefs,
			"Molden basis walk found " + to_string(file_idx) + "/" + to_string(internal) + " of " + to_string(expected_coefs) + " basis functions", file);
		dMatrix2 reordered(nmo, expected_coefs);
		//MO_sph in read_gbw's layout (file shell order, libcint m order), for -ibo
		MO_sph = dMatrix2(expected_coefs, nmo);
		for (int mo = 0; mo < nmo; mo++)
			for (int i = 0; i < expected_coefs; i++) {
				reordered(mo, perm[i]) = m_coefs(mo, i);
				MO_sph(sph_row[i], mo) = m_coefs(mo, i);
			}
		m_coefs = reordered;
	}
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


bool WFN::read_gbw(const std::filesystem::path &filename, std::ostream &file, const bool debug, const bool _has_ECPs, const std::pair<vec, vec> *natural)
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
		//file atom -> WFN atom, -1 for an embedding atom ("Na> 1.0 x y z NewECP ..." in the input: a capped ECP or
		//point charge, written as Z 119 with its charge in ch and no basis). Not part of the molecule, so dropped here.
		ivec idx(at, -1);
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
			if (geo_ints[0] == 119)
				continue;
			idx[a] = ncen;
			string temp = constants::atnr2letter(static_cast<int>(geo_ints[0]));
			err_checkf(temp != "PROBLEM", "Problem identifying atoms!", std::cout);
			err_checkf(push_back_atom(temp,
				geo_vals[0],
				geo_vals[1],
				geo_vals[2],
				static_cast<int>(geo_ints[0])),
				"Error pushing back atom", file);
		}
		err_checkf(ncen > 0, "gbw file without a single basis-carrying atom", file);
		if (ncen < at)
			file << "\nThe gbw file holds " << at - ncen << " embedding atoms (capped ECPs or point charges) next to the "
				<< ncen << " atoms of the molecule: reading the molecule only.\n";
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
			err_checkf(atom_based >= 0 && atom_based < at, "Basis set for atom " + to_string(atom_based + 1) + " of " + to_string(at), file);
			rd(&nr_shells, constants::soi, "nr_shells");
			if (idx[atom_based] < 0)
			{
				err_checkf(nr_shells == 0, "Embedding atom " + to_string(atom_based + 1) + " carries basis functions", file);
				continue;
			}
			atom_based = idx[atom_based];
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
		std::pair<vec, vec> from_density;
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
			if (!natural && operators == 1 && filename.extension() == ".gbw"
				&& std::any_of(occupations[0].begin(), occupations[0].end(), [](const double o) { return std::min({ std::abs(o), std::abs(o - 1), std::abs(o - 2) }) > 1e-6; })
				&& orca_densities_natural(filename, coefficients[0], occupations[0], dimension, from_density, file))
				natural = &from_density;
			if (natural)
			{
				err_checkf(operators == 1 && natural->second.size() == size_t(dimension) && natural->first.size() == coef_nr,
					"natural orbitals need a restricted gbw with " + to_string(dimension) + " MOs", file);
				//C(p, j) = coefficients[p * dimension + j]; natural orbitals C U, no orbital energies
				vec rotated(coef_nr, 0.0);
				for (int p = 0; p < dimension; p++)
					for (int k = 0; k < dimension; k++)
					{
						const double c = coefficients[i][size_t(p) * dimension + k];
						for (int j = 0; j < dimension; j++) rotated[size_t(p) * dimension + j] += c * natural->first[size_t(k) * dimension + j];
					}
				coefficients[i] = rotated;
				occupations[i] = natural->second;
				natural_rotation = natural->first;
				std::fill(energies[i].begin(), energies[i].end(), 0.0);
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
		//Kept before the type sort below, for -eqc to seed occ (layout at get_MO_sph)
		MO_sph = dMatrix2(operators * dimension, dimension);
		for (int s = 0; s < operators; s++)
		{
			const dMatrix2 &src = s == 0 ? reorderd_coefs_s1 : reorderd_coefs_s2;
			for (int r = 0; r < dimension; r++)
				for (int c = 0; c < dimension; c++)
					MO_sph(s * dimension + r, c) = src(r, c);
		}

		//Map to collect the end index of every type
		index = 0;
		int atom_index = 0;
		for (atom &_atom : atoms) {
			std::map<int, int> type_end;
			std::vector<basis_set_entry> basis = _atom.get_basis_set();
			int temp_bas_idx = 0;
			const int atom_start = index;
			for (unsigned int shell = 0; shell < _atom.get_shellcount_size(); shell++) {
				int type = basis[temp_bas_idx].get_type() - 1;
				temp_bas_idx += _atom.get_shellcount(shell);
				const auto higher = type_end.upper_bound(type);
				if (higher == type_end.end()) {
					type_end[type] = index + 2 * type + 1;
				}
				else {
					if (debug) {
						file << "GBW reader: moving out-of-order angular momentum block l=" << type
							 << " for atom " << atom_index + 1 << ", shell " << shell + 1
							 << " into the internal coefficient order.\n";
					}
					//behind the last block of l or below; every block above l moves up by 2l + 1 (ORCA's def2-*D put s, p, d
					//after g, so a second moved shell lands wrong if only the next l is shifted)
					const int to = higher == type_end.begin() ? atom_start : std::prev(higher)->second;
					move_columns(reorderd_coefs_s1.container(), dimension, index, 2 * type + 1, to);

					if (operators == 2) {
						move_columns(reorderd_coefs_s2.container(), dimension, index, 2 * type + 1, to);
					}

					for (auto it = higher; it != type_end.end(); ++it) it->second += 2 * type + 1;
					type_end[type] = to + 2 * type + 1;
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
				//natural orbitals of a correlated density (.mdci.nat) carry small negative occupations, which belong in the DM
				if (occupations[0][oc] != 0.0) {
					coeff_mo_s1(i, oc) = reorderd_coefs_s1(i, oc) * occupations[0][oc];
					coeff_small_s1(i, oc) = reorderd_coefs_s1(i, oc);
				}
				if (operators == 2 && occupations[1][oc] != 0.0) {
					coeff_mo_s2(i, oc) = reorderd_coefs_s2(i, oc) * occupations[1][oc];
					coeff_small_s2(i, oc) = reorderd_coefs_s2(i, oc);
				}
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
		// The ECP pointer follows the MO pointer (byte 32, or 56 with magic -1) and is 0 in an all-electron file. Read whenever
		// present, not only under -ECP: the electron-count fallback needs an exact def2 match and misses charged molecules.
		rf.seekg(ECP_start_bit, ios::beg);
		int64_t ECP_start = 0;
		rd(&ECP_start, sizeof(ECP_start), "ECP pointer");
		if (_has_ECPs)
			err_checkf(ECP_start != 0, "Could not read ECP information location from GBW file!", file);
		if (ECP_start != 0)
		{
			vector<ECP_primitive> ECP_prims;
			bool molecule_ecp = false; //only embedding atoms carry an ECP: the molecule is all-electron
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
			err_checkf(i1 <= at, "ECP block lists more atoms than the geometry", file);
			if (debug)
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
				const int wa = idx[i]; //an embedding atom's record is read past, it is not the molecule's ECP
				if (wa >= 0)
				{
					atoms[wa].set_ECP_electrons(nr_core);
					molecule_ecp = true;
				}
				rd(&max_contract, soi, "max_contract");
				err_checkf(max_contract > 0, "Error reading max_contract in ECPs", file);
				rd(&max_angular, soi, "max_angular");
				err_checkf(max_angular > 0, "Error reading max_angular in ECPs", file);
				rd(&center, soi, "center");
				err_checkf(center >= 0 && center < at, "Error reading center in ECPs", file); //0-based file atom index
				if (debug)
					file << "I read " << Z << " " << temp_0 << " " << nr_core << " " << max_contract << " " << max_angular << "\n";
				for (int l = 0; l < max_angular; l++)
				{
					rd(&exps, soi, "exps");
					err_checkf(exps > 0, "Error reading exps in ECPs", file);
					rd(&type, soi, "type");
					err_checkf(type >= 0, "Error reading type in ECPs", file);
					err_checkf(type < 200, "This type will give me a headache...", file);
					if (debug)
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
						if (debug)
							file << fun << " " << c << " " << e << " " << n << "\n";
						if (wa >= 0)
							ECP_prims.push_back(ECP_primitive(idx[center], type, e, c, static_cast<int>(n)));
					}
				}
				if (wa >= 0)
					for (int _i = 0; _i < ncen; _i++)
						if (atoms[_i].get_charge() == Z)
							atoms[_i].set_ECP_electrons(nr_core);
			}
			if (debug)
			{
				file << "Ended reading" << endl;
			}
			if (molecule_ecp)
			{
				has_ECPs = true;
				if (ECP_m == 0)
					ECP_m = 1; //the def2 core densities, as the electron-count fallback chose
				//The newline first: this runs between the caller's "Reading: x" and its " done!"
				file << "\nThe gbw file's ECP block removes " << get_nr_ECP_electrons()
					<< " core electrons from the orbitals: treating them as ECP atoms, as -ECP would.\n";
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
	//1-based like push_back_atom gives every other reader; RGBI indexes its per-atom blocks by get_nr() - 1
	for (int i = 0; i < ncen; i++)
	{
		atoms[i].set_label(constants::atnr2letter(atnbrs[i]));
		atoms[i].set_nr(i + 1);
	}
	vec charges, frozen;
	err_checkf(read_fchk_double_block(fchk, "Nuclear charges", charges), "Error reading charges", log);
	//"Nuclear charges" is Z minus an ECP core: the atom keeps Z and the core goes to its ECP count, as in the gbw reader
	if (go_get_string(fchk, "ECP-RNFroz", true) != "")
		err_checkf(read_fchk_double_block(fchk, "ECP-RNFroz", frozen), "Error reading ECP-RNFroz", log);
	int ecp_core = 0;
	for (int i = 0; i < charges.size(); i++)
	{
		const int core = i < static_cast<int>(frozen.size()) ? static_cast<int>(std::lround(frozen[i])) : 0;
		atoms[i].set_charge(static_cast<int>(charges[i]) + core);
		if (core > 0)
			atoms[i].set_ECP_electrons(core);
		ecp_core += core;
	}
	if (ecp_core > 0)
	{
		has_ECPs = true;
		if (ECP_m == 0)
			ECP_m = 1;
		log << "\nThe fchk's ECP block removes " << ecp_core << " core electrons from the orbitals: treating them as ECP atoms, as -ECP would.\n";
	}
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
	//An fchk may declare cartesian d/f/g shells and the primitives follow them, but Int_Params rebuilds a spherical libcint
	//basis from l alone, so NPA/NBO, RGBI and Mulliken would work in an AO space the file never used, which no
	//normalisation can repair. Say so at read time, where both counts are known.
	int spherical_nbf = 0, cartesian_shells = 0;
	for (size_t a = 0; a < shell_types.size(); a++)
	{
		const int l = abs(shell_types[a]);
		spherical_nbf += constants::n_spher(l);
		if (shell_types[a] > 0 && l > 1)
			cartesian_shells++;
	}
	if (spherical_nbf != nbas)
		log << "WARNING: this fchk declares " << nbas << " basis functions and " << cartesian_shells
			<< " of its shells are cartesian, so the basis kept for integrals - " << spherical_nbf
			<< " spherical functions - is not the space the file's MO coefficients are written in."
			<< " Cubes, the density fit and .tsc read the primitives and are unaffected; NPA/NBO, RGBI"
			<< " and Mulliken would not be, so use the .gbw or .molden of the same calculation for those."
			<< std::endl;
	nex = 0;
	for (int i = 0; i < 2; i++) {
		if (MOocc[i].size() == 0)
			break;
		for (int j = 0; j < nbas; j++) {
			push_back_MO(i * nbas + j + 1, MOocc[i][j], MOene[i][j], i);
			int coef_run = 0, basis_run = 0;
			for (int p = 0; p < nr_prims_shell.size(); p++)
			{
				const int l = abs(shell_types[p]), size = nr_prims_shell[p], n = shell_types[p] < 0 ? constants::n_spher(l) : constants::n_cart(l);
				vec2 shell(n, vec(size));
				for (int m = 0; m < n; m++)
				{
					//pure fchk functions m = 0, +1, -1, ... carry the Gaussian/libcint phase, the sph2cart tables ORCA's
					const double phase = shell_types[p] < 0 && orca_pure_sign_flips((m + 1) / 2) ? -1.0 : 1.0;
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
	//DM over Int_Params' basis (per atom all s, then p, ...; libcint m order), which RGBI, NPA/NBO and Mulliken pair with
	//their own overlap. Pure fchk components already carry the libcint phase, so only the order changes, as in the molden
	//reader. A cartesian d/f/g shell has no place in that basis, so those files keep DM empty (warned above).
	if (spherical_nbf == nbas)
	{
		ivec file_off(shell_types.size());
		std::vector<ivec> atom_shells(ncen);
		for (int s = 0, o = 0; s < static_cast<int>(shell_types.size()); s++)
		{
			file_off[s] = o;
			o += constants::n_spher(abs(shell_types[s]));
			atom_shells[shell2atom[s] - 1].push_back(s);
		}
		ivec perm(nbas, -1);
		for (int a = 0, internal = 0; a < ncen; a++)
			for (int l = 0; l <= 10; l++)
				for (const int s : atom_shells[a])
				{
					if (abs(shell_types[s]) != l)
						continue;
					//p shells are x, y, z in an fchk, libcint (PYPZPX) wants y, z, x
					for (int idx = 0; idx <= 2 * l; idx++)
						perm[file_off[s] + idx] = internal + static_cast<int>(l == 1 ? std::array<int, 3>{2, 0, 1}[idx] : constants::orca_2_pySCF(l, idx).value());
					internal += 2 * l + 1;
				}
		for (int i = 0; i < 2 && MOocc[i].size() > 0; i++)
		{
			const int nmo = static_cast<int>(MOocc[i].size());
			dMatrix2 C(nmo, nbas);
			for (int j = 0; j < nmo; j++)
				for (int k = 0; k < nbas; k++)
					C(j, perm[k]) = coef[i][j * nbas + k];
			dMatrix2 D = dot(diag_dot(C, MOocc[i], true), C);
			if (i == 0)
				DM = D;
			else
			{
				DM_beta = D;
				for (int r = 0; r < nbas; r++)
					for (int c = 0; c < nbas; c++)
						DM(r, c) += D(r, c);
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
	citations::queue(citations::Method::PTB); //queued: the caller's line is still open
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
		const int core = constants::ECP_core_electrons(constants::ECP_electrons_pTB, get_atom_charge(i));
		elcount += get_atom_charge(i);
		elcount -= core;
		atoms[i].set_ECP_electrons(core);
	}
	if (get_nr_ECP_electrons() > 0)
		set_has_ECPs(true, false, 3);
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
