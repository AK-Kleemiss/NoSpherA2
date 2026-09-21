#include "pch.h"
//QCT interactive menu (-qct): read, convert, edit and analyse wavefunctions and
//cubes from the terminal. Every analysis routes through the same functions as
//the command-line flags, so the menu carries no numerics of its own.
#include "cube.h"
#include "convenience.h"
#include "fchk.h"
#include "tsc_block.h"
#include "cif.h"
#include "basis_set.h"
#include "wfn_class.h"
#include "b2c.h"
#include "bondwise_analysis.h"
#include "properties.h"

namespace {
constexpr int kWidth = 76;
constexpr int kCol = 38;
// Set once std::cin fails or hits EOF; every loop in the menu ends on it so a
// closed pipe (or a script that ran out of lines) can never spin.
bool input_closed = false;

std::string rule(char c) { return "+" + std::string(kWidth - 2, c) + "+"; }
std::string boxed(std::string text) {
	if ((int)text.size() > kWidth - 4) text.resize(kWidth - 4);
	return "| " + text + std::string(kWidth - 4 - text.size(), ' ') + " |";
}
std::string centred(const std::string& text) {
	const int pad = std::max(0, (kWidth - 4 - (int)text.size()) / 2);
	return boxed(std::string(pad, ' ') + text);
}
std::string two_col(std::string left, const std::string& right) {
	if ((int)left.size() > kCol - 1) left.resize(kCol - 1);
	left.resize(kCol, ' ');
	return boxed(left + right);
}
void banner(const std::string& title) {
	std::cout << rule('=') << '\n' << centred(title) << '\n' << rule('=') << std::endl;
}
void notice(const std::string& text) { std::cout << rule('-') << '\n' << boxed(text) << '\n' << rule('-') << std::endl; }

std::string trim(std::string s) {
	const auto b = s.find_first_not_of(" \t\r\n");
	if (b == std::string::npos) return "";
	return s.substr(b, s.find_last_not_of(" \t\r\n") - b + 1);
}
// One prompt, one line. False means the input is gone, never "empty answer".
bool ask_line(const std::string& prompt, std::string& out) {
	if (input_closed) return false;
	std::cout << "  " << prompt << " > " << std::flush;
	if (!std::getline(std::cin, out)) {
		input_closed = true;
		std::cout << "\n" << std::flush;
		return false;
	}
	out = trim(out);
	return true;
}
template <class T> bool ask(const std::string& prompt, T& value) {
	std::string s;
	if (!ask_line(prompt, s)) return false;
	std::istringstream in(s);
	if (!(in >> value)) {
		std::cout << "  '" << s << "' is not a valid number." << std::endl;
		return false;
	}
	return true;
}
bool confirm(const std::string& prompt) {
	std::string s;
	return ask_line(prompt + " (y/n)", s) && !s.empty() && (s[0] == 'y' || s[0] == 'Y');
}
// Some WFN editors (change_center & co) still read with std::cin >>; drop the tail of their line.
void drop_rest_of_line() {
	if (std::cin.good()) std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	else input_closed = true;
}
char key_of(const std::string& s) { return s.empty() ? '\0' : (char)std::toupper((unsigned char)s[0]); }

std::string origin_name(const e_origin o) {
	switch (o) {
	case e_origin::CRYSTAL: return "crystal";
	case e_origin::wfn: return "wfn";
	case e_origin::cub: return "cube";
	case e_origin::ffn: return "ffn";
	case e_origin::fchk: return "fchk";
	case e_origin::wfx: return "wfx";
	case e_origin::xyz: return "xyz";
	case e_origin::molden: return "molden";
	case e_origin::gbw: return "gbw";
	case e_origin::tonto: return "tonto";
	case e_origin::xtb: return "xtb";
	case e_origin::ptb: return "ptb";
	case e_origin::OCC: return "occ";
	case e_origin::XCW_fit: return "xcw-fit";
	default: return "undefined";
	}
}
const std::vector<std::string> kWfnExt = { ".wfn", ".ffn", ".wfx", ".fch", ".fchk", ".FCh", ".FChK", ".FChk",
	".xyz", ".molden", ".gbw", ".xtb", ".stda", ".orbital_energies,restricted", ".MO_energies,r",
	".molecular_orbitals,restricted", ".MOs,r" };
bool is_cube_ext(const std::filesystem::path& p) { return p.extension() == ".cube" || p.extension() == ".cub"; }
bool is_wfn_ext(const std::filesystem::path& p) {
	return std::find(kWfnExt.begin(), kWfnExt.end(), p.extension().string()) != kWfnExt.end();
}

std::string yes_no(bool b) { return b ? "on" : "off"; }

void print_menu(const std::vector<WFN>& wavy, const int active, const bool expert, const options& opt) {
	banner("NoSpherA2 -- QCT interactive menu");
	if (wavy.empty())
		std::cout << boxed("No wavefunction loaded - press R to read one.") << '\n';
	else {
		const WFN& w = wavy[active];
		std::string status = "Active [" + std::to_string(active) + "] " + w.get_path().filename().string()
			+ "  (" + origin_name(w.get_origin()) + ")";
		if (w.get_modified()) status += "  MODIFIED";
		std::cout << boxed(status) << '\n';
		std::ostringstream o;
		o << "atoms " << w.get_ncen() << "  MOs " << w.get_nmo() << "  primitives " << w.get_nex()
			<< "  cubes " << w.get_cube_count() << "  charge " << w.get_charge() << "  mult " << w.get_multi();
		std::cout << boxed(o.str()) << '\n';
		const std::string basis = trim(w.get_basis_set_name()); //the WFN constructor seeds the name with a single space
		std::cout << boxed("basis: " + (basis.empty() ? std::string("<none>") : basis)
			+ "   loaded: " + std::to_string(wavy.size()) + " wavefunction(s)") << '\n';
	}
	std::ostringstream s;
	s << "threads " << opt.threads << "  grid radius " << opt.properties.radius << " A  resolution "
		<< opt.properties.resolution << " A  expert " << yes_no(expert) << "  debug " << yes_no(opt.debug);
	std::cout << boxed(s.str()) << '\n' << rule('-') << '\n'
		<< two_col("FILES", "EDIT") << '\n'
		<< two_col("  R  read a file (wfn / cube)", "  M  modify atoms and primitives") << '\n'
		<< two_col("  X  close active wavefunction", "  O  sort primitives (.wfn/.ffn)") << '\n'
		<< two_col("  S  save / convert", "  B  load a basis set") << '\n'
		<< two_col("  A  activate another wavefunction", "  C  cube arithmetic and basins") << '\n'
		<< rule('-') << '\n'
		<< two_col("ANALYSIS", "SETTINGS") << '\n'
		<< two_col("  P  property cubes (rho, ESP, ...)", "  L  threads (-1 = all)") << '\n'
		<< two_col("  F  Fukui / reactivity table", "  T  grid radius and resolution") << '\n'
		<< two_col("  I  ELI / QTAIM basins", "  E  expert mode toggle") << '\n'
		<< two_col("  G  Roby-Gould bond indices", "  D  debug output toggle") << '\n'
		<< two_col("  N  bonding analysis (planes, NCI)", "  Q  quit") << '\n'
		<< two_col("  U  units check (bohr / Angstrom)", "") << '\n'
		<< rule('-') << '\n'
		<< boxed("Hirshfeld surfaces, tsc/XCW and SALTED stay on the command line (see -h).") << '\n'
		<< rule('=') << std::endl;
}

bool need_wfn(const std::vector<WFN>& wavy) {
	if (!wavy.empty()) return true;
	notice("First you need to read a wavefunction!");
	return false;
}

// Flat numbered list of every cube of every wavefunction; returns (wfn, cube) or false.
bool pick_cube(std::vector<WFN>& wavy, int& w, int& c, const std::string& what) {
	std::vector<std::pair<int, int>> all;
	for (int i = 0; i < (int)wavy.size(); i++)
		for (int j = 0; j < wavy[i].get_cube_count(); j++) all.emplace_back(i, j);
	if (all.empty()) {
		notice("No cubes loaded - read one with R or compute some with P.");
		return false;
	}
	for (size_t n = 0; n < all.size(); n++)
		std::cout << "  " << std::setw(3) << n + 1 << "  [" << all[n].first << "] " << wavy[all[n].first].get_cube_path(all[n].second).filename().string() << '\n';
	int n = 0;
	if (!ask(what + " (number)", n) || n < 1 || n > (int)all.size()) return false;
	w = all[n - 1].first;
	c = all[n - 1].second;
	return true;
}
const cube* loaded_cube(std::vector<WFN>& wavy, int w, int c, bool expert) {
	if (!wavy[w].get_cube_loaded(c) && !wavy[w].read_cube(c, true, false, expert)) {
		notice("Could not read " + wavy[w].get_cube_path(c).string());
		return nullptr;
	}
	return wavy[w].get_cube_ptr(c);
}
std::filesystem::path ask_output(const std::filesystem::path& def) {
	std::string s;
	if (!ask_line("Output file [" + def.string() + "]", s) || s.empty()) return def;
	return s;
}
bool ok_to_write(const std::filesystem::path& p) {
	return !std::filesystem::exists(p) || confirm(p.string() + " exists, overwrite");
}

void read_file(options& opt, std::vector<WFN>& wavy, int& active, bool expert) {
	std::cout << boxed("wavefunctions: " + std::string(".wfn .ffn .wfx .fchk .molden .gbw .xtb .stda .xyz, Tonto MO files")) << '\n'
		<< boxed("cubes: .cube .cub (attached to the active wavefunction)") << std::endl;
	std::string s;
	if (!ask_line("File", s) || s.empty()) return;
	const std::filesystem::path p(s);
	if (!std::filesystem::exists(p)) {
		notice("No such file: " + p.string());
		return;
	}
	if (is_cube_ext(p)) {
		if (wavy.empty()) {
			wavy.emplace_back(e_origin::cub);
			active = 0;
		}
		if (wavy[active].push_back_cube(p.string(), true, expert))
			notice("Attached " + p.filename().string() + " as cube " + std::to_string(wavy[active].get_cube_count() - 1));
		else
			notice("Could not read " + p.string());
		return;
	}
	if (!is_wfn_ext(p)) {
		notice("Unknown extension '" + p.extension().string() + "'.");
		return;
	}
	// ponytail: a malformed file still leaves through err_checkf, as on the command line.
	wavy.emplace_back(p, opt.debug);
	active = (int)wavy.size() - 1;
	notice("Read " + p.filename().string() + ": " + std::to_string(wavy[active].get_ncen()) + " atoms, "
		+ std::to_string(wavy[active].get_nmo()) + " MOs, " + std::to_string(wavy[active].get_nex()) + " primitives");
}

void save_menu(options& opt, std::vector<WFN>& wavy, int active, bool expert) {
	WFN& w = wavy[active];
	std::cout << rule('-') << '\n' << boxed("SAVE / CONVERT   active: " + w.get_path().filename().string()) << '\n' << rule('-') << '\n'
		<< two_col("  1  .wfn  (all MOs)", "  6  .xyz  (coordinates, Angstrom)") << '\n'
		<< two_col("  2  .wfn  (occupied MOs only)", "  7  .cif  (atoms + wavefunction)") << '\n'
		<< two_col("  3  .wfx", "  8  cube -> .cube") << '\n'
		<< two_col("  4  .fchk (needs a basis set)", "  9  cube -> .dgrid") << '\n'
		<< two_col("  5  .47   (NBO archive)", " 10  cube -> .xdgrid") << '\n'
		<< two_col("  0  back", "") << '\n' << rule('-') << std::endl;
	int sel = 0;
	if (!ask("Choice", sel) || sel == 0) return;
	std::filesystem::path def = w.get_path();
	bool ok = false;
	switch (sel) {
	case 1:
	case 2:
	case 3:
	case 5:
	case 6:
	case 7: {
		static const std::map<int, std::string> ext = { {1, ".wfn"}, {2, ".wfn"}, {3, ".wfx"}, {5, ".47"}, {6, ".xyz"}, {7, ".cif"} };
		const std::filesystem::path out = ask_output(def.replace_extension(ext.at(sel)));
		if (input_closed || !ok_to_write(out)) return;
		if (sel == 1) ok = w.write_wfn(out, opt.debug, false);
		else if (sel == 2) ok = w.write_wfn(out, opt.debug, true);
		else if (sel == 3) ok = w.write_wfx(out, false);
		else if (sel == 5) ok = w.write_nbo(out, opt.debug, &std::cout);
		else if (sel == 6) ok = w.write_xyz(out);
		else {
			write_wfn_CIF(w, out);
			ok = std::filesystem::exists(out);
		}
		notice(ok ? "Written " + out.string() : "Writing " + out.string() + " failed.");
		return;
	}
	case 4: {
		std::string basis;
		if (!ask_line("Basis set name [" + (w.get_basis_set_name().empty() ? "def2-SVP" : w.get_basis_set_name()) + "]", basis)) return;
		if (!basis.empty()) w.set_basis_set_name(basis);
		else if (w.get_basis_set_name().empty()) w.set_basis_set_name("def2-SVP");
		const std::filesystem::path out = ask_output(def.replace_extension(".fchk"));
		if (input_closed || !ok_to_write(out)) return;
		w.assign_charge(w.calculate_charge());
		if (w.get_multi() == 0 && !w.guess_multiplicity(std::cout)) {
			notice("Could not guess the multiplicity.");
			return;
		}
		ok = free_fchk(std::cout, out, opt.basis_set_path, w, opt.debug, true);
		notice(ok ? "Written " + out.string() : "Writing " + out.string() + " failed.");
		return;
	}
	case 8:
	case 9:
	case 10: {
		int cw, cc;
		if (!pick_cube(wavy, cw, cc, "Which cube")) return;
		if (!loaded_cube(wavy, cw, cc, expert)) return;
		std::filesystem::path cdef = wavy[cw].get_cube_path(cc);
		cdef.replace_extension(sel == 8 ? ".cube" : sel == 9 ? ".dgrid" : ".xdgrid");
		const std::filesystem::path out = ask_output(cdef);
		if (input_closed || !ok_to_write(out)) return;
		if (sel == 10) wavy[cw].write_cube_xdgraph(cc, out, opt.debug);
		else wavy[cw].write_cube_file(cc, out, opt.debug);
		ok = std::filesystem::exists(out);
		notice(ok ? "Written " + out.string() : "Writing " + out.string() + " failed.");
		return;
	}
	default:
		notice("Sorry, I did not get that.");
	}
}

void modify_menu(std::vector<WFN>& wavy, int active) {
	WFN& w = wavy[active];
	std::cout << rule('-') << '\n' << boxed("MODIFY   active: " + w.get_path().filename().string()) << '\n' << rule('-') << '\n'
		<< two_col("  1  list centres", "  6  add a primitive") << '\n'
		<< two_col("  2  list primitives", "  7  change centre of a primitive") << '\n'
		<< two_col("  3  delete a centre", "  8  change type of a primitive") << '\n'
		<< two_col("  4  delete a primitive", "  9  change exponent of a primitive") << '\n'
		<< two_col("  5  add an atom", " 10  set one MO coefficient") << '\n'
		<< two_col("  0  back", "") << '\n' << rule('-') << std::endl;
	int sel = 0;
	if (!ask("Choice", sel) || sel == 0) return;
	int nr = 0;
	switch (sel) {
	case 1: w.list_centers(); return;
	case 2: w.list_primitives(); return;
	case 3:
		w.list_centers();
		if (!ask("Centre to delete", nr) || nr < 0 || nr >= w.get_ncen()) return notice("No such centre.");
		if (w.remove_center(nr)) { w.set_modified(); notice("Deleted centre " + std::to_string(nr)); }
		else notice("Could not delete centre " + std::to_string(nr));
		return;
	case 4:
		if (!ask("Primitive to delete (0-" + std::to_string(w.get_nex() - 1) + ")", nr) || nr < 0 || nr >= w.get_nex()) return notice("No such primitive.");
		if (w.remove_primitive(nr)) { w.set_modified(); notice("Deleted primitive " + std::to_string(nr)); }
		else notice("Could not delete primitive " + std::to_string(nr));
		return;
	case 5: {
		std::string label;
		double x, y, z;
		int charge;
		if (!ask_line("Label", label) || !ask("x (bohr)", x) || !ask("y (bohr)", y) || !ask("z (bohr)", z) || !ask("nuclear charge", charge)) return;
		if (charge < 1 || charge > 118) return notice("That element is not discovered yet.");
		if (w.push_back_atom(label, x, y, z, charge)) { w.set_modified(); notice("Added " + label); }
		else notice("Could not add the atom.");
		return;
	}
	case 6: {
		int cen, type;
		double e;
		if (!ask("Centre (1-" + std::to_string(w.get_ncen()) + ")", cen) || !ask("Type (1=s 2-4=p ...)", type) || !ask("Exponent", e)) return;
		std::string line;
		if (!ask_line("MO coefficients (" + std::to_string(w.get_nmo()) + " numbers, empty = zeros)", line)) return;
		vec coef(w.get_nmo(), 0.0);
		std::istringstream in(line);
		for (int i = 0; i < w.get_nmo() && (in >> coef[i]); i++);
		if (w.add_primitive(cen, type, e, coef.data())) { w.set_modified(); notice("Added primitive " + std::to_string(w.get_nex() - 1)); }
		else notice("Could not add the primitive.");
		return;
	}
	case 7:
	case 8:
	case 9:
		if (!ask("Primitive (0-" + std::to_string(w.get_nex() - 1) + ")", nr) || nr < 0 || nr >= w.get_nex()) return notice("No such primitive.");
		w.print_primitive(nr);
		if (sel == 7) w.change_center(nr);
		else if (sel == 8) w.change_type(nr);
		else w.change_exponent(nr);
		drop_rest_of_line();
		w.set_modified();
		return;
	case 10: {
		int mo, prim;
		double val;
		if (!ask("MO (0-" + std::to_string(w.get_nmo() - 1) + ")", mo) || mo < 0 || mo >= w.get_nmo()) return notice("No such MO.");
		if (!ask("Primitive (0-" + std::to_string(w.get_nex() - 1) + ")", prim) || prim < 0 || prim >= w.get_nex()) return notice("No such primitive.");
		if (!ask("New coefficient", val)) return;
		if (w.set_MO_coef(mo, prim, val)) { w.set_modified(); notice("Coefficient set."); }
		else notice("Could not set the coefficient.");
		return;
	}
	default:
		notice("Sorry, I did not get that.");
	}
}

void cube_menu(options& opt, std::vector<WFN>& wavy, bool expert) {
	std::cout << rule('-') << '\n' << boxed("CUBES") << '\n' << rule('-') << '\n'
		<< two_col("  1  integrate", "  8  A / B   -> new cube") << '\n'
		<< two_col("  2  integrate |value|", "  9  real-space R value") << '\n'
		<< two_col("  3  write 2x2x2 super cube", " 10  weighted Jaccard similarity") << '\n'
		<< two_col("  4  threshold in place", " 11  mask A where B != 0") << '\n'
		<< two_col("  5  A + B   -> new cube", " 12  mask A where B == 0") << '\n'
		<< two_col("  6  A - B   -> new cube", " 13  mask A where B < threshold") << '\n'
		<< two_col("  7  A * B   -> new cube", " 14  basins (b2c)     15  BCPs") << '\n'
		<< two_col("  0  back", expert ? " 16-19  A += -= *= /= B in place" : "") << '\n' << rule('-') << std::endl;
	int sel = 0;
	if (!ask("Choice", sel) || sel == 0) return;
	if (sel < 1 || sel > 19 || (sel > 15 && !expert)) return notice("Sorry, I did not get that.");
	int wa, ca;
	if (!pick_cube(wavy, wa, ca, sel >= 5 && sel <= 13 || sel >= 16 ? "Cube A" : "Which cube")) return;
	const cube* a = loaded_cube(wavy, wa, ca, expert);
	if (!a) return;
	std::ostringstream o;
	o << std::scientific << std::setprecision(8);
	if (sel <= 4 || sel >= 14) {
		switch (sel) {
		case 1: o << "Integrated value: " << a->sum(); break;
		case 2: o << "Integrated absolute value: " << a->diff_sum(); break;
		case 3: o << "Super cube written: " << wavy[wa].make_super_cube(ca).string(); break;
		case 4: {
			double t;
			if (!ask("Threshold", t)) return;
			o << (wavy[wa].apply_cube_thresh(ca, t) ? "Threshold applied." : "Threshold failed.");
			break;
		}
		case 14: o << (b2c(a, wavy[wa].get_atoms(), opt.debug, false) ? "Basin analysis done." : "Basin analysis failed."); break;
		case 15: o << (b2c(a, wavy[wa].get_atoms(), opt.debug, true) ? "BCP analysis done." : "BCP analysis failed."); break;
		}
		notice(o.str());
		return;
	}
	int wb, cb;
	if (!pick_cube(wavy, wb, cb, "Cube B")) return;
	const cube* b = loaded_cube(wavy, wb, cb, expert);
	if (!b) return;
	bool ok = true;
	switch (sel) {
	case 5: wavy[wa].push_back_cube(*a + *b); break;
	case 6: wavy[wa].push_back_cube(*a - *b); break;
	case 7: wavy[wa].push_back_cube(*a * *b); break;
	case 8: wavy[wa].push_back_cube(*a / *b); break;
	case 9: o << "Real-space R value: " << a->rrs(*b); break;
	case 10: o << "Weighted Jaccard similarity: " << a->jaccard(*b); break;
	case 11: ok = wavy[wa].apply_cube_mask(ca, *b); break;
	case 12: ok = wavy[wa].apply_cube_negative_mask(ca, *b); break;
	case 13: {
		double t;
		if (!ask("Threshold", t)) return;
		ok = wavy[wa].apply_cube_thresh(ca, *b, t);
		break;
	}
	case 16: ok = wavy[wa].cube_add(ca, *b); break;
	case 17: ok = wavy[wa].cube_subtract(ca, *b); break;
	case 18: ok = wavy[wa].cube_multiply(ca, *b); break;
	case 19: ok = wavy[wa].cube_divide(ca, *b); break;
	}
	if (sel >= 5 && sel <= 8) {
		if (wavy[wa].get_cube_ptr(wavy[wa].get_cube_count() - 1)->get_size(0) == 0) {
			wavy[wa].pop_back_cube();
			ok = false;
		}
		else o << "New cube " << wavy[wa].get_cube_count() - 1 << " attached to wavefunction " << wa;
	}
	if (o.str().empty()) o << (ok ? "Operation successful!" : "Operation failed - do the grids match?");
	notice(o.str());
}

// Property cubes go through properties_calculation, which re-reads the file on disk.
void property_menu(options& opt, std::vector<WFN>& wavy, int active) {
	const WFN& w = wavy[active];
	if (w.get_modified()) notice("The wavefunction is modified in memory - save it first, cubes use the file on disk.");
	std::cout << rule('-') << '\n' << boxed("PROPERTY CUBES   written next to " + w.get_path().filename().string()) << '\n' << rule('-') << '\n'
		<< boxed("Space-separated keywords, any of:") << '\n'
		<< boxed("  rho rdg elf eli lap esp def hdef hirsh srho fukui mos mo=<n>[,<m>] hirsh=<atom>") << '\n'
		<< boxed("Example: rho rdg   (NCI plot)      rho esp   (ESP on the density)") << '\n' << rule('-') << std::endl;
	std::string line;
	if (!ask_line("Properties", line) || line.empty()) return;
	properties_options p;
	p.radius = opt.properties.radius;
	p.resolution = opt.properties.resolution;
	std::istringstream in(line);
	std::string tok;
	bool any = false;
	while (in >> tok) {
		any = true;
		if (tok == "rho") p.rho = true;
		else if (tok == "rdg") p.rdg = true;
		else if (tok == "elf") p.elf = true;
		else if (tok == "eli") p.eli = true;
		else if (tok == "lap") p.lap = true;
		else if (tok == "esp") p.esp = true;
		else if (tok == "def") p.def = true;
		else if (tok == "hdef") p.hdef = true;
		else if (tok == "hirsh") p.hirsh = true;
		else if (tok == "srho") p.s_rho = true;
		else if (tok == "fukui") p.fukui = true;
		else if (tok == "mos") p.all_mos = true;
		else if (tok.rfind("mo=", 0) == 0) {
			std::istringstream nums(tok.substr(3));
			std::string n;
			while (std::getline(nums, n, ',')) p.MO_numbers.push_back(std::stoi(n));
		}
		else if (tok.rfind("hirsh=", 0) == 0) { p.hirsh = true; p.hirsh_number = std::stoi(tok.substr(6)); }
		else return notice("Unknown keyword '" + tok + "'.");
	}
	if (!any) return;
	const properties_options saved = opt.properties;
	const std::filesystem::path saved_wfn = opt.wfn;
	opt.properties = p;
	opt.wfn = w.get_path();
	std::cout << "  Calculating (details in NoSpherA2_cube.log) ..." << std::flush;
	properties_calculation(opt);
	std::cout << " done!" << std::endl;
	opt.properties = saved;
	opt.wfn = saved_wfn;
	notice("Cubes written next to " + w.get_path().string());
}

void bonding_menu(options& opt, std::vector<WFN>& wavy, int active) {
	std::cout << rule('-') << '\n' << boxed("BONDING ANALYSIS") << '\n' << rule('-') << '\n'
		<< boxed("  1  cubes in a bond plane (rho / RDG / ELI-D / Laplacian around 2-3 atoms)") << '\n'
		<< boxed("  2  the same for every bond listed in a definition file (-bonds format)") << '\n'
		<< boxed("  3  Laplacian profile along every bond (bond_<i>_<j>_lap.dat)") << '\n'
		<< boxed("  4  ELI-D basins in QTAIM basins of chosen atoms (masked cube)") << '\n'
		<< boxed("  5  Roby-Gould bond indices (RGBI)      6  ELI-D / QTAIM basin integration") << '\n'
		<< boxed("  7  promolecular NCI between fragments (.xyz files)") << '\n'
		<< boxed("  8  dipole moment of the active wavefunction") << '\n'
		<< boxed("  0  back") << '\n' << rule('-') << std::endl;
	int sel = 0;
	if (!ask("Choice", sel) || sel == 0) return;
	if (sel != 7 && !need_wfn(wavy)) return;
	const std::filesystem::path saved_wfn = opt.wfn;
	switch (sel) {
	case 1: {
		WFN& w = wavy[active];
		w.list_centers();
		int mode, a1, a2, a3;
		std::cout << boxed("orientation: 1 = around atom 1   2 = bond 1-2   3 = plane 1-2-3   4 = ring centroid") << std::endl;
		if (!ask("Orientation (1-4)", mode) || !ask("Atom 1 (1-based)", a1) || !ask("Atom 2", a2) || !ask("Atom 3", a3)) return;
		std::string props;
		if (!ask_line("Properties [rho rdg eli lap]", props)) return;
		if (props.empty()) props = "rho rdg eli lap";
		const bool rho = props.find("rho") != std::string::npos, rdg = props.find("rdg") != std::string::npos,
			eli = props.find("eli") != std::string::npos, lap = props.find("lap") != std::string::npos;
		double res[3] = { opt.properties.resolution, opt.properties.resolution, opt.properties.resolution };
		double box[3] = { 0.0, 0.0, 0.0 }; // 0 = default extent from the bond length
		const bond b = do_bonds(w, mode, true, true, res, true, box, a1, a2, a3, opt.debug, false, 1, rho, rdg, eli, lap);
		if (!b.success) return notice("Bond plane calculation failed - see the messages above.");
		for (const auto& [on, suffix] : { std::pair{rho, "_rho"}, {rdg, "_rdg"}, {eli, "_eli"}, {lap, "_lap"} })
			if (on) w.push_back_cube(b.filename + suffix + ".cube", false, false);
		notice("Cubes " + b.filename + "_*.cube attached to the active wavefunction.");
		return;
	}
	case 2: {
		std::string f;
		if (!ask_line("Bond definition file", f) || f.empty()) return;
		if (!std::filesystem::exists(f)) return notice("No such file: " + f);
		notice(autobonds(opt.debug, wavy[active], f, false) == 1 ? "Bondwise analysis done." : "Bondwise analysis failed.");
		return;
	}
	case 3: {
		std::filesystem::path p = wavy[active].get_path();
		bondwise_laplacian_plots(p);
		notice("Laplacian profiles written next to " + p.string());
		return;
	}
	case 4: {
		std::string line;
		if (!ask_line("Atoms (0-based, comma separated)", line) || line.empty()) return;
		ivec idx;
		std::istringstream in(line);
		for (std::string n; std::getline(in, n, ',');) idx.push_back(std::stoi(trim(n)));
		double bg = 0.0;
		if (!ask("Background value outside the basins", bg)) return;
		run_QTAIM_ELI_mask(wavy[active].get_path(), {}, idx, bg, opt, std::cout);
		return;
	}
	case 5:
		if (wavy[active].get_nmo() == 0) return notice("Roby-Gould bond indices need molecular orbitals.");
		{ Roby_information roby(wavy[active], opt.rgbi_group_sets, !opt.rgbi_no_sym, opt.rgbi_orbital_basis == RGBIOrbitalBasis::ANO, opt.rgbi_EVs, opt.rgbi_theta); }
		return;
	case 6: ELI_analysis(wavy[active], opt); return;
	case 7: {
		std::string line;
		if (!ask_line("Fragment .xyz files (space separated)", line) || line.empty()) return;
		pathvec frags;
		std::istringstream in(line);
		for (std::string f; in >> f;) {
			if (!std::filesystem::exists(f)) return notice("No such file: " + f);
			frags.emplace_back(f);
		}
		if (frags.size() < 2) return notice("Need at least two fragments.");
		promolecular_nci_analysis(frags, opt.properties, std::cout);
		notice("Promolecular NCI done.");
		return;
	}
	case 8:
		opt.wfn = wavy[active].get_path();
		dipole_moments(opt, std::cout);
		opt.wfn = saved_wfn;
		return;
	default:
		notice("Sorry, I did not get that.");
	}
}
} // namespace

int QCT(options& opt, std::vector<WFN>& wavy)
{
	using namespace std;
	input_closed = false;
	bool expert = false;
	int active = wavy.empty() ? 0 : (int)wavy.size() - 1;
	while (!input_closed) {
		print_menu(wavy, active, expert, opt);
		string s;
		if (!ask_line("Choice", s)) break;
		switch (key_of(s)) {
		case 'R': read_file(opt, wavy, active, expert); break;
		case 'X':
			if (!need_wfn(wavy)) break;
			if (wavy[active].get_modified() && !confirm("Unsaved changes - close anyway")) break;
			wavy.erase(wavy.begin() + active);
			active = wavy.empty() ? 0 : min(active, (int)wavy.size() - 1);
			notice("Closed.");
			break;
		case 'S': if (need_wfn(wavy)) save_menu(opt, wavy, active, expert); break;
		case 'A': {
			if (!need_wfn(wavy)) break;
			for (int i = 0; i < (int)wavy.size(); i++)
				cout << "  " << setw(3) << i << "  " << wavy[i].get_path().filename().string() << " (" << origin_name(wavy[i].get_origin()) << ")\n";
			int n;
			if (ask("Activate", n) && n >= 0 && n < (int)wavy.size()) active = n;
			else notice("No such wavefunction.");
			break;
		}
		case 'M': if (need_wfn(wavy)) modify_menu(wavy, active); break;
		case 'O':
			if (!need_wfn(wavy)) break;
			if (wavy[active].get_origin() != e_origin::wfn && wavy[active].get_origin() != e_origin::ffn) { notice("I can only sort .wfn/.ffn files!"); break; }
			cout << "  Sorting wavefunction ..." << flush;
			wavy[active].sort_wfn(wavy[active].check_order(opt.debug), opt.debug);
			cout << " done!" << endl;
			break;
		case 'B': {
			if (!need_wfn(wavy)) break;
			string name;
			if (!ask_line("Basis set name [" + (wavy[active].get_basis_set_name().empty() ? "def2-SVP" : wavy[active].get_basis_set_name()) + "]", name)) break;
			if (!name.empty()) wavy[active].set_basis_set_name(name);
			else if (wavy[active].get_basis_set_name().empty()) wavy[active].set_basis_set_name("def2-SVP");
			notice(BasisSetLibrary::read_basis_set_vanilla(opt.basis_set_path, wavy[active], opt.debug) ? "Basis set loaded." : "Could not load the basis set.");
			break;
		}
		case 'C': cube_menu(opt, wavy, expert); break;
		case 'P': if (need_wfn(wavy)) property_menu(opt, wavy, active); break;
		case 'F': {
			if (!need_wfn(wavy)) break;
			const filesystem::path saved = opt.wfn;
			opt.wfn = wavy[active].get_path();
			fukui_analysis(opt, cout);
			opt.wfn = saved;
			break;
		}
		case 'I': if (need_wfn(wavy)) ELI_analysis(wavy[active], opt); break;
		case 'G':
			if (!need_wfn(wavy)) break;
			if (wavy[active].get_nmo() == 0) { notice("Roby-Gould bond indices need molecular orbitals."); break; }
			{ Roby_information roby(wavy[active], opt.rgbi_group_sets, !opt.rgbi_no_sym, opt.rgbi_orbital_basis == RGBIOrbitalBasis::ANO, opt.rgbi_EVs, opt.rgbi_theta); }
			break;
		case 'N': bonding_menu(opt, wavy, active); break;
		case 'U':
			if (!need_wfn(wavy)) break;
			notice(check_bohr(wavy[active], opt.debug) ? "Appears to be in bohr!" : "Appears to be in Angstrom!");
			break;
		case 'L': {
			int t;
			if (!ask("Number of threads (-1 = all)", t)) break;
			if (t < -1 || t == 0) { notice("Invalid value, keeping " + to_string(opt.threads)); break; }
			opt.threads = t;
			if (t > 0) omp_set_num_threads(t);
			notice("Number of threads set to " + to_string(t));
			break;
		}
		case 'T': {
			double r, res;
			if (!ask("Grid radius around the molecule (A)", r) || !ask("Grid resolution (A)", res)) break;
			if (r <= 0 || res <= 0) { notice("Both must be positive."); break; }
			opt.properties.radius = r;
			opt.properties.resolution = res;
			break;
		}
		case 'E': expert = !expert; notice(expert ? "EXPERT MODE!" : "Expert mode off."); break;
		case 'D': opt.debug = !opt.debug; notice(opt.debug ? "Debug output on." : "Debug output off."); break;
		case 'Q':
			if ((unsaved_files(wavy) || any_of(wavy.begin(), wavy.end(), [](const WFN& w) { return w.get_modified(); }))
				&& !confirm("There are unsaved wavefunctions - quit anyway")) break;
			if (!input_closed) { banner("Bye!"); return 0; }
			break;
		default:
			notice("Sorry, I did not get that, could you try it again?");
		}
	}
	notice("Input closed - leaving the QCT menu.");
	return 0;
}
