#include "pch.h"
#include "nbo_run.h"
#include "convenience.h"
#include "wfn_class.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <regex>
#include <sstream>

/*
 * The NBO output is a fixed-width report, so every table here is matched on the shape of its
 * numbers rather than on column positions: an occupancy always carries 5 decimals, an E(2)
 * value 2 and an F(L,NL) 3. That survives the one-space shifts between NBO builds and between
 * restricted and unrestricted output, which column slicing does not.
 */

namespace {
	const std::regex re_npa(R"(^\s*([A-Za-z]{1,2})\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)(?:\s+(-?\d+\.\d+))?\s*$)");
	const std::regex re_nao(R"(^\s*(\d+)\s+([A-Za-z]{1,2})\s+(\d+)\s+(\S+)\s+(\S+)\(\s*([^)]*)\)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$)");
	const std::regex re_nbo_head(R"(^\s*(\d+)\.\s*\(\s*(-?\d+\.\d+)\)\s+(.*\S)\s*$)");
	const std::regex re_nbo_kind(R"(^([A-Z]{2}\*?)\s*\(\s*(\d+)\)\s*(.*)$)");
	const std::regex re_center(R"(([A-Z][a-z]?)\s+(\d+))");
	const std::regex re_comp(R"(([spdfg])\s*[\d.]*\(\s*([\d.]+)%\))");
	const std::regex re_polar(R"(^\s*\(\s*([\d.]+)%\)\s+(-?[\d.]+)\*\s*([A-Z][a-z]?)\s+(\d+)\s+(.*)$)");
	const std::regex re_summary(R"(^\s*(\d+)\.\s+(\S.*?)\s+(-?\d+\.\d{5})\s+(-?\d+\.\d{5}).*$)");
	const std::regex re_e2(R"(^\s*(\d+)\.\s+(\S.*?)\s+(\d+)\.\s+(\S.*?)\s+(-?\d+\.\d{2})\s+(-?\d+\.\d{2})\s+(-?\d+\.\d{3})\s*$)");
	const std::regex re_bo_head(R"(^\s*Atom\s+((?:\d+\s*)+)$)");
	const std::regex re_bo_total(R"(^\s*(\d+)\.\s+([A-Za-z]{1,2})\s+t\s+(.*)$)");
	const std::regex re_bo_cov(R"(^\s*c\s+(.*)$)");
	const std::regex re_bo_ion(R"(^\s*i\s+(.*)$)");
	const std::regex re_weight(R"(^\s*(\d+)\s+(\d+\.\d+)\s*(.*)$)");
	const std::regex re_qpnrt(R"(QPNRT\((\d+)/(\d+)\):\s*D\(0\)=([\d.eE+-]+);\s*D\(w\)=([\d.eE+-]+))");
	const std::regex re_timing(R"(Timing\(sec\):\s*search=([\d.]+);\s*Gram matrix=([\d.]+);\s*minimize=([\d.]+);\s*other=([\d.]+))");
	const std::regex re_version(R"(Cite this program \[(.+?)\])");
	const std::regex re_e2_thresh(R"(^\s*Threshold for printing:\s*([\d.]+) kcal/mol)");
	const std::regex re_e2_inter(R"(^\s*\(Intermolecular threshold:\s*([\d.]+) kcal/mol)");
	const std::regex re_parent_thresh(R"(Parent structure threshold:\s*([\d.]+)% of leading weight)");
	const std::regex re_deloc_thresh(R"(Delocalization list threshold:\s*([\d.]+) kcal/mol)");

	double to_d(const std::string& s) { return s.empty() ? 0.0 : std::stod(s); }

	std::vector<std::string> split_numbers(const std::string& s) {
		std::vector<std::string> out;
		std::istringstream iss(s);
		std::string tok;
		while (iss >> tok) out.push_back(tok);
		return out;
	}

	//Trim and squeeze the inner whitespace of an NBO label so two runs that pad differently
	//still produce the same key.
	std::string normalize(const std::string& s) {
		std::string out;
		bool space = false;
		for (const char c : s) {
			if (std::isspace(static_cast<unsigned char>(c))) { space = !out.empty(); continue; }
			if (space) out.push_back(' ');
			space = false;
			out.push_back(c);
		}
		return out;
	}

	std::vector<int> centers_of(const std::string& text) {
		std::vector<int> out;
		for (std::sregex_iterator it(text.begin(), text.end(), re_center), end; it != end; ++it)
			out.push_back(std::stoi((*it)[2].str()));
		return out;
	}

	void fill_composition(const std::string& text, NboHybrid& h) {
		for (std::sregex_iterator it(text.begin(), text.end(), re_comp), end; it != end; ++it) {
			const double v = to_d((*it)[2].str());
			switch ((*it)[1].str()[0]) {
			case 's': h.s = v; break;
			case 'p': h.p = v; break;
			case 'd': h.d = v; break;
			case 'f': h.f = v; break;
			default: break;
			}
		}
	}
}

static void parse_bond_orders(const std::filesystem::path& nbo_file, NboNrt& nrt);

NboResults parse_nbo_output(const std::filesystem::path& nbo_file) {
	std::ifstream in(nbo_file);
	err_checkf(in.good(), "Could not open NBO output " + nbo_file.string(), std::cout);

	NboResults r;
	r.name = nbo_file.stem().string();

	enum class Section { none, npa, nao, hybrids, summary, e2, weights };
	Section section = Section::none;
	std::string spin;             //"", "alpha", "beta"
	std::map<std::pair<std::string, int>, size_t> orbital_index;  //(spin, NBO number) -> position

	std::string line;
	while (std::getline(in, line)) {
		std::smatch m;
		if (r.version.empty() && std::regex_search(line, m, re_version)) r.version = m[1].str();
		if (line.find("Alpha spin orbitals") != std::string::npos) { spin = "alpha"; r.open_shell = true; section = Section::none; continue; }
		if (line.find("Beta  spin orbitals") != std::string::npos || line.find("Beta spin orbitals") != std::string::npos) { spin = "beta"; r.open_shell = true; section = Section::none; continue; }

		if (line.find("Summary of Natural Population Analysis") != std::string::npos) { section = Section::npa; continue; }
		if (line.find("NATURAL POPULATIONS:") != std::string::npos) { section = Section::nao; continue; }
		if (line.find("Bond orbital / Coefficients / Hybrids") != std::string::npos) { section = Section::hybrids; continue; }
		if (line.find("NATURAL BOND ORBITALS (Summary)") != std::string::npos) { section = Section::summary; continue; }
		if (line.find("SECOND ORDER PERTURBATION THEORY ANALYSIS") != std::string::npos) { section = Section::e2; continue; }
		//The bond-order tables are read by parse_bond_orders, which can hold their three
		//interleaved rows together; here they only end the previous section.
		if (line.find("Natural Bond Order") != std::string::npos) { section = Section::none; r.nrt.present = true; continue; }
		if (line.find("RS   Weight(%)") != std::string::npos) { section = Section::weights; r.nrt.present = true; continue; }
		if (line.find("NATURAL RESONANCE THEORY ANALYSIS") != std::string::npos) { r.nrt.present = true; section = Section::none; continue; }
		if (line.find("TOPO matrix") != std::string::npos || line.find("Natural Atomic Valencies") != std::string::npos) { section = Section::none; continue; }

		if (std::regex_search(line, m, re_qpnrt)) {
			r.nrt.structures_used = std::stoi(m[1].str());
			r.nrt.structures_found = std::stoi(m[2].str());
			r.nrt.d_0 = to_d(m[3].str());
			r.nrt.d_w = to_d(m[4].str());
			continue;
		}
		if (std::regex_search(line, m, re_timing)) {
			r.nrt.search_seconds += to_d(m[1].str());
			r.nrt.gram_seconds += to_d(m[2].str());
			r.nrt.minimize_seconds += to_d(m[3].str());
			r.nrt.other_seconds += to_d(m[4].str());
			continue;
		}
		//An open-shell run prints these once per spin with the same values; keep the first.
		if (r.e2_threshold_kcal == 0.0 && std::regex_search(line, m, re_e2_thresh)) { r.e2_threshold_kcal = to_d(m[1].str()); continue; }
		if (r.e2_intermolecular_threshold_kcal == 0.0 && std::regex_search(line, m, re_e2_inter)) { r.e2_intermolecular_threshold_kcal = to_d(m[1].str()); continue; }
		if (r.nrt.parent_threshold_percent == 0.0 && std::regex_search(line, m, re_parent_thresh)) { r.nrt.parent_threshold_percent = to_d(m[1].str()); continue; }
		if (r.nrt.deloc_threshold_kcal == 0.0 && std::regex_search(line, m, re_deloc_thresh)) { r.nrt.deloc_threshold_kcal = to_d(m[1].str()); continue; }

		switch (section) {
		case Section::npa: {
			if (line.find("* Total *") != std::string::npos) { section = Section::none; break; }
			//Only the first table is kept; an open-shell run repeats it per spin afterwards.
			if (!r.npa.empty() && line.find_first_not_of(" \t") == std::string::npos) break;
			if (!std::regex_match(line, m, re_npa)) break;
			if (!spin.empty()) break;
			NboAtomPopulation a;
			a.element = m[1].str();
			a.index = std::stoi(m[2].str());
			a.charge = to_d(m[3].str());
			a.core = to_d(m[4].str());
			a.valence = to_d(m[5].str());
			a.rydberg = to_d(m[6].str());
			a.total = to_d(m[7].str());
			if (m[8].matched) { a.spin_density = to_d(m[8].str()); a.has_spin_density = true; r.open_shell = true; }
			r.npa.push_back(a);
			break;
		}
		case Section::nao: {
			if (!std::regex_match(line, m, re_nao)) {
				if (!r.nao.empty() && line.find("---") == std::string::npos && line.find_first_not_of(" \t") != std::string::npos) section = Section::none;
				break;
			}
			if (!spin.empty()) break;
			NboNao n;
			n.index = std::stoi(m[1].str());
			n.element = m[2].str();
			n.center = std::stoi(m[3].str());
			n.lang = m[4].str();
			n.type = m[5].str();
			n.shell = normalize(m[6].str());
			n.occupancy = to_d(m[7].str());
			n.energy = to_d(m[8].str());
			r.nao.push_back(n);
			break;
		}
		case Section::hybrids: {
			if (std::regex_match(line, m, re_nbo_head)) {
				std::string rest = m[3].str();
				std::smatch k;
				if (!std::regex_match(rest, k, re_nbo_kind)) break;
				NboOrbital o;
				o.index = std::stoi(m[1].str());
				o.occupancy = to_d(m[2].str());
				o.type = k[1].str();
				o.spin = spin;
				const std::string body = k[3].str();
				//Everything up to the first composition is the atom list; a one-centre NBO
				//carries its composition on this same line.
				const size_t cut = body.find_first_of("(");
				o.description = normalize(k[1].str() + " (" + k[2].str() + ") " + (cut == std::string::npos ? body : body.substr(0, cut)));
				o.centers = centers_of(cut == std::string::npos ? body : body.substr(0, cut));
				if (cut != std::string::npos && !o.centers.empty()) {
					NboHybrid h;
					h.center = o.centers.front();
					h.coefficient = 1.0;
					h.weight_percent = 100.0;
					std::smatch e;
					const std::string head = cut == std::string::npos ? body : body.substr(0, cut);
					if (std::regex_search(head, e, re_center)) h.element = e[1].str();
					fill_composition(body, h);
					o.hybrids.push_back(h);
				}
				orbital_index[{spin, o.index}] = r.orbitals.size();
				r.orbitals.push_back(o);
				break;
			}
			if (std::regex_match(line, m, re_polar) && !r.orbitals.empty()) {
				NboHybrid h;
				h.weight_percent = to_d(m[1].str());
				h.coefficient = to_d(m[2].str());
				h.element = m[3].str();
				h.center = std::stoi(m[4].str());
				fill_composition(m[5].str(), h);
				r.orbitals.back().hybrids.push_back(h);
			}
			break;
		}
		case Section::summary: {
			if (!std::regex_match(line, m, re_summary)) break;
			const auto it = orbital_index.find({ spin, std::stoi(m[1].str()) });
			if (it == orbital_index.end()) break;
			NboOrbital& o = r.orbitals[it->second];
			//The occupancy repeats here; a mismatch means the two tables were not paired up.
			if (std::abs(o.occupancy - to_d(m[3].str())) > 1e-4) break;
			o.energy = to_d(m[4].str());
			break;
		}
		case Section::e2: {
			if (!std::regex_match(line, m, re_e2)) break;
			NboE2Entry e;
			e.donor_index = std::stoi(m[1].str());
			e.donor = normalize(m[2].str());
			e.acceptor_index = std::stoi(m[3].str());
			e.acceptor = normalize(m[4].str());
			e.energy_kcal = to_d(m[5].str());
			e.e_diff = to_d(m[6].str());
			e.fij = to_d(m[7].str());
			e.spin = spin;
			r.e2.push_back(e);
			break;
		}
		case Section::weights: {
			if (line.find("* Total *") != std::string::npos || line.find("---") != std::string::npos) break;
			if (!std::regex_match(line, m, re_weight)) { if (line.find_first_not_of(" \t") != std::string::npos) section = Section::none; break; }
			NboResonanceWeight w;
			w.structure = std::stoi(m[1].str());
			w.weight_percent = to_d(m[2].str());
			w.spin = spin;
			r.nrt.weights.push_back(w);
			break;
		}
		default: break;
		}
	}
	parse_bond_orders(nbo_file, r.nrt);
	return r;
}

/*
 * The bond-order table interleaves three rows per atom (t/c/i) across column blocks, so it is
 * read in a second pass where the three rows can be held together.
 */
static void parse_bond_orders(const std::filesystem::path& nbo_file, NboNrt& nrt) {
	std::ifstream in(nbo_file);
	if (!in.good()) return;
	std::string line, spin;
	bool active = false;
	std::vector<int> columns;
	int row = 0;
	std::vector<std::string> totals, covs;
	std::smatch m;
	auto flush = [&](const std::vector<std::string>& ions) {
		const size_t n = std::min({ columns.size(), totals.size(), covs.size(), ions.size() });
		for (size_t c = 0; c < n; c++) {
			if (columns[c] <= row) continue;         //upper triangle only; "---" marks the diagonal
			if (totals[c] == "---" || covs[c] == "---" || ions[c] == "---") continue;
			NboBondOrder b;
			b.atom1 = row;
			b.atom2 = columns[c];
			b.total = to_d(totals[c]);
			b.covalent = to_d(covs[c]);
			b.ionic = to_d(ions[c]);
			b.spin = spin;
			if (b.total != 0.0 || b.covalent != 0.0 || b.ionic != 0.0) nrt.bond_orders.push_back(b);
		}
		totals.clear(); covs.clear();
	};
	while (std::getline(in, line)) {
		//Closed shell prints "Natural Bond Order:  (total/covalent/ionic)", open shell
		//"Natural Bond Order (alpha spin):  ..." and a composite table, so match the stem only.
		if (line.find("Natural Bond Order") != std::string::npos) {
			active = true;
			columns.clear(); totals.clear(); covs.clear();
			spin = line.find("composite") != std::string::npos ? "composite"
				: line.find("alpha") != std::string::npos ? "alpha"
				: line.find("beta") != std::string::npos ? "beta" : "";
			continue;
		}
		if (!active) continue;
		if (line.find("Natural Atomic Valencies") != std::string::npos || line.find("NATURAL") != std::string::npos) { active = false; continue; }
		if (std::regex_match(line, m, re_bo_head)) {
			columns.clear();
			for (const auto& t : split_numbers(m[1].str())) columns.push_back(std::stoi(t));
			continue;
		}
		if (std::regex_match(line, m, re_bo_total)) { row = std::stoi(m[1].str()); totals = split_numbers(m[3].str()); covs.clear(); continue; }
		if (std::regex_match(line, m, re_bo_cov)) { covs = split_numbers(m[1].str()); continue; }
		if (std::regex_match(line, m, re_bo_ion)) { flush(split_numbers(m[1].str())); continue; }
	}
}

//---------------------------------------------------------------------------------------------
//JSON export; hand-written like BasisSet::write_occ_json, the tree carries no JSON library.
//---------------------------------------------------------------------------------------------
namespace {
	std::string jstr(const std::string& s) {
		std::string out = "\"";
		for (const char c : s) {
			if (c == '"' || c == '\\') { out.push_back('\\'); out.push_back(c); }
			else if (c == '\n') out += "\\n";
			else out.push_back(c);
		}
		return out + "\"";
	}
	std::string jnum(const double v) {
		std::ostringstream o;
		o << std::setprecision(10) << v;
		return o.str();
	}
}

void write_nbo_json(const NboResults& r, const std::filesystem::path& json_file) {
	std::ofstream f(json_file);
	err_checkf(f.good(), "Could not write " + json_file.string(), std::cout);
	f << "{\n";
	f << "  \"name\": " << jstr(r.name) << ",\n";
	f << "  \"source\": " << jstr(r.source) << ",\n";
	f << "  \"nbo_version\": " << jstr(r.version) << ",\n";
	f << "  \"keywords\": " << jstr(r.keywords) << ",\n";
	f << "  \"open_shell\": " << (r.open_shell ? "true" : "false") << ",\n";
	f << "  \"timings\": {\"file47_seconds\": " << jnum(r.file47_seconds) << ", \"nbo_seconds\": " << jnum(r.nbo_seconds) << "},\n";
	f << "  \"thresholds\": {\"e2_kcal\": " << jnum(r.e2_threshold_kcal)
		<< ", \"e2_intermolecular_kcal\": " << jnum(r.e2_intermolecular_threshold_kcal)
		<< ", \"nrt_parent_percent\": " << jnum(r.nrt.parent_threshold_percent)
		<< ", \"nrt_deloc_kcal\": " << jnum(r.nrt.deloc_threshold_kcal) << "},\n";

	f << "  \"npa\": [\n";
	for (size_t i = 0; i < r.npa.size(); i++) {
		const auto& a = r.npa[i];
		f << "    {\"atom\": " << a.index << ", \"element\": " << jstr(a.element)
			<< ", \"charge\": " << jnum(a.charge) << ", \"core\": " << jnum(a.core)
			<< ", \"valence\": " << jnum(a.valence) << ", \"rydberg\": " << jnum(a.rydberg)
			<< ", \"total\": " << jnum(a.total);
		if (a.has_spin_density) f << ", \"spin_density\": " << jnum(a.spin_density);
		f << "}" << (i + 1 < r.npa.size() ? "," : "") << "\n";
	}
	f << "  ],\n";

	f << "  \"nao\": [\n";
	for (size_t i = 0; i < r.nao.size(); i++) {
		const auto& n = r.nao[i];
		f << "    {\"index\": " << n.index << ", \"element\": " << jstr(n.element) << ", \"atom\": " << n.center
			<< ", \"lang\": " << jstr(n.lang) << ", \"type\": " << jstr(n.type) << ", \"shell\": " << jstr(n.shell)
			<< ", \"occupancy\": " << jnum(n.occupancy) << ", \"energy\": " << jnum(n.energy)
			<< "}" << (i + 1 < r.nao.size() ? "," : "") << "\n";
	}
	f << "  ],\n";

	f << "  \"nbos\": [\n";
	for (size_t i = 0; i < r.orbitals.size(); i++) {
		const auto& o = r.orbitals[i];
		f << "    {\"index\": " << o.index << ", \"type\": " << jstr(o.type) << ", \"description\": " << jstr(o.description)
			<< ", \"spin\": " << jstr(o.spin) << ", \"occupancy\": " << jnum(o.occupancy) << ", \"energy\": " << jnum(o.energy)
			<< ", \"centers\": [";
		for (size_t c = 0; c < o.centers.size(); c++) f << (c ? ", " : "") << o.centers[c];
		f << "], \"hybrids\": [";
		for (size_t h = 0; h < o.hybrids.size(); h++) {
			const auto& y = o.hybrids[h];
			f << (h ? ", " : "") << "{\"atom\": " << y.center << ", \"element\": " << jstr(y.element)
				<< ", \"weight_percent\": " << jnum(y.weight_percent) << ", \"coefficient\": " << jnum(y.coefficient)
				<< ", \"s\": " << jnum(y.s) << ", \"p\": " << jnum(y.p) << ", \"d\": " << jnum(y.d) << ", \"f\": " << jnum(y.f) << "}";
		}
		f << "]}" << (i + 1 < r.orbitals.size() ? "," : "") << "\n";
	}
	f << "  ],\n";

	f << "  \"e2\": [\n";
	for (size_t i = 0; i < r.e2.size(); i++) {
		const auto& e = r.e2[i];
		f << "    {\"donor_index\": " << e.donor_index << ", \"donor\": " << jstr(e.donor)
			<< ", \"acceptor_index\": " << e.acceptor_index << ", \"acceptor\": " << jstr(e.acceptor)
			<< ", \"spin\": " << jstr(e.spin) << ", \"energy_kcal\": " << jnum(e.energy_kcal)
			<< ", \"e_diff\": " << jnum(e.e_diff) << ", \"fij\": " << jnum(e.fij)
			<< "}" << (i + 1 < r.e2.size() ? "," : "") << "\n";
	}
	f << "  ],\n";

	f << "  \"nrt\": {\"present\": " << (r.nrt.present ? "true" : "false")
		<< ", \"structures_used\": " << r.nrt.structures_used << ", \"structures_found\": " << r.nrt.structures_found
		<< ", \"d_0\": " << jnum(r.nrt.d_0) << ", \"d_w\": " << jnum(r.nrt.d_w)
		<< ", \"seconds\": {\"search\": " << jnum(r.nrt.search_seconds) << ", \"gram\": " << jnum(r.nrt.gram_seconds)
		<< ", \"minimize\": " << jnum(r.nrt.minimize_seconds) << ", \"other\": " << jnum(r.nrt.other_seconds) << "},\n";
	f << "    \"weights\": [\n";
	for (size_t i = 0; i < r.nrt.weights.size(); i++) {
		const auto& w = r.nrt.weights[i];
		f << "      {\"structure\": " << w.structure << ", \"weight_percent\": " << jnum(w.weight_percent)
			<< ", \"spin\": " << jstr(w.spin) << "}" << (i + 1 < r.nrt.weights.size() ? "," : "") << "\n";
	}
	f << "    ],\n    \"bond_orders\": [\n";
	for (size_t i = 0; i < r.nrt.bond_orders.size(); i++) {
		const auto& b = r.nrt.bond_orders[i];
		f << "      {\"atom1\": " << b.atom1 << ", \"atom2\": " << b.atom2 << ", \"total\": " << jnum(b.total)
			<< ", \"covalent\": " << jnum(b.covalent) << ", \"ionic\": " << jnum(b.ionic)
			<< ", \"spin\": " << jstr(b.spin) << "}" << (i + 1 < r.nrt.bond_orders.size() ? "," : "") << "\n";
	}
	f << "    ]\n  }\n}\n";
}

//---------------------------------------------------------------------------------------------
//Comparison
//---------------------------------------------------------------------------------------------
namespace {
	//Every quantity is compared the same way: build a keyed map of the candidate, walk the
	//reference, record the worst deviation and count what is missing.
	template <typename T, typename KeyFn, typename ValFn>
	NboQuantityAgreement compare_quantity(const std::string& name, const std::vector<T>& ref, const std::vector<T>& cand,
		KeyFn key, ValFn value, const double tol) {
		NboQuantityAgreement q;
		q.quantity = name;
		std::map<std::string, const T*> lookup;
		for (const auto& c : cand) lookup[key(c)] = &c;
		for (const auto& e : ref) {
			const auto it = lookup.find(key(e));
			if (it == lookup.end()) { q.missing++; continue; }
			const double d = std::abs(value(e) - value(*it->second));
			q.compared++;
			if (d > q.max_deviation) { q.max_deviation = d; q.worst = key(e); }
			if (d > tol) q.failed++;
		}
		return q;
	}

	std::string orbital_key(const NboOrbital& o) { return o.spin + "|" + o.description; }
	std::string e2_key(const NboE2Entry& e) { return e.spin + "|" + e.donor + "->" + e.acceptor; }
	std::string bo_key(const NboBondOrder& b) { return b.spin + "|" + std::to_string(b.atom1) + "-" + std::to_string(b.atom2); }
}

std::string NboComparison::report() const {
	std::ostringstream o;
	o << (ok ? "NBO comparison PASSED" : "NBO comparison FAILED") << "\n";
	o << std::left << std::setw(22) << "quantity" << std::right << std::setw(9) << "compared"
		<< std::setw(9) << "missing" << std::setw(8) << "failed" << std::setw(14) << "max dev" << "  worst\n";
	for (const auto& q : quantities) {
		o << std::left << std::setw(22) << q.quantity << std::right << std::setw(9) << q.compared
			<< std::setw(9) << q.missing << std::setw(8) << q.failed
			<< std::setw(14) << std::scientific << std::setprecision(3) << q.max_deviation
			<< "  " << q.worst << "\n";
		o.unsetf(std::ios::scientific);
	}
	for (const auto& msg : messages) o << "  " << msg << "\n";
	return o.str();
}

NboComparison compare_nbo_results(const NboResults& ref, const NboResults& cand, const NboTolerances& tol) {
	NboComparison c;
	if (ref.open_shell != cand.open_shell) c.messages.push_back("open/closed shell disagree");

	c.quantities.push_back(compare_quantity<NboAtomPopulation>("NPA charge", ref.npa, cand.npa,
		[](const NboAtomPopulation& a) { return a.element + std::to_string(a.index); },
		[](const NboAtomPopulation& a) { return a.charge; }, tol.charge));
	if (ref.open_shell)
		c.quantities.push_back(compare_quantity<NboAtomPopulation>("NPA spin density", ref.npa, cand.npa,
			[](const NboAtomPopulation& a) { return a.element + std::to_string(a.index); },
			[](const NboAtomPopulation& a) { return a.spin_density; }, tol.charge));
	c.quantities.push_back(compare_quantity<NboNao>("NAO occupancy", ref.nao, cand.nao,
		[](const NboNao& n) { return std::to_string(n.center) + n.element + "|" + n.lang + "|" + n.type + "|" + n.shell + "|" + std::to_string(n.index); },
		[](const NboNao& n) { return n.occupancy; }, tol.occupancy));
	c.quantities.push_back(compare_quantity<NboOrbital>("NBO occupancy", ref.orbitals, cand.orbitals,
		orbital_key, [](const NboOrbital& o) { return o.occupancy; }, tol.occupancy));
	c.quantities.push_back(compare_quantity<NboOrbital>("NBO energy", ref.orbitals, cand.orbitals,
		orbital_key, [](const NboOrbital& o) { return o.energy; }, tol.energy));
	c.quantities.push_back(compare_quantity<NboOrbital>("NBO %s (hybrid 1)", ref.orbitals, cand.orbitals,
		orbital_key, [](const NboOrbital& o) { return o.hybrids.empty() ? 0.0 : o.hybrids.front().s; }, tol.hybrid_percent));
	c.quantities.push_back(compare_quantity<NboOrbital>("NBO %p (hybrid 1)", ref.orbitals, cand.orbitals,
		orbital_key, [](const NboOrbital& o) { return o.hybrids.empty() ? 0.0 : o.hybrids.front().p; }, tol.hybrid_percent));
	c.quantities.push_back(compare_quantity<NboOrbital>("NBO polarisation", ref.orbitals, cand.orbitals,
		orbital_key, [](const NboOrbital& o) { return o.hybrids.empty() ? 0.0 : o.hybrids.front().weight_percent; }, tol.hybrid_percent));
	c.quantities.push_back(compare_quantity<NboE2Entry>("E2 energy", ref.e2, cand.e2,
		e2_key, [](const NboE2Entry& e) { return e.energy_kcal; }, tol.e2_kcal));
	c.quantities.push_back(compare_quantity<NboBondOrder>("NRT bond order", ref.nrt.bond_orders, cand.nrt.bond_orders,
		bo_key, [](const NboBondOrder& b) { return b.total; }, tol.bond_order));
	c.quantities.push_back(compare_quantity<NboResonanceWeight>("NRT weight", ref.nrt.weights, cand.nrt.weights,
		[](const NboResonanceWeight& w) { return w.spin + "|" + std::to_string(w.structure); },
		[](const NboResonanceWeight& w) { return w.weight_percent; }, tol.weight_percent));

	c.ok = c.messages.empty();
	for (const auto& q : c.quantities) if (q.failed || q.missing) c.ok = false;
	return c;
}

//---------------------------------------------------------------------------------------------
//Driving the external NBO
//---------------------------------------------------------------------------------------------
namespace {
	//gennbo wants a stem in its own working directory. On Windows the licensed binary lives in
	//WSL, so the whole call is handed over with the path translated.
	std::string gennbo_command(const std::filesystem::path& dir, const std::string& stem, const std::string& exe) {
#ifdef _WIN32
		if (exe.empty() || exe.rfind("wsl", 0) == 0) {
			std::string wsl_dir = std::filesystem::absolute(dir).string();
			std::replace(wsl_dir.begin(), wsl_dir.end(), '\\', '/');
			if (wsl_dir.size() > 1 && wsl_dir[1] == ':')
				wsl_dir = "/mnt/" + std::string(1, static_cast<char>(std::tolower(wsl_dir[0]))) + wsl_dir.substr(2);
			const std::string bin = exe.empty() ? "~/nbo7/gennbo" : exe;
			return "wsl bash -lc \"cd '" + wsl_dir + "' && " + bin + " " + stem + "\"";
		}
		return "\"" + exe + "\" " + (dir / stem).string();
#else
		const std::string bin = exe.empty() ? "gennbo" : exe;
		return "cd \"" + std::filesystem::absolute(dir).string() + "\" && " + bin + " " + stem;
#endif
	}
}

int run_nbo(const NboRunOptions& o, std::ostream& log) {
	err_checkf(std::filesystem::exists(o.wavefunction), "No such wavefunction: " + o.wavefunction.string(), log);
	const std::filesystem::path dir = o.work_dir.empty() ? std::filesystem::absolute(o.wavefunction).parent_path() : o.work_dir;
	std::filesystem::create_directories(dir);
	const std::string stem = o.wavefunction.stem().string();
	const std::filesystem::path f47 = dir / (stem + ".47");
	const std::filesystem::path out = dir / (stem + ".nbo");

	WFN wavy(e_origin::NOT_YET_DEFINED);
	wavy.read_known_wavefunction_format(o.wavefunction, log, o.debug);

	auto t0 = std::chrono::steady_clock::now();
	err_checkf(wavy.write_nbo(f47, o.debug, &log, o.keywords), "Could not write " + f47.string(), log);
	const double file47_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

	const std::string command = gennbo_command(dir, stem, o.executable);
	log << "Running: " << command << std::endl;
	t0 = std::chrono::steady_clock::now();
	const int rc = std::system(command.c_str()); /* Flawfinder: ignore - the NBO binary and a path we built */
	const double nbo_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	if (rc != 0) log << "NBO returned " << rc << "; parsing whatever output exists" << std::endl;
	err_checkf(std::filesystem::exists(out), "NBO produced no " + out.string(), log);

	NboResults results = parse_nbo_output(out);
	results.source = o.wavefunction.string();
	results.keywords = o.keywords;
	results.file47_seconds = file47_seconds;
	results.nbo_seconds = nbo_seconds;

	const std::filesystem::path json = o.json_out.empty() ? dir / (stem + ".nbo.json") : o.json_out;
	write_nbo_json(results, json);
	log << "Parsed " << results.npa.size() << " atoms, " << results.nao.size() << " NAOs, "
		<< results.orbitals.size() << " NBOs, " << results.e2.size() << " E2 entries, "
		<< results.nrt.bond_orders.size() << " NRT bond orders -> " << json.string() << std::endl;
	return 0;
}
