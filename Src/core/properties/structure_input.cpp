#include "pch.h"
#include "convenience.h"
#include "cell.h"
#include "wfn_class.h"
#include "atoms.h"

void readxyzMinMax_fromWFN(
	const WFN &wavy,
	properties_options &opts,
	const bool even_steps)
{
	vec2 PosAtoms;
	PosAtoms.resize(3);
	for (int i = 0; i < 3; i++)
		PosAtoms[i].resize(wavy.get_ncen());

	for (int j = 0; j < wavy.get_ncen(); j++)
	{
		PosAtoms[0][j] = wavy.get_atom_coordinate(j, 0);
		PosAtoms[1][j] = wavy.get_atom_coordinate(j, 1);
		PosAtoms[2][j] = wavy.get_atom_coordinate(j, 2);
	}
	opts.MinMax[0] = *std::min_element(PosAtoms[0].begin(), PosAtoms[0].end());
	opts.MinMax[3] = *std::max_element(PosAtoms[0].begin(), PosAtoms[0].end());
	opts.MinMax[1] = *std::min_element(PosAtoms[1].begin(), PosAtoms[1].end());
	opts.MinMax[4] = *std::max_element(PosAtoms[1].begin(), PosAtoms[1].end());
	opts.MinMax[2] = *std::min_element(PosAtoms[2].begin(), PosAtoms[2].end());
	opts.MinMax[5] = *std::max_element(PosAtoms[2].begin(), PosAtoms[2].end());

	const double temp_rad = constants::ang2bohr(opts.radius);
	opts.MinMax[0] -= temp_rad;
	opts.MinMax[3] += temp_rad;
	opts.MinMax[1] -= temp_rad;
	opts.MinMax[4] += temp_rad;
	opts.MinMax[2] -= temp_rad;
	opts.MinMax[5] += temp_rad;

	opts.NbSteps[0] = (int)ceil(constants::bohr2ang(opts.MinMax[3] - opts.MinMax[0]) / opts.resolution);
	opts.NbSteps[1] = (int)ceil(constants::bohr2ang(opts.MinMax[4] - opts.MinMax[1]) / opts.resolution);
	opts.NbSteps[2] = (int)ceil(constants::bohr2ang(opts.MinMax[5] - opts.MinMax[2]) / opts.resolution);

	//Callers that step by span/NbSteps put the box centre, the molecule's symmetry centre, at index NbSteps/2:
	//on a grid plane only for an even count, so only then can the grid carry the point group (an odd count leaves
	//mirror planes and axes between planes). ceil() makes the parity arbitrary (9.0/0.15 gives 61), hence forced here.
	//The property-cube path steps by resolution, where parity does not move the centre; it passes even_steps = false.
	if (even_steps)
		for (int i = 0; i < 3; i++)
			if (opts.NbSteps[i] % 2 != 0)
				opts.NbSteps[i]++;
}

void readxyzMinMax_fromCIF(
	std::filesystem::path cif,
	properties_options &opts,
	vec2 &cm)
{
	using namespace std;
	cell cell(cif);

	cm[0][0] = constants::ang2bohr(cell.get_a());
	cm[0][1] = constants::ang2bohr(cell.get_b() * cell.get_cg());
	cm[0][2] = constants::ang2bohr(cell.get_c() * cell.get_cb());
	cm[1][1] = constants::ang2bohr(cell.get_b() * cell.get_sg());
	cm[1][2] = constants::ang2bohr(cell.get_c() * (cell.get_ca() - cell.get_cb() * cell.get_cg()) / cell.get_sg());
	cm[2][2] = constants::ang2bohr(cell.get_V() / (cell.get_a() * cell.get_b() * cell.get_sg()));

	opts.MinMax[0] = 0.0;
	opts.MinMax[1] = 0.0;
	opts.MinMax[2] = 0.0;

	opts.MinMax[3] = constants::ang2bohr(cell.get_a() + cell.get_b() * cell.get_cg() + cell.get_c() * cell.get_cb());
	opts.MinMax[4] = constants::ang2bohr(cell.get_b() * cell.get_sg() + cell.get_c() * (cell.get_ca() - cell.get_cb() * cell.get_cg()) / cell.get_sg());
	opts.MinMax[5] = cm[2][2];

	opts.NbSteps[0] = (int)ceil(cell.get_a() / opts.resolution);
	opts.NbSteps[1] = (int)ceil(cell.get_b() / opts.resolution);
	opts.NbSteps[2] = (int)ceil(cell.get_c() / opts.resolution);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			cm[i][j] /= opts.NbSteps[j];
}

double read_CIF(std::istream& cif_input, const cell& unit_cell, int& ncen,
	std::vector<asym_atom>& asym_atoms, std::vector<vec2>& ADPs, const bool debug)
{
	if (!cif_input)
		throw std::runtime_error("Could not open CIF file.");
	const double eight_pi2 = 78.95683520871486; // B -> U
	auto is_space = [](char c) { return std::isspace(static_cast<unsigned char>(c)) != 0; };
	// Number without its s.u. "(n)"; '.', '?' and text give false. The classic locale keeps
	// a decimal-comma locale from breaking the parse.
	auto parse_number = [](const std::string& s, double& out) {
		if (s.empty() || s == "." || s == "?")
			return false;
		std::istringstream iss(s.substr(0, s.find('(')));
		iss.imbue(std::locale::classic());
		double v;
		char rest;
		if (!(iss >> v) || (iss >> rest))
			return false;
		out = v;
		return true;
		};

	// ======================= tokenizer =======================
	std::string line;
	size_t line_no = 0, pos = 0;
	bool have_line = false;
	// getline accepting \n, \r\n and bare \r
	auto read_line = [&]() {
		line.clear();
		std::istream::sentry se(cif_input, true);
		if (!se)
			return false;
		std::streambuf* sb = cif_input.rdbuf();
		while (true)
		{
			const int c = sb->sbumpc();
			if (c == '\n')
				break;
			if (c == '\r')
			{
				if (sb->sgetc() == '\n')
					sb->sbumpc();
				break;
			}
			if (c == std::streambuf::traits_type::eof())
			{
				cif_input.setstate(std::ios::eofbit);
				if (line.empty())
					return false;
				break;
			}
			line += static_cast<char>(c);
		}
		++line_no;
		pos = 0;
		have_line = true;
		return true;
		};

	enum TokType { TAG, VALUE, LOOP, DATA, OTHER };
	TokType tok_type = VALUE;
	std::string tok;
	size_t tok_line = 0;
	auto next_token = [&]() {
		while (true)
		{
			if (!have_line)
			{
				if (!read_line())
					return false;
				// ';' in column 1 opens a text field up to the next line starting with ';'.
				// Olex2 embeds the whole FCF and .res this way; it is one value, never parsed.
				if (!line.empty() && line[0] == ';')
				{
					tok_line = line_no;
					tok = line.substr(1);
					bool closed = false;
					while (read_line())
					{
						if (!line.empty() && line[0] == ';')
						{
							closed = true;
							break;
						}
						tok += '\n';
						tok += line;
					}
					if (!closed)
						throw std::runtime_error("CIF: unterminated text field starting at line " + std::to_string(tok_line));
					pos = 1; // anything after the closing ';' is still read
					tok_type = VALUE;
					return true;
				}
			}
			while (pos < line.size() && is_space(line[pos]))
				++pos;
			if (pos >= line.size() || line[pos] == '#')
			{
				have_line = false;
				continue;
			}
			tok_line = line_no;
			const char q = line[pos];
			if (q == '\'' || q == '"')
			{
				// a quote only closes the string if followed by whitespace or EOL ('O'Brien' is one value)
				size_t j = pos + 1;
				while (j < line.size() && !(line[j] == q && (j + 1 == line.size() || is_space(line[j + 1]))))
					++j;
				if (j >= line.size())
					throw std::runtime_error("CIF: unterminated quoted string in line " + std::to_string(line_no));
				tok = line.substr(pos + 1, j - pos - 1);
				pos = j + 1;
				tok_type = VALUE;
				return true;
			}
			size_t j = pos;
			while (j < line.size() && !is_space(line[j]))
				++j;
			tok = line.substr(pos, j - pos);
			pos = j;
			const std::string lw = to_lower(tok);
			if (tok[0] == '_')
			{
				// data names are case-insensitive and "_atom_site.fract_x" (CIF2) == "_atom_site_fract_x"
				tok_type = TAG;
				tok = lw;
				std::replace(tok.begin(), tok.end(), '.', '_');
			}
			else if (lw == "loop_")
				tok_type = LOOP;
			else if (lw.rfind("data_", 0) == 0)
			{
				tok_type = DATA;
				tok = tok.substr(5);
			}
			else if (lw == "global_" || lw == "stop_" || lw.rfind("save_", 0) == 0)
				tok_type = OTHER;
			else
				tok_type = VALUE;
			return true;
		}
		};

	// ======================= parse the first data block =======================
	struct Loop
	{
		std::vector<std::string> tags, values;
		size_t first_line = 0, n_values = 0;
		int col(const std::string& t) const
		{
			const auto it = std::find(tags.begin(), tags.end(), t);
			return it == tags.end() ? -1 : static_cast<int>(it - tags.begin());
		}
		size_t rows() const { return values.size() / tags.size(); }
		const std::string& at(size_t r, int c) const { return values[r * tags.size() + static_cast<size_t>(c)]; }
	};
	std::unordered_map<std::string, std::string> items;
	std::vector<Loop> loops;
	std::string block_name;
	{
		enum State { ITEMS, HEADER, BODY } state = ITEMS;
		Loop cur;
		bool keep = false, in_block = false, done = false;
		std::string pending;
		auto finish_loop = [&]() {
			if (cur.n_values % cur.tags.size() != 0)
				throw std::runtime_error("CIF: loop at line " + std::to_string(cur.first_line) + " (" + cur.tags.front() +
					") has " + std::to_string(cur.n_values) + " values, not a multiple of its " +
					std::to_string(cur.tags.size()) + " columns");
			if (keep)
				loops.push_back(std::move(cur));
			cur = Loop();
			keep = false;
			state = ITEMS;
			};

		bool have_tok = next_token();
		while (have_tok && !done)
		{
			bool consumed = true;
			if (state == HEADER)
			{
				if (tok_type == TAG)
				{
					// Only atom-related loops are stored; reflection lists etc. are just counted.
					keep = keep || tok.rfind("_atom_site_", 0) == 0 || tok.rfind("_atom_type_", 0) == 0
						|| tok.rfind("_diffrn_radiation_", 0) == 0;
					cur.tags.push_back(tok);
				}
				else if (tok_type == VALUE && !cur.tags.empty())
				{
					state = BODY;
					consumed = false;
				}
				else
					throw std::runtime_error("CIF: malformed loop_ at line " + std::to_string(cur.first_line));
			}
			else if (state == BODY)
			{
				// Rows may span several lines (Olex2 wraps anharmonic rows), so values are counted.
				if (tok_type == VALUE)
				{
					++cur.n_values;
					if (keep)
						cur.values.push_back(tok);
				}
				else
				{
					finish_loop();
					consumed = false;
				}
			}
			else if (!pending.empty())
			{
				if (tok_type != VALUE)
					throw std::runtime_error("CIF: data item " + pending + " has no value (line " + std::to_string(tok_line) + ")");
				items[pending] = tok;
				pending.clear();
			}
			else if (tok_type == TAG)
				pending = tok;
			else if (tok_type == LOOP)
			{
				state = HEADER;
				cur.first_line = tok_line;
			}
			else if (tok_type == DATA)
			{
				if (in_block)
					done = true; // only the first data block is read
				else
				{
					in_block = true;
					block_name = tok;
				}
			}
			else if (tok_type == VALUE)
				throw std::runtime_error("CIF: value without a data name in line " + std::to_string(tok_line));
			if (consumed && !done)
				have_tok = next_token();
		}
		if (!pending.empty())
			throw std::runtime_error("CIF: data item " + pending + " has no value at end of file");
		if (state == HEADER)
			throw std::runtime_error("CIF: loop_ without values at line " + std::to_string(cur.first_line));
		if (state == BODY)
			finish_loop();
	}

	// ======================= wavelength =======================
	double wavelength = std::numeric_limits<double>::quiet_NaN();
	for (const char* t : { "_diffrn_radiation_wavelength", "_diffrn_radiation_wavelength_value" })
	{
		const auto it = items.find(t);
		if (it != items.end())
			parse_number(it->second, wavelength);
		else
			for (const Loop& lp : loops)
			{
				const int c = lp.col(t);
				if (c < 0 || lp.rows() == 0)
					continue;
				parse_number(lp.at(0, c), wavelength);
				if (lp.rows() > 1)
					std::cout << "Warning: CIF lists " << lp.rows() << " wavelengths, using the first one (" << lp.at(0, c) << ")\n";
				break;
			}
		if (wavelength == wavelength) // not NaN
			break;
	}
	if (wavelength != wavelength)
		std::cout << "Warning: no usable _diffrn_radiation_wavelength in CIF\n";

	// ======================= per-element anomalous dispersion =======================
	std::unordered_map<std::string, std::complex<double>> type_dispersion; // key: lower-case type symbol
	for (const Loop& lp : loops)
	{
		const int c_sym = lp.col("_atom_type_symbol");
		if (c_sym < 0)
			continue;
		const int c_re = lp.col("_atom_type_scat_dispersion_real");
		const int c_im = lp.col("_atom_type_scat_dispersion_imag");
		for (size_t r = 0; r < lp.rows(); ++r)
		{
			double re = 0.0, im = 0.0;
			if (c_re >= 0) parse_number(lp.at(r, c_re), re);
			if (c_im >= 0) parse_number(lp.at(r, c_im), im);
			type_dispersion[to_lower(lp.at(r, c_sym))] = { re, im };
		}
	}

	// ======================= atom list =======================
	const Loop* site = nullptr;
	for (const Loop& lp : loops)
		if (lp.col("_atom_site_label") >= 0 && lp.col("_atom_site_fract_x") >= 0)
		{
			site = &lp;
			break;
		}
	if (site == nullptr)
		throw std::runtime_error("CIF contains no _atom_site loop with labels and fractional coordinates.");
	const int c_label = site->col("_atom_site_label");
	const int c_type = site->col("_atom_site_type_symbol");
	const int c_x = site->col("_atom_site_fract_x");
	const int c_y = site->col("_atom_site_fract_y");
	const int c_z = site->col("_atom_site_fract_z");
	const int c_uiso = site->col("_atom_site_u_iso_or_equiv");
	const int c_biso = site->col("_atom_site_b_iso_or_equiv");
	if (c_type < 0 || c_y < 0 || c_z < 0)
		throw std::runtime_error("The _atom_site loop needs _atom_site_type_symbol and _atom_site_fract_x/y/z.");

	asym_atoms.clear();
	std::unordered_map<std::string, size_t> index_of;
	for (size_t r = 0; r < site->rows(); ++r)
	{
		asym_atom temp_atom;
		temp_atom.label = site->at(r, c_label);
		const std::string& type_str = site->at(r, c_type);
		temp_atom.type = constants::get_Z_from_label(type_str.c_str()) + 1;

		double fx, fy, fz;
		if (!parse_number(site->at(r, c_x), fx) || !parse_number(site->at(r, c_y), fy) || !parse_number(site->at(r, c_z), fz))
			throw std::runtime_error("Atom " + temp_atom.label + " has no numeric fractional coordinates in the CIF.");
		temp_atom.frac_pos = { fx, fy, fz };
		auto cart = unit_cell.get_coords_cartesian(fx, fy, fz, true);
		temp_atom.pos = { cart[0], cart[1], cart[2] };

		double u = 0.0;
		if (c_uiso >= 0)
			parse_number(site->at(r, c_uiso), u);
		else if (c_biso >= 0 && parse_number(site->at(r, c_biso), u))
			u /= eight_pi2;
		temp_atom.U_iso = u;

		const auto disp = type_dispersion.find(to_lower(type_str));
		temp_atom.anom = disp != type_dispersion.end() ? disp->second : std::complex<double>(0.0, 0.0);
		if (debug && disp == type_dispersion.end())
			std::cout << "No _atom_type dispersion for type " << type_str << " (atom " << temp_atom.label << "), using 0\n";

		if (!index_of.emplace(temp_atom.label, asym_atoms.size()).second)
			throw std::runtime_error("Atom label " + temp_atom.label + " appears twice in the _atom_site loop.");
		asym_atoms.push_back(std::move(temp_atom));
	}
	ADPs.assign(asym_atoms.size(), vec2(3));

	auto atom_index = [&](const std::string& label, const std::string& what) {
		const auto it = index_of.find(label);
		if (it == index_of.end())
			throw std::runtime_error(what + " found for atom '" + label + "' that is not in the _atom_site loop!");
		return it->second;
		};

	// ======================= site-specific anomalous dispersion =======================
	for (const Loop& lp : loops)
	{
		const int c_l = lp.col("_atom_site_dispersion_label");
		if (c_l < 0)
			continue;
		const int c_re = lp.col("_atom_site_dispersion_real");
		const int c_im = lp.col("_atom_site_dispersion_imag");
		for (size_t r = 0; r < lp.rows(); ++r)
		{
			asym_atom& a = asym_atoms[atom_index(lp.at(r, c_l), "Anomalous dispersion")];
			double re = a.anom.real(), im = a.anom.imag();
			if (c_re >= 0) parse_number(lp.at(r, c_re), re);
			if (c_im >= 0) parse_number(lp.at(r, c_im), im);
			a.anom = { re, im };
		}
	}

	// ======================= displacement tensors =======================
	// Fills ADPs[atom][slot] from every loop keyed by label_tag. variants are column prefixes
	// with a scale factor (U_ij as is, B_ij / 8pi^2); the first complete set is used.
	auto read_tensor = [&](const std::string& label_tag, const std::vector<std::pair<std::string, double>>& variants,
		const std::vector<std::string>& suffixes, size_t slot, const std::string& what) {
			for (const Loop& lp : loops)
			{
				const int c_l = lp.col(label_tag);
				if (c_l < 0)
					continue;
				std::vector<int> cols;
				double scale = 1.0;
				for (const auto& v : variants)
				{
					cols.clear();
					for (const auto& s : suffixes)
						cols.push_back(lp.col(v.first + s));
					scale = v.second;
					if (std::find(cols.begin(), cols.end(), -1) == cols.end())
						break;
				}
				if (std::find(cols.begin(), cols.end(), -1) != cols.end())
					throw std::runtime_error("Incomplete " + what + " loop in CIF (line " + std::to_string(lp.first_line) + ")");
				for (size_t r = 0; r < lp.rows(); ++r)
				{
					const std::string& label = lp.at(r, c_l);
					vec& target = ADPs[atom_index(label, what)][slot];
					target.assign(suffixes.size(), 0.0);
					for (size_t k = 0; k < suffixes.size(); ++k)
					{
						if (!parse_number(lp.at(r, cols[k]), target[k]))
							throw std::runtime_error(what + " of atom " + label + " is not a number: " + lp.at(r, cols[k]));
						target[k] *= scale;
					}
				}
			}
		};
	read_tensor("_atom_site_aniso_label",
		{ { "_atom_site_aniso_u_", 1.0 }, { "_atom_site_aniso_b_", 1.0 / eight_pi2 } },
		{ "11", "22", "33", "12", "13", "23" }, 0, "Anisotropic displacement parameters");
	read_tensor("_atom_site_anharm_gc_c_label",
		{ { "_atom_site_anharm_gc_c_", 1.0 } },
		{ "111", "112", "113", "122", "123", "133", "222", "223", "233", "333" }, 1, "Gram-Charlier C coefficients");
	read_tensor("_atom_site_anharm_gc_d_label",
		{ { "_atom_site_anharm_gc_d_", 1.0 } },
		{ "1111", "1112", "1113", "1122", "1123", "1133", "1222", "1223", "1233", "1333",
		  "2222", "2223", "2233", "2333", "3333" }, 2, "Gram-Charlier D coefficients");

	ncen = static_cast<int>(asym_atoms.size());

	if (debug)
	{
		std::cout << "CIF data block: " << block_name << "\nWavelength: " << wavelength << " Angstrom\n";
		for (size_t i = 0; i < asym_atoms.size(); ++i)
		{
			const asym_atom& a = asym_atoms[i];
			std::cout << "  " << a.label << "  type " << a.type
				<< "  frac (" << a.frac_pos[0] << ", " << a.frac_pos[1] << ", " << a.frac_pos[2] << ")"
				<< "  U_iso " << a.U_iso << "  f' " << a.anom.real() << "  f'' " << a.anom.imag()
				<< (ADPs[i][0].empty() ? "" : "  Uij") << (ADPs[i][1].empty() ? "" : "  C")
				<< (ADPs[i][2].empty() ? "" : "  D") << "\n";
		}
		std::cout << "Total atoms parsed: " << ncen << "\n";
	}
	return wavelength;
}
