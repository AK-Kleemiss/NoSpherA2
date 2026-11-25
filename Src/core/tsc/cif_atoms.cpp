#include "pch.h"
#include "scattering_factors.h"
#include "cell.h"

/**
 * Reads atoms from a CIF file and performs necessary operations.
 *
 * @param cif_input The input stream for the CIF file.
 * @param input_groups The vector of input groups.
 * @param unit_cell The cell object representing the unit cell.
 * @param wave The WFN object.
 * @param known_atoms The vector of known atoms.
 * @param atom_type_list The vector of atom type list.
 * @param asym_atom_to_type_list The vector of asymmetric atom to type list.
 * @param asym_atom_list The vector of asymmetric atom list.
 * @param needs_grid The vector indicating whether each atom needs a grid.
 * @param file The output stream for the file.
 * @param debug A boolean indicating whether to enable debug mode.
 */
// Reads exactly n values of one CIF loop row: a row may wrap over several lines and a quoted value ('x, y, z')
// may contain spaces. False if the file ends mid-row; line keeps the last line consumed, as the callers expect.
static bool read_cif_loop_row(std::istream &input, std::string &line, int n, svec &fields)
{
	fields.assign(n, "");
	int got = 0;
	std::string rest = line;
	while (true)
	{
		size_t i = 0;
		while (got < n && i < rest.size())
		{
			while (i < rest.size() && std::isspace(static_cast<unsigned char>(rest[i]))) i++;
			if (i >= rest.size()) break;
			std::string value;
			const char c = rest[i];
			if (c == 0x27 || c == '"')
			{
				const char quote = c;
				i++;
				while (i < rest.size() && rest[i] != quote) value.push_back(rest[i++]);
				if (i < rest.size()) i++;              // closing quote
			}
			else
			{
				while (i < rest.size() && !std::isspace(static_cast<unsigned char>(rest[i])))
					value.push_back(rest[i++]);
			}
			fields[got++] = value;
		}
		if (got >= n) return true;
		if (!getline_universal(input, rest)) return false;
		line = rest;
	}
}


std::vector<std::pair<std::string, atomID>> cif_scatterer_keys(const pathvec& cifs)
{
	using namespace std;
	vector<pair<string, atomID>> out;
	unordered_set<string> seen;
	for (const filesystem::path& cif : cifs)
	{
		ifstream in(cif);
		err_checkf(in.good(), "Could not open " + cif.string(), std::cout);
		string line;
		bool done = false;
		while (!done && getline_universal(in, line))
		{
			if (line.find("loop_") == string::npos)
				continue;
			int n = 0, label = -1, type = -1, group = -1, pos[3] = { -1, -1, -1 };
			getline_universal(in, line);
			while (trim(line).find("_") == 0)
			{
				if (line.find("_atom_site_label") != string::npos) label = n;
				else if (line.find("type_symbol") != string::npos) type = n;
				else if (line.find("disorder_group") != string::npos) group = n;
				else if (line.find("fract_x") != string::npos) pos[0] = n;
				else if (line.find("fract_y") != string::npos) pos[1] = n;
				else if (line.find("fract_z") != string::npos) pos[2] = n;
				getline_universal(in, line);
				n++;
			}
			if (label < 0 || type < 0 || pos[0] < 0 || pos[1] < 0 || pos[2] < 0)
				continue;
			//the rows as read_atoms_from_CIF reads them, so the atomID is the one the table and the .tshb carry
			while (trim(line).find("_") > 0 && trim(line).find("loop_") != 0 && line.length() > 3)
			{
				svec f;
				err_checkf(read_cif_loop_row(in, line, n, f), "CIF ends inside an atom loop row, last line: '" + line + "'", std::cout);
				for (string* s : { &f[label], &f[type] })
					s->erase(remove_if(s->begin(), s->end(), ::isspace), s->end());
				int part = 0;
				if (group >= 0 && f[group] != "." && f[group] != "?" && !f[group].empty())
					try { part = stoi(f[group]); }
					catch (const exception&) { part = 0; }
				const int Z = constants::get_Z_from_label(f[type].c_str()) + 1;
				if (Z > 0)
				{
					const atomID id(stod(f[pos[0]]), stod(f[pos[1]]), stod(f[pos[2]]), part, Z);
					if (seen.insert(id.to_hex_string()).second)
						out.emplace_back(f[label], id);
				}
				getline_universal(in, line);
			}
			done = true;
		}
	}
	return out;
}

svec read_atoms_from_CIF(std::ifstream& cif_input,
	const ivec& input_groups,
	const cell& unit_cell,
	WFN& wave,
	const svec& known_atoms,
	ivec& atom_type_list,
	ivec& asym_atom_to_type_list,
	ivec& asym_atom_list,
	bvec& needs_grid,
	std::ostream& file,
	const bool debug,
	const bool allow_empty)
{
	using namespace std;
	if (debug)
		file << "start working on cif" << endl;
	bool atoms_read = false;
	int count_fields = 0;
	int group_field = -1;
	int type_field = -1;
	int position_field[3] = { -1, -1, -1 };
	int label_field = 1000;
	string line;
	cif_input.clear();
	cif_input.seekg(0, cif_input.beg);
	svec labels(wave.get_ncen(), "");
	if (debug && input_groups.size() > 0)
		file << "Group size: " << input_groups.size() << endl;
	else if (debug)
		file << "Starting search loop" << endl;
	while (!cif_input.eof() && !atoms_read)
	{
		getline_universal(cif_input, line);
		// if (debug)
		//     file << "line: " << line << endl;
		if (line.find("loop_") != string::npos)
		{
			count_fields = 0;
			getline_universal(cif_input, line);
			if (debug)
				file << "line in loop field definition: " << trim(line) << "\n";
			while (trim(line).find("_") == 0)
			{
				if (line.find("label") != string::npos)
					label_field = count_fields;
				else if (line.find("type_symbol") != string::npos)
					type_field = count_fields;
				else if (line.find("disorder_group") != string::npos)
					group_field = count_fields;
				else if (line.find("fract_x") != string::npos)
					position_field[0] = count_fields;
				else if (line.find("fract_y") != string::npos)
					position_field[1] = count_fields;
				else if (line.find("fract_z") != string::npos)
					position_field[2] = count_fields;
				else if (label_field == 1000)
				{
					if (debug)
						file << "I don't think this is the atom block.. moving on!\n";
					break;
				}
				getline_universal(cif_input, line);
				count_fields++;
			}
			if (label_field != 1000) {
				err_checkf(position_field[0] != -1, "No x position found, impossible to continue!", std::cout);
				err_checkf(position_field[1] != -1, "No y position found, impossible to continue!", std::cout);
				err_checkf(position_field[2] != -1, "No z position found, impossible to continue!", std::cout);
				err_checkf(type_field != -1, "No type found, impossible to continue!", std::cout);
			}
			// the next loop_ or data item ends the atom rows
			while (trim(line).find("_") > 0 && trim(line).find("loop_") != 0 && line.length() > 3)
			{
				atoms_read = true;
				svec fields;
				int nr = -1;
				err_checkf(read_cif_loop_row(cif_input, line, count_fields, fields), "CIF ends inside an atom loop row, last line: '" + line + "'", file);
				fields[label_field].erase(remove_if(fields[label_field].begin(), fields[label_field].end(), ::isspace), fields[label_field].end());
				fields[type_field].erase(remove_if(fields[type_field].begin(), fields[type_field].end(), ::isspace), fields[type_field].end());
				if (debug)
					file << "label: " << setw(8) << fields[label_field] << " type: " << fields[type_field] << " frac. pos: "
					<< setw(6) << fixed << setprecision(3) << stod(fields[position_field[0]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[0]]) << " "
					<< setw(6) << fixed << setprecision(3) << stod(fields[position_field[1]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[1]]) << " "
					<< setw(6) << fixed << setprecision(3) << stod(fields[position_field[2]]) << "+/-" << get_decimal_precision_from_CIF_number(fields[position_field[2]]) << " ";
				vec position = unit_cell.get_coords_cartesian(
					stod(fields[position_field[0]]),
					stod(fields[position_field[1]]),
					stod(fields[position_field[2]]));
				vec precisions = unit_cell.get_coords_cartesian(
					min(0.01, get_decimal_precision_from_CIF_number(fields[position_field[0]])),
					min(0.01, get_decimal_precision_from_CIF_number(fields[position_field[1]])),
					min(0.01, get_decimal_precision_from_CIF_number(fields[position_field[2]])));
				for (int i = 0; i < 3; i++)
				{
					precisions[i] = abs(precisions[i]);
				}
				if (debug)
					file << " cart. pos.: " << setw(8) << position[0] << "+/-" << precisions[0] << " " << setw(8) << position[1] << "+/-" << precisions[1] << " " << setw(8) << position[2] << "+/-" << precisions[2] << "\n";

				int group_nr = 0;
				if (group_field != -1 && fields[group_field] != "." && fields[group_field] != "?"
					&& !fields[group_field].empty()) {
					// An unreadable disorder group must not abort the run.
					try { group_nr = std::stoi(fields[group_field]); }
					catch (const std::exception &) {
						file << "Could not read disorder group for atom " << fields[label_field]
							 << "; treating it as 0.\n";
						group_nr = 0;
					}
				}
				// Filter PART before matching this CIF row to a WFN atom.  Disorder
				// positions can be close enough to match an atom in another PART;
				// doing this after the match overwrites that atom's coordinates/ID and
				// corrupts the duplicate check of a later -mtc fragment.
				if (!input_groups.empty() && group_field != -1 &&
					std::find(input_groups.begin(), input_groups.end(), group_nr) == input_groups.end())
				{
					if (debug)
						file << "Wrong part!\n";
					getline_universal(cif_input, line);
					continue;
				}
				bool old_atom = false;
				const int row_Z = constants::get_Z_from_label(fields[type_field].c_str()) + 1; //0 when the type is no plain symbol
				const atomID cif_atom_id(
					stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]),
					group_nr,
					row_Z);
				const std::string atom_ID = cif_atom_id.to_hex_string();
				//serial: a parallel region per CIF row costs a thread wake-up per row, more than the comparisons themselves
				for (int run = 0; run < known_atoms.size(); run++)
				{
					if (fields[label_field] == known_atoms[run] || atom_ID == known_atoms[run])
					{
						old_atom = true;
						if (debug)
							file << "I already know this one! " << fields[label_field] << " " << known_atoms[run] << "\n";
					}
				}
				if (old_atom)
				{
					getline_universal(cif_input, line);
					continue;
				}
				vec tolerances(3);
				int wrong_element = -1; //free wfn atom at the row's position, of another element
				for (int i = 0; i < wave.get_ncen(); i++)
				{
					for (int j = 0; j < 3; j++)
					{
						/* The cap has to be well inside a bond length. It is derived
						from the fractional precision, so it grows with the cell: on a
						55-68 A protein cell it reached 1.3 A per axis and matched
						HB_A:53 to a carbon 0.83 A away, after which the element test
						dropped the atom as a mismatch. The cif and the xyz describe
						the same model, so they agree far better than this.
						*/
						tolerances[j] = 2 * max(min(abs(precisions[j]), 0.1), 0.01);
					}
					if (is_similar_abs(position[0], wave.get_atom_coordinate(i, 0), tolerances[0]) && is_similar_abs(position[1], wave.get_atom_coordinate(i, 1), tolerances[1]) && is_similar_abs(position[2], wave.get_atom_coordinate(i, 2), tolerances[2]))
					{
						/* Each wfn atom belongs to one row, and to a row of its own element.
						A mixed site (EXYZ Se1 S1) puts two atoms on one position and the label
						test below takes "se1" for sulfur, so Se1 claimed the S atom, S1 claimed
						it again and the table got a column too many (COD 7236927). Two H 0.06 A
						apart in different PARTs did the same (2238091).
						ponytail: first free atom in wfn order, right while Olex2 writes the xyz in
						row order; pick the nearest free atom if that ever stops holding. */
						if (!labels[i].empty())
							continue;
						if (row_Z > 0 && wave.get_atom_charge(i) != row_Z)
						{
							wrong_element = i;
							continue;
						}
						wave.set_atom_frac_coords(i, { stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]) });
						wave.set_atom_group_nr(i, group_nr);
						// Store exactly the identifier used above for the MTC duplicate check.
						wave.set_id_for_atom(i, cif_atom_id);
						string element = constants::atnr2letter(wave.get_atom_charge(i));
						err_checkf(element != "PROBLEM", "Problem identifying atoms!", std::cout);
						string label = fields[label_field];
						string type = fields[type_field];
						transform(element.begin(), element.end(), element.begin(), asciitolower);
						transform(label.begin(), label.end(), label.begin(), asciitolower);
						transform(type.begin(), type.end(), type.begin(), asciitolower);
						if (debug)
						{
							file << "ASYM:  " << setw(8) << element << " charge: " << setw(17) << wave.get_atom_charge(i) << "                             wfn cart. pos: "
								<< fixed << setprecision(3) << setw(16) << wave.get_atom_coordinate(i, 0) << " "
								<< fixed << setprecision(3) << setw(16) << wave.get_atom_coordinate(i, 1) << " "
								<< fixed << setprecision(3) << setw(16) << wave.get_atom_coordinate(i, 2);
							if (input_groups.size() > 0 && group_field != -1)
							{
								file << " checking disorder group: " << fields[group_field] << " vs. ";
								for (int g = 0; g < input_groups.size(); g++)
									file << input_groups[g] << ",";
							}
						}
						if (label.find(element) == string::npos || label.find(element) > 2)
						{
							if (element != "h")
							{
								if (debug)
								{
									file << "\nElement symbol not found in label, this is a problem!\n checking type...";
									if (type.find(element) != 0)
									{
										file << " ALSO FAILED! WILL IGNORE ATOM!\n";
										continue;
									}
								}
								else
								{
									if (type.find(element) != 0)
									{
										file << "\nAtom " << label << " was not matching by element determined by label reduction or type field, skipping!\n";
										continue;
									}
								}
							}
							else if (label.find("d") == string::npos && label.find("t") == string::npos)
							{
								if (debug)
								{
									file << "\nElement symbol not found in label, this is a problem!\n will check type...";
									if (type.find(element) != 0)
									{
										file << " ALSO FAILED! WILL IGNORE ATOM!\n";
										continue;
									}
								}
								else
								{
									if (type.find(element) != 0)
									{
										file << "\nAtom " << label << " was not matching by element determined by label reduction or type field, skipping!\n";
										continue;
									}
								}
							}
						}

						labels[i] = fields[label_field];
						asym_atom_list.push_back(i);
						needs_grid[i] = true;
						nr = i;
						break;
					}
				}
				if (debug)
					file << " nr= " << nr << "\n";
				if (nr != -1)
				{
					bool already_there = false;
					for (int i = 0; i < atom_type_list.size(); i++)
						if (atom_type_list[i] == wave.get_atom_charge(nr))
						{
							already_there = true;
							asym_atom_to_type_list.push_back(i);
							break;
						}
					if (already_there == false && wave.get_atom_charge(nr) != 119)
					{
						asym_atom_to_type_list.push_back((int)atom_type_list.size());
						atom_type_list.push_back(wave.get_atom_charge(nr));
					}
				}
				else if (!old_atom)
				{
					if (wrong_element != -1)
						file << "\nCIF row " << fields[label_field] << " has type " << fields[type_field] << " but the wavefunction atom at its position is "
							<< constants::atnr2letter(wave.get_atom_charge(wrong_element)) << ", skipping!\n";
					if (debug)
					{
						file << "I did not find this atom! Tolerances were: ";
						for (int j = 0; j < 3; j++)
						{
							file << setw(12) << fixed << setprecision(8) << tolerances[j];
						}
						file << "\n";
					}
				}
				getline_universal(cif_input, line);
			}
		}
	}

	// Add missing atom types to be able to calc sphericals correctly
	for (int nr = 0; nr < wave.get_ncen(); nr++)
	{
		bool already_there = false;
		for (int i = 0; i < atom_type_list.size(); i++)
		{
			if (atom_type_list[i] == wave.get_atom_charge(nr))
			{
				already_there = true;
				break;
			}
		}
		if (already_there == false && wave.get_atom_charge(nr) != 119)
		{
			atom_type_list.push_back(wave.get_atom_charge(nr));
		}
	}

	err_checkf(asym_atom_list.size() <= wave.get_ncen(), "More asymmetric unit atoms detected than in the wavefunction! Aborting!", file);
	err_checkf(allow_empty || asym_atom_list.size() != 0, "0 asym atoms is imposible! something is wrong with reading the CIF!", file);

	for (int i = 0; i < atom_type_list.size(); i++)
		err_checkf((atom_type_list[i] <= 113 || atom_type_list[i] == 119) && atom_type_list[i] > 0, "Unreasonable atom type detected: " + toString(atom_type_list[i]) + " (Happens if Atoms were not identified correctly)", file);
	// the "Reading: <cif> done!" line is closed by the cell ctor (cell.h)
	if (debug)
	{
		file << "There are " << atom_type_list.size() << " types of atoms" << endl;
		for (int i = 0; i < atom_type_list.size(); i++)
			file << setw(4) << atom_type_list[i];
		file << endl
			<< "asym_atoms_to_type_list: " << endl;
		for (int i = 0; i < asym_atom_to_type_list.size(); i++)
			file << setw(4) << asym_atom_to_type_list[i];
		file << endl;
		file << "Charges of atoms:" << endl;
		for (int i = 0; i < wave.get_ncen(); i++)
			file << setw(4) << wave.get_atom_charge(i);
		file << endl;
	}
	int size = static_cast<int>(asym_atom_list.size());
	svec labels2;
	for (int i = 0; i < size; i++)
		labels2.emplace_back(labels[asym_atom_list[i]]);
	return labels2;
}

/**
 * Reads atoms from a CIF file and performs necessary operations. Works without wavefunction for XCW routine.
 */
void read_atoms_from_CIF(std::ifstream& cif_input, const cell& unit_cell, int& ncen, bvec& needs_grid, std::vector<asym_atom>& asym_atoms, const bool debug){
	if (!cif_input)
	{
		throw std::runtime_error("Could not open CIF file.");
	}
	std::string line;
	bool in_atom_loop = false;
	bool reading_headers = false;
	std::vector<std::string> headers;
	// Column indices
	int idx_label = -1;
	int idx_type = -1;
	int idx_x = -1;
	int idx_y = -1;
	int idx_z = -1;
	while (getline_universal(cif_input, line))
	{
		if (line.empty())
			continue;
		// Trim both ends: the header comparisons below are exact, and a CIF from Windows (Olex2) ends lines with \r.
		std::string trimmed = trim(line);
		if (trimmed.empty())
			continue;
		// Detect loop_
		if (trimmed == "loop_")
		{
			in_atom_loop = false;
			reading_headers = true;
			headers.clear();
			idx_label = -1;
			idx_type = -1;
			idx_x = -1;
			idx_y = -1;
			idx_z = -1;
			continue;
		}
		// Read headers after loop_
		if (reading_headers && !trimmed.empty() && trimmed[0] == '_')
		{
			headers.push_back(trimmed);
			continue;
		}
		// End of headers -> determine whether this is the atom loop
		if (reading_headers)
		{
			reading_headers = false;
			for (size_t i = 0; i < headers.size(); ++i)
			{
				if (headers[i] == "_atom_site_label")
					idx_label = static_cast<int>(i);
				else if (headers[i] == "_atom_site_type_symbol")
					idx_type = static_cast<int>(i);
				else if (headers[i] == "_atom_site_fract_x")
					idx_x = static_cast<int>(i);
				else if (headers[i] == "_atom_site_fract_y")
					idx_y = static_cast<int>(i);
				else if (headers[i] == "_atom_site_fract_z")
					idx_z = static_cast<int>(i);
			}
			// Check if this is the desired atom loop
			if (idx_label >= 0 &&
				idx_type >= 0 &&
				idx_x >= 0 &&
				idx_y >= 0 &&
				idx_z >= 0)
			{
				in_atom_loop = true;
				if (debug)
				{
					std::cout << "Found atom_site loop\n";
					std::cout << "label idx = " << idx_label << "\n";
					std::cout << "type  idx = " << idx_type << "\n";
					std::cout << "x     idx = " << idx_x << "\n";
					std::cout << "y     idx = " << idx_y << "\n";
					std::cout << "z     idx = " << idx_z << "\n";
				}
			}
		}
		// Read atom rows
		if (in_atom_loop)
		{
			// Another loop or new data item ends this loop
			if (trimmed == "loop_" || trimmed[0] == '_')
			{
				in_atom_loop = false;
				// Reprocess line in outer loop logic
				if (trimmed == "loop_")
				{
					reading_headers = true;
					headers.clear();
				}
				continue;
			}
			std::istringstream iss(trimmed);
			std::vector<std::string> tokens;
			std::string token;
			while (iss >> token)
				tokens.push_back(token);
			const int required_cols =
				std::max({ idx_label, idx_type, idx_x, idx_y, idx_z }) + 1;
			if (static_cast<int>(tokens.size()) < required_cols)
			{
				if (debug)
				{
					std::cout << "Skipping malformed atom line:\n";
					std::cout << trimmed << "\n";
				}
				continue;
			}
			asym_atom temp_atom;
			temp_atom.label = tokens[idx_label];
			const std::string& type_str = tokens[idx_type];
			temp_atom.type = constants::get_Z_from_label(type_str.c_str()) + 1;
			const double fx = std::stod(tokens[idx_x]);
			const double fy = std::stod(tokens[idx_y]);
			const double fz = std::stod(tokens[idx_z]);
			temp_atom.frac_pos = { fx, fy, fz };
			auto cart =
				unit_cell.get_coords_cartesian(fx, fy, fz, true);
			temp_atom.pos = { cart[0], cart[1], cart[2] };
			asym_atoms.push_back(temp_atom);
			if (debug)
			{
				std::cout << "Parsed atom:\n";
				std::cout << "  label = " << temp_atom.label << "\n";
				std::cout << "  type  = " << temp_atom.type << "\n";
				std::cout << "  frac  = ("
					<< fx << ", "
					<< fy << ", "
					<< fz << ")\n";
			}
		}
	}
	ncen = static_cast<int>(asym_atoms.size());
	needs_grid.resize(ncen, true);
	if (debug)
	{
		std::cout << "\nTotal atoms parsed: " << ncen << "\n";
	}
}
