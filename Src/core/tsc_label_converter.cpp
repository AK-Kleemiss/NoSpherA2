#include "pch.h"
#include "tsc_label_converter.h"

#include "constants.h"
#include "tsc_block.h"

#include <cctype>

namespace
{
struct cif_atom { atomID id; std::string label; };

std::string uppercase(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
        [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return value;
}

std::string without_esd(std::string value)
{
    if (const auto pos = value.find('('); pos != std::string::npos)
        value.resize(pos);
    return value;
}

std::vector<std::string> cif_tokens(const std::string& line)
{
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (input >> std::quoted(value)) {
        if (value.size() >= 2 && value.front() == '\'' && value.back() == '\'')
            value = value.substr(1, value.size() - 2);
        result.push_back(value);
    }
    return result;
}

std::vector<cif_atom> read_cif_atoms(const std::filesystem::path& cif_file)
{
    std::ifstream input(cif_file);
    if (!input) throw std::runtime_error("Failed to open CIF file: " + cif_file.string());
    std::string line;
    while (std::getline(input, line)) {
        if (trim(line) != "loop_") continue;
        std::vector<std::string> columns;
        while (std::getline(input, line) && trim(line).starts_with('_'))
            columns.push_back(uppercase(trim(line)));
        const auto column = [&columns](const std::string& name) {
            const auto found = std::find(columns.begin(), columns.end(), name);
            return found == columns.end() ? -1 : static_cast<int>(found - columns.begin());
        };
        const int label = column("_ATOM_SITE_LABEL"), type = column("_ATOM_SITE_TYPE_SYMBOL");
        const int x = column("_ATOM_SITE_FRACT_X"), y = column("_ATOM_SITE_FRACT_Y");
        const int z = column("_ATOM_SITE_FRACT_Z"), group = column("_ATOM_SITE_DISORDER_GROUP");
        if (label < 0 || type < 0 || x < 0 || y < 0 || z < 0) continue;
        std::vector<cif_atom> result;
        do {
            const auto fields = cif_tokens(line);
            if (fields.size() != columns.size()) break;
            const int part = group < 0 || fields[group] == "." || fields[group] == "?" ? 0 : std::stoi(fields[group]);
            const int atomic_number = constants::get_Z_from_label(fields[type].c_str()) + 1;
            if (atomic_number <= 0) throw std::runtime_error("Unknown atom type in CIF: " + fields[type]);
            result.push_back({ atomID(std::stod(without_esd(fields[x])), std::stod(without_esd(fields[y])),
                std::stod(without_esd(fields[z])), part, atomic_number), fields[label] });
        } while (std::getline(input, line) && !trim(line).empty());
        if (!result.empty()) return result;
    }
    throw std::runtime_error("CIF contains no usable _atom_site loop");
}

tsc_block<int, cdouble> read_text_tsc(const std::filesystem::path& table_file)
{
    std::ifstream input(table_file);
    if (!input) throw std::runtime_error("Failed to open TSC file: " + table_file.string());
    ScattererLabels scatterers;
    ivec2 indices(3);
    cvec2 form_factors;
    const auto header = tsc_merge_detail::read_header(input, scatterers, true);
    form_factors.resize(scatterers.size());
    tsc_merge_detail::read_unchecked_data(input, 0, header, indices, form_factors);
    return tsc_block<int, cdouble>(form_factors, scatterers, indices, header.header);
}

tsc_block<int, cdouble> read_tsc(const std::filesystem::path& table_file)
{
    const std::string extension = uppercase(table_file.extension().string());
    if (extension == ".TSCB") return tsc_block<int, cdouble>(table_file);
    if (extension == ".TSC") return read_text_tsc(table_file);
    throw std::runtime_error("Table file must end in .tsc or .tscb");
}
}

bool convert_tsc_ids_to_labels(const std::filesystem::path& table_file,
    const std::filesystem::path& cif_file, const std::filesystem::path& output_file,
    std::ostream& log)
{
    try {
        const tsc_block<int, cdouble> input = read_tsc(table_file);
        std::unordered_map<atomID, std::string> labels;
        for (const auto& atom : read_cif_atoms(cif_file)) {
            if (!labels.emplace(atom.id, atom.label).second)
                throw std::runtime_error("CIF contains duplicate SCATTERER_IDS");
        }
        std::vector<std::string> table_labels;
        cvec2 form_factors;
        table_labels.reserve(input.scatterer_size());
        form_factors.reserve(input.scatterer_size());
        for (std::size_t i = 0; i < input.scatterer_size(); ++i) {
            const auto& scatterer = input.get_scatterer(i);
            if (!std::holds_alternative<atomID>(scatterer))
                throw std::runtime_error("Input table does not use SCATTERER_IDS");
            const atomID& id = std::get<atomID>(scatterer);
            const auto match = labels.find(id);
            if (match == labels.end())
                throw std::runtime_error("No CIF atom matches SCATTERER_ID " + id.to_hex_string());
            table_labels.push_back(match->second);
            form_factors.push_back(input.get_sf_for_scatterer(i));
        }
        tsc_block<int, cdouble> output(form_factors, table_labels, input.get_index_vector());
        output.set_AD(input.get_AD());
        output.write_tsc_file(cif_file, output_file);
        log << "Matched " << table_labels.size() << " SCATTERER_IDS and wrote " << output_file << '\n';
        return true;
    } catch (const std::exception& error) {
        log << "Could not convert SCATTERER_IDS to labels: " << error.what() << '\n';
        return false;
    }
}
