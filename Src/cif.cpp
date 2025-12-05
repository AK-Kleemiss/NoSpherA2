#include "convenience.h"
#include "constants.h"
#include "atoms.h"
#include "wfn_class.h"
#include "tsc_block.h"

const std::string get_basis_loop_header() {
    return "loop_\n_basis.id\n_basis.name\n_basis.dict\n";
}

const std::string get_wavefunction_loop_header() {
    return "loop_\n_wavefunction.id\n_wavefunction.type\n_wavefunction.radial_type\n_wavefunction.angular_type\n_wavefunction.dict\n";
}

const std::string get_basis_set_CIF(const int nr, const WFN& wavy)
{
    // Make list of unique atom types:
    ivec atom_types;
    ivec atoms_with_type;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        if (find(atom_types.begin(), atom_types.end(), wavy.get_atom_charge(i)) == atom_types.end())
        {
            atom_types.push_back(wavy.get_atom_charge(i));
            atoms_with_type.push_back(i);
        }
    }
    int _nr;
    if (nr == 0)
        _nr = 1;
    else
        _nr = nr;
    std::stringstream ss;
    ss << _nr << " '" << wavy.get_basis_set_name() << "' [\n";
    for (int i = 0; i < atom_types.size(); i++)
    {
        ss << "  {\n";
        ss << "    'atom_site_label': '" << wavy.get_atom_label(i) << "'\n";
        ss << "    'Z': " << atom_types[i] << "\n";
        ss << "    'atom_type': " << constants::atnr2letter(atom_types[i]) << "\n";
        ss << "    'nr_shells': " << wavy.get_atom_shell_count(atoms_with_type[i]) << "\n";
        ss << "    'shell_sizes': [";
        for (int j = 0; j < wavy.get_atom_shell_count(atoms_with_type[i]); j++)
        {
            ss << wavy.get_atom_shell_primitives(atoms_with_type[i], j);
        }
        ss << "]\n";
        ss << "    'shell_types': [";
        for (int j = 0; j < wavy.get_atom_shell_count(atoms_with_type[i]); j++)
        {
            ss << wavy.get_shell_type(atoms_with_type[i], j);
        }
        ss << "]\n";
        ss << "    'exponent_unit': 'a.u.'\n";
        ss << "    'primitive_exponents': [";
        for (unsigned int j = 0; j < wavy.get_atom_basis_set_size(atoms_with_type[i]); j++)
        {
            ss << wavy.get_atom_basis_set_exponent(atoms_with_type[i], j);
            if (j < wavy.get_atom_basis_set_size(atoms_with_type[i]) - 1)
            {
                ss << " ";
            }
        }
        ss << "]\n";
        ss << "    'primitive_coefficients': [";
        for (unsigned int j = 0; j < wavy.get_atom_basis_set_size(atoms_with_type[i]); j++)
        {
            ss << wavy.get_atom_basis_set_coefficient(atoms_with_type[i], j);
            if (j < wavy.get_atom_basis_set_size(atoms_with_type[i]) - 1)
            {
                ss << " ";
            }
        }
        ss << "]\n";
        ss << "  }\n";
    }
    ss << "]\n";
    return ss.str();
}

const std::string get_CIF_table(const int nr, const WFN& wavy)
{
    std::stringstream ss;
    int _nr;
    if (nr == 0)
        _nr = 1;
    else
        _nr = nr;
    ss << _nr << " 'Molecular' 'GTO' 'Cartesian' "
        << "{\n";
    ss << "  'atoms': [\n";
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        ss << "    {\n";
        ss << "      'id': " << i << "\n";
        ss << "      'atom_site_label': '" << wavy.get_atom_label(i) << "'\n";
        ss << "      'cartesian_position': [" << wavy.get_atom_coordinate(i, 0) << " " << wavy.get_atom_coordinate(i, 1) << " " << wavy.get_atom_coordinate(i, 2) << "]\n";
        ss << "      'sym_code': '.'\n";
        ss << "      'Z': " << wavy.get_atom_charge(i) << "\n";
        ss << "      'basis_set_id': " << wavy.get_atom_basis_set_id(i) << "\n";
        ss << "    }\n";
    }
    ss << "  ]\n";
    ss << "  'MOs': {\n";
    ss << "    'spins': [";
    for (int i = 0; i < wavy.get_nmo(); i++)
    {
        int spin = wavy.get_MO_op(i);
        if (spin == 0)
        {
            ss << "alpha";
        }
        else if (spin == 1)
        {
            ss << "beta";
        }
        else
        {
            ss << "unknown";
        }
        if (i < wavy.get_nmo() - 1)
        {
            ss << " ";
        }
    }
    ss << "]\n";
    ss << "    'energies': [";
    for (int i = 0; i < wavy.get_nmo(); i++)
    {
        ss << wavy.get_MO_energy(i);
        if (i < wavy.get_nmo() - 1)
        {
            ss << " ";
        }
    }
    ss << "]\n";
    ss << "    'occupancies': [";
    for (int i = 0; i < wavy.get_nmo(); i++)
    {
        ss << wavy.get_MO_occ(i);
        if (i < wavy.get_nmo() - 1)
        {
            ss << " ";
        }
    }
    ss << "]\n";
    ss << "    'coefficients': [\n";
    for (int i = 0; i < wavy.get_nmo(); i++)
    {
        ss << "      [";
        for (int j = 0; j < wavy.get_nex(); j++)
        {
            ss << wavy.get_MO_coef_f(i, j);
            if (j < wavy.get_nex() - 1)
            {
                ss << " ";
            }
        }
        ss << "]";
        if (i < wavy.get_nmo() - 1)
        {
            ss << "\n";
        }
    }
    ss << "    ]\n";
    ss << "  }\n";
    ss << "}\n";
    return ss.str();
}

void write_wfn_CIF(WFN& wavy, const std::filesystem::path& fileName)
{
    err_checkf(wavy.get_basis_set_name() != " ", "Please load a basis set before writing things to a .cif file!", std::cout);
    std::ofstream file(fileName);
    std::stringstream ss;
    ss << get_basis_loop_header();
    ss << get_basis_set_CIF(0, wavy);
    ss << "\n\n";
    ss << get_wavefunction_loop_header();
    ss << get_CIF_table(0, wavy);
    file << ss.str();
    file.close();
}

void write_wfn_CIF(WFN& wavy, const std::filesystem::path& fileName, tsc_block<int, cdouble>& tsc, options& opt)
{
    err_checkf(wavy.get_basis_set_name() != " ", "Please load a basis set before writing things to a .cif file!", std::cout);
    std::ofstream file(fileName);
    std::stringstream ss;
    ss << get_basis_loop_header();
    ss << get_basis_set_CIF(0, wavy);
    ss << "\n\n";
    ss << get_wavefunction_loop_header();
    ss << get_CIF_table(0, wavy);
    if (!tsc.is_empty()) {
        ss << "\n\n";
        ss << tsc.get_tsc_cif_block(opt);
    }
    file << ss.str();
    file.close();
}

void write_wfn_CIF(std::vector<WFN>& wavy, const std::filesystem::path& fileName, tsc_block<int, cdouble>& tsc, options& opt)
{
    std::ofstream file(fileName);
    std::stringstream ss;
    ss << get_basis_loop_header();
    for (int i = 0; i < wavy.size(); i++) {
        err_checkf(wavy[i].get_basis_set_name() != " ", "Please load a basis set before writing things to a .cif file!", std::cout);
        ss << get_basis_set_CIF(i+1, wavy[i]);
    }
    ss << "\n\n";
    ss << get_wavefunction_loop_header();
    for (int i = 0; i < wavy.size(); i++) {
        ss << get_CIF_table(i+1, wavy[i]);
    }
    if (!tsc.is_empty()) {
        ss << "\n\n";
        ss << tsc.get_tsc_cif_block(opt);
    }
    file << ss.str();
    file.close();
}