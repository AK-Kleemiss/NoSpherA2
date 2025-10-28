#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "mo_class.h"
#include "cube.h"
#include "constants.h"
#include "fchk.h"
#include "integrator.h"
#include "basis_set.h"
#include "nos_math.h"
#include "libCintMain.h"
#include "integration_params.h"

void WFN::fill_pre()
{
    for (int j = 0; j < 9; j++)
        for (int l = 0; l < 5; l++)
            for (int m = 0; m < 5; m++)
            {
                int imax = std::min(j, l);
                int imin = std::max(0, j - m);
                for (int i = imin; i <= imax; i++)
                    pre[j][l][m][i] = constants::ft[j] * constants::ft[l] / constants::ft[l - i] / constants::ft[i] * constants::ft[m] / constants::ft[m - j + i] / constants::ft[j - i];
            }
}

void WFN::fill_Afac_pre()
{
    for (int l = 0; l < 9; l++)
        for (int r = 0; r <= l / 2; r++)
            for (int s = 0; s <= (l - 2 * r) / 2; s++)
                Afac_pre[l][r][s] = constants::ft[r] * constants::ft[s] * constants::ft[l - 2 * r - 2 * s];
}

WFN::WFN()
{
    ncen = 0;
    nfunc = 0;
    nmo = 0;
    nex = 0;
    charge = 0;
    multi = 0;
    origin = e_origin::NOT_YET_DEFINED;
    ECP_m = 0;
    total_energy = 0.0;
    virial_ratio = 0.0;
    d_f_switch = false;
    modified = false;
    distance_switch = false;
    basis_set_name = " ";
    has_ECPs = false;
    comment = "Test";
    basis_set = NULL;
    fill_pre();
    fill_Afac_pre();
};

WFN::WFN(e_origin given_origin)
{
    ncen = 0;
    nfunc = 0;
    nmo = 0;
    nex = 0;
    charge = 0;
    multi = 0;
    ECP_m = 0;
    total_energy = 0.0;
    origin = given_origin;
    d_f_switch = false;
    modified = false;
    distance_switch = false;
    has_ECPs = false;
    basis_set_name = " ";
    comment = "Test";
    basis_set = NULL;
    fill_pre();
    fill_Afac_pre();
};

WFN::WFN(const std::filesystem::path & filename, const bool& debug)
{
    ncen = 0;
    nfunc = 0;
    nmo = 0;
    nex = 0;
    charge = 0;
    multi = 0;
    ECP_m = 0;
    total_energy = 0.0;
    d_f_switch = false;
    modified = false;
    distance_switch = false;
    has_ECPs = false;
    basis_set_name = " ";
    comment = "Test";
    basis_set = NULL;
    fill_pre();
    fill_Afac_pre();
    read_known_wavefunction_format(filename, std::cout, debug);
};

WFN::WFN(const std::filesystem::path& filename, const int g_charge, const int g_mult, const bool& debug) {
    ncen = 0;
    nfunc = 0;
    nmo = 0;
    nex = 0;
    charge = g_charge;
    multi = g_mult;
    ECP_m = 0;
    total_energy = 0.0;
    d_f_switch = false;
    modified = false;
    distance_switch = false;
    has_ECPs = false;
    basis_set_name = " ";
    comment = "Test";
    basis_set = NULL;
    fill_pre();
    fill_Afac_pre();
    read_known_wavefunction_format(filename, std::cout, debug);
};

bool WFN::push_back_atom(const std::string &label, const double &x, const double &y, const double &z, const int &_charge, const std::string& ID)
{
    ncen++;
    if (_charge >= 1)
        atoms.emplace_back(label, ID, ncen, x, y, z, _charge);
    else
    {
        atoms.emplace_back();
        return false;
    }
    return true;
};

bool WFN::push_back_atom(const atom &given)
{
    ncen++;
    atoms.push_back(given);
    return true;
};

bool WFN::erase_atom(const int &nr)
{
    err_checkf(nr < ncen, "unreasonable atom number", std::cout);
    removeElement(atoms, atoms[nr]);
    ncen--;
    return true;
};

bool WFN::push_back_MO(const int &nr, const double &occ, const double &ener)
{
    nmo++;
    err_checkf(nr <= nmo, "unreasonable MO number", std::cout);
    MOs.push_back(MO(nr, occ, ener));
    return true;
};

bool WFN::push_back_MO(const int &nr, const double &occ, const double &ener, const int &oper)
{
    nmo++;
    MOs.push_back(MO(nr, occ, ener, oper));
    return true;
};

bool WFN::push_back_MO(const MO &given)
{
    nmo++;
    MOs.push_back(given);
    return true;
};

void WFN::push_back_MO_coef(const int &nr, const double &value)
{
    err_checkf(nr < nmo, "not enough MOs", std::cout);
    MOs[nr].push_back_coef(value);
};

void WFN::assign_MO_coefs(const int &nr, vec &values)
{
    err_checkf(nr < nmo, "not enough MOs", std::cout);
    MOs[nr].assign_coefs(values);
};

const double& WFN::get_MO_energy(const int &mo) const
{
    err_checkf(mo < nmo, "not enough MOs", std::cout);
    return MOs[mo].get_energy();
}

const void WFN::clear_MOs()
{
    MOs.clear();
    MOs.shrink_to_fit();
    nmo = 0;
}

bool WFN::push_back_center(const int &cent)
{
    if (cent <= ncen && cent > 0)
        centers.push_back(cent);
    else
        return false;
    return true;
};

bool WFN::erase_center(const int &g_nr)
{
    centers.erase(centers.begin() + g_nr - 1);
    return true;
};

const std::string WFN::get_centers(const bool &bohr) const
{
    std::string temp;
    for (int i = 0; i < ncen; i++)
    {
        temp.append(atoms[i].get_label());
        temp.append(" ");
        if (bohr)
            temp.append(std::to_string(get_atom_coordinate(i,0)));
        else
            temp.append(std::to_string(constants::bohr2ang(get_atom_coordinate(i,0))));
        temp.append(" ");
        if (bohr)
            temp.append(std::to_string(get_atom_coordinate(i,1)));
        else
            temp.append(std::to_string(constants::bohr2ang(get_atom_coordinate(i,1))));
        temp.append(" ");
        if (bohr)
            temp.append(std::to_string(get_atom_coordinate(i,2)));
        else
            temp.append(std::to_string(constants::bohr2ang(get_atom_coordinate(i,2))));
        temp.append("\n");
    }
    return temp;
};

const void WFN::list_centers() const
{
    for (int i = 0; i < ncen; i++)
    {
        std::cout << atoms[i].get_nr() << " " << atoms[i].get_label() << " "
             << get_atom_coordinate(i,0) << " " << get_atom_coordinate(i,1) << " "
             << get_atom_coordinate(i,2) << " " << get_atom_charge(i) << std::endl;
    }
};

const MO &WFN::get_MO(const int &n) const
{
    if (n < nmo)
        return MOs[n];
    else
    {
        err_not_impl_f("Wrong MO number", std::cout);
        return MOs[0];
    }
}
const int WFN::get_MO_op_count(const int &op) const
{
    int count = 0;
#pragma omp parallel for reduction(+ : count)
    for (int i = 0; i < nmo; i++)
        if (MOs[i].get_op() == op)
            count++;
    return count;
};

void WFN::delete_MO(const int &nr)
{
    err_checkf(nr < nmo, "not enough MOs", std::cout);
    MOs.erase(MOs.begin() + nr);
    nmo--;
};

bool WFN::push_back_type(const int &type)
{
    types.push_back(type);
    return true;
};

bool WFN::erase_type(const int &nr)
{
    err_checkf(nr >= 1, "Wrong type to erase!", std::cout);
    types.erase(types.begin() + nr - 1);
    return true;
};

bool WFN::push_back_exponent(const double &e)
{
    exponents.push_back(e);
    return true;
};

bool WFN::erase_exponent(const int &nr)
{
    if (nr < 1)
        return false;
    exponents.erase(exponents.begin() + (nr - 1));
    return true;
};

bool WFN::remove_primitive(const int &nr)
{
    nex--;
    if (erase_center(nr) && erase_exponent(nr) && erase_type(nr))
    {
        for (int n = 0; n < nmo; n++)
            MOs[n].erase_coef(nr, nex);
        return true;
    }
    else
        return false;
};

bool WFN::add_primitive(const int &cent, const int &type, const double &e, double *values)
{
    nex++;
    if (push_back_center(cent) && push_back_type(type) && push_back_exponent(e))
        for (int n = 0; n < nmo; n++)
            MOs[n].push_back_coef(values[n]);
    else
        return false;
    return true;
};

void WFN::change_type(const int &nr)
{
    err_checkf(nr < nex, "Wrong input", std::cout);
    bool end = false;
    while (!end)
    {
        std::cout << "Please enter the new type you want to assign: ";
        int new_type = 0;
        std::cin >> new_type;
        if (new_type > ncen || new_type < 0)
        {
            std::cout << "Sorry, wrong input, try again!\n";
            continue;
        }
        types[nr - 1] = new_type;
        end = true;
        set_modified();
    }
};

void WFN::change_exponent(const int &nr)
{
    err_checkf(nr < nex, "Wrong input", std::cout);
    bool end = false;
    while (!end)
    {
        std::cout << "Please enter the new exponent you want to assign: ";
        int new_exp = 0;
        std::cin >> new_exp;
        if (new_exp > ncen || new_exp < 0)
        {
            std::cout << "Sorry, wrong input, try again!\n";
            continue;
        }
        exponents[nr - 1] = new_exp;
        end = true;
        set_modified();
    }
};

void WFN::change_center(const int &nr)
{
    bool end = false;
    while (!end)
    {
        std::cout << "Please enter the new center you want to assign: ";
        int new_center = 0;
        std::cin >> new_center;
        if (new_center > ncen || new_center < 0)
        {
            std::cout << "Sorry, wrong input, try again!\n";
            continue;
        }
        centers[nr - 1] = new_center;
        end = true;
        set_modified();
    }
};

bool WFN::set_MO_coef(const int &nr_mo, const int &nr_primitive, const double &value)
{
    err_checkf(nr_mo <= MOs.size(), "MO doesn't exist!", std::cout);
    return MOs[nr_mo].set_coefficient(nr_primitive, value);
};

const void WFN::list_primitives() const
{
    for (int i = 0; i < nex; i++)
    {
        std::cout << i << " center: " << centers[i] << " type: " << types[i] << " exponent: " << exponents[i] << std::endl;
    }
};

bool WFN::remove_center(const int &nr)
{
    erase_center(nr);
    try
    {
        for (int i = 0; i < nex; i++)
            if (centers[i] == nr)
                remove_primitive(i);
    }
    catch (...)
    {
        err_checkf(false, "Problem removing center", std::cout);
        return false;
    }
    set_modified();
    return true;
}

bool WFN::add_exp(const int &cent, const int &type, const double &e)
{
    nex++;
    if (!push_back_center(cent) || !push_back_type(type) || !push_back_exponent(e))
        return false;
    else
        return true;
};

const double& WFN::get_MO_coef(const int &nr_mo, const int &nr_primitive) const
{
    err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
    return MOs[nr_mo].get_coefficient(nr_primitive);
};

const double& WFN::get_MO_coef_f(const int &nr_mo, const int &nr_primitive) const
{
    err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
    return MOs[nr_mo].get_coefficient_f(nr_primitive);
};

const double* WFN::get_MO_coef_ptr(const int &nr_mo)
{
    err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
    return MOs[nr_mo].get_coefficient_ptr();
};

const int WFN::get_MO_primitive_count(const int &nr_mo) const
{
    err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
    return MOs[nr_mo].get_primitive_count();
};

const std::string WFN::hdr(const bool &occupied) const
{
    std::string temp = "GAUSSIAN            ";
    if (!occupied)
    {
        if (nmo > 100)
            temp.append(std::to_string(nmo));
        else if (nmo < 100 && nmo > 10)
        {
            temp.append(" ");
            temp.append(std::to_string(nmo));
        }
        else if (nmo < 10 && nmo > 0)
        {
            temp.append("  ");
            temp.append(std::to_string(nmo));
        }
        else if (nmo < 0)
            return "PROBLEM";
    }
    else
    {
        const int occupied_mos = get_nmo(true);
        if (occupied_mos > 100)
            temp.append(std::to_string(occupied_mos));
        else if (occupied_mos < 100 && occupied_mos > 10)
        {
            temp.append(" ");
            temp.append(std::to_string(occupied_mos));
        }
        else if (occupied_mos < 10 && occupied_mos > 0)
        {
            temp.append("  ");
            temp.append(std::to_string(occupied_mos));
        }
        else if (occupied_mos < 0)
            return "PROBLEM";
    }
    temp.append(" MOL ORBITALS    ");
    if (nex > 100)
        temp.append(std::to_string(nex));
    else if (nex < 100 && nex > 10)
    {
        temp.append(" ");
        temp.append(std::to_string(nex));
    }
    else if (nex < 10 && nex > 0)
    {
        temp.append("  ");
        temp.append(std::to_string(nex));
    }
    else if (nex < 0)
        return "PROBLEM";
    temp.append(" PRIMITIVES      ");
    if (ncen > 100)
        temp.append(std::to_string(ncen));
    else if (ncen < 100 && ncen > 10)
    {
        temp.append(" ");
        temp.append(std::to_string(ncen));
    }
    else if (ncen < 10 && ncen > 0)
    {
        temp.append("  ");
        temp.append(std::to_string(ncen));
    }
    else if (ncen < 0)
        return "PROBLEM";
    temp.append(" NUCLEI\n");
    return temp;
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
        // file << "There is already a wavefunction loaded, do you want to continue and possibly overwrite the existing wavefunction?" << endl;
        // if (!yesno()) return false;
        // else file << "okay, carrying on then..." << endl;
        file << "There is already a wavefunction loaded, aborting!" << endl;
        return false;
    }
    origin = e_origin::wfn;
    err_checkf(std::filesystem::exists(fileName), "Couldn't open or find " + fileName.string() + ", leaving", file);
    ifstream rf(fileName);
    if (rf.good())
        path = fileName;
    string line;
    rf.seekg(0);
    getline(rf, line);
    comment = line;
    getline(rf, line);
    stringstream stream(line);
    string header_tmp;
    int e_nmo, e_nex, e_nuc = 0; // number of expected MOs, Exponents and nuclei
    stream >> header_tmp >> e_nmo >> header_tmp >> header_tmp >> e_nex >> header_tmp >> e_nuc;
    if (debug)
    {
        file << "e_nmo: " << e_nmo << ", e_nex: " << e_nex << ", e_nuc : " << e_nuc << endl;
    }
    //----------------------------- Read Atoms ------------------------------------------------------------
    ivec dum_nr, dum_ch;
    dum_nr.resize(e_nuc);
    dum_ch.resize(e_nuc);
    svec dum_label;
    vec dum_x, dum_y, dum_z;
    dum_x.resize(e_nuc);
    dum_y.resize(e_nuc);
    dum_z.resize(e_nuc);
    char tempchar[20];
    size_t length;
    for (int i = 0; i < e_nuc; i++)
    {
        // int dump = 0;
        getline(rf, line);
        if (debug)
            file << i << ".run, line:" << line << endl;
        length = line.copy(tempchar, 4, 0);
        tempchar[length] = '\0';
        string temp;
        temp = tempchar;
        length = line.copy(tempchar, 4, 5);
        tempchar[length] = '\0';
        dum_nr[i] = i;
        length = line.copy(tempchar, 12, 24);
        tempchar[length] = '\0';
        dum_x[i] = stod(tempchar);
        length = line.copy(tempchar, 12, 36);
        tempchar[length] = '\0';
        dum_y[i] = stod(tempchar);
        length = line.copy(tempchar, 12, 48);
        tempchar[length] = '\0';
        dum_z[i] = stod(tempchar);
        length = line.copy(tempchar, 3, 70);
        tempchar[length] = '\0';
        dum_ch[i] = stoi(tempchar);
        dum_label.push_back(shrink_string_to_atom(temp, dum_ch[i]));
    }
    //------------------------------------ Read center assignements -------------------------------------------
    ivec dum_center;
    dum_center.resize(e_nex);
    if (debug)
        for (int i = 0; i < e_nex; i++)
            dum_center[i] = 99;
    // int run = 0;
    int exnum = 0;
    // int dump = 0;
    getline(rf, line);
    while (line.compare(0, 6, "CENTRE") == 0 && !rf.eof())
    {
        if (exnum + 20 <= e_nex)
        {
            for (int i = 0; i < 20; i++)
            {
                length = line.copy(tempchar, 3, 20 + 3 * i);
                tempchar[length] = '\0';
                dum_center[exnum] = stoi(tempchar);
                if (dum_center[exnum] > e_nuc)
                {
                   std::cout << "this center doesn't exist.. some weird problem!\n";
                    return false;
                }
                exnum++;
            }
        }
        else
        {
            if (exnum < e_nex)
            {
                for (int i = 0; i < e_nex % 20; i++)
                {
                    length = line.copy(tempchar, 3, 20 + 3 * i);
                    tempchar[length] = '\0';
                    dum_center[exnum] = stoi(tempchar);
                    if (dum_center[exnum] > e_nuc)
                    {
                        file << "this center doesn't exist.. some weird problem!\n";
                        return false;
                    }
                    exnum++;
                }
            }
            else
            {
                getline(rf, line);
                continue;
            }
        }
        getline(rf, line);
        if (exnum > e_nex)
        {
            file << "run went higher than expected values in center reading, thats suspicius, lets stop here...\n";
            return false;
        }
        // run++;
    }
    if (debug)
    {
        file << exnum << endl;
    }
    if (exnum < e_nex)
    {
        file << "We have a problem adding center assignements!\n";
        return false;
    }
    if (debug)
        file << "finished with centers, moving to types...\n";
    //------------------------------------ Read Types ---------------------------------------------------------
    vector<unsigned int> dum_type;
    dum_type.resize(e_nex);
    // run = 0;
    exnum = 0;
    while (line.compare(0, 4, "TYPE") == 0 && !rf.eof())
    {
        if (exnum + 20 <= e_nex)
        {
            for (int i = 0; i < 20; i++)
            {
                length = line.copy(tempchar, 2, size_t(21 + 3 * i));
                tempchar[length] = '\0';
                dum_type[exnum] = stoi(tempchar);
                exnum++;
            }
        }
        else if (exnum < e_nex)
        {
            for (int i = 0; i < e_nex % 20; i++)
            {
                length = line.copy(tempchar, 2, 21 + 3 * i);
                tempchar[length] = '\0';
                dum_type[exnum] = stoi(tempchar);
                exnum++;
            }
        }
        else
        {
            getline(rf, line);
            continue;
        }
        getline(rf, line);
        if (exnum > e_nex)
        {
            file << "exnum went higher than expected values in type reading, thats suspicius, lets stop here...\n";
            return false;
        }
        // run++;
    }
    if (exnum < e_nex)
    {
        file << "We have a problem adding type assignements!\n";
        return false;
    }
    if (debug)
        file << "finished with types, reading exponents now...\n";
    //----------------------------- Read exponents -------------------------------
    vec dum_exp;
    dum_exp.resize(e_nex);
    // run = 0;
    exnum = 0;
    string replace = "E";
    bool three_exponents = false;
    while (line.compare(0, 9, "EXPONENTS") == 0 && !rf.eof())
    {
        if (exnum + 5 <= e_nex)
        {
            if (exnum == 0)
            {
                const char test = line.at(10);
                string temp_str(" ");
                const char empty = temp_str.at(0);
                if (test != empty)
                    three_exponents = true;
            }
            for (int i = 0; i < 5; i++)
            {

                if (!three_exponents)
                {
                    line.replace(20 + i * 14, 1, replace);
                    length = line.copy(tempchar, 13, 11 + 14 * i);
                    tempchar[length] = '\0';
                    dum_exp[exnum] = stod(tempchar);
                    exnum++;
                }
                else
                {
                    line.replace(19 + i * 14, 1, replace);
                    length = line.copy(tempchar, 14, 10 + 14 * i);
                    tempchar[length] = '\0';
                    dum_exp[exnum] = stod(tempchar);
                    exnum++;
                }
            }
        }
        else if (exnum < e_nex)
        {
            for (int i = 0; i < e_nex % 5; i++)
            {
                if (!three_exponents)
                {
                    line.replace(20 + i * 14, 1, replace);
                    length = line.copy(tempchar, 13, 11 + 14 * i);
                    tempchar[length] = '\0';
                    dum_exp[exnum] = stod(tempchar);
                    exnum++;
                }
                else
                {
                    line.replace(19 + i * 14, 1, replace);
                    length = line.copy(tempchar, 14, 10 + 14 * i);
                    tempchar[length] = '\0';
                    dum_exp[exnum] = stod(tempchar);
                    exnum++;
                }
            }
        }
        else
        {
            getline(rf, line);
            continue;
        }
        getline(rf, line);
        if (exnum > e_nex)
        {
            file << "exnum went higher than expected values in exponent reading, thats suspicius, lets stop here...\n";
            return false;
        }
        // run++;
    }
    if (exnum < e_nex)
    {
        file << "We have a problem adding exponents!\n";
        return false;
    }
    if (debug)
    {
        file << "finished with exponents, reading MOs now...\n";
        file << "line: " << line << endl;
    }
    for (int i = 0; i < e_nuc; i++)
        err_checkf(push_back_atom(dum_label[i], dum_x[i], dum_y[i], dum_z[i], dum_ch[i]), "Error while making atoms!!\n", file);
    for (int j = 0; j < e_nex; j++)
    {
        err_checkf(add_exp(dum_center[j], dum_type[j], dum_exp[j]), "Error while writing MO coefficients...\n", file);
    }
    isBohr = true;
    int linecount = 0;
    int monum = 0;
    vec2 temp_val;
    temp_val.resize(e_nmo);
    for (int i = 0; i < e_nmo; i++)
        temp_val[i].resize(e_nex);
    //-------------------------------- Read MOs --------------------------------------
    // bool orca_switch = false;
    // int temp_orca = check_order(debug),
    int temp_nr = 0;
    int oper = 0;
    double temp_occ = -1.0, temp_ener = 0.0, last_ener = -DBL_MAX;
    // if (temp_orca % 10 == 3)
    //   orca_switch = true;
    while (!(line.compare(0, 3, "END") == 0) && !rf.eof())
    {
        if (monum == e_nmo)
        {            
            file << "monum went higher than expected values in MO reading, thats suspicius, lets stop here...\n";
            break;
        }
        stringstream stream2(line);
        string tmp;
        temp_nr = 0;
        temp_occ = -1.0;
        temp_ener = 0.0;
        if (temp_nr == 0)
        {
            length = line.copy(tempchar, 6, 2);
            tempchar[length] = '\0';
            temp_nr = stoi(tempchar);
        }
        if (temp_occ == -1.0)
        {
            length = line.copy(tempchar, 12, 36);
            tempchar[length] = '\0';
            temp_occ = stod(tempchar);
        }
        if (temp_ener == 0)
        {
            length = line.copy(tempchar, 12, 61);
            tempchar[length] = '\0';
            // if we have a "=" in the line, we have to make it a space
            if (tempchar[0] == '=') tempchar[0] = ' ';
            temp_ener = stod(tempchar);
        }
        if (temp_ener > last_ener)
        {
            last_ener = temp_ener;
        }
        else
        {
            last_ener = -DBL_MAX;
            oper++;
        }
        push_back_MO(temp_nr, temp_occ, temp_ener, oper);
        //---------------------------Start reading MO coefficients-----------------------
        getline(rf, line);
        linecount = 0;
        exnum = 0;
        while (!(line.compare(0, 2, "MO") == 0) && !rf.eof())
        {
            if (exnum + 5 <= e_nex)
            {
                for (int i = 0; i < 5; i++)
                {
                    if (!three_exponents)
                        line.replace(12 + i * 16, 1, replace);
                    else
                        line.replace(11 + i * 16, 1, replace);
                    length = line.copy(tempchar, 16, 16 * i);
                    tempchar[length] = '\0';
                    temp_val[monum][exnum] = stod(tempchar);
                    exnum++;
                }
            }
            else if (exnum < e_nex)
            {
                for (int i = 0; i < (e_nex % 5); i++)
                {
                    if (!three_exponents)
                    {
                        line.replace(12 + i * 16, 1, replace);
                        length = line.copy(tempchar, 15, 1 + 16 * i);
                        tempchar[length] = '\0';
                        temp_val[monum][exnum] = stod(tempchar);
                    }
                    else
                    {
                        line.replace(11 + i * 16, 1, replace);
                        length = line.copy(tempchar, 16, 16 * i);
                        tempchar[length] = '\0';
                        temp_val[monum][exnum] = stod(tempchar);
                    }
                    exnum++;
                }
            }
            else
            {
                getline(rf, line);
                continue;
            }
            getline(rf, line);
            if (linecount * 5 > e_nex + 1)
            {
                file << "linecount went higher than expected values in exponent reading, thats suspicius, lets stop here...\n";
                return false;
            }
            // run++;
        }
        monum++;
    }
    err_checkf(monum + 1 >= e_nmo, "less MOs than expected, quitting...\nmonum: " + to_string(monum) + " e_nmo : " + to_string(e_nmo), file);
    for (int i = 0; i < e_nmo; i++)
    {
        for (int j = 0; j < e_nex; j++)
        {
            MOs[i].push_back_coef(temp_val[i][j]);
        }
    }
    constants::exp_cutoff = std::log(constants::density_accuracy / get_maximum_MO_coefficient());
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

    getline(rf, line);
    stringstream stream(line);
    int e_nuc = 0; // number of expected MOs, Exponents and nuclei
    stream >> e_nuc;
    if (debug)
        file << "e_nuc: " << e_nuc << endl;
    getline(rf, line);
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
        svec temp;
        getline(rf, line);
        stream.str(line);
        if (debug)
            file << i << ".run, line:" << line << endl;
        dum_nr[i] = i;
        temp = split_string<string>(line, " ");
        remove_empty_elements(temp);
        dum_label[i] = temp[0];
        dum_x[i] = constants::ang2bohr(stod(temp[1]));
        dum_y[i] = constants::ang2bohr(stod(temp[2]));
        dum_z[i] = constants::ang2bohr(stod(temp[3]));
        dum_ch[i] = constants::get_Z_from_label(dum_label[i].c_str()) + 1;
        if (debug)
        {
            file << "label:" << dum_label[i]
                 << " nr: " << dum_nr[i]
                 << " x: " << dum_x[i]
                 << " y: " << dum_y[i]
                 << " z: " << dum_z[i]
                 << " charge: " << dum_ch[i] << endl;
        }
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
    using namespace std;
    err_checkf(std::filesystem::exists(fileName), "Couldn't open or find " + fileName.string() + ", leaving", file);
    ifstream rf(fileName.c_str());
    path = fileName;
    string line;
    rf.seekg(0);
    getline(rf, line);
    while (line.find("<Title>") == string::npos)
        getline(rf, line);
    getline(rf, line);
    if (debug)
        file << "comment line " << line << endl;
    comment = line;
    rf.seekg(0);
    while (line.find("<Number of Nuclei>") == string::npos)
        getline(rf, line);
    getline(rf, line);
    if (debug)
        file << line << endl;
    int temp_ncen = stoi(line);
    rf.seekg(0);
    while (line.find("<Number of Primitives>") == string::npos)
        getline(rf, line);
    getline(rf, line);
    if (debug)
        file << "nex line: " << line << endl;
    int temp_nex = stoi(line);
    rf.seekg(0);
    while (line.find("<Number of Occupied Molecular Orbitals>") == string::npos)
        getline(rf, line);
    getline(rf, line);
    if (debug)
        file << "nmo line: " << line << endl;
    int temp_nmo = stoi(line);
    rf.seekg(0);
    while (line.find("<Atomic Numbers>") == string::npos)
        getline(rf, line);
    ivec nrs;
    while (true)
    {
        getline(rf, line);
        // if (debug) file << "atom number line: " << line << endl;
        if (line.find("</Atomic Numbers>") != string::npos)
            break;
        nrs.push_back(stoi(line));
    }
    err_checkf(nrs.size() == temp_ncen, "Mismatch in atom number size", file);
    rf.seekg(0);
    while (line.find("<Nuclear Cartesian Coordinates>") == string::npos)
        getline(rf, line);
    vec2 pos;
    pos.resize(3);
    double temp[3]{0, 0, 0};
    while (true)
    {
        getline(rf, line);
        if (line.find("</Nuclear Cartesian Coordinates>") != string::npos)
            break;
        istringstream is(line);
        is >> temp[0] >> temp[1] >> temp[2];
        for (int i = 0; i < 3; i++)
            pos[i].push_back(temp[i]);
    }
    err_checkf(pos[0].size() == temp_ncen, "Mismatch in atom position size", file);
    for (int i = 0; i < temp_ncen; i++)
        push_back_atom(constants::atnr2letter(nrs[i]) + to_string(i + 1), pos[0][i], pos[1][i], pos[2][i], nrs[i]);
    err_checkf(ncen == temp_ncen, "Mismatch in atom position size", file);
    isBohr = true;
    for (int i = 0; i < 3; i++)
        pos[i].resize(0);
    pos.resize(0);
    rf.seekg(0);
    while (line.find("<Net Charge>") == string::npos)
        getline(rf, line);
    getline(rf, line);
    charge = stoi(line);
    rf.seekg(0);
    while (line.find("<Electronic Spin Multiplicity>") == string::npos)
        getline(rf, line);
    getline(rf, line);
    multi = stoi(line);
    rf.seekg(0);
    while (line.find("<Primitive Centers>") == string::npos)
        getline(rf, line);
    while (true)
    {
        getline(rf, line);
        if (line.find("</Primitive Centers>") != string::npos)
            break;
        int number = CountWords(line.c_str());
        istringstream is(line);
        int _temp;
        for (int i = 0; i < number; i++)
        {
            is >> _temp;
            centers.push_back(_temp);
        }
    }
    rf.seekg(0);
    while (line.find("<Primitive Types>") == string::npos)
        getline(rf, line);
    while (true)
    {
        getline(rf, line);
        if (line.find("</Primitive Types>") != string::npos)
            break;
        int number = CountWords(line.c_str());
        istringstream is(line);
        int _temp;
        for (int i = 0; i < number; i++)
        {
            is >> _temp;
            types.push_back(_temp);
        }
    }
    rf.seekg(0);
    while (line.find("<Primitive Exponents>") == string::npos)
        getline(rf, line);
    while (true)
    {
        getline(rf, line);
        bool please_break = false;
        if (line.find("</Primitive Exponents>") != string::npos)
        {
            if (CountWords(line.c_str()) != 2)
                please_break = true;
            else
                break;
        }
        int number = CountWords(line.c_str());
        if (please_break)
            number -= 2;
        istringstream is(line);
        double _temp;
        for (int i = 0; i < number; i++)
        {
            is >> _temp;
            exponents.push_back(_temp);
        }
        if (please_break)
            break;
    }
    err_checkf(exponents.size() == temp_nex && centers.size() == temp_nex && types.size() == temp_nex, "Mismatch in numbers! aborting!", file);
    nex = temp_nex;
    rf.seekg(0);
    while (line.find("<Molecular Orbital Occupation Numbers>") == string::npos)
        getline(rf, line);
    vec occ;
    while (true)
    {
        getline(rf, line);
        bool please_break = false;
        if (line.find("</Molecular Orbital Occupation Numbers>") != string::npos)
        {
            if (CountWords(line.c_str()) != 4)
                please_break = true;
            else
                break;
        }
        int number = CountWords(line.c_str());
        istringstream is(line);
        double _temp;
        for (int i = 0; i < number; i++)
        {
            is >> _temp;
            occ.push_back(_temp);
        }
        if (please_break)
            break;
    }
    rf.seekg(0);
    while (line.find("<Molecular Orbital Energies>") == string::npos)
        getline(rf, line);
    vec ener;
    while (true)
    {
        getline(rf, line);
        bool please_break = false;
        if (line.find("</Molecular Orbital Energies>") != string::npos)
        {
            if (CountWords(line.c_str()) != 3)
                please_break = true;
            else
                break;
        }
        int number = CountWords(line.c_str());
        istringstream is(line);
        double _temp;
        for (int i = 0; i < number; i++)
        {
            is >> _temp;
            ener.push_back(_temp);
        }
        if (please_break)
            break;
    }
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
        }
        err_checkf(push_back_MO(i + 1, occ[i], ener[i], oper), "Error poshing back MO! MO: " + to_string(i), file);
    }
    occ.resize(0);
    ener.resize(0);
    rf.seekg(0);
    while (line.find("<Molecular Orbital Primitive Coefficients>") == string::npos)
        getline(rf, line);
    vec coef;
    while (line.find("</Molecular Orbital Primitive Coefficients>") == string::npos)
    {
        while (line.find("<MO Number>") == string::npos)
            getline(rf, line);
        getline(rf, line);
        // if (debug) file << "mo Nr line: " << line << endl;
        int nr = stoi(line);
        nr--;
        while (line.find("</MO Number>") == string::npos)
            getline(rf, line);
        while (coef.size() != nex)
        {
            getline(rf, line);
            // if (nr == 1 && debug) file << "first MO Coef lines: " << line << endl;
            int number = CountWords(line.c_str());
            istringstream is(line);
            double _temp;
            for (int i = 0; i < number; i++)
            {
                is >> _temp;
                coef.push_back(_temp);
            }
            err_checkf(coef.size() <= nex, "Error reading coefficients! MO: " + to_string(MOs.size()), file);
        }
        for (int i = 0; i < nex; i++)
            MOs[nr].push_back_coef(coef[i]);
        coef.resize(0);
        getline(rf, line);
    }

    //Trying to actually read in all the information where it belongs
  //  int n_occ = 0;
  //  vec2 MOs_mat;
  //  vec2 reordered_MOs_mat(MOs[0].get_primitive_count(), vec(MOs.size()));
  //  for (int i = 0; i < MOs.size(); i++) {
  //      if (MOs[i].get_occ() <= 0.0)continue;
  //      MOs_mat.push_back(MOs[i].get_coefficients());
  //      n_occ++;
  //  }
  //  centers;

  //  MOs_mat = transpose(MOs_mat);

  //  for (int type_idx = 0; type_idx < MOs_mat.size(); type_idx++) {
  //      int type = types[type_idx];
  //      if (type == 1) {//s-Type
  //          atoms[centers[type_idx] - 1].push_back_basis_set(coef[type_idx], exponents[type_idx], 1, atoms[centers[type_idx] - 1].get_shellcount_size());
  //          reordered_MOs_mat[type_idx] = MOs_mat[type_idx];
  //          continue;
  //      }else if (type >1 && type < 5)
        //{
        //    type = 1;
        //}
        //else if (type >= 5 && type < 10){
        //    type = 2;
        //}
        //else if (type >= 10 && type < 17) {
        //    type = 3;
        //}
        //else if (type >= 17 && type < 26) {
        //    type = 4;
        //}
        //else if (type >= 26 && type < 37) {
        //    type = 5;
        //}
  //      else{
        //    file << "Higher angular momentum basis functions than G, not supported!" << endl;
  //          exit(1);
        //}

  //      //Pretending to not know about contraction....
  //      int shell = atoms[centers[type_idx] - 1].get_shellcount_size();
  //      for (int m = -type; m <= type; m++) {
  //          atoms[centers[type_idx] - 1].push_back_basis_set(coef[type_idx + m + type], exponents[type_idx + m + type], type + 1, shell);
  //          reordered_MOs_mat[type_idx + constants::orca_2_pySCF[type][m]] = MOs_mat[type_idx + m + type];
  //      }
  //      type_idx += 2 * type;
  //  }


  //  vec coeff_mo(n_occ * MOs_mat.size(), 0.0);
  //  vec coeff_small(n_occ * MOs_mat.size(), 0.0);
  //  for (int i = 0; i < MOs_mat.size(); i++) {
  //      for (int oc = 0; oc < MOs.size(); oc++) {
        //    if (MOs[oc].get_occ() <= 0.0)continue;
  //          coeff_mo[i * n_occ + oc] = MOs_mat[i][oc] * MOs[oc].get_occ();
  //          coeff_small[i * n_occ + oc] = MOs_mat[i][oc];
  //      }
  //  }

  //  DM = dot(coeff_mo, coeff_small, (int)MOs_mat.size(), (int)n_occ, (int)MOs_mat.size(), (int)n_occ, false, true);

    

    while (line.find("<Energy =") == string::npos)
        getline(rf, line);
    getline(rf, line);
    total_energy = stod(line);
    while (line.find("<Virial Ratio") == string::npos)
        getline(rf, line);
    getline(rf, line);
    virial_ratio = stod(line);
    rf.close();
    constants::exp_cutoff = std::log(constants::density_accuracy / get_maximum_MO_coefficient());
    return true;
};

const double WFN::get_maximum_MO_coefficient(bool occu) const {
    double max_coef = 0.0;
    for (int i = 0; i < nmo; i++) {
        if (occu && MOs[i].get_occ() == 0.0) 
            continue;
        for(int j=0; j<nex; j++){
            if (std::abs(MOs[i].get_coefficients()[j]) > max_coef) {
                max_coef  = std::abs(MOs[i].get_coefficients()[j]);
            }
        }
    }
    return max_coef;
};

bool WFN::read_molden(const std::filesystem::path &filename, std::ostream &file, const bool debug)
{
    using namespace std;
    err_checkf(std::filesystem::exists(filename), "couldn't open or find " + filename.string() + ", leaving", file);
    if (debug)
        file << "File is valid, continuing...\n"
             << GetCurrentDir << endl;
    origin = e_origin::molden;
    ifstream rf(filename.c_str());
    if (rf.good())
        path = filename;
    string line;
    rf.seekg(0);
    // d_f_switch = true;

    getline(rf, line);
    err_checkf(line.find("Molden Format") != string::npos, "Does not look like proper molden format file!", file);
    getline(rf, line);
    comment = split_string<string>(line, "]")[1];
    bool au_bohr = false; // au = false, angs = true;
    while (line.find("[Atoms]") == string::npos)
    {
        getline(rf, line);
    }
    if (split_string<string>(line, "]")[1].find("angs") != string::npos)
        au_bohr = true;
    else if (split_string<string>(line, "]")[1].find("Angs") != string::npos)
        au_bohr = true;
    getline(rf, line);
    svec temp;
    while (line.find("]") == string::npos)
    {
        temp = split_string<string>(line, " ");
        remove_empty_elements(temp);
        if (au_bohr)
            err_checkf(push_back_atom(temp[0],
                                      constants::ang2bohr(stod(temp[3])),
                                      constants::ang2bohr(stod(temp[4])),
                                      constants::ang2bohr(stod(temp[5])),
                                      stoi(temp[2])),
                       "Error pushing back atom", file);
        else
            err_checkf(push_back_atom(temp[0],
                                      stod(temp[3]),
                                      stod(temp[4]),
                                      stod(temp[5]),
                                      stoi(temp[2])),
                       "Error pushing back atom", file);
        getline(rf, line);
    }
    err_checkf(line.find("[STO]") == string::npos, "ERROR: STOs are not yet suupported!", file);
    getline(rf, line);
    int atoms_with_basis = 0;
    while (atoms_with_basis < ncen && line.find("[") == string::npos)
    {
        svec line_digest = split_string<string>(line, " ");
        remove_empty_elements(line_digest);
        const int atom_based = stoi(line_digest[0]) - 1;
        getline(rf, line);
        int shell = 0;
        while (line.size() > 2)
        {
            line_digest = split_string<string>(line, " ");
            remove_empty_elements(line_digest);
            int shell_type;
            if (line_digest[0] == "s" || line_digest[0] == "S")
                shell_type = 1;
            else if (line_digest[0] == "p" || line_digest[0] == "P")
                shell_type = 2;
            else if (line_digest[0] == "d" || line_digest[0] == "D")
                shell_type = 3;
            else if (line_digest[0] == "f" || line_digest[0] == "F")
                shell_type = 4;
            else if (line_digest[0] == "g" || line_digest[0] == "G")
                shell_type = 5;
            else if (line_digest[0] == "h" || line_digest[0] == "H" || line_digest[0] == "i" || line_digest[0] == "I")
                err_not_impl_f("Higher angular momentum basis functions than G", file);
            getline(rf, line);
            const int number_of_functions = stoi(line_digest[1]);
            for (int i = 0; i < number_of_functions; i++)
            {
                line_digest = split_string<string>(line, " ");
                remove_empty_elements(line_digest);
                err_checkf(atoms[atom_based].push_back_basis_set(stod(line_digest[0]), stod(line_digest[1]), shell_type, shell), "Error pushing back basis", file);
                getline(rf, line);
            }
            shell++;
        }
        atoms_with_basis++;
        getline(rf, line);
    }
    bool d5 = false;
    bool f7 = false;
    bool g9 = false;
    while (line.find("[MO]") == string::npos)
    {
        if (line.find("[5D]") != string::npos || line.find("[5d]") != string::npos)
        {
            d5 = true;
        }
        if (line.find("[7F]") != string::npos || line.find("[7f]") != string::npos)
        {
            f7 = true;
        }
        if (line.find("[9G]") != string::npos || line.find("[9g]") != string::npos)
        {
            g9 = true;
        }
        if (line.find("[5D7F]") != string::npos || line.find("[5d7f]") != string::npos)
        {
            f7 = true;
            d5 = true;
        }
        if (line.find("[5D7F9G]") != string::npos || line.find("[5d7f9g]") != string::npos)
        {
            f7 = true;
            d5 = true;
            g9 = true;
        }
        getline(rf, line); // Read more lines until we reach MO block
    }
    vec3 coefficients(2);
    vec occ;
    if (d5 && f7 && g9)
    {
        int run = 0;
        string sym;
        bool spin; // alpha = false, beta = true
        double ene, occup;
        int expected_coefs = 0;
        vector<primitive> prims;
        ivec temp_shellsizes;
        for (int a = 0; a < ncen; a++)
        {
            int current_shell = -1;
            // int l = 0;
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
                        expected_coefs += 5;
                    }
                    else if (atoms[a].get_basis_set_type(s) == 4)
                    {
                        expected_coefs += 7;
                    }
                    else if (atoms[a].get_basis_set_type(s) == 5)
                    {
                        expected_coefs += 9;
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
        getline(rf, line);
        int MO_run = 0;
        vec2 p_pure_2_cart;
        vec2 d_pure_2_cart;
        vec2 f_pure_2_cart;
        vec2 g_pure_2_cart;
        err_checkf(generate_sph2cart_mat(p_pure_2_cart, d_pure_2_cart, f_pure_2_cart, g_pure_2_cart), "Error creating the conversion matrix", file);
        while (!rf.eof() && rf.good() && line.size() > 2 && line.find("[") == string::npos)
        {
            run++;
            temp = split_string<string>(line, " ");
            remove_empty_elements(temp);
            sym = temp[1];
            getline(rf, line);
            temp = split_string<string>(line, " ");
            remove_empty_elements(temp);
            ene = stod(temp[1]);
            getline(rf, line);
            temp = split_string<string>(line, " ");
            remove_empty_elements(temp);
            if (temp[1] == "Alpha" || temp[1] == "alpha")
                spin = false;
            else
                spin = true;
            getline(rf, line);
            temp = split_string<string>(line, " ");
            remove_empty_elements(temp);
            occup = stod(temp[1]);
            push_back_MO(run, occup, ene, spin);
            occ.push_back(occup);
            coefficients[spin].push_back(vec());
            // int run_coef = 0;
            int p_run = 0;
            vec2 p_temp(3);
            int d_run = 0;
            vec2 d_temp(5);
            int f_run = 0;
            vec2 f_temp(7);
            int g_run = 0;
            vec2 g_temp(9);
            int basis_run = 0;
            for (int i = 0; i < expected_coefs; i++)
            {
                getline(rf, line);
                temp = split_string<string>(line, " ");
                remove_empty_elements(temp);
                coefficients[spin][MO_run].push_back(stod(temp[1]));
                // err_checkf(temp_shellsizes[basis_run] == 1, "Please do not feed me contracted basis sets yet...", file);
                switch (prims[basis_run].get_type())
                {
                case 1:
                {
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        double t = stod(temp[1]) * prims[basis_run + s].get_coef();
                        if (abs(t) < 1E-10)
                            t = 0;
                        push_back_MO_coef(MO_run, t);
                        if (MO_run == 0)
                        {
                            push_back_exponent(prims[basis_run + s].get_exp());
                            push_back_center(prims[basis_run].get_center());
                            push_back_type(prims[basis_run].get_type());
                            nex++;
                        }
                    }
                    basis_run += temp_shellsizes[basis_run];
                    break;
                }
                case 2:
                {
                    if (p_run == 0)
                    {
                        for (int _i = 0; _i < 3; _i++)
                        {
                            p_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                        }
                    }
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        p_temp[p_run][s] = stod(temp[1]) * prims[basis_run + s].get_coef();
                    }
                    p_run++;
                    if (p_run == 3)
                    {
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            double temp_coef = 0;
                            for (int cart = 0; cart < 3; cart++)
                            {
                                temp_coef = p_temp[cart][s];
                                if (abs(temp_coef) < 1E-10)
                                    temp_coef = 0;
                                push_back_MO_coef(MO_run, temp_coef);
                                if (MO_run == 0)
                                {
                                    push_back_exponent(prims[basis_run + s].get_exp());
                                    push_back_center(prims[basis_run].get_center());
                                    if (cart == 0)
                                        push_back_type(prims[basis_run].get_type() + 2);
                                    else if (cart == 1)
                                        push_back_type(prims[basis_run].get_type());
                                    else if (cart == 2)
                                        push_back_type(prims[basis_run].get_type() + 1);
                                    nex++;
                                }
                            }
                        }
                        p_run = 0;
                        basis_run += temp_shellsizes[basis_run];
                    }
                    break;
                }
                case 3:
                {
                    if (d_run == 0)
                    {
                        for (int _i = 0; _i < 5; _i++)
                        {
                            d_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                        }
                    }
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        d_temp[d_run][s] = stod(temp[1]) * prims[basis_run + s].get_coef();
                    }
                    d_run++;
                    if (d_run == 5)
                    {
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            double temp_coef = 0;
                            for (int cart = 0; cart < 6; cart++)
                            {
                                temp_coef = 0;
                                for (int spher = 0; spher < 5; spher++)
                                {
                                    temp_coef += d_pure_2_cart[cart][spher] * d_temp[spher][s];
                                }
                                if (abs(temp_coef) < 1E-10)
                                    temp_coef = 0;
                                push_back_MO_coef(MO_run, temp_coef);
                                if (MO_run == 0)
                                {
                                    push_back_exponent(prims[basis_run + s].get_exp());
                                    push_back_center(prims[basis_run].get_center());
                                    push_back_type(5 + cart);
                                    nex++;
                                }
                            }
                        }
                        d_run = 0;
                        basis_run += temp_shellsizes[basis_run];
                    }
                    break;
                }
                case 4:
                {
                    if (f_run == 0)
                    {
                        for (int _i = 0; _i < 7; _i++)
                        {
                            f_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                        }
                    }
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        f_temp[f_run][s] = stod(temp[1]) * prims[basis_run + s].get_coef();
                    }
                    f_run++;
                    if (f_run == 7)
                    {
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            double temp_coef = 0;
                            for (int cart = 0; cart < 10; cart++)
                            {
                                temp_coef = 0;
                                for (int spher = 0; spher < 7; spher++)
                                {
                                    temp_coef += f_pure_2_cart[cart][spher] * f_temp[spher][s];
                                }
                                if (abs(temp_coef) < 1E-10)
                                    temp_coef = 0;
                                push_back_MO_coef(MO_run, temp_coef);
                                if (MO_run == 0)
                                {
                                    push_back_exponent(prims[basis_run + s].get_exp());
                                    push_back_center(prims[basis_run].get_center());
                                    push_back_type(11 + cart);
                                    nex++;
                                }
                            }
                        }
                        f_run = 0;
                        basis_run += temp_shellsizes[basis_run];
                    }
                    break;
                }
                case 5:
                {
                    if (g_run == 0)
                    {
                        for (int _i = 0; _i < 9; _i++)
                        {
                            g_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                        }
                    }
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        g_temp[g_run][s] = stod(temp[1]) * prims[basis_run + s].get_coef();
                    }
                    g_run++;
                    if (g_run == 9)
                    {
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            double temp_coef = 0;
                            for (int cart = 0; cart < 15; cart++)
                            {
                                temp_coef = 0;
                                for (int spher = 0; spher < 9; spher++)
                                {
                                    temp_coef += g_pure_2_cart[cart][spher] * g_temp[spher][s];
                                }
                                if (abs(temp_coef) < 1E-10)
                                    temp_coef = 0;
                                push_back_MO_coef(MO_run, temp_coef);
                                if (MO_run == 0)
                                {
                                    push_back_exponent(prims[basis_run].get_exp());
                                    push_back_center(prims[basis_run].get_center());
                                    push_back_type(21 + cart);
                                    nex++;
                                }
                            }
                        }
                        g_run = 0;
                        basis_run += temp_shellsizes[basis_run];
                    }
                    break;
                }
                }
            }
            err_checkf(p_run == 0 && d_run == 0 && f_run == 0 && g_run == 0, "There should not be any unfinished shells! Aborting reading molden file after MO " + to_string(MO_run) + "!", file);
            MO_run++;
            getline(rf, line);
        }
    }
    else if (!d5 && !f7 && !g9)
    {
        int run = 0;
        string sym;
        bool spin; // alpha = false, beta = true
        double ene, occup;
        int expected_coefs = 0;
        vector<primitive> prims;
        ivec temp_shellsizes;
        for (int a = 0; a < ncen; a++)
        {
            int current_shell = -1;
            // int l = 0;
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
        getline(rf, line);
        int MO_run = 0;
        while (!rf.eof() && rf.good() && line.size() > 2 && line.find("[") == string::npos)
        {
            run++;
            temp = split_string<string>(line, " ");
            remove_empty_elements(temp);
            sym = temp[1];
            getline(rf, line);
            temp = split_string<string>(line, " ");
            remove_empty_elements(temp);
            ene = stod(temp[1]);
            getline(rf, line);
            temp = split_string<string>(line, " ");
            remove_empty_elements(temp);
            if (temp[1] == "Alpha" || temp[1] == "alpha")
                spin = false;
            else
                spin = true;
            getline(rf, line);
            temp = split_string<string>(line, " ");
            remove_empty_elements(temp);
            occup = stod(temp[1]);
            push_back_MO(run, occup, ene, spin);
            occ.push_back(occup);
            coefficients[spin].push_back(vec());
            // int run_coef = 0;
            int p_run = 0;
            vec2 p_temp(3);
            int d_run = 0;
            vec2 d_temp(6);
            int f_run = 0;
            vec2 f_temp(10);
            int g_run = 0;
            vec2 g_temp(15);
            int basis_run = 0;
            for (int i = 0; i < expected_coefs; i++)
            {
                getline(rf, line);
                temp = split_string<string>(line, " ");
                remove_empty_elements(temp);
                coefficients[spin][MO_run].push_back(stod(temp[1]));
                switch (prims[basis_run].get_type())
                {
                case 1:
                {
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        double t = stod(temp[1]) * prims[basis_run + s].get_coef();
                        if (abs(t) < 1E-10)
                            t = 0;
                        push_back_MO_coef(MO_run, t);
                        if (MO_run == 0)
                        {
                            push_back_exponent(prims[basis_run + s].get_exp());
                            push_back_center(prims[basis_run].get_center());
                            push_back_type(prims[basis_run].get_type());
                            nex++;
                        }
                    }
                    basis_run += temp_shellsizes[basis_run];
                    break;
                }
                case 2:
                {
                    if (p_run == 0)
                    {
                        for (int _i = 0; _i < 3; _i++)
                        {
                            p_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                        }
                    }
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        p_temp[p_run][s] = stod(temp[1]) * prims[basis_run + s].get_coef();
                    }
                    p_run++;
                    if (p_run == 3)
                    {
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            double temp_coef = 0;
                            for (int cart = 0; cart < 3; cart++)
                            {
                                temp_coef = p_temp[cart][s];
                                if (abs(temp_coef) < 1E-10)
                                    temp_coef = 0;
                                push_back_MO_coef(MO_run, temp_coef);
                                if (MO_run == 0)
                                {
                                    push_back_exponent(prims[basis_run + s].get_exp());
                                    push_back_center(prims[basis_run].get_center());
                                    push_back_type(prims[basis_run].get_type() + cart);
                                    nex++;
                                }
                            }
                        }
                        p_run = 0;
                        basis_run += temp_shellsizes[basis_run];
                    }
                    break;
                }
                case 3:
                {
                    if (d_run == 0)
                    {
                        for (int _i = 0; _i < 6; _i++)
                        {
                            d_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                        }
                    }
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        d_temp[d_run][s] = stod(temp[1]) * prims[basis_run + s].get_coef();
                    }
                    d_run++;
                    if (d_run == 6)
                    {
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            for (int _i = 0; _i < 3; _i++)
                                push_back_MO_coef(MO_run, d_temp[_i][s]);
                            for (int _i = 3; _i < 6; _i++)
                                push_back_MO_coef(MO_run, d_temp[_i][s] * sqrt(3));
                            for (int cart = 0; cart < 6; cart++)
                            {
                                if (MO_run == 0)
                                {
                                    push_back_exponent(prims[basis_run + s].get_exp());
                                    push_back_center(prims[basis_run].get_center());
                                    push_back_type(5 + cart);
                                    nex++;
                                }
                            }
                        }
                        d_run = 0;
                        basis_run += temp_shellsizes[basis_run];
                    }
                    break;
                }
                case 4:
                {
                    if (f_run == 0)
                    {
                        for (int _i = 0; _i < 10; _i++)
                        {
                            f_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                        }
                    }
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        f_temp[f_run][s] = stod(temp[1]) * prims[basis_run + s].get_coef();
                    }
                    f_run++;
                    if (f_run == 10)
                    {
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            for (int cart = 0; cart < 10; cart++)
                            {
                                // THIS IS A MESS AND NEEDS REEVALUATION; THEY ARE CERTAINLY NOT CORRECT!
                                if (cart < 3 || cart == 9)
                                {
                                    if (cart == 9)
                                        push_back_MO_coef(MO_run, f_temp[cart][s] * sqrt(15));
                                    else
                                        push_back_MO_coef(MO_run, f_temp[cart][s]);
                                }
                                else if (cart == 3 || cart == 4)
                                    push_back_MO_coef(MO_run, f_temp[cart + 1][s] * sqrt(15));
                                else if (cart == 5)
                                    push_back_MO_coef(MO_run, f_temp[cart + 3][s] * sqrt(15));
                                else if (cart == 6)
                                    push_back_MO_coef(MO_run, f_temp[cart - 3][s] * sqrt(5));
                                else if (cart == 7)
                                    push_back_MO_coef(MO_run, f_temp[cart - 1][s] * sqrt(5));
                                else if (cart == 8)
                                    push_back_MO_coef(MO_run, f_temp[cart - 1][s] * sqrt(15));
                                if (MO_run == 0)
                                {
                                    push_back_exponent(prims[basis_run + s].get_exp());
                                    push_back_center(prims[basis_run].get_center());
                                    push_back_type(11 + cart);
                                    nex++;
                                }
                            }
                        }
                        f_run = 0;
                        basis_run += temp_shellsizes[basis_run];
                    }
                    break;
                }
                case 5:
                {
                    if (g_run == 0)
                    {
                        for (int _i = 0; _i < 15; _i++)
                        {
                            g_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                        }
                    }
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        g_temp[g_run][s] = stod(temp[1]) * prims[basis_run + s].get_coef();
                    }
                    g_run++;
                    if (g_run == 15)
                    {
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            for (int cart = 0; cart < 15; cart++)
                            {
                                push_back_MO_coef(MO_run, g_temp[cart][s]);
                                if (MO_run == 0)
                                {
                                    push_back_exponent(prims[basis_run].get_exp());
                                    push_back_center(prims[basis_run].get_center());
                                    push_back_type(21 + cart);
                                    nex++;
                                }
                            }
                        }
                        g_run = 0;
                        basis_run += temp_shellsizes[basis_run];
                    }
                    break;
                }
                }
            }
            err_checkf(p_run == 0 && d_run == 0 && f_run == 0 && g_run == 0, "There should not be any unfinished shells! Aborting reading molden file after MO " + to_string(MO_run) + "!", file);
            MO_run++;
            getline(rf, line);
        }
    }
    else
    {
        err_not_impl_f("PLEASE DONT MIX CARTESIAN AND SPERHICAL HARMINICS; THAT IS ANNOYING!", std::cout);
    }
    //Make the matrix symmetric
    while (coefficients[0].size() < coefficients[0][0].size()) {
        coefficients[0].push_back(vec(coefficients[0][0].size(), 0.0));
        occ.push_back(0);
    }
    while (coefficients[0][0].size() < coefficients[0].size()) {
        for (int i = 0; i < coefficients[0].size(); i++) {
            coefficients[0][i].push_back(0.0);
        }
        occ.push_back(0);
    }
    vec _coefficients = flatten<double>(coefficients);
    dMatrix2 m_coefs = reshape<dMatrix2>(_coefficients, Shape2D((int)coefficients[0].size(), (int)coefficients[0].size()));
    dMatrix2 temp_co = diag_dot(m_coefs, occ, true);
    DM = dot(temp_co, m_coefs);
    constants::exp_cutoff = std::log(constants::density_accuracy / get_maximum_MO_coefficient());
    return true;
};

bool WFN::read_tonto(const std::filesystem::path& filename, std::ostream& file, const bool debug, const std::filesystem::path& energies_filename, const std::filesystem::path& orbitals_filename)
{
    using namespace std;
    err_checkf(std::filesystem::exists(filename), "couldn't open or find " + filename.string() + ", leaving", file);
    std::filesystem::path energies_file, orbitals_file, stdout_file;
    ifstream rf;
    string line;
    if (energies_filename == "" && orbitals_filename == "") {
        if (filename.string().find("stdout") != std::string::npos) {
            string jobname;
            stdout_file = filename;
            rf.open(stdout_file.string().c_str(), ios::in);
            rf.seekg(0);
            while (rf.good() && line.find("Name ...") == string::npos) {
                getline(rf, line);
            }
            jobname = split_string<string>(line, " ")[2];
            energies_file = filename.parent_path() / (jobname + ".orbital_energies,restricted");
            if (!std::filesystem::exists(energies_file))
                energies_file = filename.parent_path() / (jobname + ".MO_energies,r");
            orbitals_file = filename.parent_path() / (jobname + ".molecular_orbitals,restricted");
            if (!std::filesystem::exists(orbitals_file))
                orbitals_file = filename.parent_path() / (jobname + ".MOs,r");
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
        stdout_file = filename;
        err_checkf(std::filesystem::exists(orbitals_file), "couldn't open or find " + orbitals_file.string() + ", leaving", file);
        err_checkf(std::filesystem::exists(energies_file), "couldn't open or find " + energies_file.string() + ", leaving", file);
        err_checkf(std::filesystem::exists(stdout_file), "couldn't open or find " + stdout_file.string() + ", leaving", file);
        rf.open(stdout_file.string().c_str(), ios::in);
	}


    if (debug)
        file << "File is valid, continuing...\n" << GetCurrentDir << endl;
    origin = e_origin::tonto;
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

    rf_e.close();
    rf_o.close();

    //Now read the stdout file to get the atomic positions and basis set
    rf.seekg(0);
    err_checkf(rf.good(), "couldn't open " + stdout_file.string() + ", leaving", file);
    //fast forward to the atom section
    while (rf.good() && line.find("Molecule information") == string::npos) {
        getline(rf, line);
    }
    err_checkf(rf.good(), "Couldn't find molecule information in " + stdout_file.string(), file);
    //skip 8 lines to get to the charge
    for (int i = 0; i < 8; i++) {
        getline(rf, line);
    }
    svec line_digest = split_string<string>(line, " ");
    //get the charge form the last element in the line
    charge = stoi(line_digest[line_digest.size() - 1]);
    getline(rf, line);
    line_digest = split_string<string>(line, " ");
    multi = stoi(line_digest[line_digest.size() - 1]);
    getline(rf, line);
    getline(rf, line);
    line_digest = split_string<string>(line, " ");
    const int expected_atoms = stoi(line_digest[line_digest.size() - 1]);
    getline(rf, line);
    line_digest = split_string<string>(line, " ");
    const int expected_electrons = stoi(line_digest[line_digest.size() - 1]);
    while (rf.good() && line.find("Atom coordinates") == string::npos) {
        getline(rf, line);
    }
    //skip 11 lines to get to the atom list
    for (int i = 0; i < 11; i++) {
        getline(rf, line);
    }
    line_digest = split_string<string>(line, " ");
    remove_empty_elements(line_digest);
    err_checkf(line_digest[0] == "1", "Atom list does nor start with one? let me stop...", std::cout);

    while (!line.empty() && line.front() != '_')
    {
        line_digest = split_string<string>(line, " ");
        remove_empty_elements(line_digest);
        string label = line_digest[1];
        int atomic_number = static_cast<int>(stod(line_digest[2]));
        double x = constants::ang2bohr(stod(line_digest[3]));
        double y = constants::ang2bohr(stod(line_digest[4]));
        double z = constants::ang2bohr(stod(line_digest[5]));
        err_checkf(push_back_atom(label, x, y, z, atomic_number), "Error pushing back an atom!", std::cout);
        getline(rf, line);
    }
    err_checkf(ncen == expected_atoms, "Did not read expected numebr of atoms!", std::cout);
    while (rf.good() && line.find("Gaussian basis sets") == string::npos) {
        getline(rf, line);
    }
    //get 3 lines down to the basis set name
    for (int i = 0; i < 3; i++) {
        getline(rf, line);
    }
    //read basis set
    line_digest = split_string<string>(line, " ");
    remove_empty_elements(line_digest);
    basis_set_name = line_digest[line_digest.size() - 1];
    getline(rf, line);//emtpy line
    getline(rf, line);//number of basis sets that will be printed below
    line_digest = split_string<string>(line, " ");
    remove_empty_elements(line_digest);
    const int no_basis_sets = stoi(line_digest[line_digest.size() - 1]);
    getline(rf, line);//number of shells

    line_digest = split_string<string>(line, " ");
    remove_empty_elements(line_digest);
    const int no_shells = stoi(line_digest[line_digest.size() - 1]);

    getline(rf, line);//number of shell pair (we do not consider this)
    getline(rf, line);//No of basis functions

    line_digest = split_string<string>(line, " ");
    remove_empty_elements(line_digest);
    const int no_bf = stoi(line_digest[line_digest.size() - 1]);

    getline(rf, line);//No of primitives

    line_digest = split_string<string>(line, " ");
    remove_empty_elements(line_digest);
    const int no_prim = stoi(line_digest[line_digest.size() - 1]);
    std::vector<std::pair<std::string, atom>> basis_set_data;
    std::map<char, int> l_map = { {'S',0}, {'P',1}, {'D',2}, {'F',3}, {'G',4}, {'H',5}, {'I',6}, {'s',0}, {'p',1}, {'d',2}, {'f',3}, {'g',4}, {'h',5}, {'i',6} };
    for (int nbs = 0; nbs < no_basis_sets; nbs++)
    {
        //two empty lines
        getline(rf, line);
        getline(rf, line);
        getline(rf, line);//looks like "Basis set H:3-21G"
        line_digest = split_string<string>(line, " ");
        const string atom_type = split_string<string>(line_digest[2], ":")[0];
        getline(rf, line); // empty line
        getline(rf, line);//looks like "No. of shells .... N"
        line_digest = split_string<string>(line, " ");
        const int shells_local = stoi(line_digest[4]);
        getline(rf, line);//looks like "No. of basis functions .... N"
        line_digest = split_string<string>(line, " ");
        const int bfs_local = stoi(line_digest[5]);
        getline(rf, line);//looks like "No. of primitives .... N"
        line_digest = split_string<string>(line, " ");
        const int prims_local = stoi(line_digest[4]);
        /*
__________________________________

 -L-   Fn    Exponent  Contraction
        #         /au          /au
__________________________________

*/
        for (int i = 0; i < 7; i++) getline(rf, line); //skip 6 lines to get to the shells
        atom temp_at(atom_type, "0000000000000", 0, 0, 0, 0, constants::get_Z_from_label(atom_type.c_str()));
        for (int s = 0; s < shells_local; s++)
        {
            //getline(rf, line); //get shell line
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
                getline(rf, line);
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
        for (const auto& [atom_type, atom_template] : basis_set_data)
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
    dMatrix2 coefficients = reshape<dMatrix2>(orbitals, Shape2D(expected_coefs, expected_coefs));
    vec occ(expected_coefs, 0.0);
    int alpha_els = 0, beta_els = 0, temp_els = get_nr_electrons();
    while (temp_els > 1)
    {
        alpha_els++;
        beta_els++;
        temp_els -= 2;
        if (debug)
            file << temp_els << std::endl;
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
    if (mult != 1) {
        for (int i = 0; i < alpha_els; i++)
            occ[i] = 1.0;
        for (int i = alpha_els; i < alpha_els + beta_els; i++)
            occ[i] = 1.0;
    }
    else {
        for (int i = 0; i < alpha_els; i++)
            occ[i] = 2.0;
    }
	for (int MO_run = 0; MO_run < expected_coefs; MO_run++)
    {
        push_back_MO(MO_run, occ[MO_run], energies[MO_run], 0);
        int p_run = 0;
        vec2 p_temp(3);
        int d_run = 0;
        vec2 d_temp(6);
        int f_run = 0;
        vec2 f_temp(10);
        int g_run = 0;
        vec2 g_temp(15);
        int basis_run = 0;
        for (int i = 0; i < expected_coefs; i++)
        {
            switch (prims[basis_run].get_type())
            {
            case 1:
            {
                for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                {
                    double t = coefficients(MO_run,i) * prims[basis_run + s].get_coef();
                    if (abs(t) < 1E-10)
                        t = 0;
                    push_back_MO_coef(MO_run, t);
                    if (MO_run == 0)
                    {
                        push_back_exponent(prims[basis_run + s].get_exp());
                        push_back_center(prims[basis_run].get_center());
                        push_back_type(prims[basis_run].get_type());
                        nex++;
                    }
                }
                basis_run += temp_shellsizes[basis_run];
                break;
            }
            case 2:
            {
                if (p_run == 0)
                {
                    for (int _i = 0; _i < 3; _i++)
                    {
                        p_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                    }
                }
                for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                {
                    p_temp[p_run][s] = coefficients(MO_run, i) * prims[basis_run + s].get_coef();
                }
                p_run++;
                if (p_run == 3)
                {
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        double temp_coef = 0;
                        for (int cart = 0; cart < 3; cart++)
                        {
                            temp_coef = p_temp[cart][s];
                            if (abs(temp_coef) < 1E-10)
                                temp_coef = 0;
                            push_back_MO_coef(MO_run, temp_coef);
                            if (MO_run == 0)
                            {
                                push_back_exponent(prims[basis_run + s].get_exp());
                                push_back_center(prims[basis_run].get_center());
                                push_back_type(prims[basis_run].get_type() + cart);
                                nex++;
                            }
                        }
                    }
                    p_run = 0;
                    basis_run += temp_shellsizes[basis_run];
                }
                break;
            }
            case 3:
            {
                if (d_run == 0)
                {
                    for (int _i = 0; _i < 6; _i++)
                    {
                        d_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                    }
                }
                for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                {
                    d_temp[d_run][s] = coefficients(MO_run, i) * prims[basis_run + s].get_coef() / sqrt(1.5);
                }
                d_run++;
                if (d_run == 6)
                {
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        for (int _i = 0; _i < 6; _i++)
                            push_back_MO_coef(MO_run, d_temp[_i][s]);
                        for (int cart = 0; cart < 6; cart++)
                        {
                            if (MO_run == 0)
                            {
                                push_back_exponent(prims[basis_run + s].get_exp());
                                push_back_center(prims[basis_run].get_center());
                                push_back_type(5 + cart);
                                nex++;
                            }
                        }
                    }
                    d_run = 0;
                    basis_run += temp_shellsizes[basis_run];
                }
                break;
            }
            case 4:
            {
                if (f_run == 0)
                {
                    for (int _i = 0; _i < 10; _i++)
                    {
                        f_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                    }
                }
                for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                {
                    f_temp[f_run][s] = coefficients(MO_run, i) * prims[basis_run + s].get_coef() / sqrt(5.0);
                }
                f_run++;
                if (f_run == 10)
                {
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        for (int cart = 0; cart < 10; cart++)
                        {
                            // tonto swaps type 17 and 16
                            if (cart != 5 && cart != 6)
                                push_back_MO_coef(MO_run, f_temp[cart][s]);
                            else if (cart == 5)
                                push_back_MO_coef(MO_run, f_temp[cart + 1][s]);
							else if (cart == 6)
								push_back_MO_coef(MO_run, f_temp[cart - 1][s]);
                            if (MO_run == 0)
                            {
                                push_back_exponent(prims[basis_run + s].get_exp());
                                push_back_center(prims[basis_run].get_center());
                                push_back_type(11 + cart);
                                nex++;
                            }
                        }
                    }
                    f_run = 0;
                    basis_run += temp_shellsizes[basis_run];
                }
                break;
            }
            case 5:
            {
                if (g_run == 0)
                {
                    for (int _i = 0; _i < 15; _i++)
                    {
                        g_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                    }
                }
                for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                {
                    g_temp[g_run][s] = coefficients(MO_run, i) * prims[basis_run + s].get_coef() / sqrt(13.125);
                }
                g_run++;
                if (g_run == 15)
                {
                    for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                    {
                        for (int cart = 0; cart < 15; cart++)
                        {
                            push_back_MO_coef(MO_run, g_temp[cart][s]);
                            if (MO_run == 0)
                            {
                                push_back_exponent(prims[basis_run].get_exp());
                                push_back_center(prims[basis_run].get_center());
                                push_back_type(21 + cart);
                                nex++;
                            }
                        }
                    }
                    g_run = 0;
                    basis_run += temp_shellsizes[basis_run];
                }
                break;
            }
            }
        }
        err_checkf(p_run == 0 && d_run == 0 && f_run == 0 && g_run == 0, "There should not be any unfinished shells! Aborting reading molden file after MO " + to_string(MO_run) + "!", file);
    }

    dMatrix2 temp_co = diag_dot(coefficients, occ, true);
    DM = dot(temp_co, coefficients);
    constants::exp_cutoff = std::log(constants::density_accuracy / get_maximum_MO_coefficient());
    return true;
};

template<typename T>
void move_columns(T& matrix, int dimension, int from_col, int num_cols, int to_col) {
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
    ifstream rf(filename.c_str(), ios::binary);
    if (rf.good())
        path = filename;
    string line;
    int geo_start_bit = 8;
    int basis_start_bit = 16;
    int MO_start_bit = 24;
    int soi = constants::soi; 
    int geo_int_lim = 5;

    try
    {
        rf.seekg(0, ios::beg);
        int64_t magic = 0;
        rf.read((char *)&magic, sizeof(magic));
        if (magic == -1){
            geo_start_bit += 24;
            basis_start_bit += 24;
            MO_start_bit += 24;
            soi = 8;
            geo_int_lim = 1;
        }
        // Reading geometry
        rf.seekg(geo_start_bit, ios::beg);
        int64_t geo_start = 0;
        rf.read((char *)&geo_start, sizeof(geo_start));
        err_checkf(geo_start != 0, "Could not read geometry information location from GBW file!", file);
        if (debug)
            file << "I read the pointer of geometry successfully" << endl;
        rf.seekg(geo_start, ios::beg);
        int at = 0;
        rf.read((char *)&at, constants::soi);
        double geo_vals[6]{0, 0, 0, 0, 0, 0}; // x,y,z, ch, exp_fin_nuc, mass
        int geo_ints[5]{0, 0, 0, 0, 0};
        for (int a = 0; a < at; a++)
        {
            for (int i = 0; i < 6; i++)
            {
                rf.read((char *)&(geo_vals[i]), constants::sod);
                err_checkf(rf.good(), "Error reading geo_val", file);
            }
            for (int i = 0; i < geo_int_lim; i++)
            {
                rf.read((char *)&(geo_ints[i]), soi);
                err_checkf(rf.good(), "Error reading geo_int", file);
            }
            string temp = constants::atnr2letter(geo_ints[0]);
            err_checkf(temp != "PROBLEM", "Problem identifying atoms!", std::cout);
            err_checkf(push_back_atom(temp,
                                      geo_vals[0],
                                      geo_vals[1],
                                      geo_vals[2],
                                      geo_ints[0]),
                       "Error pushing back atom", file);
        }
        if (debug)
            file << "I read the geometry of " << at << " atoms successfully" << endl;

        rf.seekg(basis_start_bit, ios::beg);
        int64_t basis_start = 0;
        rf.read((char *)&basis_start, constants::soli);
        err_checkf(basis_start != 0, "Could not read beasis information location from GBW file!", file);
        if (debug)
            file << "I read the pointer of basis set successfully" << endl;
        rf.seekg(basis_start, ios::beg);
        int atoms2 = 0, temp = 0;
        rf.read((char *)&temp, constants::soi);
        rf.read((char *)&atoms2, constants::soi);
        // long unsigned int atoms_with_basis = 0;
        vec exp(37, 0);
        vec con(37, 0);
        for (int a = 0; a < atoms2; a++)
        {
            int atom_based = 0, nr_shells = 0;
            rf.read((char *)&atom_based, constants::soi);
            err_checkf(rf.good(), "Error reading atom_based", file);
            rf.read((char *)&nr_shells, constants::soi);
            err_checkf(rf.good(), "Error reading nr_shells", file);
            int shell = 0;
            for (int p = 0; p < nr_shells; p++)
            {
                int ang_mom = 0, coeff_ind = 0, nr_funct = 0, center = 0;
                rf.read((char *)&ang_mom, constants::soi);
                err_checkf(rf.good(), "Error reading ang_mom", file);
                if (ang_mom >= 5)
                    err_not_impl_f("Higher angular momentum basis functions than G", file);
                rf.read((char *)&coeff_ind, constants::soi);
                err_checkf(rf.good(), "Error reading ceof_ind", file);
                rf.read((char *)&nr_funct, constants::soi);
                err_checkf(rf.good(), "Error reading nr_func", file);
                rf.read((char *)&center, constants::soi);
                err_checkf(rf.good(), "Error reading center", file);
                for (int b = 0; b < 37; b++)
                {
                    rf.read((char *)&(exp[b]), constants::sod);
                    err_checkf(rf.good(), "Error reading exp", file);
                }
                for (int b = 0; b < 37; b++)
                {
                    rf.read((char *)&(con[b]), constants::sod);
                    err_checkf(rf.good(), "Error reading con", file);
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
                        expected_coefs += 5;
                    }
                    else if (atoms[a].get_basis_set_type(s) == 4)
                    {
                        expected_coefs += 7;
                    }
                    else if (atoms[a].get_basis_set_type(s) == 5)
                    {
                        expected_coefs += 9;
                    }
                    else if (atoms[a].get_basis_set_type(s) == 6)
                    {
                        expected_coefs += 11;
                    }
                    else if (atoms[a].get_basis_set_type(s) == 7)
                    {
                        expected_coefs += 13;
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
        basis_set_name = "GBW read basis set";
        // int norm_const_run = 0;
        int MO_run = 0;
        vec2 p_pure_2_cart;
        vec2 d_pure_2_cart;
        vec2 f_pure_2_cart;
        vec2 g_pure_2_cart;
        err_checkf(generate_sph2cart_mat(p_pure_2_cart, d_pure_2_cart, f_pure_2_cart, g_pure_2_cart), "Error creating the conversion matrix", file);
        if (debug)
            file << "I read the basis of " << atoms2 << " atoms successfully" << endl;

        rf.seekg(MO_start_bit, ios::beg);
        int64_t MOs_start = 0;
        rf.read((char *)&MOs_start, constants::soli);
        err_checkf(rf.good(), "Error reading MO_start", file);
        err_checkf(MOs_start != 0, "Could not read MO information location from GBW file!", file);
        if (debug)
            file << "I read the pointer of MOs successfully" << endl;
        rf.seekg(MOs_start, ios::beg);
        int operators = 0, dimension = 0;
        rf.read((char *)&operators, constants::soi);
        err_checkf(rf.good(), "Error reading operators", file);
        rf.read((char *)&dimension, soi);
        err_checkf(rf.good(), "Error reading dimnesion", file);
        size_t coef_nr = size_t(dimension) * size_t(dimension);
        vec2 coefficients(operators);
        vec2 occupations(operators);
        vec2 energies(operators);
        ivec2 irreps(operators);
        ivec2 cores(operators);
        for (int i = 0; i < operators; i++)
        {
            coefficients[i].resize(coef_nr, 0);
            occupations[i].resize(dimension, 0);
            energies[i].resize(dimension, 0);
            irreps[i].resize(dimension, 0);
            cores[i].resize(dimension, 0);
            if (debug)
                file << "operators: " << operators << " coef_nr: " << coef_nr << " dimension: " << dimension << endl;
            rf.read((char *)coefficients[i].data(), constants::sod * coef_nr);
            err_checkf(rf.good(), "Error reading coefficients", file);
            if (debug)
                file << "I read the coefficients successfully" << endl;
            rf.read((char *)occupations[i].data(), constants::sod * dimension);
            err_checkf(rf.good(), "Error reading occupations", file);
            if (debug)
                file << "I read the occupations successfully" << endl;
            rf.read((char *)energies[i].data(), constants::sod * dimension);
            err_checkf(rf.good(), "Error reading energies", file);
            if (debug)
                file << "I read the energies successfully" << endl;
            rf.read((char *)irreps[i].data(), constants::soi * dimension);
            err_checkf(rf.good(), "Error reading irreps", file);
            if (debug)
                file << "I read the irreps successfully" << endl;
            rf.read((char *)cores[i].data(), constants::soi * dimension);
            err_checkf(rf.good(), "Error reading cores", file);
            if (debug)
            {
                file << "I read the cores successfully\nI am expecting " << expected_coefs << " coefficients per MO" << endl;
            }
            for (int j = 0; j < dimension; j++)
            {
                push_back_MO(i * dimension + j + 1, occupations[i][j], energies[i][j], i);
                // int run_coef = 0;
                int p_run = 0;
                vec2 p_temp(3);
                int d_run = 0;
                vec2 d_temp(5);
                int f_run = 0;
                vec2 f_temp(7);
                int g_run = 0;
                vec2 g_temp(9);
                int basis_run = 0;
                // if (debug) {
                //   file << "Starting the " << j << ". loop... wish me luck... " << endl;
                // }
                for (int p = 0; p < expected_coefs; p++)
                {
                    switch (prims[basis_run].get_type())
                    {
                    case 1:
                    {
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            double t = coefficients[i][j + p * dimension] * prims[basis_run + s].get_coef();
                            if (abs(t) < 1E-10)
                                t = 0;
                            push_back_MO_coef(MO_run, t);
                            if (MO_run == 0)
                            {
                                push_back_exponent(prims[basis_run + s].get_exp());
                                push_back_center(prims[basis_run].get_center());
                                push_back_type(prims[basis_run].get_type());
                                nex++;
                            }
                        }
                        basis_run += temp_shellsizes[basis_run];
                        break;
                    }
                    case 2:
                    {
                        if (p_run == 0)
                        {
                            for (int _i = 0; _i < 3; _i++)
                            {
                                p_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                            }
                        }
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            p_temp[p_run][s] = coefficients[i][j + p * dimension] * prims[basis_run + s].get_coef();
                        }
                        p_run++;
                        if (p_run == 3)
                        {
                            for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                            {
                                double temp_coef = 0;
                                for (int cart = 0; cart < 3; cart++)
                                {
                                    temp_coef = p_temp[cart][s];
                                    if (abs(temp_coef) < 1E-10)
                                        temp_coef = 0;
                                    push_back_MO_coef(MO_run, temp_coef);
                                    if (MO_run == 0)
                                    {
                                        push_back_exponent(prims[basis_run + s].get_exp());
                                        push_back_center(prims[basis_run].get_center());
                                        if (cart == 0)
                                            push_back_type(prims[basis_run].get_type() + 2);
                                        else if (cart == 1)
                                            push_back_type(prims[basis_run].get_type());
                                        else if (cart == 2)
                                            push_back_type(prims[basis_run].get_type() + 1);
                                        nex++;
                                    }
                                }
                            }
                            p_run = 0;
                            basis_run += temp_shellsizes[basis_run];
                        }
                        break;
                    }
                    case 3:
                    {
                        if (d_run == 0)
                        {
                            for (int _i = 0; _i < 5; _i++)
                            {
                                d_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                            }
                        }
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            d_temp[d_run][s] = coefficients[i][j + p * dimension] * prims[basis_run + s].get_coef();
                        }
                        d_run++;
                        if (d_run == 5)
                        {
                            for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                            {
                                double temp_coef = 0;
                                for (int cart = 0; cart < 6; cart++)
                                {
                                    temp_coef = 0;
                                    for (int spher = 0; spher < 5; spher++)
                                    {
                                        temp_coef += d_pure_2_cart[cart][spher] * d_temp[spher][s];
                                    }
                                    if (abs(temp_coef) < 1E-10)
                                        temp_coef = 0;
                                    push_back_MO_coef(MO_run, temp_coef);
                                    if (MO_run == 0)
                                    {
                                        push_back_exponent(prims[basis_run + s].get_exp());
                                        push_back_center(prims[basis_run].get_center());
                                        push_back_type(5 + cart);
                                        nex++;
                                    }
                                }
                            }
                            d_run = 0;
                            basis_run += temp_shellsizes[basis_run];
                        }
                        break;
                    }
                    case 4:
                    {
                        if (f_run == 0)
                        {
                            for (int _i = 0; _i < 7; _i++)
                            {
                                f_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                            }
                        }
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            f_temp[f_run][s] = coefficients[i][j + p * dimension] * prims[basis_run + s].get_coef();
                        }
                        f_run++;
                        if (f_run == 7)
                        {
                            for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                            {
                                double temp_coef = 0;
                                for (int cart = 0; cart < 10; cart++)
                                {
                                    temp_coef = 0;
                                    for (int spher = 0; spher < 7; spher++)
                                    {
                                        temp_coef += f_pure_2_cart[cart][spher] * f_temp[spher][s];
                                    }
                                    if (abs(temp_coef) < 1E-10)
                                        temp_coef = 0;
                                    push_back_MO_coef(MO_run, temp_coef);
                                    if (MO_run == 0)
                                    {
                                        push_back_exponent(prims[basis_run + s].get_exp());
                                        push_back_center(prims[basis_run].get_center());
                                        push_back_type(11 + cart);
                                        nex++;
                                    }
                                }
                            }
                            f_run = 0;
                            basis_run += temp_shellsizes[basis_run];
                        }
                        break;
                    }
                    case 5:
                    {
                        if (g_run == 0)
                        {
                            for (int _i = 0; _i < 9; _i++)
                            {
                                g_temp[_i].resize(temp_shellsizes[basis_run], 0.0);
                            }
                        }
                        for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                        {
                            g_temp[g_run][s] = coefficients[i][j + p * dimension] * prims[basis_run + s].get_coef();
                        }
                        g_run++;
                        if (g_run == 9)
                        {
                            for (int s = 0; s < temp_shellsizes[basis_run]; s++)
                            {
                                double temp_coef = 0;
                                for (int cart = 0; cart < 15; cart++)
                                {
                                    temp_coef = 0;
                                    for (int spher = 0; spher < 9; spher++)
                                    {
                                        temp_coef += g_pure_2_cart[cart][spher] * g_temp[spher][s];
                                    }
                                    if (abs(temp_coef) < 1E-10)
                                        temp_coef = 0;
                                    push_back_MO_coef(MO_run, temp_coef);
                                    if (MO_run == 0)
                                    {
                                        push_back_exponent(prims[basis_run].get_exp());
                                        push_back_center(prims[basis_run].get_center());
                                        push_back_type(21 + cart);
                                        nex++;
                                    }
                                }
                            }
                            g_run = 0;
                            basis_run += temp_shellsizes[basis_run];
                        }
                        break;
                    }
                    default:
                    {
                        if (debug)
                            file << "This is not supposed to happen!" << endl;
                        err_not_impl_f("Types higher than g type in gbws", file);
                        break;
                    }
                    }
                }
                err_checkf(p_run == 0 && d_run == 0 && f_run == 0 && g_run == 0, "There should not be any unfinished shells! Aborting reading gbw file after MO " + to_string(MO_run) + "!\nStatus (p,d,f,g): " + to_string(p_run) + " " + to_string(d_run) + " " + to_string(f_run) + " " + to_string(g_run), file);
                MO_run++;
            }
        }


        dMatrix2 reorderd_coefs_s1(dimension,dimension), reorderd_coefs_s2;
        if (operators == 2) reorderd_coefs_s2 = dMatrix2(dimension, dimension);

        dMatrixRef2 coefs_2D_s1_span(coefficients[0].data(), dimension, dimension);
        dMatrixRef2 coefs_2D_s2_span(coefficients[1].data(), dimension, dimension);

        int index = 0;
        for (const atom& _atom : atoms) {
            std::vector<basis_set_entry> basis = _atom.get_basis_set();
            int temp_bas_idx = 0;
            for (unsigned int shell = 0; shell < _atom.get_shellcount_size(); shell++) {
                int type = basis[temp_bas_idx].get_type() - 1;
                temp_bas_idx += _atom.get_shellcount(shell);
                for (int m = -type; m <= type; m++) {
                    auto coefs_2D_s1_slice = Kokkos::submdspan(coefs_2D_s1_span, index + m + type, Kokkos::full_extent);
                    auto reord_coefs_slice = Kokkos::submdspan(reorderd_coefs_s1.to_mdspan(), index + constants::orca_2_pySCF(type,m), Kokkos::full_extent);
                    std::copy(coefs_2D_s1_slice.data_handle(), coefs_2D_s1_slice.data_handle() + dimension, reord_coefs_slice.data_handle());
                    if (operators == 2) {
                        auto coefs_2D_s2_slice = Kokkos::submdspan(coefs_2D_s2_span, index + m + type, Kokkos::full_extent);
                        reord_coefs_slice = Kokkos::submdspan(reorderd_coefs_s2.to_mdspan(), index + constants::orca_2_pySCF(type,m), Kokkos::full_extent);
                        std::copy(coefs_2D_s2_slice.data_handle(), coefs_2D_s2_slice.data_handle() + dimension, reord_coefs_slice.data_handle());
                    } 
                }
                index += 2 * type + 1;
            }
        }

        //Map to collect the end index of every type
        index = 0;
        for (atom& _atom : atoms) {
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
                    std::cout << "function of type " << type << " is out of line!" << std::endl;
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
        }



        int n_occ = 0;
        for (int i = 0; i < occupations[0].size(); i++) {if (occupations[0][i] > 0.0) n_occ++;}
        
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
            // Reading ECPs?
            rf.seekg(32, ios::beg);
            long int ECP_start = 0;
            rf.read((char *)&ECP_start, sizeof(ECP_start));
            err_checkf(rf.good(), "Error reading center in ECPs", file);
            err_checkf(ECP_start != 0, "Could not read ECP information location from GBW file!", file);
            if (debug)
                file << "I read the pointer of ECP successfully" << endl;
            rf.seekg(ECP_start, ios::beg);
            long int i1 = 0;
            int i2 = 0;
            const int soi = 4;
            const int sod = 8;
            rf.read((char *)&i1, 8);
            err_checkf(rf.good(), "Error reading center in ECPs", file);
            file << "First line: " << i1 << endl;
            for (int i = 0; i < i1; i++)
            {
                rf.read((char *)&i2, 1);
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
                rf.read((char *)&Z, soi);
                err_checkf(Z > 0, "Error reading Z in ECPs", file);
                err_checkf(rf.good(), "Error reading Z in ECPs", file);
                rf.read((char *)&temp_0, soi);
                err_checkf(temp_0 > 0, "Error reading temp_0 in ECPs", file);
                err_checkf(rf.good(), "Error reading temp_0 in ECPs", file);
                char *temp_c = new char[temp_0];
                rf.read(temp_c, temp_0);
                err_checkf(rf.good(), "Error reading temp_c in ECPs", file);
                rf.read((char *)&nr_core, soi);
                err_checkf(nr_core >= 0, "Error reading nr_core in ECPs", file);
                err_checkf(rf.good(), "Error reading nr_core in ECPs", file);
                atoms[i].set_ECP_electrons(nr_core);
                rf.read((char *)&max_contract, soi);
                err_checkf(max_contract > 0, "Error reading max_contract in ECPs", file);
                err_checkf(rf.good(), "Error reading max_contract in ECPs", file);
                rf.read((char *)&max_angular, soi);
                err_checkf(max_angular > 0, "Error reading max_angular in ECPs", file);
                err_checkf(rf.good(), "Error reading max_angular in ECPs", file);
                rf.read((char *)&center, soi);
                err_checkf(center > 0, "Error reading center in ECPs", file);
                err_checkf(rf.good(), "Error reading center in ECPs", file);
                file << "I read " << Z << " " << temp_0 << " " << nr_core << " " << max_contract << " " << max_angular << endl;
                for (int l = 0; l < max_angular; l++)
                {
                    rf.read((char *)&exps, soi);
                    err_checkf(exps > 0, "Error reading exps in ECPs", file);
                    err_checkf(rf.good(), "Error reading center in ECPs", file);
                    rf.read((char *)&type, soi);
                    err_checkf(type >= 0, "Error reading type in ECPs", file);
                    err_checkf(rf.good(), "Error reading center in ECPs", file);
                    err_checkf(type < 200, "This type will give me a headache...", file);
                    err_checkf(rf.good(), "Error reading center in ECPs", file);
                    file << "There are " << exps << " exponents of type " << type << " for angular momentum " << l << endl;
                    for (int fun = 0; fun < exps; fun++)
                    {

                        rf.read((char *)&n, sod);
                        err_checkf(n < 200, "This Exponent will give me a headache...", file);
                        err_checkf(rf.good(), "Error reading center in ECPs", file);
                        rf.read((char *)&c, sod);
                        err_checkf(c < 200, "This Coefficient will give me a headache...", file);
                        err_checkf(rf.good(), "Error reading center in ECPs", file);
                        rf.read((char *)&e, sod);
                        err_checkf(e < 200, "This Exponent will give me a headache...", file);
                        err_checkf(rf.good(), "Error reading center in ECPs", file);
                        file << fun << " " << c << " " << e << " " << n << endl;
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
    constants::exp_cutoff = std::log(constants::density_accuracy / get_maximum_MO_coefficient());
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
                std::cout << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy..." << std::endl;
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
                std::cout << "ERROR in type assignement!!" << std::endl;
            }
            if (debug)
            {
                std::cout << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << std::endl
                 << "start: " << get_shell_start(a, s) << std::flush
                 << " stop: " << get_shell_end(a, s) << std::flush << std::endl
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
                    std::cout << factor << std::endl;
                for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
                {
                    if (debug)
                    {
                        std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << std::endl
                         << "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << std::endl;
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
                    std::cout << factor << std::endl;
                for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
                {
                    if (debug)
                    {
                        std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << std::endl
                         << "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << std::endl;
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
                    std::cout << factor << std::endl;
                for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
                {
                    if (debug)
                    {
                        std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << std::endl
                         << "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << std::endl;
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
                    std::cout << factor << std::endl;
                for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
                {
                    if (debug)
                    {
                        std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << std::endl
                         << "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << std::endl;
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
                std::cout << "This shell has: " << get_shell_end(a, s) - get_shell_start(a, s) + 1 << " primitives" << std::endl;
        }
    }
    return norm_const;
}

const double WFN::get_atom_coordinate(const unsigned int &nr, const unsigned int &axis) const
{
    err_checkf(!((int)nr >= ncen || axis > 2), "This input is invalid for get_atom_coordinate!", std::cout);
    return atoms[nr].get_coordinate(axis);
};

bool WFN::write_wfn(const std::filesystem::path &fileName, const bool &debug, const bool occupied)
{
    using namespace std;
    if (debug)
    {
        if (std::filesystem::exists(fileName))
        {
           std::cout << "File already existed!";
            return false;
        }
        else
        {
            if (debug)
               std::cout << "File didn't exist before, writing comment to it now." << endl;
        }
    }

    ofstream rf(fileName, ios::out);
    string line;
    if (!rf.is_open())
    {
       std::cout << "Sorry, can't open the file...\n";
        return false;
    }
    rf << comment << endl;
    if (debug)
       std::cout << "comment written, now for the header..\n";
    rf << hdr(occupied);
    if (debug)
    {
       std::cout << "header written, now for the centers..\n";
       std::cout << "this is the header: \n"
             << hdr(occupied);
    }
    rf.flush();
    for (int i = 0; i < ncen; i++)
    {
        rf << setw(5) << atoms[i].get_label();
        rf << setw(3) << i + 1 << "    (CENTRE ";
        if (i < 9)
            rf << ' ';
        rf << i + 1 << ") ";
        rf << fixed << showpoint << setprecision(8);
        rf << setw(12) << get_atom_coordinate(i,0);
        rf << setw(12) << get_atom_coordinate(i,1);
        rf << setw(12) << get_atom_coordinate(i,2);
        rf << "  CHARGE = ";
        rf << fixed << showpoint << setprecision(1) << setw(2) << get_atom_charge(i);
        rf << ".0";
        rf << '\n';
    }
    if (debug)
       std::cout << "centers written, now for the center_assignement..\n";
    if (debug)
       std::cout << "ncen: " << ncen << " nex: " << nex << " nmo: " << nmo << endl;
    int run = 0;
    int exnum = 0;
    for (int i = 0; i < nex / 20; i++)
    {
        rf << "CENTRE ASSIGNMENTS  ";
        for (int j = 0; j < 20; j++)
        {
            rf << setw(3) << centers[exnum];
            if (exnum > nex)
            {
               std::cout << "run is too big in center writing";
                if (debug)
                   std::cout << "in 20er-lines...\n";
                return false;
            }
            exnum++;
        }
        run++;
        rf << '\n';
    }
    if (debug)
       std::cout << "this should be the last line... \n";
    if (exnum < nex)
    {
        rf << "CENTRE ASSIGNMENTS  ";
        for (int j = 0; j < nex % 20; j++)
        {
            rf << setw(3) << centers[exnum];
            if (exnum > nex)
            {
               std::cout << "run is too big in center writing";
                if (debug)
                   std::cout << " in last line... trying to access # " << exnum << "\n";
                return false;
            }
            exnum++;
        }
        rf << '\n';
    }
    if (run * 20 < nex / 20 - 1)
    {
       std::cout << "Problem during writing of Centre assignments... stopping...\n";
        return false;
    }
    if (debug)
       std::cout << "center assignements written, now for the types..\n";
    run = 0;
    exnum = 0;
    for (int i = 0; i < nex / 20; i++)
    {
        rf << "TYPE ASSIGNMENTS    ";
        for (int j = 0; j < 20; j++)
        {
            rf << setw(3) << types[exnum];
            if (exnum > nex)
            {
               std::cout << "run is too big in types writing\n";
                return false;
            }
            exnum++;
        }
        run++;
        rf << '\n';
    }
    if (exnum < nex)
    {
        rf << "TYPE ASSIGNMENTS    ";
        int final_j = 0;
        for (int j = 0; j < nex % 20; j++)
        {
            rf << setw(3) << types[exnum];
            if (exnum > nex)
            {
               std::cout << "run is too big in types writing";
                return false;
            }
            final_j = j;
            exnum++;
        }
        if (debug)
           std::cout << "final_j: " << final_j << endl;
        rf << '\n';
    }
    if (run * 20 < nex / 20 - 1)
    {
       std::cout << "Problem during writing of Type assignments... stopping...";
        return false;
    }
    if (debug)
       std::cout << "types assignements written, now for the exponents..\n";
    run = 0;
    exnum = 0;
    for (int i = 0; i < nex / 5; i++)
    {
        rf << "EXPONENTS ";
        for (int j = 0; j < 5; j++)
        {
            stringstream stream;
            string temp;
            stream << uppercase << scientific << setw(14) << setprecision(7) << exponents[exnum];
            temp = stream.str();
            rf << temp;
            if (exnum > nex)
            {
               std::cout << "run is too big in exponents writing";
                return false;
            }
            exnum++;
        }
        run++;
        rf << '\n';
    }
    if (exnum < nex)
    {
        rf << "EXPONENTS ";
        for (int j = 0; j < nex % 5; j++)
        {
            stringstream stream;
            string temp;
            stream << uppercase << scientific << setw(14) << setprecision(7) << exponents[exnum];
            temp = stream.str();
            rf << temp;
            if (run > nex)
            {
               std::cout << "run is too big in exponents writing";
                return false;
            }
            exnum++;
        }
        rf << '\n';
    }
    if (run * 5 < nex / 5 - 1)
    {
       std::cout << "Problem during writing of Exponents... stopping...";
        return false;
    }
    if (debug)
       std::cout << "exponents assignements written, now for the MOs.." << endl
             << "For informational purposes: ncen "
             << ncen << " nmo " << nmo << " nex " << nex << endl;
    int mo_run = 1;
    for (int mo_counter = 0; mo_counter < nmo; mo_counter++)
    {
        if (occupied && MOs[mo_counter].get_occ() == 0)
            continue;
        // rf << MOs[mo_counter].hdr();
        rf << "MO" << setw(3) << mo_run << setw(29) << "OCC NO =" << setw(13) << fixed << setprecision(8) << MOs[mo_counter].get_occ()
           << setw(14) << "ORB. ENERGY =" << setw(13) << fixed << setprecision(8) << MOs[mo_counter].get_energy() << endl;
        run = 0;
        for (int i = 0; i < nex / 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                stringstream stream;
                string temp;
                stream << uppercase << scientific << showpoint << setprecision(8) << setw(16) << MOs[mo_counter].get_coefficient(run);
                temp = stream.str();
                rf << temp;
                if (run > nex)
                {
                   std::cout << "run (" << run << ") is too big in MO ceofficients writing" << endl;
                    return false;
                }
                run++;
            }
            rf << '\n';
        }
        if (run < nex)
        {
            if (debug)
               std::cout << "Still some left to write... going in % for loop...." << endl;
            for (int j = 0; j < nex % 5; j++)
            {
                stringstream stream;
                string temp;
                stream << uppercase << scientific << showpoint << setprecision(8) << setw(16) << MOs[mo_counter].get_coefficient(run);
                temp = stream.str();
                rf << temp;
                if (run > nex)
                {
                   std::cout << "run (" << run << ") is too big in MO ceofficients writing" << endl;
                    return false;
                }
                run++;
            }
            rf << '\n';
        }
        mo_run++;
    }
    if (run != nex)
    {
       std::cout << "Problem during writing of MOs... stopping...";
        if (debug)
           std::cout << "run: " << run << endl;
        return false;
    }
    rf << "END DATA" << endl;
    rf.flush();
    rf.close();
    return true;
};

bool WFN::write_nbo(const std::filesystem::path& fileName, const bool& debug)
{
    //We want to write a .47 file that looks like this according to the NBO manual:
    /*
 $GENNBO NATOMS=7 NBAS=28 UPPER BODM FORMAT=PRECISE $END 
 $NBO $END
 $COORD
 Methylamine in 3-21G basis set
 6 6 0.745914 0.011106 0.000000
 7 7 -0.721743 -0.071848 0.000000
 1 1 1.042059 1.060105 0.000000
 1 1 1.129298 -0.483355 0.892539
 1 1 1.129298 -0.483355 -0.892539
 1 1 -1.076988 0.386322 -0.827032
 1 1 -1.076988 0.386322 0.827032
 $END
 $BASIS
    CENTER = 1 1 1 1 1 1 1 1 1 2 2 2 2
             2 2 2 2 2 3 3 4 4 5 5 6 6
             7 7 
    LABEL = 1 1 101 102 103 1 101 102 103 1 1 101 102 
            103 1 101 102 103 1 1 1 1 1 1 1 1 
            1 1 
 $END 
 $CONTRACT
  NSHELL = 16 
    NEXP = 27 
   NCOMP = 1 4 4 1 4 4 1 1 1 1 1 1 1 
           1 1 1 
   NPRIM = 3 2 1 3 2 1 2 1 2 1 2 1 2 
           1 2 1 
    NPTR = 1 4 6 7 10 12 13 15 16 18 19 21 22 24 25 27 
    EXP = 0.172256000000E+03 0.259109000000E+02 0.553335000000E+01 
          0.366498000000E+01 0.770545000000E+00 0.195857000000E+00 
          0.242766000000E+03 0.364851000000E+02 0.781449000000E+01 
          0.542522000000E+01 0.114915000000E+01 0.283205000000E+00 
          0.544717800000E+01 0.824547240000E+00 0.183191580000E+00 
          0.544717800000E+01 0.824547240000E+00 0.183191580000E+00 
          0.544717800000E+01 0.824547240000E+00 0.183191580000E+00 
          0.544717800000E+01 0.824547240000E+00 0.183191580000E+00 
          0.544717800000E+01 0.824547240000E+00 0.183191580000E+00 
     CS = 0.617669074000E-01 0.358794043000E+00 0.700713084000E+00 
         -0.395895162000E+00 0.121583436000E+01 0.100000000000E+01 
          0.598657005000E-01 0.352955003000E+00 0.706513006000E+00 
         -0.413300077000E+00 0.122441727000E+01 0.100000000000E+01 
          0.156284979000E+00 0.904690877000E+00 0.100000000000E+01 
          0.156284979000E+00 0.904690877000E+00 0.100000000000E+01 
          0.156284979000E+00 0.904690877000E+00 0.100000000000E+01 
          0.156284979000E+00 0.904690877000E+00 0.100000000000E+01 
          0.156284979000E+00 0.904690877000E+00 0.100000000000E+01 
     CP = 0.000000000000E+00 0.000000000000E+00 0.000000000000E+00 
          0.236459947000E+00 0.860618806000E+00 0.100000000000E+01 
          0.000000000000E+00 0.000000000000E+00 0.000000000000E+00 
          0.237972016000E+00 0.858953059000E+00 0.100000000000E+01 
          0.000000000000E+00 0.000000000000E+00 0.000000000000E+00 
          0.000000000000E+00 0.000000000000E+00 0.000000000000E+00 
          0.000000000000E+00 0.000000000000E+00 0.000000000000E+00 
          0.000000000000E+00 0.000000000000E+00 0.000000000000E+00 
          0.000000000000E+00 0.000000000000E+00 0.000000000000E+00 
 $END 
 $OVERLAP 
          0.100000000000E+01 0.191447444408E+00 0.100000000000E+01 
 $END 
 $DENSITY 
          0.203642496554E+01 0.110916720865E+00 0.103889621321E+00 
 $END 
 $LCAOMO 
         -0.581395484288E-03 -0.241638924924E-02 -0.179639931958E-02 
 $END
    */
    using namespace std;

	//make a copy of this WFN to compute the overlap matrix
    WFN copy(*this);
    copy.delete_basis_set();
    
    //Add the decontracted basis to the copy
    for (int e = 0; e < nex; e++){
        int type = 0;
        if (types[e] == 1)
            type = 1;
        else if (types[e] == 2)
            type = 2;
        else if (types[e] == 5)
            type = 3;
        else if (types[e] == 11)
            type = 4;
        else if (types[e] == 21)
            type = 5;
        else
            continue;
        unsigned int shellsize = std::max(copy.get_atom_shell_count(centers[e] - 1), 0);
		copy.push_back_atom_basis_set(centers[e] - 1, exponents[e], 1.0, type, shellsize);

	}

    vec OVLP_matrix = {};
	Int_Params int_params(copy);

	compute2c_Overlap_Cart(int_params, OVLP_matrix);
	//We have the overlap matrix, now write it to file

    ofstream rf(fileName, ios::out);
    string line;
	stringstream stream;
    if (!rf.is_open())
    {
        std::cout << "Sorry, can't open the file...\n";
        return false;
    }

    

	rf << " $GENNBO NATOMS=" << ncen << " NBAS=" << nex << " UPPER BODM FORMAT=PRECISE $END" << endl;
	rf << " $NBO NBO NRT $END" << endl;
	rf << " $COORD" << endl;
	rf << " .47 file generated by NoSpherA2 based on " << path << endl;
	for (int i = 0; i < ncen; i++)
        //To-do: Fix second charge mention for ECPs depending on the mode of the wfn / maybe by file type?!
		rf << " " << setw(2) << get_atom_charge(i) << " " << setw(2) << get_atom_charge(i) << " " << fixed << setprecision(6) << setw(8) << get_atom_coordinate(i, 0) << " " << fixed << setprecision(6) << setw(8) << get_atom_coordinate(i, 1) << " " << fixed << setprecision(6) << setw(8) << get_atom_coordinate(i, 2) << endl;
	rf << " $END" << endl;
	rf << " $BASIS" << endl;
	rf << "    CENTER = ";
    for (int i = 0; i < nex; i++) {
		rf << centers.at(i) << " ";
        if (i%13 == 0 && i!= 0)
			rf << "\n             ";
    }
    rf.flush();
	rf << "\n    LABEL = ";

    int highest_angular = -1;
    for (int i = 0; i < nex; i++) {
        rf << constants::type_2_nbo(types.at(i)) << " ";
        if (i % 13 == 0 && i != 0)
            rf << "\n             ";
    }
    rf.flush();
	rf << "\n $END" << endl;

    string comp_string = "";
	string exp_string = "";
    string nprim_string = "";
    string nptr_string = "";
    svec cx_string(5);
    int nshell = 0;
    for (int i = 0; i < nex; i++) {
        string nr = to_string(constants::type_2_nbo(types.at(i)));
        if (nr.back() == '1') {
            nshell++;
            if (nr.front() == '1') {
                comp_string += "  1 ";
                cx_string[0] += "0.100000000000E+01 ";
                cx_string[1] += "0.000000000000E+00 ";
                cx_string[2] += "0.000000000000E+00 ";
                cx_string[3] += "0.000000000000E+00 ";
                cx_string[4] += "0.000000000000E+00 ";
                highest_angular = std::max(highest_angular, 0);
            }
            else if (nr.front() == '2') {
                comp_string += "  3 ";
                cx_string[0] += "0.000000000000E+00 ";
                cx_string[1] += "0.100000000000E+01 ";
                cx_string[2] += "0.000000000000E+00 ";
                cx_string[3] += "0.000000000000E+00 ";
                cx_string[4] += "0.000000000000E+00 ";
                highest_angular = std::max(highest_angular, 1);
            }
            else if (nr.front() == '3') {
                comp_string += "  6 ";
                cx_string[0] += "0.000000000000E+00 ";
                cx_string[1] += "0.000000000000E+00 ";
                cx_string[2] += "0.100000000000E+01 ";
                cx_string[3] += "0.000000000000E+00 ";
                cx_string[4] += "0.000000000000E+00 ";
                highest_angular = std::max(highest_angular, 2);
            }
            else if (nr.front() == '4') {
                comp_string += " 10 ";
                cx_string[0] += "0.000000000000E+00 ";
                cx_string[1] += "0.000000000000E+00 ";
                cx_string[2] += "0.000000000000E+00 ";
                cx_string[3] += "0.100000000000E+01 ";
                cx_string[4] += "0.000000000000E+00 ";
                highest_angular = std::max(highest_angular, 3);
            }
            else if (nr.front() == '5') {
                comp_string += " 15 ";
                cx_string[0] += "0.000000000000E+00 ";
                cx_string[1] += "0.000000000000E+00 ";
                cx_string[2] += "0.000000000000E+00 ";
                cx_string[3] += "0.000000000000E+00 ";
                cx_string[4] += "0.100000000000E+01 ";
                highest_angular = std::max(highest_angular, 4);
            }
			stream << setw(18) << scientific << setprecision(11) << exponents.at(i);
            exp_string += stream.str();
            stream.str("");
            nprim_string += "  1 ";
            stream << setw(4) << to_string(i + 1);
            nptr_string += stream.str();
            stream.str("");
            if ((nshell % 13) == 0 && i != 0) {
                comp_string += "\n            ";
                nprim_string += "\n            ";
                nptr_string += "\n           ";
            }
            if ((nshell % 3) == 0 && i != 0) {
                exp_string += "\n          ";
                for (int j = 0; j < 5; j++)
                    cx_string[j] += "\n           ";
            }
        }
    }
	rf << " $CONTRACT" << endl;
	rf << "  NSHELL = " << nshell << endl;
	rf << "    NEXP = " << nshell << endl;
    rf << "   NCOMP =  " << comp_string << endl;
	rf << "   NPRIM =  " << nprim_string << endl;
    rf << "    NPTR = " << nptr_string << endl;
	rf << "     EXP =" << exp_string << endl;
	if (highest_angular >= 0) rf << "      CS = " << cx_string[0] << endl;
	if (highest_angular >= 1) rf << "      CP = " << cx_string[1] << endl;
	if (highest_angular >= 2) rf << "      CD = " << cx_string[2] << endl;
	if (highest_angular >= 3) rf << "      CF = " << cx_string[3] << endl;
	if (highest_angular >= 4) rf << "      CG = " << cx_string[4] << endl;
	rf << " $END" << endl;
	rf << " $OVERLAP " << endl;
	dMatrixRef2 OVLP_mat(OVLP_matrix.data(), nex, nex);
	int runner = 0;
    for (int i = 0; i < nex; i++) {
        for (int j = i; j < nex; j++) {
			runner++;
            stream.str("");
            stream << uppercase << scientific << showpoint << setprecision(12) << setw(19) << OVLP_mat(i, j);
            rf << stream.str() << " ";
            if (runner % 3 == 0)
				rf << "\n";
        }
    }
	rf << " $END" << endl;
	rf << " $DENSITY " << endl;
    runner = 0;
	//Build the density matrix in the decontracted basis
    vec coef;
    vec occ(nex, 0.0);
    coef.resize(nex*nex);
    for(int mo = 0; mo < nex; mo++){
        if (mo >= nmo) {
            for (int i = 0; i < nex; i++) {
                coef[mo * nex + i] = 0.0;
            }
            continue;
        }
        for(int i = 0; i < nex; i++){
			coef[mo*nex + i] = MOs[mo].get_coefficient_f(i);
        }
		occ[mo] = MOs[mo].get_occ();
    }
    dMatrix2 coefficients = reshape<dMatrix2>(coef, Shape2D(nex, nex));
	dMatrix2 coef_occ = diag_dot(coefficients, occ);
    dMatrix2 DM_new = dot(coef_occ, coefficients);

    for (int i = 0; i < nex; i++) {
        for (int j = i; j < nex; j++) {
            runner++;
            stream.str("");
            stream << uppercase << scientific << showpoint << setprecision(12) << setw(19) << DM_new(i, j);
            rf << stream.str() << " ";
            if (runner % 3 == 0)
                rf << "\n";
        }
    }
	rf << " $END" << endl;
	rf << " $LCAOMO " << endl;
    for (int mo_counter = 0; mo_counter < nex; mo_counter++)
    {
        if (debug)
            std::cout << "Writing MO #" << mo_counter + 1 << "...\n";
        for (int i = 0; i < nex; i++) {
            stream.str("");
            stream << uppercase << scientific << showpoint << setprecision(12) << setw(19) << coefficients(mo_counter, i);
            rf << stream.str() << " ";
            if ((i + 1) % 3 == 0)
                rf << "\n";
        }
	}
	rf << " $END" << endl;
    rf.close();
    return true;
};

bool WFN::write_xyz(const std::filesystem::path& fileName)
{
    using namespace std;
    try {
        ofstream f(fileName, ios::out);
        f << ncen << endl;
        f << "XYZ File written by NoSpherA2 based on " << path << endl;
        for (int i = 0; i < ncen; i++)
            if (atoms[i].get_label() == "")
                f << constants::atnr2letter(get_atom_charge(i)) << setw(14) << setprecision(8) << get_atom_coordinate(i,0) << setw(14) << setprecision(8) << get_atom_coordinate(i,1) << setw(14) << setprecision(8) << get_atom_coordinate(i,2) << endl;
            else{
                if (isBohr) {
                    f << atoms[i].get_label() << setw(14) << setprecision(8) << constants::bohr2ang(get_atom_coordinate(i,0)) << setw(14) << setprecision(8) << constants::bohr2ang(get_atom_coordinate(i,1)) << setw(14) << setprecision(8) << constants::bohr2ang(get_atom_coordinate(i,2)) << endl;
                }else
                {
                    f << atoms[i].get_label() << setw(14) << setprecision(8) << get_atom_coordinate(i,0) << setw(14) << setprecision(8) << get_atom_coordinate(i,1) << setw(14) << setprecision(8) << get_atom_coordinate(i,2) << endl;
                }
            }
        f.flush();
        f.close();
    }
    catch (exception) {
        err("Error writing the xyz file! Aborting!",std::cout);
        return false;
    }
    return true;
};

void WFN::print_primitive(const int &nr) const
{
    std::cout << "center assignement: " << centers[nr] << " type: " << types[nr]
         << " exponent: " << exponents[nr] << std::endl
         << "MO coefficients:";
    for (int i = 0; i < nmo; i++)
    {
        std::cout << MOs[nr].get_coefficient(i) << "   ";
        if (i % 5 == 0)
            std::cout << std::endl;
    }
};

const int WFN::get_nmo(const bool &only_occ) const
{
    if (!only_occ)
        return nmo;
    else
    {
        int count = 0;
#pragma omp parallel for reduction(+ : count)
        for (int i = 0; i < MOs.size(); i++)
        {
            if (MOs[i].get_occ() != 0.0)
            {
                count++;
            }
        }
        return count;
    }
};

const unsigned int WFN::get_nr_electrons() const
{
    unsigned int count = 0;
    for (int i = 0; i < ncen; i++)
        count += get_atom_charge(i);
    count -= charge;
    return count;
};

const unsigned int WFN::get_nr_ECP_electrons() const
{
    unsigned int count = 0;
    for (int i = 0; i < ncen; i++)
        count += atoms[i].get_ECP_electrons();
    return count;
}

double WFN::count_nr_electrons(void) const 
{
    double count = 0;
    for (int i = 0; i < nmo; i++)
        count += MOs[i].get_occ();
    return count;
};

const double WFN::get_atom_basis_set_exponent(const int &nr_atom, const int &nr_prim) const
{
    if (nr_atom <= ncen && nr_atom >= 0 && (int)atoms[nr_atom].get_basis_set_size() >= nr_prim && nr_prim >= 0)
        return atoms[nr_atom].get_basis_set_exponent(nr_prim);
    else
        return -1;
};

const double WFN::get_atom_basis_set_coefficient(const int &nr_atom, const int &nr_prim) const
{
    if (nr_atom <= ncen && nr_atom >= 0 && (int)atoms[nr_atom].get_basis_set_size() >= nr_prim && nr_prim >= 0)
        return atoms[nr_atom].get_basis_set_coefficient(nr_prim);
    else
        return -1;
};

bool WFN::change_atom_basis_set_exponent(const int &nr_atom, const int &nr_prim, const double &value)
{
    if (nr_atom <= ncen && nr_atom >= 0 && (int)atoms[nr_atom].get_basis_set_size() >= nr_prim && nr_prim >= 0)
    {
        atoms[nr_atom].set_basis_set_exponent(nr_prim, value);
        set_modified();
        return true;
    }
    else
        return false;
};

bool WFN::change_atom_basis_set_coefficient(const int &nr_atom, const int &nr_prim, const double &value)
{
    err_checkf(nr_atom <= ncen && nr_atom >= 0 && (int)atoms[nr_atom].get_basis_set_size() >= nr_prim && nr_prim >= 0, "Wrong input!", std::cout);
    atoms[nr_atom].set_basis_set_coefficient(nr_prim, value);
    set_modified();
    return true;
};

const int WFN::get_atom_primitive_count(const int &nr) const
{
    if (nr <= ncen && nr >= 0)
        return (int)atoms[nr].get_basis_set_size();
    else
        return -1;
};

const int WFN::get_basis_set_shell(const unsigned int &nr_atom, const unsigned int &nr_prim) const
{
    if ((int)nr_atom <= ncen && atoms[nr_atom].get_basis_set_size() >= (int)nr_prim)
    {
        return atoms[nr_atom].get_basis_set_shell(nr_prim);
    }
    else
        return -1;
};

const int WFN::get_atom_primitive_type(const int& nr_atom, const int& nr_prim) const
{
    if (nr_atom < atoms.size() && nr_atom >= 0 && nr_prim < (int)atoms[nr_atom].get_basis_set_size() && nr_prim >= 0)
        return atoms[nr_atom].get_basis_set_type(nr_prim);
    else
        return -1;
};

const int WFN::get_atom_shell_count(const unsigned int &nr) const
{
    if ((int)nr <= ncen)
        return (int)atoms[nr].get_shellcount_size();
    else
        return -1;
};

const int WFN::get_atom_shell_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
    if ((int)nr_atom <= ncen && (int)nr_shell < atoms[nr_atom].get_shellcount_size())
        return atoms[nr_atom].get_shellcount(nr_shell);
    else
        return -1;
};

const int WFN::get_shell_type(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
    if (static_cast<int>(nr_atom) <= ncen && nr_shell <= atoms[nr_atom].get_shellcount_size())
    {
        int primitive_counter = 0;
        while (atoms[nr_atom].get_basis_set_shell(primitive_counter) != nr_shell)
            primitive_counter++;
        return atoms[nr_atom].get_basis_set_type(primitive_counter);
    }
    else
        return -1;
};

const int WFN::get_shell_center(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
    if (static_cast<int>(nr_atom) <= ncen && nr_shell <= atoms[nr_atom].get_shellcount_size())
        return centers[get_shell_start_in_primitives(nr_atom, nr_shell)];
    else
        return -1;
};

const int WFN::get_shell_start(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
    if (static_cast<int>(nr_atom) <= ncen && nr_shell <= atoms[nr_atom].get_shellcount_size() - 1)
    {
        int primitive_counter = 0;
#pragma loop(no_vector)
        for (int s = 0; s < static_cast<int>(nr_shell); s++)
            primitive_counter += atoms[nr_atom].get_shellcount(s);
        return primitive_counter;
    }
    else
        return -1;
};

const int WFN::get_shell_start_in_primitives(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
    if (static_cast<int>(nr_atom) <= ncen && nr_shell <= atoms[nr_atom].get_shellcount_size() - 1)
    {
        int primitive_counter = 0;
        for (unsigned int a = 0; a < nr_atom; a++)
            for (unsigned int s = 0; s < atoms[a].get_shellcount_size(); s++)
                switch (get_shell_type(a, s))
                {
                case 1:
                    primitive_counter += atoms[a].get_shellcount(s);
                    break;
                case 2:
                    primitive_counter += (3 * atoms[a].get_shellcount(s));
                    break;
                case 3:
                    primitive_counter += (6 * atoms[a].get_shellcount(s));
                    break;
                case 4:
                    primitive_counter += (10 * atoms[a].get_shellcount(s));
                    break;
                }
        for (unsigned int s = 0; s < nr_shell; s++)
        {
            switch (get_shell_type(nr_atom, s))
            {
            case 1:
                primitive_counter += atoms[nr_atom].get_shellcount(s);
                break;
            case 2:
                primitive_counter += (3 * atoms[nr_atom].get_shellcount(s));
                break;
            case 3:
                primitive_counter += (6 * atoms[nr_atom].get_shellcount(s));
                break;
            case 4:
                primitive_counter += (10 * atoms[nr_atom].get_shellcount(s));
                break;
            }
        }
        return primitive_counter;
    }
    else
        return -1;
};

const int WFN::get_shell_end(const unsigned int &nr_atom, const unsigned int &nr_shell) const
{
    if (static_cast<int>(nr_atom) <= ncen && nr_atom >= 0 && nr_shell <= atoms[nr_atom].get_shellcount_size() && static_cast<int>(nr_atom) >= 0)
    {
        if (nr_shell == atoms[nr_atom].get_shellcount_size() - 1)
            return (int)atoms[nr_atom].get_basis_set_size() - 1;
        int primitive_counter = 0;
        while (atoms[nr_atom].get_basis_set_shell(primitive_counter) != (nr_shell + 1))
            primitive_counter++;
        return primitive_counter - 1;
    }
    else
        return -1;
};

const std::string WFN::get_atom_label(const unsigned int &nr) const
{
    std::string error_return{ '?' };
    if (nr < static_cast<unsigned int>(ncen))
        return atoms[nr].get_label();
    else
        return error_return;
};

const int WFN::get_nr_basis_set_loaded() const
{
    int count = 0;
    for (int a = 0; a < ncen; a++)
        if (atoms[a].get_basis_set_loaded())
            count++;
    return count;
};

const bool WFN::get_atom_basis_set_loaded(const int &nr) const
{
    if (nr <= ncen && nr >= 0)
        return atoms[nr].get_basis_set_loaded();
    else
    {
        std::cout << "invalid atom choice in atom_basis_set_loaded!" << std::endl;
        return false;
    }
};

const int WFN::get_atom_charge(const int &nr) const
{
    if (nr <= ncen && nr >= 0)
        return atoms[nr].get_charge();
    else
    {
        std::cout << "invalid atom choice in atom_basis_set_loaded!" << std::endl;
        return -1;
    }
};

void WFN::push_back_DM(const double &value)
{
    UT_DensityMatrix.push_back(value);
};

void WFN::resize_DM(const int &size, const double &value)
{
    UT_DensityMatrix.resize(size, value);
};

const double WFN::get_DM(const int &nr) const
{
    if (nr >= 0 && nr < UT_DensityMatrix.size())
        return UT_DensityMatrix[nr];
    else
    {
        std::cout << "Requested nr out of range! Size: " << UT_DensityMatrix.size() << " nr: " << nr << std::endl;
        return -1;
    }
};

bool WFN::set_DM(const int &nr, const double &value)
{
    if (nr >= 0 && nr < UT_DensityMatrix.size())
    {
        UT_DensityMatrix[nr] = value;
        return true;
    }
    else
    {
        std::cout << "invalid arguments for set_DM! Input was: " << nr << ";" << value << std::endl;
        return false;
    }
};

void WFN::push_back_SDM(const double &value)
{
    UT_SpinDensityMatrix.push_back(value);
};

void WFN::resize_SDM(const int &size, const double &value)
{
    UT_SpinDensityMatrix.resize(size, value);
};

const double WFN::get_SDM(const int &nr) const
{
    if (nr >= 0 && nr < UT_SpinDensityMatrix.size())
        return UT_SpinDensityMatrix[nr];
    else
    {
        std::cout << "Requested nr out of range! Size: " << UT_SpinDensityMatrix.size() << " nr: " << nr << std::endl;
        return -1;
    }
};

bool WFN::set_SDM(const int &nr, const double &value)
{
    if (nr >= 0 && nr < UT_SpinDensityMatrix.size())
    {
        UT_SpinDensityMatrix[nr] = value;
        return true;
    }
    else
    {
        std::cout << "invalid arguments for set_SDM! Input was: " << nr << ";" << value << std::endl;
        return false;
    }
};

bool WFN::build_DM(std::string basis_set_path, bool debug) {
    using namespace std;
    int elcount = -get_charge();
    if (debug)
       std::cout << "elcount: " << elcount << std::endl;
    for (int i = 0; i < ncen; i++)
    {
        elcount += get_atom_charge(i);
        elcount -= constants::ECP_electrons_pTB[get_atom_charge(i)];
    }
    if (debug)
       std::cout << "elcount after: " << elcount << std::endl;
    int alpha_els = 0, beta_els = 0, temp_els = elcount;
    while (temp_els > 1)
    {
        alpha_els++;
        beta_els++;
        temp_els -= 2;
        if (debug)
           std::cout << temp_els << std::endl;
        err_checkf(alpha_els >= 0 && beta_els >= 0, "Error setting alpha and beta electrons! a or b are negative!",std::cout);
        err_checkf(alpha_els + beta_els <= elcount, "Error setting alpha and beta electrons! Sum a + b > elcount!",std::cout);
        err_checkf(temp_els > -elcount, "Error setting alpha and beta electrons! Ran below -elcount!",std::cout);
    }
    alpha_els += temp_els;
    if (debug)
       std::cout << "al/be els:" << alpha_els << " " << beta_els << std::endl;
    const int mult = get_multi();
    int diff = 0;
    if (mult != 0)
        diff = get_multi() - 1;
    if (debug)
       std::cout << "diff: " << diff << std::endl;
    while (alpha_els - beta_els != diff)
    {
        alpha_els++;
        beta_els--;
        err_checkf(alpha_els >= 0 && beta_els >= 0, "Error setting alpha and beta electrons!",std::cout);
    }
    if (debug)
    {
       std::cout << "alpha, beta, elcount: " << setw(5) << alpha_els << setw(5) << beta_els << setw(5) << elcount << endl;
    }
    if (get_nr_basis_set_loaded() == 0)
    {
        if (debug)
           std::cout << "No basis set loaded, will load a complete basis set now!" << endl;
        err_checkf(read_basis_set_vanilla(basis_set_path, *this, debug), "ERROR during reading of missing basis set!",std::cout);
    }
    else if (get_nr_basis_set_loaded() < get_ncen())
    {
       std::cout << "Not all atoms have a basis set loaded!\nLaoding the missing atoms..." << flush;
        err_checkf(read_basis_set_missing(basis_set_path, *this, debug), "ERROR during reading of missing basis set!",std::cout);
    }
    else if (get_nr_basis_set_loaded() > get_ncen())
    {
        err_checkf(false, "# of loaded > # atoms\nSorry, this should not happen... aborting!!!",std::cout);
    }
    // set_modified();
    vec CMO;
    vec CMO_beta;
    if (debug)
    {
       std::cout << "Origin: " << get_origin() << endl;
    }
    if (get_origin() == 2 || get_origin() == 4 || get_origin() == 9 || get_origin() == 8)
    {
        //-----------------------check ordering and order accordingly----------------------
        sort_wfn(check_order(debug), debug);
        //---------------normalize basis set---------------------------------
        if (debug)
           std::cout << "starting to normalize the basis set" << endl;
        vec norm_const;
        //-----------debug output---------------------------------------------------------
        if (debug)
        {
           std::cout << "exemplary output before norm_const of the first atom with all it's properties: " << endl;
            print_atom_long(0);
           std::cout << "ended normalizing the basis set, now for the MO_coeffs" << endl;
           std::cout << "Status report:" << endl;
           std::cout << "size of norm_const: " << norm_const.size() << endl;
           std::cout << "WFN MO counter: " << get_nmo() << endl;
           std::cout << "Number of atoms: " << get_ncen() << endl;
           std::cout << "Primitive count of zero MO: " << get_MO_primitive_count(0) << endl;
           std::cout << "Primitive count of first MO: " << get_MO_primitive_count(1) << endl;
        }

        //-------------------normalize the basis set shell wise into a copy vector---------
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
                   std::cout << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy..." << endl;
                    break;
                }
                temp_c = pow(temp_c, 0.25) * get_atom_basis_set_coefficient(a, p);
                if (debug)
                   std::cout << "temp_c:" << temp_c << std::endl;
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
                err_chkf(type_temp != -1, "ERROR in type assignement!!",std::cout);
                if (debug)
                {
                   std::cout << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << endl
                        << "start: " << get_shell_start(a, s)
                        << " stop: " << get_shell_end(a, s) << endl
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
                            double term = constants::PI3 / pow(aiaj, 3);
                            term = pow(term, 0.5);
                            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
                        }
                    }
                    if (factor == 0)
                        return false;
                    factor = pow(factor, -0.5);
                    if (debug)
                       std::cout << factor << endl;
                    for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
                    {
                        if (debug)
                        {
                           std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i)
                                << " Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << endl;
                        }
                        // contraction_coefficients[a][i] = factor * get_atom_basis_set_coefficient(a, i);
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
                            aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
                            double term = constants::PI3 / (4 * pow(aiaj, 5));
                            term = pow(term, 0.5);
                            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
                        }
                    }
                    if (factor == 0)
                        return false;
                    factor = pow(factor, -0.5);
                    if (debug)
                       std::cout << factor << endl;
                    for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
                    {
                        if (debug)
                        {
                           std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i)
                                << " Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << endl;
                        }
                        // contraction_coefficients[a][i] = factor * get_atom_basis_set_coefficient(a, i);
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
                            aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
                            double term = constants::PI3 / (16 * pow(aiaj, 7));
                            term = pow(term, 0.5);
                            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
                        }
                    }
                    if (factor == 0)
                        return false;
                    factor = (pow(factor, -0.5)) / sqrt(3);
                    if (debug)
                       std::cout << factor << endl;
                    for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
                    {
                        if (debug)
                        {
                           std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i)
                                << " Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << endl;
                        }
                        // contraction_coefficients[a][i] = factor * get_atom_basis_set_coefficient(a, i);
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
                            aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
                            double term = constants::PI3 / (64 * pow((aiaj), 9));
                            term = pow(term, 0.5);
                            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
                        }
                    }
                    if (factor == 0)
                        return false;
                    factor = pow(factor, -0.5) / sqrt(15);
                    if (debug)
                       std::cout << factor << endl;
                    for (int i = get_shell_start(a, s); i <= get_shell_end(a, s); i++)
                    {
                        if (debug)
                        {
                           std::cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i)
                                << " Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << endl;
                        }
                        // contraction_coefficients[a][i] = factor * get_atom_basis_set_coefficient(a, i);
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
                   std::cout << "This shell has: " << get_shell_end(a, s) - get_shell_start(a, s) + 1 << " primitives" << endl;
            }
        }
        //-----------debug output---------------------------------------------------------
        if (debug)
        {
           std::cout << "exemplary output of the first atom with all it's properties: " << endl;
            print_atom_long(0);
           std::cout << "ended normalizing the basis set, now for the norm_cprims" << endl;
           std::cout << "Status report:" << endl;
           std::cout << "size of norm_const: " << norm_const.size() << endl;
           std::cout << "WFN MO counter: " << get_nmo() << endl;
           std::cout << "Number of atoms: " << get_ncen() << endl;
           std::cout << "Primitive count of zero MO: " << get_MO_primitive_count(0) << endl;
           std::cout << "Primitive count of first MO: " << get_MO_primitive_count(1) << endl;
        }
        //---------------------To not mix up anything start normalizing WFN_matrix now--------------------------
        int run = 0;
        vec2 changed_coefs;
        changed_coefs.resize(get_nmo());
        if (debug)
        {
           std::cout << "Opening norm_cprim!" << endl;
            ofstream norm_cprim("norm_prim.debug", ofstream::out);
            for (int m = 0; m < get_nmo(); m++)
            {
                norm_cprim << m << ". MO:" << endl;
                changed_coefs[m].resize(get_MO_primitive_count(m), 0.0);
                for (int p = 0; p < get_MO_primitive_count(m); p++)
                {
                    changed_coefs[m][p] = get_MO_coef(m, p) / norm_const[p];
                    if (m == 0)
                       std::cout << p << ". primitive; " << m << ". MO "
                        << "norm nonst: " << norm_const[p]
                        << " temp after normalization: " << changed_coefs[m][p] << "\n";
                    norm_cprim << " " << changed_coefs[m][p] << endl;
                    run++;
                }
            }
            norm_cprim.flush();
            norm_cprim.close();
           std::cout << "See norm_cprim.debug for the CPRIM vectors" << endl;
           std::cout << "Total count in CPRIM: " << run << endl;
        }
        else
        {
#pragma omp parallel for
            for (int m = 0; m < get_nmo(); m++)
            {
                changed_coefs[m].resize(get_MO_primitive_count(m), 0.0);
                for (int p = 0; p < get_MO_primitive_count(m); p++)
                {
                    changed_coefs[m][p] = get_MO_coef(m, p) / norm_const[p];
                }
            }
        }
        //--------------Build CMO of alessandro from the first elements of each shell-------------
        int nao = 0;
        for (int a = 0; a < get_ncen(); a++)
        {
            for (int s = 0; s < get_atom_shell_count(a); s++)
            {
                switch (get_shell_type(a, s))
                {
                case 1:
                    nao++;
                    break;
                case 2:
                    nao += 3;
                    break;
                case 3:
                    nao += 6;
                    break;
                case 4:
                    nao += 10;
                    break;
                }
            }
        }
        int nshell = 0;
        for (int m = 0; m < get_nmo(); m++)
        {
            int run_2 = 0;
            for (int a = 0; a < get_ncen(); a++)
            {
                for (int s = 0; s < get_atom_shell_count(a); s++)
                {
                    // if (debug)std::cout << "Going to load the " << get_shell_start_in_primitives(a, s) << ". value\n"l;
                    switch (get_shell_type(a, s))
                    {
                    case 1:
                        CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s)]);
                        if (debug && get_atom_shell_primitives(a, s) != 1 && m == 0)
                           std::cout << "Pushing back 1 coefficient for S shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives! Shell start is: " << get_shell_start(a, s) << endl;
                        break;
                    case 2:
                        for (int i = 0; i < 3; i++)
                            CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
                        if (debug && get_atom_shell_primitives(a, s) != 1 && m == 0)
                           std::cout << "Pushing back 3 coefficients for P shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!" << endl;
                        break;
                    case 3:
                        for (int i = 0; i < 6; i++)
                            CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
                        if (debug && get_atom_shell_primitives(a, s) != 1 && m == 0)
                           std::cout << "Pushing back 6 coefficient for D shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!" << endl;
                        break;
                    case 4:
                        // this hardcoded piece is due to the order of f-type functions in the fchk
                        for (int i = 0; i < 3; i++)
                            CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
                        CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 6]);
                        for (int i = 0; i < 2; i++)
                            CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i + 3]);
                        for (int i = 0; i < 2; i++)
                            CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i + 7]);
                        CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 5]);
                        CMO.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 9]);
                        if (debug && get_atom_shell_primitives(a, s) != 1 && m == 0)
                           std::cout << "Pushing back 10 coefficient for F shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!" << endl;
                        break;
                    }
                    run_2++;
                }
                if (debug && m == 0)
                   std::cout << "finished with atom!" << endl;
            }
            if (debug)
               std::cout << "finished with MO!" << endl;
            if (nshell != run_2)
                nshell = run_2;
        }
        if (alpha_els != beta_els)
        {
            for (int m = alpha_els; m < alpha_els + beta_els; m++)
            {
                int run_2 = 0;
                for (int a = 0; a < get_ncen(); a++)
                {
                    for (int s = 0; s < get_atom_shell_count(a); s++)
                    {
                        if (debug)
                           std::cout << "Going to load the " << get_shell_start_in_primitives(a, s) << ". value" << endl;
                        switch (get_shell_type(a, s))
                        {
                        case 1:
                            CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s)]);
                            if (m == 0)
                                nao++;
                            if (debug && get_atom_shell_primitives(a, s) != 1)
                               std::cout << "Pushing back 1 coefficient for S shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives! Shell start is: " << get_shell_start(a, s) << endl;
                            break;
                        case 2:
                            for (int i = 0; i < 3; i++)
                                CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
                            if (debug && get_atom_shell_primitives(a, s) != 1)
                               std::cout << "Pushing back 3 coefficients for P shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!" << endl;
                            if (m == 0)
                                nao += 3;
                            break;
                        case 3:
                            for (int i = 0; i < 6; i++)
                                CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
                            if (debug && get_atom_shell_primitives(a, s) != 1)
                               std::cout << "Pushing back 6 coefficient for D shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!" << endl;
                            if (m == 0)
                                nao += 6;
                            break;
                        case 4:
                            // this hardcoded piece is due to the order of f-type functions in the fchk
                            for (int i = 0; i < 3; i++)
                                CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i]);
                            CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 6]);
                            for (int i = 0; i < 2; i++)
                                CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i + 3]);
                            for (int i = 0; i < 2; i++)
                                CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + i + 7]);
                            CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 5]);
                            CMO_beta.push_back(changed_coefs[m][get_shell_start_in_primitives(a, s) + 9]);
                            if (debug && get_atom_shell_primitives(a, s) != 1)
                               std::cout << "Pushing back 10 coefficient for F shell, this shell has " << get_atom_shell_primitives(a, s) << " primitives!" << endl;
                            if (m == 0)
                                nao += 10;
                            break;
                        }
                        run_2++;
                    }
                    if (debug)
                       std::cout << "finished with atom!" << endl;
                }
                if (debug)
                   std::cout << "finished with MO!" << endl;
                if (nshell != run_2)
                    nshell = run_2;
            }
        }

        if (debug)
        {
            ofstream cmo("cmo.debug", ofstream::out);
            for (int p = 0; p < CMO.size(); p++)
            {
                for (int i = 0; i < 5; i++)
                {
                    cmo << scientific << setw(14) << setprecision(7) << CMO[p + i] << " ";
                }
                p += 4;
                cmo << endl;
            }
            cmo.flush();
            cmo.close();
           std::cout << CMO.size() << " Elements in CMO" << endl;
           std::cout << norm_const.size() << " = nprim" << endl;
           std::cout << nao << " = nao" << endl;
           std::cout << nshell << " = nshell" << endl;
        }
        //------------------ make the DM -----------------------------
        int naotr = nao * (nao + 1) / 2;
        vec kp;
        resize_DM(naotr, 0.0);
        if (alpha_els != beta_els)
            resize_SDM(naotr, 0.0);
        if (debug)
        {
           std::cout << "I made kp!" << endl
                << nao << " is the maximum for iu" << endl;
           std::cout << "Making DM now!" << endl;
        }
        for (int iu = 0; iu < nao; iu++)
        {
#pragma omp parallel for
            for (int iv = 0; iv <= iu; iv++)
            {
                const int iuv = (iu * (iu + 1) / 2) + iv;
                // if (debug)std::cout << "iu: " << iu << " iv: " << iv << " iuv: " << iuv << " kp(iu): " << iu * (iu + 1) / 2 << endl;
                double temp;
                // if (debug)std::cout << "Working on MO: ";
                for (int m = 0; m < get_nmo(); m++)
                {
                    // if (debug && m == 0)std::cout << m << " " << flush;
                    // else if (debug && m != get_nmo() - 1)std::cout << "." << flush;
                    // elsestd::cout << get_nmo() - 1 << flush;
                    if (alpha_els != beta_els)
                    {
                        if (m < alpha_els)
                        {
                            temp = get_MO_occ(m) * CMO[iu + (m * nao)] * CMO[iv + (m * nao)];
                            err_checkf(set_SDM(iuv, get_SDM(iuv) + temp), "Something went wrong while writing the SDM! iuv=" + to_string(iuv), std::cout);
                            err_checkf(set_DM(iuv, get_DM(iuv) + temp), "Something went wrong while writing the DM! iuv=" + to_string(iuv), std::cout);
                        }
                        else
                        {
                            temp = get_MO_occ(m) * CMO_beta[iu + ((m - alpha_els) * nao)] * CMO_beta[iv + ((m - alpha_els) * nao)];
                            err_checkf(set_SDM(iuv, get_SDM(iuv) - temp), "Something went wrong while writing the SDM! iuv=" + to_string(iuv), std::cout);
                            err_checkf(set_DM(iuv, get_DM(iuv) + temp), "Something went wrong while writing the DM! iuv=" + to_string(iuv), std::cout);
                        }
                    }
                    else
                    {
                        if (get_MO_occ(m) == 0.0)
                            continue;
                        temp = get_MO_occ(m) * CMO[iu + (m * nao)] * CMO[iv + (m * nao)];
                        err_checkf(set_DM(iuv, get_DM(iuv) + temp), "Something went wrong while writing the DM!",std::cout);
                    }
                    // else if (debug)std::cout << "DM after: " << get_DM(iuv) << endl;
                }
                // if (debug)std::cout << endl;
            }
        }
    }
    else
    {
       std::cout << "Sorry, this origin is not supported yet!" << endl;
        return false;
    }
    return true;
};

int WFN::check_order(const bool &debug) const
{
    for (int i = 0; i < ncen; i++)
    {
        if (!get_atom_basis_set_loaded(i))
        {
            std::cout << "Sorry, consistency check only works if basis set is loaded for all atoms!" << std::endl
                 << "Failing atom: " << i << " " << get_atom_label(i) << std::endl;
            return -1;
        }
    }
    //---------------------------check if the wfn is in the right order----------------------
    int order = 0;   // 1=gaussian (P=2222 3333 4444) 2=tonto (234 234 234 234) 3=ORCA (423 423 423 423)
    int f_order = 0; // 1=gaussian (F=11 12 13 17 14 15 18 19 16 20) 2=tonto=ORCA 3=ORCA (11 12 13 14 15 17 16 18 19 20) 4=natural (11 12 13 14 15 16 17 18 19 20)
    int primcounter = 0;
    bool order_found = false;
    for (int a = 0; a < get_ncen(); a++)
    {
        for (int s = 0; s < get_atom_shell_count(a); s++)
        {
            int type = get_shell_type(a, s);
            switch (type)
            {
            case 1: // S Orbital
            {
                for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                {
                    if (types[primcounter] != 1)
                    {
                        order = -1;
                        if (debug)
                        {
                            std::cout << "This should not happen, the order of your file is not ok for S-types! Checked #:" << primcounter << std::endl;
                        }
                    }
                    else
                        primcounter++;
                }
                break;
            }
            case 2: // P Orbital
            {
                if (order_found)
                {
                    if (order == 1)
                    {
                        for (int r = 0; r < get_atom_shell_primitives(a, s); r++)
                        {
                            if (types[primcounter] != 2 || types[primcounter + get_atom_shell_primitives(a, s)] != 3 || types[primcounter + 2 * get_atom_shell_primitives(a, s)] != 4)
                            {
                                if (debug)
                                {
                                    std::cout << "The found order does not match all type entries! primcounter: " << primcounter << std::endl;
                                }
                            }
                            primcounter++;
                        }
                        primcounter += 2 * get_atom_shell_primitives(a, s);
                    }
                    else if (order == 2)
                    {
                        for (int r = 0; r < get_atom_shell_primitives(a, s); r++)
                        {
                            if (types[primcounter] != 2 || types[primcounter + 1] != 3 || types[primcounter + 2] != 4)
                            {
                                if (debug)
                                {
                                    std::cout << "The found order does not match all type entries! primcounter: " << primcounter << std::endl;
                                }
                                return -1;
                            }
                            primcounter += 3;
                        }
                    }
                    else if (order == 3)
                    {
                        for (int r = 0; r < get_atom_shell_primitives(a, s); r++)
                        {
                            if (types[primcounter] != 4 || types[primcounter + 1] != 2 || types[primcounter + 2] != 3)
                            {
                                if (debug)
                                {
                                    std::cout << "The found order does not match all type entries! primcounter: " << primcounter << std::endl;
                                }
                                return -1;
                            }
                            primcounter += 3;
                        }
                    }
                }
                else
                {
                    if (types[primcounter] == 2)
                    {
                        if (debug && a == 0)
                        {
                            std::cout << "Seems to be either tonto or gaussian file..." << std::endl;
                        }
                        if (types[primcounter + 1] == 3)
                        {
                            order = 2;
                            order_found = true;
                            if (debug)
                            {
                                std::cout << "This wfn file is in tonto order!" << std::endl;
                            }
                        }
                        else if (types[primcounter + 1] == 2 && get_atom_shell_primitives(a, s) > 1)
                        {
                            order = 1;
                            order_found = true;
                            if (debug)
                            {
                                std::cout << "This wfn file is in gaussian order!" << std::endl;
                            }
                        }
                        else
                        {
                            order = 1;
                            order_found = true;
                            if (debug)
                            {
                                std::cout << "Either some error or this shell only has 1 p-primitive and "
                                    << "i didn't find any order yet... assuming gaussian" << std::endl;
                            }
                        }
                    }
                    else if (types[primcounter] == 4)
                    {
                        if (debug && a == 0)
                        {
                            std::cout << "Seems as if this file was ORCA ordered..." << std::endl;
                        }
                        order = 3;
                        if (types[primcounter + 1] == 2)
                        {
                            if (debug)
                            {
                                std::cout << "Yep, that's right! making it permanent now!" << std::endl;
                            }
                            order_found = true;
                        }
                    }
                    else
                    {
                        std::cout << "I can't recognize this order of the .wfn file..." << std::endl;
                        return -1;
                    }
                    s--;
                }
                break;
            }
            case 3: // D Orbital
            {
                if (order_found)
                {
                    switch (order)
                    {
                    case 1:
                    {
                        for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                        {
                            if (types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] != 5 || types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] != 6 || types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] != 7 || types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] != 8 || types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] != 9 || types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] != 10)
                            {
                                order = -1;
                                if (debug)
                                {
                                    std::cout << "The checked types are 6 from #" << primcounter << " and are:" << std::endl
                                         << types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] << std::endl;
                                }
                                std::cout << "Something seems to be wrong in the order of your D-Types..." << std::endl;
                            }
                            else
                                primcounter += 6;
                        }
                    }
                    break;
                    case 2:
                    case 3:
                    {
                        for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                        {
                            if (types[get_shell_start_in_primitives(a, s) + 0 + 6 * i] != 5 || types[get_shell_start_in_primitives(a, s) + 1 + 6 * i] != 6 || types[get_shell_start_in_primitives(a, s) + 2 + 6 * i] != 7 || types[get_shell_start_in_primitives(a, s) + 3 + 6 * i] != 8 || types[get_shell_start_in_primitives(a, s) + 4 + 6 * i] != 9 || types[get_shell_start_in_primitives(a, s) + 5 + 6 * i] != 10)
                            {
                                order = -1;
                                if (debug)
                                {
                                    std::cout << "The checked types are 6 from #" << primcounter << " and are:" << std::endl
                                         << types[get_shell_start_in_primitives(a, s) + 0 + 6 * i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 1 + 6 * i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 2 + 6 * i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 3 + 6 * i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 4 + 6 * i] << " "
                                         << types[get_shell_start_in_primitives(a, s) + 5 + 6 * i] << std::endl;
                                }
                                std::cout << "Something seems to be wrong in the order of your D-Types..." << std::endl;
                            }
                            else
                                primcounter += 6;
                        }
                    }
                    break;
                    }
                }
                else
                {
                    std::cout << "That's highly suspicious, no order for P-type but D-type functions? Let me stop before i do something stupid!" << std::endl;
                    return -1;
                }
            }
            break;
            case 4:
                for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                {
                    if (types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] != 11 || types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] != 12 || types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] != 13 || types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] != 17 || types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] != 14 || types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] != 15 || types[get_shell_start_in_primitives(a, s) + 6 * get_atom_shell_primitives(a, s) + i] != 18 || types[get_shell_start_in_primitives(a, s) + 7 * get_atom_shell_primitives(a, s) + i] != 19 || types[get_shell_start_in_primitives(a, s) + 8 * get_atom_shell_primitives(a, s) + i] != 16 || types[get_shell_start_in_primitives(a, s) + 9 * get_atom_shell_primitives(a, s) + i] != 20)
                    {
                        if (types[primcounter] == 11 || types[primcounter + 1] == 12 || types[primcounter + 2] == 13 || types[primcounter + 3] == 14 || types[primcounter + 4] == 15 || types[primcounter + 5] == 16 || types[primcounter + 6] == 17 || types[primcounter + 7] == 18 || types[primcounter + 8] == 19 || types[primcounter + 9] == 20)
                        {
                            f_order = 4;
                            if (debug)
                            {
                                std::cout << "The checked types are 10 from #" << primcounter << " and are:\n"
                                     << types[primcounter] << " " << types[primcounter + 1] << " " << types[primcounter + 2] << " "
                                     << types[primcounter + 3] << " " << types[primcounter + 4] << " " << types[primcounter + 5] << " "
                                     << types[primcounter + 6] << " " << types[primcounter + 7] << " " << types[primcounter + 8] << " "
                                     << types[primcounter + 9] << std::endl
                                     << "Appears to be already okay..." << std::endl;
                            }
                            primcounter += 10;
                        }
                        else if (types[primcounter] != 11 || types[primcounter + 1] != 12 || types[primcounter + 2] != 13 || types[primcounter + 3] != 14 || types[primcounter + 4] != 15 || types[primcounter + 5] != 17 || types[primcounter + 6] != 16 || types[primcounter + 7] != 18 || types[primcounter + 8] != 19 || types[primcounter + 9] != 20)
                        {
                            order = -1;
                            if (debug)
                            {
                                std::cout << "The checked types are 10 from #" << primcounter << " and are:\n"
                                    << types[primcounter] << " " << types[primcounter + 1] << " " << types[primcounter + 2] << " "
                                    << types[primcounter + 3] << " " << types[primcounter + 4] << " " << types[primcounter + 5] << " "
                                    << types[primcounter + 6] << " " << types[primcounter + 7] << " " << types[primcounter + 8] << " "
                                    << types[primcounter + 9] << "\n"
                                    << "Something seems to be wrong in the order of your F-Types..." << std::endl;
                            }
                        }
                        else
                        {
                            f_order = 3;
                            primcounter += 10;
                        }
                    }
                    else
                    {
                        f_order = 1;
                        primcounter += 10;
                    }
                }
                break;
            default:
                std::cout << "ERROR in type assignement!!" << std::endl;
                return -1;
                break;
            }
        }
    }
    /*if(debug){
       std::cout << "Going to return " << f_order*10+order << endl;
        Enter();
    }*/
    return (f_order * 10 + order);
};

bool WFN::sort_wfn(const int &g_order, const bool &debug)
{
    set_modified();
    int primcounter = 0;
    int f_order = 0;
    int order = g_order;
    // Sorry for this way of forwarding the order, i think 2 switches would have been more nicely
    while (order >= 10)
    {
        f_order++;
        order -= 10;
    }
    switch (order)
    {
    case 1:
        for (int a = 0; a < get_ncen(); a++)
            for (int s = 0; s < get_atom_shell_count(a); s++)
                switch (get_shell_type(a, s))
                {
                case 1:
                    primcounter += get_atom_shell_primitives(a, s);
                    break;
                case 2:
                {
                    if (get_atom_shell_primitives(a, s) > 1)
                    {
                        ivec temp_centers;
                        temp_centers.resize(3 * get_atom_shell_primitives(a, s));
                        ivec temp_types;
                        temp_types.resize(3 * get_atom_shell_primitives(a, s));
                        vec temp_exponents;
                        temp_exponents.resize(3 * get_atom_shell_primitives(a, s));
                        vec2 temp_MO_coefficients;
                        temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 3);
                        for (int i = 0; i < get_atom_shell_primitives(a, s) * 3; i++)
                            temp_MO_coefficients[i].resize(nmo);
                        for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                            for (int c = 0; c < 3; c++)
                            {
                                temp_centers[3 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
                                temp_types[3 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
                                temp_exponents[3 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
                                for (int m = 0; m < nmo; m++)
                                    temp_MO_coefficients[3 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
                            }
                        for (int i = 0; i < 3 * get_atom_shell_primitives(a, s); i++)
                        {
                            centers[primcounter + i] = temp_centers[i];
                            types[primcounter + i] = temp_types[i];
                            exponents[primcounter + i] = temp_exponents[i];
                            for (int m = 0; m < nmo; m++)
                                err_checkf(set_MO_coef(m, primcounter + i, temp_MO_coefficients[i][m]), "Error while assigning new MO coefficient!", std::cout);
                        }
                    }
                    primcounter += get_atom_shell_primitives(a, s) * 3;
                    break;
                }
                case 3:
                {
                    if (get_atom_shell_primitives(a, s) > 1)
                    {
                        ivec temp_centers;
                        temp_centers.resize(6 * get_atom_shell_primitives(a, s));
                        ivec temp_types;
                        temp_types.resize(6 * get_atom_shell_primitives(a, s));
                        vec temp_exponents;
                        temp_exponents.resize(get_atom_shell_primitives(a, s) * 6);
                        vec2 temp_MO_coefficients;
                        temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 6);
                        for (int i = 0; i < get_atom_shell_primitives(a, s) * 6; i++)
                            temp_MO_coefficients[i].resize(nmo);
                        for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                            for (int c = 0; c < 6; c++)
                            {
                                temp_centers[6 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
                                temp_types[6 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
                                temp_exponents[6 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
                                for (int m = 0; m < nmo; m++)
                                    temp_MO_coefficients[6 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
                            }
                        for (int i = 0; i < 6 * get_atom_shell_primitives(a, s); i++)
                        {
                            centers[primcounter + i] = temp_centers[i];
                            types[primcounter + i] = temp_types[i];
                            exponents[primcounter + i] = temp_exponents[i];
                            for (int m = 0; m < nmo; m++)
                                if (!set_MO_coef(m, primcounter + i, temp_MO_coefficients[i][m]))
                                {
                                    std::cout << "Error while assigning new MO coefficient!" << std::endl;
                                    return false;
                                }
                        }
                    }
                    primcounter += get_atom_shell_primitives(a, s) * 6;
                    break;
                }
                case 4:
                {
                    if (get_atom_shell_primitives(a, s) > 1)
                    {
                        ivec temp_centers;
                        temp_centers.resize(10 * get_atom_shell_primitives(a, s));
                        ivec temp_types;
                        temp_types.resize(10 * get_atom_shell_primitives(a, s));
                        vec temp_exponents;
                        temp_exponents.resize(get_atom_shell_primitives(a, s) * 10);
                        vec2 temp_MO_coefficients;
                        temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 10);
                        for (int i = 0; i < get_atom_shell_primitives(a, s) * 10; i++)
                            temp_MO_coefficients[i].resize(nmo);
                        for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                            for (int c = 0; c < 10; c++)
                            {
                                temp_centers[10 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
                                temp_types[10 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
                                temp_exponents[10 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
                                for (int m = 0; m < nmo; m++)
                                    temp_MO_coefficients[10 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
                            }
                        for (int i = 0; i < 10 * get_atom_shell_primitives(a, s); i++)
                        {
                            centers[primcounter + i] = temp_centers[i];
                            types[primcounter + i] = temp_types[i];
                            exponents[primcounter + i] = temp_exponents[i];
                            for (int m = 0; m < nmo; m++)
                                if (!set_MO_coef(m, primcounter + i, temp_MO_coefficients[i][m]))
                                {
                                    std::cout << "Error while assigning new MO coefficient!" << std::endl;
                                    return false;
                                }
                        }
                    }
                    primcounter += get_atom_shell_primitives(a, s) * 10;
                    break;
                }
                }
        break;
    case 2:
        if (debug)
        {
            std::cout << "Nothing to be done here, i like tonto type..." << std::endl;
        }
        break;
    case 3:
        for (int a = 0; a < get_ncen(); a++)
        {
            for (int s = 0; s < get_atom_shell_count(a); s++)
            {
                switch (get_shell_type(a, s))
                {
                case 1:
                    primcounter += get_atom_shell_primitives(a, s);
                    break;
                case 2:
                {
                    int temp_center;
                    int temp_type;
                    double temp_exponent;
                    vec temp_MO_coefficients;
                    temp_MO_coefficients.resize(nmo);
                    for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                    {
                        temp_center = centers[primcounter];
                        temp_type = types[primcounter];
                        temp_exponent = exponents[primcounter];
                        for (int m = 0; m < nmo; m++)
                            temp_MO_coefficients[m] = MOs[m].get_coefficient(primcounter);
                        for (int j = 0; j < 2; j++)
                        {
                            centers[primcounter + j] = centers[primcounter + 1 + j];
                            types[primcounter + j] = types[primcounter + 1 + j];
                            exponents[primcounter + j] = exponents[primcounter + 1 + j];
                            for (int m = 0; m < nmo; m++)
                            {
                                err_checkf(set_MO_coef(m, primcounter + j, MOs[m].get_coefficient(primcounter + 1 + j)), "Error while assigning new MO coefficient!", std::cout);
                            }
                        }
                        centers[primcounter + 2] = temp_center;
                        types[primcounter + 2] = temp_type;
                        exponents[primcounter + 2] = temp_exponent;
                        for (int m = 0; m < nmo; m++)
                            if (!set_MO_coef(m, primcounter + 2, temp_MO_coefficients[m]))
                            {
                                std::cout << "Error while assigning new MO coefficient!" << std::endl;
                                return false;
                            }
                        primcounter += 3;
                    }
                    break;
                }
                case 3:
                    primcounter += get_atom_shell_primitives(a, s) * 6;
                    break;
                case 4:
                    primcounter += get_atom_shell_primitives(a, s) * 10;
                    break;
                }
            }
        }
        break;
    default:
        std::cout << "order type: " << f_order << " " << order << " not supported!" << std::endl;
        return false;
        break;
    }
    primcounter = 0;
    switch (f_order)
    {
    case 1:
        for (int a = 0; a < get_ncen(); a++)
            for (int s = 0; s < get_atom_shell_count(a); s++)
                switch (get_shell_type(a, s))
                {
                case 1:
                    primcounter += get_atom_shell_primitives(a, s);
                    break;
                case 2:
                    primcounter += get_atom_shell_primitives(a, s) * 3;
                    break;
                case 3:
                    primcounter += get_atom_shell_primitives(a, s) * 6;
                    break;
                case 4:
                {
                    ivec temp_center;
                    temp_center.resize(10);
                    ivec temp_type;
                    temp_type.resize(10);
                    vec temp_exponent;
                    temp_exponent.resize(10);
                    vec2 temp_MO_coefficients;
                    temp_MO_coefficients.resize(nmo);
                    for (int m = 0; m < nmo; m++)
                        temp_MO_coefficients[m].resize(10);
                    for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                    {
                        for (int j = 0; j < 10; j++)
                        {
                            temp_center[j] = centers[get_shell_start_in_primitives(a, s) + 10 * i + j];
                            temp_type[j] = types[get_shell_start_in_primitives(a, s) + 10 * i + j];
                            temp_exponent[j] = exponents[get_shell_start_in_primitives(a, s) + 10 * i + j];
                            for (int m = 0; m < nmo; m++)
                                temp_MO_coefficients[m][j] = MOs[m].get_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + j);
                        }
                        ivec mask{0, 1, 2, 6, 3, 4, 7, 8, 5, 9};
                        for (int j = 0; j < 10; j++)
                        {
                            centers[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_center[mask[j]];
                            types[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_type[mask[j]];
                            exponents[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_exponent[mask[j]];
                            for (int m = 0; m < nmo; m++)
                                set_MO_coef(m, get_shell_start_in_primitives(a, s) + 10 * i + mask[j], temp_MO_coefficients[m][j]);
                        }
                    }
                    primcounter += 10;
                }
                break;
                }
        break;
    case 2:
    case 3:
        for (int a = 0; a < get_ncen(); a++)
            for (int s = 0; s < get_atom_shell_count(a); s++)
                switch (get_shell_type(a, s))
                {
                case 1:
                    primcounter += get_atom_shell_primitives(a, s);
                    break;
                case 2:
                    primcounter += get_atom_shell_primitives(a, s) * 3;
                    break;
                case 3:
                    primcounter += get_atom_shell_primitives(a, s) * 6;
                    break;
                case 4:
                {
                    ivec temp_center;
                    temp_center.resize(10);
                    ivec temp_type;
                    temp_type.resize(10);
                    vec temp_exponent;
                    temp_exponent.resize(10);
                    vec2 temp_MO_coefficients;
                    temp_MO_coefficients.resize(nmo);
                    for (int m = 0; m < nmo; m++)
                        temp_MO_coefficients[m].resize(10);
                    for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
                    {
                        for (int j = 0; j < 10; j++)
                        {
                            temp_center[j] = centers[get_shell_start_in_primitives(a, s) + 10 * i + j];
                            temp_type[j] = types[get_shell_start_in_primitives(a, s) + 10 * i + j];
                            temp_exponent[j] = exponents[get_shell_start_in_primitives(a, s) + 10 * i + j];
                            for (int m = 0; m < nmo; m++)
                                temp_MO_coefficients[m][j] = MOs[m].get_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + j);
                        }
                        ivec mask{0, 1, 2, 3, 4, 6, 5, 7, 8, 9};
                        for (int j = 0; j < 10; j++)
                        {
                            centers[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_center[mask[j]];
                            types[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_type[mask[j]];
                            exponents[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_exponent[mask[j]];
                            for (int m = 0; m < nmo; m++)
                                set_MO_coef(m, get_shell_start_in_primitives(a, s) + 10 * i + mask[j], temp_MO_coefficients[m][j]);
                        }
                    }
                    primcounter += 10;
                }
                break;
                }
        break;
    case 4:
        if (debug)
        {
            std::cout << "This is fine, i like them well ordered!" << std::endl;
        }
        break;
    case 0:
        if (debug)
        {
            std::cout << "There seems to be no f-functions!" << std::endl;
        }
        break;
    default:
        std::cout << "f-order type: " << f_order << " not supported!" << std::endl;
        return false;
        break;
    }
    return true;
};

void WFN::set_has_ECPs(const bool &in, const bool &apply_to_atoms, const int &ECP_mode)
{
    has_ECPs = in;
    ECP_m = ECP_mode;
    if (apply_to_atoms && ECP_mode == 1) // This is the def2 ECPs
    {
#pragma omp parallel for
        for (int i = 0; i < ncen; i++)
        {
            if (constants::ECP_electrons[get_atom_charge(i)] != 0)
            {
                atoms[i].set_ECP_electrons(constants::ECP_electrons[get_atom_charge(i)]);
            }
        }
    }
    if (apply_to_atoms && ECP_mode == 2) // xTB ECPs
    {
#pragma omp parallel for
        for (int i = 0; i < ncen; i++)
        {
            if (constants::ECP_electrons_xTB[get_atom_charge(i)] != 0)
            {
                atoms[i].set_ECP_electrons(constants::ECP_electrons_xTB[get_atom_charge(i)]);
            }
        }
    }
    if (apply_to_atoms && ECP_mode == 3) // pTB ECPs
    {
#pragma omp parallel for
        for (int i = 0; i < ncen; i++)
        {
            if (constants::ECP_electrons_pTB[get_atom_charge(i)] != 0)
            {
                atoms[i].set_ECP_electrons(constants::ECP_electrons_pTB[get_atom_charge(i)]);
            }
        }
    }
};

void WFN::set_ECPs(ivec &nr, ivec &elcount)
{
    has_ECPs = true;
    err_chkf(nr.size() == elcount.size(), "mismatch in size of atoms and ECP electrons!", std::cout);
#pragma omp parallel for
    for (int i = 0; i < ncen; i++)
    {
        for (int j = 0; j < nr.size(); j++)
        {
            std::cout << "checking " << get_atom_charge(i) << " against " << nr[j] << std::endl;
            if (get_atom_charge(i) == nr[j])
            {
                std::cout << "Adding " << elcount[j] << " electron to atom " << i << "with charge " << get_atom_charge(i) << std::endl;
                atoms[i].set_ECP_electrons(elcount[j]);
                break;
            }
        }
    }
};

void WFN::operator=(const WFN &right)
{
    isBohr = right.isBohr;
    ncen = right.get_ncen();
    origin = right.get_origin();
    nex = right.get_nex();
    multi = right.get_multi();
    charge = right.get_charge();
    nfunc = right.get_nfunc();
    has_ECPs = right.get_has_ECPs();
    basis_set = right.get_basis_set_ptr();
    method = right.get_method();
    ECP_m = right.get_ECP_mode();
    comment = right.get_comment();
    basis_set_name = right.get_basis_set_name();
    path = right.get_path();
    multi = right.get_multi();
    total_energy = right.get_total_energy();
    virial_ratio = right.get_virial_ratio();
    UT_DensityMatrix = right.get_DensityMatrix();
    UT_SpinDensityMatrix = right.get_SpinDensityMatrix();
    for (int i = 0; i < nex; i++)
    {
        push_back_center(right.get_center(i));
        push_back_type(right.get_type(i));
        push_back_exponent(right.get_exponent(i));
    }
    for (int i = 0; i < right.get_nmo(); i++)
    {
        push_back_MO(i, right.get_MO_primitive_count(i), right.get_MO_energy(i));
        for (int j = 0; j < right.get_MO_primitive_count(i); j++)
            push_back_MO_coef(i, right.get_MO_coef(i, j));
    }
    for (int c = 0; c < right.cub.size(); c++)
    {
        cub.push_back(right.cub[c]);
    }
    for (int a = 0; a < right.atoms.size(); a++)
    {
        atoms.push_back(right.atoms[a]);
    }
    fill_pre();
    fill_Afac_pre();

};

int WFN::calculate_charge()
{
    int atomic_charges = 0;
    double mo_charges = 0;
    for (int a = 0; a < ncen; a++)
    {
        int nr = get_atom_charge(a);
        if (nr == 0)
        {
            std::cout << "ERROR: Atomtype misunderstanding!" << std::endl;
            return -1000;
        }
        atomic_charges += nr;
    }
    for (int mo = 0; mo < nmo; mo++)
    {
        mo_charges += get_MO_occ(mo);
    }
    return atomic_charges - (int)mo_charges;
};

int WFN::calculate_charge(std::ostream &file)
{
    int atomic_charges = 0;
    double mo_charges = 0;
    for (int a = 0; a < ncen; a++)
    {
        int nr = get_atom_charge(a);
        if (nr == 0)
        {
            file << "ERROR: Atomtype misunderstanding!" << std::endl;
            return -1000;
        }
        atomic_charges += nr;
    }
    for (int mo = 0; mo < nmo; mo++)
    {
        mo_charges += get_MO_occ(mo);
    }
    return atomic_charges - (int)mo_charges;
};

bool WFN::guess_multiplicity(std::ostream &file)
{
    if (get_nr_electrons() % 2 == 0)
    {
        file << "With " << get_nr_electrons() << " electrons your system appears to have multiplicity 1." << std::endl;
        assign_multi(1);
    }
    else if (get_nr_electrons() % 2 == 1)
    {
        file << "With " << get_nr_electrons() % 2 << " electron this seems to be open shell, assuming multiplicity 2." << std::endl;
        assign_multi(2);
    }
    else
    {
        file << "This is awkward... i dind't think of this case yet, contact Florian to implement it!" << std::endl;
        return false;
    }
    return true;
};

bool WFN::push_back_cube(const std::string &filepath, const bool &full, const bool &expert)
{
    cub.emplace_back(filepath, full, *this, std::cout, expert);
    return true;
};

void WFN::pop_back_cube()
{
    cub.pop_back();
}

const unsigned int WFN::get_atom_integer_mass(const unsigned int &atomnr) const
{
    if (get_atom_charge(atomnr) > 86)
    {
        std::cout << "Sorry, only implemented until Rn yet, ask Florian for increases!" << std::endl;
        return 0;
    }
    if (get_atom_charge(atomnr) <= 0)
    {
        std::cout << "sorry, something seems wrong with the atoms you requested!" << std::endl;
        return 0;
    }
    return constants::integer_masses[get_atom_charge(atomnr) - 1];
};

const double WFN::get_atom_real_mass(const int &atomnr) const
{
    if (get_atom_charge(atomnr) > 86)
    {
        std::cout << "Sorry, only implemented until Xe yet, ask Florian for increases!" << std::endl;
        return 0;
    }
    if (get_atom_charge(atomnr) <= 0)
    {
        std::cout << "sorry, something seems wrong with the atoms you requested!" << std::endl;
        return 0;
    }
    return constants::real_masses[get_atom_charge(atomnr) - 1];
}

const double& WFN::get_MO_occ(const int &nr) const
{
    return MOs[nr].get_occ();
};

const int& WFN::get_MO_op(const int &nr) const
{
    return MOs[nr].get_op();
};

void WFN::delete_unoccupied_MOs()
{
    for (int i = static_cast<int>(MOs.size()) - 1; i >= 0; i--)
    {
        if (get_MO_occ(i) == 0.0)
        {
            MOs.erase(MOs.begin() + i);
            nmo--;
        }
    }
};
void WFN::delete_Qs() {
    const int old_ncen = ncen;
    for (int i = static_cast<int>(atoms.size()) - 1; i >= 0; i--) {
        if (atoms[i].get_charge() == 119) {
            atoms.erase(atoms.begin() + i);
            ncen--;
            for (int j = 0; j < centers.size(); j++)
                if (centers[j] >= i) 
                    centers[j]--;
        }
    }
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
    std::string line;
    getline(fchk, line);
    std::string title = line;
    getline(fchk, line);
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
    const int el = read_fchk_integer(fchk, "Number of electrons", false);
    getline(fchk, line);
    const int ael = read_fchk_integer(line);
    getline(fchk, line);
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
    if (!read_fchk_integer_block(fchk, "Atomic numbers", atnbrs))
    {
        log << "Error reading atnbrs" << std::endl;
        return false;
    }
    ncen = static_cast<int>(atnbrs.size());
    atoms.resize(ncen);
    for (int i = 0; i < ncen; i++)
        atoms[i].set_label(constants::atnr2letter(atnbrs[i]));
    vec charges;
    if (!read_fchk_double_block(fchk, "Nuclear charges", charges))
    {
        log << "Error reading charges" << std::endl;
        return false;
    }
    for (int i = 0; i < charges.size(); i++)
        atoms[i].set_charge(static_cast<int>(charges[i]));
    vec coords;
    if (!read_fchk_double_block(fchk, "Current cartesian coordinates", coords))
    {
        log << "Error reading coordinates" << std::endl;
        return false;
    }
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
    if (!read_fchk_integer_block(fchk, "Shell types", shell_types, false))
    {
        log << "Error reading shell types" << std::endl;
        return false;
    }
    bool is_spherical = false;
    for (int i = 0; i < shell_types.size(); i++)
        if (shell_types[i] < -1)
            is_spherical = true;
    if (debug)
        log << "This fchk contains spherical harmonics, which will be transformed into cartesian functions!" << std::endl
            << "Loading basis set information..." << std::endl;
    ivec nr_prims_shell;
    if (!read_fchk_integer_block(fchk, "Number of primitives per shell", nr_prims_shell))
    {
        log << "Error reading primitives per shell" << std::endl;
        return false;
    }
    ivec shell2atom;
    if (!read_fchk_integer_block(fchk, "Shell to atom map", shell2atom))
    {
        log << "Error reading shell2atom" << std::endl;
        return false;
    }
    vec exp;
    if (!read_fchk_double_block(fchk, "Primitive exponents", exp))
    {
        log << "Error reading Primitive exponents" << std::endl;
        return false;
    }
    vec con;
    if (!read_fchk_double_block(fchk, "Contraction coefficients", con))
    {
        log << "Error reading Contraction coefficients" << std::endl;
        return false;
    }
    vec2 coef(2);
    vec2 MOocc(2), MOene(2);
    if (r_u_ro_switch == 0 || r_u_ro_switch == 2)
    { // Restricted or Restricted-Open-Shell
        if (!read_fchk_double_block(fchk, "Alpha Orbital Energies", MOene[0]))
        {
            log << "Error during reading of Alpha Energies" << std::endl;
            return false;
        }
        if (!read_fchk_double_block(fchk, "MO coefficients", coef[0]))
        {
            log << "Error during reading of Alpha MOs" << std::endl;
            return false;
        }
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
        if (!read_fchk_double_block(fchk, "Alpha Orbital Energies", MOene[0]))
        {
            log << "Error during reading of Alpha Energies" << std::endl;
            return false;
        }
        if (!read_fchk_double_block(fchk, "Beta Orbital Energies", MOene[1]))
        {
            log << "Error during reading of Beta Energies" << std::endl;
            return false;
        }
        if (!read_fchk_double_block(fchk, "Alpha MO coefficients", coef[0]))
        {
            log << "Error during reading of Alpha MOs" << std::endl;
            return false;
        }
        if (!read_fchk_double_block(fchk, "Beta MO coefficients", coef[1]))
        {
            log << "Error during reading of Beta MOs" << std::endl;
            return false;
        }
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

    int nprims = 0;
    int nshell = (int)shell_types.size();
    for (int i = 0; i < nshell; i++) {
        nprims += sht2nbas(abs(shell_types[i])) * nr_prims_shell[i];
    }
    nex = nprims;
    vec con_coefs;
    int exp_run = 0;
    for (int a = 0; a < shell_types.size(); a++)
    {
         double confac = 1.0;
         if (abs(shell_types[a]) == 0)
         {
             for (int i = 0; i < nr_prims_shell[a]; i++)
             {
                 confac = pow(8 * pow(exp[exp_run], 3) / constants::PI3, 0.25);
                 con_coefs.push_back(con[exp_run] * confac);
                 push_back_exponent(exp[exp_run]);
                 push_back_center(shell2atom[a]);
                 push_back_type(abs(shell_types[a])+1);
                 exp_run++;
             }
         }
         else if (abs(shell_types[a]) == 1)
         {
             for (int cart = 0; cart < 3; cart++) {
                 for (int i = 0; i < nr_prims_shell[a]; i++)
                 {
                     confac = pow(128 * pow(exp[exp_run + i], 5) / constants::PI3, 0.25);
                     con_coefs.push_back(con[exp_run + i] * confac);
                     push_back_exponent(exp[exp_run + i]);
                     push_back_center(shell2atom[a]);
                     push_back_type(2 + cart);
                 }
             }
             exp_run += nr_prims_shell[a];
         }
         else if (abs(shell_types[a]) == 2)
         {
             for (int cart = 0; cart < 6; cart++) {
                 for (int i = 0; i < nr_prims_shell[a]; i++)
                 {
                     confac = pow(2048 * pow(exp[exp_run + i], 7) / constants::PI3, 0.25);
                     con_coefs.push_back(con[exp_run + i] * confac);
                     push_back_exponent(exp[exp_run + i]);
                     push_back_center(shell2atom[a]);
                     push_back_type(5 + cart);
                 }
             }
             exp_run += nr_prims_shell[a];
         }
         else if (abs(shell_types[a]) == 3)
         {
             for (int cart = 0; cart < 10; cart++) {
                 for (int i = 0; i < nr_prims_shell[a]; i++)
                 {
                     confac = pow(32768 * pow(exp[exp_run + i], 9) / constants::PI3, 0.25);
                     con_coefs.push_back(con[exp_run + i] * confac);
                     push_back_exponent(exp[exp_run + i]);
                     push_back_center(shell2atom[a]);
                     push_back_type(11 + cart);
                 }
             }
             exp_run += nr_prims_shell[a];
         }
         else if (abs(shell_types[a]) == 4)
         {
             for (int cart = 0; cart < 15; cart++) {
                 for (int i = 0; i < nr_prims_shell[a]; i++)
                 {
                     confac = pow(524288 * pow(exp[exp_run + i], 11) / constants::PI3, 0.25);
                     con_coefs.push_back(con[exp_run + i] * confac);
                     push_back_exponent(exp[exp_run + i]);
                     push_back_center(shell2atom[a]);
                     push_back_type(21 + cart);
                 }
             }
             exp_run += nr_prims_shell[a];
         }
         else if (abs(shell_types[a]) == 5)
         {
             for (int cart = 0; cart < 21; cart++) {
                 for (int i = 0; i < nr_prims_shell[a]; i++)
                 {
                     confac = pow(8388608 * pow(exp[exp_run + i], 13) / constants::PI3, 0.25);
                     con_coefs.push_back(con[exp_run + i] * confac);
                     push_back_exponent(exp[exp_run + i]);
                     push_back_center(shell2atom[a]);
                     push_back_type(36 + cart);
                 }
             }
             exp_run += nr_prims_shell[a];
         }
         else if (abs(shell_types[a]) == 6)
         {
             //to-do: Have to calcualte confac for higher l
         }
         
    }
    vec2 p_pure_2_cart;
    vec2 d_pure_2_cart;
    vec2 f_pure_2_cart;
    vec2 g_pure_2_cart;
    err_checkf(generate_sph2cart_mat(p_pure_2_cart, d_pure_2_cart, f_pure_2_cart, g_pure_2_cart), "Error creating the conversion matrix", log);
    if (debug)
        log << "I read the basis of " << ncen << " atoms successfully" << std::endl;
    for (int i = 0; i < 2; i++) {
        if (MOocc[i].size() == 0)
            break;
        for (int j = 0; j < nbas; j++) {
            push_back_MO(i * nbas + j + 1, MOocc[i][j], MOene[i][j], 0);
            int cc_run = 0, coef_run = 0;
            for (int p = 0; p < nr_prims_shell.size(); p++)
            {
                int sw = abs(shell_types[p]);
                switch (sw)
                {
                case 0:
                {
                    for (int s = 0; s < nr_prims_shell[p]; s++)
                    {
                        push_back_MO_coef(j, coef[i][j * nbas + coef_run] * con_coefs[cc_run + s]);
                    }
                    coef_run++;
                    cc_run += nr_prims_shell[p];
                    break;
                }
                case 1:
                {
                    for (int cart = 0; cart < 3; cart++)
                    {
                        for (int s = 0; s < nr_prims_shell[p]; s++)
                        {
                            push_back_MO_coef(j, coef[i][j * nbas + coef_run + cart] * con_coefs[cc_run + s]);
                        }
                    }
                    coef_run += 3;
                    cc_run += 3 * nr_prims_shell[p];
                    break;
                }
                case 2:
                {
                    double temp_coef = 0;
                    for (int cart = 0; cart < 6; cart++)
                    {
                        temp_coef = 0;
                        for (int spher = 0; spher < 5; spher++)
                        {
                            temp_coef += d_pure_2_cart[cart][spher] * coef[i][j * nbas + coef_run + spher];
                        }
                        for (int s = 0; s < nr_prims_shell[p]; s++)
                        {
                            push_back_MO_coef(j, temp_coef * con_coefs[cc_run + s]);
                        }
                    }
                    coef_run += 5;
                    cc_run += 6 * nr_prims_shell[p];
                    break;
                }
                case 3:
                {
                    double temp_coef = 0;
                    for (int cart = 0; cart < 10; cart++)
                    {
                        temp_coef = 0;
                        for (int spher = 0; spher < 7; spher++)
                        {
                            temp_coef += f_pure_2_cart[cart][spher] * coef[i][j * nbas + coef_run + spher];
                        }
                        for (int s = 0; s < nr_prims_shell[p]; s++)
                        {
                            push_back_MO_coef(j, temp_coef * con_coefs[cc_run + s]);
                        }
                    }
                    coef_run += 7;
                    cc_run += 10 * nr_prims_shell[p];
                    break;
                }
                case 4:
                {
                    double temp_coef = 0;
                    for (int cart = 0; cart < 15; cart++)
                    {
                        temp_coef = 0;
                        for (int spher = 0; spher < 9; spher++)
                        {
                            temp_coef += g_pure_2_cart[cart][spher] * coef[i][j * nbas + coef_run + spher];
                        }
                        for (int s = 0; s < nr_prims_shell[p]; s++)
                        {
                            push_back_MO_coef(j, temp_coef * con_coefs[cc_run + s]);
                        }
                    }
                    coef_run += 9;
                    cc_run += 15 * nr_prims_shell[p];
                    break;
                }
                default:
                {
                    if (debug)
                        log << "This is not supposed to happen!" << std::endl;
                    err_not_impl_f("Types higher than g type in fchk", log);
                    break;
                }
                }
            }


        }
    }
    constants::exp_cutoff = std::log(constants::density_accuracy / get_maximum_MO_coefficient());
    return true;
};

const double WFN::compute_dens(
    const double &Pos1,
    const double &Pos2,
    const double &Pos3,
    vec2 &d,
    vec &phi) const
{
    if (d_f_switch)
    {
        err_checkf(d.size() >= 5, "d is too small!", std::cout);
        err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
        return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
    }
    else
    {
        err_checkf(d.size() >= 16, "d is too small!", std::cout);
        err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
        return compute_dens_cartesian(Pos1, Pos2, Pos3, d, phi);
    }
};

const double WFN::compute_dens(
    const double &Pos1,
    const double &Pos2,
    const double &Pos3) const
{
    vec2 d;
    vec phi(nmo, 0.0);

    if (d_f_switch)
    {
        d.resize(5);
        for (int i = 0; i < 5; i++)
            d[i].resize(ncen, 0.0);
        err_not_impl_f("Nah.. not yet implemented correctly", std::cout);
        return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
    }
    else
    {
        d.resize(16);
        for (int i = 0; i < 16; i++)
            d[i].resize(ncen, 0.0);
        return compute_dens_cartesian(Pos1, Pos2, Pos3, d, phi);
    }
};

const double WFN::compute_spin_dens(
    const double &Pos1,
    const double &Pos2,
    const double &Pos3,
    vec2 &d,
    vec &phi) const
{
    if (d_f_switch)
    {
        err_checkf(d.size() >= 5, "d is too small!", std::cout);
        err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
        err_not_impl_f("Nah.. not yet implemented correctly", std::cout);
        return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
    }
    else
    {
        err_checkf(d.size() >= 4, "d is too small!", std::cout);
        err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
        return compute_spin_dens_cartesian(Pos1, Pos2, Pos3, d, phi);
    }
};

const double WFN::compute_spin_dens(
    const double &Pos1,
    const double &Pos2,
    const double &Pos3) const
{
    vec2 d;
    vec phi(nmo, 0.0);

    if (d_f_switch)
    {
        d.resize(5);
        for (int i = 0; i < 5; i++)
            d[i].resize(ncen, 0.0);
        err_not_impl_f("Nah.. not yet implemented correctly", std::cout);
        return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
    }
    else
    {
        d.resize(16);
        for (int i = 0; i < 16; i++)
            d[i].resize(ncen, 0.0);
        return compute_spin_dens_cartesian(Pos1, Pos2, Pos3, d, phi);
    }
};

const double WFN::compute_dens_cartesian(
    const double &Pos1,
    const double &Pos2,
    const double &Pos3,
    vec2 &d,
    vec &phi) const
{
    std::fill(phi.begin(), phi.end(), 0.0);
    double Rho = 0.0;
    int iat, j, k;
    int l[3];
    double ex;
    int mo;

    // precalculate some distances and powers of distances for faster computation
    for (iat = 0; iat < ncen; iat++)
    {
        d[0][iat] = Pos1 - atoms[iat].get_coordinate(0);
        d[1][iat] = Pos2 - atoms[iat].get_coordinate(1);
        d[2][iat] = Pos3 - atoms[iat].get_coordinate(2);
        d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
        d[4][iat] = d[0][iat] * d[0][iat];
        d[5][iat] = d[1][iat] * d[1][iat];
        d[6][iat] = d[2][iat] * d[2][iat];
        d[7][iat] = d[0][iat] * d[4][iat];
        d[8][iat] = d[1][iat] * d[5][iat];
        d[9][iat] = d[2][iat] * d[6][iat];
        d[10][iat] = d[0][iat] * d[7][iat];
        d[11][iat] = d[1][iat] * d[8][iat];
        d[12][iat] = d[2][iat] * d[9][iat];
        d[13][iat] = d[0][iat] * d[10][iat];
        d[14][iat] = d[1][iat] * d[11][iat];
        d[15][iat] = d[2][iat] * d[12][iat];
    }

    for (j = 0; j < nex; j++)
    {
        iat = centers[j] - 1;
        constants::type2vector(types[j], l);
        ex = -exponents[j] * d[3][iat];
        if (ex < constants::exp_cutoff)
        { // corresponds to cutoff of maximum density contribution of 1E-5
            continue;
        }
        ex = exp(ex);
        for (k = 0; k < 3; k++)
        {
            if (l[k] == 0)
                continue;
            else if (l[k] == 1)
                ex *= d[k][iat];
            else if (l[k] == 2)
                ex *= d[k + 4][iat];
            else if (l[k] == 3)
                ex *= d[k + 7][iat];
            else if (l[k] == 4)
                ex *= d[k + 10][iat];
            else if (l[k] == 5)
                ex *= d[k + 13][iat];
        }
        double *run = phi.data();
        const MO *run2 = MOs.data();
        for (mo = 0; mo < nmo; mo++)
        {
            *run += (*run2).get_coefficient_f(j) * ex; // build MO values at this point
            //*run += *(coef_ptrs[mo]) * ex;
            run2++, run++;
            //(coef_ptrs[mo])++, run++;
        }
    }

    double* run = phi.data();
    const MO* run2 = MOs.data();
    for (mo = 0; mo < nmo; mo++)
    {
        Rho += (*run2).get_occ() * pow(*run, 2);
        run++, run2++;
    }

    return Rho;
}

const double WFN::compute_spin_dens_cartesian(
    const double &Pos1,
    const double &Pos2,
    const double &Pos3,
    vec2 &d,
    vec &phi) const
{
    std::fill(phi.begin(), phi.end(), 0.0);
    double alpha = 0.0, beta = 0.0;
    int iat, j, k;
    int l[3];
    double ex;
    int mo;

    for (iat = 0; iat < ncen; iat++)
    {
        d[0][iat] = Pos1 - atoms[iat].get_coordinate(0);
        d[1][iat] = Pos2 - atoms[iat].get_coordinate(1);
        d[2][iat] = Pos3 - atoms[iat].get_coordinate(2);
        d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
        d[4][iat] = d[0][iat] * d[0][iat];
        d[5][iat] = d[1][iat] * d[1][iat];
        d[6][iat] = d[2][iat] * d[2][iat];
        d[7][iat] = d[0][iat] * d[4][iat];
        d[8][iat] = d[1][iat] * d[5][iat];
        d[9][iat] = d[2][iat] * d[6][iat];
        d[10][iat] = d[0][iat] * d[7][iat];
        d[11][iat] = d[1][iat] * d[8][iat];
        d[12][iat] = d[2][iat] * d[9][iat];
        d[13][iat] = d[0][iat] * d[10][iat];
        d[14][iat] = d[1][iat] * d[11][iat];
        d[15][iat] = d[2][iat] * d[12][iat];
    }

    for (j = 0; j < nex; j++)
    {
        iat = centers[j] - 1;
        constants::type2vector(types[j], l);
        ex = -exponents[j] * d[3][iat];
        if (ex < constants::exp_cutoff)
        { // corresponds to cutoff of ex ~< 1E-20
            continue;
        }
        ex = exp(ex);
        for (k = 0; k < 3; k++)
        {
            if (l[k] == 0)
                continue;
            else if (l[k] == 1)
                ex *= d[k][iat];
            else if (l[k] == 2)
                ex *= d[k + 4][iat];
            else if (l[k] == 3)
                ex *= d[k + 7][iat];
            else if (l[k] == 4)
                ex *= d[k + 10][iat];
            else if (l[k] == 5)
                ex *= d[k + 13][iat];
        }
        double *run = phi.data();
        const MO *run2 = MOs.data();
        for (mo = 0; mo < nmo; mo++)
        {
            *run += (*run2).get_coefficient_f(j) * ex; // build MO values at this point
            run2++, run++;
        }
    }

    double *run = phi.data();
    const MO *run2 = MOs.data();
    for (mo = 0; mo < nmo; mo++)
    {
        if ((*run2).get_op())
            beta += (*run2).get_occ() * pow(*run, 2);
        else
            alpha += (*run2).get_occ() * pow(*run, 2);
        run++, run2++;
    }

    return alpha - beta;
}

const double WFN::compute_MO_spherical(
    const double &Pos1,
    const double &Pos2,
    const double &Pos3,
    const int &MO) const
{
    err_not_impl_f("This one is not tested an will most likely not work, therefore aborting!", std::cout);
    return 0.0;
    (void)Pos1;
    (void)Pos2;
    (void)Pos3;
    (void)MO;
    //err_checkf(d_f_switch, "Only works for spheriacl wavefunctions!", std::cout);
    //int iat;
    //int l = 0;
    //// ex will carry information about radial function
    //double ex;
    //vec2 d(5);
    //for (int i = 0; i < 5; i++)
    //    d[i].resize(ncen);
    //double phi(0.0);

    //for (iat = 0; iat < ncen; iat++)
    //{
    //    d[0][iat] = Pos1 - atoms[iat].x;
    //    d[1][iat] = Pos2 - atoms[iat].y;
    //    d[2][iat] = Pos3 - atoms[iat].z;
    //    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
    //    d[4][iat] = sqrt(d[3][iat]);
    //}
    ///*Here d[0] = x
    //             d[1] = y
    //             d[2] = z
    //             d[3] = r^2
    //             d[4] = r
    //             */

    //for (int j = 0; j < nex; j++)
    //{
    //    iat = centers[j] - 1;
    //    ex = -exponents[j] * d[3][iat];
    //    if (ex < -46.0517)
    //    { // corresponds to cutoff of ex ~< 1E-20
    //        continue;
    //    }
    //    // apply radial function:
    //    ex = exp(ex);
    //    int lam = 0, m = 0;
    //    if (l > 1)
    //    {
    //        if (l <= 4)
    //        {
    //            ex *= d[4][iat];
    //            lam = 1;
    //            if (l == 2)
    //                m = 0;
    //            else if (l == 3)
    //                m = 1;
    //            else if (l == 4)
    //                m = -1;
    //        }
    //        else if (l <= 9)
    //        {
    //            ex *= d[3][iat];
    //            lam = 2;
    //            if (l == 5)
    //                m = 0;
    //            else if (l == 6)
    //                m = 1;
    //            else if (l == 7)
    //                m = -1;
    //            else if (l == 8)
    //                m = 2;
    //            else if (l == 9)
    //                m = -2;
    //        }
    //        else if (l <= 16)
    //        {
    //            ex *= d[3][iat] * d[4][iat];
    //            lam = 3;
    //            if (l == 10)
    //                m = 0;
    //            else if (l == 11)
    //                m = 1;
    //            else if (l == 12)
    //                m = -1;
    //            else if (l == 13)
    //                m = 2;
    //            else if (l == 14)
    //                m = -2;
    //            else if (l == 15)
    //                m = 3;
    //            else if (l == 16)
    //                m = -3;
    //        }
    //        else if (l <= 25)
    //        {
    //            ex *= d[3][iat] * d[3][iat];
    //            lam = 4;
    //            if (l == 17)
    //                m = 0;
    //            else if (l == 18)
    //                m = 1;
    //            else if (l == 19)
    //                m = -1;
    //            else if (l == 20)
    //                m = 2;
    //            else if (l == 21)
    //                m = -2;
    //            else if (l == 22)
    //                m = 3;
    //            else if (l == 23)
    //                m = -3;
    //            else if (l == 24)
    //                m = 4;
    //            else if (l == 25)
    //                m = -4;
    //        }
    //        else if (l <= 36)
    //        {
    //            ex *= pow(d[4][iat], 5);
    //            lam = 5;
    //            if (l == 26)
    //                m = 0;
    //            else if (l == 27)
    //                m = 1;
    //            else if (l == 28)
    //                m = -1;
    //            else if (l == 29)
    //                m = 2;
    //            else if (l == 30)
    //                m = -2;
    //            else if (l == 31)
    //                m = 3;
    //            else if (l == 32)
    //                m = -3;
    //            else if (l == 33)
    //                m = 4;
    //            else if (l == 34)
    //                m = -4;
    //            else if (l == 35)
    //                m = 5;
    //            else if (l == 36)
    //                m = -5;
    //        }
    //    }
    //    // calc spherical harmonic
    //    double d_t[]{d[0][iat], d[1][iat], d[2][iat], d[3][iat], d[4][iat]};
    //    double SH = spherical_harmonic(lam, m, d_t);
    //    SH *= ex;                                 // multiply radial part with spherical harmonic
    //    phi += MOs[MO].get_coefficient_f(j) * SH; // build MO values at this point
    //}
    //shrink_vector<vec>(d);

    //return phi;
}

const double WFN::compute_dens_spherical(
    const double &Pos1,
    const double &Pos2,
    const double &Pos3,
    vec2 &d,
    vec &phi) const
{
    err_not_impl_f("This one is not tested an will most likely not work, therefore aborting!", std::cout);
    return 0.0;
    (void)Pos1;
    (void)Pos2;
    (void)Pos3;
    (void)d;
    (void)phi;
    //err_checkf(d_f_switch, "Only works for spheriacl wavefunctions!", std::cout);
    //std::fill(phi.begin(), phi.end(), 0.0);
    //double Rho = 0.0;
    //int iat;
    //int l;
    //// ex will carry information about radial function
    //double ex;
    //int mo;

    //for (iat = 0; iat < ncen; iat++)
    //{
    //    d[0][iat] = Pos1 - atoms[iat].x;
    //    d[1][iat] = Pos2 - atoms[iat].y;
    //    d[2][iat] = Pos3 - atoms[iat].z;
    //    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
    //    d[4][iat] = sqrt(d[3][iat]);
    //}
    ///*Here d[0] = x
    //             d[1] = y
    //             d[2] = z
    //             d[3] = r^2
    //             d[4] = r
    //             */
    //for (int j = 0; j < nex; j++)
    //{
    //    iat = centers[j] - 1;
    //    ex = -exponents[j] * d[3][iat];
    //    if (ex < -46.0517)
    //    { // corresponds to cutoff of ex ~< 1E-20
    //        continue;
    //    }
    //    ex = exp(ex);
    //    // apply radial function:
    //    l = types[j];
    //    int lam = 0, m = 0;
    //    if (l > 1)
    //    {
    //        if (l <= 4)
    //        {
    //            ex *= d[4][iat];
    //            lam = 1;
    //            if (l == 2)
    //                m = 0;
    //            else if (l == 3)
    //                m = 1;
    //            else if (l == 4)
    //                m = -1;
    //        }
    //        else if (l <= 9)
    //        {
    //            ex *= d[3][iat];
    //            lam = 2;
    //            if (l == 5)
    //                m = 0;
    //            else if (l == 6)
    //                m = 1;
    //            else if (l == 7)
    //                m = -1;
    //            else if (l == 8)
    //                m = 2;
    //            else if (l == 9)
    //                m = -2;
    //        }
    //        else if (l <= 16)
    //        {
    //            ex *= d[3][iat] * d[4][iat];
    //            lam = 3;
    //            if (l == 10)
    //                m = 0;
    //            else if (l == 11)
    //                m = 1;
    //            else if (l == 12)
    //                m = -1;
    //            else if (l == 13)
    //                m = 2;
    //            else if (l == 14)
    //                m = -2;
    //            else if (l == 15)
    //                m = 3;
    //            else if (l == 16)
    //                m = -3;
    //        }
    //        else if (l <= 25)
    //        {
    //            ex *= d[3][iat] * d[3][iat];
    //            lam = 4;
    //            if (l == 17)
    //                m = 0;
    //            else if (l == 18)
    //                m = 1;
    //            else if (l == 19)
    //                m = -1;
    //            else if (l == 20)
    //                m = 2;
    //            else if (l == 21)
    //                m = -2;
    //            else if (l == 22)
    //                m = 3;
    //            else if (l == 23)
    //                m = -3;
    //            else if (l == 24)
    //                m = 4;
    //            else if (l == 25)
    //                m = -4;
    //        }
    //        else if (l <= 36)
    //        {
    //            ex *= pow(d[4][iat], 5);
    //            lam = 5;
    //            if (l == 26)
    //                m = 0;
    //            else if (l == 27)
    //                m = 1;
    //            else if (l == 28)
    //                m = -1;
    //            else if (l == 29)
    //                m = 2;
    //            else if (l == 30)
    //                m = -2;
    //            else if (l == 31)
    //                m = 3;
    //            else if (l == 32)
    //                m = -3;
    //            else if (l == 33)
    //                m = 4;
    //            else if (l == 34)
    //                m = -4;
    //            else if (l == 35)
    //                m = 5;
    //            else if (l == 36)
    //                m = -5;
    //        }
    //    }
    //    double d_t[]{d[0][iat], d[1][iat], d[2][iat], d[3][iat], d[4][iat]};
    //    double SH = spherical_harmonic(lam, m, d_t);
    //    SH *= ex; // multiply radial part with spherical harmonic
    //    auto run = phi.data();
    //    auto run2 = MOs.data();
    //    for (mo = 0; mo < phi.size(); mo++)
    //    {
    //        *run += (*run2).get_coefficient_f(j) * SH; // build MO values at this point
    //        run++, run2++;
    //    }
    //}

    //auto run = phi.data();
    //auto run2 = MOs.data();
    //for (mo = 0; mo < phi.size(); mo++)
    //{
    //    Rho += (*run2).get_occ() * pow(*run, 2);
    //    run++, run2++;
    //}

    //return Rho;
}

void WFN::pop_back_MO()
{
    MOs.pop_back();
    nmo--;
}

const void WFN::computeValues(
    const std::array<double,3>& PosGrid, // [3] vector with current position on te grid
    double &Rho,           // Value of Electron Density
    double &normGrad,      // Gradiant Vector
    double *Hess,          // Hessian Matrix, later used to determine lambda2
    double &Elf,           // Value of the ELF
    double &Eli,           // Value of the ELI
    double &Lap            // Value for the Laplacian
) const
{
    const int _nmo = get_nmo(false);
    vec phi(10 * _nmo, 0.0);
    double *phi_temp;
    double chi[10]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double d[4]{0, 0, 0, 0};
    int iat = 0;
    int l[3]{0, 0, 0};
    double ex = 0;
    double xl[3][3]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    double Grad[3]{0, 0, 0};
    double tau = 0;

    Rho = 0;
    Elf = 0;
    Hess[0] = 0;
    Hess[1] = 0;
    Hess[2] = 0;
    Hess[8] = 0;
    Hess[4] = 0;
    Hess[5] = 0;

    for (int j = 0; j < nex; j++)
    {
        iat = get_center(j) - 1;

        constants::type2vector(get_type(j), l);
        d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
        d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
        d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
        d[3] = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        double temp = -get_exponent(j) * (d[3]);
        if (temp < constants::exp_cutoff)
            continue;
        ex = exp(temp);
        for (int k = 0; k < 3; k++)
        {
            if (l[k] == 0)
            {
                xl[k][0] = 1.0;
                xl[k][1] = 0.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 1)
            {
                xl[k][0] = d[k];
                xl[k][1] = 1.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 2)
            {
                xl[k][0] = d[k] * d[k];
                xl[k][1] = 2 * d[k];
                xl[k][2] = 2;
            }
            else if (l[k] == 3)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d[k];
                xl[k][1] = 3 * d2;
                xl[k][2] = 6 * d[k];
            }
            else if (l[k] == 4)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d2;
                xl[k][1] = 4 * d2 * d[k];
                xl[k][2] = 12 * d2;
            }
            else
            {
                return;
            }
        }
        double ex2 = 2 * get_exponent(j);
        chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
        chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
        chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
        chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
        double temp_ex = pow(ex2, 2);
        chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
        chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
        chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;
        chi[7] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[2][0] * ex;
        chi[8] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[1][0] * ex;
        chi[9] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * ex;
        for (int mo = 0; mo < _nmo; mo++)
        {
            phi_temp = &phi[mo * 10];
            for (int i = 0; i < 10; i++)
                //                if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
                phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i]; // build MO values at this point
        }
    }

    for (int mo = 0; mo < _nmo; mo++)
    {
        const double occ = get_MO_occ(mo);
        const double docc = 2 * occ;
        if (occ != 0)
        {
            phi_temp = &phi[mo * 10];
            Rho += occ * pow(*phi_temp, 2);
            Grad[0] += docc * *phi_temp * phi_temp[1];
            Grad[1] += docc * *phi_temp * phi_temp[2];
            Grad[2] += docc * *phi_temp * phi_temp[3];
            Hess[0] += docc * (*phi_temp * phi_temp[4] + pow(phi_temp[1], 2));
            Hess[4] += docc * (*phi_temp * phi_temp[5] + pow(phi_temp[2], 2));
            Hess[8] += docc * (*phi_temp * phi_temp[6] + pow(phi_temp[3], 2));
            Hess[1] += docc * (*phi_temp * phi_temp[7] + phi_temp[1] * phi_temp[2]);
            Hess[2] += docc * (*phi_temp * phi_temp[8] + phi_temp[1] * phi_temp[3]);
            Hess[5] += docc * (*phi_temp * phi_temp[9] + phi_temp[2] * phi_temp[3]);
            tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
        }
    }

    Hess[3] = Hess[1];
    Hess[6] = Hess[2];
    Hess[7] = Hess[5];
    if (Rho > 0)
    {
        normGrad = constants::alpha_coef * sqrt(Grad[0] * Grad[0] + Grad[1] * Grad[1] + Grad[2] * Grad[2]) / pow(Rho, constants::c_43);
        Elf = 1 / (1 + pow(constants::ctelf * pow(Rho, constants::c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
        Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
    }
    Lap = Hess[0] + Hess[4] + Hess[8];
};

const void WFN::computeELIELF(
    const std::array<double, 3>& PosGrid, // [3] vector with current position on te grid
    double &Elf,           // Value of the ELF
    double &Eli            // Value of the ELI
) const
{
    const int _nmo = get_nmo(false);
    vec phi(4 * _nmo, 0.0);
    double *phi_temp;
    double chi[4]{0, 0, 0, 0};
    double d[3]{0, 0, 0};
    int iat = 0;
    int l[3]{0, 0, 0};
    double ex = 0;
    double xl[3][3]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int j = 0; j < nex; j++)
    {
        iat = get_center(j) - 1;

        constants::type2vector(get_type(j), l);
        d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
        d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
        d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
        double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (temp < constants::exp_cutoff)
            continue;
        ex = exp(temp);
        for (int k = 0; k < 3; k++)
        {
            if (l[k] == 0)
            {
                xl[k][0] = 1.0;
                xl[k][1] = 0.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 1)
            {
                xl[k][0] = d[k];
                xl[k][1] = 1.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 2)
            {
                xl[k][0] = d[k] * d[k];
                xl[k][1] = 2 * d[k];
                xl[k][2] = 2;
            }
            else if (l[k] == 3)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d[k];
                xl[k][1] = 3 * d2;
                xl[k][2] = 6 * d[k];
            }
            else if (l[k] == 4)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d2;
                xl[k][1] = 4 * d2 * d[k];
                xl[k][2] = 12 * d2;
            }
            else
            {
                return;
            }
        }
        double ex2 = 2 * get_exponent(j);
        chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
        chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
        chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
        chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
        for (int mo = 0; mo < _nmo; mo++)
        {
            phi_temp = &phi[mo * 4];
            for (int i = 0; i < 4; i++)
                //                if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
                phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i]; // build MO values at this point
        }
    }

    double Grad[3]{0, 0, 0};
    double tau = 0;
    double Rho = 0;

    for (int mo = 0; mo < _nmo; mo++)
    {
        const double occ = get_MO_occ(mo);
        const double docc = 2 * occ;
        if (occ != 0)
        {
            phi_temp = &phi[mo * 4];
            Rho += occ * pow(*phi_temp, 2);
            Grad[0] += docc * *phi_temp * phi_temp[1];
            Grad[1] += docc * *phi_temp * phi_temp[2];
            Grad[2] += docc * *phi_temp * phi_temp[3];
            tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
        }
    }
    if (Rho > 0)
    {
        Elf = 1 / (1 + pow(constants::ctelf * pow(Rho, constants::c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
        Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
    }
};

const void WFN::computeELI(
    const std::array<double, 3>& PosGrid, // [3] vector with current position on te grid
    double &Eli            // Value of the ELI
) const
{
    const int _nmo = get_nmo(false);
    vec phi(4 * _nmo, 0.0);
    double *phi_temp;
    double chi[4]{0, 0, 0, 0};
    double d[3]{0, 0, 0};
    int iat = 0;
    int l[3]{0, 0, 0};
    double ex = 0;
    double xl[3][3]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int j = 0; j < nex; j++)
    {
        iat = get_center(j) - 1;

        constants::type2vector(get_type(j), l);
        d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
        d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
        d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
        double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (temp < constants::exp_cutoff)
            continue;
        ex = exp(temp);
        for (int k = 0; k < 3; k++)
        {
            if (l[k] == 0)
            {
                xl[k][0] = 1.0;
                xl[k][1] = 0.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 1)
            {
                xl[k][0] = d[k];
                xl[k][1] = 1.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 2)
            {
                xl[k][0] = d[k] * d[k];
                xl[k][1] = 2 * d[k];
                xl[k][2] = 2;
            }
            else if (l[k] == 3)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d[k];
                xl[k][1] = 3 * d2;
                xl[k][2] = 6 * d[k];
            }
            else if (l[k] == 4)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d2;
                xl[k][1] = 4 * d2 * d[k];
                xl[k][2] = 12 * d2;
            }
            else
            {
                return;
            }
        }
        double ex2 = 2 * get_exponent(j);
        chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
        chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
        chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
        chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
        for (int mo = 0; mo < _nmo; mo++)
        {
            phi_temp = &phi[mo * 4];
            for (int i = 0; i < 4; i++)
                //                if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
                phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i]; // build MO values at this point
        }
    }

    double Grad[3]{0, 0, 0};
    double tau = 0;
    double Rho = 0;

    for (int mo = 0; mo < _nmo; mo++)
    {
        const double occ = get_MO_occ(mo);
        const double docc = 2 * occ;
        if (occ != 0)
        {
            phi_temp = &phi[mo * 4];
            Rho += occ * pow(*phi_temp, 2);
            Grad[0] += docc * *phi_temp * phi_temp[1];
            Grad[1] += docc * *phi_temp * phi_temp[2];
            Grad[2] += docc * *phi_temp * phi_temp[3];
            tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
        }
    }
    Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
};

const void WFN::computeELF(
    const std::array<double, 3>& PosGrid, // [3] vector with current position on te grid
    double &Elf            // Value of the ELF
) const
{
    const int _nmo = get_nmo(false);
    vec phi(4 * _nmo, 0.0);
    double *phi_temp;
    double chi[4]{0, 0, 0, 0};
    double d[3]{0, 0, 0};
    int iat = 0;
    int l[3]{0, 0, 0};
    double ex = 0;
    double xl[3][3]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int j = 0; j < nex; j++)
    {
        iat = get_center(j) - 1;

        constants::type2vector(get_type(j), l);
        d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
        d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
        d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
        double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (temp < constants::exp_cutoff)
            continue;
        ex = exp(temp);
        for (int k = 0; k < 3; k++)
        {
            if (l[k] == 0)
            {
                xl[k][0] = 1.0;
                xl[k][1] = 0.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 1)
            {
                xl[k][0] = d[k];
                xl[k][1] = 1.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 2)
            {
                xl[k][0] = d[k] * d[k];
                xl[k][1] = 2 * d[k];
                xl[k][2] = 2;
            }
            else if (l[k] == 3)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d[k];
                xl[k][1] = 3 * d2;
                xl[k][2] = 6 * d[k];
            }
            else if (l[k] == 4)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d2;
                xl[k][1] = 4 * d2 * d[k];
                xl[k][2] = 12 * d2;
            }
            else
            {
                return;
            }
        }
        double ex2 = 2 * get_exponent(j);
        chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
        chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
        chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
        chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
        for (int mo = 0; mo < _nmo; mo++)
        {
            phi_temp = &phi[mo * 4];
            for (int i = 0; i < 4; i++)
                //                if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
                phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i]; // build MO values at this point
        }
    }

    double Grad[3]{0, 0, 0};
    double tau = 0;
    double Rho = 0;

    for (int mo = 0; mo < _nmo; mo++)
    {
        const double occ = get_MO_occ(mo);
        const double docc = 2 * occ;
        if (occ != 0)
        {
            phi_temp = &phi[mo * 4];
            Rho += occ * pow(*phi_temp, 2);
            Grad[0] += docc * *phi_temp * phi_temp[1];
            Grad[1] += docc * *phi_temp * phi_temp[2];
            Grad[2] += docc * *phi_temp * phi_temp[3];
            tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
        }
    }
    Elf = 1 / (1 + pow(constants::ctelf * pow(Rho, constants::c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
};

const void WFN::computeLapELIELF(
    const std::array<double, 3>& PosGrid, // [3] vector with current position on te grid
    double &Elf,           // Value of the ELF
    double &Eli,           // Value of the ELI
    double &Lap            // Value for the Laplacian
) const
{
    const int _nmo = get_nmo(false);
    vec phi(7 * _nmo, 0.0);
    double *phi_temp;
    double chi[7]{0, 0, 0, 0, 0, 0, 0};
    double d[3]{0, 0, 0};
    int iat = 0;
    int l[3]{0, 0, 0};
    double ex = 0;
    double xl[3][3]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int j = 0; j < nex; j++)
    {
        iat = get_center(j) - 1;

        constants::type2vector(get_type(j), l);
        d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
        d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
        d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
        double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (temp < constants::exp_cutoff)
            continue;
        ex = exp(temp);
        for (int k = 0; k < 3; k++)
        {
            if (l[k] == 0)
            {
                xl[k][0] = 1.0;
                xl[k][1] = 0.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 1)
            {
                xl[k][0] = d[k];
                xl[k][1] = 1.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 2)
            {
                xl[k][0] = d[k] * d[k];
                xl[k][1] = 2 * d[k];
                xl[k][2] = 2;
            }
            else if (l[k] == 3)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d[k];
                xl[k][1] = 3 * d2;
                xl[k][2] = 6 * d[k];
            }
            else if (l[k] == 4)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d2;
                xl[k][1] = 4 * d2 * d[k];
                xl[k][2] = 12 * d2;
            }
            else
            {
                return;
            }
        }
        double ex2 = 2 * get_exponent(j);
        chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
        chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
        chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
        chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
        double temp_ex = pow(ex2, 2);
        chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
        chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
        chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;
        for (int mo = 0; mo < _nmo; mo++)
        {
            phi_temp = &phi[mo * 7];
            for (int i = 0; i < 7; i++)
                //                if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
                phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i]; // build MO values at this point
        }
    }

    double Grad[3]{0, 0, 0};
    double Hess[3]{0, 0, 0};
    double tau = 0;
    double Rho = 0;

    Elf = 0;

    for (int mo = 0; mo < _nmo; mo++)
    {
        const double occ = get_MO_occ(mo);
        const double docc = 2 * occ;
        if (occ != 0)
        {
            phi_temp = &phi[mo * 7];
            Rho += occ * pow(*phi_temp, 2);
            Grad[0] += docc * *phi_temp * phi_temp[1];
            Grad[1] += docc * *phi_temp * phi_temp[2];
            Grad[2] += docc * *phi_temp * phi_temp[3];
            Hess[0] += docc * (*phi_temp * phi_temp[4] + pow(phi_temp[1], 2));
            Hess[1] += docc * (*phi_temp * phi_temp[5] + pow(phi_temp[2], 2));
            Hess[2] += docc * (*phi_temp * phi_temp[6] + pow(phi_temp[3], 2));
            tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
        }
    }
    Elf = 1 / (1 + pow(constants::ctelf * pow(Rho, constants::c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
    Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
    Lap = Hess[0] + Hess[1] + Hess[2];
};

const void WFN::computeLapELI(
    const std::array<double, 3>& PosGrid, // [3] vector with current position on te grid
    double &Eli,           // Value of the ELI
    double &Lap            // Value for the Laplacian
) const
{
    const int _nmo = get_nmo(false);
    vec phi(7 * _nmo, 0.0);
    double *phi_temp;
    double chi[7]{0, 0, 0, 0, 0, 0, 0};
    double d[3]{0, 0, 0};
    int iat = 0;
    int l[3]{0, 0, 0};
    double ex = 0;
    double xl[3][3]{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int j = 0; j < nex; j++)
    {
        iat = get_center(j) - 1;

        constants::type2vector(get_type(j), l);
        d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
        d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
        d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
        double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (temp < constants::exp_cutoff)
            continue;
        ex = exp(temp);
        for (int k = 0; k < 3; k++)
        {
            if (l[k] == 0)
            {
                xl[k][0] = 1.0;
                xl[k][1] = 0.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 1)
            {
                xl[k][0] = d[k];
                xl[k][1] = 1.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 2)
            {
                xl[k][0] = d[k] * d[k];
                xl[k][1] = 2 * d[k];
                xl[k][2] = 2;
            }
            else if (l[k] == 3)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d[k];
                xl[k][1] = 3 * d2;
                xl[k][2] = 6 * d[k];
            }
            else if (l[k] == 4)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d2;
                xl[k][1] = 4 * d2 * d[k];
                xl[k][2] = 12 * d2;
            }
            else
            {
                return;
            }
        }
        double ex2 = 2 * get_exponent(j);
        chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
        chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
        chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
        chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
        double temp_ex = pow(ex2, 2);
        chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
        chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
        chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;
        for (int mo = 0; mo < _nmo; mo++)
        {
            phi_temp = &phi[mo * 7];
            for (int i = 0; i < 7; i++)
                //                if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
                phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i]; // build MO values at this point
        }
    }

    double Grad[3]{0, 0, 0};
    double Hess[3]{0, 0, 0};
    double tau = 0;
    double Rho = 0;

    for (int mo = 0; mo < _nmo; mo++)
    {
        const double occ = get_MO_occ(mo);
        const double docc = 2 * occ;
        if (occ != 0)
        {
            phi_temp = &phi[mo * 7];
            Rho += occ * pow(*phi_temp, 2);
            Grad[0] += docc * *phi_temp * phi_temp[1];
            Grad[1] += docc * *phi_temp * phi_temp[2];
            Grad[2] += docc * *phi_temp * phi_temp[3];
            Hess[0] += docc * (*phi_temp * phi_temp[4] + pow(phi_temp[1], 2));
            Hess[1] += docc * (*phi_temp * phi_temp[5] + pow(phi_temp[2], 2));
            Hess[2] += docc * (*phi_temp * phi_temp[6] + pow(phi_temp[3], 2));
            tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
        }
    }
    Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
    Lap = Hess[0] + Hess[1] + Hess[2];
};

const double WFN::computeLap(
    const std::array<double, 3>& PosGrid // [3] vector with current position on te grid
) const
{
    const int _nmo = get_nmo(false);
    vec phi(7 * _nmo, 0.0);
    double* phi_temp;
    double chi[7]{ 0, 0, 0, 0, 0, 0, 0 };
    double d[3]{ 0, 0, 0 };
    int iat = 0;
    int l[3]{ 0, 0, 0 };
    double ex = 0;
    double xl[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

    for (int j = 0; j < nex; j++)
    {
        iat = get_center(j) - 1;

        constants::type2vector(get_type(j), l);
        d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
        d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
        d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
        double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (temp < constants::exp_cutoff)
            continue;
        ex = exp(temp);
        for (int k = 0; k < 3; k++)
        {
            if (l[k] == 0)
            {
                xl[k][0] = 1.0;
                xl[k][1] = 0.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 1)
            {
                xl[k][0] = d[k];
                xl[k][1] = 1.0;
                xl[k][2] = 0.0;
            }
            else if (l[k] == 2)
            {
                xl[k][0] = d[k] * d[k];
                xl[k][1] = 2 * d[k];
                xl[k][2] = 2;
            }
            else if (l[k] == 3)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d[k];
                xl[k][1] = 3 * d2;
                xl[k][2] = 6 * d[k];
            }
            else if (l[k] == 4)
            {
                double d2 = d[k] * d[k];
                xl[k][0] = d2 * d2;
                xl[k][1] = 4 * d2 * d[k];
                xl[k][2] = 12 * d2;
            }
            else
            {
                return -100;
            }
        }
        double ex2 = 2 * get_exponent(j);
        chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
        chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
        chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
        chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
        double temp_ex = pow(ex2, 2);
        chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
        chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
        chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;
        for (int mo = 0; mo < _nmo; mo++)
        {
            phi_temp = &phi[mo * 7];
            for (int i = 0; i < 7; i++)
                phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i]; // build MO values at this point
        }
    }

    double Hess[3]{ 0, 0, 0 };

    for (int mo = 0; mo < _nmo; mo++)
    {
        const double occ = get_MO_occ(mo);
        const double docc = 2 * occ;
        if (occ != 0)
        {
            phi_temp = &phi[mo * 7];
            Hess[0] += docc * (*phi_temp * phi_temp[4] + pow(phi_temp[1], 2));
            Hess[1] += docc * (*phi_temp * phi_temp[5] + pow(phi_temp[2], 2));
            Hess[2] += docc * (*phi_temp * phi_temp[6] + pow(phi_temp[3], 2));
        }
    }
    return Hess[0] + Hess[1] + Hess[2];
};

const double WFN::computeMO(
    const std::array<double, 3>& PosGrid, // [3] array with current position on the grid
    const int &mo) const
{
    double result = 0.0;
    int iat = 0;
    int l[3]{0, 0, 0};
    double ex = 0;
    double temp = 0;

    // x, y, z and dsqd
    vec2 d(4);
    for (int i = 0; i < 4; i++)
        d[i].resize(ncen);

    for (iat = 0; iat < ncen; iat++)
    {
        d[0][iat] = PosGrid[0] - atoms[iat].get_coordinate(0);
        d[1][iat] = PosGrid[1] - atoms[iat].get_coordinate(1);
        d[2][iat] = PosGrid[2] - atoms[iat].get_coordinate(2);
        d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
    }

    for (int j = 0; j < nex; j++)
    {
        iat = get_center(j) - 1;
        // if (iat != atom) continue;
        constants::type2vector(get_type(j), l);
        temp = -get_exponent(j) * d[3][iat];
        if (temp < constants::exp_cutoff)
            continue;
        ex = exp(temp);
        for (int k = 0; k < 3; k++)
        {
            if (l[k] == 0)
                continue;
            else if (l[k] == 1)
                ex *= d[k][iat];
            else if (l[k] == 2)
                ex *= d[k][iat] * d[k][iat];
            else if (l[k] == 3)
                ex *= pow(d[k][iat], 3);
            else if (l[k] == 4)
                ex *= pow(d[k][iat], 4);
            else if (l[k] == 5)
                ex *= pow(d[k][iat], 5);
        }
        result += MOs[mo].get_coefficient_f(j) * ex; // build MO values at this point
    }
    shrink_vector<vec>(d);
    return result;
}

double Integrate(int &m, double i, double &expn)
{
    int x;
    if (i <= 10)
    {
        if (expn == 0.0)
            return 0.0;
        double a = m + 0.5;
        double term = 1 / a;
        double partsum = term;
        for (x = 2; x < 50; x++)
        {
            a++;
            term *= (i / a);
            partsum += term;
            if (term / partsum < 1E-8)
                return 0.5 * partsum * expn;
        }
    }
    else
    {
        double a = m;
        double b = a + 0.5;
        a -= 0.5;
        double id = 1 / i;
        double approx = 0.88622692 * sqrt(id) * pow(id, m);
        for (x = 0; x < m; x++)
        {
            b--;
            approx *= b;
        }
        double mult = 0.5 * expn * id;
        if (mult == 0)
            return approx;
        double prop = mult / approx;
        double term = 1;
        double partsum = 1;
        for (x = 1; x < i + m; x++)
        {
            term *= a * id;
            partsum += term;
            if (abs(term * prop / partsum) < 1E-8)
                return approx - mult * partsum;
            a--;
        }
    }
    return -1;
};
const double WFN::fj(int &j, int &l, int &m, double &aa, double &bb) const
{
    double temp = 0.0;
    double temp2 = 0.0;
    int a = 0, b = 0;
    for (int i = std::max(0, j - m); i <= std::min(j, l); i++)
    {
        // pre = factorial[l] / factorial[l - i] / factorial[i] * factorial[m] / factorial[m - j + i] / factorial[j - i];
        temp2 = static_cast<double>(pre[j][l][m][i]);
        a = l - i;
        b = m + i - j;
        if (a != 0)
            temp2 *= pow(aa, a);
        if (b != 0)
            temp2 *= pow(bb, b);
        temp += temp2;
    }
    return temp;
};

const double WFN::Afac(int &l, int &r, int &i, double &PC, double &gamma, double &fjtmp) const
{
    double temp = fjtmp * pow(0.25 / gamma, r + i) / Afac_pre[l][r][i];
    int num = l - 2 * r - 2 * i;
    if (num != 0)
        temp *= pow(PC, num);
    if (i % 2 == 1)
        return -temp;
    else
        return temp;
}

bool WFN::read_ptb(const std::filesystem::path &filename, std::ostream &file, const bool debug)
{
    origin = e_origin::ptb;
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
    err_checkf(read_block_from_fortran_binary(inFile, &one), "Error reading initial number", std::cout);
    err_checkf(one != 2, "Error reading first number in the xtb file!", std::cout);

    int infos[4] = { 0, 0, 0, 0 }; //ncent nbf nmomax nprims
    err_checkf(read_block_from_fortran_binary(inFile, infos), "Error reading sizes of data", std::cout);
    int ncent = infos[0];
    int nbf = infos[1];
    int nmomax = infos[2];
    int nprims = infos[3];

    svec atyp(ncent);
    char temp[3]{0, 0, '\0'};
    for (int i = 0; i < ncent; ++i)
    {
        err_checkf(read_block_from_fortran_binary(inFile, temp), "Error reading atom label " + std::to_string(i), std::cout);
        atyp[i] = temp;
        atyp[i].erase(remove(atyp[i].begin(), atyp[i].end(), ' '), atyp[i].end());
    }

    vec x(ncent), y(ncent), z(ncent);
    ivec _charge(ncent);
    for (int i = 0; i < ncent; ++i)
    {
        err_checkf(read_block_from_fortran_binary(inFile, &x[i]), "Error reading atom data for atom " + std::to_string(i), std::cout);
        err_checkf(read_block_from_fortran_binary(inFile, &y[i]), "Error reading atom data for atom " + std::to_string(i), std::cout);
        err_checkf(read_block_from_fortran_binary(inFile, &z[i]), "Error reading atom data for atom " + std::to_string(i), std::cout);
        err_checkf(read_block_from_fortran_binary(inFile, &_charge[i]), "Error reading atom data for atom " + std::to_string(i), std::cout);
    }

    ivec lao(nprims), aoatcart(nprims), ipao(nprims);
    for (int i = 0; i < nprims; ++i) err_checkf(read_block_from_fortran_binary(inFile, &lao[i]), "Error reading basis set information lao of primitive " + std::to_string(i), std::cout);
    for (int i = 0; i < nprims; ++i) err_checkf(read_block_from_fortran_binary(inFile, &aoatcart[i]), "Error reading basis set information aotcart of primitive " + std::to_string(i), std::cout);
    for (int i = 0; i < nprims; ++i) err_checkf(read_block_from_fortran_binary(inFile, &ipao[i]), "Error reading basis set information ipao of primitive " + std::to_string(i), std::cout);

    vec exps(nprims), contr(nprims);
    err_checkf(read_block_from_fortran_binary(inFile, exps.data()), "Error reading exponents!", std::cout);
    err_checkf(read_block_from_fortran_binary(inFile, contr.data()), "Error reading contraction coefs!", std::cout);
    vec occ(nmomax), eval(nmomax);
    err_checkf(read_block_from_fortran_binary(inFile, occ.data()), "Error reading occupancies!", std::cout);
    err_checkf(read_block_from_fortran_binary(inFile, eval.data()), "Error reading energies!", std::cout);
    vec tempvec((size_t)nbf* (size_t)nmomax);
    err_checkf(read_block_from_fortran_binary(inFile, tempvec.data()), "Error reading MO coefficients!", std::cout);
    dMatrix2 momat = reshape<dMatrix2>(tempvec, Shape2D(nmomax, nbf));

    // making it into the wavefunction data
    for (int i = 0; i < ncent; i++)
    {
        err_checkf(push_back_atom(atom(atyp[i], "0000000000000", i, x[i], y[i], z[i], _charge[i])), "Error adding atom to WFN!", file);
    }
    err_checkf(ncen == ncent, "Error adding atoms to WFN!", file);

    int elcount = -get_charge();
    if (debug)
        file << "elcount: " << elcount << std::endl;
    for (int i = 0; i < ncen; i++)
    {
        elcount += get_atom_charge(i);
        elcount -= constants::ECP_electrons_pTB[get_atom_charge(i)];
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
            file << temp_els << std::endl;
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
    //for (int i = beta_els; i < alpha_els; i++)
    //{
    //    occ[i] = 1.0;
    //}
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
        else
            err_checkf(push_back_MO(MO(i, occ[i], eval[i], 1)), "Error adding MO to WFN!", file);
    }
    err_checkf(nmo == nmomax, "Error adding MOs to WFN!", file);

    // we need to generate the primitive coefficients from the contr and exp from the momat, then we cann add them MO-wise
    for (int i = 0; i < nprims; i++)
    {
        vec values;
        for (int j = 0; j < nmomax; j++)
        {
            values.push_back(momat(j, ipao[i] - 1) * contr[i]);
        }
        add_primitive(aoatcart[i], lao[i], exps[i], values.data());
    }
    // To Do: This needs reviewing in light of mdspan!
    //while (momat.size() < nbf) {
    //    momat.push_back(vec(nbf, 0.0));
    //    occ.push_back(0);
    //}
    // build density matrix
    dMatrix2 temp_co = diag_dot(momat, occ, true);
    DM = dot(temp_co, momat, false, false);

    err_checkf(nprims == nex, "Error adding primitives to WFN!", file);
    inFile.close();
    if(debug)
        this->write_wfn("test_convert_from_xtb.wfn", false, false);
    constants::exp_cutoff = std::log(constants::density_accuracy / get_maximum_MO_coefficient());
    return true;
}

const double WFN::computeESP(const std::array<double, 3>& PosGrid, const vec2 &d2) const
{
    double ESP = 0;
    double P[3]{0, 0, 0};
    double Pi[3]{0, 0, 0};
    double Pj[3]{0, 0, 0};
    double PC[3]{0, 0, 0};
    double Fn[11]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double Al[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double Am[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double An[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int maplrsl[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int maplrsm[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int maplrsn[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int l_i[3]{0, 0, 0};
    int l_j[3]{0, 0, 0};
    int iat = 0, jat = 0, MaxFn = 0;
    double ex_sum = 0,
           sqd = 0,
           sqpc = 0,
           prefac = 0,
           expc = 0,
           term = 0,
           addesp = 0,
           fjtmp = 0,
           twoexpc = 0,
           iex = 0,
           jex = 0;

    const int MO = get_nmo(true);
    const int nprim = get_nex();

    double temp;
    int maxl, maxm, maxn;
    vec2 pos(3);
    for (int i = 0; i < 3; i++)
        pos[i].resize(get_ncen());

    for (iat = 0; iat < get_ncen(); iat++)
    {
        pos[0][iat] = atoms[iat].get_coordinate(0);
        pos[1][iat] = atoms[iat].get_coordinate(1);
        pos[2][iat] = atoms[iat].get_coordinate(2);
        ESP += get_atom_charge(iat) * pow(sqrt(pow(PosGrid[0] - pos[0][iat], 2) + pow(PosGrid[1] - pos[1][iat], 2) + pow(PosGrid[2] - pos[2][iat], 2)), -1);
    }

    for (int iprim = 0; iprim < nprim; iprim++)
    {
        iat = get_center(iprim) - 1;
        constants::type2vector(get_type(iprim), l_i);
        iex = get_exponent(iprim);
        for (int jprim = iprim; jprim < nprim; jprim++)
        {
            jat = get_center(jprim) - 1;
            constants::type2vector(get_type(jprim), l_j);
            ex_sum = get_exponent(iprim) + get_exponent(jprim);
            jex = get_exponent(jprim);

            sqd = d2[iat][jat];

            prefac = constants::TWO_PI / ex_sum * exp(-iex * jex * sqd / ex_sum);
            if (prefac < 1E-10)
                continue;

            for (int i = 0; i < 3; i++)
            {
                P[i] = (pos[i][iat] * iex + pos[i][jat] * jex) / ex_sum;
                Pi[i] = P[i] - pos[i][iat];
                Pj[i] = P[i] - pos[i][jat];
                PC[i] = P[i] - PosGrid[i];
            }

            sqpc = pow(PC[0], 2) + pow(PC[1], 2) + pow(PC[2], 2);

            expc = exp(-ex_sum * sqpc);
            MaxFn = 0;
            for (int i = 0; i < 3; i++)
                MaxFn += l_i[i] + l_j[i];
            temp = Integrate(MaxFn, ex_sum * sqpc, expc);
            Fn[MaxFn] = temp;
            twoexpc = 2 * ex_sum * sqpc;
            for (int nu = MaxFn - 1; nu >= 0; nu--)
                Fn[nu] = (expc + twoexpc * Fn[nu + 1]) / (2 * (nu + 1) - 1);

            maxl = -1;
            for (int l = 0; l <= l_i[0] + l_j[0]; l++)
            {
                if (l % 2 != 1)
                    fjtmp = fj(l, l_i[0], l_j[0], Pi[0], Pj[0]); // *factorial[l];
                else
                    fjtmp = -fj(l, l_i[0], l_j[0], Pi[0], Pj[0]); // * factorial[l];
                for (int r = 0; r <= l / 2; r++)
                    for (int s = 0; s <= (l - 2 * r) / 2; s++)
                    {
                        maxl++;
                        Al[maxl] = Afac(l, r, s, PC[0], ex_sum, fjtmp);
                        maplrsl[maxl] = l - 2 * r - s;
                    }
            }
            maxm = -1;
            for (int l = 0; l <= l_i[1] + l_j[1]; l++)
            {
                if (l % 2 != 1)
                    fjtmp = fj(l, l_i[1], l_j[1], Pi[1], Pj[1]); // *factorial[l];
                else
                    fjtmp = -fj(l, l_i[1], l_j[1], Pi[1], Pj[1]); // * factorial[l];
                for (int r = 0; r <= l / 2; r++)
                    for (int s = 0; s <= (l - 2 * r) / 2; s++)
                    {
                        maxm++;
                        Am[maxm] = Afac(l, r, s, PC[1], ex_sum, fjtmp);
                        maplrsm[maxm] = l - 2 * r - s;
                    }
            }
            maxn = -1;
            for (int l = 0; l <= l_i[2] + l_j[2]; l++)
            {
                if (l % 2 != 1)
                    fjtmp = fj(l, l_i[2], l_j[2], Pi[2], Pj[2]); // *factorial[l];
                else
                    fjtmp = -fj(l, l_i[2], l_j[2], Pi[2], Pj[2]); // * factorial[l];
                for (int r = 0; r <= l / 2; r++)
                    for (int s = 0; s <= (l - 2 * r) / 2; s++)
                    {
                        maxn++;
                        An[maxn] = Afac(l, r, s, PC[2], ex_sum, fjtmp);
                        maplrsn[maxn] = l - 2 * r - s;
                    }
            }

            term = 0.0;
            for (int l = 0; l <= maxl; l++)
            {
                if (Al[l] == 0)
                    continue;
                for (int m = 0; m <= maxm; m++)
                {
                    if (Am[m] == 0)
                        continue;
                    for (int n = 0; n <= maxn; n++)
                    {
                        if (An[n] == 0)
                            continue;
                        term += Al[l] * Am[m] * An[n] * Fn[maplrsl[l] + maplrsm[m] + maplrsn[n]];
                    }
                }
            }

            if (term == 0)
                continue;

            if (iprim != jprim)
                term *= 2.0;

            term *= prefac;
            addesp = 0.0;
            for (int mo = 0; mo < MO; mo++)
                addesp += get_MO_occ(mo) * get_MO_coef_f(mo, iprim) * get_MO_coef_f(mo, jprim);
            ESP -= addesp * term;
        }
    }
    return ESP;
};

const double WFN::computeESP_noCore(const std::array<double, 3>& PosGrid, const vec2 &d2) const
{
    double ESP = 0;
    double P[3]{0, 0, 0};
    double Pi[3]{0, 0, 0};
    double Pj[3]{0, 0, 0};
    double PC[3]{0, 0, 0};
    double Fn[11]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double Al[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double Am[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double An[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int maplrsl[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int maplrsm[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int maplrsn[54]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int l_i[3]{0, 0, 0};
    int l_j[3]{0, 0, 0};
    int iat = 0, jat = 0, MaxFn = 0;
    double ex_sum = 0,
           sqd = 0,
           sqpc = 0,
           prefac = 0,
           expc = 0,
           term = 0,
           addesp = 0,
           fjtmp = 0,
           twoexpc = 0,
           iex = 0,
           jex = 0;

    const int MO = get_nmo(true);
    const int nprim = get_nex();

    double temp;
    int maxl, maxm, maxn;
    vec2 pos(3);
    for (int i = 0; i < 3; i++)
        pos[i].resize(get_ncen());

    for (int iprim = 0; iprim < nprim; iprim++)
    {
        iat = get_center(iprim) - 1;
        constants::type2vector(get_type(iprim), l_i);
        iex = get_exponent(iprim);
        for (int jprim = iprim; jprim < nprim; jprim++)
        {
            jat = get_center(jprim) - 1;
            constants::type2vector(get_type(jprim), l_j);
            ex_sum = get_exponent(iprim) + get_exponent(jprim);
            jex = get_exponent(jprim);

            sqd = d2[iat][jat];

            prefac = constants::TWO_PI / ex_sum * exp(-iex * jex * sqd / ex_sum);
            if (prefac < 1E-10)
                continue;

            for (int i = 0; i < 3; i++)
            {
                P[i] = (pos[i][iat] * iex + pos[i][jat] * jex) / ex_sum;
                Pi[i] = P[i] - pos[i][iat];
                Pj[i] = P[i] - pos[i][jat];
                PC[i] = P[i] - PosGrid[i];
            }

            sqpc = pow(PC[0], 2) + pow(PC[1], 2) + pow(PC[2], 2);

            expc = exp(-ex_sum * sqpc);
            MaxFn = 0;
            for (int i = 0; i < 3; i++)
                MaxFn += l_i[i] + l_j[i];
            temp = Integrate(MaxFn, ex_sum * sqpc, expc);
            Fn[MaxFn] = temp;
            twoexpc = 2 * ex_sum * sqpc;
            for (int nu = MaxFn - 1; nu >= 0; nu--)
                Fn[nu] = (expc + twoexpc * Fn[nu + 1]) / (2 * (nu + 1) - 1);

            maxl = -1;
            for (int l = 0; l <= l_i[0] + l_j[0]; l++)
            {
                if (l % 2 != 1)
                    fjtmp = fj(l, l_i[0], l_j[0], Pi[0], Pj[0]);
                else
                    fjtmp = -fj(l, l_i[0], l_j[0], Pi[0], Pj[0]);
                for (int r = 0; r <= l / 2; r++)
                    for (int s = 0; s <= (l - 2 * r) / 2; s++)
                    {
                        maxl++;
                        Al[maxl] = Afac(l, r, s, PC[0], ex_sum, fjtmp);
                        maplrsl[maxl] = l - 2 * r - s;
                    }
            }
            maxm = -1;
            for (int l = 0; l <= l_i[1] + l_j[1]; l++)
            {
                if (l % 2 != 1)
                    fjtmp = fj(l, l_i[1], l_j[1], Pi[1], Pj[1]);
                else
                    fjtmp = -fj(l, l_i[1], l_j[1], Pi[1], Pj[1]);
                for (int r = 0; r <= l / 2; r++)
                    for (int s = 0; s <= (l - 2 * r) / 2; s++)
                    {
                        maxm++;
                        Am[maxm] = Afac(l, r, s, PC[1], ex_sum, fjtmp);
                        maplrsm[maxm] = l - 2 * r - s;
                    }
            }
            maxn = -1;
            for (int l = 0; l <= l_i[2] + l_j[2]; l++)
            {
                if (l % 2 != 1)
                    fjtmp = fj(l, l_i[2], l_j[2], Pi[2], Pj[2]);
                else
                    fjtmp = -fj(l, l_i[2], l_j[2], Pi[2], Pj[2]);
                for (int r = 0; r <= l / 2; r++)
                    for (int s = 0; s <= (l - 2 * r) / 2; s++)
                    {
                        maxn++;
                        An[maxn] = Afac(l, r, s, PC[2], ex_sum, fjtmp);
                        maplrsn[maxn] = l - 2 * r - s;
                    }
            }

            term = 0.0;
            for (int l = 0; l <= maxl; l++)
            {
                if (Al[l] == 0)
                    continue;
                for (int m = 0; m <= maxm; m++)
                {
                    if (Am[m] == 0)
                        continue;
                    for (int n = 0; n <= maxn; n++)
                    {
                        if (An[n] == 0)
                            continue;
                        term += Al[l] * Am[m] * An[n] * Fn[maplrsl[l] + maplrsm[m] + maplrsn[n]];
                    }
                }
            }

            if (term == 0)
                continue;

            if (iprim != jprim)
                term *= 2.0;

            term *= prefac;
            addesp = 0.0;
            for (int mo = 0; mo < MO; mo++)
                addesp += get_MO_occ(mo) * get_MO_coef_f(mo, iprim) * get_MO_coef_f(mo, jprim);
            ESP -= addesp * term;
        }
    }
    return ESP;
};

bool WFN::delete_basis_set()
{
    for (int a = 0; a < get_ncen(); a++)
    {
        atoms[a].clear_shellcount();
        int nr_prim = get_atom_primitive_count(a);
        for (int p = 0; p < nr_prim; p++)
        {
            bool succes = erase_atom_primitive(a, 0);
            if (!succes)
                return false;
        }
    }
    return true;
};

std::string WFN::get_atom_label(const int& nr) const
{
    return atoms[nr].get_label();
};

int WFN::get_atom_ECP_electrons(const int& nr) const
{
    return atoms[nr].get_ECP_electrons();
};

basis_set_entry WFN::get_atom_basis_set_entry(const int& nr, const int& bs) const
{
    return atoms[nr].get_basis_set_entry(bs);
};

bool WFN::erase_atom_primitive(const unsigned int& nr, const unsigned int& nr_prim)
{
    if ((int)nr <= ncen && (int)nr_prim < atoms[nr].get_basis_set_size())
    {
        atoms[nr].erase_basis_set(nr_prim);
        return true;
    }
    else
        return false;
};

std::filesystem::path WFN::get_cube_path(const int& nr) const
{
    return cub[nr].get_path();
};

void WFN::write_cube_file(const int& nr, const std::filesystem::path& filename, const bool& debug) {
    err_checkf(nr < cub.size(), "Wrong cube selected!", std::cout);
    if (cub[nr].get_path() != filename) {
        cub[nr].set_path(filename);
    }
    if (debug) {
        std::cout << "Writing cube file to: " << filename << std::endl;
    }
    cub[nr].write_file(filename);
    if (debug) {
        std::cout << "Cube file written!" << std::endl;
    }
}
void WFN::write_cube_dgrid(const int& nr, const std::filesystem::path& filename, const bool& debug) {
    err_checkf(nr < cub.size(), "Wrong cube selected!", std::cout);
    if (cub[nr].get_path() != filename) {
        cub[nr].set_path(filename);
    }
    if (debug) {
        std::cout << "Writing cube file to: " << filename << std::endl;
    }
    cub[nr].write_file(filename);
    if (debug) {
        std::cout << "Cube file written!" << std::endl;
    }
}
void WFN::write_cube_xdgraph(const int& nr, const std::filesystem::path& filename, const bool& debug) {
    err_checkf(nr < cub.size(), "Wrong cube selected!", std::cout);
    if (cub[nr].get_path() != filename) {
        cub[nr].set_path(filename);
    }
    if (debug) {
        std::cout << "Writing cube file to: " << filename << std::endl;
    }
    cub[nr].write_xdgraph(filename);
    if (debug) {
        std::cout << "Cube file written!" << std::endl;
    }
}
