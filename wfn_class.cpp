#include "wfn_class.h"
#include "convenience.h"
#include "mo_class.h"
#include "cube.h"

using namespace std;


const int ECP_electrons[] = { 0,0,                                                                                             0,
 0,  0,                                                                                                  0,  0,  0,  0,  0,  0,
 0,  0,                                                                                                  0,  0,  0,  0,  0,  0,
 0,  0,                                                          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
28, 28,                                                         28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
46, 46, 46, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60 };

bool debug_wfn = false;
bool debug_wfn_deep = false;

void WFN::fill_pre()
{
  for (int j = 0; j < 9; j++)
    for (int l = 0; l < 5; l++)
      for (int m = 0; m < 5; m++) {
        int imax = min(j, l);
        int imin = max(0, j - m);
        for (int i = imin; i <= imax; i++)
          pre[j][l][m][i] = ft[j] * ft[l] / ft[l - i] / ft[i] * ft[m] / ft[m - j + i] / ft[j - i];
      }
}

void WFN::fill_Afac_pre()
{
  for (int l = 0; l < 9; l++)
    for (int r = 0; r <= l / 2; r++)
      for (int s = 0; s <= (l - 2 * r) / 2; s++)
        Afac_pre[l][r][s] = ft[r] * ft[s] * ft[l - 2 * r - 2 * s];
}

WFN::WFN()
{
  ncen = 0;
  nfunc = 0;
  nmo = 0;
  nex = 0;
  charge = 0;
  multi = 0;
  origin = 0;
  total_energy = 0.0;
  virial_ratio = 0.0;
  d_f_switch = false;
  modified = false;
  distance_switch = false;
  has_ECPs = false;
  fill_pre();
  fill_Afac_pre();
};

WFN::WFN(int given_origin)
{
  ncen = 0;
  nfunc = 0;
  nmo = 0;
  nex = 0;
  charge = 0;
  multi = 0;
  total_energy = 0.0;
  origin = given_origin;
  d_f_switch = false;
  modified = false;
  distance_switch = false;
  has_ECPs = false;
  basis_set_name = " ";
  comment = "Test";
  fill_pre();
  fill_Afac_pre();
};

bool WFN::push_back_atom(const string& label, const double& x, const double& y, const double& z, const int& charge)
{
  ncen++;
  if (charge >= 1) atoms.push_back(atom(label, ncen, x, y, z, charge));
  else {
    atoms.push_back(atom());
    return false;
  }
  return true;
};

bool WFN::push_back_atom(const atom& given)
{
  ncen++;
  atoms.push_back(given);
  return true;
};

bool WFN::erase_atom(const int& nr)
{
  if (ncen < 1) return false;
  int n = nr - 1;
  atoms.erase(atoms.begin() + n);
  ncen--;
  return true;
};

bool WFN::push_back_MO(const int& nr, const double& occ, const double& ener)
{
  nmo++;
  //hier fehlen noch sinnhaftigkeitsfragen
  MOs.push_back(MO(nr, occ, ener));
  return true;
};

bool WFN::push_back_MO(MO& given)
{
  nmo++;
  //hier fehlen noch sinnhaftigkeitsfragen
  MOs.push_back(given);
  return true;
};

void WFN::push_back_MO_coef(const int& nr, const double& value)
{
  err_checkf(nr < nmo, "not enough MOs", std::cout);
  MOs[nr].push_back_coef(value);
};

void WFN::assign_MO_coefs(const int& nr, vector<double>& values)
{
  err_checkf(nr < nmo, "not enough MOs", std::cout);
  MOs[nr].assign_coefs(values);
};

double WFN::get_MO_energy(const int& mo)const
{
  if (mo > nmo) return -1;
  else return MOs[mo].get_energy();
}

bool WFN::push_back_center(const int& cent)
{
  if (cent <= ncen && cent > 0) centers.push_back(cent);
  else return false;
  return true;
};

bool WFN::erase_center(const int& g_nr)
{
  centers.erase(centers.begin() + g_nr - 1);
  return true;
};

string WFN::get_centers(const bool& bohr)
{
  string temp;
  for (int i = 0; i < ncen; i++) {
    temp.append(atoms[i].label);
    temp.append(" ");
    if (bohr) temp.append(to_string(atoms[i].x));
    else temp.append(to_string(atoms[i].x * 0.52917721067));
    temp.append(" ");
    if (bohr) temp.append(to_string(atoms[i].y));
    else temp.append(to_string(atoms[i].y * 0.52917721067));
    temp.append(" ");
    if (bohr) temp.append(to_string(atoms[i].z));
    else temp.append(to_string(atoms[i].z * 0.52917721067));
    temp.append("\n");
  }
  return temp;
};

void WFN::list_centers()
{
  for (int i = 0; i < ncen; i++) {
    cout << atoms[i].nr << " " << atoms[i].label << " "
      << atoms[i].x << " " << atoms[i].y << " "
      << atoms[i].z << " " << atoms[i].charge << endl;
  }
};

bool WFN::push_back_type(const int& type)
{
  types.push_back(type);
  return true;
};

bool WFN::erase_type(const int& nr)
{
  err_checkf(nr >= 1, "Wrong type to erase!", std::cout);
  types.erase(types.begin() + nr-1);
  return true;
};

bool WFN::push_back_exponent(const double& e)
{
  exponents.push_back(e);
  return true;

};

bool WFN::erase_exponent(const int& nr)
{
  if (nr < 1) return false;
  exponents.erase(exponents.begin() + (nr - 1));
  return true;
};

bool WFN::remove_primitive(const int& nr)
{
  nex--;
  if (erase_center(nr) && erase_exponent(nr) && erase_type(nr)) {
    for (int n = 0; n < nmo; n++) MOs[n].erase_coef(nr, nex);
    return true;
  }
  else return false;
};

bool WFN::add_primitive(const int& cent, const int& type, const double& e, double* values)
{
  nex++;
  if (push_back_center(cent) && push_back_type(type) && push_back_exponent(e))
    for (int n = 0; n < nmo; n++) MOs[n].push_back_coef(values[n]);
  else return false;
  return true;
};

void WFN::change_type(const int& nr)
{
  bool end = false;
  while (!end) {
    cout << "Please enter the new type you want to assign: ";
    int new_type = 0;
    cin >> new_type;
    if (new_type > ncen || new_type < 0) {
      cout << "Sorry, wrong input, try again!\n";
      continue;
    }
    types[nr - 1] = new_type;
    end = true;
    set_modified();
  }
};

void WFN::change_exponent(const int& nr)
{
  bool end = false;
  while (!end) {
    cout << "Please enter the new exponent you want to assign: ";
    int new_exp = 0;
    cin >> new_exp;
    if (new_exp > ncen || new_exp < 0) {
      cout << "Sorry, wrong input, try again!\n";
      continue;
    }
    exponents[nr - 1] = new_exp;
    end = true;
    set_modified();
  }
};

void WFN::change_center(const int& nr)
{
  bool end = false;
  while (!end) {
    cout << "Please enter the new center you want to assign: ";
    int new_center = 0;
    cin >> new_center;
    if (new_center > ncen || new_center < 0) {
      cout << "Sorry, wrong input, try again!\n";
      continue;
    }
    centers[nr - 1] = new_center;
    end = true;
    set_modified();
  }
};

void WFN::change_MO_coef(const int& nr)
{
  cout << "Which coefficient out of " << nex << " do you want to change?\n";
  int sel = 0;
  bool end = false;
  while (!end) {
    cin >> sel;
    if (sel > nex || sel <= 0) {
      cout << "sorry this is not inside the possible range\n";
      continue;
    }
    MOs[nr - 1].change_coefficient(sel);
    end = true;
  }
  set_modified();
};

bool WFN::change_MO_coef(const int& nr_mo, const int& nr_primitive, const double& value)
{
  err_checkf(nr_mo >= MOs.size(), "MO doesn't exist!", std::cout);
  return MOs[nr_mo].change_coefficient(nr_primitive, value);
};

void WFN::list_primitives()
{
  for (int i = 0; i < nex; i++) {
    cout << i << " center: " << centers[i] << " type: " << types[i] << " exponent: " << exponents[i] << endl;
  }
};

bool WFN::remove_center(const int& nr)
{
  erase_center(nr);
  try {
    for (int i = 0; i < nex; i++) if (centers[i] == nr) remove_primitive(i);
  }
  catch (...) {
    return false;
  }
  set_modified();
  return true;
}

bool WFN::add_exp(const int& cent, const int& type, const double& e)
{
  nex++;
  if (!push_back_center(cent) || !push_back_type(type) || !push_back_exponent(e)) return false;
  else return true;
};

double WFN::get_MO_coef(const int& nr_mo, const int& nr_primitive) const
{
  err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
  return MOs[nr_mo].get_coefficient(nr_primitive);
};

double WFN::get_MO_coef_f(const int& nr_mo, const int& nr_primitive) const
{
  return MOs[nr_mo].get_coefficient_f(nr_primitive);
};

double* WFN::get_MO_coef_ptr(const int& nr_mo)
{
  return MOs[nr_mo].get_coefficient_ptr();
};

int WFN::get_MO_primitive_count(const int& nr_mo) const
{
  err_checkf(nr_mo < MOs.size() && nr_mo >= 0, "WRONG INPUT!", std::cout);
  return MOs[nr_mo].get_primitive_count();
};

string WFN::hdr(const bool& occupied)
{
  string temp = "GAUSSIAN            ";
  if (!occupied) {
    if (nmo > 100) temp.append(to_string(nmo));
    else if (nmo < 100 && nmo>10) {
      temp.append(" ");
      temp.append(to_string(nmo));
    }
    else if (nmo < 10 && nmo>0) {
      temp.append("  ");
      temp.append(to_string(nmo));
    }
    else if (nmo < 0) return "PROBLEM";
  }
  else {
    unsigned long int occupied_mos = 0;
#pragma omp parallel for reduction(+:occupied_mos)
    for (int i = 0; i < nmo; i++) {
      if (MOs[i].get_occ() > 0) occupied_mos++;
    }
    if (occupied_mos > 100) temp.append(to_string(occupied_mos));
    else if (occupied_mos < 100 && occupied_mos>10) {
      temp.append(" ");
      temp.append(to_string(occupied_mos));
    }
    else if (occupied_mos < 10 && occupied_mos>0) {
      temp.append("  ");
      temp.append(to_string(occupied_mos));
    }
    else if (occupied_mos < 0) return "PROBLEM";
  }
  temp.append(" MOL ORBITALS    ");
  if (nex > 100) temp.append(to_string(nex));
  else if (nex < 100 && nex>10) {
    temp.append(" ");
    temp.append(to_string(nex));
  }
  else if (nex < 10 && nex>0) {
    temp.append("  ");
    temp.append(to_string(nex));
  }
  else if (nex < 0) return "PROBLEM";
  temp.append(" PRIMITIVES      ");
  if (ncen > 100) temp.append(to_string(ncen));
  else if (ncen < 100 && ncen>10) {
    temp.append(" ");
    temp.append(to_string(ncen));
  }
  else if (ncen < 10 && ncen>0) {
    temp.append("  ");
    temp.append(to_string(ncen));
  }
  else if (ncen < 0) return "PROBLEM";
  temp.append(" NUCLEI\n");
  return temp;
};

void WFN::read_known_wavefunction_format(const string& fileName, ostream& file, const bool debug)
{
  if (fileName.find(".wfn") != string::npos) origin = 2, err_checkf(read_wfn(fileName, debug, file), "Problem reading wfn", file);
  else if (fileName.find(".ffn") != string::npos) origin = 4, err_checkf(read_wfn(fileName, debug, file), "Problem reading ffn", file);
  else if (fileName.find(".wfx") != string::npos) origin = 6, err_checkf(read_wfx(fileName, debug, file), "Problem reading wfx", file);
  else if (fileName.find(".fch") != string::npos) {
    origin = 4, err_checkf(read_fchk(fileName, file, debug), "Problem reading fchk", file);
    if (debug) err_checkf(write_wfn("test.wfn", debug, false), "Problem writing test.wfn", file);
  }
  else if (fileName.find(".xyz") != string::npos) origin = 7, err_checkf(read_xyz(fileName, file, debug), "Problem reading xyz", file);
  else if (fileName.find(".molden") != string::npos) origin = 8, err_checkf(read_molden(fileName, file, debug), "Problem reading molden file", file);
  else if (fileName.find(".gbw") != string::npos) origin = 9, err_checkf(read_gbw(fileName, file, debug), "Problem reading gbw file", file);
  else err_checkf(false, "Unknown filetype!", file);
};

bool WFN::read_wfn(const string& fileName, const bool& debug, ostream& file)
{
  if (ncen > 0) {
    file << "There is already a wavefunction loaded, do you want to continue and possibly overwrite the existing wavefunction?" << endl;
    if (!yesno()) return false;
    else file << "okay, carrying on then..." << endl;
  }
  if (debug) {
    debug_wfn = true;
    debug_wfn_deep = true;
  }
  err_checkf(exists(fileName), "Couldn't open or find " + fileName + ", leaving", file);
  ifstream rf(fileName.c_str());
  if (rf.good())
    path = fileName;
  string line;
  rf.seekg(0);
  getline(rf, line);
  comment = line;
  getline(rf, line);
  stringstream stream(line);
  string header_tmp;
  int e_nmo, e_nex, e_nuc = 0; //number of expected MOs, Exponents and nuclei
  stream >> header_tmp >> e_nmo >> header_tmp >> header_tmp >> e_nex >> header_tmp >> e_nuc;
  if (debug) {
    file << "e_nmo: " << e_nmo << ", e_nex: " << e_nex << ", e_nuc : " << e_nuc << endl;
  }
  //----------------------------- Read Atoms ------------------------------------------------------------
  vector<unsigned int> dum_nr, dum_ch;
  dum_nr.resize(e_nuc); dum_ch.resize(e_nuc);
  vector<string> dum_label;
  vector<double> dum_x, dum_y, dum_z;
  dum_x.resize(e_nuc); dum_y.resize(e_nuc); dum_z.resize(e_nuc);
  char tempchar[20];
  size_t length;
  for (int i = 0; i < e_nuc; i++) {
    int dump;
    getline(rf, line);
    if (debug) file << i << ".run, line:" << line << endl;
    length = line.copy(tempchar, 4, 0);
    tempchar[length] = '\0';
    if (debug) file << "First line-copy succesfull" << endl;
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
    //if (debug) {
    //  file << "label:" << dum_label[i]
    //    << " nr: " << dum_nr[i]
    //    << " x: " << dum_x[i]
    //    << " y: " << dum_y[i]
    //    << " z: " << dum_z[i]
    //    << " charge: " << dum_ch[i] << endl;
    //}
  }
  //------------------------------------ Read center assignements -------------------------------------------
  vector<unsigned int> dum_center;
  dum_center.resize(e_nex);
  if (debug) for (int i = 0; i < e_nex; i++) dum_center[i] = 99;
  int run = 0;
  int exnum = 0;
  int dump;
  getline(rf, line);
  while (line.compare(0, 6, "CENTRE") == 0 && !rf.eof()) {
    //if (debug) file << "run: " << run << " exnum: " << exnum << " line: " << line << endl;
    if (exnum + 20 <= e_nex) {
      //if (debug) file << "dum_center: ";
      for (int i = 0; i < 20; i++) {
        length = line.copy(tempchar, 3, 20 + 3 * i);
        tempchar[length] = '\0';
        //				if(debug) file << "tempchar: " << tempchar << endl;
        dum_center[exnum] = stoi(tempchar);
        //if (debug) file << dum_center[exnum] << " ";
        if (dum_center[exnum] > e_nuc) {
          cout << "this center doesn't exist.. some weird problem!\n";
          return false;
        }
        exnum++;
      }
      //if (debug) file << endl << "run: " << run * 20 << "   center: " << dum_center[exnum - 1] << endl;
    }
    else {
      //if (debug) file << "exnum+20>e_nex...\n";
      if (exnum < e_nex) {
        for (int i = 0; i < e_nex % 20; i++) {
          length = line.copy(tempchar, 3, 20 + 3 * i);
          tempchar[length] = '\0';
          //					if(debug) file << "tempchar: " << tempchar << endl;
          dum_center[exnum] = stoi(tempchar);
          //if (debug) file << dum_center[exnum] << endl;
          if (dum_center[exnum] > e_nuc) {
            file << "this center doesn't exist.. some weird problem!\n";
            return false;
          }
          exnum++;
        }
      }
      else {
        getline(rf, line);
        continue;
      }
    }
    getline(rf, line);
    if (exnum > e_nex) {
      file << "run went higher than expected values in center reading, thats suspicius, lets stop here...\n";
      return false;
    }
    run++;
  }
  if (debug) {
    file << exnum << endl;
  }
  if (exnum < e_nex) {
    file << "We have a problem adding center assignements!\n";
    return false;
  }
  if (debug_wfn) file << "finished with centers, moving to types...\n";
  //------------------------------------ Read Types ---------------------------------------------------------
  vector<unsigned int> dum_type;
  dum_type.resize(e_nex);
  run = 0;
  exnum = 0;
  while (line.compare(0, 4, "TYPE") == 0 && !rf.eof()) {
    //if (debug_wfn_deep) file << "run: "
    //  << run << " exnum: "
    //  << exnum << " line: "
    //  << line << endl;
    if (exnum + 20 <= e_nex) {
      //if (debug_wfn_deep) file << "dum_type: ";
      for (int i = 0; i < 20; i++) {
        length = line.copy(tempchar, 2, size_t(21 + 3 * i));
        tempchar[length] = '\0';
        //				if(debug_wfn_deep) file << "tempchar: " << tempchar << endl;
        dum_type[exnum] = stoi(tempchar);
        //if (debug_wfn_deep) file << dum_type[exnum] << " ";
        exnum++;
      }
      //if (debug_wfn_deep) file << endl << "run: " << run * 20 << "   type: " << dum_type[exnum - 1] << endl;
    }
    else if (exnum < e_nex) {
      for (int i = 0; i < e_nex % 20; i++) {
        length = line.copy(tempchar, 2, 21 + 3 * i);
        tempchar[length] = '\0';
        dum_type[exnum] = stoi(tempchar);
        exnum++;
      }
    }
    else {
      getline(rf, line);
      continue;
    }
    getline(rf, line);
    if (exnum > e_nex) {
      file << "exnum went higher than expected values in type reading, thats suspicius, lets stop here...\n";
      return false;
    }
    run++;
  }
  if (exnum < e_nex) {
    file << "We have a problem adding type assignements!\n";
    return false;
  }
  if (debug_wfn) file << "finished with types, reading exponents now...\n";
  //----------------------------- Read exponents -------------------------------
  vector<double> dum_exp;
  dum_exp.resize(e_nex);
  run = 0;
  exnum = 0;
  string replace = "E";
  bool three_exponents = false;
  while (line.compare(0, 9, "EXPONENTS") == 0 && !rf.eof()) {
    //if (debug_wfn_deep) file << "run: "
    //  << run << " exnum: "
    //  << exnum << " line: "
    //  << line << endl;
    if (exnum + 5 <= e_nex) {
      if (exnum == 0) {
        const char test = line.at(10);
        string temp_str(" ");
        const char empty = temp_str.at(0);
        if (test != empty)
          three_exponents = true;
        //if (debug_wfn) file << "Three exp: " << three_exponents << endl;
      }
      //if (debug_wfn_deep) file << "dum_exp:" << endl;
      for (int i = 0; i < 5; i++) {

        if (!three_exponents) {
          line.replace(20 + i * 14, 1, replace);
          length = line.copy(tempchar, 13, 11 + 14 * i);
          tempchar[length] = '\0';
          //if (debug_wfn_deep) file << "	tempchar: " << tempchar << " ";
          dum_exp[exnum] = stod(tempchar);
          //if (debug_wfn_deep) file << dum_exp[exnum] << " ";
          exnum++;
        }
        else {
          line.replace(19 + i * 14, 1, replace);
          length = line.copy(tempchar, 14, 10 + 14 * i);
          tempchar[length] = '\0';
          //if (debug_wfn_deep) file << "	tempchar: " << tempchar << endl;
          dum_exp[exnum] = stod(tempchar);
          //if (debug_wfn_deep) file << dum_exp[exnum] << " ";
          exnum++;
        }
      }
      //if (debug_wfn_deep) file << endl;
      //if (debug_wfn_deep) file << "run: " << run * 5 << "   exponent: " << dum_exp[exnum - 1] << endl;
    }
    else if (exnum < e_nex) {
      for (int i = 0; i < e_nex % 5; i++) {
        if (!three_exponents) {
          line.replace(20 + i * 14, 1, replace);
          length = line.copy(tempchar, 13, 11 + 14 * i);
          tempchar[length] = '\0';
          //if (debug_wfn_deep) file << "	tempchar: " << tempchar << endl;
          dum_exp[exnum] = stod(tempchar);
          //if (debug_wfn_deep) file << dum_exp[exnum] << " ";
          exnum++;
        }
        else {
          line.replace(19 + i * 14, 1, replace);
          length = line.copy(tempchar, 14, 10 + 14 * i);
          tempchar[length] = '\0';
          //if (debug_wfn_deep) file << "	tempchar: " << tempchar << endl;
          dum_exp[exnum] = stod(tempchar);
          //if (debug_wfn_deep) file << dum_exp[exnum] << " ";
          exnum++;
        }
      }
    }
    else {
      getline(rf, line);
      continue;
    }
    getline(rf, line);
    if (exnum > e_nex) {
      file << "exnum went higher than expected values in exponent reading, thats suspicius, lets stop here...\n";
      return false;
    }
    run++;
  }
  if (exnum < e_nex) {
    file << "We have a problem adding exponents!\n";
    return false;
  }
  if (debug_wfn) {
    file << "finished with exponents, reading MOs now...\n";
    file << "line: " << line << endl;
  }
  int linecount = 0;
  exnum = 0;
  int monum = 0;
  vector< vector<double> > temp_val;
  temp_val.resize(e_nmo);
  for (int i = 0; i < e_nmo; i++)
    temp_val[i].resize(e_nex);
  //-------------------------------- Read MOs --------------------------------------
  bool orca_switch;
  int temp_orca = check_order(debug_wfn);
  if (temp_orca % 10 == 3)
    orca_switch = true;
  while (!(line.compare(0, 3, "END") == 0) && !rf.eof()) {
    bool b = 0;
    if (monum == e_nmo) {
      //if (debug_wfn) file << "read all MOs I expected, finishing read...." << endl;
      b = true;
      break;
    }
    stringstream stream(line);
    string tmp;
    int temp_nr = 0;
    double temp_occ = -1.0;
    double temp_ener = 0.0;
    /*if(!orca_switch)
      stream >> tmp >> temp_nr >> tmp >> tmp >> tmp >> tmp >> tmp >> temp_occ >> tmp >> tmp >> tmp >> temp_ener;
    else
      stream >> tmp >> temp_nr >> tmp >> tmp >> tmp >> temp_occ >> tmp >> tmp >> tmp >> temp_ener;
    */
    if (temp_nr == 0) {
      length = line.copy(tempchar, 6, 2);
      tempchar[length] = '\0';
      temp_nr = stoi(tempchar);
    }
    if (temp_occ == -1.0) {
      length = line.copy(tempchar, 12, 36);
      tempchar[length] = '\0';
      temp_occ = stod(tempchar);
    }
    if (temp_ener == 0) {
      if (temp_orca % 10 == 3)
        length = line.copy(tempchar, 12, 63);
      else if (temp_orca % 10 == 1)
        length = line.copy(tempchar, 12, 62);
      else
        length = line.copy(tempchar, 12, 62);
      tempchar[length] = '\0';
      temp_ener = stod(tempchar);
    }
    //if (debug_wfn_deep)
    //  file << endl << "tempchar: " << tempchar << " temp_ener: " << temp_ener << endl;
    push_back_MO(temp_nr, temp_occ, temp_ener);
    //if (debug_wfn_deep) {
    //  file << endl << "This is the header for the new MO: " << MOs[monum].hdr() << endl;
    //}
    //---------------------------Start reading MO coefficients-----------------------
    getline(rf, line);
    linecount = 0;
    exnum = 0;
    while (!(line.compare(0, 2, "MO") == 0) && !rf.eof()) {
      //if (debug_wfn_deep) file << "linecount: "
      //  << linecount << " exnum: "
      //  << exnum << " monum: "
      //  << monum << " line: "
      //  << line << endl;
      if (exnum + 5 <= e_nex) {
        for (int i = 0; i < 5; i++) {
          if (!three_exponents)
            line.replace(12 + i * 16, 1, replace);
          else
            line.replace(11 + i * 16, 1, replace);
          length = line.copy(tempchar, 16, 16 * i);
          tempchar[length] = '\0';
          //if (debug_wfn_deep)
          //	file << "	tempchar: " << tempchar << " ";
          temp_val[monum][exnum] = stod(tempchar);
          //if (debug_wfn_deep)
          //  file << temp_val[monum][exnum];
          exnum++;
        }
      }
      else if (exnum < e_nex) {
        if (debug_wfn) file << "This should be the last line! e_nex: " << e_nex << " exnum: " << exnum << endl;
        for (int i = 0; i < (e_nex % 5); i++) {
          if (!three_exponents) {
            line.replace(12 + i * 16, 1, replace);
            length = line.copy(tempchar, 15, 1 + 16 * i);
            tempchar[length] = '\0';
            //if (debug_wfn_deep)
              //file << "	tempchar: " << tempchar << " ";
            temp_val[monum][exnum] = stod(tempchar);
            //if (debug_wfn_deep)
            //  file << temp_val[monum][exnum];
          }
          else {
            line.replace(11 + i * 16, 1, replace);
            length = line.copy(tempchar, 16, 16 * i);
            tempchar[length] = '\0';
            //if (debug_wfn_deep)
              //file << "	tempchar: " << tempchar << " ";
            temp_val[monum][exnum] = stod(tempchar);
            //if (debug_wfn_deep)
            //  file << temp_val[monum][exnum];
          }
          exnum++;
        }
      }
      else {
        getline(rf, line);
        continue;
      }
      //if (debug_wfn_deep) file << endl;
      getline(rf, line);
      if (linecount * 5 > e_nex + 1) {
        file << "linecount went higher than expected values in exponent reading, thats suspicius, lets stop here...\n";
        return false;
      }
      run++;
    }
    monum++;
    //if (debug_wfn_deep) file << "monum after: " << monum << endl;
    if (b == true) break;
  }
  err_checkf(monum + 1 >= e_nmo, "less MOs than expected, quitting...\nmonum: " + to_string(monum) + " e_nmo : " + to_string(e_nmo), file);
  //---------------------Start writing everything from the temp arrays into wave ---------------------

  for (int i = 0; i < e_nuc; i++)
    err_checkf(push_back_atom(dum_label[i], dum_x[i], dum_y[i], dum_z[i], dum_ch[i]), "Error while making atoms!!\n", file);
  for (int j = 0; j < e_nex; j++) {
    err_checkf(add_exp(dum_center[j], dum_type[j], dum_exp[j]), "Error while writing MO coefficients...\n", file);
  }
  for (int i = 0; i < e_nmo; i++) {
    for (int j = 0; j < e_nex; j++) {
      MOs[i].push_back_coef(temp_val[i][j]);
    }
  }
  return true;
};

bool WFN::read_xyz(const string& filename, ostream& file, const bool debug)
{
  err_checkf(exists(filename), "Couldn't open or find " + filename + ", leaving", file);
  origin = 7;
  ifstream rf(filename.c_str());
  if (rf.good())
    path = filename;
  string line;
  rf.seekg(0);

  getline(rf, line);
  stringstream stream(line);
  int e_nuc = 0; //number of expected MOs, Exponents and nuclei
  stream >> e_nuc;
  if (debug)
    file << "e_nuc: " << e_nuc << endl;
  getline(rf, line);
  comment = line;
  //----------------------------- Read Atoms ------------------------------------------------------------
  vector<unsigned int> dum_nr, dum_ch;
  dum_nr.resize(e_nuc); dum_ch.resize(e_nuc);
  vector<string> dum_label;
  vector<double> dum_x, dum_y, dum_z;
  dum_x.resize(e_nuc); dum_y.resize(e_nuc); dum_z.resize(e_nuc); dum_label.resize(e_nuc);
  for (int i = 0; i < e_nuc; i++) {
    vector <string> temp;
    getline(rf, line);
    stream.str(line);
    if (debug) file << i << ".run, line:" << line << endl;
    dum_nr[i] = i;
    temp = split_string<string>(line, " ");
    remove_empty_elements(temp);
    dum_label[i] = temp[0];
    dum_x[i] = ang2bohr(stod(temp[1]));
    dum_y[i] = ang2bohr(stod(temp[2]));
    dum_z[i] = ang2bohr(stod(temp[3]));
    dum_ch[i] = get_Z_from_label(dum_label[i].c_str()) + 1;
    if (debug) {
      file << "label:" << dum_label[i]
        << " nr: " << dum_nr[i]
        << " x: " << dum_x[i]
        << " y: " << dum_y[i]
        << " z: " << dum_z[i]
        << " charge: " << dum_ch[i] << endl;
    }
  }
  //---------------------Start writing everything from the temp arrays into wave ---------------------
  if (debug) file << "finished with reading the file, now i'm going to make everything permantent in the wavefunction...\n";

  for (int i = 0; i < e_nuc; i++)
    err_checkf(push_back_atom(dum_label[i], dum_x[i], dum_y[i], dum_z[i], dum_ch[i]), "Error while making atoms!!", file);
  return true;
};

bool WFN::read_wfx(const string& fileName, const bool& debug, ostream& file)
{
  err_checkf(exists(fileName), "Couldn't open or find " + fileName + ", leaving", file);
  ifstream rf(fileName.c_str());
  path = fileName;
  string line;
  rf.seekg(0);
  getline(rf, line);
  while (line.find("<Title>") == string::npos)
    getline(rf, line);
  getline(rf, line);
  if (debug) file << "comment line " << line << endl;
  comment = line;
  rf.seekg(0);
  while (line.find("<Number of Nuclei>") == string::npos)
    getline(rf, line);
  getline(rf, line);
  if (debug) file << line << endl;
  int temp_ncen = stoi(line);
  rf.seekg(0);
  while (line.find("<Number of Primitives>") == string::npos)
    getline(rf, line);
  getline(rf, line);
  if (debug) file << "nex line: " << line << endl;
  int temp_nex = stoi(line);
  rf.seekg(0);
  while (line.find("<Number of Occupied Molecular Orbitals>") == string::npos)
    getline(rf, line);
  getline(rf, line);
  if (debug) file << "nmo line: " << line << endl;
  int temp_nmo = stoi(line);
  rf.seekg(0);
  while (line.find("<Atomic Numbers>") == string::npos)
    getline(rf, line);
  vector <int> nrs;
  while (true) {
    getline(rf, line);
    //if (debug) file << "atom number line: " << line << endl;
    if (line.find("</Atomic Numbers>") != string::npos)
      break;
    nrs.push_back(stoi(line));
  }
  err_checkf(nrs.size() == temp_ncen, "Mismatch in atom number size", file);
  rf.seekg(0);
  while (line.find("<Nuclear Cartesian Coordinates>") == string::npos)
    getline(rf, line);
  vector < vector <double> > pos;
  pos.resize(3);
  double temp[3];
  while (true) {
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
    push_back_atom(atnr2letter(nrs[i]) + to_string(i+1), pos[0][i], pos[1][i], pos[2][i], nrs[i]);
  err_checkf(ncen == temp_ncen, "Mismatch in atom position size", file);
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
  while (true) {
    getline(rf, line);
    if (line.find("</Primitive Centers>") != string::npos)
      break;
    int number = CountWords(line.c_str());
    istringstream is(line);
    int temp;
    for (int i = 0; i < number; i++) {
      is >> temp;
      centers.push_back(temp);
    }
  }
  rf.seekg(0);
  while (line.find("<Primitive Types>") == string::npos)
    getline(rf, line);
  while (true) {
    getline(rf, line);
    if (line.find("</Primitive Types>") != string::npos)
      break;
    int number = CountWords(line.c_str());
    istringstream is(line);
    int temp;
    for (int i = 0; i < number; i++) {
      is >> temp;
      types.push_back(temp);
    }
  }
  rf.seekg(0);
  while (line.find("<Primitive Exponents>") == string::npos)
    getline(rf, line);
  while (true) {
    getline(rf, line);
    bool please_break = false;
    if (line.find("</Primitive Exponents>") != string::npos) {
      if (CountWords(line.c_str()) != 2)
        please_break = true;
      else
        break;
    }
    int number = CountWords(line.c_str());
    if (please_break)
      number -= 2;
    istringstream is(line);
    double temp;
    for (int i = 0; i < number; i++) {
      is >> temp;
      exponents.push_back(temp);
    }
    if (please_break)
      break;
  }
  err_checkf(exponents.size() == temp_nex && centers.size() == temp_nex && types.size() == temp_nex, "Mismatch in numbers! aborting!", file);
  nex = temp_nex;
  rf.seekg(0);
  while (line.find("<Molecular Orbital Occupation Numbers>") == string::npos)
    getline(rf, line);
  vector <double> occ;
  while (true) {
    getline(rf, line);
    bool please_break = false;
    if (line.find("</Molecular Orbital Occupation Numbers>") != string::npos) {
      if (CountWords(line.c_str()) != 4)
        please_break = true;
      else
        break;
    }
    int number = CountWords(line.c_str());
    istringstream is(line);
    double temp;
    for (int i = 0; i < number; i++) {
      is >> temp;
      occ.push_back(temp);
    }
    if (please_break)
      break;
  }
  rf.seekg(0);
  while (line.find("<Molecular Orbital Energies>") == string::npos)
    getline(rf, line);
  vector <double> ener;
  while (true) {
    getline(rf, line);
    bool please_break = false;
    if (line.find("</Molecular Orbital Energies>") != string::npos) {
      if (CountWords(line.c_str()) != 3)
        please_break = true;
      else
        break;
    }
    int number = CountWords(line.c_str());
    istringstream is(line);
    double temp;
    for (int i = 0; i < number; i++) {
      is >> temp;
      ener.push_back(temp);
    }
    if (please_break)
      break;
  }
  for (int i = 0; i < temp_nmo; i++) {
    err_checkf(push_back_MO(i + 1, occ[i], ener[i]), "Error poshing back MO! MO: " + to_string(i), file);
  }
  occ.resize(0);
  ener.resize(0);
  rf.seekg(0);
  while (line.find("<Molecular Orbital Primitive Coefficients>") == string::npos)
    getline(rf, line);
  vector <double> coef;
  while (line.find("</Molecular Orbital Primitive Coefficients>") == string::npos) {
    while (line.find("<MO Number>") == string::npos)
      getline(rf, line);
    getline(rf, line);
    //if (debug) file << "mo Nr line: " << line << endl;
    int nr = stoi(line);
    while (line.find("</MO Number>") == string::npos)
      getline(rf, line);
    while (coef.size() != nex) {
      getline(rf, line);
      //if (nr == 1 && debug) file << "first MO Coef lines: " << line << endl;
      int number = CountWords(line.c_str());
      istringstream is(line);
      double temp;
      for (int i = 0; i < number; i++) {
        is >> temp;
        coef.push_back(temp);
      }
      err_checkf(coef.size() <= nex, "Error reading coefficients! MO: " + to_string(MOs.size()), file);
    }
    for (int i = 0; i < nex; i++)
      MOs[nr - 1].push_back_coef(coef[i]);
    coef.resize(0);
    getline(rf, line);
  }
  while (line.find("<Energy =") == string::npos)
    getline(rf, line);
  getline(rf, line);
  total_energy = stod(line);
  while (line.find("<Virial Ratio") == string::npos)
    getline(rf, line);
  getline(rf, line);
  virial_ratio = stod(line);
  rf.close();
  return true;
};

bool WFN::read_molden(const string& filename, ostream& file, const bool debug)
{
  err_checkf(exists(filename), "couldn't open or find " + filename + ", leaving", file);
  if (debug)
    file << "File is valid, continuing...\n" << GetCurrentDir << endl;
  origin = 8;
  ifstream rf(filename.c_str());
  if (rf.good())
    path = filename;
  string line;
  rf.seekg(0);
  //d_f_switch = true;

  getline(rf, line);
  err_checkf(line.find("Molden Format") != string::npos, "Does not look like proper molden format file!", file);
  getline(rf, line);
  comment = split_string<string>(line, "]")[1];
  bool au_bohr = false; //au = false, angs = true;
  while (line.find("[Atoms]") == string::npos) {
    getline(rf, line);
  }
  if (split_string<string>(line, "]")[1].find("angs") != string::npos) au_bohr = true;
  else if (split_string<string>(line, "]")[1].find("Angs") != string::npos) au_bohr = true;
  getline(rf, line);
  vector<string> temp;
  while (line.find("]") == string::npos) {
    temp = split_string<string>(line, " ");
    remove_empty_elements(temp);
    if (au_bohr)
      err_checkf(push_back_atom(temp[0],
        ang2bohr(stod(temp[3])),
        ang2bohr(stod(temp[4])),
        ang2bohr(stod(temp[5])),
        stoi(temp[2])), "Error pushing back atom", file);
    else
      err_checkf(push_back_atom(temp[0],
        stod(temp[3]),
        stod(temp[4]),
        stod(temp[5]),
        stoi(temp[2])), "Error pushing back atom", file);
    getline(rf, line);
  }
  err_checkf(line.find("[STO]") == string::npos, "ERROR: STOs are not yet suupported!", file);
  getline(rf, line);
  unsigned int atoms_with_basis = 0;
  while (atoms_with_basis < ncen && line.find("[") == string::npos) {
    vector <string> line_digest = split_string<string>(line, " ");
    remove_empty_elements(line_digest);
    const int atom_based = stoi(line_digest[0]);
    getline(rf, line);
    int shell = 0;
    while (line.size() > 2) {
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
      for (int i = 0; i < number_of_functions; i++) {
        line_digest = split_string<string>(line, " ");
        remove_empty_elements(line_digest);
        err_checkf(atoms[atom_based - 1].push_back_basis_set(stod(line_digest[0]), stod(line_digest[1]), shell_type, shell), "Error pushing back basis", file);
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
  while (line.find("[MO]") == string::npos) {
    if (line.find("[5D]") != string::npos) {
      d5 = true;
    }
    if (line.find("[7F]") != string::npos) {
      f7 = true;
    }
    if (line.find("[9G]") != string::npos) {
      g9 = true;
    }
    getline(rf, line); //Read more lines until we reach MO block
  }
  err_checkf(d5 && f7 && g9, "Molden files without spheircal harmonics not yet implemented!", file);
  int run = 0;
  string sym;
  bool spin; //alpha = false, beta = true
  double ene, occup;
  int expected_coefs = 0;
  vector<primitive> prims;
  vector<int> temp_shellsizes;
  for (int a = 0; a < ncen; a++) {
    int current_shell = -1;
    int l = 0;
    for (int s = 0; s < atoms[a].basis_set.size(); s++) {
      if (atoms[a].basis_set[s].shell != current_shell) {
        if (atoms[a].basis_set[s].type == 1) {
          expected_coefs++;
          l = 0;
        }
        else if (atoms[a].basis_set[s].type == 2) {
          expected_coefs += 3;
          l = 1;
        }
        else if (atoms[a].basis_set[s].type == 3) {
          expected_coefs += 5;
          l = 2;
        }
        else if (atoms[a].basis_set[s].type == 4) {
          expected_coefs += 7;
          l = 3;
        }
        else if (atoms[a].basis_set[s].type == 5) {
          expected_coefs += 9;
          l = 4;
        }
        current_shell++;
      }
      temp_shellsizes.push_back(atoms[a].shellcount[current_shell]);
      prims.push_back(primitive(a + 1,
        atoms[a].basis_set[s].type,
        atoms[a].basis_set[s].exponent,
        atoms[a].basis_set[s].coefficient));
      //prims.back().unnormalize();
    }
  }
  getline(rf, line);
  //vector<double> norm_const = get_norm_const(file);
  int norm_const_run = 0;
  int MO_run = 0;
  vector<vector<double>> d_pure_2_cart;
  vector<vector<double>> f_pure_2_cart;
  vector<vector<double>> g_pure_2_cart;
  err_checkf(generate_sph2cart_mat(d_pure_2_cart, f_pure_2_cart, g_pure_2_cart), "Error creating the conversion matrix", file);
  while (!rf.eof() && rf.good() && line.size() > 2 && line.find("[") == string::npos) {
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
    if (temp[1] == "Alpha" || temp[1] == "alpha") spin = false;
    else spin = true;
    getline(rf, line);
    temp = split_string<string>(line, " ");
    remove_empty_elements(temp);
    occup = stod(temp[1]);
    push_back_MO(run, occup, ene);
    int run_coef = 0;
    int p_run = 0;
    int d_run = 0;
    vector<vector<double>> d_temp(5);
    int f_run = 0;
    vector<vector<double>> f_temp(7);
    int g_run = 0;
    vector<vector<double>> g_temp(9);
    int basis_run = 0;
    for (int i = 0; i < expected_coefs; i++) {
      getline(rf, line);
      temp = split_string<string>(line, " ");
      remove_empty_elements(temp);
      //err_checkf(temp_shellsizes[basis_run] == 1, "Please do not feed me contracted basis sets yet...", file);
      switch (prims[basis_run].type) {
      case 1: {
        for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
          double t = stod(temp[1]) * prims[basis_run + s].coefficient;
          if (abs(t) < 1E-10) t = 0;
          push_back_MO_coef(MO_run, t);
          if (MO_run == 0) {
            push_back_exponent(prims[basis_run + s].exp);
            push_back_center(prims[basis_run].center);
            push_back_type(prims[basis_run].type);
            nex++;
          }
        }
        basis_run += temp_shellsizes[basis_run];
        break;
      }
      case 2: {
        for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
          double t = stod(temp[1]) * prims[basis_run + s].coefficient;
          if (abs(t) < 1E-10) t = 0;
          push_back_MO_coef(MO_run, t);
          if (MO_run == 0) {
            push_back_exponent(prims[basis_run + s].exp);
            push_back_center(prims[basis_run].center);
            push_back_type(prims[basis_run].type + p_run);
            nex++;
          }
        }
        p_run++;
        if (p_run == 3) {
          p_run = 0;
          basis_run += temp_shellsizes[basis_run];
        }
        break;
      }
      case 3: {
        if (d_run == 0) {
          for (int i = 0; i < 5; i++) {
            d_temp[i].resize(temp_shellsizes[basis_run], 0.0);
            fill(d_temp[i].begin(), d_temp[i].end(), 0.0);
          }
        }
        for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
          d_temp[d_run][s] = stod(temp[1]) * prims[basis_run + s].coefficient;
        }
        d_run++;
        if (d_run == 5) {
          for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
            double temp_coef = 0;
            for (int cart = 0; cart < 6; cart++) {
              temp_coef = 0;
              for (int spher = 0; spher < 5; spher++) {
                temp_coef += d_pure_2_cart[cart][spher] * d_temp[spher][s];
              }
              if (abs(temp_coef) < 1E-10) temp_coef = 0;
              push_back_MO_coef(MO_run, temp_coef);
              if (MO_run == 0) {
                push_back_exponent(prims[basis_run + s].exp);
                push_back_center(prims[basis_run].center);
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
      case 4: {
        if (f_run == 0) {
          for (int i = 0; i < 7; i++) {
            f_temp[i].resize(temp_shellsizes[basis_run], 0.0);
            fill(f_temp[i].begin(), f_temp[i].end(), 0.0);
          }
        }
        for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
          f_temp[f_run][s] = stod(temp[1]) * prims[basis_run + s].coefficient;
        }
        f_run++;
        if (f_run == 7) {
          for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
            double temp_coef = 0;
            for (int cart = 0; cart < 10; cart++) {
              temp_coef = 0;
              for (int spher = 0; spher < 7; spher++) {
                temp_coef += f_pure_2_cart[cart][spher] * f_temp[spher][s];
              }
              if (abs(temp_coef) < 1E-10) temp_coef = 0;
              push_back_MO_coef(MO_run, temp_coef);
              if (MO_run == 0) {
                push_back_exponent(prims[basis_run + s].exp);
                push_back_center(prims[basis_run].center);
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
      case 5: {
        if (g_run == 0) {
          for (int i = 0; i < 9; i++) {
            g_temp[i].resize(temp_shellsizes[basis_run], 0.0);
            fill(g_temp[i].begin(), g_temp[i].end(), 0.0);
          }
        }
        for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
          g_temp[g_run][s] = stod(temp[1]) * prims[basis_run + s].coefficient;
        }
        g_run++;
        if (g_run == 9) {
          for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
            double temp_coef = 0;
            for (int cart = 0; cart < 15; cart++) {
              temp_coef = 0;
              for (int spher = 0; spher < 9; spher++) {
                temp_coef += g_pure_2_cart[cart][spher] * g_temp[spher][s];
              }
              if (abs(temp_coef) < 1E-10) temp_coef = 0;
              push_back_MO_coef(MO_run, temp_coef);
              if (MO_run == 0) {
                push_back_exponent(prims[basis_run].exp);
                push_back_center(prims[basis_run].center);
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


  return true;
};

bool WFN::read_gbw(const string& filename, ostream& file, const bool debug)
{
  err_checkf(exists(filename), "couldn't open or find " + filename + ", leaving", file);
  if (debug)
    file << "File is valid, continuing...\n" << GetCurrentDir << endl;
  origin = 9;
  ifstream rf(filename.c_str(), ios::binary);
  if (rf.good())
    path = filename;
  string line;

  //Reading geometry
  rf.seekg(8, ios::beg);
  long int geo_start;
  rf.read((char*)&geo_start, sizeof(geo_start));
  err_checkf(geo_start != 0, "Could not read geometry information location from GBW file!", file);
  rf.seekg(geo_start, ios::beg);
  int at;
  rf.read((char*)&at, 4);
  double geo_vals[6]; //x,y,z, ch, exp_fin_nuc, mass
  int geo_ints[5];
  for (int a = 0; a < at; a++) {
    for (int i = 0; i < 6; i++) {
      rf.read((char*)&(geo_vals[i]), 8);
      err_checkf(rf.good(), "Error reading geo_val", file);
    }
    for (int i = 0; i < 5; i++) {
      rf.read((char*)&(geo_ints[i]), 4);
      err_checkf(rf.good(), "Error reading geo_int", file);
    }
    err_checkf(push_back_atom(atnr2letter(geo_ints[0]),
      geo_vals[0],
      geo_vals[1],
      geo_vals[2],
      geo_ints[0]), "Error pushing back atom", file);
  }

  rf.seekg(16, ios::beg);
  long int basis_start;
  rf.read((char*)&basis_start, 8);
  err_checkf(basis_start != 0, "Could not read beasis information location from GBW file!", file);
  rf.seekg(basis_start, ios::beg);
  int atoms2, temp;
  rf.read((char*)&temp, 4);
  rf.read((char*)&atoms2, 4);
  long unsigned int atoms_with_basis = 0;
  double exp[37];
  double con[37];
  for (int a = 0; a < atoms2; a++) {
    int atom_based, nr_shells;
    rf.read((char*)&atom_based, 4);
    err_checkf(rf.good(), "Error reading atom_based", file);
    rf.read((char*)&nr_shells, 4);
    err_checkf(rf.good(), "Error reading nr_shells", file);
    int shell = 0;
    for (int p = 0; p < nr_shells; p++) {
      int ang_mom, coeff_ind, nr_funct, center;
      rf.read((char*)&ang_mom, 4);
      err_checkf(rf.good(), "Error reading ang_mom", file);
      if (ang_mom >= 5) err_not_impl_f("Higher angular momentum basis functions than G", file);
      rf.read((char*)&coeff_ind, 4);
      err_checkf(rf.good(), "Error reading ceof_ind", file);
      rf.read((char*)&nr_funct, 4);
      err_checkf(rf.good(), "Error reading nr_func", file);
      rf.read((char*)&center, 4);
      err_checkf(rf.good(), "Error reading center", file);
      for (int b = 0; b < 37; b++) {
        rf.read((char*)&(exp[b]), 8);
        err_checkf(rf.good(), "Error reading exp", file);
      }
      for (int b = 0; b < 37; b++) {
        rf.read((char*)&(con[b]), 8);
        err_checkf(rf.good(), "Error reading con", file);
        if (exp[b] != 0 && con[b] != 0) {
          err_checkf(atoms[atom_based].push_back_basis_set(exp[b], con[b], ang_mom + 1, shell), "Error pushing back basis", file);
        }
      }
      shell++;
    }
    atoms_with_basis++;
  }
  int expected_coefs = 0;
  vector<primitive> prims;
  vector<int> temp_shellsizes;
  for (int a = 0; a < ncen; a++) {
    int current_shell = -1;
    int l = 0;
    for (int s = 0; s < atoms[a].basis_set.size(); s++) {
      if (atoms[a].basis_set[s].shell != current_shell) {
        if (atoms[a].basis_set[s].type == 1) {
          expected_coefs++;
          l = 0;
        }
        else if (atoms[a].basis_set[s].type == 2) {
          expected_coefs += 3;
          l = 1;
        }
        else if (atoms[a].basis_set[s].type == 3) {
          expected_coefs += 5;
          l = 2;
        }
        else if (atoms[a].basis_set[s].type == 4) {
          expected_coefs += 7;
          l = 3;
        }
        else if (atoms[a].basis_set[s].type == 5) {
          expected_coefs += 9;
          l = 4;
        }
        current_shell++;
      }
      temp_shellsizes.push_back(atoms[a].shellcount[current_shell]);
      prims.push_back(primitive(a + 1,
        atoms[a].basis_set[s].type,
        atoms[a].basis_set[s].exponent,
        atoms[a].basis_set[s].coefficient));
    }
  }
  int norm_const_run = 0;
  int MO_run = 0;
  vector<vector<double>> d_pure_2_cart;
  vector<vector<double>> f_pure_2_cart;
  vector<vector<double>> g_pure_2_cart;
  err_checkf(generate_sph2cart_mat(d_pure_2_cart, f_pure_2_cart, g_pure_2_cart), "Error creating the conversion matrix", file);

  rf.seekg(24, ios::beg);
  long int MOs_start;
  rf.read((char*)&MOs_start, 8);
  err_checkf(rf.good(), "Error reading MO_start", file);
  err_checkf(MOs_start != 0, "Could not read MO information location from GBW file!", file);
  rf.seekg(MOs_start, ios::beg);
  int operators, dimension;
  rf.read((char*)&operators, 4);
  err_checkf(rf.good(), "Error reading operators", file);
  rf.read((char*)&dimension, 4);
  err_checkf(rf.good(), "Error reading dimnesion", file);
  size_t coef_nr = size_t(dimension) * size_t(dimension);
  vector<vector<double>> coefficients(operators);
  vector<vector<double>> occupations(operators);
  vector<vector<double>> energies(operators);
  vector<vector<int>> irreps(operators);
  vector<vector<int>> cores(operators);
  for (int i = 0; i < operators; i++) {
    coefficients[i].resize(coef_nr);
    occupations[i].resize(dimension);
    energies[i].resize(dimension);
    irreps[i].resize(dimension);
    cores[i].resize(dimension);
    rf.read((char*)coefficients[i].data(), 8 * coef_nr);
    err_checkf(rf.good(), "Error reading coefficients", file);
    rf.read((char*)occupations[i].data(), 8 * dimension);
    err_checkf(rf.good(), "Error reading occupations", file);
    rf.read((char*)energies[i].data(), 8 * dimension);
    err_checkf(rf.good(), "Error reading energies", file);
    rf.read((char*)irreps[i].data(), 4 * dimension);
    err_checkf(rf.good(), "Error reading irreps", file);
    rf.read((char*)cores[i].data(), 4 * dimension);
    err_checkf(rf.good(), "Error reading cores", file);
    for (int j = 0; j < dimension; j++) {
      push_back_MO(i * dimension + j + 1, occupations[i][j], energies[i][j]);
      int run_coef = 0;
      int p_run = 0;
      int d_run = 0;
      vector<vector<double>> d_temp(5);
      int f_run = 0;
      vector<vector<double>> f_temp(7);
      int g_run = 0;
      vector<vector<double>> g_temp(9);
      int basis_run = 0;
      for (int p = 0; p < expected_coefs; p++) {
        switch (prims[basis_run].type) {
        case 1: {
          for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
            double t = coefficients[i][j + p * dimension] * prims[basis_run + s].coefficient;
            if (abs(t) < 1E-10) t = 0;
            push_back_MO_coef(MO_run, t);
            if (MO_run == 0) {
              push_back_exponent(prims[basis_run + s].exp);
              push_back_center(prims[basis_run].center);
              push_back_type(prims[basis_run].type);
              nex++;
            }
          }
          basis_run += temp_shellsizes[basis_run];
          break;
        }
        case 2: {
          for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
            double t = coefficients[i][j + p * dimension] * prims[basis_run + s].coefficient;
            if (abs(t) < 1E-10) t = 0;
            push_back_MO_coef(MO_run, t);
            if (MO_run == 0) {
              push_back_exponent(prims[basis_run + s].exp);
              push_back_center(prims[basis_run].center);
              if (p_run == 0)
                push_back_type(prims[basis_run].type + 2);
              else if (p_run == 1)
                push_back_type(prims[basis_run].type);
              else if (p_run == 2)
                push_back_type(prims[basis_run].type + 1);
              nex++;
            }
          }
          p_run++;
          if (p_run == 3) {
            p_run = 0;
            basis_run += temp_shellsizes[basis_run];
          }
          break;
        }
        case 3: {
          if (d_run == 0) {
            for (int i = 0; i < 5; i++) {
              d_temp[i].resize(temp_shellsizes[basis_run], 0.0);
              fill(d_temp[i].begin(), d_temp[i].end(), 0.0);
            }
          }
          for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
            d_temp[d_run][s] = coefficients[i][j + p * dimension] * prims[basis_run + s].coefficient;
          }
          d_run++;
          if (d_run == 5) {
            for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
              double temp_coef = 0;
              for (int cart = 0; cart < 6; cart++) {
                temp_coef = 0;
                for (int spher = 0; spher < 5; spher++) {
                  temp_coef += d_pure_2_cart[cart][spher] * d_temp[spher][s];
                }
                if (abs(temp_coef) < 1E-10) temp_coef = 0;
                push_back_MO_coef(MO_run, temp_coef);
                if (MO_run == 0) {
                  push_back_exponent(prims[basis_run + s].exp);
                  push_back_center(prims[basis_run].center);
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
        case 4: {
          if (f_run == 0) {
            for (int i = 0; i < 7; i++) {
              f_temp[i].resize(temp_shellsizes[basis_run], 0.0);
              fill(f_temp[i].begin(), f_temp[i].end(), 0.0);
            }
          }
          for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
            f_temp[f_run][s] = coefficients[i][j + p * dimension] * prims[basis_run + s].coefficient;
          }
          f_run++;
          if (f_run == 7) {
            for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
              double temp_coef = 0;
              for (int cart = 0; cart < 10; cart++) {
                temp_coef = 0;
                for (int spher = 0; spher < 7; spher++) {
                  temp_coef += f_pure_2_cart[cart][spher] * f_temp[spher][s];
                }
                if (abs(temp_coef) < 1E-10) temp_coef = 0;
                push_back_MO_coef(MO_run, temp_coef);
                if (MO_run == 0) {
                  push_back_exponent(prims[basis_run + s].exp);
                  push_back_center(prims[basis_run].center);
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
        case 5: {
          if (g_run == 0) {
            for (int i = 0; i < 9; i++) {
              g_temp[i].resize(temp_shellsizes[basis_run], 0.0);
              fill(g_temp[i].begin(), g_temp[i].end(), 0.0);
            }
          }
          for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
            g_temp[g_run][s] = coefficients[i][j + p * dimension] * prims[basis_run + s].coefficient;
          }
          g_run++;
          if (g_run == 9) {
            for (int s = 0; s < temp_shellsizes[basis_run]; s++) {
              double temp_coef = 0;
              for (int cart = 0; cart < 15; cart++) {
                temp_coef = 0;
                for (int spher = 0; spher < 9; spher++) {
                  temp_coef += g_pure_2_cart[cart][spher] * g_temp[spher][s];
                }
                if (abs(temp_coef) < 1E-10) temp_coef = 0;
                push_back_MO_coef(MO_run, temp_coef);
                if (MO_run == 0) {
                  push_back_exponent(prims[basis_run].exp);
                  push_back_center(prims[basis_run].center);
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
    }
  }
  return true;
};

vector<double> WFN::get_norm_const(ofstream& file, bool debug)
{
  err_checkf(get_nr_basis_set_loaded() != 0, "No basis set loaded!", file);
  err_checkf(get_nr_basis_set_loaded() == get_ncen(), "Not all atoms have a basis set loaded!", file);
  vector <double> norm_const;
  //-------------------normalize the basis set shell wise into a copy vector---------
  vector <vector <double> > basis_coefficients;
  basis_coefficients.resize(ncen);
  for (int a = 0; a < ncen; a++)
    for (int p = 0; p < get_atom_primitive_count(a); p++) {
      double temp = get_atom_basis_set_exponent(a, p);
      switch (get_atom_primitive_type(a, p)) {
      case 1:
        temp = 2 * temp / PI;
        temp = pow(temp, 0.75);
        temp = temp * get_atom_basis_set_coefficient(a, p);
        basis_coefficients[a].push_back(temp);
        break;
      case 2:
        temp = 128 * pow(temp, 5);
        temp = temp / pow(PI, 3);
        temp = pow(temp, 0.25);
        temp = get_atom_basis_set_coefficient(a, p) * temp;
        basis_coefficients[a].push_back(temp);
        break;
      case 3:
        temp = 2048 * pow(temp, 7);
        temp = temp / (9 * pow(PI, 3));
        temp = pow(temp, 0.25);
        temp = get_atom_basis_set_coefficient(a, p) * temp;
        basis_coefficients[a].push_back(temp);
        break;
      case 4:
        temp = 32768 * pow(temp, 9);
        temp = temp / (225 * pow(PI, 3));
        temp = pow(temp, 0.25);
        temp = get_atom_basis_set_coefficient(a, p) * temp;
        basis_coefficients[a].push_back(temp);
        break;
      case -1:
        cout << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy..." << endl;
        break;
      }
    }
  for (int a = 0; a < ncen; a++) {
    double aiaj = 0.0;
    double factor = 0.0;
    for (int s = 0; s < get_atom_shell_count(a); s++) {
      int type_temp = get_shell_type(a, s);
      if (type_temp == -1) {
        cout << "ERROR in type assignement!!" << endl;
      }
      if (debug) {
        cout << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << endl;
        int testcount = 0;
        cout << "start: " << get_shell_start(a, s, false) << flush;
        cout << " stop: " << get_shell_end(a, s, false) << flush << endl;
        cout << "factor: ";
      }
      switch (type_temp) {
      case 1:
        factor = 0;
        for (int i = get_shell_start(a, s, false); i <= get_shell_end(a, s, false); i++) {
          for (int j = get_shell_start(a, s, false); j <= get_shell_end(a, s, false); j++) {
            double aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
            double term = (PI / aiaj);
            term = pow(term, 1.5);
            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
          }
        }
        err_checkf(factor != 0, "Factor of 0 is unphysical!", file);
        factor = pow(factor, -0.5);
        if (debug) cout << factor << endl;
        for (int i = get_shell_start(a, s, false); i <= get_shell_end(a, s, false); i++) {
          if (debug) {
            cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << endl;
            cout << "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << endl;
          }
          basis_coefficients[a][i] *= factor;
          norm_const.push_back(basis_coefficients[a][i]);
        }
        break;
      case 2:
        factor = 0;
        for (int i = get_shell_start(a, s, false); i <= get_shell_end(a, s, false); i++) {
          for (int j = get_shell_start(a, s, false); j <= get_shell_end(a, s, false); j++) {
            double aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
            double term = 4 * pow(aiaj, 5);
            term = pow(PI, 3) / term;
            term = pow(term, 0.5);
            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
          }
        }
        err_checkf(factor != 0, "Factor of 0 is unphysical!", file);
        factor = pow(factor, -0.5);
        if (debug) cout << factor << endl;
        for (int i = get_shell_start(a, s, false); i <= get_shell_end(a, s, false); i++) {
          if (debug) {
            cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << endl;
            cout << "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << endl;
          }
          basis_coefficients[a][i] *= factor;
          for (int k = 0; k < 3; k++) norm_const.push_back(basis_coefficients[a][i]);
        }
        break;
      case 3:
        factor = 0;
        for (int i = get_shell_start(a, s, false); i <= get_shell_end(a, s, false); i++) {
          for (int j = get_shell_start(a, s, false); j <= get_shell_end(a, s, false); j++) {
            double aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
            double term = 16 * pow(aiaj, 7);
            term = pow(PI, 3) / term;
            term = pow(term, 0.5);
            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
          }
        }
        err_checkf(factor != 0, "Factor of 0 is unphysical!", file);
        factor = (pow(factor, -0.5)) / sqrt(3);
        if (debug) cout << factor << endl;
        for (int i = get_shell_start(a, s, false); i <= get_shell_end(a, s, false); i++) {
          if (debug) {
            cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << endl;
            cout << "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << endl;
          }
          basis_coefficients[a][i] *= factor;
          for (int k = 0; k < 3; k++) norm_const.push_back(basis_coefficients[a][i]);
          for (int k = 0; k < 3; k++) norm_const.push_back(sqrt(3) * basis_coefficients[a][i]);
        }
        break;
      case 4:
        factor = 0;
        for (int i = get_shell_start(a, s, false); i <= get_shell_end(a, s, false); i++) {
          for (int j = get_shell_start(a, s, false); j <= get_shell_end(a, s, false); j++) {
            double aiaj = get_atom_basis_set_exponent(a, i) + get_atom_basis_set_exponent(a, j);
            double term = 64 * pow((aiaj), 9);
            term = pow(PI, 3) / term;
            term = pow(term, 0.5);
            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
          }
        }
        err_checkf(factor != 0, "Factor of 0 is unphysical!", file);
        factor = pow(factor, -0.5) / sqrt(15);
        if (debug) cout << factor << endl;
        for (int i = get_shell_start(a, s, false); i <= get_shell_end(a, s, false); i++) {
          if (debug) {
            cout << "Contraction coefficient before: " << get_atom_basis_set_coefficient(a, i) << endl;
            cout << "Contraction coefficient after:  " << factor * get_atom_basis_set_coefficient(a, i) << endl;
          }
          basis_coefficients[a][i] *= factor;
          for (int l = 0; l < 3; l++) 	norm_const.push_back(basis_coefficients[a][i]);
          for (int l = 0; l < 6; l++) 	norm_const.push_back(sqrt(5) * basis_coefficients[a][i]);
          norm_const.push_back(sqrt(15) * basis_coefficients[a][i]);
        }
        break;
      }
      if (debug)	cout << "This shell has: " << get_shell_end(a, s, false) - get_shell_start(a, s, false) + 1 << " primitives" << endl;
    }
  }
  return norm_const;
}

double WFN::get_atom_coordinate(const int& nr, const int& axis)
{
  err_checkf(!(nr < 0 || nr >= ncen || axis < 0 || axis>2), "This input is invalid for get_atom_coordinate!", std::cout);
  switch (axis) {
  case 0: return atoms[nr].x; break;
  case 1: return atoms[nr].y; break;
  case 2: return atoms[nr].z; break;
  default: return -2;
  }
};

bool WFN::write_wfn(const string& fileName, const bool& debug, const bool occupied)
{
  if (debug) {
    debug_wfn = true;
    debug_wfn_deep = true;
    if (exists(fileName)) {
      cout << "File already existed, do you want to overwrite it?";
      if (!yesno()) return false;
    }
    else {
      if (debug_wfn) cout << "File didn't exist before, writing comment to it now." << endl;
    }
  }

  ofstream rf(fileName.c_str(), ios::out);
  string line;
  char choice;
  if (!rf.is_open()) {
    cout << "Sorry, can't open the file...\n";
    return false;
  }
  rf << comment << endl;
  if (debug_wfn) cout << "comment written, now for the header..\n";
  rf << hdr(occupied);
  if (debug_wfn) {
    cout << "header written, now for the centers..\n";
    cout << "this is the header: \n" << hdr(occupied);
  }
  rf.flush();
  for (int i = 0; i < ncen; i++) {
    rf << atoms[i].label;
    if (i < 9) rf << "     ";
    else rf << "    ";
    rf << i + 1 << "    (CENTRE ";
    if (i < 9) rf << ' ';
    rf << i + 1 << ") ";
    rf << fixed << showpoint << setprecision(8);
    rf << setw(12) << atoms[i].x;
    rf << setw(12) << atoms[i].y;
    rf << setw(12) << atoms[i].z;
    rf << "  CHARGE = ";
    rf << fixed << showpoint << setprecision(1) << setw(2) << atoms[i].charge;
    rf << ".0";
    rf << '\n';
  }
  if (debug_wfn) cout << "centers written, now for the center_assignement..\n";
  if (debug_wfn) cout << "ncen: " << ncen << " nex: " << nex << " nmo: " << nmo << endl;
  int run = 0;
  int exnum = 0;
  for (int i = 0; i < nex / 20; i++) {
    rf << "CENTRE ASSIGNMENTS  ";
    for (int j = 0; j < 20; j++) {
      rf << setw(3) << centers[exnum];
      if (exnum > nex) {
        cout << "run is too big in center writing";
        if (debug_wfn) cout << "in 20er-lines...\n";
        return false;
      }
      exnum++;
    }
    run++;
    rf << '\n';
  }
  if (debug_wfn) cout << "this should be the last line... \n";
  if (exnum < nex) {
    rf << "CENTRE ASSIGNMENTS  ";
    for (int j = 0; j < nex % 20; j++) {
      rf << setw(3) << centers[exnum];
      if (exnum > nex) {
        cout << "run is too big in center writing";
        if (debug_wfn) cout << " in last line... trying to access # " << exnum << "\n";
        return false;
      }
      exnum++;
    }
    rf << '\n';
  }
  if (run * 20 < nex / 20 - 1) {
    cout << "Problem during writing of Centre assignments... stopping...\n";
    return false;
  }
  if (debug_wfn) cout << "center assignements written, now for the types..\n";
  run = 0;
  exnum = 0;
  for (int i = 0; i < nex / 20; i++) {
    rf << "TYPE ASSIGNMENTS    ";
    for (int j = 0; j < 20; j++) {
      rf << setw(3) << types[exnum];
      if (exnum > nex) {
        cout << "run is too big in types writing\n";
        return false;
      }
      exnum++;
    }
    run++;
    rf << '\n';
  }
  if (exnum < nex) {
    rf << "TYPE ASSIGNMENTS    ";
    int final_j = 0;
    for (int j = 0; j < nex % 20; j++) {
      rf << setw(3) << types[exnum];
      if (exnum > nex) {
        cout << "run is too big in types writing";
        return false;
      }
      final_j = j;
      exnum++;
    }
    if (debug_wfn) cout << "final_j: " << final_j << endl;
    rf << '\n';
  }
  if (run * 20 < nex / 20 - 1) {
    cout << "Problem during writing of Type assignments... stopping...";
    return false;
  }
  if (debug_wfn) cout << "types assignements written, now for the exponents..\n";
  run = 0;
  exnum = 0;
  for (int i = 0; i < nex / 5; i++) {
    rf << "EXPONENTS ";
    for (int j = 0; j < 5; j++) {
      stringstream stream;
      string temp;
      stream << uppercase << scientific << setw(14) << setprecision(7) << exponents[exnum];
      temp = stream.str();
      rf << temp;
      if (exnum > nex) {
        cout << "run is too big in exponents writing";
        return false;
      }
      exnum++;
    }
    run++;
    rf << '\n';
  }
  if (exnum < nex) {
    rf << "EXPONENTS ";
    for (int j = 0; j < nex % 5; j++) {
      stringstream stream;
      string temp;
      stream << uppercase << scientific << setw(14) << setprecision(7) << exponents[exnum];
      temp = stream.str();
      rf << temp;
      if (run > nex) {
        cout << "run is too big in exponents writing";
        return false;
      }
      exnum++;
    }
    rf << '\n';
  }
  if (run * 5 < nex / 5 - 1) {
    cout << "Problem during writing of Exponents... stopping...";
    return false;
  }
  if (debug_wfn) cout << "exponents assignements written, now for the MOs.." << endl << "For informational purposes: ncen "
    << ncen << " nmo " << nmo << " nex " << nex << endl;
  int mo_run = 1;
  for (int mo_counter = 0; mo_counter < nmo; mo_counter++) {
    if (occupied && MOs[mo_counter].get_occ() == 0) continue;
    //rf << MOs[mo_counter].hdr();
    rf << "MO" << setw(3) << mo_run << setw(29) << "OCC NO =" << setw(13) << fixed << setprecision(8) << MOs[mo_counter].get_occ()
      << setw(14) << "ORB. ENERGY =" << setw(13) << fixed << setprecision(8) << MOs[mo_counter].get_energy() << endl;
    run = 0;
    for (int i = 0; i < nex / 5; i++) {
      for (int j = 0; j < 5; j++) {
        stringstream stream;
        string temp;
        stream << uppercase << scientific << showpoint << setprecision(8) << setw(16) << MOs[mo_counter].get_coefficient(run);
        temp = stream.str();
        rf << temp;
        if (run > nex) {
          cout << "run (" << run << ") is too big in MO ceofficients writing" << endl;
          return false;
        }
        run++;
      }
      rf << '\n';
    }
    if (run < nex) {
      if (debug_wfn_deep) cout << "Still some left to write... going in % for loop...." << endl;
      for (int j = 0; j < nex % 5; j++) {
        stringstream stream;
        string temp;
        stream << uppercase << scientific << showpoint << setprecision(8) << setw(16) << MOs[mo_counter].get_coefficient(run);
        temp = stream.str();
        rf << temp;
        if (run > nex) {
          cout << "run (" << run << ") is too big in MO ceofficients writing" << endl;
          return false;
        }
        run++;
      }
      rf << '\n';
    }
    mo_run++;
  }
  if (run != nex) {
    cout << "Problem during writing of MOs... stopping...";
    if (debug_wfn_deep) cout << "run: " << run << endl;
    return false;
  }
  rf << "END DATA" << endl;
  rf.flush();
  rf.close();
  return true;
};

void WFN::print_primitive(const int& nr)
{
  cout << "center assignement: " << centers[nr] << " type: " << types[nr]
    << " exponent: " << exponents[nr] << endl << "MO coefficients:";
  for (int i = 0; i < nmo; i++) {
    cout << MOs[nr].get_coefficient(i) << "   ";
    if (i % 5 == 0) cout << endl;
  }
};

int WFN::get_nmo(const bool& only_occ) const
{
  if (!only_occ) return nmo;
  else {
    int count = 0;
#pragma omp parallel for reduction(+:count)
    for (int i = 0; i < MOs.size(); i++) {
      if (MOs[i].get_occ() != 0.0) {
        count++;
      }
    }
    return count;
  }
};

unsigned int WFN::get_nr_electrons(const bool& debug)
{
  unsigned int count = 0;
  for (int i = 0; i < ncen; i++)
    count += atoms[i].charge;
  count -= charge;
  return count;
};

unsigned int WFN::get_nr_ECP_electrons() {
  unsigned int count = 0;
  for (int i = 0; i < ncen; i++)
    count += atoms[i].ECP_electrons;
  return count;
}

double WFN::count_nr_electrons(void)
{
  double count = 0;
  for (int i = 0; i < nmo; i++)
    count += MOs[i].get_occ();
  return count;
};

double WFN::get_atom_basis_set_exponent(const int& nr_atom, const int& nr_prim)
{
  if (nr_atom <= ncen && nr_atom >= 0 && atoms[nr_atom].basis_set.size() >= nr_prim && nr_prim >= 0)
    return atoms[nr_atom].basis_set[nr_prim].exponent;
  else return -1;
};

double WFN::get_atom_basis_set_coefficient(const int& nr_atom, const int& nr_prim)
{
  if (nr_atom <= ncen && nr_atom >= 0 && atoms[nr_atom].basis_set.size() >= nr_prim && nr_prim >= 0)
    return atoms[nr_atom].basis_set[nr_prim].coefficient;
  else return -1;
};

bool WFN::change_atom_basis_set_exponent(const int& nr_atom, const int& nr_prim, const double& value)
{
  if (nr_atom <= ncen && nr_atom >= 0 && atoms[nr_atom].basis_set.size() >= nr_prim && nr_prim >= 0) {
    atoms[nr_atom].basis_set[nr_prim].exponent = value;
    set_modified();
    return true;
  }
  else return false;
};

bool WFN::change_atom_basis_set_coefficient(const int& nr_atom, const int& nr_prim, const double& value)
{
  err_checkf(nr_atom <= ncen && nr_atom >= 0 && atoms[nr_atom].basis_set.size() >= nr_prim && nr_prim >= 0, "Wrong input!", cout);
  atoms[nr_atom].basis_set[nr_prim].coefficient = value;
  set_modified();
  return true;
};

int WFN::get_atom_primitive_count(const int& nr)
{
  if (nr <= ncen && nr >= 0) return atoms[nr].basis_set.size();
  else return -1;
};

bool WFN::erase_atom_primitive(const int& nr, const int& nr_prim)
{
  if (nr <= ncen && nr >= 0 && nr_prim >= 0 && nr_prim < atoms[nr].basis_set.size()) {
    atoms[nr].basis_set.erase(atoms[nr].basis_set.begin() + nr_prim);
    return true;
  }
  else return false;
};

int WFN::get_basis_set_shell(const int& nr_atom, const int& nr_prim)
{
  if (nr_atom <= ncen && nr_atom >= 0 && atoms[nr_atom].basis_set.size() >= nr_prim && nr_prim >= 0) {
    return atoms[nr_atom].basis_set[nr_prim].shell;
  }
  else return -1;
};

int WFN::get_atom_shell_count(const int& nr)
{
  if (nr <= ncen && nr >= 0) return atoms[nr].shellcount.size();
  else return -1;
};

int WFN::get_atom_shell_primitives(const int& nr_atom, const int& nr_shell)
{
  if (nr_atom <= ncen && nr_atom >= 0 && nr_shell < atoms[nr_atom].shellcount.size() && nr_shell >= 0)
    return atoms[nr_atom].shellcount[nr_shell];
  else return -1;
};

int WFN::get_shell_type(const int& nr_atom, const int& nr_shell)
{
  if (nr_atom <= ncen && nr_atom >= 0 && nr_shell <= atoms[nr_atom].shellcount.size() && nr_shell >= 0) {
    int primitive_counter = 0;
    while (atoms[nr_atom].basis_set[primitive_counter].shell != nr_shell)
      primitive_counter++;
    return atoms[nr_atom].basis_set[primitive_counter].type;
  }
  else return -1;
};

int WFN::get_shell_center(const int& nr_atom, const int& nr_shell)
{
  if (nr_atom <= ncen && nr_atom >= 0 && nr_shell <= atoms[nr_atom].shellcount.size() && nr_shell >= 0)
    return centers[get_shell_start_in_primitives(nr_atom, nr_shell)];
  else return -1;
};

int WFN::get_shell_start(const int& nr_atom, const int& nr_shell, const bool& debug)
{
  if (nr_atom <= ncen && nr_atom >= 0 && nr_shell <= atoms[nr_atom].shellcount.size() - 1 && nr_shell >= 0) {
    int primitive_counter = 0;
#pragma loop(no_vector)
    for (int s = 0; s < nr_shell; s++)
      primitive_counter += atoms[nr_atom].shellcount[s];
    return primitive_counter;
  }
  else return -1;
};

int WFN::get_shell_start_in_primitives(const int& nr_atom, const int& nr_shell)
{
  if (nr_atom <= ncen && nr_atom >= 0 && nr_shell <= atoms[nr_atom].shellcount.size() - 1 && nr_shell >= 0) {
    int primitive_counter = 0;
    for (int a = 0; a < nr_atom; a++)
      for (int s = 0; s < atoms[a].shellcount.size(); s++)
        switch (get_shell_type(a, s)) {
        case 1:
          primitive_counter += atoms[a].shellcount[s];
          break;
        case 2:
          primitive_counter += (3 * atoms[a].shellcount[s]);
          break;
        case 3:
          primitive_counter += (6 * atoms[a].shellcount[s]);
          break;
        case 4:
          primitive_counter += (10 * atoms[a].shellcount[s]);
          break;
        }
    for (int s = 0; s < nr_shell; s++) {
      switch (get_shell_type(nr_atom, s)) {
      case 1:
        primitive_counter += atoms[nr_atom].shellcount[s];
        break;
      case 2:
        primitive_counter += (3 * atoms[nr_atom].shellcount[s]);
        break;
      case 3:
        primitive_counter += (6 * atoms[nr_atom].shellcount[s]);
        break;
      case 4:
        primitive_counter += (10 * atoms[nr_atom].shellcount[s]);
        break;
      }
    }
    return primitive_counter;
  }
  else return -1;
};

int WFN::get_shell_end(const int& nr_atom, const int& nr_shell, const bool& debug)
{
  if (nr_atom <= ncen && nr_atom >= 0 && nr_shell <= atoms[nr_atom].shellcount.size() && nr_shell >= 0) {
    if (nr_shell == atoms[nr_atom].shellcount.size() - 1) return atoms[nr_atom].basis_set.size() - 1;
    int primitive_counter = 0;
    while (atoms[nr_atom].basis_set[primitive_counter].shell != (nr_shell + 1))	primitive_counter++;
    return primitive_counter - 1;
  }
  else return -1;
};

string WFN::get_atom_label(const int& nr)
{
  string error_return;
  error_return = '?';
  if (nr < ncen && nr >= 0) return atoms[nr].label;
  else return error_return;
};

int WFN::get_nr_basis_set_loaded()
{
  int count = 0;
  for (int a = 0; a < ncen; a++) if (atoms[a].get_basis_set_loaded()) count++;
  return count;
};

bool WFN::get_atom_basis_set_loaded(const int& nr)
{
  if (nr <= ncen && nr >= 0) return atoms[nr].get_basis_set_loaded();
  else {
    cout << "invalid atom choice in atom_basis_set_loaded!" << endl;
    return false;
  }
};

int WFN::get_atom_charge(const int& nr) const
{
  if (nr <= ncen && nr >= 0) return atoms[nr].charge;
  else {
    cout << "invalid atom choice in atom_basis_set_loaded!" << endl;
    return false;
  }
};

void WFN::push_back_DM(const double& value)
{
  DensityMatrix.push_back(value);
};

double WFN::get_DM(const int& nr)
{
  if (nr >= 0 && nr < DensityMatrix.size()) return DensityMatrix[nr];
  else {
    cout << "Requested nr out of range! Size: " << DensityMatrix.size() << " nr: " << nr << endl;
    return -1;
  }
};

bool WFN::set_DM(const int& nr, const double& value)
{
  if (nr >= 0 && nr < DensityMatrix.size()) {
    DensityMatrix[nr] = value;
    return true;
  }
  else {
    cout << "invalid arguments for set_DM! Input was: " << nr << ";" << value << endl;
    return false;
  }
};

void WFN::push_back_SDM(const double& value)
{
  SpinDensityMatrix.push_back(value);
};

double WFN::get_SDM(const int& nr)
{
  if (nr >= 0 && nr < SpinDensityMatrix.size()) return SpinDensityMatrix[nr];
  else {
    cout << "Requested nr out of range! Size: " << SpinDensityMatrix.size() << " nr: " << nr << endl;
    return -1;
  }
};

bool WFN::set_SDM(const int& nr, const double& value)
{
  if (nr >= 0 && nr < SpinDensityMatrix.size()) {
    SpinDensityMatrix[nr] = value;
    return true;
  }
  else {
    cout << "invalid arguments for set_SDM! Input was: " << nr << ";" << value << endl;
    return false;
  }
};

int WFN::check_order(const bool& debug)
{
  for (int i = 0; i < ncen; i++) {
    if (!get_atom_basis_set_loaded(i)) {
      cout << "Sorry, consistency check only works if basis set is loaded for all atoms!" << endl
        << "Failing atom: " << i << " " << get_atom_label(i) << endl;
      return -1;
    }
  }
  //---------------------------check if the wfn is in the right order----------------------
  int order = 0; //1= gaussian (P=2222 3333 4444) 2= tonto (234 234 234 234) 3= ORCA (423 423 423 423)
  int f_order = 0; //1=gaussian (F=11 12 13 17 14 15 18 19 16 20) 2=tonto=ORCA 3=ORCA (11 12 13 14 15 17 16 18 19 20)
  int primcounter = 0;
  bool order_found = false;
  for (int a = 0; a < get_ncen(); a++) {
    for (int s = 0; s < get_atom_shell_count(a); s++) {
      switch (get_shell_type(a, s)) {
      case 1:
        for (int i = 0; i < get_atom_shell_primitives(a, s); i++) {
          if (types[primcounter] != 1) {
            order = -1;
            if (debug) {
              cout << "This should not happen, the order of your file is not ok for S-types! Checked #:" << primcounter << endl;
              return -1;
            }
          }
          else primcounter++;
        }
        break;
      case 2:
        if (order_found) {
          if (order == 1) {
            for (int r = 0; r < get_atom_shell_primitives(a, s); r++) {
              if (types[primcounter] != 2 || types[primcounter + get_atom_shell_primitives(a, s)] != 3
                || types[primcounter + 2 * get_atom_shell_primitives(a, s)] != 4) {
                if (debug) {
                  cout << "The found order does not match all type entries! primcounter: " << primcounter << endl;
                }
                return -1;
              }
              primcounter++;
            }
            primcounter += 2 * get_atom_shell_primitives(a, s);
          }
          else if (order == 2) {
            for (int r = 0; r < get_atom_shell_primitives(a, s); r++) {
              if (types[primcounter] != 2 || types[primcounter + 1] != 3 || types[primcounter + 2] != 4) {
                if (debug) {
                  cout << "The found order does not match all type entries! primcounter: " << primcounter << endl;
                }
                return -1;
              }
              primcounter += 3;
            }
          }
          else if (order == 3)
          {
            for (int r = 0; r < get_atom_shell_primitives(a, s); r++) {
              if (types[primcounter] != 4 || types[primcounter + 1] != 2 || types[primcounter + 2] != 3) {
                if (debug) {
                  cout << "The found order does not match all type entries! primcounter: " << primcounter << endl;
                }
                return -1;
              }
              primcounter += 3;
            }
          }
        }
        else {
          if (types[primcounter] == 2) {
            if (debug && a == 0) {
              cout << "Seems to be either tonto or gaussian file..." << endl;
            }
            order = 1;
            if (types[primcounter + 1] == 3) {
              order = 2;
              order_found = true;
              if (debug) {
                cout << "This wfn file is in tonto order!" << endl;
              }
            }
            else if (types[primcounter + 1] == 2 && get_atom_shell_primitives(a, s) > 1) {
              order = 1;
              order_found = true;
              if (debug) {
                cout << "This wfn file is in gaussian order!" << endl;
              }
            }
            else {
              order = 1;
              order_found = true;
              if (debug) {
                cout << "Either some error or this shell only has 1 p-primitive and "
                  << "i didn't find any order yet... assuming gaussian" << endl;
              }
            }
          }
          else if (types[primcounter] == 4) {
            if (debug && a == 0) {
              cout << "Seems as if this file was ORCA ordered..." << endl;
            }
            order = 3;
            if (types[primcounter + 1] == 2) {
              if (debug) {
                cout << "Yep, that's right! making it permanent now!" << endl;
              }
              order_found = true;
            }
          }
          else {
            order = -1;
            cout << "I can't recognize this order of the .wfn file..." << endl;
            return -1;
          }
          s--;
        }
        break;
      case 3: {
        if (order_found) {
          switch (order) {
          case 1: {
            for (int i = 0; i < get_atom_shell_primitives(a, s); i++) {
              if (types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] != 5 || types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] != 6
                || types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] != 7 || types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] != 8
                || types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] != 9 || types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] != 10) {
                order = -1;
                if (debug) {
                  cout << "The checked types are 6 from #" << primcounter << " and are:" << endl;
                  cout << types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] << endl;
                  return -1;
                }
                cout << "Something seems to be wrong in the order of your D-Types..." << endl;
              }
              else primcounter += 6;
            }
          }
                break;
          case 2:
          case 3: {
            for (int i = 0; i < get_atom_shell_primitives(a, s); i++) {
              if (types[get_shell_start_in_primitives(a, s) + 0 + 6 * i] != 5 || types[get_shell_start_in_primitives(a, s) + 1 + 6 * i] != 6
                || types[get_shell_start_in_primitives(a, s) + 2 + 6 * i] != 7 || types[get_shell_start_in_primitives(a, s) + 3 + 6 * i] != 8
                || types[get_shell_start_in_primitives(a, s) + 4 + 6 * i] != 9 || types[get_shell_start_in_primitives(a, s) + 5 + 6 * i] != 10) {
                order = -1;
                if (debug) {
                  cout << "The checked types are 6 from #" << primcounter << " and are:" << endl;
                  cout << types[get_shell_start_in_primitives(a, s) + 0 + 6 * i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 1 + 6 * i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 2 + 6 * i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 3 + 6 * i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 4 + 6 * i] << " "
                    << types[get_shell_start_in_primitives(a, s) + 5 + 6 * i] << endl;
                  return -1;
                }
                cout << "Something seems to be wrong in the order of your D-Types..." << endl;
              }
              else primcounter += 6;
            }
          }
                break;
          }
        }
        else {
          cout << "That's highly suspicious, no order for P-type but D-type functions? Let me stop before i do something stupid!" << endl;
          return -1;
        }
      }
            break;
      case 4:
        for (int i = 0; i < get_atom_shell_primitives(a, s); i++) {
          if (types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] != 11
            || types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] != 12
            || types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] != 13
            || types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] != 17
            || types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] != 14
            || types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] != 15
            || types[get_shell_start_in_primitives(a, s) + 6 * get_atom_shell_primitives(a, s) + i] != 18
            || types[get_shell_start_in_primitives(a, s) + 7 * get_atom_shell_primitives(a, s) + i] != 19
            || types[get_shell_start_in_primitives(a, s) + 8 * get_atom_shell_primitives(a, s) + i] != 16
            || types[get_shell_start_in_primitives(a, s) + 9 * get_atom_shell_primitives(a, s) + i] != 20) {
            if (types[primcounter] != 11 || types[primcounter + 1] != 12 || types[primcounter + 2] != 13
              || types[primcounter + 3] != 14 || types[primcounter + 4] != 15 || types[primcounter + 5] != 17
              || types[primcounter + 6] != 16 || types[primcounter + 7] != 18 || types[primcounter + 8] != 19 || types[primcounter + 9] != 20) {
              order = -1;
              if (debug) {
                cout << "The checked types are 10 from #" << primcounter << " and are:" << endl;
                cout << types[primcounter] << " " << types[primcounter + 1] << " " << types[primcounter + 2] << " "
                  << types[primcounter + 3] << " " << types[primcounter + 4] << " " << types[primcounter + 5] << " "
                  << types[primcounter + 6] << " " << types[primcounter + 7] << " " << types[primcounter + 8] << " "
                  << types[primcounter + 9] << endl;
                cout << "Something seems to be wrong in the order of your F-Types..." << endl;
                return -1;
              }
            }
            else {
              f_order = 3;
              primcounter += 10;
            }
          }
          else {
            f_order = 1;
            primcounter += 10;
          }
        }
        break;
      default:
        cout << "ERROR in type assignement!!" << endl;
        return -1;
        break;

      }
    }
  }
  /*if(debug){
    cout << "Going to return " << f_order*10+order << endl;
    Enter();
  }*/
  return (f_order * 10 + order);
};

bool WFN::sort_wfn(const int& g_order, const bool& debug)
{
  set_modified();
  int primcounter = 0;
  int f_order = 0;
  int order = g_order;
  //Sorry for this way of forwarding the order, i think 2 switches would have been more nicely
  if (order >= 10 && order < 20) {
    f_order = 1;
    order -= 10;
  }
  else if (order >= 20 && order < 30) {
    f_order = 2;
    order -= 20;
  }
  else if (order >= 30 && order < 40) {
    f_order = 3;
    order -= 30;
  }
  switch (order) {
  case 1:
    for (int a = 0; a < get_ncen(); a++)
      for (int s = 0; s < get_atom_shell_count(a); s++)
        switch (get_shell_type(a, s)) {
        case 1:
          primcounter += get_atom_shell_primitives(a, s);
          break;
        case 2: {
          if (get_atom_shell_primitives(a, s) > 1) {
            vector<int> temp_centers;
            temp_centers.resize(3 * get_atom_shell_primitives(a, s));
            vector<int> temp_types;
            temp_types.resize(3 * get_atom_shell_primitives(a, s));
            vector<double> temp_exponents;
            temp_exponents.resize(3 * get_atom_shell_primitives(a, s));
            vector <vector<double> > temp_MO_coefficients;
            temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 3);
            for (int i = 0; i < get_atom_shell_primitives(a, s) * 3; i++) temp_MO_coefficients[i].resize(nmo);
            for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
              for (int c = 0; c < 3; c++) {
                temp_centers[3 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
                temp_types[3 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
                temp_exponents[3 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
                for (int m = 0; m < nmo; m++)
                  temp_MO_coefficients[3 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
              }
            for (int i = 0; i < 3 * get_atom_shell_primitives(a, s); i++) {
              centers[primcounter + i] = temp_centers[i];
              types[primcounter + i] = temp_types[i];
              exponents[primcounter + i] = temp_exponents[i];
              for (int m = 0; m < nmo; m++)
                err_checkf(MOs[m].change_coefficient(primcounter + i, temp_MO_coefficients[i][m]), "Error while assigning new MO coefficient!", std::cout);
            }
          }
          primcounter += get_atom_shell_primitives(a, s) * 3;
          break;
        }
        case 3: {
          if (get_atom_shell_primitives(a, s) > 1) {
            vector<int> temp_centers;
            temp_centers.resize(6 * get_atom_shell_primitives(a, s));
            vector<int> temp_types;
            temp_types.resize(6 * get_atom_shell_primitives(a, s));
            vector<double> temp_exponents;
            temp_exponents.resize(get_atom_shell_primitives(a, s) * 6);
            vector <vector<double> > temp_MO_coefficients;
            temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 6);
            for (int i = 0; i < get_atom_shell_primitives(a, s) * 6; i++) temp_MO_coefficients[i].resize(nmo);
            for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
              for (int c = 0; c < 6; c++) {
                temp_centers[6 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
                temp_types[6 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
                temp_exponents[6 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
                for (int m = 0; m < nmo; m++)
                  temp_MO_coefficients[6 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
              }
            for (int i = 0; i < 6 * get_atom_shell_primitives(a, s); i++) {
              centers[primcounter + i] = temp_centers[i];
              types[primcounter + i] = temp_types[i];
              exponents[primcounter + i] = temp_exponents[i];
              for (int m = 0; m < nmo; m++)
                if (!MOs[m].change_coefficient(primcounter + i, temp_MO_coefficients[i][m])) {
                  cout << "Error while assigning new MO coefficient!" << endl;
                  return false;
                }
            }
          }
          primcounter += get_atom_shell_primitives(a, s) * 6;
          break;
        }
        case 4: {
          if (get_atom_shell_primitives(a, s) > 1) {
            vector<int> temp_centers;
            temp_centers.resize(10 * get_atom_shell_primitives(a, s));
            vector<int> temp_types;
            temp_types.resize(10 * get_atom_shell_primitives(a, s));
            vector<double> temp_exponents;
            temp_exponents.resize(get_atom_shell_primitives(a, s) * 10);
            vector <vector<double> > temp_MO_coefficients;
            temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 10);
            for (int i = 0; i < get_atom_shell_primitives(a, s) * 10; i++) temp_MO_coefficients[i].resize(nmo);
            for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
              for (int c = 0; c < 10; c++) {
                temp_centers[10 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
                temp_types[10 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
                temp_exponents[10 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
                for (int m = 0; m < nmo; m++)
                  temp_MO_coefficients[10 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
              }
            for (int i = 0; i < 10 * get_atom_shell_primitives(a, s); i++) {
              centers[primcounter + i] = temp_centers[i];
              types[primcounter + i] = temp_types[i];
              exponents[primcounter + i] = temp_exponents[i];
              for (int m = 0; m < nmo; m++)
                if (!MOs[m].change_coefficient(primcounter + i, temp_MO_coefficients[i][m])) {
                  cout << "Error while assigning new MO coefficient!" << endl;
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
    if (debug) {
      cout << "Nothing to be done here, i like tonto type..." << endl;
    }
    break;
  case 3:
    for (int a = 0; a < get_ncen(); a++) {
      for (int s = 0; s < get_atom_shell_count(a); s++) {
        switch (get_shell_type(a, s)) {
        case 1:
          primcounter += get_atom_shell_primitives(a, s);
          break;
        case 2: {
          int temp_center;
          int temp_type;
          double temp_exponent;
          vector<double> temp_MO_coefficients;
          temp_MO_coefficients.resize(nmo);
          for (int i = 0; i < get_atom_shell_primitives(a, s); i++) {
            temp_center = centers[primcounter];
            temp_type = types[primcounter];
            temp_exponent = exponents[primcounter];
            for (int m = 0; m < nmo; m++)
              temp_MO_coefficients[m] = MOs[m].get_coefficient(primcounter);
            for (int j = 0; j < 2; j++) {
              centers[primcounter + j] = centers[primcounter + 1 + j];
              types[primcounter + j] = types[primcounter + 1 + j];
              exponents[primcounter + j] = exponents[primcounter + 1 + j];
              for (int m = 0; m < nmo; m++) {
                err_checkf(MOs[m].change_coefficient(primcounter + j, MOs[m].get_coefficient(primcounter + 1 + j)), "Error while assigning new MO coefficient!", std::cout);
              }
            }
            centers[primcounter + 2] = temp_center;
            types[primcounter + 2] = temp_type;
            exponents[primcounter + 2] = temp_exponent;
            for (int m = 0; m < nmo; m++)
              if (!MOs[m].change_coefficient(primcounter + 2, temp_MO_coefficients[m])) {
                cout << "Error while assigning new MO coefficient!" << endl;
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
    cout << "order type: " << f_order << " " << order << " not supported!" << endl;
    return false;
    break;
  }
  primcounter = 0;
  switch (f_order) {
  case 1:
    for (int a = 0; a < get_ncen(); a++)
      for (int s = 0; s < get_atom_shell_count(a); s++)
        switch (get_shell_type(a, s)) {
        case 1:
          primcounter += get_atom_shell_primitives(a, s);
          break;
        case 2:
          primcounter += get_atom_shell_primitives(a, s) * 3;
          break;
        case 3:
          primcounter += get_atom_shell_primitives(a, s) * 6;
          break;
        case 4: {
          vector <int> temp_center;
          temp_center.resize(10);
          vector <int> temp_type;
          temp_type.resize(10);
          vector <double> temp_exponent;
          temp_exponent.resize(10);
          vector < vector <double> > temp_MO_coefficients;
          temp_MO_coefficients.resize(nmo);
          for (int m = 0; m < nmo; m++)
            temp_MO_coefficients[m].resize(10);
          for (int i = 0; i < get_atom_shell_primitives(a, s); i++) {
            for (int j = 0; j < 10; j++) {
              temp_center[j] = centers[get_shell_start_in_primitives(a, s) + 10 * i + j];
              temp_type[j] = types[get_shell_start_in_primitives(a, s) + 10 * i + j];
              temp_exponent[j] = exponents[get_shell_start_in_primitives(a, s) + 10 * i + j];
              for (int m = 0; m < nmo; m++)
                temp_MO_coefficients[m][j] = MOs[m].get_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + j);
            }
            vector <int> mask{ 0,1,2,6,3,4,7,8,5,9 };
            for (int j = 0; j < 10; j++) {
              centers[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_center[mask[j]];
              types[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_type[mask[j]];
              exponents[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_exponent[mask[j]];
              for (int m = 0; m < nmo; m++)
                MOs[m].change_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + mask[j], temp_MO_coefficients[m][j]);
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
        switch (get_shell_type(a, s)) {
        case 1:
          primcounter += get_atom_shell_primitives(a, s);
          break;
        case 2:
          primcounter += get_atom_shell_primitives(a, s) * 3;
          break;
        case 3:
          primcounter += get_atom_shell_primitives(a, s) * 6;
          break;
        case 4: {
          vector <int> temp_center;
          temp_center.resize(10);
          vector <int> temp_type;
          temp_type.resize(10);
          vector <double> temp_exponent;
          temp_exponent.resize(10);
          vector < vector <double> > temp_MO_coefficients;
          temp_MO_coefficients.resize(nmo);
          for (int m = 0; m < nmo; m++) temp_MO_coefficients[m].resize(10);
          for (int i = 0; i < get_atom_shell_primitives(a, s); i++) {
            for (int j = 0; j < 10; j++) {
              temp_center[j] = centers[get_shell_start_in_primitives(a, s) + 10 * i + j];
              temp_type[j] = types[get_shell_start_in_primitives(a, s) + 10 * i + j];
              temp_exponent[j] = exponents[get_shell_start_in_primitives(a, s) + 10 * i + j];
              for (int m = 0; m < nmo; m++)
                temp_MO_coefficients[m][j] = MOs[m].get_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + j);
            }
            vector <int> mask{ 0,1,2,3,4,6,5,7,8,9 };
            for (int j = 0; j < 10; j++) {
              centers[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_center[mask[j]];
              types[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_type[mask[j]];
              exponents[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_exponent[mask[j]];
              for (int m = 0; m < nmo; m++)
                MOs[m].change_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + mask[j], temp_MO_coefficients[m][j]);
            }
          }
          primcounter += 10;
        }
              break;
        }
    break;
  case 0:
    if (debug) {
      cout << "There seems to be no f-functions!" << endl;
    }
    break;
  default:
    cout << "f-order type: " << f_order << " not supported!" << endl;
    return false;
    break;
  }
  return true;
};

void WFN::set_has_ECPs(const bool& in, const bool& apply_to_atoms) {
  has_ECPs = in;
  if (apply_to_atoms) {
#pragma omp parallel for
    for (int i = 0; i < ncen; i++) {
      if (ECP_electrons[atoms[i].charge] != 0) {
        atoms[i].ECP_electrons = ECP_electrons[atoms[i].charge];
      }
    }
  }
};

void WFN::operator=(const WFN& right)
{
  ncen = right.get_ncen();
  nmo = right.get_nmo();
  origin = right.get_origin();
  nex = right.get_nex();
  multi = right.get_multi();
  charge = right.get_charge();
  for (int i = 0; i < ncen; i++) {
    push_back_center(right.get_center(i));
    push_back_type(right.get_type(i));
    push_back_exponent(right.get_exponent(i));
  }
  for (int i = 0; i < nmo; i++) {
    push_back_MO(i, right.get_MO_primitive_count(i), right.get_MO_energy(i));
    for (int j = 0; j < right.get_MO_primitive_count(i); j++) 
      push_back_MO_coef(i, right.get_MO_coef(i, j));
  }
  for (int c = 0; c < right.cub.size(); c++) {
    cub.push_back(right.cub[c]);
  }
  for (int a = 0; a < right.atoms.size(); a++) {
    atoms.push_back(right.atoms[a]);
  }
};

int WFN::calculate_charge()
{
  int atomic_charges = 0;
  double mo_charges = 0;
  for (int a = 0; a < ncen; a++) {
    int nr = atoms[a].charge;
    if (nr == 0) {
      cout << "ERROR: Atomtype misunderstanding!" << endl;
      return -1000;
    }
    atomic_charges += nr;
  }
  for (int mo = 0; mo < nmo; mo++) {
    mo_charges += get_MO_occ(mo);
  }
  return atomic_charges - mo_charges;
};

int WFN::calculate_charge(ofstream& file)
{
  int atomic_charges = 0;
  double mo_charges = 0;
  for (int a = 0; a < ncen; a++) {
    int nr = atoms[a].charge;
    if (nr == 0) {
      file << "ERROR: Atomtype misunderstanding!" << endl;
      return -1000;
    }
    atomic_charges += nr;
  }
  for (int mo = 0; mo < nmo; mo++) {
    mo_charges += get_MO_occ(mo);
  }
  return atomic_charges - mo_charges;
};

bool WFN::guess_multiplicity(bool expert)
{
  bool f = false;
  if (get_nr_electrons(f) % 2 == 0 && !expert) {
    bool temp = true;
    cout << "With " << get_nr_electrons(f) << " electrons your system appears to have multiplicity 1." << endl;
    assign_multi(1);
  }
  else if (get_nr_electrons(f) % 2 == 0) {
    cout << "With " << get_nr_electrons(f) << " electrons your system appears to have multiplicity 1. Do you agree?" << endl;
    if (yesno()) assign_multi(1);
    else {
      bool end = false;
      while (!end) {
        cout << "What is the multiplicity of your molecule? ";
        int temp = 0;
        cin >> temp;
        if (temp % 2 != 1) {
          cout << "This multiplicity does not match your number of electrons! Try again!" << endl;
          continue;
        }
        else {
          assign_multi(temp);
          end = true;
        }
      }
    }
  }
  else if (get_nr_electrons(f) % 2 == 1 && !expert) {
    cout << "With " << get_nr_electrons(f) % 2 << " electron this seems to be open shell, assuming multiplicity 2." << endl;
    assign_multi(2);
  }
  else if (get_nr_electrons(f) % 2 == 1) {
    cout << "With " << get_nr_electrons(f) % 2 << " electron this seems to be open shell, assuming multiplicity 2. Do you agree?" << endl;
    if (yesno()) assign_multi(2);
    else {
      bool end = false;
      while (!end) {
        cout << "What is the multiplicity of your molecule? ";
        int tmp = 0;
        cin >> tmp;
        if (tmp % 2 != 1) {
          cout << "This Multiplicity and electron count does not match! Try again!" << endl;
          continue;
        }
        else {
          assign_multi(tmp);
          end = true;
        }
      }
    }
  }
  else {
    cout << "This is awkward... i dind't think of this case yet, contact florian to implement it!" << endl;
    return false;
  }
  return true;
};

bool WFN::guess_multiplicity(ofstream& file, bool expert)
{
  bool f = false;
  if (get_nr_electrons(f) % 2 == 0 && !expert) {
    bool temp = true;
    file << "With " << get_nr_electrons(f) << " electrons your system appears to have multiplicity 1." << endl;
    assign_multi(1);
  }
  else if (get_nr_electrons(f) % 2 == 0) {
    file << "With " << get_nr_electrons(f) << " electrons your system appears to have multiplicity 1. Do you agree?" << endl;
    if (yesno()) assign_multi(1);
    else {
      bool end = false;
      while (!end) {
        file << "What is the multiplicity of your molecule? ";
        int temp = 0;
        cin >> temp;
        if (temp % 2 != 1) {
          file << "This multiplicity does not match your number of electrons! Try again!" << endl;
          continue;
        }
        else {
          assign_multi(temp);
          end = true;
        }
      }
    }
  }
  else if (get_nr_electrons(f) % 2 == 1 && !expert) {
    file << "With " << get_nr_electrons(f) % 2 << " electron this seems to be open shell, assuming multiplicity 2." << endl;
    assign_multi(2);
  }
  else if (get_nr_electrons(f) % 2 == 1) {
    file << "With " << get_nr_electrons(f) % 2 << " electron this seems to be open shell, assuming multiplicity 2. Do you agree?" << endl;
    if (yesno()) assign_multi(2);
    else {
      bool end = false;
      while (!end) {
        file << "What is the multiplicity of your molecule? ";
        int tmp = 0;
        cin >> tmp;
        if (tmp % 2 != 1) {
          file << "This Multiplicity and electron count does not match! Try again!" << endl;
          continue;
        }
        else {
          assign_multi(tmp);
          end = true;
        }
      }
    }
  }
  else {
    file << "This is awkward... i dind't think of this case yet, contact florian to implement it!" << endl;
    return false;
  }
  return true;
};

/*
bool WFN::change_center(int nr, int value){
  if(nr>centers.size()||nr<0||value <=0){
    cout << "This is an imposible choice" << endl;
    Enter();
    return false;
  }
  else centers[nr]=value;
};

bool WFN::change_type(int nr, int value){
  if(nr>types.size()||nr<0||value <=0||value>20){
    cout << "This is an imposible choice" << endl;								NOT NEEDED AT THIS POINT
    Enter();
    return false;
  }
  else types[nr]=value;
};

bool WFN::change_exponent(int nr, double value){
  if(nr>exponents.size()||nr<0||value <=0){
    cout << "This is an imposible choice" << endl;
    Enter();
    return false;
  }
  else exponents[nr]=value;
};
*/

bool WFN::push_back_cube(const string& filepath, const bool& full, const bool& expert)
{
  cub.push_back(cube(filepath, full, *this, cout, expert));
  return true;
};

void WFN::pop_back_cube()
{
  cub.pop_back();
}

double* WFN::get_ptr_mo_coefficients(const int& mo)
{
  return MOs[mo].get_ptr_coefficients();
};

unsigned int WFN::get_atom_integer_mass(const unsigned int& atomnr)
{
  if (get_atom_charge(atomnr) > 86) {
    cout << "Sorry, only implemented until Rn yet, ask Florian for increases!" << endl;
    return 0;
  }
  if (get_atom_charge(atomnr) <= 0) {
    cout << "sorry, something seems wrong with the atoms you requested!" << endl;
    return 0;
  }
  return integer_masses[get_atom_charge(atomnr) - 1];
};

double WFN::get_atom_real_mass(const int& atomnr)
{
  if (get_atom_charge(atomnr) > 86) {
    cout << "Sorry, only implemented until Xe yet, ask Florian for increases!" << endl;
    return 0;
  }
  if (get_atom_charge(atomnr) <= 0) {
    cout << "sorry, something seems wrong with the atoms you requested!" << endl;
    return 0;
  }
  return real_masses[get_atom_charge(atomnr) - 1];
}

double WFN::get_MO_occ(const int& nr) const
{
  return MOs[nr].get_occ();
};

void WFN::delete_unoccupied_MOs()
{
  for (int i = MOs.size() - 1; i >= 0; i--) {
    if (get_MO_occ(i) == 0.0) {
      MOs.erase(MOs.begin() + i);
      nmo--;
    }
  }
};

bool WFN::read_fchk(const string& filename, ostream& log, const bool debug)
{
  vector<vector<double>> mat_5d6d, mat_7f10f, mat_9g15g, mat_11h21h;
  if (!generate_cart2sph_mat(mat_5d6d, mat_7f10f, mat_9g15g, mat_11h21h)) log << "Error during geenration of matrix" << endl;
  int r_u_ro_switch = 0;
  ifstream fchk(filename, ios::in);
  if (!fchk.is_open()) {
    log << "ERROR while opening .fchk file!" << endl;
    return false;
  }
  string line;
  getline(fchk, line);
  string title = line;
  getline(fchk, line);
  string calculation_level = line;
  if (line[10] == 'R') {
    if (line[11] == 'O') {
      if (line[12] == '3')
        r_u_ro_switch = 0;
      else
        r_u_ro_switch = 2;
    }
  }
  else if (line[10] == 'U') //Unrestricted
    r_u_ro_switch = 1;
  int el, ael, bel;
  el = read_fchk_integer(fchk, "Number of electrons", false);
  getline(fchk, line);
  ael = read_fchk_integer(line);
  getline(fchk, line);
  bel = read_fchk_integer(line);
  if (ael != bel && r_u_ro_switch == 0) r_u_ro_switch = 1; //If U was not correctly recognized
  if (calculation_level.find("CASSCF") != string::npos && ael != bel) r_u_ro_switch = 2; // CASSCF requires open shell treatment
  int nbas = read_fchk_integer(fchk, "Number of basis functions");
  line = go_get_string(fchk, "Number of independent functions");
  int indbas;
  if (line == "") indbas = read_fchk_integer(fchk, "Number of independant functions");
  else indbas = read_fchk_integer(line);
  line = go_get_string(fchk, "Virial Ratio");
  if (line != "")
    virial_ratio = read_fchk_double(line);
  line = go_get_string(fchk, "Total Energy");
  if (line != "")
    total_energy = read_fchk_double(line);
  vector<int> atnbrs;
  if (!read_fchk_integer_block(fchk, "Atomic numbers", atnbrs)) {
    log << "Error reading atnbrs" << endl;
    return false;
  }
  ncen = atnbrs.size();
  atoms.resize(ncen);
  for (int i = 0; i < ncen; i++)
    atoms[i].label = atnr2letter(atnbrs[i]);
  vector<double> charges;
  if (!read_fchk_double_block(fchk, "Nuclear charges", charges)) {
    log << "Error reading charges" << endl;
    return false;
  }
  for (int i = 0; i < charges.size(); i++)
    atoms[i].charge = charges[i];
  vector<double> coords;
  if (!read_fchk_double_block(fchk, "Current cartesian coordinates", coords)) {
    log << "Error reading coordinates" << endl;
    return false;
  }
  if (coords.size() != ncen * 3) {
    log << "Inconsistant number of atoms and coordinates" << endl;
    return false;
  }
  for (int i = 0; i < ncen; i++) {
    atoms[i].x = coords[3 * i];
    atoms[i].y = coords[3 * i + 1];
    atoms[i].z = coords[3 * i + 2];
  }
  vector <int> shell_types;
  if (!read_fchk_integer_block(fchk, "Shell types", shell_types, false)) {
    log << "Error reading shell types" << endl;
    return false;
  }
  bool is_spherical = false;
  for (int i = 0; i < shell_types.size(); i++) if (shell_types[i] < -1) is_spherical = true;
  if (debug) log << "This fchk contains spherical harmonics, which will be transformed into cartesian functions!" << endl
    << "Loading basis set information..." << endl;
  vector <int> nr_prims_shell;
  if (!read_fchk_integer_block(fchk, "Number of primitives per shell", nr_prims_shell)) {
    log << "Error reading primitives per shell" << endl;
    return false;
  }
  vector <int> shell2atom;
  if (!read_fchk_integer_block(fchk, "Shell to atom map", shell2atom)) {
    log << "Error reading shell2atom" << endl;
    return false;
  }
  vector <double> exp;
  if (!read_fchk_double_block(fchk, "Primitive exponents", exp)) {
    log << "Error reading Primitive exponents" << endl;
    return false;
  }
  vector <double> contraction;
  if (!read_fchk_double_block(fchk, "Contraction coefficients", contraction)) {
    log << "Error reading Contraction coefficients" << endl;
    return false;
  }
  vector<double> acoef, bcoef;
  vector<double> MOocc, aMOene, bMOene;
  int l_nmo;
  if (r_u_ro_switch == 0 || r_u_ro_switch == 2) { // Restricted or Restricted-Open-Shell
    if (!read_fchk_double_block(fchk, "Alpha Orbital Energies", aMOene)) {
      log << "Error during reading of Alpha Energies" << endl;
      return false;
    }
    if (!read_fchk_double_block(fchk, "MO coefficients", acoef)) {
      log << "Error during reading of Alpha MOs" << endl;
      return false;
    }
    MOocc.resize(aMOene.size());
    if (r_u_ro_switch == 0)
#pragma omp parallel for
      for (int i = 0; i < MOocc.size(); i++) {
        if (i < ael) MOocc[i] = 2.0;
        else MOocc[i] = 0.0;
      }
    else
#pragma omp parallel for
      for (int i = 0; i < MOocc.size(); i++) {
        if (i < bel) MOocc[i] = 2.0;
        else if (i < ael) MOocc[i] = 1.0;
        else MOocc[i] = 0.0;
      }
  }
  else { // Unrestricted
    l_nmo = 2 * nbas;
    if (!read_fchk_double_block(fchk, "Alpha Orbital Energies", aMOene)) {
      log << "Error during reading of Alpha Energies" << endl;
      return false;
    }
    if (!read_fchk_double_block(fchk, "Beta Orbital Energies", bMOene)) {
      log << "Error during reading of Beta Energies" << endl;
      return false;
    }
    if (!read_fchk_double_block(fchk, "Alpha MO coefficients", acoef)) {
      log << "Error during reading of Alpha MOs" << endl;
      return false;
    }
    if (!read_fchk_double_block(fchk, "Beta MO coefficients", bcoef)) {
      log << "Error during reading of Beta MOs" << endl;
      return false;
    }
    MOocc.resize(aMOene.size() + bMOene.size());
#pragma omp parallel for
    for (int i = 0; i < aMOene.size(); i++) {
      if (i < ael) MOocc[i] = 1.0;
      else MOocc[i] = 0.0;
    }
#pragma omp parallel for
    for (int i = aMOene.size(); i < aMOene.size() + bMOene.size(); i++) {
      if (i - aMOene.size() < bel) MOocc[i] = 1.0;
      else MOocc[i] = 0.0;
    }
  }
  if (debug) log << "Finished reading the file! Transferring to WFN object!" << endl;
  vector<int> shelltypesspherical;
  int nbas5D;
  nmo = nbas;
  int nshell = shell_types.size();
  if (is_spherical) {
    shelltypesspherical.resize(shell_types.size());
    if (debug) {
      log << "shelltype:" << endl;
      for (int i = 0; i < nshell; i++) log << setw(3) << shell_types[i] << endl;
    }
    shelltypesspherical = shell_types;
    for (int i = 0; i < nshell; i++) if (shell_types[i] <= -2) shell_types[i] = -shell_types[i];
    nbas5D = nbas;
    nbas = 0;
    for (int i = 0; i < nshell; i++) nbas += sht2nbas(shell_types[i]);
  }
  if (debug) {
    log << "sizes" << endl;
    log << setw(3) << nshell << endl;
    log << setw(3) << nbas5D << endl;
    log << setw(3) << nbas << endl;
  }
  vector<int> shelltypescartesian(size_t(shell_types.size()), 0);
  shelltypescartesian = shell_types;
  int nbasCart = nbas;
  int nprims = 0;
  for (int i = 0; i < nshell; i++) nprims += sht2nbas(shell_types[i]) * nr_prims_shell[i];
  if (debug) {
    log << "nprim" << endl;
    log << setw(3) << nprims << endl;
    log << "amocoeff of MO:1";
    for (int i = 0; i < 4; i++) {
      if (i % 2 == 0)
        log << endl;
      log << setprecision(8) << setw(16) << scientific << acoef[i];
    }
    log << endl << "First of MOs 1-4:" << endl;
    int temp = 0;
    for (int i = 0; i < 4; i++) {
      if (i % 2 == 0)
        log << endl;
      if (is_spherical)
        temp = i * nbas5D;
      else
        temp = i * nbas;
      log << setprecision(8) << setw(16) << scientific << acoef[temp];
    }
  }
  //vector<primitive> basis;
  vector<vector<double>> COa, COb, CObasa, CObasb, CObasa_spherical, CObasb_spherical;
  //vector<int> basshell, bascen, bastype, basstart, basend, primstart, primend;
  vector<double> primconnorm;
  //create arrays
  //basshell.resize(nbas);
  //bascen.resize(nbas);
  //bastype.resize(nbas);
  //primstart.resize(nbas);
  //primend.resize(nbas);
  primconnorm.resize(nprims);
  //basstart.resize(ncen);
  //basend.resize(nbas);
  exponents.resize(nprims);
  centers.resize(nprims);
  types.resize(nprims);
  COa.resize(nmo);
#pragma omp parallel for
  for (int i = 0; i < nmo; i++)
    COa[i].resize(nprims);
  CObasa.resize(nbas);
#pragma omp parallel for
  for (int i = 0; i < nbas; i++)
    CObasa[i].resize(nbas);
  if (r_u_ro_switch == 1) {
    COb.resize(nmo);
#pragma omp parallel for
    for (int i = 0; i < nmo; i++)
      COb[i].resize(nprims);
    CObasb.resize(nbas);
#pragma omp parallel for
    for (int i = 0; i < nbas; i++)
      CObasb[i].resize(nbas);
  }

  if (is_spherical) {
    // NEEEEEEEEDS TO BE DONE!!!!!!!!
    CObasa_spherical.resize(nbas5D);
#pragma omp parallel for
    for (int mo = 0; mo < nbas5D; mo++)
      CObasa_spherical[mo].resize(nbas);
#pragma omp parallel for
    for (int mo = 0; mo < nbas5D; mo++)
      for (int b = 0; b < nbas5D; b++)
        CObasa_spherical[mo][b] = acoef[nbas5D * b + mo];
    if (debug) {
      log << endl << "CObasa5d";
      int run = 0;
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < nbas5D; j++) {
          if (run % 2 == 0) log << endl;
          log << setprecision(8) << setw(16) << scientific << CObasa_spherical[i][j];
          run++;
        }
        log << endl;
      }
      log << endl;
    }
    if (r_u_ro_switch == 1) {
      CObasb_spherical.resize(nbas5D);
#pragma omp parallel for
      for (int mo = 0; mo < nbas5D; mo++)
        CObasb_spherical[mo].resize(nbas5D);
#pragma omp parallel for
      for (int mo = 0; mo < nbas5D; mo++)
        for (int b = 0; b < nbas5D; b++)
          CObasb_spherical[mo][b] = bcoef[nbas5D * b + mo];
    }
    int ipos_spher = 0, ipos_cart = 0;
    for (int shell = 0; shell < nshell; shell++) {
      int temp_typ5D = shelltypesspherical[shell];
      int temp_typ6D = shell_types[shell];
      int shell_size5D = sht2nbas(temp_typ5D);
      int shell_size6D = sht2nbas(temp_typ6D);
      if (debug) {
        log << setw(3) << ipos_spher;
        log << setw(3) << ipos_cart;
        log << setw(3) << temp_typ5D;
        log << setw(3) << temp_typ6D;
        log << setw(3) << shell_size5D;
        log << setw(3) << shell_size6D;
        log << endl;
      }
      if (temp_typ5D >= -1) {// S and P shells are fine!
#pragma omp parallel for
        for (int i = 0; i < nbas5D; i++)
          for (int j = 0; j < shell_size6D; j++)
            CObasa[ipos_cart + j][i] = CObasa_spherical[ipos_spher + j][i];
        if (debug) {
          int run = 0;
          for (int i = 0; i < nbas5D; i++) {
            if (run % 2 == 0) log << endl;
            log << setprecision(8) << setw(16) << scientific << CObasa_spherical[ipos_spher][i];
            run++;
          }
          log << endl;
          for (int i = 0; i < nbas; i++) {
            if (run % 2 == 0) log << endl;
            log << setprecision(8) << setw(16) << scientific << CObasa[ipos_cart][i];
            run++;
          }
          log << endl;
        }
        //if (debug) {
        //	int run = 0;
        //	for (int i = 0; i < nbas5D; i++) {
        //		if (run % 2 == 0) log << endl;
        //		log << setprecision(8) << setw(16) << scientific << CObasa_spherical[ipos_spher][i];
        //		run++;
        //	}
        //	log << endl;
        //}
      }
      else if (temp_typ5D == -2) {// 5D -> 6D
#pragma omp parallel for
        for (int i = 0; i < nbas5D; i++)
          for (int j = 0; j < 6; j++)
            CObasa[ipos_cart + j][i] =
            mat_5d6d[j][0] * CObasa_spherical[ipos_spher][i]
            + mat_5d6d[j][1] * CObasa_spherical[ipos_spher + 1][i]
            + mat_5d6d[j][2] * CObasa_spherical[ipos_spher + 2][i]
            + mat_5d6d[j][3] * CObasa_spherical[ipos_spher + 3][i]
            + mat_5d6d[j][4] * CObasa_spherical[ipos_spher + 4][i];
        if (debug) {
          int run = 0;
          for (int j = 0; j < 5; j++)
            for (int i = 0; i < nbas5D; i++) {
              if (run % 2 == 0) log << endl;
              log << setprecision(8) << setw(16) << scientific << CObasa_spherical[ipos_spher + j][i];
              run++;
            }
          log << endl;
        }
        if (debug) {
          int run = 0;
          for (int j = 0; j < 6; j++) {
            for (int i = 0; i < nbas5D; i++) {
              if (run % 2 == 0) log << endl;
              log << setprecision(8) << setw(16) << scientific << CObasa[ipos_cart + j][i];
              run++;
            }
            log << endl;
          }
          log << endl;
        }
      }
      else if (temp_typ5D == -3) // 7F -> 10F
#pragma omp parallel for
        for (int i = 0; i < nbas5D; i++)
          for (int j = 0; j < 10; j++)
            CObasa[ipos_cart + j][i] =
            mat_7f10f[j][0] * CObasa_spherical[ipos_spher][i]
            + mat_7f10f[j][1] * CObasa_spherical[ipos_spher + 1][i]
            + mat_7f10f[j][2] * CObasa_spherical[ipos_spher + 2][i]
            + mat_7f10f[j][3] * CObasa_spherical[ipos_spher + 3][i]
            + mat_7f10f[j][4] * CObasa_spherical[ipos_spher + 4][i]
            + mat_7f10f[j][5] * CObasa_spherical[ipos_spher + 5][i]
            + mat_7f10f[j][6] * CObasa_spherical[ipos_spher + 6][i];
      else if (temp_typ5D == -4) // 9G -> 15G
#pragma omp parallel for
        for (int i = 0; i < nbas5D; i++)
          for (int j = 0; j < 15; j++)
            CObasa[ipos_cart + j][i] =
            mat_9g15g[j][0] * CObasa_spherical[ipos_spher][i]
            + mat_9g15g[j][1] * CObasa_spherical[ipos_spher + 1][i]
            + mat_9g15g[j][2] * CObasa_spherical[ipos_spher + 2][i]
            + mat_9g15g[j][3] * CObasa_spherical[ipos_spher + 3][i]
            + mat_9g15g[j][4] * CObasa_spherical[ipos_spher + 4][i]
            + mat_9g15g[j][5] * CObasa_spherical[ipos_spher + 5][i]
            + mat_9g15g[j][6] * CObasa_spherical[ipos_spher + 6][i]
            + mat_9g15g[j][7] * CObasa_spherical[ipos_spher + 7][i]
            + mat_9g15g[j][8] * CObasa_spherical[ipos_spher + 8][i];
      else if (temp_typ5D == -5) // 11H -> 21H
#pragma omp parallel for
        for (int i = 0; i < nbas5D; i++)
          for (int j = 0; j < 21; j++)
            CObasa[ipos_cart + j][i] =
            mat_11h21h[j][0] * CObasa_spherical[ipos_spher][i]
            + mat_11h21h[j][1] * CObasa_spherical[ipos_spher + 1][i]
            + mat_11h21h[j][2] * CObasa_spherical[ipos_spher + 2][i]
            + mat_11h21h[j][3] * CObasa_spherical[ipos_spher + 3][i]
            + mat_11h21h[j][4] * CObasa_spherical[ipos_spher + 4][i]
            + mat_11h21h[j][5] * CObasa_spherical[ipos_spher + 5][i]
            + mat_11h21h[j][6] * CObasa_spherical[ipos_spher + 6][i]
            + mat_11h21h[j][7] * CObasa_spherical[ipos_spher + 7][i]
            + mat_11h21h[j][8] * CObasa_spherical[ipos_spher + 8][i]
            + mat_11h21h[j][9] * CObasa_spherical[ipos_spher + 9][i]
            + mat_11h21h[j][10] * CObasa_spherical[ipos_spher + 10][i];
      if (r_u_ro_switch == 1) {
        if (temp_typ5D >= -1) // S and P shells are fine!
#pragma omp parallel for
          for (int i = 0; i < nbas5D; i++)
            for (int j = 0; j < shell_size6D; j++)
              CObasb[ipos_cart + j][i] = CObasb_spherical[ipos_spher + j][i];
        else if (temp_typ5D == -2) // 5D -> 6D
#pragma omp parallel for
          for (int i = 0; i < nbas5D; i++)
            for (int j = 0; j < 6; j++)
              CObasb[ipos_cart + j][i] =
              mat_5d6d[j][0] * CObasb_spherical[ipos_spher][i]
              + mat_5d6d[j][1] * CObasb_spherical[ipos_spher + 1][i]
              + mat_5d6d[j][2] * CObasb_spherical[ipos_spher + 2][i]
              + mat_5d6d[j][3] * CObasb_spherical[ipos_spher + 3][i]
              + mat_5d6d[j][4] * CObasb_spherical[ipos_spher + 4][i];
        else if (temp_typ5D == -3) // 7F -> 10F
#pragma omp parallel for
          for (int i = 0; i < nbas5D; i++)
            for (int j = 0; j < 10; j++)
              CObasb[ipos_cart + j][i] =
              mat_7f10f[j][0] * CObasb_spherical[ipos_spher][i]
              + mat_7f10f[j][1] * CObasb_spherical[ipos_spher + 1][i]
              + mat_7f10f[j][2] * CObasb_spherical[ipos_spher + 2][i]
              + mat_7f10f[j][3] * CObasb_spherical[ipos_spher + 3][i]
              + mat_7f10f[j][4] * CObasb_spherical[ipos_spher + 4][i]
              + mat_7f10f[j][5] * CObasb_spherical[ipos_spher + 5][i]
              + mat_7f10f[j][6] * CObasb_spherical[ipos_spher + 6][i];
        else if (temp_typ5D == -4) // 9G -> 15G
#pragma omp parallel for
          for (int i = 0; i < nbas5D; i++)
            for (int j = 0; j < 15; j++)
              CObasb[ipos_cart + j][i] =
              mat_9g15g[j][0] * CObasb_spherical[ipos_spher][i]
              + mat_9g15g[j][1] * CObasb_spherical[ipos_spher + 1][i]
              + mat_9g15g[j][2] * CObasb_spherical[ipos_spher + 2][i]
              + mat_9g15g[j][3] * CObasb_spherical[ipos_spher + 3][i]
              + mat_9g15g[j][4] * CObasb_spherical[ipos_spher + 4][i]
              + mat_9g15g[j][5] * CObasb_spherical[ipos_spher + 5][i]
              + mat_9g15g[j][6] * CObasb_spherical[ipos_spher + 6][i]
              + mat_9g15g[j][7] * CObasb_spherical[ipos_spher + 7][i]
              + mat_9g15g[j][8] * CObasb_spherical[ipos_spher + 8][i];
        else if (temp_typ5D == -5) // 11H -> 21H
#pragma omp parallel for
          for (int i = 0; i < nbas5D; i++)
            for (int j = 0; j < 21; j++)
              CObasb[ipos_cart + j][i] =
              mat_11h21h[j][0] * CObasb_spherical[ipos_spher][i]
              + mat_11h21h[j][1] * CObasb_spherical[ipos_spher + 1][i]
              + mat_11h21h[j][2] * CObasb_spherical[ipos_spher + 2][i]
              + mat_11h21h[j][3] * CObasb_spherical[ipos_spher + 3][i]
              + mat_11h21h[j][4] * CObasb_spherical[ipos_spher + 4][i]
              + mat_11h21h[j][5] * CObasb_spherical[ipos_spher + 5][i]
              + mat_11h21h[j][6] * CObasb_spherical[ipos_spher + 6][i]
              + mat_11h21h[j][7] * CObasb_spherical[ipos_spher + 7][i]
              + mat_11h21h[j][8] * CObasb_spherical[ipos_spher + 8][i]
              + mat_11h21h[j][9] * CObasb_spherical[ipos_spher + 9][i]
              + mat_11h21h[j][10] * CObasb_spherical[ipos_spher + 10][i];
      }
      ipos_cart += shell_size6D;
      ipos_spher += shell_size5D;
    }
    if (debug) {
      log << endl << "CObasa";
      int run = 0;
      for (int i = 0; i < 10; i++) {
        run = 0;
        for (int j = 0; j < nbas; j++) {
          if (run % 2 == 0)
            log << endl;
          log << setprecision(8) << setw(16) << scientific << CObasa[i][j];
          run++;
        }
        log << endl;
      }
      log << endl;
    }
  }
  else {
#pragma omp parallel for
    for (int mo = 0; mo < nbas; mo++)
      for (int b = 0; b < nbas; b++)
        CObasa[b][mo] = acoef[nbas * b + mo];
    if (r_u_ro_switch == 1)
#pragma omp parallel for
      for (int mo = 0; mo < nbas; mo++)
        for (int b = 0; b < nbas; b++)
          CObasb[b][mo] = bcoef[nbas * b + mo];
  }
  if (debug) {
    log << endl;
    int run = 0;
    for (int i = 0; i < contraction.size(); i++) {
      if (run % 2 == 0) log << endl;
      log << setprecision(8) << setw(16) << scientific << contraction[i];
      run++;
    }
    log << flush;
  }
  int k = 0, iexp = 0, ibasis = 0;
  //double tnormgau;
  for (int i = 0; i < nshell; i++) {
    int j;
    for (j = 0; j < nr_prims_shell[i] * sht2nbas(shell_types[i]); j++)
      centers[k + j] = shell2atom[i];
    int lim = sht2nbas(shell_types[i]);
    for (j = 0; j < lim; j++) {
      int temp = shell2function(shell_types[i], j);
#pragma omp parallel for
      for (int l = 0; l < nr_prims_shell[i]; l++)
        types[k + l] = temp;
      for (int l = 0; l < nr_prims_shell[i]; l++) {
        exponents[k] = exp[iexp + l];
        primconnorm[k] = contraction[iexp + l] * normgauss(types[k], exponents[k]);
        if (debug) log << setw(22) << setprecision(16) << normgauss(types[k], exponents[k]) << endl;
#pragma omp parallel for
        for (int mo = 0; mo < nmo; mo++) {
          COa[mo][k] = CObasa[ibasis][mo] * primconnorm[k];
          if (r_u_ro_switch == 1)//R or RO
            COb[mo][k] = CObasb[ibasis][mo] * primconnorm[k];
        }
        k++;
      }
      ibasis++;
    }
    iexp += nr_prims_shell[i];
  }
  nex = nprims;
  if (r_u_ro_switch != 1)
    MOs.resize(aMOene.size());
  else
    MOs.resize(aMOene.size() + bMOene.size());
#pragma omp parallel for
  for (int i = 0; i < MOs.size(); i++) {
    MOs[i].set_ener(aMOene[i]);
    MOs[i].set_nr(i + 1);
    MOs[i].set_occ(MOocc[i]);
    MOs[i].assign_coefs(COa[i]);
  }
  if (r_u_ro_switch == 1) {
#pragma omp parallel for
    for (int i = aMOene.size(); i < aMOene.size() + bMOene.size(); i++) {
      const int b_step = i - aMOene.size();
      MOs[i].set_ener(bMOene[b_step]);
      MOs[i].set_nr(i + 1);
      MOs[i].set_occ(MOocc[i]);
      MOs[i].assign_coefs(COb[b_step]);
    }
  }
  return true;
};

double WFN::compute_dens(
  const double& Pos1,
  const double& Pos2,
  const double& Pos3,
  vector<vector<double>>& d,
  vector<double>& phi,
  const bool& add_ECP_dens
)
{
  if (d_f_switch) {
    err_checkf(d.size() >= 5, "d is too small!", std::cout);
    err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
    return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi, add_ECP_dens);
  }
  else {
    err_checkf(d.size() >= 4, "d is too small!", std::cout);
    err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
    return compute_dens_cartesian(Pos1, Pos2, Pos3, d, phi, add_ECP_dens);
  }
};

double WFN::compute_dens(
  const double& Pos1,
  const double& Pos2,
  const double& Pos3,
  const bool& add_ECP_dens
)
{
  vector<vector<double>> d;
  vector<double> phi(nmo, 0.0);
  const int n = get_nmo(true);
  if (d_f_switch) {
    d.resize(5);
    for (int i = 0; i < 5; i++)
      d[i].resize(ncen, 0.0);
    err_not_impl_f("Nah.. not yet implemented correctly", std::cout);
    return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi, add_ECP_dens);
  }
  else {
    d.resize(4);
    for (int i = 0; i < 4; i++)
      d[i].resize(ncen, 0.0);
    return compute_dens_cartesian(Pos1, Pos2, Pos3, d, phi, add_ECP_dens);
  }
};

double WFN::compute_dens_cartesian(
  const double& Pos1,
  const double& Pos2,
  const double& Pos3,
  vector<vector<double>>& d,
  vector<double>& phi,
  const bool& add_ECP_dens
)
{
  std::fill(phi.begin(), phi.end(), 0.0);
  double Rho = 0.0;
  int iat;
  int l[3];
  double ex;
  int mo;

  for (iat = 0; iat < ncen; iat++) {
    d[0][iat] = Pos1 - atoms[iat].x;
    d[1][iat] = Pos2 - atoms[iat].y;
    d[2][iat] = Pos3 - atoms[iat].z;
    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
    if (add_ECP_dens && has_ECPs && atoms[iat].ECP_electrons != 0) {  // This adds a tight core density based on 
      Rho += 8 * atoms[iat].ECP_electrons * exp(-4 * PI * d[3][iat]); // a spherical gaussian to fill in the gap 
    }                                                                 // left by ECPs integrating to the correct number of electrons
  }

  for (int j = 0; j < nex; j++) {
    iat = centers[j] - 1;
    //if (iat != atom) continue;
    type2vector(types[j], l);
    ex = -exponents[j] * d[3][iat];
    if (ex < -46.0517) { //corresponds to cutoff of ex ~< 1E-20
      continue;
    }
    ex = exp(ex);
    for (int k = 0; k < 3; k++) {
      if (l[k] == 0)		continue;
      else if (l[k] == 1)	ex *= d[k][iat];
      else if (l[k] == 2)	ex *= d[k][iat] * d[k][iat];
      else if (l[k] == 3)	ex *= pow(d[k][iat], 3);
      else if (l[k] == 4)	ex *= pow(d[k][iat], 4);
      else if (l[k] == 5)	ex *= pow(d[k][iat], 5);
    }
    auto run = phi.data();
    auto run2 = MOs.data();
    for (mo = 0; mo < nmo; mo++) {
      *run += (*run2).get_coefficient_f(j) * ex;      //build MO values at this point
      run2++, run++;
    }
  }

  auto run = phi.data();
  auto run2 = MOs.data();
  for (mo = 0; mo < nmo; mo++) {
    Rho += (*run2).get_occ() * pow(*run, 2);
    run++, run2++;
  }

  return Rho;
}

double WFN::compute_MO_spherical(
  const double& Pos1,
  const double& Pos2,
  const double& Pos3,
  const int& MO
)
{
  err_not_impl_f("This one is not tested an will most likely not work, therefore aborting!", cout);
  err_checkf(d_f_switch, "Only works for spheriacl wavefunctions!", std::cout);
  int iat;
  int l;
  //ex will carry information about radial function
  double ex;
  vector<vector<double>> d(5);
  for (int i = 0; i < 5; i++)
    d[i].resize(ncen);
  double phi(0.0);

  for (iat = 0; iat < ncen; iat++) {
    d[0][iat] = Pos1 - atoms[iat].x;
    d[1][iat] = Pos2 - atoms[iat].y;
    d[2][iat] = Pos3 - atoms[iat].z;
    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
    d[4][iat] = sqrt(d[3][iat]);
  }
  /*Here d[0] = x
         d[1] = y
         d[2] = z
         d[3] = r^2
         d[4] = r
         */

  for (int j = 0; j < nex; j++) {
    iat = centers[j] - 1;
    ex = -exponents[j] * d[3][iat];
    if (ex < -46.0517) { //corresponds to cutoff of ex ~< 1E-20
      continue;
    }
    //apply radial function:
    ex = exp(ex);
    l = types[j];
    if (l > 1) {
      if (l <= 4)	ex *= d[4][iat];
      else if (l <= 9)	ex *= d[3][iat];
      else if (l <= 16)	ex *= d[3][iat] * d[4][iat];
      else if (l <= 25)	ex *= d[3][iat] * d[3][iat];
      else if (l <= 36)	ex *= pow(d[4][iat], 5);
    }
    //calc spherical harmonic
    double SH = 0;
    switch (l) {
    case 1: { //S
      SH = c_1_4p; break;
    }
    case 4: { //P 0 Z
      SH = c_3_4p * d[2][iat] / d[4][iat]; break;
    }
    case 2: { //P 1 X
      SH = c_3_4p * d[0][iat] / d[4][iat]; break;
    }
    case 3: { //P -1 Y
      SH = c_3_4p * d[1][iat] / d[4][iat]; break;
    }
    case 5: { //D 0 Z2
      SH = c_5_16p / d[3][iat] * (3 * pow(d[2][iat], 2) - d[3][iat]); break;
    }
    case 6: { //D 1 XZ
      SH = c_15_4p * d[0][iat] * d[2][iat] / d[3][iat]; break;
    }
    case 7: { //D -1 YZ
      SH = c_15_4p * d[1][iat] * d[2][iat] / d[3][iat]; break;
    }
    case 8: { //D 2 X2-Y2
      SH = c_15_16p * (pow(d[0][iat], 2) - pow(d[1][iat], 2)) / d[3][iat]; break;
    }
    case 9: { //D -2 XY
      SH = c_15_4p * d[1][iat] * d[0][iat] / d[3][iat]; break;
    }
    case 10: { //F 0 Z3
      SH = c_7_16p / pow(d[4][iat], 3) * (5 * pow(d[2][iat], 3) - 3 * d[3][iat] * d[2][iat]); break;
    }
    case 11: { //F 1 XZZ
      SH = c_21_32p / pow(d[4][iat], 3) * d[0][iat] * (5 * pow(d[2][iat], 2) - d[3][iat]); break;
    }
    case 12: { //F -1 YZZ
      SH = c_21_32p / pow(d[4][iat], 3) * d[1][iat] * (5 * pow(d[2][iat], 2) - d[3][iat]); break;
    }
    case 13: { //F 2 Z(X2-Y2)
      SH = c_105_16p / pow(d[4][iat], 3) * ((pow(d[0][iat], 2) - pow(d[1][iat], 2)) * d[2][iat]); break;
    }
    case 14: { //F -2 XYZ
      SH = c_105_4p / pow(d[4][iat], 3) * d[0][iat] * d[1][iat] * d[2][iat]; break;
    }
    case 15: { //F 3 X(X^2-3Y^2)
      SH = c_35_32p / pow(d[4][iat], 3) * d[0][iat] * (pow(d[0][iat], 2) - 3 * pow(d[1][iat], 2)); break;
    }
    case 16: { //F -3 Y(3X^2-Y^2)
      SH = c_35_32p / pow(d[4][iat], 3) * d[1][iat] * (3 * pow(d[0][iat], 2) - pow(d[1][iat], 2)); break;
    }
    case 17: { //G 0 Z^4
      SH = c_9_256p / pow(d[4][iat], 4) * (35 * pow(d[2][iat], 4) - 30 * pow(d[2][iat] * d[4][iat], 2) + 3.0 * pow(d[4][iat], 4)); break;
    }
    case 18: { //G 1 X(7Z^3-3ZR^2)
      SH = c_45_32p / pow(d[4][iat], 4) * d[0][iat] * (7 * pow(d[2][iat], 3) - 3 * d[2][iat] * d[3][iat]); break;
    }
    case 19: { //G -1 Y(7Z^3-3ZR^2)
      SH = c_45_32p / pow(d[4][iat], 4) * d[1][iat] * (7 * pow(d[2][iat], 3) - 3 * d[2][iat] * d[3][iat]); break;
    }
    case 20: { //G 2
      SH = c_45_64p / pow(d[4][iat], 4) * (pow(d[0][iat], 2) - pow(d[1][iat], 2)) * (7 * pow(d[2][iat], 2) - d[3][iat]); break;
    }
    case 21: { //G -2
      SH = c_45_64p / pow(d[4][iat], 4) * d[0][iat] * d[1][iat] * (7 * pow(d[2][iat], 2) - d[3][iat]); break;
    }
    case 22: { //G 3 XZ(X^2-3Y^2)
      SH = c_315_32p / pow(d[4][iat], 4) * d[0][iat] * (pow(d[0][iat], 2) - 3 * pow(d[1][iat], 2)) * d[2][iat]; break;
    }
    case 23: { //G -3 XZ(3X^2-Y^2)
      SH = c_315_32p / pow(d[4][iat], 4) * d[1][iat] * (3 * pow(d[0][iat], 2) - pow(d[1][iat], 2)) * d[2][iat]; break;
    }
    case 24: { //G 4 X^2(X^-3Y^2)-Y^2(3X^2-Y^2)
      SH = c_315_256p / pow(d[4][iat], 4) * ((pow(d[0][iat], 2) * (pow(d[0][iat], 2) - 3 * pow(d[1][iat], 2))) -
        (pow(d[1][iat], 2) * (3 * pow(d[0][iat], 2) - pow(d[1][iat], 2)))); break;
    }
    case 25: { //G -4 XY(X^2-Y^2)
      SH = c_315_16p / pow(d[4][iat], 4) * d[0][iat] * d[1][iat] * (3 * pow(d[0][iat], 2) - pow(d[1][iat], 2)); break;
    }
    default: {
      err_not_impl_f("Higher than G types", std::cout);
    }
    };
    SH *= ex; // multiply radial part with spherical harmonic
    phi += MOs[MO].get_coefficient_f(j) * SH;      //build MO values at this point
  }
  shrink_vector<vector<double>>(d);

  return phi;
}

double WFN::compute_dens_spherical(
  const double& Pos1,
  const double& Pos2,
  const double& Pos3,
  vector<vector<double>>& d,
  vector<double>& phi,
  const bool& add_ECP_dens
)
{
  err_not_impl_f("This one is not tested an will most likely not work, therefore aborting!", cout);
  err_checkf(d_f_switch, "Only works for spheriacl wavefunctions!", std::cout);
  std::fill(phi.begin(), phi.end(), 0.0);
  double Rho = 0.0;
  int iat;
  int l;
  //ex will carry information about radial function
  double ex;
  int mo;

  for (iat = 0; iat < ncen; iat++) {
    d[0][iat] = Pos1 - atoms[iat].x;
    d[1][iat] = Pos2 - atoms[iat].y;
    d[2][iat] = Pos3 - atoms[iat].z;
    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
    if (add_ECP_dens && has_ECPs && atoms[iat].ECP_electrons != 0) {
      Rho += 8 * atoms[iat].ECP_electrons * exp(-4 * PI * d[3][iat]);
    }
    d[4][iat] = sqrt(d[3][iat]);
  }
  /*Here d[0] = x
         d[1] = y
         d[2] = z
         d[3] = r^2
         d[4] = r
         */

  for (int j = 0; j < nex; j++) {
    iat = centers[j] - 1;
    ex = -exponents[j] * d[3][iat];
    if (ex < -46.0517) { //corresponds to cutoff of ex ~< 1E-20
      continue;
    }
    ex = exp(ex);
    //apply radial function:
    l = types[j];
    if (l > 1) {
      if (l <= 4)	ex *= d[4][iat];
      else if (l <= 9)	ex *= d[3][iat];
      else if (l <= 16)	ex *= d[3][iat] * d[4][iat];
      else if (l <= 25)	ex *= d[3][iat] * d[3][iat];
      else if (l <= 36)	ex *= pow(d[4][iat], 5);
    }
    //calc spherical harmonic
    double SH = 0;
    switch (l) {
    case 1: { //S
      SH = c_1_4p; break;
    }
    case 4: { //P 0 Z
      SH = c_3_4p * d[2][iat] / d[4][iat]; break;
    }
    case 2: { //P 1 X
      SH = c_3_4p * d[0][iat] / d[4][iat]; break;
    }
    case 3: { //P -1 Y
      SH = c_3_4p * d[1][iat] / d[4][iat]; break;
    }
    case 5: { //D 0 Z2
      SH = c_5_16p * (3 / d[3][iat] * pow(d[2][iat], 2) - 1.0); break;
    }
    case 6: { //D 1 XZ
      SH = c_15_4p * d[0][iat] * d[2][iat] / d[3][iat]; break;
    }
    case 7: { //D -1 YZ
      SH = c_15_4p * d[1][iat] * d[2][iat] / d[3][iat]; break;
    }
    case 8: { //D 2 X2-Y2
      SH = c_15_16p * (pow(d[0][iat], 2) - pow(d[1][iat], 2)) / d[3][iat]; break;
    }
    case 9: { //D -2 XY
      SH = c_15_4p * d[1][iat] * d[0][iat] / d[3][iat]; break;
    }
    case 10: { //F 0 Z3
      SH = c_7_16p / pow(d[4][iat], 3) * (5 * pow(d[2][iat], 3) - 3 * d[3][iat] * d[2][iat]); break;
    }
    case 11: { //F 1 XZZ
      SH = c_21_32p / pow(d[4][iat], 3) * d[0][iat] * (5 * pow(d[2][iat], 2) - d[3][iat]); break;
    }
    case 12: { //F -1 YZZ
      SH = c_21_32p / pow(d[4][iat], 3) * d[1][iat] * (5 * pow(d[2][iat], 2) - d[3][iat]); break;
    }
    case 13: { //F 2 Z(X2-Y2)
      SH = c_105_16p / pow(d[4][iat], 3) * ((pow(d[0][iat], 2) - pow(d[1][iat], 2)) * d[2][iat]); break;
    }
    case 14: { //F -2 XYZ
      SH = c_105_4p / pow(d[4][iat], 3) * d[0][iat] * d[1][iat] * d[2][iat]; break;
    }
    case 15: { //F 3 X(X^2-3Y^2)
      SH = c_35_32p / pow(d[4][iat], 3) * d[0][iat] * (pow(d[0][iat], 2) - 3 * pow(d[1][iat], 2)); break;
    }
    case 16: { //F -3 Y(3X^2-Y^2)
      SH = c_35_32p / pow(d[4][iat], 3) * d[1][iat] * (3 * pow(d[0][iat], 2) - pow(d[1][iat], 2)); break;
    }
    case 17: { //G 0 Z^4
      SH = c_9_256p / pow(d[4][iat], 4) * (35 * pow(d[2][iat], 4) - 30 * pow(d[2][iat] * d[4][iat], 2) + 3.0 * pow(d[4][iat], 4)); break;
    }
    case 18: { //G 1 X(7Z^3-3ZR^2)
      SH = c_45_32p / pow(d[4][iat], 4) * d[0][iat] * (7 * pow(d[2][iat], 3) - 3 * d[2][iat] * d[3][iat]); break;
    }
    case 19: { //G -1 Y(7Z^2-3ZR^2)
      SH = c_45_32p / pow(d[4][iat], 4) * d[1][iat] * (7 * pow(d[2][iat], 3) - 3 * d[2][iat] * d[3][iat]); break;
    }
    case 20: { //G 2
      SH = c_45_64p / pow(d[4][iat], 4) * (pow(d[0][iat], 2) - pow(d[1][iat], 2)) * (7 * pow(d[2][iat], 2) - d[3][iat]); break;
    }
    case 21: { //G -2
      SH = c_45_64p / pow(d[4][iat], 4) * d[0][iat] * d[1][iat] * (7 * pow(d[2][iat], 2) - d[3][iat]); break;
    }
    case 22: { //G 3 XZ(X^2-3Y^2)
      SH = c_315_32p / pow(d[4][iat], 4) * d[0][iat] * (pow(d[0][iat], 2) - 3 * pow(d[1][iat], 2)) * d[2][iat]; break;
    }
    case 23: { //G -3 XZ(3X^2-Y^2)
      SH = c_315_32p / pow(d[4][iat], 4) * d[1][iat] * (3 * pow(d[0][iat], 2) - pow(d[1][iat], 2)) * d[2][iat]; break;
    }
    case 24: { //G 4 X^2(X^-3Y^2)-Y^2(3X^2-Y^2)
      SH = c_315_256p / pow(d[4][iat], 4) * ((pow(d[0][iat], 2) * (pow(d[0][iat], 2) - 3 * pow(d[1][iat], 2))) -
        (pow(d[1][iat], 2) * (3 * pow(d[0][iat], 2) - pow(d[1][iat], 2)))); break;
    }
    case 25: { //G -4 XY(X^2-Y^2)
      SH = c_315_16p / pow(d[4][iat], 4) * d[0][iat] * d[1][iat] * (3 * pow(d[0][iat], 2) - pow(d[1][iat], 2)); break;
    }
    default: {
      err_checkf(false, "Not yet implemented!", std::cout);
    }
    };
    SH *= ex; // multiply radial part with spherical harmonic
    auto run = phi.data();
    auto run2 = MOs.data();
    for (mo = 0; mo < phi.size(); mo++) {
      *run += (*run2).get_coefficient_f(j) * SH;      //build MO values at this point
      run++, run2++;
    }
  }

  auto run = phi.data();
  auto run2 = MOs.data();
  for (mo = 0; mo < phi.size(); mo++) {
    Rho += (*run2).get_occ() * pow(*run, 2);
    run++, run2++;
  }

  return Rho;
}

void WFN::computeValues(
  const double* PosGrid,		// [3] vector with current position on te grid
  double& Rho,			// Value of Electron Density
  double& normGrad,		// Gradiant Vector
  double* Hess,			// Hessian Matrix, later used to determine lambda2
  double& Elf,			// Value of the ELF
  double& Eli,			// Value of the ELI
  double& Lap,			// Value for the Laplacian
  const bool& add_ECP_dens
)
{
  const int nmo = get_nmo(false);
  vector<double> phi(10 * nmo, 0.0);
  double* phi_temp;
  double chi[10];
  double d[4];
  int iat;
  int l[3];
  double ex;
  double xl[3][3];
  double Grad[3];
  double tau;

  Rho = 0; tau = 0; Elf = 0;
  Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;
  Hess[0] = 0; Hess[1] = 0; Hess[2] = 0;
  Hess[8] = 0; Hess[4] = 0; Hess[5] = 0;

  for (int j = 0; j < nex; j++) {
    iat = get_center(j) - 1;

    type2vector(get_type(j), l);
    d[0] = PosGrid[0] - atoms[iat].x;
    d[1] = PosGrid[1] - atoms[iat].y;
    d[2] = PosGrid[2] - atoms[iat].z;
    d[3] = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    if (add_ECP_dens && has_ECPs && atoms[iat].ECP_electrons != 0) {
      Rho += 8 * atoms[iat].ECP_electrons * exp(-4 * PI * d[3]);
    }
    double temp = -get_exponent(j) * (d[3]);
    if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
      continue;
    ex = exp(temp);
    for (int k = 0; k < 3; k++) {
      if (l[k] == 0) {
        xl[k][0] = 1.0;
        xl[k][1] = 0.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 1) {
        xl[k][0] = d[k];
        xl[k][1] = 1.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 2) {
        xl[k][0] = d[k] * d[k];
        xl[k][1] = 2 * d[k];
        xl[k][2] = 2;
      }
      else if (l[k] == 3) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d[k];
        xl[k][1] = 3 * d2;
        xl[k][2] = 6 * d[k];
      }
      else if (l[k] == 4) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d2;
        xl[k][1] = 4 * d2 * d[k];
        xl[k][2] = 12 * d2;
      }
      else {
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
    for (int mo = 0; mo < nmo; mo++) {
      phi_temp = &phi[mo * 10];
      for (int i = 0; i < 10; i++)
        //				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
        phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i];      //build MO values at this point
    }
  }

  for (int mo = 0; mo < nmo; mo++) {
    const double occ = get_MO_occ(mo);
    const double docc = 2 * occ;
    if (occ != 0) {
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
  if (Rho > 0) {
    normGrad = alpha * sqrt(Grad[0] * Grad[0] + Grad[1] * Grad[1] + Grad[2] * Grad[2]) / pow(Rho, c_43);
    Elf = 1 / (1 + pow(ctelf * pow(Rho, c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
    Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);
  }
  Lap = Hess[0] + Hess[4] + Hess[8];
};

void WFN::computeELIELF(
  const double* PosGrid,		// [3] vector with current position on te grid
  double& Elf,			// Value of the ELF
  double& Eli 			// Value of the ELI
)
{
  const int nmo = get_nmo(false);
  vector<double> phi(4 * nmo, 0.0);
  double* phi_temp;
  double chi[4];
  double d[3];
  int iat;
  int l[3];
  double ex;
  double xl[3][3];

  for (int j = 0; j < nex; j++) {
    iat = get_center(j) - 1;

    type2vector(get_type(j), l);
    d[0] = PosGrid[0] - atoms[iat].x;
    d[1] = PosGrid[1] - atoms[iat].y;
    d[2] = PosGrid[2] - atoms[iat].z;
    double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
      continue;
    ex = exp(temp);
    for (int k = 0; k < 3; k++) {
      if (l[k] == 0) {
        xl[k][0] = 1.0;
        xl[k][1] = 0.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 1) {
        xl[k][0] = d[k];
        xl[k][1] = 1.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 2) {
        xl[k][0] = d[k] * d[k];
        xl[k][1] = 2 * d[k];
        xl[k][2] = 2;
      }
      else if (l[k] == 3) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d[k];
        xl[k][1] = 3 * d2;
        xl[k][2] = 6 * d[k];
      }
      else if (l[k] == 4) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d2;
        xl[k][1] = 4 * d2 * d[k];
        xl[k][2] = 12 * d2;
      }
      else {
        return;
      }
    }
    double ex2 = 2 * get_exponent(j);
    chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
    chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
    chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
    chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
    for (int mo = 0; mo < nmo; mo++) {
      phi_temp = &phi[mo * 4];
      for (int i = 0; i < 4; i++)
        //				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
        phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i];      //build MO values at this point
    }
  }

  double Grad[3];
  double tau;
  double Rho;

  Rho = 0; tau = 0;
  Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;

  for (int mo = 0; mo < nmo; mo++) {
    const double occ = get_MO_occ(mo);
    const double docc = 2 * occ;
    if (occ != 0) {
      phi_temp = &phi[mo * 4];
      Rho += occ * pow(*phi_temp, 2);
      Grad[0] += docc * *phi_temp * phi_temp[1];
      Grad[1] += docc * *phi_temp * phi_temp[2];
      Grad[2] += docc * *phi_temp * phi_temp[3];
      tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
    }
  }
  if (Rho > 0) {
    Elf = 1 / (1 + pow(ctelf * pow(Rho, c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
    Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);
  }
};

void WFN::computeELI(
  const double* PosGrid,		// [3] vector with current position on te grid
  double& Eli 			// Value of the ELI
)
{
  const int nmo = get_nmo(false);
  vector<double> phi(4 * nmo, 0.0);
  double* phi_temp;
  double chi[4];
  double d[3];
  int iat;
  int l[3];
  double ex;
  double xl[3][3];

  for (int j = 0; j < nex; j++) {
    iat = get_center(j) - 1;

    type2vector(get_type(j), l);
    d[0] = PosGrid[0] - atoms[iat].x;
    d[1] = PosGrid[1] - atoms[iat].y;
    d[2] = PosGrid[2] - atoms[iat].z;
    double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
      continue;
    ex = exp(temp);
    for (int k = 0; k < 3; k++) {
      if (l[k] == 0) {
        xl[k][0] = 1.0;
        xl[k][1] = 0.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 1) {
        xl[k][0] = d[k];
        xl[k][1] = 1.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 2) {
        xl[k][0] = d[k] * d[k];
        xl[k][1] = 2 * d[k];
        xl[k][2] = 2;
      }
      else if (l[k] == 3) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d[k];
        xl[k][1] = 3 * d2;
        xl[k][2] = 6 * d[k];
      }
      else if (l[k] == 4) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d2;
        xl[k][1] = 4 * d2 * d[k];
        xl[k][2] = 12 * d2;
      }
      else {
        return;
      }
    }
    double ex2 = 2 * get_exponent(j);
    chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
    chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
    chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
    chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
    for (int mo = 0; mo < nmo; mo++) {
      phi_temp = &phi[mo * 4];
      for (int i = 0; i < 4; i++)
        //				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
        phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i];      //build MO values at this point
    }
  }

  double Grad[3];
  double tau;
  double Rho;

  Rho = 0; tau = 0;
  Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;

  for (int mo = 0; mo < nmo; mo++) {
    const double occ = get_MO_occ(mo);
    const double docc = 2 * occ;
    if (occ != 0) {
      phi_temp = &phi[mo * 4];
      Rho += occ * pow(*phi_temp, 2);
      Grad[0] += docc * *phi_temp * phi_temp[1];
      Grad[1] += docc * *phi_temp * phi_temp[2];
      Grad[2] += docc * *phi_temp * phi_temp[3];
      tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
    }
  }
  Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);
};

void WFN::computeELF(
  const double* PosGrid,		// [3] vector with current position on te grid
  double& Elf			// Value of the ELF
)
{
  const int nmo = get_nmo(false);
  vector<double> phi(4 * nmo, 0.0);
  double* phi_temp;
  double chi[4];
  double d[3];
  int iat;
  int l[3];
  double ex;
  double xl[3][3];

  for (int j = 0; j < nex; j++) {
    iat = get_center(j) - 1;

    type2vector(get_type(j), l);
    d[0] = PosGrid[0] - atoms[iat].x;
    d[1] = PosGrid[1] - atoms[iat].y;
    d[2] = PosGrid[2] - atoms[iat].z;
    double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
      continue;
    ex = exp(temp);
    for (int k = 0; k < 3; k++) {
      if (l[k] == 0) {
        xl[k][0] = 1.0;
        xl[k][1] = 0.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 1) {
        xl[k][0] = d[k];
        xl[k][1] = 1.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 2) {
        xl[k][0] = d[k] * d[k];
        xl[k][1] = 2 * d[k];
        xl[k][2] = 2;
      }
      else if (l[k] == 3) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d[k];
        xl[k][1] = 3 * d2;
        xl[k][2] = 6 * d[k];
      }
      else if (l[k] == 4) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d2;
        xl[k][1] = 4 * d2 * d[k];
        xl[k][2] = 12 * d2;
      }
      else {
        return;
      }
    }
    double ex2 = 2 * get_exponent(j);
    chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
    chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
    chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
    chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
    for (int mo = 0; mo < nmo; mo++) {
      phi_temp = &phi[mo * 4];
      for (int i = 0; i < 4; i++)
        //				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
        phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i];      //build MO values at this point
    }
  }

  double Grad[3];
  double tau;
  double Rho;

  Rho = 0; tau = 0;
  Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;

  for (int mo = 0; mo < nmo; mo++) {
    const double occ = get_MO_occ(mo);
    const double docc = 2 * occ;
    if (occ != 0) {
      phi_temp = &phi[mo * 4];
      Rho += occ * pow(*phi_temp, 2);
      Grad[0] += docc * *phi_temp * phi_temp[1];
      Grad[1] += docc * *phi_temp * phi_temp[2];
      Grad[2] += docc * *phi_temp * phi_temp[3];
      tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
    }
  }
  Elf = 1 / (1 + pow(ctelf * pow(Rho, c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
};

void WFN::computeLapELIELF(
  const double* PosGrid,		// [3] vector with current position on te grid
  double& Elf,			// Value of the ELF
  double& Eli,			// Value of the ELI
  double& Lap				// Value for the Laplacian
)
{
  const int nmo = get_nmo(false);
  vector<double> phi(7 * nmo, 0.0);
  double* phi_temp;
  double chi[7];
  double d[3];
  int iat;
  int l[3];
  double ex;
  double xl[3][3];

  for (int j = 0; j < nex; j++) {
    iat = get_center(j) - 1;

    type2vector(get_type(j), l);
    d[0] = PosGrid[0] - atoms[iat].x;
    d[1] = PosGrid[1] - atoms[iat].y;
    d[2] = PosGrid[2] - atoms[iat].z;
    double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
      continue;
    ex = exp(temp);
    for (int k = 0; k < 3; k++) {
      if (l[k] == 0) {
        xl[k][0] = 1.0;
        xl[k][1] = 0.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 1) {
        xl[k][0] = d[k];
        xl[k][1] = 1.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 2) {
        xl[k][0] = d[k] * d[k];
        xl[k][1] = 2 * d[k];
        xl[k][2] = 2;
      }
      else if (l[k] == 3) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d[k];
        xl[k][1] = 3 * d2;
        xl[k][2] = 6 * d[k];
      }
      else if (l[k] == 4) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d2;
        xl[k][1] = 4 * d2 * d[k];
        xl[k][2] = 12 * d2;
      }
      else {
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
    for (int mo = 0; mo < nmo; mo++) {
      phi_temp = &phi[mo * 7];
      for (int i = 0; i < 7; i++)
        //				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
        phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i];      //build MO values at this point
    }
  }

  double Grad[3];
  double Hess[3];
  double tau;
  double Rho;

  Rho = 0; tau = 0; Elf = 0;
  Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;
  Hess[0] = 0; Hess[1] = 0; Hess[2] = 0;

  for (int mo = 0; mo < nmo; mo++) {
    const double occ = get_MO_occ(mo);
    const double docc = 2 * occ;
    if (occ != 0) {
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
  Elf = 1 / (1 + pow(ctelf * pow(Rho, c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
  Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);
  Lap = Hess[0] + Hess[1] + Hess[2];
};

void WFN::computeLapELI(
  const double* PosGrid,		// [3] vector with current position on te grid
  double& Eli,			// Value of the ELI
  double& Lap				// Value for the Laplacian
)
{
  const int nmo = get_nmo(false);
  vector<double> phi(7 * nmo, 0.0);
  double* phi_temp;
  double chi[7];
  double d[3];
  int iat;
  int l[3];
  double ex;
  double xl[3][3];

  for (int j = 0; j < nex; j++) {
    iat = get_center(j) - 1;

    type2vector(get_type(j), l);
    d[0] = PosGrid[0] - atoms[iat].x;
    d[1] = PosGrid[1] - atoms[iat].y;
    d[2] = PosGrid[2] - atoms[iat].z;
    double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
      continue;
    ex = exp(temp);
    for (int k = 0; k < 3; k++) {
      if (l[k] == 0) {
        xl[k][0] = 1.0;
        xl[k][1] = 0.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 1) {
        xl[k][0] = d[k];
        xl[k][1] = 1.0;
        xl[k][2] = 0.0;
      }
      else if (l[k] == 2) {
        xl[k][0] = d[k] * d[k];
        xl[k][1] = 2 * d[k];
        xl[k][2] = 2;
      }
      else if (l[k] == 3) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d[k];
        xl[k][1] = 3 * d2;
        xl[k][2] = 6 * d[k];
      }
      else if (l[k] == 4) {
        double d2 = d[k] * d[k];
        xl[k][0] = d2 * d2;
        xl[k][1] = 4 * d2 * d[k];
        xl[k][2] = 12 * d2;
      }
      else {
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
    for (int mo = 0; mo < nmo; mo++) {
      phi_temp = &phi[mo * 7];
      for (int i = 0; i < 7; i++)
        //				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
        phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i];      //build MO values at this point
    }
  }

  double Grad[3];
  double Hess[3];
  double tau;
  double Rho;

  Rho = 0; tau = 0;
  Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;
  Hess[0] = 0; Hess[1] = 0; Hess[2] = 0;

  for (int mo = 0; mo < nmo; mo++) {
    const double occ = get_MO_occ(mo);
    const double docc = 2 * occ;
    if (occ != 0) {
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
  Eli = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);
  Lap = Hess[0] + Hess[1] + Hess[2];
};

double WFN::computeMO(
  const double* PosGrid,			// [3] array with current position on the grid
  const int& mo
)
{
  double result = 0.0;
  int iat;
  int l[3];
  double ex;
  double temp;

  // x, y, z and dsqd
  vector<vector<double> > d(4);
  for (int i = 0; i < 4; i++)
    d[i].resize(ncen);

  for (iat = 0; iat < ncen; iat++) {
    d[0][iat] = PosGrid[0] - atoms[iat].x;
    d[1][iat] = PosGrid[1] - atoms[iat].y;
    d[2][iat] = PosGrid[2] - atoms[iat].z;
    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
  }

  for (int j = 0; j < nex; j++) {
    iat = get_center(j) - 1;
    //if (iat != atom) continue;
    type2vector(get_type(j), l);
    temp = -get_exponent(j) * d[3][iat];
    if (temp < -46.0517) //corresponds to cutoff of ex ~< 1E-20
      continue;
    ex = exp(temp);
    for (int k = 0; k < 3; k++) {
      if (l[k] == 0)		continue;
      else if (l[k] == 1)	ex *= d[k][iat];
      else if (l[k] == 2)	ex *= d[k][iat] * d[k][iat];
      else if (l[k] == 3)	ex *= pow(d[k][iat], 3);
      else if (l[k] == 4)	ex *= pow(d[k][iat], 4);
      else if (l[k] == 5)	ex *= pow(d[k][iat], 5);
    }
    result += MOs[mo].get_coefficient_f(j) * ex;      //build MO values at this point
  }
  shrink_vector<vector<double>>(d);
  return result;
}

double Integrate(int& m, double i, double& expn)
{
  int x;
  if (i <= 10) {
    if (expn == 0.0)
      return 0.0;
    double a = m + 0.5;
    double term = 1 / a;
    double partsum = term;
    for (x = 2; x < 50; x++) {
      a++;
      term *= (i / a);
      partsum += term;
      if (term / partsum < 1E-8)
        return 0.5 * partsum * expn;
    }
  }
  else {
    double a = m;
    double b = a + 0.5;
    a -= 0.5;
    double id = 1 / i;
    double approx = 0.88622692 * sqrt(id) * pow(id, m);
    for (x = 0; x < m; x++) {
      b--;
      approx *= b;
    }
    double mult = 0.5 * expn * id;
    if (mult == 0)
      return approx;
    double prop = mult / approx;
    double term = 1;
    double partsum = 1;
    for (x = 1; x < i + m; x++) {
      term *= a * id;
      partsum += term;
      if (abs(term * prop / partsum) < 1E-8)
        return approx - mult * partsum;
      a--;
    }
  }
  return -1;
};
double WFN::fj(int& j, int& l, int& m, double& aa, double& bb)
{
  double temp = 0.0;
  double temp2;
  int a, b;
  for (int i = max(0, j - m); i <= min(j, l); i++) {
    //pre = factorial[l] / factorial[l - i] / factorial[i] * factorial[m] / factorial[m - j + i] / factorial[j - i];
    temp2 = pre[j][l][m][i];
    a = l - i;
    b = m + i - j;
    if (a != 0)	temp2 *= pow(aa, a);
    if (b != 0)	temp2 *= pow(bb, b);
    temp += temp2;
  }
  return temp;
};

double WFN::Afac(int& l, int& r, int& i, double& PC, double& gamma, double& fjtmp)
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

double WFN::computeESP(double* PosGrid, vector<vector<double> >& d2)
{
  double ESP = 0;
  double P[3];
  double Pi[3];
  double Pj[3];
  double PC[3];
  double Fn[11];
  double Al[54];
  double Am[54];
  double An[54];
  int maplrsl[54];
  int maplrsm[54];
  int maplrsn[54];
  int l_i[3];
  int l_j[3];
  int iat, jat, MaxFn;
  double ex_sum, sqd, sqpc, prefac, expc, term, addesp, fjtmp, c, twoexpc, iex, jex;

  const int MO = get_nmo();
  const int nprim = get_nex();

  double temp;
  int maxl, maxm, maxn;
  vector<vector<double> > pos(3);
  for (int i = 0; i < 3; i++) pos[i].resize(get_ncen());

  for (int iat = 0; iat < get_ncen(); iat++) {
    pos[0][iat] = atoms[iat].x;
    pos[1][iat] = atoms[iat].y;
    pos[2][iat] = atoms[iat].z;
    ESP += get_atom_charge(iat) * pow(sqrt(pow(PosGrid[0] - pos[0][iat], 2) + pow(PosGrid[1] - pos[1][iat], 2) + pow(PosGrid[2] - pos[2][iat], 2)), -1);
  }

  for (int iprim = 0; iprim < nprim; iprim++) {
    iat = get_center(iprim) - 1;
    type2vector(get_type(iprim), l_i);
    iex = get_exponent(iprim);
    for (int jprim = iprim; jprim < nprim; jprim++) {
      jat = get_center(jprim) - 1;
      type2vector(get_type(jprim), l_j);
      ex_sum = get_exponent(iprim) + get_exponent(jprim);
      jex = get_exponent(jprim);

      sqd = d2[iat][jat];

      prefac = TWO_PI / ex_sum * exp(-iex * jex * sqd / ex_sum);
      if (prefac < 1E-10) continue;

      for (int i = 0; i < 3; i++) {
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
      for (int l = 0; l <= l_i[0] + l_j[0]; l++) {
        if (l % 2 != 1)
          fjtmp = fj(l, l_i[0], l_j[0], Pi[0], Pj[0]);// *factorial[l];
        else
          fjtmp = -fj(l, l_i[0], l_j[0], Pi[0], Pj[0]);// * factorial[l];
        for (int r = 0; r <= l / 2; r++)
          for (int s = 0; s <= (l - 2 * r) / 2; s++) {
            maxl++;
            Al[maxl] = Afac(l, r, s, PC[0], ex_sum, fjtmp);
            maplrsl[maxl] = l - 2 * r - s;
          }
      }
      maxm = -1;
      for (int l = 0; l <= l_i[1] + l_j[1]; l++) {
        if (l % 2 != 1)
          fjtmp = fj(l, l_i[1], l_j[1], Pi[1], Pj[1]);// *factorial[l];
        else
          fjtmp = -fj(l, l_i[1], l_j[1], Pi[1], Pj[1]);// * factorial[l];
        for (int r = 0; r <= l / 2; r++)
          for (int s = 0; s <= (l - 2 * r) / 2; s++) {
            maxm++;
            Am[maxm] = Afac(l, r, s, PC[1], ex_sum, fjtmp);
            maplrsm[maxm] = l - 2 * r - s;
          }
      }
      maxn = -1;
      for (int l = 0; l <= l_i[2] + l_j[2]; l++) {
        if (l % 2 != 1)
          fjtmp = fj(l, l_i[2], l_j[2], Pi[2], Pj[2]);// *factorial[l];
        else
          fjtmp = -fj(l, l_i[2], l_j[2], Pi[2], Pj[2]);// * factorial[l];
        for (int r = 0; r <= l / 2; r++)
          for (int s = 0; s <= (l - 2 * r) / 2; s++) {
            maxn++;
            An[maxn] = Afac(l, r, s, PC[2], ex_sum, fjtmp);
            maplrsn[maxn] = l - 2 * r - s;
          }
      }

      term = 0.0;
      for (int l = 0; l < maxl; l++) {
        if (Al[l] == 0)
          continue;
        for (int m = 0; m < maxm; m++) {
          if (Am[m] == 0)
            continue;
          for (int n = 0; n < maxn; n++) {
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