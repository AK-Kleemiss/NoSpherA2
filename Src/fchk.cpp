#include "pch.h"
#include "fchk.h"
#include "convenience.h"
#include "basis_set.h"
#include "constants.h"
//----------------------------FCHK Preparation and Gaussian--------------------------------------
/*
KEPT AS AN EXAMPLE HOW TO CALL G09 FROM WITHIN C++
bool gaussian(const string &programPath, const bool &debug){
  if(debug){
    debug_fchk=true;
  }
  string workpath;
  if(programPath.size()<3) workpath="/usr/local/g09/g09";
  else workpath=programPath;
  bool success=false;
  string envbsd=workpath;
  envbsd.erase(envbsd.find("/g09/g09")+4,4);
  envbsd.append("/bsd");
  string envlocal=workpath;
  envlocal.erase(envlocal.find("/g09/g09")+4,4);
  envlocal.append("/local");
  string envextras=workpath;
  envextras.erase(envextras.find("/g09/g09")+4,4);
  envextras.append("/extras");
  string envshort=workpath;
  envshort.erase(envshort.find("/g09/g09")+4,4);
  if(debug_fchk){
   std::cout << envbsd << endl;
   std::cout << envlocal << endl;
   std::cout << envextras << endl;
   std::cout << envshort << endl;
  }
  string EXEDIR=("GAUSS_EXEDIR=");
  EXEDIR.append(envbsd);
  EXEDIR.append(":");
  EXEDIR.append(envlocal);
  EXEDIR.append(":");
  EXEDIR.append(envextras);
  EXEDIR.append(":");
  EXEDIR.append(envshort);
  if(debug_fchk)std::cout << EXEDIR << endl;
  string BSDDIR=("GAUSS_BSDDIR=");
  BSDDIR.append(envbsd);
  if(debug_fchk)std::cout << BSDDIR << endl;
  char *exedir = &EXEDIR[0u];
  char *bsddir = &BSDDIR[0u];
  cls();
 std::cout << "Running gaussian... Please wait... " << endl;
  string execute_command;
  execute_command=workpath;
  execute_command.append(" gaussian.com");
  if(system(execute_command.c_str()))std::cout << "Well system returned non 0... No idea what this means now..." << endl;
  if(!debug)std::cout << "Finished gaussian!!" << endl;
  ifstream ifile("gaussian.log",ios::in);
  string line;
  string preline;
  while (!ifile.eof()){
    preline=line;
    getline(ifile,line);
  }
  if(line.find("Normal termination")==-1&&preline.find("Normal termination")==-1) success=false;
  else success=true;
  return success;
};

string prepare_gaussian(const string& basis_set_path, const string& fchkname, WFN& wave, const int& ncpus, const float& mem, bool debug) {
  if (exists("gaussian.com")) {
   std::cout << "gaussian.com already exists, do you want me to overwrite it?";
    if (!yesno()) return "WRONG";
  }
  if (wave.get_modified()) {
   std::cout << "The wavefunction has been modified after reading. Please make sure that you know what you are doing!" << endl;
   std::cout << "Do you want to continue?" << flush;
    if (!yesno()) return "WRONG";
  }
  //-----------------READING BASIS SET--------------------------------
  if (wave.get_nr_basis_set_loaded() == 0) {
   std::cout << "No basis set loaded, will load a complete basis set now!" << endl;
    if (!read_basis_set_vanilla(basis_set_path, wave, debug, false)) {
     std::cout << "Problem during reading of the basis set!" << endl;
      return "WRONG";
    }
  }
  else if (wave.get_nr_basis_set_loaded() < wave.get_ncen()) {
   std::cout << "Not all atoms have a basis set loaded!" << endl
      << "Do you want to laod the missing atoms?" << flush;
    if (!yesno()) {
     std::cout << "Do you want to load a new basis set?" << flush;
      if (!yesno()) {
       std::cout << "Okay, aborting then!!!" << endl;
        return "WRONG";
      }
     std::cout << "deleting old one..." << endl;
      if (!wave.delete_basis_set()) {
       std::cout << "ERROR while deleting a basis set!" << endl;
        return "WRONG";
      }
      else if (!read_basis_set_vanilla(basis_set_path, wave, debug, true)) {
       std::cout << "ERROR during reading of the new basis set!" << endl;
        return "WRONG";
      }
    }
    else {
     std::cout << "okay, loading the missing atoms.." << endl;
      if (!read_basis_set_missing(basis_set_path, wave, debug)) {
       std::cout << "ERROR during reading of missing basis set!" << endl;
        return "WRONG";
      }
    }
  }
  else if (wave.get_nr_basis_set_loaded() == wave.get_ncen()) {
   std::cout << "There already is a basis set loaded!" << endl
      << "Do you want me to delete the old one and load a new one?" << flush;
    if (!yesno()) {
     std::cout << "Okay, continuing with the old one..." << endl;
    }
    else {
     std::cout << "Deleting the old basis set!" << endl;
      if (!wave.delete_basis_set()) {
       std::cout << "ERROR during deleting of the basis set!";
        return "WRONG";
      }
     std::cout << "Going to load a new one now!" << endl;
      if (!read_basis_set_vanilla(basis_set_path, wave, debug, true)) {
       std::cout << "Problem during reading of the basis set!" << endl;
        return "WRONG";
      }
    }
  }
  else {
   std::cout << "# of loaded > # atoms" << endl
      << "Sorry, this should not happen... aborting!!!" << endl;
    return "WRONG";
  }
  //-----------------------check ordering and order accordingly----------------------
  if (!wave.sort_wfn(wave.check_order(debug), debug)) {
   std::cout << "Could not order the wavefunction, aborting!" << endl;
    return "WRONG";
  }
  //------------setting up the .com file with basis set---------------
  ofstream com;
  com.open("gaussian.com", ios::out);
  string temp;
  temp = "gaussian.chk";
  int run = 0;
  if (exists(temp)) {
   std::cout << "Do you want me to overwrite gaussian.chk?";
    if (!yesno()) {
     std::cout << "Then pelase give a new .chk filename: ";
      string filename;
      cin >> filename;
      if (filename.find(".chk") == -1) {
       std::cout << "this does not look like a .chk file! Please try again!";
        return "WRONG";
      }
    }
  }
  com << "%chk=" << temp << endl;
  com << "%mem=" << mem << "GB" << endl;
  com << "%nproc=" << ncpus << endl;
  if(wave.get_d_f_switch()) com << "# SCF=(MaxCycle=3000,Conver=-1) b3lyp/gen 5D 7F nosymm IOp(3/32=2)" << endl << endl;
    else
  com << "# SCF=(MaxCycle=3000,Conver=-1) rhf/gen 6D 10F nosymm IOp(3/32=2)" << endl << endl;
  com << "TITLE" << endl << endl;
  com << wave.get_charge() << " " << wave.get_multi() << endl;
  com << wave.get_centers(check_bohr(wave, true, debug));
  //	com << wave.get_centers (false);
  com << endl;
  svec elements_list;
  if (debug)std::cout << "elements_list.size()= " << elements_list.size() << endl;
  for (int a = 0; a < wave.get_ncen(); a++) {
    string label_temp;
    label_temp = wave.get_atom_label(a);
    bool found = false;
    if (elements_list.size() != 0) {
      for (int i = 0; i < elements_list.size(); i++) {
        if (elements_list[i].compare(label_temp) == 0) {
          found = true;
          if (debug)std::cout << "Found " << label_temp << " in the elements_list!" << endl;
        }
      }
    }
    if (!found) {
      elements_list.push_back(label_temp);
      if (debug)std::cout << "elements_list.size()= " << elements_list.size() << endl;
      com << label_temp << " 0" << endl;
      int working_shell = -1;
      for (int p = 0; p < wave.get_atom_primitive_count(a); p++) {
        int temp_shell = wave.get_basis_set_shell(a, p);
        if (temp_shell != working_shell) {
          working_shell = temp_shell;
          switch (wave.get_shell_type(a, temp_shell)) {
          case 1:
            com << "S " << wave.get_atom_shell_primitives(a, temp_shell) << " 1.00" << endl;
            break;
          case 2:
            com << "P " << wave.get_atom_shell_primitives(a, temp_shell) << " 1.00" << endl;
            break;
          case 3:
            com << "D " << wave.get_atom_shell_primitives(a, temp_shell) << " 1.00" << endl;
            break;
          case 4:
            com << "F " << wave.get_atom_shell_primitives(a, temp_shell) << " 1.00" << endl;
            break;
          }
        }
        com << scientific << setw(17) << setprecision(10) << wave.get_atom_basis_set_exponent(a, p);
        com << " ";
        com << scientific << setw(17) << setprecision(10) << wave.get_atom_basis_set_coefficient(a, p);
        com << endl;
      }
      com << "****" << endl;
    }
  }
  com << endl;
  com.flush();
  com.close();
  if (debug) {
   std::cout << "Wrote the gaussian Input! Please check it before i continue..." << endl;
  }
  return temp;
};

bool modify_fchk(const string& fchk_name, const string& basis_set_path, WFN& wave, bool& debug, const bool& read) {
  wave.set_modified();
  if (debug) debug_fchk = true;
  vec CMO;
  int nao = 0;
  int nshell = 0;
  int naotr = 0;
  if (debug) {
   std::cout << "Origin: " << wave.get_origin() << endl;
  }
  if (wave.get_origin() == 2 || wave.get_origin() == 4) {
    //-----------------READING BASIS SET--------------------------------
    if (read) {
      if (!read_basis_set_vanilla(basis_set_path, wave, debug, false)) {
       std::cout << "Problem during reading of the basis set!" << endl;
        return false;
      }
    }
    double pi = 3.14159265358979;
    //---------------normalize basis set---------------------------------
    if (debug)std::cout << "starting to normalize the basis set" << endl;
    vec norm_const;
    //-----------debug output---------------------------------------------------------
    if (debug) {
     std::cout << "exemplary output before norm_const of the first atom with all it's properties: " << endl;
      wave.print_atom_long(0);
     std::cout << "ended normalizing the basis set, now for the MO_coeffs" << endl;
     std::cout << "Status report:" << endl;
     std::cout << "size of norm_const: " << norm_const.size() << endl;
     std::cout << "WFN MO counter: " << wave.get_nmo() << endl;
     std::cout << "Number of atoms: " << wave.get_ncen() << endl;
     std::cout << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
     std::cout << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
    }
    //-----------------------check ordering and order accordingly----------------------
    wave.sort_wfn(wave.check_order(debug), debug);
    //-------------------normalize the basis set shell wise---------------------
    for (int a = 0; a < wave.get_ncen(); a++) {
      for (int p = 0; p < wave.get_atom_primitive_count(a); p++) {
        double temp = wave.get_atom_basis_set_exponent(a, p);
        switch (wave.get_atom_primitive_type(a, p)) {
        case 1:
          temp = 2 * temp / pi;
          temp = pow(temp, 0.75);
          temp = temp * wave.get_atom_basis_set_coefficient(a, p);
          wave.change_atom_basis_set_coefficient(a, p, temp);
          break;
        case 2:
          temp = 128 * pow(temp, 5);
          temp = temp / pow(pi, 3);
          temp = pow(temp, 0.25);
          temp = wave.get_atom_basis_set_coefficient(a, p) * temp;
          wave.change_atom_basis_set_coefficient(a, p, temp);
          break;
        case 3:
          temp = 2048 * pow(temp, 7);
          temp = temp / (9 * pow(pi, 3));
          temp = pow(temp, 0.25);
          temp = wave.get_atom_basis_set_coefficient(a, p) * temp;
          wave.change_atom_basis_set_coefficient(a, p, temp);
          break;
        case 4:
          temp = 32768 * pow(temp, 9);
          temp = temp / (225 * pow(pi, 3));
          temp = pow(temp, 0.25);
          temp = wave.get_atom_basis_set_coefficient(a, p) * temp;
          wave.change_atom_basis_set_coefficient(a, p, temp);
          break;
        case -1:
         std::cout << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy..." << endl;
          break;
        }
      }
    }
    for (int a = 0; a < wave.get_ncen(); a++) {
      double aiaj = 0.0;
      double factor = 0.0;
      for (int s = 0; s < wave.get_atom_shell_count(a); s++) {
        int type_temp = wave.get_shell_type(a, s);
        if (type_temp == -1) {
         std::cout << "ERROR in type assignement!!" << endl;
        }
        if (debug) {
         std::cout << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << endl
          << "start: " << wave.get_shell_start(a, s, false) << flush
          << " stop: " << wave.get_shell_end(a, s, false) << flush << endl
          << "factor: ";
        }
        switch (type_temp) {
        case 1:
          factor = 0;
          for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
            for (int j = wave.get_shell_start(a, s, false); j <= wave.get_shell_end(a, s, false); j++) {
              aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
              double term = (pi / aiaj);
              term = pow(term, 1.5);
              factor += wave.get_atom_basis_set_coefficient(a, i) * wave.get_atom_basis_set_coefficient(a, j) * term;
            }
          }
          if (factor == 0) return false;
          factor = pow(factor, -0.5);
          if (debug)std::cout << factor << endl;
          for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
            if (debug) {
             std::cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i)
                << " after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
            }
            wave.change_atom_basis_set_coefficient(a, i, factor * wave.get_atom_basis_set_coefficient(a, i));
            norm_const.push_back(wave.get_atom_basis_set_coefficient(a, i));
          }
          break;
        case 2:
          factor = 0;
          for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
            for (int j = wave.get_shell_start(a, s, false); j <= wave.get_shell_end(a, s, false); j++) {
              aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
              double term = 4 * pow(aiaj, 5);
              term = pow(pi, 3) / term;
              term = pow(term, 0.5);
              factor += wave.get_atom_basis_set_coefficient(a, i) * wave.get_atom_basis_set_coefficient(a, j) * term;
            }
          }
          if (factor == 0) return false;
          factor = pow(factor, -0.5);
          if (debug)std::cout << factor << endl;
          for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
            if (debug) {
             std::cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i)
                << " after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
            }
            wave.change_atom_basis_set_coefficient(a, i, factor * wave.get_atom_basis_set_coefficient(a, i));
            for (int k = 0; k < 3; k++) norm_const.push_back(wave.get_atom_basis_set_coefficient(a, i));
          }
          break;
        case 3:
          factor = 0;
          for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
            for (int j = wave.get_shell_start(a, s, false); j <= wave.get_shell_end(a, s, false); j++) {
              aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
              double term = 16 * pow(aiaj, 7);
              term = pow(pi, 3) / term;
              term = pow(term, 0.5);
              factor += wave.get_atom_basis_set_coefficient(a, i) * wave.get_atom_basis_set_coefficient(a, j) * term;
            }
          }
          if (factor == 0) return false;
          factor = (pow(factor, -0.5)) / sqrt(3);
          if (debug)std::cout << factor << endl;
          for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
            if (debug) {
             std::cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i)
                << " after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
            }
            wave.change_atom_basis_set_coefficient(a, i, factor * wave.get_atom_basis_set_coefficient(a, i));
            for (int k = 0; k < 3; k++) norm_const.push_back(wave.get_atom_basis_set_coefficient(a, i));
            for (int k = 0; k < 3; k++) norm_const.push_back(sqrt(3) * wave.get_atom_basis_set_coefficient(a, i));
          }
          break;
        case 4:
          factor = 0;
          for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
            for (int j = wave.get_shell_start(a, s, false); j <= wave.get_shell_end(a, s, false); j++) {
              aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
              double term = 64 * pow((aiaj), 9);
              term = pow(pi, 3) / term;
              term = pow(term, 0.5);
              factor += wave.get_atom_basis_set_coefficient(a, i) * wave.get_atom_basis_set_coefficient(a, j) * term;
            }
          }
          if (factor == 0) return false;
          factor = pow(factor, -0.5) / sqrt(15);
          if (debug)std::cout << factor << endl;
          for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
            if (debug) {
             std::cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i)
                << " after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
            }
            wave.change_atom_basis_set_coefficient(a, i, factor * wave.get_atom_basis_set_coefficient(a, i));
            for (int l = 0; l < 3; l++) norm_const.push_back(wave.get_atom_basis_set_coefficient(a, i));
            for (int l = 0; l < 6; l++) norm_const.push_back(sqrt(5) * wave.get_atom_basis_set_coefficient(a, i));
            norm_const.push_back(sqrt(15) * wave.get_atom_basis_set_coefficient(a, i));
          }
          break;
        }
        if (debug)	std::cout << "This shell has: " << wave.get_shell_end(a, s, false) - wave.get_shell_start(a, s, false) + 1 << " primitives" << endl;
      }
    }
    //-----------debug output---------------------------------------------------------
    if (debug) {
     std::cout << "exemplary output of the first atom with all it's properties: " << endl;
      wave.print_atom_long(0);
     std::cout << "ended normalizing the basis setnow for norm_cprim" << endl;
     std::cout << "Status report:" << endl;
     std::cout << "size of norm_const: " << norm_const.size() << endl;
     std::cout << "WFN MO counter: " << wave.get_nmo() << endl;
     std::cout << "Number of atoms: " << wave.get_ncen() << endl;
     std::cout << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
     std::cout << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
    }
    //---------------------To not mix up anything start normalizing WFN_matrix now--------------------------
    int run = 0;
    ofstream norm_cprim;
    if (debug) norm_cprim.open("norm_prim.debug", ofstream::out);
    for (int m = 0; m < wave.get_nmo(); m++) {
      if (debug) norm_cprim << m << ". MO:" << endl;
      for (int p = 0; p < wave.get_MO_primitive_count(m); p++) {
        if (debug)std::cout << p << ". primitive; " << m << ". MO " << "norm nonst: " << norm_const[p] << endl;
        double temp = wave.get_MO_coef(m, p) / norm_const[p];
        if (debug) {
         std::cout << " temp after normalization: " << temp << endl;
          norm_cprim << " " << temp << endl;
        }
        run++;
        if (!wave.change_MO_coef(m, p, temp)) {
         std::cout << "ERROR in changing the coefficients after normalising!";
          if (debug)std::cout << "m:" << m << " p: " << p << " temp:" << temp;
         std::cout << endl;
          return false;
        }
      }
    }
    if (debug) {
      norm_cprim.flush();
      norm_cprim.close();
     std::cout << "See norm_cprim.debug for the CPRIM vectors" << endl;
     std::cout << "Total count in CPRIM: " << run << endl;
    }
    //--------------Build CMO of alessandro from the first elements of each shell-------------
    for (int m = 0; m < wave.get_nmo(); m++) {
      int run_2 = 0;
      for (int a = 0; a < wave.get_ncen(); a++) {
        for (int s = 0; s < wave.get_atom_shell_count(a); s++) {
          if (debug)std::cout << "Going to load the " << wave.get_shell_start_in_primitives(a, s) << ". value" << endl;
          switch (wave.get_shell_type(a, s)) {
          case 1:
            CMO.push_back(wave.get_MO_coef(m, wave.get_shell_start_in_primitives(a, s)));
            if (m == 0) nao++;
            if (debug && wave.get_atom_shell_primitives(a, s) != 1)
             std::cout << "Pushing back 1 coefficient for S shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives! Shell start is: " << wave.get_shell_start(a, s, false) << endl;
            break;
          case 2:
            for (int i = 0; i < 3; i++) CMO.push_back(wave.get_MO_coef(m, wave.get_shell_start_in_primitives(a, s) + i));
            if (debug && wave.get_atom_shell_primitives(a, s) != 1)
             std::cout << "Pushing back 3 coefficients for P shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
            if (m == 0) nao += 3;
            break;
          case 3:
            for (int i = 0; i < 6; i++) CMO.push_back(wave.get_MO_coef(m, wave.get_shell_start_in_primitives(a, s) + i));
            if (debug && wave.get_atom_shell_primitives(a, s) != 1)
             std::cout << "Pushing back 6 coefficient for D shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
            if (m == 0) nao += 6;
            break;
          case 4:
            //this hardcoded piece is due to the order of f-type functions in the fchk
            for (int i = 0; i < 3; i++) CMO.push_back(wave.get_MO_coef(m, wave.get_shell_start_in_primitives(a, s) + i));
            CMO.push_back(wave.get_MO_coef(m, wave.get_shell_start_in_primitives(a, s) + 6));
            for (int i = 0; i < 2; i++) CMO.push_back(wave.get_MO_coef(m, wave.get_shell_start_in_primitives(a, s) + i + 3));
            for (int i = 0; i < 2; i++) CMO.push_back(wave.get_MO_coef(m, wave.get_shell_start_in_primitives(a, s) + i + 7));
            CMO.push_back(wave.get_MO_coef(m, wave.get_shell_start_in_primitives(a, s) + 5));
            CMO.push_back(wave.get_MO_coef(m, wave.get_shell_start_in_primitives(a, s) + 9));
            if (debug && wave.get_atom_shell_primitives(a, s) != 1)
             std::cout << "Pushing back 10 coefficient for F shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
            if (m == 0) nao += 10;
            break;
          }
          run_2++;
        }
        if (debug)std::cout << "finished with atom!" << endl;
      }
      if (debug)std::cout << "finished with MO!" << endl;
      if (nshell != run_2) nshell = run_2;
    }
    if (debug) {
      ofstream cmo("cmo.debug", ofstream::out);
      for (int p = 0; p < CMO.size(); p++) {
        string temp;
        for (int i = 0; i < 5; i++) {
          stringstream stream;
          stream << scientific << setw(14) << setprecision(7) << CMO[p + i] << " ";
          temp += stream.str();
        }
        p += 4;
        cmo << temp << endl;
      }
      cmo.flush();
      cmo.close();
     std::cout << CMO.size() << " Elements in CMO" << endl;
     std::cout << norm_const.size() << " = nprim" << endl;
     std::cout << nao << " = nao" << endl;
     std::cout << nshell << " = nshell" << endl;
    }
    //------------------ make the DM -----------------------------
    naotr = nao * (nao + 1) / 2;
    vec kp;
    for (int i = 0; i < naotr; i++) wave.push_back_DM(0.0);
    if (debug)std::cout << "I made kp!" << endl << nao << " is the maximum for iu" << endl;
    for (int iu = 0; iu < nao; iu++) {
      for (int iv = 0; iv <= iu; iv++) {
        int iuv = (iu * (iu + 1) / 2) + iv;
        if (debug)std::cout << "iu: " << iu << " iv: " << iv << " iuv: " << iuv << " kp(iu): " << iu * (iu + 1) / 2 << endl;
        for (int m = 0; m < wave.get_nmo(); m++) {
          //if(debug)std::cout << "DM before: " << wave.get_DM(iuv) << endl;
          const int temp1 = iu + (m * nao);
          const int temp2 = iv + (m * nao);
          err_checkf(
            wave.set_DM(iuv, wave.get_DM(iuv) + 2 * CMO[temp1] * CMO[temp2]),
            "Something went wrong while writing the DM! iuv=" +
            to_string(iuv) + " while CMO has size: " +
            to_string(CMO.size()) + "\nCMO[" + to_string(temp1) + "] = " +
            to_string(CMO[temp1]) + "; CMO[" + to_string(temp2) + "] = " +
            to_string(CMO[temp2]) + " DM: " + to_string(wave.get_DM(iuv)) +
            "\nDM_size: " + to_string(wave.get_DM_size()), std::cout);
          //else if (debug)std::cout << "DM after: " << wave.get_DM(iuv) << endl;
        }
      }
    }
    if (debug) {
     std::cout << "DM is:" << endl;
      for (int p = 0; p < wave.get_DM_size(); p++) {
        string temp;
        for (int i = 0; i < 5; i++) {
          stringstream stream;
          stream << scientific << setw(14) << setprecision(7) << wave.get_DM(i + p) << " ";
          temp += stream.str();
        }
        p += 4;
       std::cout << temp << endl;
      }
     std::cout << wave.get_DM_size() << " Elements in DM" << endl;
    }
  }
  // open fchk for copying
  string temp_fchk = fchk_name;
  temp_fchk.append(".fchk");
  if (!exists(temp_fchk)) {
   std::cout << "This is worrysome... The fchk should be there.." << endl;
   std::cout << "fchk_name: " << temp_fchk << endl;
    return false;
  }
  ifstream ifchk(temp_fchk.c_str());
  if (!ifchk.is_open()) {
   std::cout << "ERROR while opening .fchk ifile!" << endl;
    return false;
  }
  string tfn("temp.fchk");
  ofstream ofchk;
  ofchk.open(tfn, ofstream::out);
  string line;
  int dum_nao = 0;
  if (wave.get_origin() == 2 || wave.get_origin() == 4) {
    while (line.find("Alpha Orbital Energies") == -1 && !ifchk.eof()) {
      getline(ifchk, line);
      if (debug)std::cout << "line: " << line << endl;
      if (line.find("Alpha Orbital Energies") == -1) ofchk << line << endl;
      else {
        char tempchar[100];
        size_t length;
        length = line.copy(tempchar, 11, 50);
        tempchar[length] = '\0';
        dum_nao = stoi(tempchar);
        if (debug) {
         std::cout << "nao read from fchk: " << dum_nao << " and from basis set: " << nao << endl;
        }
        ofchk << line << endl;
      }
    }
    int counter = 0;
    while (line.find("Alpha MO") == -1 && !ifchk.eof()) {
      getline(ifchk, line);
      string temp = " ";
      for (int j = 0; j < 5; j++) {
        if (counter + j < wave.get_nmo()) {
          stringstream stream;
          stream << scientific << setw(15) << setprecision(8) << wave.get_MO_energy(counter + j);
          if (j < 4) stream << " ";
          temp += stream.str();
          temp.replace(12 + j * 16, 1, "E");
        }
        else if (counter + j < dum_nao) {
          char tempchar[100];
          size_t length;
          length = line.copy(tempchar, 15, 1 + j * 16);
          tempchar[length] = '\0';
          stringstream stream;
          stream << scientific << setw(15) << setprecision(8) << stod(tempchar);
          if (j < 4) stream << " ";
          temp += stream.str();
          temp.replace(12 + j * 16, 1, "E");
        }
      }
      counter += 5;
      temp += '\n';
      if (temp.size() > 3) ofchk << temp;
    }
    ofchk.flush();
    ofchk << "Alpha MO coefficients                      R   N=" << setw(12) << nao * dum_nao << endl;
    //now write the CMO and skip lines in IFCHK
    for (int i = 0; i < nao * dum_nao; i++) {
      string temp = " ";
      for (int j = 0; j < 5; j++) {
        if (i + j < CMO.size()) {
          stringstream stream;
          stream << scientific << setw(15) << setprecision(8) << CMO[i + j] << " ";
          temp += stream.str();
          temp.replace(12 + j * 16, 1, "E");
        }
        else if (i + j < nao * nao) {
          stringstream stream;
          stream << scientific << setw(15) << setprecision(8) << 0.0 << " ";
          temp += stream.str();
          temp.replace(12 + j * 16, 1, "E");
        }
      }
      i += 4;
      temp += '\n';
      ofchk << temp;
      getline(ifchk, line);
    }
  }
  if (wave.get_origin() == 1) {
    while (line.find("Total SCF Density") == -1 && !ifchk.eof()) {
      getline(ifchk, line);
      if (debug)std::cout << "line: " << line << endl;
      if (line.find("Total SCF Density") == -1) ofchk << line << endl;
    }
    ofchk.flush();
  }
  ofchk << "Total SCF Density                          R   N=" << setw(12) << wave.get_DM_size() << endl;
  getline(ifchk, line);
  //now write the DM and skip lines in IFCHK
  for (int i = 0; i < wave.get_DM_size(); i++) {
    string temp = " ";
    if (debug)std::cout << "i: " << i << " DM_size= " << wave.get_DM_size() << " Element ";
    for (int j = 0; j < 5; j++) {
      if (i + j < wave.get_DM_size()) {
        stringstream stream;
        if (debug)std::cout << i + j << " ";
        stream << scientific << setw(15) << setprecision(8) << wave.get_DM(i + j) << " ";
        if (i + j < wave.get_DM_size()) {
          temp += stream.str();
          temp.replace(12 + j * 16, 1, "E");
        }
      }
    }
    i += 4;
    if (debug)std::cout << endl;
    temp += '\n';
    ofchk << temp;
    getline(ifchk, line);
  }
  while (!ifchk.eof()) {
    getline(ifchk, line);
    ofchk << line;
    if (!ifchk.eof()) ofchk << endl;
  }
  ofchk.flush();
  ofchk.close();
  ifchk.close();
  if (debug) {
   std::cout << "Keeping the old fchk file and the new temp.fchk" << endl;
    copy_file(tfn, temp_fchk);
  }
  else copy_file(tfn, temp_fchk);
  if (remove("temp.fchk") != 0)std::cout << "error deleting temp.fchk!" << endl;
  cls();
  return true;
};
*/
//------------------Functions to read from .fchk files----------------------------------
int read_fchk_integer(const std::string& in)
{
    return std::stoi(in.substr(49, in.length() - 49));
};
double read_fchk_double(const std::string& in)
{
    return std::stod(in.substr(49, in.length() - 49));
};
bool read_fchk_integer_block(std::ifstream& in, const char* heading, ivec& result, bool rewind)
{
    if (result.size() != 0)
        result.clear();
    std::string line = go_get_string(in, heading, rewind);
    int limit = read_fchk_integer(line);
    int run = 0;
    int temp;
    getline(in, line);
    while (run < limit)
    {
        if (in.eof())
            return false;
        temp = stoi(line.substr(12 * (run % 6), 12 * (run % 6 + 1)));
        result.push_back(temp);
        run++;
        if (run % 6 == 0)
            getline(in, line);
    }
    return true;
};
bool read_fchk_double_block(std::ifstream& in, const char* heading, vec& result, bool rewind)
{
    if (result.size() != 0)
        result.clear();
    std::string line = go_get_string(in, heading, rewind);
    int limit = read_fchk_integer(line);
    int run = 0;
    double temp;
    getline(in, line);
    while (run < limit)
    {
        if (in.eof())
            return false;
        temp = stod(line.substr(16 * (run % 5), 16 * (run % 5 + 1)));
        result.push_back(temp);
        run++;
        if (run % 5 == 0)
            getline(in, line);
    }
    return true;
};
int read_fchk_integer(std::ifstream& in, const char* search, bool rewind)
{
    std::string temp = go_get_string(in, search, rewind);
    return std::stoi(temp.substr(49, temp.length() - 49));
};
double read_fchk_double(std::ifstream& in, const char* search, bool rewind)
{
    std::string temp = go_get_string(in, search, rewind);
    return std::stod(temp.substr(49, temp.length() - 49));
};

bool free_fchk(std::ostream &file, const std::filesystem::path &fchk_name, const std::filesystem::path &basis_set_path, WFN &wave, bool &debug, bool force_overwrite)
{
    using namespace std;
    int elcount = 0;
    elcount -= wave.get_charge();
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        elcount += wave.get_atom_charge(i);
    }
    int alpha_els = 0, beta_els = 0, temp_els = elcount;
    while (temp_els > 1)
    {
        alpha_els++;
        beta_els++;
        temp_els -= 2;
    }
    alpha_els += temp_els;
    int diff = wave.get_multi() - 1;
    while (alpha_els - beta_els != diff)
    {
        alpha_els++;
        beta_els--;
    }
    if (debug)
    {
        file << "alpha, beta, elcount: " << setw(5) << alpha_els << setw(5) << beta_els << setw(5) << elcount << endl;
    }
    if (wave.get_nr_basis_set_loaded() == 0)
    {
        if (debug)
            file << "No basis set loaded, will load a complete basis set now!" << endl;
        err_checkf(read_basis_set_vanilla(basis_set_path, wave, debug), "ERROR during reading of missing basis set!", file);
    }
    else if (wave.get_nr_basis_set_loaded() < wave.get_ncen())
    {
        file << "Not all atoms have a basis set loaded!\nLaoding the missing atoms..." << flush;
        err_checkf(read_basis_set_missing(basis_set_path, wave, debug), "ERROR during reading of missing basis set!", file);
    }
    else if (wave.get_nr_basis_set_loaded() > wave.get_ncen())
    {
        err_checkf(false, "# of loaded > # atoms\nSorry, this should not happen... aborting!!!", file);
    }
    // wave.set_modified();
    vec CMO;
    vec CMO_beta;
    if (debug)
    {
        file << "Origin: " << wave.get_origin() << endl;
    }
    if (wave.get_origin() == 2 || wave.get_origin() == 4 || wave.get_origin() == 9 || wave.get_origin() == 8)
    {
        //-----------------------check ordering and order accordingly----------------------
        wave.sort_wfn(wave.check_order(debug), debug);
        //---------------normalize basis set---------------------------------
        if (debug)
            file << "starting to normalize the basis set" << endl;
        vec norm_const;
        //-----------debug output---------------------------------------------------------
        if (debug)
        {
            file << "exemplary output before norm_const of the first atom with all it's properties: " << endl;
            wave.print_atom_long(0);
            file << "ended normalizing the basis set, now for the MO_coeffs" << endl;
            file << "Status report:" << endl;
            file << "size of norm_const: " << norm_const.size() << endl;
            file << "WFN MO counter: " << wave.get_nmo() << endl;
            file << "Number of atoms: " << wave.get_ncen() << endl;
            file << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
            file << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
        }

        //-------------------normalize the basis set shell wise into a copy vector---------
        vec2 basis_coefficients(wave.get_ncen());
#pragma omp parallel for
        for (int a = 0; a < wave.get_ncen(); a++)
        {
            for (int p = 0; p < wave.get_atom_primitive_count(a); p++)
            {
                double temp_c = wave.get_atom_basis_set_exponent(a, p);
                switch (wave.get_atom_primitive_type(a, p))
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
                    file << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy..." << endl;
                    break;
                }
                temp_c = pow(temp_c, 0.25) * wave.get_atom_basis_set_coefficient(a, p);
                basis_coefficients[a].push_back(temp_c);
            }
        }
        for (int a = 0; a < wave.get_ncen(); a++)
        {
            double aiaj = 0.0;
            double factor = 0.0;
            for (int s = 0; s < wave.get_atom_shell_count(a); s++)
            {
                int type_temp = wave.get_shell_type(a, s);
                err_chkf(type_temp != -1, "ERROR in type assignement!!", file);
                if (debug)
                {
                    file << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << endl
                         << "start: " << wave.get_shell_start(a, s)
                         << " stop: " << wave.get_shell_end(a, s) << endl
                         << "factor: ";
                }
                switch (type_temp)
                {
                case 1:
                    factor = 0;
                    for (int i = wave.get_shell_start(a, s); i <= wave.get_shell_end(a, s); i++)
                    {
                        for (int j = wave.get_shell_start(a, s); j <= wave.get_shell_end(a, s); j++)
                        {
                            aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
                            double term = constants::PI3 / pow(aiaj, 3);
                            term = pow(term, 0.5);
                            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
                        }
                    }
                    if (factor == 0)
                        return false;
                    factor = pow(factor, -0.5);
                    if (debug)
                        file << factor << endl;
                    for (int i = wave.get_shell_start(a, s); i <= wave.get_shell_end(a, s); i++)
                    {
                        if (debug)
                        {
                            file << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i)
                                 << " Contraction coefficient after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
                        }
                        // contraction_coefficients[a][i] = factor * wave.get_atom_basis_set_coefficient(a, i);
                        basis_coefficients[a][i] *= factor;
                        norm_const.emplace_back(basis_coefficients[a][i]);
                    }
                    break;
                case 2:
                    factor = 0;
                    for (int i = wave.get_shell_start(a, s); i <= wave.get_shell_end(a, s); i++)
                    {
                        for (int j = wave.get_shell_start(a, s); j <= wave.get_shell_end(a, s); j++)
                        {
                            aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
                            double term = constants::PI3 / (4 * pow(aiaj, 5));
                            term = pow(term, 0.5);
                            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
                        }
                    }
                    if (factor == 0)
                        return false;
                    factor = pow(factor, -0.5);
                    if (debug)
                        file << factor << endl;
                    for (int i = wave.get_shell_start(a, s); i <= wave.get_shell_end(a, s); i++)
                    {
                        if (debug)
                        {
                            file << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i)
                                 << " Contraction coefficient after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
                        }
                        // contraction_coefficients[a][i] = factor * wave.get_atom_basis_set_coefficient(a, i);
                        basis_coefficients[a][i] *= factor;
                        for (int k = 0; k < 3; k++)
                            norm_const.emplace_back(basis_coefficients[a][i]);
                    }
                    break;
                case 3:
                    factor = 0;
                    for (int i = wave.get_shell_start(a, s); i <= wave.get_shell_end(a, s); i++)
                    {
                        for (int j = wave.get_shell_start(a, s); j <= wave.get_shell_end(a, s); j++)
                        {
                            aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
                            double term = constants::PI3 / (16 * pow(aiaj, 7));
                            term = pow(term, 0.5);
                            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
                        }
                    }
                    if (factor == 0)
                        return false;
                    factor = (pow(factor, -0.5)) / sqrt(3);
                    if (debug)
                        file << factor << endl;
                    for (int i = wave.get_shell_start(a, s); i <= wave.get_shell_end(a, s); i++)
                    {
                        if (debug)
                        {
                            file << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i)
                                 << " Contraction coefficient after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
                        }
                        // contraction_coefficients[a][i] = factor * wave.get_atom_basis_set_coefficient(a, i);
                        basis_coefficients[a][i] *= factor;
                        for (int k = 0; k < 3; k++)
                            norm_const.emplace_back(basis_coefficients[a][i]);
                        for (int k = 0; k < 3; k++)
                            norm_const.emplace_back(sqrt(3) * basis_coefficients[a][i]);
                    }
                    break;
                case 4:
                    factor = 0;
                    for (int i = wave.get_shell_start(a, s); i <= wave.get_shell_end(a, s); i++)
                    {
                        for (int j = wave.get_shell_start(a, s); j <= wave.get_shell_end(a, s); j++)
                        {
                            aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
                            double term = constants::PI3 / (64 * pow((aiaj), 9));
                            term = pow(term, 0.5);
                            factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
                        }
                    }
                    if (factor == 0)
                        return false;
                    factor = pow(factor, -0.5) / sqrt(15);
                    if (debug)
                        file << factor << endl;
                    for (int i = wave.get_shell_start(a, s); i <= wave.get_shell_end(a, s); i++)
                    {
                        if (debug)
                        {
                            file << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i)
                                 << " Contraction coefficient after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
                        }
                        // contraction_coefficients[a][i] = factor * wave.get_atom_basis_set_coefficient(a, i);
                        basis_coefficients[a][i] *= factor;
                        for (int l = 0; l < 3; l++)
                            norm_const.emplace_back(basis_coefficients[a][i]);
                        for (int l = 0; l < 6; l++)
                            norm_const.emplace_back(sqrt(5) * basis_coefficients[a][i]);
                        norm_const.emplace_back(sqrt(15) * basis_coefficients[a][i]);
                    }
                    break;
                }
                if (debug)
                    file << "This shell has: " << wave.get_shell_end(a, s) - wave.get_shell_start(a, s) + 1 << " primitives" << endl;
            }
        }
        //-----------debug output---------------------------------------------------------
        if (debug)
        {
            file << "exemplary output of the first atom with all it's properties: " << endl;
            wave.print_atom_long(0);
            file << "ended normalizing the basis set, now for the norm_cprims" << endl;
            file << "Status report:" << endl;
            file << "size of norm_const: " << norm_const.size() << endl;
            file << "WFN MO counter: " << wave.get_nmo() << endl;
            file << "Number of atoms: " << wave.get_ncen() << endl;
            file << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
            file << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
        }
        //---------------------To not mix up anything start normalizing WFN_matrix now--------------------------
        int run = 0;
        vec2 changed_coefs;
        changed_coefs.resize(wave.get_nmo());
        if (debug)
        {
            file << "Opening norm_cprim!" << endl;
            ofstream norm_cprim("norm_prim.debug", ofstream::out);
            for (int m = 0; m < wave.get_nmo(); m++)
            {
                norm_cprim << m << ". MO:" << endl;
                changed_coefs[m].resize(wave.get_MO_primitive_count(m), 0.0);
                for (int p = 0; p < wave.get_MO_primitive_count(m); p++)
                {
                    changed_coefs[m][p] = wave.get_MO_coef(m, p) / norm_const[p];
                    if (m == 0)
                        file << p << ". primitive; " << m << ". MO "
                             << "norm nonst: " << norm_const[p]
                             << " temp after normalization: " << changed_coefs[m][p] << "\n";
                    norm_cprim << " " << changed_coefs[m][p] << endl;
                    run++;
                }
            }
            norm_cprim.flush();
            norm_cprim.close();
            file << "See norm_cprim.debug for the CPRIM vectors" << endl;
            file << "Total count in CPRIM: " << run << endl;
        }
        else
        {
#pragma omp parallel for
            for (int m = 0; m < wave.get_nmo(); m++)
            {
                changed_coefs[m].resize(wave.get_MO_primitive_count(m), 0.0);
                for (int p = 0; p < wave.get_MO_primitive_count(m); p++)
                {
                    changed_coefs[m][p] = wave.get_MO_coef(m, p) / norm_const[p];
                }
            }
        }
        //--------------Build CMO of alessandro from the first elements of each shell-------------
        int nao = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
        {
            for (int s = 0; s < wave.get_atom_shell_count(a); s++)
            {
                switch (wave.get_shell_type(a, s))
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
        for (int m = 0; m < wave.get_nmo(); m++)
        {
            int run_2 = 0;
            for (int a = 0; a < wave.get_ncen(); a++)
            {
                for (int s = 0; s < wave.get_atom_shell_count(a); s++)
                {
                    // if (debug) file << "Going to load the " << wave.get_shell_start_in_primitives(a, s) << ". value\n"l;
                    switch (wave.get_shell_type(a, s))
                    {
                    case 1:
                        CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s)]);
                        if (debug && wave.get_atom_shell_primitives(a, s) != 1 && m == 0)
                            file << "Pushing back 1 coefficient for S shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives! Shell start is: " << wave.get_shell_start(a, s) << endl;
                        break;
                    case 2:
                        for (int i = 0; i < 3; i++)
                            CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i]);
                        if (debug && wave.get_atom_shell_primitives(a, s) != 1 && m == 0)
                            file << "Pushing back 3 coefficients for P shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
                        break;
                    case 3:
                        for (int i = 0; i < 6; i++)
                            CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i]);
                        if (debug && wave.get_atom_shell_primitives(a, s) != 1 && m == 0)
                            file << "Pushing back 6 coefficient for D shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
                        break;
                    case 4:
                        // this hardcoded piece is due to the order of f-type functions in the fchk
                        for (int i = 0; i < 3; i++)
                            CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i]);
                        CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + 6]);
                        for (int i = 0; i < 2; i++)
                            CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i + 3]);
                        for (int i = 0; i < 2; i++)
                            CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i + 7]);
                        CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + 5]);
                        CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + 9]);
                        if (debug && wave.get_atom_shell_primitives(a, s) != 1 && m == 0)
                            file << "Pushing back 10 coefficient for F shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
                        break;
                    }
                    run_2++;
                }
                if (debug && m == 0)
                    file << "finished with atom!" << endl;
            }
            if (debug)
                file << "finished with MO!" << endl;
            if (nshell != run_2)
                nshell = run_2;
        }
        if (alpha_els != beta_els)
        {
            for (int m = alpha_els; m < alpha_els + beta_els; m++)
            {
                int run_2 = 0;
                for (int a = 0; a < wave.get_ncen(); a++)
                {
                    for (int s = 0; s < wave.get_atom_shell_count(a); s++)
                    {
                        if (debug)
                            file << "Going to load the " << wave.get_shell_start_in_primitives(a, s) << ". value" << endl;
                        switch (wave.get_shell_type(a, s))
                        {
                        case 1:
                            CMO_beta.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s)]);
                            if (m == 0)
                                nao++;
                            if (debug && wave.get_atom_shell_primitives(a, s) != 1)
                                file << "Pushing back 1 coefficient for S shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives! Shell start is: " << wave.get_shell_start(a, s) << endl;
                            break;
                        case 2:
                            for (int i = 0; i < 3; i++)
                                CMO_beta.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i]);
                            if (debug && wave.get_atom_shell_primitives(a, s) != 1)
                                file << "Pushing back 3 coefficients for P shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
                            if (m == 0)
                                nao += 3;
                            break;
                        case 3:
                            for (int i = 0; i < 6; i++)
                                CMO_beta.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i]);
                            if (debug && wave.get_atom_shell_primitives(a, s) != 1)
                                file << "Pushing back 6 coefficient for D shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
                            if (m == 0)
                                nao += 6;
                            break;
                        case 4:
                            // this hardcoded piece is due to the order of f-type functions in the fchk
                            for (int i = 0; i < 3; i++)
                                CMO_beta.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i]);
                            CMO_beta.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + 6]);
                            for (int i = 0; i < 2; i++)
                                CMO_beta.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i + 3]);
                            for (int i = 0; i < 2; i++)
                                CMO_beta.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i + 7]);
                            CMO_beta.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + 5]);
                            CMO_beta.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + 9]);
                            if (debug && wave.get_atom_shell_primitives(a, s) != 1)
                                file << "Pushing back 10 coefficient for F shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
                            if (m == 0)
                                nao += 10;
                            break;
                        }
                        run_2++;
                    }
                    if (debug)
                        file << "finished with atom!" << endl;
                }
                if (debug)
                    file << "finished with MO!" << endl;
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
            file << CMO.size() << " Elements in CMO" << endl;
            file << norm_const.size() << " = nprim" << endl;
            file << nao << " = nao" << endl;
            file << nshell << " = nshell" << endl;
        }
        //------------------ make the DM -----------------------------
        int naotr = nao * (nao + 1) / 2;
        vec kp;
        wave.resize_DM(naotr, 0.0);
        if (alpha_els != beta_els)
            wave.resize_SDM(naotr, 0.0);
        if (debug)
        {
            file << "I made kp!" << endl
                 << nao << " is the maximum for iu" << endl;
           std::cout << "Making DM now!" << endl;
        }
        for (int iu = 0; iu < nao; iu++)
        {
#pragma omp parallel for
            for (int iv = 0; iv <= iu; iv++)
            {
                const int iuv = (iu * (iu + 1) / 2) + iv;
                // if (debug) file << "iu: " << iu << " iv: " << iv << " iuv: " << iuv << " kp(iu): " << iu * (iu + 1) / 2 << endl;
                double temp;
                // if (debug) file << "Working on MO: ";
                for (int m = 0; m < wave.get_nmo(); m++)
                {
                    // if (debug && m == 0) file << m << " " << flush;
                    // else if (debug && m != wave.get_nmo() - 1) file << "." << flush;
                    // else file << wave.get_nmo() - 1 << flush;
                    if (alpha_els != beta_els)
                    {
                        if (m < alpha_els)
                        {
                            temp = wave.get_MO_occ(m) * CMO[iu + (m * nao)] * CMO[iv + (m * nao)];
                            err_checkf(wave.set_SDM(iuv, wave.get_SDM(iuv) + temp), "Something went wrong while writing the SDM! iuv=" + to_string(iuv), file);
                            err_checkf(wave.set_DM(iuv, wave.get_DM(iuv) + temp), "Something went wrong while writing the DM! iuv=" + to_string(iuv), file);
                        }
                        else
                        {
                            temp = wave.get_MO_occ(m) * CMO_beta[iu + ((m - alpha_els) * nao)] * CMO_beta[iv + ((m - alpha_els) * nao)];
                            err_checkf(wave.set_SDM(iuv, wave.get_SDM(iuv) - temp), "Something went wrong while writing the SDM! iuv=" + to_string(iuv), file);
                            err_checkf(wave.set_DM(iuv, wave.get_DM(iuv) + temp), "Something went wrong while writing the DM! iuv=" + to_string(iuv), file);
                        }
                    }
                    else
                    {
                        if (wave.get_MO_occ(m) == 0.0)
                            continue;
                        temp = wave.get_MO_occ(m) * CMO[iu + (m * nao)] * CMO[iv + (m * nao)];
                        err_checkf(wave.set_DM(iuv, wave.get_DM(iuv) + temp), "Something went wrong while writing the DM!", file);
                    }
                    // else if (debug) file << "DM after: " << wave.get_DM(iuv) << endl;
                }
                // if (debug) file << endl;
            }
        }
        if (debug)
        {
           std::cout << "Done with DM!" << endl;
            ofstream dm("dm.debug", ofstream::out);
            file << "DM is in dm.debug" << endl;
            for (int p = 0; p < wave.get_DM_size(); p++)
            {
                for (int i = 0; i < 5; i++)
                {
                    if (i + p >= wave.get_DM_size())
                        continue;
                    dm << scientific << setw(14) << setprecision(7) << wave.get_DM(i + p) << " ";
                }
                p += 4;
                dm << endl;
            }
            file << wave.get_DM_size() << " Elements in DM" << endl;
            dm.flush();
            dm.close();
            if (alpha_els != beta_els)
            {
                ofstream sdm("sdm.debug", ofstream::out);
                file << "SDM is in sdm.debug" << endl;
                for (int p = 0; p < wave.get_SDM_size(); p++)
                {
                    for (int i = 0; i < 5; i++)
                    {
                        if (i + p >= wave.get_DM_size())
                            continue;
                        sdm << scientific << setw(14) << setprecision(7) << wave.get_SDM(i + p) << " ";
                    }
                    p += 4;
                    sdm << endl;
                }
                file << wave.get_SDM_size() << " Elements in SDM" << endl;
                sdm.flush();
                sdm.close();
            }
            file << "Starting to write fchk now!" << endl;
        }
        // open fchk for writing
        std::filesystem::path temp_fchk = fchk_name;
        temp_fchk.replace_extension(".fchk");
        file << "Writing " << temp_fchk << " ..." << flush;
        if (std::filesystem::exists(temp_fchk) && !force_overwrite)
        {
            file << "The fchk already exists!" << endl;
            return false;
        }
        ofstream fchk(temp_fchk.c_str());
        if (!fchk.is_open())
        {
            file << "ERROR while opening .fchk ifile!" << endl;
            return false;
        }
        // stringstream st_s;
        // string s;
        // st_s.str("");
        fchk << "TITLE\n";

        if (wave.get_method() == "rhf" && wave.get_multi() == 1)
            fchk << "SP        RHF                                                         Gen\n";
        else if (wave.get_method() == "rks" && wave.get_multi() == 1)
            fchk << "SP        B3LYP                                                       Gen\n";
        else if (wave.get_method() == "rhf" && wave.get_multi() > 1)
            fchk << "SP        UHF                                                         Gen\n";
        else if (wave.get_method() == "rks" && wave.get_multi() > 1)
            fchk << "SP        UB3LYP                                                      Gen\n";
        else
            fchk << "SP        B3LYP                                                       Gen\n";
        fchk << "Number of atoms                            I" << setw(17) << wave.get_ncen();
        fchk << "\nInfo1-9                                    I   N=           9\n"
             << "          53          51           0           0           0         111\n           1           1           2"; // I have NO CLUE what the hell this means...
        fchk << "\nCharge                                     I" << setw(17) << wave.get_charge()
             << "\nMultiplicity                               I" << setw(17) << wave.get_multi();
        fchk << "\nNumber of electrons                        I" << setw(17) << elcount
             << "\nNumber of alpha electrons                  I" << setw(17) << alpha_els;
        fchk << "\nNumber of beta electrons                   I" << setw(17) << beta_els
             << "\nNumber of basis functions                  I" << setw(17) << nao;
        fchk << "\nNumber of independent functions            I" << setw(17) << nao
             << "\nNumber of point charges in /Mol/           I                0\nNumber of translation vectors              I                0\n";
        fchk << "Atomic numbers                             I   N=" << setw(12) << wave.get_ncen() << "\n";
        for (int i = 0; i < wave.get_ncen(); i++)
        {
            fchk << setw(12) << wave.get_atom_charge(i);
            if (((i + 1) % 6 == 0 && i != 0) || i == wave.get_ncen() - 1)
                fchk << "\n";
        }
        fchk << "Nuclear charges                            R   N=" << setw(12) << wave.get_ncen() << "\n";
        for (int i = 0; i < wave.get_ncen(); i++)
        {
            fchk << uppercase << scientific << setw(16) << setprecision(8) << static_cast<double>(wave.get_atom_charge(i));
            if (((i + 1) % 5 == 0 && i != 0) || i == wave.get_ncen() - 1)
                fchk << "\n";
        }
        fchk << "Current cartesian coordinates              R   N=" << setw(12) << wave.get_ncen() * 3 << "\n";
        unsigned int runs = 0;
        for (int i = 0; i < wave.get_ncen(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_coordinate(i, j);
                runs++;
                if (runs % 5 == 0 || (i == wave.get_ncen() - 1 && j == 2))
                    fchk << "\n";
            }
        }
        fchk << "Integer atomic weights                     I   N=" << setw(12) << wave.get_ncen() << "\n";
        for (int i = 0; i < wave.get_ncen(); i++)
        {
            fchk << setw(12) << wave.get_atom_integer_mass((unsigned int)i);
            if (((i + 1) % 6 == 0 && i != 0) || i == wave.get_ncen() - 1)
                fchk << "\n";
        }
        fchk << "Real atomic weights                        R   N=" << setw(12) << wave.get_ncen() << "\n";
        for (int i = 0; i < wave.get_ncen(); i++)
        {
            fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_real_mass(i);
            if (((i + 1) % 5 == 0 && i != 0) || i == wave.get_ncen() - 1)
                fchk << "\n";
        }

        unsigned int n_contracted_shells = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
            for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
                n_contracted_shells++;
        fchk << "Number of contracted shells                I" << setw(17) << n_contracted_shells;
        unsigned int n_prim_shells = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
            for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
                for (int p = 0; p < wave.get_atom_shell_primitives(a, sh); p++)
                    n_prim_shells++;
        fchk << "\nNumber of primitive shells                 I" << setw(17) << n_prim_shells
             << "\nPure/Cartesian d shells                    I" << setw(17) << 1;
        int d_count = 0;
        int f_count = 0;
        int max_contraction = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
            for (int _s = 0; _s < wave.get_atom_shell_count(a); _s++)
            {
                if (wave.get_atom_shell_primitives(a, _s) > max_contraction)
                    max_contraction = wave.get_atom_shell_primitives(a, _s);
                if (wave.get_shell_type(a, _s) == 3)
                    d_count++;
                else if (wave.get_shell_type(a, _s) == 4)
                    f_count++;
            }
        fchk << "\nPure/Cartesian f shells                    I" << setw(17) << 1
             << "\nHighest angular momentum                   I" << setw(17) << 1 + (d_count > 0) + (f_count > 0)
             << "\nLargest degree of contraction              I" << setw(17) << max_contraction
             << "\nShell types                                I   N=" << setw(12) << n_contracted_shells << "\n";
        runs = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
            for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
            {
                fchk << setw(12) << wave.get_shell_type(a, sh) - 1;
                runs++;
                if ((runs % 6 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1))
                    fchk << "\n";
            }
        fchk << "Number of primitives per shell             I   N=" << setw(12) << n_contracted_shells << "\n";
        runs = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
            for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
            {
                fchk << setw(12) << wave.get_atom_shell_primitives(a, sh);
                runs++;
                if ((runs % 6 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1))
                    fchk << "\n";
            }
        fchk << "Shell to atom map                          I   N=" << setw(12) << n_contracted_shells << "\n";
        runs = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
            for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
            {
                fchk << setw(12) << wave.get_shell_center(a, sh);
                runs++;
                if ((runs % 6 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1))
                    fchk << "\n";
            }

        fchk << "Primitive exponents                        R   N=" << setw(12) << n_prim_shells << "\n";
        runs = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
        {
            int p_run = 0;
            for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
                for (int p = 0; p < wave.get_atom_shell_primitives(a, sh); p++)
                {
                    fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_basis_set_exponent(a, p_run);
                    runs++;
                    p_run++;
                    if ((runs % 5 == 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1 && p == wave.get_atom_shell_primitives(a, sh) - 1))
                        fchk << "\n";
                }
        }
        fchk << "Contraction coefficients                   R   N=" << setw(12) << n_prim_shells << "\n";
        runs = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
        {
            int p_run = 0;
            for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
                for (int p = 0; p < wave.get_atom_shell_primitives(a, sh); p++)
                {
                    fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_basis_set_coefficient(a, p_run);
                    runs++;
                    p_run++;
                    if ((runs % 5 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1 && p == wave.get_atom_shell_primitives(a, sh) - 1))
                        fchk << "\n";
                }
        }
        fchk << "Coordinates of each shell                  R   N=" << setw(12) << 3 * n_contracted_shells << "\n";
        runs = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
            for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
                for (int i = 0; i < 3; i++)
                {
                    fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_coordinate(a, i);
                    runs++;
                    if ((runs % 5 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1 && i == 2))
                        fchk << "\n";
                }
        fchk << "Constraint Structure                       R   N=" << setw(12) << 3 * wave.get_ncen() << "\n";
        runs = 0;
        for (int a = 0; a < wave.get_ncen(); a++)
            for (int i = 0; i < 3; i++)
            {
                fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_coordinate(a, i);
                runs++;
                if ((runs % 5 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && i == 2))
                    fchk << "\n";
            }
        double energy = 0;
        for (int m = 0; m < wave.get_nmo(); m++)
            if (wave.get_MO_occ(m) > 0)
                energy += wave.get_MO_energy(m);
        fchk << "Virial Ratio                               R" << uppercase << scientific << setw(27) << setprecision(15) << 2.0
             << "\nSCF Energy                                 R" << uppercase << scientific << setw(27) << setprecision(15) << energy
             << "\nTotal Energy                               R" << uppercase << scientific << setw(27) << setprecision(15) << energy
             << "\nRMS Density                                R" << uppercase << scientific << setw(27) << setprecision(15) << 2.0 * pow(10, -9)
             << "\nAlpha Orbital Energies                     R   N=" << setw(12) << nao << "\n";
        runs = 0;
        for (int m = 0; m < nao; m++)
        {
            if (m < wave.get_nmo())
                fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_MO_energy(m);
            else
                fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_MO_energy(alpha_els - 1) + m;
            runs++;
            if ((runs % 5 == 0 && runs != 0) || m == nao - 1)
                fchk << "\n";
        }

        if (alpha_els != beta_els)
        {
            fchk << "Beta Orbital Energies                      R   N=" << setw(12) << nao << "\n";
            runs = 0;
            for (int m = 0; m < nao; m++)
            {
                if (m + nao < wave.get_nmo())
                    fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_MO_energy(m + nao);
                else
                    fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_MO_energy(wave.get_nmo() - 1) + m;
                runs++;
                if ((runs % 5 == 0 && runs != 0) || m == nao - 1)
                    fchk << "\n";
            }
        }

        fchk << "Alpha MO coefficients                      R   N=" << setw(12) << nao * nao << "\n";
        runs = 0;
        for (int i = 0; i < nao * nao; i++)
        {
            if (i < nao * wave.get_nmo())
                fchk << uppercase << scientific << setw(16) << setprecision(8) << CMO[i];
            else
                fchk << uppercase << scientific << setw(16) << setprecision(8) << 0.0;
            runs++;
            if ((runs % 5 == 0 && runs != 0) || i == nao * nao - 1)
                fchk << "\n";
        }

        if (alpha_els != beta_els)
        {
            fchk << "Beta MO coefficients                       R   N=" << setw(12) << nao * nao << "\n";
            runs = 0;
            for (int i = 0; i < nao * nao; i++)
            {
                if (i < nao * beta_els)
                    fchk << uppercase << scientific << setw(16) << setprecision(8) << CMO_beta[i];
                else
                    fchk << uppercase << scientific << setw(16) << setprecision(8) << 0.0;
                runs++;
                if ((runs % 5 == 0 && runs != 0) || i == nao * nao - 1)
                    fchk << "\n";
            }
        }

        fchk << "Total SCF Density                          R   N=" << setw(12) << wave.get_DM_size() << "\n";
        runs = 0;
        for (int i = 0; i < wave.get_DM_size(); i++)
        {
            fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_DM(i);
            runs++;
            if ((runs % 5 == 0 && runs != 0) || i == wave.get_DM_size() - 1)
                fchk << "\n";
        }
        if (alpha_els != beta_els)
        {
            fchk << "Spin SCF Density                           R   N=" << setw(12) << wave.get_SDM_size() << "\n";
            runs = 0;
            for (int i = 0; i < wave.get_SDM_size(); i++)
            {
                fchk << uppercase << scientific << setw(16) << setprecision(8) << wave.get_SDM(i);
                runs++;
                if ((runs % 5 == 0 && runs != 0) || i == wave.get_SDM_size() - 1)
                    fchk << "\n";
            }
        }
        fchk.flush();
        fchk.close();
        file << " ...done!" << endl;
    }
    if (debug)
        file << "Finished writing fchk!" << endl;
    return true;
};
