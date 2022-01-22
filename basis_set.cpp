#include "basis_set.h"
#include "convenience.h"
struct basis_set_entry;
class WFN;
using namespace std;
bool debug_dens = false;

//-------------Reading basis sets and determining the density matrix from the wfn coefficients -----------------------

//test
//read the basis set and give the info to the atoms of the wavefunction
bool read_basis_set(const string& basis_set_path, WFN& wave, bool debug)
{
  if (debug) debug_dens = true;
  string basis_set_name;
  string temp;
  bool end = false;
  bool manual = false;
  while (!end) {
    if (wave.get_basis_set_name().length() < 3) {
      cout << "Which basis set do you want to load?" << endl;
      cin >> basis_set_name;
    }
    else if (!manual) basis_set_name = wave.get_basis_set_name();
    //assemble basis set name and look if file exists
    temp = basis_set_path;
    temp.append(basis_set_name);
    if (exists(temp)) {
      if (debug_dens) cout << "basis set is valid, continueing..." << endl;
      end = true;
    }
    else {
      cout << "sorry, could not find this basis set in the basis set directory specified in the programs.config file!" << endl
        << "Do you want to speify a new name?" << endl;
      if (!yesno()) return false;
      manual = true;
      cout << "What is the name of the basis set in the directory: ";
      cin >> basis_set_name;
    }
  }
  if (debug_dens) cout << "File of basis set to load: " << temp << endl;
  ifstream ifile(temp.c_str(), ios::in);
  //  Looking for all the types of atoms we need to find
  vector<string> elements_list;
  bool found = false;
  for (int i = 0; i < wave.get_ncen(); i++) {
    if (debug_dens) cout << "i: " << i << endl;
    for (int j = 0; j < elements_list.size(); j++) {
      if (elements_list[j].compare(wave.get_atom_label(i)) == 0) found = true;
      if (debug_dens) cout << "   j: " << j << " Atom label: " << wave.get_atom_label(i) << endl;
    }
    if (!found) {
      elements_list.push_back(wave.get_atom_label(i));
      if (debug_dens) cout << "Added an atom which was not there yet!" << wave.get_atom_label(i) << endl;
    }
    found = false;
  }
  if (debug_dens) {
    cout << "Number of elements in elements_list: " << elements_list.size() << endl;
    cout << "This is the elements list:" << endl;
    for (int l = 0; l < elements_list.size(); l++) cout << l << ": " << elements_list[l] << endl;
    Enter();
  }
  int found_counter = 0;
  for (int i = 0; i < elements_list.size(); i++) {
    if (debug_dens) cout << "before: " << elements_list[i] << " " << i << endl;
    while (elements_list[i].find(" ") != -1) elements_list[i].erase(elements_list[i].find(" "), 1);
    elements_list[i].append(":");
    if (debug_dens) {
      cout << "after: " << elements_list[i] << " " << i << endl;
      Enter();
    }
    //scan the tonto style basis set file for the entries we are looking or:
    string line;
    getline(ifile, line);
    int file_type = 0;
    //check if we support that type of basis set
    while (line.find("keys=") == string::npos && !ifile.eof())
      getline(ifile, line);
    if (debug_dens) {
      cout << "Line after looking for keys=: " << line << endl;
    }
    if (line.find("keys=") < line.size() && debug_dens) cout << "Found keys=!" << endl;
    if (line.find("turbomole") < line.size()) {
      file_type = 1;
      if (debug_dens) cout << "This file is written in turbomole type!" << endl;
    }
    else if (line.find("gamess-us") < line.size()) {
      file_type = 2;
      if (debug_dens) cout << "This file is written in gamess-us type!" << endl;
    }
    else if (line.find("gaussian") < line.size()) {
      file_type = 3;
      if (debug_dens) cout << "This file is written in gaussian type!" << endl;
    }
    else if (line.find("CRYSTAL") < line.size()) {
      file_type = 1;
      wave.set_d_f_switch(true);
      if (debug_dens) cout << "This file is written in CRYSTAL type!" << endl;
    }
    else {
      cout << "This type of file is not supported, please provide another basis set!" << endl;
      return false;
    }
    if (ifile.eof()) {
      cout << "Please provide a basis set in the turbomole, gaussian or gamess-us format compatible with tonto."
        << "Look at the example files \"examble.basis\" and \"examble2.basis\" in the wfn_cpp folder if you want to see how it has to look like" << endl;
      return false;
    }
    while (!(line.find(elements_list[i]) < line.size()) && !ifile.eof()) {
      getline(ifile, line);
      if (debug_dens) cout << "line while search for " << elements_list[i] << " :" << line << endl;
    }
    if (debug_dens && line.find(elements_list[i]) != -1) {
      cout << "I found an entry i know from the element list!" << endl;
      cout << "The line is: " << line << endl;
      Enter();
    }
    if (ifile.eof()) {
      cout << "Could not find the atom you were looking for in the basis set file... " << endl;
      Enter();
      return false;
    }
    unsigned int shell = 0;
    if (line.find("{") == -1) {
      getline(ifile, line);
      if (debug) {
        cout << "I read an additional line!" << endl;
        //Enter();
      }
    }
    while (line.find("}") == string::npos && !ifile.eof()) {
      getline(ifile, line);
      stringstream stream;
      stream << line;
      if (line.find("}") != string::npos)
        break;
      int count = 0;
      int nr_exp = 0;
      double pi = 3.14159265358979;
      char c_temp = '?';
      double temp_vals[2];
      int dum = 0;
      if (file_type == 1) {
        stream >> count >> c_temp;
        if (debug_dens) cout << "count: " << count << " type: " << c_temp << endl;
      }
      else if (file_type == 2) {
        stream >> c_temp >> count;
        if (debug_dens) cout << "count: " << count << " type: " << c_temp << endl;
      }
      else if (file_type == 3) {
        stream >> c_temp >> count;
        if (debug_dens) cout << "count: " << count << " type: " << c_temp << endl;
      }
      for (int j = 0; j < count; j++) {
        getline(ifile, line);
        if (debug_dens) {
          cout << "read the " << j << ". line: " << line << endl;
        }
        stringstream stream2;
        stream2 << line;
        if (file_type == 1) stream2 >> temp_vals[0] >> temp_vals[1];
        else if (file_type == 2) stream2 >> dum >> temp_vals[0] >> temp_vals[1];
        else if (file_type == 3) stream2 >> temp_vals[0] >> temp_vals[1];
        switch (c_temp) {
        case 's':
        case 'S':
          temp_vals[1] *= pow((2 * temp[0] / pi), 0.75);
          break;
        case 'p':
        case 'P':
          temp_vals[1] *= pow((128 * pow(temp[0], 5) / pow(pi, 3)), 0.25);
          break;
        case 'd':
        case 'D':
          temp_vals[1] *= (pow(2048 * pow(temp[0], 7) / (9 * pow(pi, 3)), 0.25));
          break;
        case 'f':
        case 'F':
          temp_vals[1] *= pow((32768 * pow(temp[0], 9) / (225 * pow(pi, 3))), 0.25);
          break;
        default:
          cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
          return false;
        } // end switch of types
        //this is where i started copying
        for (int h = 0; h < wave.get_ncen(); h++) {
          string temp_label;
          temp_label = wave.get_atom_label(h);
          while (temp_label.find(" ") != -1) temp_label.erase(temp_label.find(" "), 1);
          temp_label.append(":");
          if (elements_list[i].find(temp_label) != -1) {
            if (debug_dens) {
              cout << "It's a match!" << endl;
              cout << "element_label: " << elements_list[i] << " temp_label: " << temp_label << endl;
            }
            switch (c_temp) {
            case 's':
            case 'S':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 1, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type S" << endl;
              break;
            case 'p':
            case 'P':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 2, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type P" << endl;
              break;
            case 'd':
            case 'D':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 3, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type D" << endl;
              break;
            case 'f':
            case 'F':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 4, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type F" << endl;
              break;
            default:
              cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
              return false;
            } // end switch of types
          } // end if(find atom_label + : 
        }//end for h = ncen
        nr_exp++;
        if (debug_dens) cout << "recapitulation[" << j << "]... type: " << c_temp << " coef: " << temp[0] << " exp: " << temp[1] << endl;
        if (dum > count) {
          cout << "this should not happen, lets stop before i do something silly!" << endl;
          return false;
        }
      }
      shell++;
    }//end while line != }
    if (debug_dens) cout << "I found }: " << line << endl;
    shell = 0;
    ifile.seekg(0);
  }//end for element_list.size()
  if (debug_dens) {
    cout << "FINISHED WITH READING BASIS SET!" << endl;
    Enter();
  }
  ifile.close();
  return true;
};

bool read_basis_set_vanilla(const string& basis_set_path, WFN& wave, const bool& debug, bool manual)
{
  if (debug) debug_dens = true;
  string basis_set_name;
  string temp;
  bool end = false;
  while (!end) {
    if (wave.get_basis_set_name().length() < 3 || manual) {
      cout << "Which basis set do you want to load?" << endl;
      cin >> basis_set_name;
    }
    else if (!manual) basis_set_name = wave.get_basis_set_name();
    //assemble basis set name and look if file exists
    temp = basis_set_path;
    temp.append(basis_set_name);
    if (exists(temp)) {
      if (debug_dens) cout << "basis set is valid, continuing..." << endl;
      end = true;
    }
    else {
      cout << "sorry, could not find this basis set in the basis set directory specified in the ~/.cuQCT file!" << endl
        << "Do you want to specify a new name?" << endl;
      if (!yesno()) return false;
      manual = true;
      cout << "What is the name of the basis set in the directory: ";
      cin >> basis_set_name;
    }
  }
  wave.set_basis_set_name(basis_set_name);
  if (debug_dens) cout << "File of basis set to load: " << temp << endl;
  ifstream ifile(temp.c_str(), ios::in);
  //  Looking for all the types of atoms we need to find
  vector<string> elements_list;
  bool found = false;
  for (int i = 0; i < wave.get_ncen(); i++) {
    if (debug_dens) cout << "i: " << i << endl;
    for (int j = 0; j < elements_list.size(); j++) {
      if (elements_list[j].compare(wave.get_atom_label(i)) == 0) found = true;
      if (debug_dens) cout << "   j: " << j << " Atom label: " << wave.get_atom_label(i) << endl;
    }
    if (!found) {
      elements_list.push_back(wave.get_atom_label(i));
      if (debug_dens) cout << "Added an atom which was not there yet! " << wave.get_atom_label(i) << "!" << endl;
    }
    found = false;
  }
  if (debug_dens) {
    cout << "Number of elements in elements_list: " << elements_list.size() << endl;
    cout << "This is the elements list:" << endl;
    for (int l = 0; l < elements_list.size(); l++) cout << l << ": " << elements_list[l] << "," << endl;
    Enter();
  }
  int found_counter = 0;
  for (int i = 0; i < elements_list.size(); i++) {
    if (debug_dens) cout << "before: " << elements_list[i] << " " << i << endl;
    while (elements_list[i].find(" ") != -1) { elements_list[i].erase(elements_list[i].find(" "), 1); } //Cut out spaces
    elements_list[i].append(":");
    if (debug_dens) {
      cout << "after: " << elements_list[i] << " " << i << endl;
      Enter();
    }
    //scan the tonto style basis set file for the entries we are looking or:
    string line;
    getline(ifile, line);
    int file_type = 0;
    //check if we support that type of basis set
    while (line.find("keys=") == -1 && !ifile.eof())
      getline(ifile, line);
    if (debug_dens) {
      cout << "Line after looking for keys=: " << line << endl;
    }
    if (line.find("keys=") < line.size() && debug_dens) cout << "Found keys=!" << endl;
    if (line.find("turbomole") < line.size()) {
      file_type = 1;
      if (debug_dens) cout << "This file is written in turbomole type!" << endl;
    }
    else if (line.find("gamess-us") < line.size()) {
      file_type = 2;
      if (debug_dens) cout << "This file is written in gamess-us type!" << endl;
    }
    else if (line.find("gaussian") < line.size()) {
      file_type = 3;
      if (debug_dens) cout << "This file is written in gaussian type!" << endl;
    }
    else if (line.find("CRYSTAL") < line.size()) {
      file_type = 1;
      wave.set_d_f_switch(true);
      if (debug_dens) cout << "This file is written in CRYSTAL type!" << endl;
    }
    else {
      cout << "This type of file is not supported, please provide another basis set!" << endl;
      return false;
    }
    if (ifile.eof()) {
      cout << "Please provide a basis set in the turbomole, gaussian or gamess-us format compatible with tonto."
        << "Look at the example files \"examble.basis\" and \"examble2.basis\" in the wfn_cpp folder if you want to see how it has to look like" << endl;
      return false;
    }
    while (!(line.find(elements_list[i]) < line.size()) && !ifile.eof())
      getline(ifile, line);
    if (debug_dens) cout << "line while search for " << elements_list[i] << " :" << line << endl;
    if (debug_dens && line.find(elements_list[i]) != -1) {
      cout << "I found an entry i know from the element list!" << endl;
      cout << "The line is: " << line << endl;
      Enter();
    }
    if (ifile.eof()) {
      cout << "Could not find the atom you were looking for in the basis set file... " << endl;
      Enter();
      return false;
    }
    unsigned int shell = 0;
    if (line.find("{") == -1) {
      getline(ifile, line);
      if (debug) {
        cout << "I read an additional line!" << endl;
        //Enter();
      }
    }
    while (line.find("}") == string::npos && !ifile.eof()) {
      getline(ifile, line);
      stringstream stream;
      stream << line;
      if (line.find("}") != string::npos)
        break;
      int count = 0;
      int nr_exp = 0;
      char c_temp = '?';
      double temp_vals[2]{ 0,0 };
      int dum = 0;
      if (file_type == 1) {
        stream >> count >> c_temp;
        if (debug_dens) cout << "count: " << count << " type: " << c_temp << endl;
      }
      else if (file_type == 2) {
        stream >> c_temp >> count;
        if (debug_dens) cout << "count: " << count << " type: " << c_temp << endl;
      }
      else if (file_type == 3) {
        stream >> c_temp >> count;
        if (debug_dens) cout << "count: " << count << " type: " << c_temp << endl;
      }
      for (int j = 0; j < count; j++) {
        getline(ifile, line);
        if (debug_dens) {
          cout << "read the " << j << ". line: " << line << endl;
        }
        stringstream stream2;
        stream2 << line;
        if (file_type == 1) stream2 >> temp_vals[0] >> temp_vals[1];
        else if (file_type == 2) stream2 >> dum >> temp_vals[0] >> temp_vals[1];
        else if (file_type == 3) stream2 >> temp_vals[0] >> temp_vals[1];
        //this is where i started copying
        for (int h = 0; h < wave.get_ncen(); h++) {
          string temp_label;
          temp_label = wave.get_atom_label(h);
          while (temp_label.find(" ") != -1)
            temp_label.erase(temp_label.find(" "), 1);
          temp_label.append(":");
          if (elements_list[i].find(temp_label) != -1) {
            switch (c_temp) {
            case 's':
            case 'S':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 1, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type S" << endl;
              break;
            case 'p':
            case 'P':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 2, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type P" << endl;
              break;
            case 'd':
            case 'D':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 3, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type D" << endl;
              break;
            case 'f':
            case 'F':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 4, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type F" << endl;
              break;
            default:
              cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
              return false;
            } // end switch of types
          } // end if(find atom_label + : 
        }//end for h = ncen
        nr_exp++;
        if (debug_dens) cout << "recapitulation[" << j << "]... type: " << c_temp << " coef: " << temp[0] << " exp: " << temp[1] << endl;
        if (dum > count) {
          cout << "this should not happen, lets stop before i do something silly!" << endl;
          return false;
        }
      }
      shell++;
    }//end while line != }
    if (debug_dens) cout << "I found }: " << line << endl;
    shell = 0;
    ifile.seekg(0);
  }//end for element_list.size()
  if (debug_dens) {
    cout << "FINISHED WITH READING BASIS SET!" << endl;
    Enter();
  }
  ifile.close();
  return true;
};

bool read_basis_set_missing(const string& basis_set_path, WFN& wave, bool debug)
{
  if (debug) debug_dens = true;
  string basis_set_name;
  string temp;
  bool end = false;
  bool manual = false;
  while (!end) {
    if (wave.get_basis_set_name().length() < 3) {
      cout << "Which basis set do you want to load?" << endl;
      cin >> basis_set_name;
    }
    else if (!manual) basis_set_name = wave.get_basis_set_name();
    //assemble basis set name and look if file exists
    temp = basis_set_path;
    temp.append(basis_set_name);
    if (exists(temp)) {
      if (debug_dens) cout << "basis set is valid, continueing..." << endl;
      end = true;
    }
    else {
      cout << "sorry, could not find this basis set in the basis set directory specified in the programs.config file!" << endl
        << "Do you want to speify a new name?" << endl;
      if (!yesno()) return false;
      manual = true;
      cout << "What is the name of the basis set in the directory: ";
      cin >> basis_set_name;
    }
  }
  if (debug_dens) cout << "File of basis set to load: " << temp << endl;
  ifstream ifile(temp.c_str(), ios::in);
  //  Looking for all the types of atoms we need to find
  vector<string> elements_list;
  bool found = false;
  for (int i = 0; i < wave.get_ncen(); i++) {
    if (debug_dens) cout << "i: " << i << endl;
    for (int j = 0; j < elements_list.size(); j++) {
      if (elements_list[j].find(wave.get_atom_label(i)) != -1) found = true;
      if (debug_dens) cout << "   j: " << j << endl;
    }
    if (!found) {
      elements_list.push_back(wave.get_atom_label(i));
      if (debug_dens) cout << "Added an atom which was not there yet! " << wave.get_atom_label(i) << endl;
    }
    found = false;
  }
  if (debug_dens) {
    cout << "Number of elements in elements_list: " << elements_list.size() << endl;
    cout << "This is the elements list:" << endl;
    for (int l = 0; l < elements_list.size(); l++) cout << l << ": " << elements_list[l] << endl;
    Enter();
  }
  int found_counter = 0;
  for (int i = 0; i < elements_list.size(); i++) {
    if (debug_dens) cout << "before: " << elements_list[i] << " " << i << endl;
    if (elements_list[i].find(" ")) elements_list[i].erase(elements_list[i].find(" "), 1);
    elements_list[i].append(":");
    if (debug_dens) {
      cout << "after: " << elements_list[i] << " " << i << endl;
      Enter();
    }
    //scan the tonto style basis set file for the entries we are looking or:
    string line;
    getline(ifile, line);
    int file_type = 0;
    //check if we support that type of basis set
    while (line.find("keys=") == -1 && !ifile.eof()) {
      if (debug_dens) {
        cout << "line.size of first line: " << line.size() << "line.find(\"keys=\"): " << line.find("keys=") << endl;
      }
      getline(ifile, line);
    }
    if (debug_dens) {
      cout << "Line after looking for keys=: " << line << endl;
    }
    if (line.find("keys=") < line.size() && debug_dens) cout << "Found keys=!" << endl;
    if (line.find("turbomole") < line.size()) {
      file_type = 1;
      if (debug_dens) cout << "This file is written in turbomole type!" << endl;
    }
    else if (line.find("gamess-us") < line.size()) {
      file_type = 2;
      if (debug_dens) cout << "This file is written in gamess-us type!" << endl;
    }
    else if (line.find("gaussian") < line.size()) {
      file_type = 3;
      if (debug_dens) cout << "This file is written in gaussian type!" << endl;
    }
    else {
      cout << "This type of file is not supported, please provide another basis set!" << endl;
      return false;
    }
    if (ifile.eof()) {
      cout << "Please provide a basis set in the turbomole, gaussian or gamess-us format compatible with tonto."
        << "Look at the example files \"examble.basis\" and \"examble2.basis\" in the wfn_cpp folder if you want to see how it has to look like" << endl;
      return false;
    }
    while (!(line.find(elements_list[i]) < line.size()) && !ifile.eof()) {
      getline(ifile, line);
      if (debug_dens) cout << "line while search for " << elements_list[i] << " :" << line << endl;
    }
    if (debug_dens && line.find(elements_list[i]) != -1) {
      cout << "I found an entry i know from the element list!" << endl;
      cout << "The line is: " << line << endl;
      Enter();
    }
    if (ifile.eof()) {
      cout << "Could not find the atom you were looking for in the basis set file... " << endl;
      Enter();
      return false;
    }
    unsigned int shell = 0;
    if (line.find("{") == -1) {
      getline(ifile, line);
      if (debug) cout << "I read an additional line!" << endl;
    }
    while (line.find("}") == -1 && !ifile.eof()) {
      getline(ifile, line);
      stringstream stream;
      stream << line;
      int count = 0;
      int nr_exp = 0;
      double pi = 3.14159265358979;
      char c_temp = '?';
      double temp[2];
      int dum = 0;
      if (file_type == 1) {
        stream >> count >> c_temp;
        if (debug_dens) cout << "count: " << count << " type: " << c_temp << endl;
      }
      else if (file_type == 2) {
        stream >> c_temp >> count;
        if (debug_dens) cout << "count: " << count << " type: " << c_temp << endl;
      }
      else if (file_type == 3) {
        stream >> c_temp >> count;
        if (debug_dens) cout << "count: " << count << " type: " << c_temp << endl;
      }
      for (int j = 0; j < count; j++) {
        getline(ifile, line);
        if (debug_dens) {
          cout << "read the " << j << ". line: " << line << endl;
        }
        stringstream stream2;
        stream2 << line;
        if (file_type == 1) stream2 >> temp[0] >> temp[1];
        else if (file_type == 2) stream2 >> dum >> temp[0] >> temp[1];
        else if (file_type == 3) stream2 >> temp[0] >> temp[1];
        //this is where i started copying
        for (int h = 0; h < wave.get_ncen(); h++) {
          //skip atoms taht already have a basis set!
          if (wave.get_atom_basis_set_loaded(h)) continue;
          string temp_label;
          temp_label = wave.get_atom_label(h);
          if (temp_label.find(" ")) temp_label.erase(temp_label.find(" "), 1);
          temp_label.append(":");
          if (elements_list[i].find(temp_label) != -1) {
            if (debug_dens) {
              cout << "It's a match!" << endl;
              cout << "element_label: " << elements_list[i] << " temp_label: " << temp_label << endl;
            }
            switch (c_temp) {
            case 's':
            case 'S':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 1, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type S" << endl;
              break;
            case 'p':
            case 'P':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 2, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type P" << endl;
              break;
            case 'd':
            case 'D':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 3, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type D" << endl;
              break;
            case 'f':
            case 'F':
              if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 4, shell)) {
                cout << "ERROR while pushing back atoms basis set" << endl;
                Enter();
              }
              if (debug_dens) cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                << " and exp: " << temp[0] << " and type F" << endl;
              break;
            default:
              cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
              return false;
            } // end switch of types
          } // end if(find atom_label + : 
        }//end for h = ncen
        nr_exp++;
        if (debug_dens) cout << "recapitulation[" << j << "]... type: " << c_temp << " coef: " << temp[0] << " exp: " << temp[1] << endl;
        if (dum > count) {
          cout << "this should not happen, lets stop before i do something silly!" << endl;
          return false;
        }
      }
      shell++;
    }//end while line != }
    if (debug_dens) cout << "I found }!" << endl;
    shell = 0;
    ifile.seekg(0);
  }//end for element_list.size()
  if (debug_dens) {
    cout << "FINISHED WITH READING BASIS SET!" << endl;
    Enter();
  }
  ifile.close();
  return true;
};

bool delete_basis_set_vanilla(const string& basis_set_path, WFN& wave, bool debug)
{
  for (int a = 0; a < wave.get_ncen(); a++) {
    int nr_prim = wave.get_atom_primitive_count(a);
    for (int p = 0; p < nr_prim; p++) {
      bool succes = wave.erase_atom_primitive(a, 0);
      if (!succes) return false;
    }
  }
  return true;
};

#include "wfn_class.h"
#include "atoms.h"
