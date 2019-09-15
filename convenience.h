#ifndef __CONVENIENCE_H__
#define __CONVENIENCE_H__

#include <iostream>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <vector>

class WFN;

//------------------general functions for easy use of terminal input--------------------
bool yesno();
bool is_similar(double first, double second, double tolerance);
void Enter();
void cls();
std::string get_home_path(void);
void join_path(std::string &s1, std::string &s2);
template<class T> std::string toString(const T& t)
{
     std::ostringstream stream;
     stream << t;
     return stream.str();
}

template<class T> T fromString(const std::string& s)
{
     std::istringstream stream (s);
     T t;
     stream >> t;
     return t;
}
inline void QCT(){
std::cout << "  ____   _____ _______ " << std::endl;
std::cout << " / __ \\ / ____|__   __|" << std::endl;
std::cout << "| |  | | |       | |   " << std::endl;
std::cout << "| |  | | |       | |   " << std::endl;
std::cout << "| |__| | |____   | |   " << std::endl;
std::cout << " \\___\\_\\\\_____|  |_|   " << std::endl;
std::cout << "                       "<< std::endl;
};
inline void cuQCT(ofstream &file){
  file << "             ____   _____ _______ " << std::endl;
  file << "            / __ \\ / ____|__   __|" << std::endl;
  file << "  ___ _   _| |  | | |       | |   " << std::endl;
  file << " / __| | | | |  | | |       | |   " << std::endl;
  file << "| (__| |_| | |__| | |____   | |   " << std::endl;
  file << " \\___|\\__,_|\\___\\_\\\\_____|  |_|   " << std::endl;
  file << "                       "<< std::endl;
};
inline void cuQCT(){
std::cout << "             ____   _____ _______ " << std::endl;
std::cout << "            / __ \\ / ____|__   __|" << std::endl;
std::cout << "  ___ _   _| |  | | |       | |   " << std::endl;
std::cout << " / __| | | | |  | | |       | |   " << std::endl;
std::cout << "| (__| |_| | |__| | |____   | |   " << std::endl;
std::cout << " \\___|\\__,_|\\___\\_\\\\_____|  |_|   " << std::endl;
std::cout << "                       "<< std::endl;
};
inline void copyright(){
	std::cout << "This software is part of the cuQCT software suite developed by Florian Kleemiss.\nPlease give credit and cite corresponding pieces!\n";
}
inline void copyright(ofstream &file){
  file << "This software is part of the cuQCT software suite developed by Florian Kleemiss.\nPlease give credit and cite corresponding pieces!\n";
}
inline bool exists(const std::string &name){
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
};
std::string atnr2letter(const int nr);
void copy_file(std::string from, std::string to);
std::string shrink_string(std::string &input);
std::string shrink_string_to_atom(std::string &input, const int atom_number);
std::string get_filename_from_path(const std::string &input);
std::string get_foldername_from_path(const std::string &input);
std::string get_basename_without_ending(const std::string &input);
//------------------Functions to handle .wfn files--------------------------------------
bool writewfn(WFN &wavefunction, const std::string &path, bool &debug, bool occ);
bool readwfn(WFN &wavefunction, const std::string &path, bool &debug);
//------------------Functions to work with configuration files--------------------------
void write_template_confi();
int program_confi(std::string &gaussian_path, std::string &turbomole_path, 
                   std::string &basis, int &ncpus, float &mem, bool debug = false, bool expert = false, unsigned int counter = 0);
bool check_bohr(WFN &wave, bool interactive, bool debug);
int filetype_identifier(std::string &file, bool debug = false);

bool open_file_dialog(std::string &path, bool debug, std::vector <std::string> filter);
bool save_file_dialog(std::string &path, bool debug, const std::vector<std::string> &endings, const std::string &filename_given);
bool save_file_dialog(std::string &path, bool debug, const std::vector<std::string> &endings);
void select_cubes(std::vector <std::vector <unsigned int> > &selection, std::vector<WFN> &wavy, unsigned int nr_of_cubes=1, bool wfnonly=false, bool debug = false);
bool unsaved_files(std::vector<WFN> &wavy);
int get_Z_from_label(const char * tmp);

#include "wfn_class.h"

#endif
