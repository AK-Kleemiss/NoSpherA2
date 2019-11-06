#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#ifdef _WIN32
#include <direct.h>
#include <io.h>
#include <windows.h>
#include <Windows.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#include <sys/wait.h>
#endif
#include <fcntl.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <iomanip>
#ifdef __cplusplus__
 #include  <cstdlib>
#else
 #include <stdlib.h>
#endif

#include "convenience.h"

using namespace std;

//--------------------------General convenience terminal functions---------------------------------
inline void QCT(){
 cout << "  ____   _____ _______ " << std::endl;
 cout << " / __ \\ / ____|__   __|" << std::endl;
 cout << "| |  | | |       | |   " << std::endl;
 cout << "| |  | | |       | |   " << std::endl;
 cout << "| |__| | |____   | |   " << std::endl;
 cout << " \\___\\_\\\\_____|  |_|   " << std::endl;
 cout << "                       "<< std::endl;
};

inline bool cuQCT(ofstream& file){
  file << "             ____   _____ _______ " << std::endl;
  file << "            / __ \\ / ____|__   __|" << std::endl;
  file << "  ___ _   _| |  | | |       | |   " << std::endl;
  file << " / __| | | | |  | | |       | |   " << std::endl;
  file << "| (__| |_| | |__| | |____   | |   " << std::endl;
  file << " \\___|\\__,_|\\___\\_\\\\_____|  |_|   " << std::endl;
  file << "                       "<< std::endl;
  return true;
};

inline void cuQCT(){
  cout << "             ____   _____ _______ " << std::endl;
  cout << "            / __ \\ / ____|__   __|" << std::endl;
  cout << "  ___ _   _| |  | | |       | |   " << std::endl;
  cout << " / __| | | | |  | | |       | |   " << std::endl;
  cout << "| (__| |_| | |__| | |____   | |   " << std::endl;
  cout << " \\___|\\__,_|\\___\\_\\\\_____|  |_|   " << std::endl;
  cout << "                       "<< std::endl;
};

inline void copyright(){
	std::cout << "This software is part of the cuQCT software suite developed by Florian Kleemiss.\nPlease give credit and cite corresponding pieces!\n";
};

inline bool copyright(ofstream& file){
  file << "This software is part of the cuQCT software suite developed by Florian Kleemiss.\nPlease give credit and cite corresponding pieces!\n";
  return true;
};

bool yesno(){
	bool end=false;
	while (!end) {
		char dum ='?';
		cout << "(Y/N)?";
		cin >> dum;
		if(dum == 'y'||dum == 'Y'){
			cout << "Okay..." << endl;
			 return true;
		}
		else if(dum == 'N'||dum == 'n') return false;
		else cout << "Sorry, i did not understand that!" << endl;
	}
	return false;
};

void Enter(){
	//cout << "press ENTER to continue... " << flush;
	//cin.ignore();
	//cin.get();
};

bool is_similar(double first, double second, double tolerance) {
	if (first < second * (1 - tolerance) || first > second * (1 + tolerance)) 
		return false;
	else 
		return true;
};

void cls(){
//    cout << string( 100, '\n' );
#ifdef _WIN32
	if(system("CLS")) cout << "this should not happen...!" << endl;
#else
	if(system("clear")) cout << "this should not happen...!" << endl;
#endif
};

string atnr2letter(const int nr){
	vector <string> Labels{"DM","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca"
		,"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"
		,"Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe"
		,"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn"};
	if(nr>86 || nr <=0 ){
		cout << "Only yet implemented from H-Rn, ask Florian for improvements or give a reasonable number between 1-86!" << endl;
		Enter();
		return ("PROBLEM");
	}
	else return Labels[nr];
};

int get_Z_from_label(const char * tmp){
	if (strcmp(tmp, "H" ) == 0)  return 0;
else 	if (strcmp(tmp, "He") == 0)  return 1;
else 	if (strcmp(tmp, "Li") == 0)  return 2;
else 	if (strcmp(tmp, "Be") == 0)  return 3;
else 	if (strcmp(tmp, "B" ) == 0)  return 4;
else 	if (strcmp(tmp, "C" ) == 0)  return 5;
else 	if (strcmp(tmp, "N" ) == 0)  return 6;
else 	if (strcmp(tmp, "O" ) == 0)  return 7;
else 	if (strcmp(tmp, "F" ) == 0)  return 8;
else 	if (strcmp(tmp, "Ne") == 0)  return 9;
else 	if (strcmp(tmp, "Na") == 0)  return 10;
else 	if (strcmp(tmp, "Mg") == 0)  return 11;
else 	if (strcmp(tmp, "Al") == 0)  return 12;
else 	if (strcmp(tmp, "Si") == 0)  return 13;
else 	if (strcmp(tmp, "P" ) == 0)  return 14;
else 	if (strcmp(tmp, "S" ) == 0)  return 15;
else 	if (strcmp(tmp, "Cl") == 0)  return 16;
else 	if (strcmp(tmp, "Ar") == 0)  return 17;
else 	if (strcmp(tmp, "K" ) == 0)  return 18;
else 	if (strcmp(tmp, "Ca") == 0)  return 19;
else 	if (strcmp(tmp, "Sc") == 0)  return 20;
else 	if (strcmp(tmp, "Ti") == 0)  return 21;
else 	if (strcmp(tmp, "V" ) == 0)  return 22;
else 	if (strcmp(tmp, "Cr") == 0)  return 23;
else 	if (strcmp(tmp, "Mn") == 0)  return 24;
else 	if (strcmp(tmp, "Fe") == 0)  return 25;
else 	if (strcmp(tmp, "Co") == 0)  return 26;
else 	if (strcmp(tmp, "Ni") == 0)  return 27;
else 	if (strcmp(tmp, "Cu") == 0)  return 28;
else 	if (strcmp(tmp, "Zn") == 0)  return 29;
else 	if (strcmp(tmp, "Ga") == 0)  return 30;
else 	if (strcmp(tmp, "Ge") == 0)  return 31;
else 	if (strcmp(tmp, "As") == 0)  return 32;
else 	if (strcmp(tmp, "Se") == 0)  return 33;
else 	if (strcmp(tmp, "Br") == 0)  return 34;
else 	if (strcmp(tmp, "Kr") == 0)  return 35;
else 	if (strcmp(tmp, "Rb") == 0)  return 36;
else 	if (strcmp(tmp, "Sr") == 0)  return 37;
else 	if (strcmp(tmp, "Y" ) == 0)  return 38;
else 	if (strcmp(tmp, "Zr") == 0)  return 39;
else 	if (strcmp(tmp, "Nb") == 0)  return 40;
else 	if (strcmp(tmp, "Mo") == 0)  return 41;
else 	if (strcmp(tmp, "Tc") == 0)  return 42;
else 	if (strcmp(tmp, "Ru") == 0)  return 43;
else 	if (strcmp(tmp, "Rh") == 0)  return 44;
else 	if (strcmp(tmp, "Pd") == 0)  return 45;
else 	if (strcmp(tmp, "Ag") == 0)  return 46;
else 	if (strcmp(tmp, "Cd") == 0)  return 47;
else 	if (strcmp(tmp, "In") == 0)  return 48;
else 	if (strcmp(tmp, "Sn") == 0)  return 49;
else 	if (strcmp(tmp, "Sb") == 0)  return 50;
else 	if (strcmp(tmp, "Te") == 0)  return 51;
else 	if (strcmp(tmp, "I" ) == 0)  return 52;
else 	if (strcmp(tmp, "Xe") == 0)  return 53;
else 	if (strcmp(tmp, "Cs") == 0)  return 54;
else 	if (strcmp(tmp, "Ba") == 0)  return 55;
else 	if (strcmp(tmp, "La") == 0)  return 56;
else 	if (strcmp(tmp, "Ce") == 0)  return 57;
else 	if (strcmp(tmp, "Pr") == 0)  return 58;
else 	if (strcmp(tmp, "Nd") == 0)  return 59;
else 	if (strcmp(tmp, "Pm") == 0)  return 60;
else 	if (strcmp(tmp, "Sm") == 0)  return 61;
else 	if (strcmp(tmp, "Eu") == 0)  return 62;
else 	if (strcmp(tmp, "Gd") == 0)  return 63;
else 	if (strcmp(tmp, "Tb") == 0)  return 64;
else 	if (strcmp(tmp, "Dy") == 0)  return 65;
else 	if (strcmp(tmp, "Ho") == 0)  return 66;
else 	if (strcmp(tmp, "Er") == 0)  return 67;
else 	if (strcmp(tmp, "Tm") == 0)  return 68;
else 	if (strcmp(tmp, "Yb") == 0)  return 69;
else 	if (strcmp(tmp, "Lu") == 0)  return 70;
else 	if (strcmp(tmp, "Hf") == 0)  return 71;
else 	if (strcmp(tmp, "Ta") == 0)  return 72;
else 	if (strcmp(tmp, "W" ) == 0)  return 73;
else 	if (strcmp(tmp, "Re") == 0)  return 74;
else 	if (strcmp(tmp, "Os") == 0)  return 75;
else 	if (strcmp(tmp, "Ir") == 0)  return 76;
else 	if (strcmp(tmp, "Pt") == 0)  return 77;
else 	if (strcmp(tmp, "Au") == 0)  return 78;
else 	if (strcmp(tmp, "Hg") == 0)  return 79;
else 	if (strcmp(tmp, "Ti") == 0)  return 80;
else 	if (strcmp(tmp, "Pb") == 0)  return 81;
else 	if (strcmp(tmp, "Bi") == 0)  return 82;
else 	if (strcmp(tmp, "Po") == 0)  return 83;
else 	if (strcmp(tmp, "At") == 0)  return 84;
else 	if (strcmp(tmp, "Rn") == 0)  return 85;
else                                 return -1;
}

void copy_file(string from, string to){
    ifstream source(from.c_str(), ios::binary);
    ofstream dest(to.c_str(), ios::binary);

    dest << source.rdbuf();

    source.close();
    dest.close();
};

//----------------------------WFN Handling-------------------------------------------------------------

bool writewfn(WFN &wavefunction, const string &path, bool &debug,bool occ){
	try{
		if(wavefunction.write_wfn (path,debug, occ)){
			wavefunction.set_path(path);
			return true;
		}
		else{
			cout << "Sorry, something went wrong, try to look above what i might be..\n";
			throw (int)1;
		}
	}
	catch(exception& e)
	{
		cerr << "Sorry, something went wrong, try to look above what i might be.. the error code is " << e.what() << endl;
		throw (int) 3;
	}
}; 

bool readwfn(WFN &wavefunction, const string &path, bool &debug){
	try{
		if(wavefunction.read_wfn (path,debug)){
			wavefunction.set_path(path);
			return true;
		}
		else{
			cout << "Sorry, something went wrong, try to look above what i might be..\n";
			throw (int)1;
		}
	}
	catch(exception& e)
	{
		cerr << "Sorry, something went wrong, try to look above what i might be.. the error code is " << e.what() << endl;
		throw (int) 3;
	}
};

bool readwfn(WFN& wavefunction, const string& path, bool& debug, ofstream &file) {
	try {
		if (wavefunction.read_wfn(path, debug,file)) {
			wavefunction.set_path(path);
			return true;
		}
		else {
			cout << "Sorry, something went wrong, try to look above what i might be..\n";
			throw (int)1;
		}
	}
	catch (exception& e)
	{
		cerr << "Sorry, something went wrong, try to look above what i might be.. the error code is " << e.what() << endl;
		throw (int)3;
	}
};

//---------------------------Configuration files ---------------------------------------------------

string get_home_path(void) {
#ifdef _WIN32
	string temp1 = getenv("HOMEDRIVE");
	string temp2 = getenv("HOMEPATH");
	temp1.append(temp2);
	return temp1;
#else
	string home = getenv("HOME");
	return home;
#endif
}

void join_path(string &s1, string &s2) {
	char ch = s1.back();
#ifdef _WIN32
	s1.append( "\\" );
#else
	if(s1.substr(s1.length() -1) == "/")
		s1.append("/");
#endif
	s1.append(s2);
}

string get_filename_from_path(const string &input){
#ifdef _WIN32
	return input.substr(input.rfind("\\")+1);
#else
	return input.substr(input.rfind("/")+1);
#endif
}

string get_foldername_from_path(const string &input){
#ifdef _WIN32
	return input.substr(0,input.rfind("\\")+1);
#else
	return input.substr(0,input.rfind("/")+1);
#endif
}

string get_basename_without_ending(const string &input){
	return input.substr(0,input.rfind("."));
}

void write_template_confi(){
	ofstream conf;
	string line;
	string programs = get_home_path();
	string filename = ".cuQCT.conf";
	join_path(programs, filename);
	if(exists(programs)){
		cout << "File already exist, do you want to overwrite it?" << endl;
		if(!yesno()) return;
	}
	conf.open(programs.c_str(),ios::out);
#ifdef _WIN32
	conf << "gaussian=\"D:\\g09\\g09\\\"" << endl;
	conf << "turbomole=\"D:\\turbomole\\dscf7.1\\\"" << endl;
	conf << "basis=\"D:\\tonto\\basis_sets\\\"" << endl;
#else
	conf << "gaussian=\"/usr/local/g09/g09\"" << endl;
	conf << "turbomole=\"/usr/local/bin/dscf7.1\"" << endl;
	conf << "basis=\"/basis_sets/\"" << endl;
#endif
	conf << "cpu=4" << endl;
	conf << "mem=4.0" << endl;
	conf << "rho=1" << endl;
	conf << "rdg=1" << endl;
	conf << "eli=0" << endl;
	conf << "elf=0" << endl;
	conf << "lap=0" << endl;
	conf << "esp=0" << endl;
	conf << "efv=0" << endl;
	conf << "def=0" << endl;
	conf << "hir=0" << endl;
	conf.flush();
	conf.close();
#ifdef _WIN32
//	const wchar_t* fileLPCWSTR = programs.c_str();
	wstring stemp = wstring(programs.begin(), programs.end());
	int attr = GetFileAttributes(stemp.c_str());
	if ((attr & FILE_ATTRIBUTE_HIDDEN) == 0) {
		SetFileAttributes(stemp.c_str(), attr | FILE_ATTRIBUTE_HIDDEN);
	}
#endif
	return;
};

int program_confi(string &gaussian_path, string &turbomole_path, string &basis, int &ncpus, float &mem, bool debug, bool expert, unsigned int counter) {
	counter++;
	if (counter == 3) {
		cout << "Too many iterations of tries to read config file, better abort..." << endl;
		return -1;
	}
	string programs = get_home_path();
	string filename = ".cuQCT.conf";
	join_path(programs, filename);
	ifstream conf(programs.c_str());
	if(debug) cout << programs << endl;
	string line;
	if(conf.good()){
		if(debug) cout << "File is valid, continuing..." << endl;
	}
	else{
		if (expert) {
			cout << "couldn't open or find .cuQCT.conf, in your home folder: " << programs << ", do you want me to write a template for you?" << endl;
			if (yesno()) {
				write_template_confi();
				if(program_confi(gaussian_path, turbomole_path, basis, ncpus, mem, debug, expert, counter)!=1) return -1;
				cout << "Wrote a template for you, read default values!" << endl;
				return 0;
			}
			else return -1;
		}
		else {
			write_template_confi();
			program_confi(gaussian_path, turbomole_path, basis, ncpus, mem, debug);
			return 0;
		}
	}
	conf.seekg(0);
	getline(conf,line);
	size_t length;
	char tempchar[200];
	int run=0;
	while(!conf.eof()){
		switch(run){
			case 0:
				length = line.copy(tempchar,line.size()-11,10);
				tempchar[length]='\0';
				gaussian_path=tempchar;
				break;
			case 1:
				length = line.copy(tempchar,line.size()-12,11);
				tempchar[length]='\0';
				turbomole_path=tempchar;
				break;
			case 2:
				length = line.copy(tempchar,line.size()-8,7);
				tempchar[length]='\0';
				basis=tempchar;
				break;
			case 3: length = line.copy(tempchar,line.size()-3,4);
				tempchar[length]='\0';
				sscanf(tempchar, "%d", &ncpus);
				break;
			case 4: length = line.copy(tempchar,line.size()-3,4);
				tempchar[length]='\0';
				sscanf(tempchar, "%f", &mem);
				break;
			default:
				if(debug) cout << "found everything i was looking for, if you miss something check the switch" << endl;
				break;
		}
		if(debug) cout << run << ". line: " << tempchar << endl;
		run++;
		getline(conf,line);
	}
	return 1;
};

bool check_bohr(WFN &wave, bool interactive, bool debug){
	double min_length=300.0;
	for(int i=0; i<wave.get_ncen(); i++){
		double atom1[3];
		for(int x=0; x<3; x++)
			atom1[x]=wave.get_atom_coordinate (i, x, debug);
		for(int j=i+1; j<wave.get_ncen(); j++){
			double atom2[3];
			for(int x=0; x<3; x++)
				atom2[x]=wave.get_atom_coordinate (j,x,debug);
			double d[3];
			d[0]=atom1[0]-atom2[0];
			d[1]=atom1[1]-atom2[1];
			d[2]=atom1[2]-atom2[2];
			double length=sqrt(pow(d[0],2)+pow(d[1],2)+pow(d[2],2));
			if(debug) cout << "Length for: " << i << ";" << j << ": " << length << ", min_length: " << min_length << endl;
			if(length<min_length) min_length=length;
		}
	}
	if(debug || interactive){
		if(min_length<2)
			cout << "Decided it's written in Angstrom" << endl;
		else
			cout << "Decided it's written in Bohr" << endl;
	}
	if(interactive){
	cout << "Do you agree ";
	if(yesno()) return (!(min_length<2));
	else{
		wave.set_dist_switch();
		return (min_length<2);
	}}
	else return (!(min_length<2));
};

int filetype_identifier(string &file, bool debug){
	/*
	List of filetypes and correpsonding values:
					-1: unreadable keyword
	-i/o *.wfn 		2: wfn
	-i/o *.ffn 		4: ffn
	-i *.out 		1: crystal output
	-c/o *.cub(e) 	3: cube file
	-g/o *.grd      6: XDGraph grid file
	-o *.(F)fc(C)hk 5: fchk
	*/
	if (debug){
		cout << "Testing WFN:  " << file.find(".wfn") << endl
			 << "Testing out:  " << file.find(".out") << endl
			 << "Testing FFN:  " << file.find(".ffn") << endl
			 << "Testing CUB:  " << file.find(".cub") << endl
			 << "Testing CUBE: " << file.find(".cube") << endl
			 << "Testing Grid: " << file.find(".grd") << endl
			 << "Testing fchk: " << file.find(".fchk") << endl
			 << "Testing FChk: " << file.find(".FChk") << endl
			 << "Testing Fchk: " << file.find(".Fchk") << endl;
		cout << "string::npos: " << string::npos << endl;
	}
	int temp_type=0;
	size_t found,temp;
	temp=0;
	if(debug) cout << "Temp before any checks: " << temp << endl;
	vector <string> types{".out",".wfn",".ffn",".cub",".cube",".grd",".fchk",".Fchk",".FChk"};
	if(file.find(".wfn")!=string::npos){
		if (debug) cout << "Checking for" << ".wfn" << endl;
		temp_type = 2;
		found = file.rfind(".wfn");
		if(debug) cout << "Found: " << found << endl;
		for(int i=0; i<types.size();i++)
			if(file.rfind(types[i])>=found && file.rfind(types[i])!=string::npos)
				temp=file.rfind(types[i]);
		if(debug) cout << "Temp: " << temp << endl;
		if(temp==found) return temp_type;
		else {
			temp=0;
			if (debug) cout << "Moving on!" << endl;
		}
	}
	if(file.find(".out")!=string::npos){
		if (debug) cout << "Checking for" << ".out" << endl;
		temp_type = 1;
		found = file.rfind(".out");
		for(int i=0; i<types.size();i++)
			if(file.rfind(types[i])>=found && file.rfind(types[i])!=string::npos)
				temp=file.rfind(types[i]);
		if(temp==found) return temp_type;
		else temp=0;
	}
	if(file.find(".ffn")!=string::npos){
		if (debug) cout << "Checking for" << ".ffn" << endl;
		temp_type = 4;
		found = file.rfind(".ffn");
		for(int i=0; i<types.size();i++)
			if(file.rfind(types[i])>=found && file.rfind(types[i])!=string::npos)
				temp=file.rfind(types[i]);
		if(temp==found) return temp_type;
		else temp=0;
	}
	if(file.find(".cub")!=string::npos){
		if (debug) cout << "Checking for" << ".cub" << endl;
		temp_type = 3;
		found = file.rfind(".cub");
		for(int i=0; i<types.size();i++)
			if(file.rfind(types[i])>=found && file.rfind(types[i])!=string::npos)
				temp=file.rfind(types[i]);
		if(temp==found) return temp_type;
		else{
			temp=0;
			if (debug) cout << "Moving on!" << endl;
		}
	}
	if(file.find(".cube")!=string::npos){
		temp_type = 3;
		found = file.rfind(".cube");
		for(int i=0; i<types.size();i++)
			if(file.rfind(types[i])>=found && file.rfind(types[i])!=string::npos)
				temp=file.rfind(types[i]);
		if(temp==found) return temp_type;
		else temp=0;
	}
	if(file.find(".grd")!=string::npos){
		temp_type = 6;
		found = file.rfind(".grd");
		for(int i=0; i<types.size();i++)
			if(file.rfind(types[i])>=found && file.rfind(types[i])!=string::npos)
				temp=file.rfind(types[i]);
		if(temp==found) return temp_type;
		else temp=0;
	}
	if(file.find(".fchk")!=string::npos){
		temp_type = 5;
		found = file.rfind(".fchk");
		for(int i=0; i<types.size();i++)
			if(file.rfind(types[i])>=found && file.rfind(types[i])!=string::npos)
				temp=file.rfind(types[i]);
		if(temp==found) return temp_type;
		else temp=0;
	}
	if(file.find(".FChk")!=string::npos){
		temp_type = 5;
		found = file.rfind(".FChk");
		for(int i=0; i<types.size();i++)
			if(file.rfind(types[i])>=found && file.rfind(types[i])!=string::npos)
				temp=file.rfind(types[i]);
		if(temp==found) return temp_type;
		else temp=0;
	}
	if(file.find(".Fchk")!=string::npos){
		temp_type = 5;
		found = file.rfind(".Fchk");
		for(int i=0; i<types.size();i++)
			if(file.rfind(types[i])>=found && file.rfind(types[i])!=string::npos)
				temp=file.rfind(types[i]);
		if(temp==found) return temp_type;
		else temp=0;
	}
	return -1;
}

string shrink_string(string &input){
	while(input.find(" ")!=-1){input.erase(input.find(" "),1);}
	while(input.find("1")!=-1){input.erase(input.find("1"),1);}
	while(input.find("2")!=-1){input.erase(input.find("2"),1);}
	while(input.find("3")!=-1){input.erase(input.find("3"),1);}
	while(input.find("4")!=-1){input.erase(input.find("4"),1);}
	while(input.find("5")!=-1){input.erase(input.find("5"),1);}
	while(input.find("6")!=-1){input.erase(input.find("6"),1);}
	while(input.find("7")!=-1){input.erase(input.find("7"),1);}
	while(input.find("8")!=-1){input.erase(input.find("8"),1);}
	while(input.find("9")!=-1){input.erase(input.find("9"),1);}
	while(input.find("0")!=-1){input.erase(input.find("0"),1);}
	while(input.find("(")!=-1){input.erase(input.find("("),1);}
	while(input.find(")")!=-1){input.erase(input.find(")"),1);}
	return input;
};

string shrink_string_to_atom(string &input, const int atom_number){
	while(input.find(" ")!=-1){input.erase(input.find(" "),1);}
	while(input.find("1")!=-1){input.erase(input.find("1"),1);}
	while(input.find("2")!=-1){input.erase(input.find("2"),1);}
	while(input.find("3")!=-1){input.erase(input.find("3"),1);}
	while(input.find("4")!=-1){input.erase(input.find("4"),1);}
	while(input.find("5")!=-1){input.erase(input.find("5"),1);}
	while(input.find("6")!=-1){input.erase(input.find("6"),1);}
	while(input.find("7")!=-1){input.erase(input.find("7"),1);}
	while(input.find("8")!=-1){input.erase(input.find("8"),1);}
	while(input.find("9")!=-1){input.erase(input.find("9"),1);}
	while(input.find("0")!=-1){input.erase(input.find("0"),1);}
	while(input.find("(")!=-1){input.erase(input.find("("),1);}
	while(input.find(")")!=-1){input.erase(input.find(")"),1);}
	string temp=atnr2letter(atom_number);
	if(input.find(temp)!=1) return temp;
	if (temp!="PROBLEM")
	 while(input.size()>temp.size())
	  input.pop_back();
	return input;
};

bool open_file_dialog(string &path, bool debug, vector <string> filter){
	char pwd[1024];
	if(GetCurrentDir( pwd, 1024)==NULL) return false;
	string current_path(pwd);
#ifdef _WIN32
	char filename[ 1024 ];

	OPENFILENAMEA ofn;
	    ZeroMemory( &filename, sizeof( filename ) );
	    ZeroMemory( &ofn,      sizeof( ofn ) );
	    ofn.lStructSize  = sizeof( ofn );
	    ofn.hwndOwner    = NULL;  // If you have a window to center over, put its HANDLE here
	    ofn.lpstrFilter  = "wfn Files\0*.wfn\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0Any File\0*.*\0";
	    ofn.lpstrFile    = filename;
	    ofn.nMaxFile     = 1024;
	    ofn.lpstrTitle   = "Select a File";
	    ofn.Flags        = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST;

	if (GetOpenFileNameA( &ofn )){
		if(debug) cout << "You chose the file \"" << filename << "\"\n";
		if(exists(filename)){
			path=filename;
			return true;
		}
	}
	else
	{
		// All this stuff below is to tell you exactly how you messed up above.
		// Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
		switch (CommDlgExtendedError())
		{
			case CDERR_DIALOGFAILURE   : cout << "CDERR_DIALOGFAILURE\n";   break;
			case CDERR_FINDRESFAILURE  : cout << "CDERR_FINDRESFAILURE\n";  break;
			case CDERR_INITIALIZATION  : cout << "CDERR_INITIALIZATION\n";  break;
			case CDERR_LOADRESFAILURE  : cout << "CDERR_LOADRESFAILURE\n";  break;
			case CDERR_LOADSTRFAILURE  : cout << "CDERR_LOADSTRFAILURE\n";  break;
			case CDERR_LOCKRESFAILURE  : cout << "CDERR_LOCKRESFAILURE\n";  break;
			case CDERR_MEMALLOCFAILURE : cout << "CDERR_MEMALLOCFAILURE\n"; break;
			case CDERR_MEMLOCKFAILURE  : cout << "CDERR_MEMLOCKFAILURE\n";  break;
			case CDERR_NOHINSTANCE     : cout << "CDERR_NOHINSTANCE\n";     break;
			case CDERR_NOHOOK          : cout << "CDERR_NOHOOK\n";          break;
			case CDERR_NOTEMPLATE      : cout << "CDERR_NOTEMPLATE\n";      break;
			case CDERR_STRUCTSIZE      : cout << "CDERR_STRUCTSIZE\n";      break;
			case FNERR_BUFFERTOOSMALL  : cout << "FNERR_BUFFERTOOSMALL\n";  break;
			case FNERR_INVALIDFILENAME : cout << "FNERR_INVALIDFILENAME\n"; break;
			case FNERR_SUBCLASSFAILURE : cout << "FNERR_SUBCLASSFAILURE\n"; break;
			default                    : cout << "You cancelled.\n";
		}
	}
	return false;
#else
	char file[1024];
	string command;
	command="zenity --file-selection --title=\"Select a file to load\" --filename=\"";
	command += current_path;
	command += "/\"";
	for(int i=0; i<filter.size(); i++){
		command += " --file-filter=\"";
		command += filter[i];
		command += "\" ";
	}
	command += " 2> /dev/null";
	FILE *f = popen(command.c_str(), "r");
	if(!f){
		cout << "ERROR" << endl;
		return false;
	}
	if(fgets(file, 1024, f)==NULL) return false;
	if (debug) cout << "Filename: " << file << endl;
	path=file;
	stringstream ss(path);
	getline(ss, path);
	if(pclose(f)!=0) cout << "Zenity returned non zero, whatever that means..." << endl;
#endif
	return true;
};

bool save_file_dialog(string &path, bool debug, const vector <string> &endings, const string &filename_given){
	char pwd[1024];
	if(GetCurrentDir( pwd, 1024)==NULL) return false;
	string current_path(pwd);
#ifdef _WIN32
	char filename[ 1024 ];

	OPENFILENAMEA sfn;
	 ZeroMemory( &filename, sizeof( filename ) );
	 ZeroMemory( &sfn,      sizeof( sfn ) );
	 sfn.lStructSize  = sizeof( sfn );
	 sfn.hwndOwner    = NULL;  // If you have a window to center over, put its HANDLE here
	 sfn.lpstrFilter  = "wfn Files\0*.wfn\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0Any File\0*.*\0";
	 sfn.lpstrFile    = filename;
	 sfn.nMaxFile     = 1024;
	 sfn.lpstrTitle   = "Select a File, yo!";
	 sfn.Flags        = OFN_DONTADDTORECENT;
	bool end=false;
	while(!end){
		if (GetSaveFileNameA( &sfn )){
			if(debug) cout << "You chose the file \"" << filename << "\"\n";
			if(exists(filename)){
				cout << filename << " exists, do you want to overwrite it?";
				if(yesno()){
					path=filename;
					bool found=false;
					for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
					if(found) end=true;
				}
				else return false;
			}
			else{
				path=filename;
				bool found=false;
				for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
				if(found) end=true;
			}
		}
		else
		{
			// All this stuff below is to tell you exactly how you messed up above.
			// Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
			switch (CommDlgExtendedError())
			{
				case CDERR_DIALOGFAILURE   : cout << "CDERR_DIALOGFAILURE\n";   break;
				case CDERR_FINDRESFAILURE  : cout << "CDERR_FINDRESFAILURE\n";  break;
				case CDERR_INITIALIZATION  : cout << "CDERR_INITIALIZATION\n";  break;
				case CDERR_LOADRESFAILURE  : cout << "CDERR_LOADRESFAILURE\n";  break;
				case CDERR_LOADSTRFAILURE  : cout << "CDERR_LOADSTRFAILURE\n";  break;
				case CDERR_LOCKRESFAILURE  : cout << "CDERR_LOCKRESFAILURE\n";  break;
				case CDERR_MEMALLOCFAILURE : cout << "CDERR_MEMALLOCFAILURE\n"; break;
				case CDERR_MEMLOCKFAILURE  : cout << "CDERR_MEMLOCKFAILURE\n";  break;
				case CDERR_NOHINSTANCE     : cout << "CDERR_NOHINSTANCE\n";     break;
				case CDERR_NOHOOK          : cout << "CDERR_NOHOOK\n";          break;
				case CDERR_NOTEMPLATE      : cout << "CDERR_NOTEMPLATE\n";      break;
				case CDERR_STRUCTSIZE      : cout << "CDERR_STRUCTSIZE\n";      break;
				case FNERR_BUFFERTOOSMALL  : cout << "FNERR_BUFFERTOOSMALL\n";  break;
				case FNERR_INVALIDFILENAME : cout << "FNERR_INVALIDFILENAME\n"; break;
				case FNERR_SUBCLASSFAILURE : cout << "FNERR_SUBCLASSFAILURE\n"; break;
				default                    : cout << "You cancelled.\n";
			}
			return false;
		}
	}
	return true;
#else
	char file[1024];
	string command;
	command="zenity --file-selection --title=\"Select where to save\" --filename=\"";
	command += current_path;
	command += filename_given;
	command += "/\" --save --confirm-overwrite 2> /dev/null";
	bool end=false;
	while(!end){
		FILE *f = popen(command.c_str(), "r");
		if(!f){
			cout << "ERROR" << endl;
			return false;
		}
		if(fgets(file, 1024, f)==NULL) return false;
		if (debug) cout << "Filename: " << file << endl;
		path=file;
		stringstream ss(path);
		getline(ss, path);
		if(debug) cout << "Path: " << path << endl;
		if(pclose(f)!=0) cout << "Zenity returned non zero, whatever that means..." << endl;
		bool found=false;
		for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
		if(found) end=true;
	}
#endif
	return true;
};

bool save_file_dialog(string &path, bool debug, const vector <string> &endings){
	char pwd[1024];
	if(GetCurrentDir( pwd, 1024)==NULL) return false;
	string current_path(pwd);
#ifdef _WIN32
	char filename[ 1024 ];

	OPENFILENAMEA sfn;
	 ZeroMemory( &filename, sizeof( filename ) );
	 ZeroMemory( &sfn,      sizeof( sfn ) );
	 sfn.lStructSize  = sizeof( sfn );
	 sfn.hwndOwner    = NULL;  // If you have a window to center over, put its HANDLE here
	 sfn.lpstrFilter  = "wfn Files\0*.wfn\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0Any File\0*.*\0";
	 sfn.lpstrFile    = filename;
	 sfn.nMaxFile     = 1024;
	 sfn.lpstrTitle   = "Select a File, yo!";
	 sfn.Flags        = OFN_DONTADDTORECENT;
	bool end=false;
	while(!end){
		if (GetSaveFileNameA( &sfn )){
			if(debug) cout << "You chose the file \"" << filename << "\"\n";
			if(exists(filename)){
				cout << filename << " exists, do you want to overwrite it?";
				if(yesno()){
					path=filename;
					bool found=false;
					for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
					if(found) end=true;
				}
				else return false;
			}
			else{
				path=filename;
				bool found=false;
				for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
				if(found) end=true;
			}
		}
		else
		{
			// All this stuff below is to tell you exactly how you messed up above.
			// Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
			switch (CommDlgExtendedError())
			{
				case CDERR_DIALOGFAILURE   : cout << "CDERR_DIALOGFAILURE\n";   break;
				case CDERR_FINDRESFAILURE  : cout << "CDERR_FINDRESFAILURE\n";  break;
				case CDERR_INITIALIZATION  : cout << "CDERR_INITIALIZATION\n";  break;
				case CDERR_LOADRESFAILURE  : cout << "CDERR_LOADRESFAILURE\n";  break;
				case CDERR_LOADSTRFAILURE  : cout << "CDERR_LOADSTRFAILURE\n";  break;
				case CDERR_LOCKRESFAILURE  : cout << "CDERR_LOCKRESFAILURE\n";  break;
				case CDERR_MEMALLOCFAILURE : cout << "CDERR_MEMALLOCFAILURE\n"; break;
				case CDERR_MEMLOCKFAILURE  : cout << "CDERR_MEMLOCKFAILURE\n";  break;
				case CDERR_NOHINSTANCE     : cout << "CDERR_NOHINSTANCE\n";     break;
				case CDERR_NOHOOK          : cout << "CDERR_NOHOOK\n";          break;
				case CDERR_NOTEMPLATE      : cout << "CDERR_NOTEMPLATE\n";      break;
				case CDERR_STRUCTSIZE      : cout << "CDERR_STRUCTSIZE\n";      break;
				case FNERR_BUFFERTOOSMALL  : cout << "FNERR_BUFFERTOOSMALL\n";  break;
				case FNERR_INVALIDFILENAME : cout << "FNERR_INVALIDFILENAME\n"; break;
				case FNERR_SUBCLASSFAILURE : cout << "FNERR_SUBCLASSFAILURE\n"; break;
				default                    : cout << "You cancelled.\n";
			}
			return false;
		}
	}
	return true;
#else
	char file[1024];
	string command;
	command="zenity --file-selection --title=\"Select where to save\" --filename=\"";
	command += current_path;
	command += "/\" --save --confirm-overwrite 2> /dev/null";
	bool end=false;
	while(!end){
		FILE *f = popen(command.c_str(), "r");
		if(!f){
			cout << "ERROR" << endl;
			return false;
		}
		if(fgets(file, 1024, f)==NULL) return false;
		if (debug) cout << "Filename: " << file << endl;
		path=file;
		stringstream ss(path);
		getline(ss, path);
		if(debug) cout << "Path: " << path << endl;
		if(pclose(f)!=0) cout << "Zenity returned non zero, whatever that means..." << endl;
		bool found=false;
		for(int i=0; i< endings.size(); i++) if (path.find(endings[i])!=string::npos) found=true;
		if(found) end=true;
	}
#endif
	return true;
};

void select_cubes(vector <vector <unsigned int> > &selection, vector<WFN> &wavy, unsigned int nr_of_cubes, bool wfnonly, bool debug){
	//asks which wfn to use, if wfnonly is set or whcih cubes up to nr of cubes to use
	//Returns values in selection[0][i] for iths selection of wavefunction and
	// selection[1][i] for iths selection of cube
	cout << "Which of the following cubes to use? Need to select " << nr_of_cubes << " file";
	if(nr_of_cubes>1) cout << "s in total."<< endl;
	else cout << "." << endl;
	cout << endl << endl;
	for(int w=0; w<wavy.size(); w++){
		stringstream stream;
		cout << "_____________________________________________________________" << endl;
		cout << "WFN ";
		stream << setw(2) << w;
		cout << stream.str() << ") " << get_basename_without_ending(wavy[w].get_path()) << endl;
		stream.str("");
		for(int c=0; c< wavy[w].cub.size(); c++){
			if(c==0) cout << "        |" << endl << "Cube    |" << endl ;
			else cout << "        |" << endl;
			if(!wfnonly){
				cout << setw(2) << w;
				cout << ".";
				cout << setw(2) << c;
			}
			else cout << "     ";
			cout << "   |_ " << get_basename_without_ending(wavy[w].cub[c].path);
			if(!exists(wavy[w].cub[c].path))
				cout << " (MEM ONLY)";
			cout << endl;
		}
		cout << "_____________________________________________________________" << endl << endl << endl;
	}
	bool happy=false;
	unsigned int selected_cubes=0;
	do{
		cout << "Select " << selected_cubes+1 << ". ";
		if(wfnonly) cout << "WFN ";
		else cout << "cube ";
		cout << "please: ";
		string input;
		cin >> input;
		if(!wfnonly){
			if(input.find('.')==string::npos){
				cout << "no . found in input!" << endl;
				continue;
			}
		}
		else{
			if(input.find('.')==string::npos) cout << "Ignoring the .!" << endl;
			unsigned int nr_wave=fromString<unsigned int> (input);
			if(nr_wave<0 || nr_wave>=wavy.size()){
				cout << "Invalid choice!" << endl;
				continue;
			}
			selected_cubes++;
			selection[0].push_back(nr_wave);
			if(selected_cubes==nr_of_cubes) return;
			else continue;
		}
		if(debug){
			cout << "input: " << input << endl;
			cout << "with . found at: " << input.find('.') << endl;
			cout << "substr1: " << input.substr(0,input.find('.')) << endl;
			cout << "substr2: " << input.substr(input.find('.')+1) << endl;
		}
		string wave(input.substr(0,input.find('.')));
		string cube(input.substr(input.find('.')+1));
		unsigned int nr_wave=fromString<unsigned int> (wave);
		unsigned int nr_cube=fromString<unsigned int> (cube);
		if (debug) cout << "Translated: " << nr_wave << " " << nr_cube << endl;
		if(nr_wave<0 || nr_wave>=wavy.size() || nr_cube<0 || nr_cube>=wavy[nr_wave].cub.size()){
			cout << "Invalid choice!" << endl;
			continue;
		}
		selection[0][selected_cubes]=nr_wave;
		selection[1][selected_cubes]=nr_cube;
		selected_cubes++;
		if(selected_cubes==nr_of_cubes){
			if(debug) cout << "Going to return!" << endl;
			return;
		}
	}while(true);
};

bool unsaved_files(vector<WFN> &wavy){
	for(int w=0; w<wavy.size(); w++)
		for(int c=0; c< wavy[w].cub.size(); c++)
			if(!exists(wavy[w].cub[c].path))
				return true;
	return false;
};

