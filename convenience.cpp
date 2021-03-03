#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
#ifdef _WIN32
#include <direct.h>
#include <io.h>
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
#include "cell.h"

using namespace std;
const double c_pi = 3.141592653589793;

//--------------------------General convenience terminal functions---------------------------------


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

bool is_similar_rel(double first, double second, double tolerance) {
	double diff = abs(first - second);
	if ( diff > abs((first + second+0.01) * tolerance/2)) 
		return false;
	else 
		return true;
};

bool is_similar(double first, double second, double tolerance) {
	double diff = abs(first - second);
	if (diff > pow(10,tolerance))
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
		,"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn"
		,"Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"
		};
	if(nr>113 || nr <=0 ){
		cout << "Only yet implemented from H-Lr, ask Florian for improvements or give a reasonable number between 1-113!" << endl;
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

cosinus_annaeherung::cosinus_annaeherung() : mSize(0), mBase_values(nullptr), mStepwidth(1.0) {
	resize(100);
}

void cosinus_annaeherung::resize(size_t size)
{
	mSize = size;
	if (mBase_values) delete[] mBase_values;
	mBase_values = new double[mSize + 1];
#pragma omp parallel for
	for (auto i = 0; i < mSize + 1; i++)  // Fuer einen Werte mehr die Stueststellen speichern
	{
		double y = cos((MPI2 * i) / mSize);
		// cout << "resize: i="<<i<<" y=" << y << endl;
		mBase_values[i] = y;
	}
	mStepwidth = MPI2 / size;
}

double cosinus_annaeherung::calculate_error_at(double x) const
{
	return cos(x) - get(x);
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
//	wstring stemp = wstring(programs.begin(), programs.end());
//	int attr = GetFileAttributes(stemp.c_str());
//	if ((attr & FILE_ATTRIBUTE_HIDDEN) == 0) {
//		SetFileAttributes(stemp.c_str(), attr | FILE_ATTRIBUTE_HIDDEN);
//	}
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
	int dump;
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
				dump = sscanf(tempchar, "%d", &ncpus);
				break;
			case 4: length = line.copy(tempchar,line.size()-3,4);
				tempchar[length]='\0';
				dump = sscanf(tempchar, "%f", &mem);
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

/*bool open_file_dialog(string &path, bool debug, vector <string> filter){
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
*/

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


void progress_bar::write(double fraction)
{
	// clamp fraction to valid range [0,1]
	if (fraction < 0)
		fraction = 0;
	else if (fraction > 1)
		fraction = 1;

	auto width = bar_width - message.size();
	auto offset = bar_width - static_cast<unsigned>(width * round(fraction * 20) * 0.05);

	os << message;
	os.write(full_bar.data() + offset, width);
	os << " [" << std::setw(3) << static_cast<int>(100 * round(fraction*20)*0.05) << "%] " << std::flush;
}

void readxyzMinMax_fromWFN(
	WFN& wavy,
	double* CoordMinMax,
	double* NbSteps,
	double Radius,
	double Increments
)
{
	const double bohrtoa = 0.52917720859;
	int iAtoms, j;
	int NbAtoms = wavy.atoms.size();

	vector < vector <double > > PosAtoms;
	PosAtoms.resize(wavy.get_ncen());
	for (int i = 0; i < wavy.get_ncen(); i++)
		PosAtoms[i].resize(3);
	bool bohrang = !check_bohr(wavy, false, false);
	j = 0;
	for (iAtoms = 0; iAtoms < wavy.get_ncen(); iAtoms++)
	{

		PosAtoms[j][0] = wavy.atoms[iAtoms].x;
		PosAtoms[j][1] = wavy.atoms[iAtoms].y;
		PosAtoms[j][2] = wavy.atoms[iAtoms].z;
		if (!bohrang) {
			cout << "Dividing atom positions!" << endl;
			PosAtoms[j][0] /= bohrtoa;
			PosAtoms[j][1] /= bohrtoa;
			PosAtoms[j][2] /= bohrtoa;
		}
		if (iAtoms == 0)
		{
			CoordMinMax[0] = PosAtoms[j][0];
			CoordMinMax[3] = PosAtoms[j][0];
			CoordMinMax[1] = PosAtoms[j][1];
			CoordMinMax[4] = PosAtoms[j][1];
			CoordMinMax[2] = PosAtoms[j][2];
			CoordMinMax[5] = PosAtoms[j][2];
		}

		else
		{
			if (CoordMinMax[0] > PosAtoms[j][0])
				CoordMinMax[0] = PosAtoms[j][0];
			if (CoordMinMax[3] < PosAtoms[j][0])
				CoordMinMax[3] = PosAtoms[j][0];

			if (CoordMinMax[1] > PosAtoms[j][1])
				CoordMinMax[1] = PosAtoms[j][1];
			if (CoordMinMax[4] < PosAtoms[j][1])
				CoordMinMax[4] = PosAtoms[j][1];

			if (CoordMinMax[2] > PosAtoms[j][2])
				CoordMinMax[2] = PosAtoms[j][2];
			if (CoordMinMax[5] < PosAtoms[j][2])
				CoordMinMax[5] = PosAtoms[j][2];
		}
		j++;
	}

	CoordMinMax[0] -= Radius/bohrtoa;
	CoordMinMax[3] += Radius/bohrtoa;
	CoordMinMax[1] -= Radius/bohrtoa;
	CoordMinMax[4] += Radius/bohrtoa;
	CoordMinMax[2] -= Radius/bohrtoa;
	CoordMinMax[5] += Radius/bohrtoa;

	NbSteps[0] = ceil((CoordMinMax[3] - CoordMinMax[0]) / Increments * bohrtoa);
	NbSteps[1] = ceil((CoordMinMax[4] - CoordMinMax[1]) / Increments * bohrtoa);
	NbSteps[2] = ceil((CoordMinMax[5] - CoordMinMax[2]) / Increments * bohrtoa);
}

void readxyzMinMax_fromCIF(
	string cif,
	double* CoordMinMax,
	double* NbSteps,
	vector < vector < double > > &cm,
	double Resolution,
	ofstream& file,
	bool debug
)
{
	if (debug)
		file << "starting to read cif!" << endl;
	if (!exists(cif)) {
		file << "CIF does not exists!" << endl;
		return;
	}
	ifstream cif_input(cif.c_str(), ios::in);
	vector<bool> found;
	string line;
	found.resize(7);
	for (int k = 0; k < 7; k++)
		found[k] = false;
	double a = 0.0, b = 0.0, c = 0.0, v = 0.0;
	double alpha = 0.0, beta = 0.0, gamma = 0.0;
	vector <string> cell_keywords;
	cell_keywords.push_back("_cell_length_a");
	cell_keywords.push_back("_cell_length_b");
	cell_keywords.push_back("_cell_length_c");
	cell_keywords.push_back("_cell_angle_alpha");
	cell_keywords.push_back("_cell_angle_beta");
	cell_keywords.push_back("_cell_angle_gamma");
	cell_keywords.push_back("_cell_volume");
	if (debug)
		file << "Starting while !.eof()" << endl;
	while (!cif_input.eof()) {
		if (debug)
			file << "While line! " << setw(80) << line
			<< setw(10) << a << found[0]
			<< setw(10) << b << found[1]
			<< setw(10) << c << found[2]
			<< setw(10) << alpha << found[3]
			<< setw(10) << beta << found[4]
			<< setw(10) << gamma << found[5]
			<< setw(10) << v << found[6] << endl;
		getline(cif_input, line);
		for (int k = 0; k < cell_keywords.size(); k++) {
			if (line.find(cell_keywords[k]) != string::npos) {
				switch (k) {
				case 0:
					a = stod(line.substr(cell_keywords[k].length(), line.find("(")));
					break;
				case 1:
					b = stod(line.substr(cell_keywords[k].length(), line.find("(")));
					break;
				case 2:
					c = stod(line.substr(cell_keywords[k].length(), line.find("(")));
					break;
				case 3:
					alpha = stod(line.substr(cell_keywords[k].length(), line.find("(")));
					break;
				case 4:
					beta = stod(line.substr(cell_keywords[k].length(), line.find("(")));
					break;
				case 5:
					gamma = stod(line.substr(cell_keywords[k].length(), line.find("(")));
					break;
				case 6:
					v = stod(line.substr(cell_keywords[k].length(), line.find("(")));
					break;
				default:
					file << "This is weird... should never get here... aborting!" << endl;
					return;
				}
				found[k] = true;
			}
		}
		if (found[0] == true && found[1] == true && found[2] == true && found[3] == true && found[4] == true && found[5] == true && found[6] == true)
			break;
	}
	double ca = cos(0.0174532925199432944444444444444 * alpha);
	double cb = cos(0.0174532925199432944444444444444 * beta);
	double cg = cos(0.0174532925199432944444444444444 * gamma);
	double sa = sin(0.0174532925199432944444444444444 * alpha);
	double sb = sin(0.0174532925199432944444444444444 * beta);
	double sg = sin(0.0174532925199432944444444444444 * gamma);
	double Vp = sqrt((1 - ca * ca - cb * cb - cg * cg) + 2 * (ca * cb * cg));
	double V = a * b * c * Vp;

	if (debug)
		file << "Making cm" << endl
		<< a << " " << b << " " << c << " " << ca << " " << cb << " " << cg << " " << sa << " " << sb << " " << sg << " " << V  << " vs. " << v << endl;

	cm[0][0] = a / 0.529177249;
	cm[1][0] = 0.0;
	cm[2][0] = 0.0;

	cm[0][1] = b * cg / 0.529177249;
	cm[1][1] = b * sg / 0.529177249;
	cm[2][1] = 0.0;

	cm[0][2] = c * cb / 0.529177249;
	cm[1][2] = (c * (ca - cb * cg) / sg) / 0.529177249;
	cm[2][2] = V / (a * b * sg) / 0.529177249;

	CoordMinMax[0] = 0.0;
	CoordMinMax[1] = 0.0;
	CoordMinMax[2] = 0.0;

	CoordMinMax[3] = (a + b * cg + c * cb) / 0.529177249;
	CoordMinMax[4] = (b * sg + c * (ca - cb * cg) / sg) / 0.529177249;
	CoordMinMax[5] = V / (a * b * sg) / 0.529177249;

	NbSteps[0] = ceil(a / Resolution);
	NbSteps[1] = ceil(b / Resolution);
	NbSteps[2] = ceil(c / Resolution);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			cm[i][j] /= NbSteps[j];

	cif_input.close();
}

const int type_vector[168]{ 
	0, 0, 0,
	1, 0, 0,
	0, 1, 0,
	0, 0, 1,
	2, 0, 0,
	0, 2, 0,
	0, 0, 2,
	1, 1, 0,
	1, 0, 1,
	0, 1, 1,
	3, 0, 0,
	0, 3, 0,
	0, 0, 3,
	2, 1, 0,
	2, 0, 1,
	0, 2, 1,
	1, 2, 0,
	1, 0, 2,
	0, 1, 2,
	1, 1, 1,
	0, 0, 4,
	0, 1, 3,
	0, 2, 2,
	0, 3, 1,
	0, 4, 0,
	1, 0, 3,
	1, 1, 2,
	1, 2, 1,
	1, 3, 0,
	2, 0, 2,
	2, 1, 1,
	2, 2, 0,
	3, 0, 1,
	3, 1, 0,
	4, 0, 0,
	0, 0, 5,
	0, 1, 4,
	0, 2, 3,
	0, 3, 2,
	0, 4, 1,
	0, 5, 0,
	1, 0, 4,
	1, 1, 3,
	1, 2, 2,
	1, 3, 1,
	1, 4, 0,
	2, 0, 3,
	2, 1, 2,
	2, 2, 1,
	2, 3, 0,
	3, 0, 2,
	3, 1, 1,
	3, 2, 0,
	4, 0, 1,
	4, 1, 0,
	5, 0, 0};

void type2vector(
	const int index,
	int* vector) {
	if (index < 1 || index > 35) {
		vector[0] = -1;
		vector[1] = -1;
		vector[2] = -1;
		return;
	}
	const int temp = index - 1;
	vector[0] = type_vector[temp * 3]; 
	vector[1] = type_vector[temp * 3 + 1]; 
	vector[2] = type_vector[temp * 3 + 2];
}

bool read_fracs_ADPs_from_CIF(string cif, WFN& wavy, cell &unit_cell, ofstream &log3, bool debug) {
	vector<vector<double>> Uij, Cijk, Dijkl;
	ifstream asym_cif_input(cif.c_str(), std::ios::in);
	asym_cif_input.clear();
	asym_cif_input.seekg(0, asym_cif_input.beg);
	string line;
	vector <string> labels;
	int count_fields = 0;
	int position_field[3] = { 0,0,0 };
	int label_field = 100;
	vector < vector <double > > positions;
	positions.resize(wavy.get_ncen());

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < wavy.get_ncen(); i++) {
		positions[i].resize(3);
	}
	bool atoms_read = false;
	while (!asym_cif_input.eof() && !atoms_read) {
		getline(asym_cif_input, line);
		if (line.find("loop_") != string::npos) {
			while (line.find("_") != string::npos) {
				getline(asym_cif_input, line);
				if (debug) log3 << "line in loop field definition: " << line << endl;
				if (line.find("label") != string::npos)
					label_field = count_fields;
				else if (line.find("fract_x") != string::npos)
					position_field[0] = count_fields;
				else if (line.find("fract_y") != string::npos)
					position_field[1] = count_fields;
				else if (line.find("fract_z") != string::npos)
					position_field[2] = count_fields;
				else if (label_field == 100) {
					if (debug) log3 << "I don't think this is the atom block.. moving on!" << endl;
					break;
				}
				count_fields++;
			}
			while (line.find("_") == string::npos && line.length() > 3) {
				atoms_read = true;
				stringstream s(line);
				vector <string> fields;
				fields.resize(count_fields);
				for (int i = 0; i < count_fields; i++)
					s >> fields[i];
				if (debug) log3 << "label: " << fields[label_field] << " frac_position: " << stod(fields[position_field[0]]) << " " << stod(fields[position_field[1]]) << " " << stod(fields[position_field[2]]) << endl;
				positions[labels.size()] = unit_cell.get_coords_cartesian(stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]));
				bool found_this_one = false;
				if (debug) log3 << "label: " << fields[label_field] << " cartesian position: " << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
				for (int i = 0; i < wavy.get_ncen(); i++) {
					if (is_similar(positions[labels.size()][0], wavy.atoms[i].x, -1)
						&& is_similar(positions[labels.size()][1], wavy.atoms[i].y, -1)
						&& is_similar(positions[labels.size()][2], wavy.atoms[i].z, -1)) {
						if (debug) log3 << "WFN position: " << wavy.atoms[i].x << " " << wavy.atoms[i].y << " " << wavy.atoms[i].z << endl
							<< "Found an atom: " << fields[label_field] << " Corresponding to atom charge " << wavy.atoms[i].charge << endl;
						wavy.atoms[i].label = fields[label_field];
						wavy.atoms[i].frac_coords = { stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]) };
						found_this_one = true;
						break;
					}
				}
				if (!found_this_one && debug)
					log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
				labels.push_back(fields[label_field]);
				getline(asym_cif_input, line);
			}
		}
	}

	asym_cif_input.clear();
	asym_cif_input.seekg(0, asym_cif_input.beg);
	count_fields = 0;
	int ADP_field[15] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	label_field = 100;
	atoms_read = false;
	Uij.resize(wavy.get_ncen());
	while (!asym_cif_input.eof() && !atoms_read) {
		getline(asym_cif_input, line);
		if (line.find("loop_") != string::npos) {
			while (line.find("_") != string::npos) {
				getline(asym_cif_input, line);
				if (debug) log3 << "line in loop field definition: " << line << endl;
				if (line.find("aniso_label") != string::npos)
					label_field = count_fields;
				else if (line.find("aniso_U_11") != string::npos)
					ADP_field[0] = count_fields;
				else if (line.find("aniso_U_22") != string::npos)
					ADP_field[1] = count_fields;
				else if (line.find("aniso_U_33") != string::npos)
					ADP_field[2] = count_fields;
				else if (line.find("aniso_U_12") != string::npos)
					ADP_field[3] = count_fields;
				else if (line.find("aniso_U_13") != string::npos)
					ADP_field[4] = count_fields;
				else if (line.find("aniso_U_23") != string::npos)
					ADP_field[5] = count_fields;
				else if (label_field == 100) {
					if (debug) log3 << "I don't think this is the Uij block.. moving on!" << endl;
					break;
				}
				count_fields++;
			}
			while (line.find("_") == string::npos && line.length() > 3) {
				atoms_read = true;
				stringstream s(line);
				vector <string> fields;
				fields.resize(count_fields);
				for (int i = 0; i < count_fields; i++)
					s >> fields[i];
				if (debug) log3 << "label: " << fields[label_field] << endl;
				bool found_this_one = false;
				for (int i = 0; i < wavy.get_ncen(); i++) {
					if (fields[label_field] == wavy.atoms[i].label) {
						Uij[i].resize(6);
						for (int j = 0; j < 6; j++)
							Uij[i][j] = stod(fields[ADP_field[j]]);
						found_this_one = true;
						break;
					}
				}
				if (!found_this_one && debug)
					log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
				getline(asym_cif_input, line);
			}
		}
	}

	asym_cif_input.clear();
	asym_cif_input.seekg(0, asym_cif_input.beg);
	count_fields = 0;
	label_field = 100;
	atoms_read = false;
	Cijk.resize(wavy.get_ncen());
	while (!asym_cif_input.eof() && !atoms_read) {
		getline(asym_cif_input, line);
		if (line.find("loop_") != string::npos) {
			while (line.find("_") != string::npos) {
				getline(asym_cif_input, line);
				if (debug) log3 << "line in loop field definition: " << line << endl;
				if (line.find("C_label") != string::npos)
					label_field = count_fields;
				else if (line.find("C_111") != string::npos)
					ADP_field[0] = count_fields;
				else if (line.find("C_112") != string::npos)
					ADP_field[1] = count_fields;
				else if (line.find("C_113") != string::npos)
					ADP_field[2] = count_fields;
				else if (line.find("C_122") != string::npos)
					ADP_field[3] = count_fields;
				else if (line.find("C_123") != string::npos)
					ADP_field[4] = count_fields;
				else if (line.find("C_133") != string::npos)
					ADP_field[5] = count_fields;
				else if (line.find("C_222") != string::npos)
					ADP_field[6] = count_fields;
				else if (line.find("C_223") != string::npos)
					ADP_field[7] = count_fields;
				else if (line.find("C_233") != string::npos)
					ADP_field[8] = count_fields;
				else if (line.find("C_333") != string::npos)
					ADP_field[9] = count_fields;
				else if (label_field == 100) {
					if (debug) log3 << "I don't think this is the Cijk block.. moving on!" << endl;
					break;
				}
				count_fields++;
			}
			while (line.find("_") == string::npos && line.length() > 3) {
				atoms_read = true;
				stringstream s(line);
				vector <string> fields;
				fields.resize(count_fields);
				for (int i = 0; i < count_fields; i++)
					s >> fields[i];
				if (debug) log3 << "label: " << fields[label_field] << endl;
				bool found_this_one = false;
				for (int i = 0; i < wavy.get_ncen(); i++) {
					if (fields[label_field] == wavy.atoms[i].label) {
						Cijk[i].resize(10);
						for (int j = 0; j < 6; j++)
							Cijk[i][j] = stod(fields[ADP_field[j]]);
						found_this_one = true;
						break;
					}
				}
				if (!found_this_one && debug)
					log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
				getline(asym_cif_input, line);
			}
		}
	}

	asym_cif_input.clear();
	asym_cif_input.seekg(0, asym_cif_input.beg);
	count_fields = 0;
	label_field = 100;
	atoms_read = false;
	Dijkl.resize(wavy.get_ncen());
	while (!asym_cif_input.eof() && !atoms_read) {
		getline(asym_cif_input, line);
		if (line.find("loop_") != string::npos) {
			while (line.find("_") != string::npos) {
				getline(asym_cif_input, line);
				if (debug) log3 << "line in loop field definition: " << line << endl;
				if (line.find("D_label") != string::npos)
					label_field = count_fields;
				else if (line.find("D_1111") != string::npos)
					ADP_field[0] = count_fields;
				else if (line.find("D_1112") != string::npos)
					ADP_field[1] = count_fields;
				else if (line.find("D_1113") != string::npos)
					ADP_field[2] = count_fields;
				else if (line.find("D_1122") != string::npos)
					ADP_field[3] = count_fields;
				else if (line.find("D_1123") != string::npos)
					ADP_field[4] = count_fields;
				else if (line.find("D_1133") != string::npos)
					ADP_field[5] = count_fields;
				else if (line.find("D_1222") != string::npos)
					ADP_field[6] = count_fields;
				else if (line.find("D_1223") != string::npos)
					ADP_field[7] = count_fields;
				else if (line.find("D_1233") != string::npos)
					ADP_field[8] = count_fields;
				else if (line.find("D_1333") != string::npos)
					ADP_field[9] = count_fields;
				else if (line.find("D_2222") != string::npos)
					ADP_field[10] = count_fields;
				else if (line.find("D_2223") != string::npos)
					ADP_field[11] = count_fields;
				else if (line.find("D_2233") != string::npos)
					ADP_field[12] = count_fields;
				else if (line.find("D_2333") != string::npos)
					ADP_field[13] = count_fields;
				else if (line.find("D_3333") != string::npos)
					ADP_field[14] = count_fields;
				else if (label_field == 100) {
					if (debug) log3 << "I don't think this is the Dijk block.. moving on!" << endl;
					break;
				}
				count_fields++;
			}
			while (line.find("_") == string::npos && line.length() > 3) {
				atoms_read = true;
				stringstream s(line);
				vector <string> fields;
				fields.resize(count_fields);
				for (int i = 0; i < count_fields; i++)
					s >> fields[i];
				if (debug) log3 << "label: " << fields[label_field] << endl;
				bool found_this_one = false;
				for (int i = 0; i < wavy.get_ncen(); i++) {
					if (fields[label_field] == wavy.atoms[i].label) {
						Dijkl[i].resize(15);
						for (int j = 0; j < 6; j++)
							Dijkl[i][j] = stod(fields[ADP_field[j]]);
						found_this_one = true;
						break;
					}
				}
				if (!found_this_one && debug)
					log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
				getline(asym_cif_input, line);
			}
		}
	}

	for (int i = 0; i < wavy.get_ncen(); i++) 
		wavy.atoms[i].assign_ADPs(Uij[i], Cijk[i], Dijkl[i]);

	return true;
};

