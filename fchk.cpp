#include <iostream>
#include <fstream>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#include <sys/wait.h>
#endif
#include <sys/stat.h>
#include <fcntl.h>
#include <string>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <sys/types.h>
#include <vector>

#include "fchk.h"
#include "convenience.h"
#include "basis_set.h"
using namespace std;
bool debug_fchk=false;
//----------------------------FCHK Preparation and Gaussian--------------------------------------

bool chk2fchk(const string &outputname, const string &chk_name, const string &gaussian_path){
	bool success=false;
	string fchk_path=gaussian_path;
	size_t gaussian_directory=fchk_path.find("g09/g09");
	fchk_path.erase(fchk_path.find("g09/g09"),7);
	fchk_path.append("g09/formchk");
	string execute_command;
	execute_command=fchk_path;
	execute_command.append(" ");
	execute_command.append(chk_name);
	execute_command.append(" ");
	execute_command.append(outputname);
	if(system(execute_command.c_str())) cout << "Well system returned non 0...! Don't ask me what this means now..." << endl;
	return true;
};

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
		cout << envbsd << endl;
		cout << envlocal << endl;
		cout << envextras << endl;
		cout << envshort << endl;
	} 
	string EXEDIR=("GAUSS_EXEDIR=");
	EXEDIR.append(envbsd);
	EXEDIR.append(":");
	EXEDIR.append(envlocal);
	EXEDIR.append(":");
	EXEDIR.append(envextras);
	EXEDIR.append(":");
	EXEDIR.append(envshort);
	if(debug_fchk) cout << EXEDIR << endl;
	string BSDDIR=("GAUSS_BSDDIR=");
	BSDDIR.append(envbsd);
	if(debug_fchk) cout << BSDDIR << endl;
	char *exedir = &EXEDIR[0u];
	char *bsddir = &BSDDIR[0u];
	cls();
	cout << "Running gaussian... Please wait... " << endl;
	string execute_command;
	execute_command=workpath;
	execute_command.append(" gaussian.com");
	if(system(execute_command.c_str())) cout << "Well system returned non 0... No idea what this means now..." << endl;
	if(!debug) cout << "Finished gaussian!!" << endl; 
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

string prepare_gaussian(const string &basis_set_path, const string &fchkname, WFN &wave, const int &ncpus, const float &mem, bool debug){
	if(exists("gaussian.com")){
		cout << "gaussian.com already exists, do you want me to overwrite it?";
		if(!yesno()) return "WRONG";
	}
	if(wave.get_modified ()){
		cout << "The wavefunction has been modified after reading. Please make sure that you know what you are doing!" << endl;
		cout << "Do you want to continue?" << flush;
		if(!yesno()) return "WRONG";
	}
	//-----------------READING BASIS SET--------------------------------
	if(wave.get_nr_basis_set_loaded ()==0){
		cout << "No basis set loaded, will load a complete basis set now!" << endl;
		if(debug) Enter();
		if(!read_basis_set_vanilla (basis_set_path, wave, debug,false)) {
			cout << "Problem during reading of the basis set!" << endl;
			return "WRONG";
		}
	}
	else if(wave.get_nr_basis_set_loaded ()<wave.get_ncen ()){
		cout << "Not all atoms have a basis set loaded!" << endl
			<< "Do you want to laod the missing atoms?" << flush;
		if(!yesno()){
			cout << "Do you want to load a new basis set?" <<flush;
			if(!yesno()){
				cout << "Okay, aborting then!!!" << endl;
				if(debug) Enter();
				return "WRONG";
			}
			cout << "deleting old one..." << endl;
			if(debug) Enter();
			if(!delete_basis_set_vanilla (basis_set_path, wave, debug)){
				cout << "ERROR while deleting a basis set!" << endl;
				return "WRONG";
			}
			else if(!read_basis_set_vanilla (basis_set_path, wave, debug,true)){
				cout << "ERROR during reading of the new basis set!" << endl;
				return "WRONG";
			}
		}
		else{
			cout << "okay, loading the missing atoms.." << endl;
			if(debug) Enter();
			if(!read_basis_set_missing(basis_set_path, wave, debug)){
				cout << "ERROR during reading of missing basis set!" << endl;
				return "WRONG";
			}
		}
	}
	else if(wave.get_nr_basis_set_loaded ()==wave.get_ncen ()){
		cout << "There already is a basis set loaded!" << endl
			<< "Do you want me to delete the old one and load a new one?" << flush;
		if(debug) Enter();
		if(!yesno()){
			cout << "Okay, continuing with the old one..." << endl;
			if(debug) Enter();
		}
		else{
			cout << "Deleting the old basis set!" << endl;
			if(debug) Enter();
			if(!delete_basis_set_vanilla (basis_set_path, wave, debug)){
				cout << "ERROR during deleting of the basis set!";
				Enter();
				return "WRONG";
			}
			cout << "Going to load a new one now!" << endl;
			if(!read_basis_set_vanilla (basis_set_path, wave, debug,true)) {
				cout << "Problem during reading of the basis set!" << endl;
				return "WRONG";
			}
		}
	}
	else{
		cout << "# of loaded > # atoms" << endl
			<< "Sorry, this should not happen... aborting!!!" << endl;
		return "WRONG";
	}
	//-----------------------check ordering and order accordingly----------------------
	if(!wave.sort_wfn(wave.check_order(debug),debug)){
		cout << "Could not order the wavefunction, aborting!" << endl;
		if(debug) Enter();
		return "WRONG";
	}
	//------------setting up the .com file with basis set---------------
	ofstream com;
	com.open("gaussian.com",ios::out);
	string temp;
	temp="gaussian.chk";
	int run=0;
	if(exists(temp)){
		cout << "Do you want me to overwrite gaussian.chk?";
		if(!yesno()){
			cout << "Then pelase give a new .chk filename: ";
			string filename;
			cin >> filename;
			if(filename.find(".chk")==-1){
				cout << "this does not look like a .chk file! Please try again!";
				return "WRONG";
			}
		}
	}
	com << "%chk="<< temp << endl;
	com << "%mem=" << mem << "GB" << endl;
	com << "%nproc=" << ncpus << endl;
/*	if(wave.get_d_f_switch()) com << "# SCF=(MaxCycle=3000,Conver=-1) b3lyp/gen 5D 7F nosymm IOp(3/32=2)" << endl << endl;
	else*/
	com << "# SCF=(MaxCycle=3000,Conver=-1) rhf/gen 6D 10F nosymm IOp(3/32=2)" << endl << endl;
	com << "TITLE" << endl << endl;
	com << wave.get_charge() << " " << wave.get_multi() << endl;
	com << wave.get_centers (check_bohr (wave,true,debug));
//	com << wave.get_centers (false);
	com << endl;
	vector<string> elements_list;
	if(debug) cout << "elements_list.size()= " << elements_list.size() << endl;
	for(int a=0; a<wave.get_ncen(); a++){
		string label_temp;
		label_temp=wave.get_atom_label(a);
		bool found=false;
		if(elements_list.size()!=0){
			for (int i=0; i< elements_list.size(); i++) {
				if(elements_list[i].compare(label_temp)==0){
					found=true;
					if(debug) cout << "Found " << label_temp << " in the elements_list!" << endl;
				}
			}
		}
		if(!found){
			elements_list.push_back(label_temp);
			if(debug) cout << "elements_list.size()= " << elements_list.size() << endl;
			com << label_temp << " 0" << endl;
			int working_shell=-1;
			for (int p=0; p< wave.get_atom_primitive_count (a);p++){
				int temp_shell=wave.get_basis_set_shell (a,p);
				if(temp_shell!=working_shell){
					working_shell=temp_shell;
					switch (wave.get_shell_type (a,temp_shell)){
						case 1:
							com << "S " << wave.get_atom_shell_primitives (a,temp_shell) << " 1.00" << endl;
							break;
						case 2:
							com << "P " << wave.get_atom_shell_primitives (a,temp_shell) << " 1.00" << endl;
							break;
						case 3:
							com << "D " << wave.get_atom_shell_primitives (a,temp_shell) << " 1.00" << endl;
							break;
						case 4:
							com << "F " << wave.get_atom_shell_primitives (a,temp_shell) << " 1.00" << endl;
							break;
					}
				}
				com << scientific << setw(17) << setprecision(10) << wave.get_atom_basis_set_exponent (a,p);
				com << " ";
				com << scientific << setw(17) << setprecision(10) << wave.get_atom_basis_set_coefficient (a,p);
				com << endl;
		    }
			com << "****" << endl;
		}
	}
	com << endl;
	com.flush();
	com.close();
	if(debug){
		cout << "Wrote the gaussian Input! Please check it before i continue..." << endl;
		Enter();
	}
	return temp;
};

bool new_fchk(const string &basis_set_path, const std::string &fchkname, 
              const std::string &gaussian_path, WFN &wave, const int &ncpus, const float &mem, bool debug){
	if(debug) debug_fchk=true;
	if(debug_fchk) cout << "let's run a gaussian calculation to obtain a fchk file" << endl;
	string chkfile=prepare_gaussian (basis_set_path, fchkname, wave, ncpus, mem, debug);
	if(chkfile=="WRONG"){
		cout << "ERROR during Preparation of gaussian, aborting!" << endl;
		if(debug) Enter();
		return false;
	}
	bool success=false;
	success=gaussian(gaussian_path,debug);
	if(success){
		success=chk2fchk(fchkname,chkfile,gaussian_path);
		if(debug){
			cout << "made the fchk!" << endl;
			Enter();
		}
	}
	else return false;
	if(debug){
		cout << "I will keep the gaussian files which have been created temporarily!!" << endl;
		Enter();
	}
	if(success&&!debug_fchk){
		remove("gaussian.com");
		remove("gaussian.log");
		remove("gaussian_run.log");
		remove("gaussian.chk");
		remove("formchk.log");
	}
	return success;
};

bool modify_fchk(const string &fchk_name, const string &basis_set_path, WFN &wave, bool &debug, const bool &read){
	wave.set_modified();
	if(debug) debug_fchk=true;
	vector<double> CMO;
	int nao=0;
	int nshell=0;
	int naotr=0;
	if(debug){
		cout << "Origin: " << wave.get_origin() << endl;
		Enter();
	}
	if(wave.get_origin ()==2||wave.get_origin()==4){
		//-----------------READING BASIS SET--------------------------------
		if(read){
			if(!read_basis_set_vanilla (basis_set_path, wave, debug,false)) {
				cout << "Problem during reading of the basis set!" << endl;
				return false;
			}
		}
		double pi=3.14159265358979;
		//---------------normalize basis set---------------------------------
		if(debug) cout << "starting to normalize the basis set" << endl;
		vector<double> norm_const;
		//-----------debug output---------------------------------------------------------
		if(debug){
			cout << "exemplary output before norm_const of the first atom with all it's properties: "<< endl;
			wave.print_atom_long(0);
			cout << "ended normalizing the basis set, now for the MO_coeffs" << endl;
			cout << "Status report:" << endl;
			cout << "size of norm_const: " << norm_const.size() << endl;
			cout << "WFN MO counter: " << wave.get_nmo() << endl;
			cout << "Number of atoms: " << wave.get_ncen() << endl;
			cout << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
			cout << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
			Enter();
		}
		//-----------------------check ordering and order accordingly----------------------
		wave.sort_wfn(wave.check_order(debug),debug);
		//-------------------normalize the basis set shell wise---------------------
		for(int a=0; a< wave.get_ncen(); a++){
			for(int p=0; p< wave.get_atom_primitive_count(a); p++){
				double temp=wave.get_atom_basis_set_exponent(a,p);
				switch(wave.get_atom_primitive_type(a,p)){
					case 1:
						temp=2*temp/pi;
						temp=pow(temp,0.75);
						temp=temp*wave.get_atom_basis_set_coefficient(a,p);
						wave.change_atom_basis_set_coefficient(a,p,temp);
						break;
					case 2:
						temp=128*pow(temp,5);
						temp=temp/pow(pi,3);
						temp=pow(temp,0.25);
						temp=wave.get_atom_basis_set_coefficient(a,p)*temp;
						wave.change_atom_basis_set_coefficient(a,p,temp);
						break;
					case 3:
						temp=2048*pow(temp,7);
						temp=temp/(9*pow(pi,3));
						temp=pow(temp,0.25);
						temp=wave.get_atom_basis_set_coefficient(a,p)*temp;
						wave.change_atom_basis_set_coefficient(a,p,temp);
						break;
					case 4:
						temp=32768*pow(temp,9);
						temp=temp/(225*pow(pi,3));
						temp=pow(temp,0.25);
						temp=wave.get_atom_basis_set_coefficient(a,p)*temp;
						wave.change_atom_basis_set_coefficient(a,p,temp);
						break;
					case -1:
						cout << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy..." << endl;
						Enter();
						break;
				}
			}
		}
		for(int a=0; a< wave.get_ncen(); a++){
			double aiaj=0.0;
			double factor=0.0;
			for(int s=0; s< wave.get_atom_shell_count(a); s++){
				int type_temp=wave.get_shell_type(a,s);		
				if(type_temp==-1){
					cout << "ERROR in type assignement!!"<<endl;
					Enter();
				}
				if (debug){
					cout << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << endl;
					int testcount=0;
					cout << "start: " << wave.get_shell_start(a,s,false) << flush;
					cout <<" stop: " << wave.get_shell_end(a,s,false) << flush << endl;
					cout << "factor: ";
				}
				switch (type_temp){
					case 1:
						factor=0;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							for(int j=wave.get_shell_start(a,s,false); j<= wave.get_shell_end(a,s,false); j++){
								double aiaj=wave.get_atom_basis_set_exponent (a,i)+wave.get_atom_basis_set_exponent (a,j);
								double term=(pi/aiaj);
								term=pow(term,1.5);
								factor+=wave.get_atom_basis_set_coefficient(a,i)*wave.get_atom_basis_set_coefficient(a,j)*term;
							}
						}
						if(factor==0) return false;
						factor=pow(factor,-0.5);
						if(debug) cout << factor << endl;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							if(debug){
								cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a,i) << endl;
								cout << "Contraction coefficient after:  " << factor*wave.get_atom_basis_set_coefficient(a,i) << endl;
							}
							wave.change_atom_basis_set_coefficient(a,i,factor*wave.get_atom_basis_set_coefficient(a,i));
							norm_const.push_back(wave.get_atom_basis_set_coefficient(a,i));
						}
						break;
					case 2:
						factor=0;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							for(int j=wave.get_shell_start(a,s,false); j<= wave.get_shell_end(a,s,false); j++){
								double aiaj=wave.get_atom_basis_set_exponent (a,i)+wave.get_atom_basis_set_exponent (a,j);
								double term=4*pow(aiaj,5);
								term=pow(pi,3)/term;
								term=pow(term,0.5);
								factor+=wave.get_atom_basis_set_coefficient(a,i)*wave.get_atom_basis_set_coefficient(a,j)*term;
							}
						}
						if(factor==0) return false;
						factor=pow(factor,-0.5);
						if(debug) cout << factor << endl;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							if(debug){
								cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a,i) << endl;
								cout << "Contraction coefficient after:  " << factor*wave.get_atom_basis_set_coefficient(a,i) << endl;
							}
							wave.change_atom_basis_set_coefficient(a,i,factor*wave.get_atom_basis_set_coefficient(a,i));
							for(int k=0; k<3; k++) norm_const.push_back(wave.get_atom_basis_set_coefficient(a,i));
						}
						break;
					case 3:
						factor=0;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							for(int j=wave.get_shell_start(a,s,false); j<= wave.get_shell_end(a,s,false); j++){
								double aiaj=wave.get_atom_basis_set_exponent (a,i)+wave.get_atom_basis_set_exponent (a,j);
								double term=16*pow(aiaj,7);
								term=pow(pi,3)/term;
								term=pow(term,0.5);
								factor+=wave.get_atom_basis_set_coefficient(a,i)*wave.get_atom_basis_set_coefficient(a,j)*term;
							}
						}
						if(factor==0) return false;
						factor=(pow(factor,-0.5))/sqrt(3);
						if(debug) cout << factor << endl;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							if(debug){
								cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a,i) << endl;
								cout << "Contraction coefficient after:  " << factor*wave.get_atom_basis_set_coefficient(a,i) << endl;
							}
							wave.change_atom_basis_set_coefficient(a,i,factor*wave.get_atom_basis_set_coefficient(a,i));
							for(int k=0; k<3; k++) norm_const.push_back(wave.get_atom_basis_set_coefficient(a,i));
							for(int k=0; k<3; k++) norm_const.push_back(sqrt(3)*wave.get_atom_basis_set_coefficient(a,i));
						}
						break;
					case 4:
						factor=0;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							for(int j=wave.get_shell_start(a,s,false); j<= wave.get_shell_end(a,s,false); j++){
								double aiaj=wave.get_atom_basis_set_exponent (a,i)+wave.get_atom_basis_set_exponent (a,j);
								double term=64*pow((aiaj),9);
								term=pow(pi,3)/term;
								term=pow(term,0.5);
								factor+=wave.get_atom_basis_set_coefficient(a,i)*wave.get_atom_basis_set_coefficient(a,j)*term;
							}
						}
						if(factor==0) return false;
						factor=pow(factor,-0.5)/sqrt(15);
						if(debug) cout << factor << endl;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							if(debug){
								cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a,i) << endl;
								cout << "Contraction coefficient after:  " << factor*wave.get_atom_basis_set_coefficient(a,i) << endl;
							}
							wave.change_atom_basis_set_coefficient(a,i,factor*wave.get_atom_basis_set_coefficient(a,i));
							for(int l=0; l<3;l++) norm_const.push_back(wave.get_atom_basis_set_coefficient(a,i));
							for(int l=0; l<6;l++) norm_const.push_back(sqrt(5)*wave.get_atom_basis_set_coefficient(a,i));
							norm_const.push_back(sqrt(15)*wave.get_atom_basis_set_coefficient(a,i));
						} 
						break;
				}
				if (debug)	cout << "This shell has: " <<  wave.get_shell_end(a,s,false)-wave.get_shell_start(a,s,false)+1 << " primitives" << endl;
			}
		}
		//-----------debug output---------------------------------------------------------
		if(debug){
			cout << "exemplary output of the first atom with all it's properties: "<< endl;
			wave.print_atom_long(0);
			cout << "ended normalizing the basis set, now for the MO_coeffs" << endl;
			cout << "Status report:" << endl;
			cout << "size of norm_const: " << norm_const.size() << endl;
			cout << "WFN MO counter: " << wave.get_nmo() << endl;
			cout << "Number of atoms: " << wave.get_ncen() << endl;
			cout << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
			cout << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
			Enter();
			Enter();
		}
		//---------------------To not mix up anything start normalizing WFN_matrix now--------------------------
		int run=0;
		ofstream norm_cprim;
		if(debug) norm_cprim.open("norm_prim.debug", ofstream::out);
		for(int m=0; m< wave.get_nmo(); m++){
			if(debug) norm_cprim << m << ". MO:" << endl;
			for(int p=0; p<wave.get_MO_primitive_count(m);p++){
				if(debug) cout << p << ". primitive; " << m << ". MO " << "norm nonst: " << norm_const[p] << endl;
				double temp=wave.get_MO_coef(m,p, debug)/norm_const[p];
				if(debug){
					cout << " temp after normalization: " << temp << endl;
					norm_cprim << " " << temp << endl;
				}
				run++;
				if(!wave.change_MO_coef(m,p,temp, debug)){
					cout << "ERROR in changing the coefficients after normalising!";
					if(debug) cout << "m:" << m << " p: " << p << " temp:" << temp; 
					cout << endl;
					Enter();
					return false;
				}
			}
		}
		if(debug) {
			norm_cprim.flush();
			norm_cprim.close();
			cout << "See norm_cprim.debug for the CPRIM vectors" << endl;
			cout << "Total count in CPRIM: " << run << endl;
			Enter();
			Enter();
		}
		//--------------Build CMO of alessandro from the first elements of each shell-------------
		for (int m =0; m < wave.get_nmo (); m++){
			int run_2=0;
			for (int a =0; a < wave.get_ncen (); a++){
				for (int s=0; s< wave.get_atom_shell_count(a); s++){
					if(debug) cout << "Going to load the " << wave.get_shell_start_in_primitives(a,s) << ". value" << endl;
					switch (wave.get_shell_type(a,s)){
						case 1:
							CMO.push_back(wave.get_MO_coef(m,wave.get_shell_start_in_primitives(a,s), false));
							if(m==0) nao++;
							if(debug&&wave.get_atom_shell_primitives(a,s)!=1)
								cout << "Pushing back 1 coefficient for S shell, this shell has " << wave.get_atom_shell_primitives(a,s) << " primitives! Shell start is: " << wave.get_shell_start(a,s,false) << endl;
						break;  
						case 2:
							for (int i=0; i<3 ; i++) CMO.push_back(wave.get_MO_coef(m,wave.get_shell_start_in_primitives(a,s)+i, false));
							if(debug&&wave.get_atom_shell_primitives(a,s)!=1) 
								cout << "Pushing back 3 coefficients for P shell, this shell has " << wave.get_atom_shell_primitives(a,s) << " primitives!" << endl;
							if(m==0) nao+=3;
						break;
						case 3:
							for (int i=0; i<6 ; i++) CMO.push_back(wave.get_MO_coef(m,wave.get_shell_start_in_primitives(a,s)+i,false));
							if(debug&&wave.get_atom_shell_primitives(a,s)!=1) 
								cout << "Pushing back 6 coefficient for D shell, this shell has " << wave.get_atom_shell_primitives(a,s) << " primitives!" << endl;
							if(m==0) nao+=6;
						break;
						case 4:
							//this hardcoded piece is due to the order of f-type functions in the fchk
							for (int i=0; i<3 ; i++) CMO.push_back(wave.get_MO_coef(m,wave.get_shell_start_in_primitives(a,s)+i,false));
							CMO.push_back(wave.get_MO_coef(m,wave.get_shell_start_in_primitives(a,s)+6,false));
							for (int i=0; i<2 ; i++) CMO.push_back(wave.get_MO_coef(m,wave.get_shell_start_in_primitives(a,s)+i+3,false));
							for (int i=0; i<2 ; i++) CMO.push_back(wave.get_MO_coef(m,wave.get_shell_start_in_primitives(a,s)+i+7,false));
							CMO.push_back(wave.get_MO_coef(m,wave.get_shell_start_in_primitives(a,s)+5,false));
							CMO.push_back(wave.get_MO_coef(m,wave.get_shell_start_in_primitives(a,s)+9,false));
							if(debug&&wave.get_atom_shell_primitives(a,s)!=1) 
								cout << "Pushing back 10 coefficient for F shell, this shell has " << wave.get_atom_shell_primitives(a,s) << " primitives!" << endl;
							if(m==0) nao+=10;
						break;
					}
					run_2++;
				}
				if(debug) cout << "finished with atom!" << endl;
			}
			if(debug) cout << "finished with MO!" << endl;
			if(nshell!=run_2) nshell=run_2;
		}
		if(debug){
			Enter();
			Enter();
			ofstream cmo ("cmo.debug", ofstream::out);
			for(int p=0; p< CMO.size() ; p++){
				string temp;
				for( int i=0; i< 5; i++){
					stringstream stream;
					stream << scientific << setw(14) << setprecision(7) << CMO[p+i] << " ";
					temp += stream.str();
				}
				p+=4;
				cmo << temp << endl;
			}
			cmo.flush();
			cmo.close();
			cout << CMO.size() << " Elements in CMO" << endl;
			cout << norm_const.size() << " = nprim" << endl;
			cout << nao << " = nao" << endl;
			cout << nshell << " = nshell" << endl;
			Enter();
		}
		//------------------ make the DM -----------------------------
		naotr=nao*(nao+1)/2;
		vector<double> kp;
		for(int i=0; i < naotr; i++) wave.push_back_DM (0.0);
		if(debug) cout << "I made kp!" << endl << nao << " is the maximum for iu" << endl;
		for(int iu=0; iu < nao; iu++){
			for(int iv=0; iv <= iu; iv++){
				int iuv=(iu*(iu+1)/2)+iv;
				if(debug) cout << "iu: " << iu << " iv: " << iv << " iuv: " << iuv << " kp(iu): " << iu*(iu+1)/2<< endl;
				for(int m=0; m < wave.get_nmo(); m++){
					if(debug) cout << "DM before: " << wave.get_DM(iuv) << endl;
					if(!wave.set_DM(iuv,wave.get_DM(iuv)+2*CMO[iu+(m*nao)]*CMO[iv+(m*nao)])){
						cout << "Something went wrong while writing the DM! iuv=" << iuv << endl;
						cout << "accessing CMO elements: " << iu+(m*nao) << " and " << iv+(m*nao) << " while CMO has size: " << CMO.size() << endl;
						cout << "CMO[" << iu+(m*nao) << "] = " << CMO[iu*(m*nao)] << "; CMO[" << iv+(m*nao) << "] = " << CMO[iv+(m*nao)] << " DM: " << wave.get_DM(iuv) << endl;
						cout << "DM_size: " << wave.get_DM_size() << endl;
						Enter();
						return false;
					}
					else if (debug) cout << "DM after: " << wave.get_DM(iuv) << endl;
				}
			}
		}
		if(debug){
			cout << "DM is:" << endl;
			for(int p=0; p< wave.get_DM_size() ; p++){
				string temp;
				for( int i=0; i< 5; i++){
					stringstream stream;
					stream << scientific << setw(14) << setprecision(7) << wave.get_DM(i+p) << " ";
					temp += stream.str();
				}
				p+=4;
				cout << temp << endl;
			}
			cout << wave.get_DM_size() << " Elements in DM" << endl;
			Enter();
		}
	}
	// open fchk for copying
	string temp_fchk=fchk_name;
	temp_fchk.append(".fchk");
	if(!exists(temp_fchk)){
		cout << "This is worrysome... The fchk should be there.." << endl;
		cout << "fchk_name: " << temp_fchk << endl;
		Enter();
		return false;
	}
	ifstream ifchk(temp_fchk.c_str());
	if(!ifchk.is_open()){
		cout << "ERROR while opening .fchk ifile!" << endl;
		Enter();
		return false;
	}
	ofstream ofchk;
	ofchk.open("temp.fchk",ofstream::out);
	string line;
	int dum_nao=0;
	if(wave.get_origin()==2||wave.get_origin()==4){
		while(line.find("Alpha Orbital Energies")==-1 &&!ifchk.eof()){
			getline(ifchk,line);
			if(debug) cout << "line: " << line << endl;
			if(line.find("Alpha Orbital Energies")==-1) ofchk << line << endl;
			else{
				char tempchar[100];
				size_t length;
				length = line.copy(tempchar,11,50);
				tempchar[length]='\0';
				sscanf(tempchar,"%d",&dum_nao);
				if(debug){
					cout << "nao read from fchk: " << dum_nao << " and from basis set: " << nao << endl;
					Enter();
				}
				ofchk << line << endl;
			}
		}
		int counter=0;
		while(line.find("Alpha MO")==-1 && !ifchk.eof()){	
			getline(ifchk,line);
			string temp=" ";
			for(int j=0; j<5; j++){
				if(counter+j<wave.get_nmo()){
					stringstream stream;
					stream << scientific << setw(15) << setprecision(8) << wave.get_MO_energy (counter+j);
					if(j<4) stream << " ";
					temp += stream.str();
					temp.replace(12+j*16,1,"E");
				}
				else if(counter+j<dum_nao){
					double dum_ener=0.0;
					char tempchar[100];
					size_t length;
					length = line.copy(tempchar,15,1+j*16);
					tempchar[length]='\0';
					sscanf(tempchar,"%lf",&dum_ener);
					stringstream stream;
					stream << scientific << setw(15) << setprecision(8) << dum_ener; 
					if(j<4) stream << " ";
					temp += stream.str();
					temp.replace(12+j*16,1,"E");
				}
			}
			counter+=5;
			temp+='\n';
			if(temp.size()>3) ofchk << temp;
		}
		ofchk.flush();   
		ofchk << "Alpha MO coefficients                      R   N=" << setw(12) << nao*dum_nao << endl;
		//now write the CMO and skip lines in IFCHK
		for(int i=0; i<nao*dum_nao; i++){
			string temp=" ";
			for(int j=0; j<5; j++){
				if(i+j<CMO.size()){
					stringstream stream;
					stream << scientific << setw(15) << setprecision(8) << CMO[i+j] << " ";
					temp += stream.str();
					temp.replace(12+j*16,1,"E");
				}
				else if(i+j<nao*nao){
					stringstream stream;
					stream << scientific << setw(15) << setprecision(8) << 0.0 << " ";
					temp += stream.str();
					temp.replace(12+j*16,1,"E");
				}
			}
			i+=4;
			temp+='\n';
			ofchk << temp;
			getline(ifchk,line);
		}
	}
	if(wave.get_origin()==1){
		while(line.find("Total SCF Density")==-1 && !ifchk.eof()){
			getline(ifchk,line);
			if(debug) cout << "line: " << line << endl;
			if(line.find("Total SCF Density")==-1) ofchk << line << endl;
		}
		ofchk.flush();
	}
	ofchk << "Total SCF Density                          R   N=" << setw(12) << wave.get_DM_size() << endl;
	getline(ifchk,line);
	//now write the DM and skip lines in IFCHK
	for(int i=0; i<wave.get_DM_size(); i++){
		string temp=" ";
		if(debug) cout << "i: " << i << " DM_size= " << wave.get_DM_size() << " Element ";
		for(int j=0; j<5; j++){
			if(i+j<wave.get_DM_size ()){
				stringstream stream;
				if(debug) cout << i+j << " ";
				stream << scientific << setw(15) << setprecision(8) << wave.get_DM(i+j) << " ";
				if(i+j<wave.get_DM_size ()){
					temp += stream.str();
					temp.replace(12+j*16,1,"E");
				}
			}
		}
		i+=4;
		if(debug) cout << endl;
		temp+='\n';
		ofchk << temp;
		getline(ifchk,line);
	}
	if(debug) Enter();
	while(!ifchk.eof()){
		getline(ifchk,line);
		ofchk << line ;
		if(!ifchk.eof()) ofchk << endl ;
	}
	ofchk.flush();
	ofchk.close();
	ifchk.close();
	if(debug){
		cout << "Do you want me to keep the old fchk file and the new temp.fchk?";
		if(!yesno()){
			copy_file("temp.fchk",temp_fchk);
		}
	}
	else copy_file("temp.fchk",temp_fchk);
	if(remove("temp.fchk")!=0) cout << "error deleting temp.fchk!" << endl;
	cls();
	return true;
};

bool free_fchk(const string &fchk_name, const string &basis_set_path, WFN &wave, bool &debug, bool force_overwrite){
	if(wave.get_nr_basis_set_loaded ()==0){
		cout << "No basis set loaded, will load a complete basis set now!" << endl;
		if(debug) Enter();
		if(!read_basis_set_vanilla (basis_set_path, wave, debug,false)) {
			cout << "Problem during reading of the basis set!" << endl;
			return "WRONG";
		}
	}
	else if(wave.get_nr_basis_set_loaded ()<wave.get_ncen ()){
		cout << "Not all atoms have a basis set loaded!" << endl
			<< "Do you want to laod the missing atoms?" << flush;
		if(!yesno()){
			cout << "Do you want to load a new basis set?" <<flush;
			if(!yesno()){
				cout << "Okay, aborting then!!!" << endl;
				if(debug) Enter();
				return "WRONG";
			}
			cout << "deleting old one..." << endl;
			if(debug) Enter();
			if(!delete_basis_set_vanilla (basis_set_path, wave, debug)){
				cout << "ERROR while deleting a basis set!" << endl;
				return "WRONG";
			}
			else if(!read_basis_set_vanilla (basis_set_path, wave, debug,true)){
				cout << "ERROR during reading of the new basis set!" << endl;
				return "WRONG";
			}
		}
		else{
			cout << "okay, loading the missing atoms.." << endl;
			if(debug) Enter();
			if(!read_basis_set_missing(basis_set_path, wave, debug)){
				cout << "ERROR during reading of missing basis set!" << endl;
				return "WRONG";
			}
		}
	}
	else if(wave.get_nr_basis_set_loaded ()==wave.get_ncen ()){
		cout << "There already is a basis set loaded!" << endl
			<< "Do you want me to delete the old one and load a new one?" << flush;
		if(debug) Enter();
		if(!yesno()){
			cout << "Okay, continuing with the old one..." << endl;
			if(debug) Enter();
		}
		else{
			cout << "Deleting the old basis set!" << endl;
			if(debug) Enter();
			if(!delete_basis_set_vanilla (basis_set_path, wave, debug)){
				cout << "ERROR during deleting of the basis set!";
				Enter();
				return "WRONG";
			}
			cout << "Going to load a new one now!" << endl;
			if(!read_basis_set_vanilla (basis_set_path, wave, debug,true)) {
				cout << "Problem during reading of the basis set!" << endl;
				return "WRONG";
			}
		}
	}
	else{
		cout << "# of loaded > # atoms" << endl
			<< "Sorry, this should not happen... aborting!!!" << endl;
		return "WRONG";
	}
	//wave.set_modified();
	if(debug) debug_fchk=true;
	vector<double> CMO;
	int nao=0;
	int nshell=0;
	int naotr=0;
	if(debug){
		cout << "Origin: " << wave.get_origin() << endl;
		Enter();
	}
	if(wave.get_origin ()==2||wave.get_origin()==4){
		double pi=3.14159265358979;
		//---------------normalize basis set---------------------------------
		if(debug) cout << "starting to normalize the basis set" << endl;
		vector<double> norm_const;
		//-----------debug output---------------------------------------------------------
		if(debug){
			cout << "exemplary output before norm_const of the first atom with all it's properties: "<< endl;
			wave.print_atom_long(0);
			cout << "ended normalizing the basis set, now for the MO_coeffs" << endl;
			cout << "Status report:" << endl;
			cout << "size of norm_const: " << norm_const.size() << endl;
			cout << "WFN MO counter: " << wave.get_nmo() << endl;
			cout << "Number of atoms: " << wave.get_ncen() << endl;
			cout << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
			cout << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
			Enter();
		}
		//-----------------------check ordering and order accordingly----------------------
		wave.sort_wfn(wave.check_order(debug),debug);
		//-------------------normalize the basis set shell wise into a copy vector---------
		vector <vector <double> > basis_coefficients;
		basis_coefficients.resize(wave.get_ncen());
		for(int a=0; a< wave.get_ncen(); a++)
			for(int p=0; p< wave.get_atom_primitive_count(a); p++){
				double temp=wave.get_atom_basis_set_exponent(a,p);
				switch(wave.get_atom_primitive_type(a,p)){
					case 1:
						temp=2*temp/pi;
						temp=pow(temp,0.75);
						temp=temp*wave.get_atom_basis_set_coefficient(a,p);
						basis_coefficients[a].push_back(temp);
						break;
					case 2:
						temp=128*pow(temp,5);
						temp=temp/pow(pi,3);
						temp=pow(temp,0.25);
						temp=wave.get_atom_basis_set_coefficient(a,p)*temp;
						basis_coefficients[a].push_back(temp);
						break;
					case 3:
						temp=2048*pow(temp,7);
						temp=temp/(9*pow(pi,3));
						temp=pow(temp,0.25);
						temp=wave.get_atom_basis_set_coefficient(a,p)*temp;
						basis_coefficients[a].push_back(temp);
						break;
					case 4:
						temp=32768*pow(temp,9);
						temp=temp/(225*pow(pi,3));
						temp=pow(temp,0.25);
						temp=wave.get_atom_basis_set_coefficient(a,p)*temp;
						basis_coefficients[a].push_back(temp);
						break;
					case -1:
						cout << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy..." << endl;
						Enter();
						break;
				}
			}
		for(int a=0; a< wave.get_ncen(); a++){
			double aiaj=0.0;
			double factor=0.0;
			for(int s=0; s<wave.get_atom_shell_count(a); s++){
				int type_temp=wave.get_shell_type(a,s);
				if(type_temp==-1){
					cout << "ERROR in type assignement!!"<<endl;
					Enter();
				}
				if (debug){
					cout << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << endl;
					int testcount=0;
					cout << "start: " << wave.get_shell_start(a,s,false) << flush;
					cout <<" stop: " << wave.get_shell_end(a,s,false) << flush << endl;
					cout << "factor: ";
				}
				switch (type_temp){
					case 1:
						factor=0;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							for(int j=wave.get_shell_start(a,s,false); j<= wave.get_shell_end(a,s,false); j++){
								double aiaj=wave.get_atom_basis_set_exponent (a,i)+wave.get_atom_basis_set_exponent (a,j);
								double term=(pi/aiaj);
								term=pow(term,1.5);
								factor+=basis_coefficients[a][i]*basis_coefficients[a][j]*term;
							}
						}
						if(factor==0) return false;
						factor=pow(factor,-0.5);
						if(debug) cout << factor << endl;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							if(debug){
								cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a,i) << endl;
								cout << "Contraction coefficient after:  " << factor*wave.get_atom_basis_set_coefficient(a,i) << endl;
							}
							basis_coefficients[a][i]*=factor;
							norm_const.push_back(basis_coefficients[a][i]);
						}
						break;
					case 2:
						factor=0;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							for(int j=wave.get_shell_start(a,s,false); j<= wave.get_shell_end(a,s,false); j++){
								double aiaj=wave.get_atom_basis_set_exponent (a,i)+wave.get_atom_basis_set_exponent (a,j);
								double term=4*pow(aiaj,5);
								term=pow(pi,3)/term;
								term=pow(term,0.5);
								factor+=basis_coefficients[a][i]*basis_coefficients[a][j]*term;
							}
						}
						if(factor==0) return false;
						factor=pow(factor,-0.5);
						if(debug) cout << factor << endl;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							if(debug){
								cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a,i) << endl;
								cout << "Contraction coefficient after:  " << factor*wave.get_atom_basis_set_coefficient(a,i) << endl;
							}
							basis_coefficients[a][i]*=factor;
							for(int k=0; k<3; k++) norm_const.push_back(basis_coefficients[a][i]);
						}
						break;
					case 3:
						factor=0;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							for(int j=wave.get_shell_start(a,s,false); j<= wave.get_shell_end(a,s,false); j++){
								double aiaj=wave.get_atom_basis_set_exponent (a,i)+wave.get_atom_basis_set_exponent (a,j);
								double term=16*pow(aiaj,7);
								term=pow(pi,3)/term;
								term=pow(term,0.5);
								factor+=basis_coefficients[a][i]*basis_coefficients[a][j]*term;
							}
						}
						if(factor==0) return false;
						factor=(pow(factor,-0.5))/sqrt(3);
						if(debug) cout << factor << endl;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							if(debug){
								cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a,i) << endl;
								cout << "Contraction coefficient after:  " << factor*wave.get_atom_basis_set_coefficient(a,i) << endl;
							}
							basis_coefficients[a][i]*=factor;
							for(int k=0; k<3; k++) norm_const.push_back(basis_coefficients[a][i]);
							for(int k=0; k<3; k++) norm_const.push_back(sqrt(3)*basis_coefficients[a][i]);
						}
						break;
					case 4:
						factor=0;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							for(int j=wave.get_shell_start(a,s,false); j<= wave.get_shell_end(a,s,false); j++){
								double aiaj=wave.get_atom_basis_set_exponent (a,i)+wave.get_atom_basis_set_exponent (a,j);
								double term=64*pow((aiaj),9);
								term=pow(pi,3)/term;
								term=pow(term,0.5);
								factor+=basis_coefficients[a][i]*basis_coefficients[a][j]*term;
							}
						}
						if(factor==0) return false;
						factor=pow(factor,-0.5)/sqrt(15);
						if(debug) cout << factor << endl;
						for(int i=wave.get_shell_start(a,s,false); i<= wave.get_shell_end(a,s,false); i++){
							if(debug){
								cout << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a,i) << endl;
								cout << "Contraction coefficient after:  " << factor*wave.get_atom_basis_set_coefficient(a,i) << endl;
							}
							basis_coefficients[a][i]*=factor;
							for(int l=0; l<3;l++) 	norm_const.push_back(basis_coefficients[a][i]);
							for(int l=0; l<6;l++) 	norm_const.push_back(sqrt(5)*basis_coefficients[a][i]);
													norm_const.push_back(sqrt(15)*basis_coefficients[a][i]);
						}
						break;
				}
				if (debug)	cout << "This shell has: " <<  wave.get_shell_end(a,s,false)-wave.get_shell_start(a,s,false)+1 << " primitives" << endl;
			}
		}
		//-----------debug output---------------------------------------------------------
		if(debug){
			cout << "exemplary output of the first atom with all it's properties: "<< endl;
			wave.print_atom_long(0);
			cout << "ended normalizing the basis set, now for the MO_coeffs" << endl;
			cout << "Status report:" << endl;
			cout << "size of norm_const: " << norm_const.size() << endl;
			cout << "WFN MO counter: " << wave.get_nmo() << endl;
			cout << "Number of atoms: " << wave.get_ncen() << endl;
			cout << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
			cout << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
			Enter();
			Enter();
		}
		//---------------------To not mix up anything start normalizing WFN_matrix now--------------------------
		int run=0;
		vector <vector <double> > changed_coefs;
		changed_coefs.resize(wave.get_nmo());
		ofstream norm_cprim;
		if(debug) norm_cprim.open("norm_prim.debug", ofstream::out);
		for(int m=0; m< wave.get_nmo(); m++){
			if(debug) norm_cprim << m << ". MO:" << endl;
			for(int p=0; p<wave.get_MO_primitive_count(m);p++){
				if(debug) cout << p << ". primitive; " << m << ". MO " << "norm nonst: " << norm_const[p] << endl;
				changed_coefs[m].push_back(wave.get_MO_coef(m,p, debug)/norm_const[p]);
				if(debug){
					cout << " temp after normalization: " << changed_coefs[m][p] << endl;
					norm_cprim << " " << changed_coefs[m][p] << endl;
				}
				run++;
			}
		}
		if(debug) {
			norm_cprim.flush();
			norm_cprim.close();
			cout << "See norm_cprim.debug for the CPRIM vectors" << endl;
			cout << "Total count in CPRIM: " << run << endl;
			Enter();
			Enter();
		}
		//--------------Build CMO of alessandro from the first elements of each shell-------------
		for (int m=0; m < wave.get_nmo (); m++){
			int run_2=0;
			for (int a =0; a < wave.get_ncen (); a++){
				for (int s=0; s< wave.get_atom_shell_count(a); s++){
					if(debug) cout << "Going to load the " << wave.get_shell_start_in_primitives(a,s) << ". value" << endl;
					switch (wave.get_shell_type(a,s)){
						case 1:
							CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a,s)]);
							if(m==0) nao++;
							if(debug&&wave.get_atom_shell_primitives(a,s)!=1)
								cout << "Pushing back 1 coefficient for S shell, this shell has " << wave.get_atom_shell_primitives(a,s) << " primitives! Shell start is: " << wave.get_shell_start(a,s,false) << endl;
						break;
						case 2:
							for (int i=0; i<3 ; i++) CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a,s)+i]);
							if(debug&&wave.get_atom_shell_primitives(a,s)!=1)
								cout << "Pushing back 3 coefficients for P shell, this shell has " << wave.get_atom_shell_primitives(a,s) << " primitives!" << endl;
							if(m==0) nao+=3;
						break;
						case 3:
							for (int i=0; i<6 ; i++) CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a,s)+i]);
							if(debug&&wave.get_atom_shell_primitives(a,s)!=1)
								cout << "Pushing back 6 coefficient for D shell, this shell has " << wave.get_atom_shell_primitives(a,s) << " primitives!" << endl;
							if(m==0) nao+=6;
						break;
						case 4:
							//this hardcoded piece is due to the order of f-type functions in the fchk
							for (int i=0; i<3 ; i++) 	CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a,s)+i]);
														CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a,s)+6]);
							for (int i=0; i<2 ; i++) 	CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a,s)+i+3]);
							for (int i=0; i<2 ; i++) 	CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a,s)+i+7]);
														CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a,s)+5]);
														CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a,s)+9]);
							if(debug&&wave.get_atom_shell_primitives(a,s)!=1)
								cout << "Pushing back 10 coefficient for F shell, this shell has " << wave.get_atom_shell_primitives(a,s) << " primitives!" << endl;
							if(m==0) nao+=10;
						break;
					}
					run_2++;
				}
				if(debug) cout << "finished with atom!" << endl;
			}
			if(debug) cout << "finished with MO!" << endl;
			if(nshell!=run_2) nshell=run_2;
		}
		if(debug){
			Enter();
			Enter();
			ofstream cmo ("cmo.debug", ofstream::out);
			for(int p=0; p< CMO.size() ; p++){
				string temp;
				for( int i=0; i< 5; i++){
					stringstream stream;
					stream << scientific << setw(14) << setprecision(7) << CMO[p+i] << " ";
					temp += stream.str();
				}
				p+=4;
				cmo << temp << endl;
			}
			cmo.flush();
			cmo.close();
			cout << CMO.size() << " Elements in CMO" << endl;
			cout << norm_const.size() << " = nprim" << endl;
			cout << nao << " = nao" << endl;
			cout << nshell << " = nshell" << endl;
			Enter();
		}
		//------------------ make the DM -----------------------------
		naotr=nao*(nao+1)/2;
		vector<double> kp;
		for(int i=0; i < naotr; i++) wave.push_back_DM (0.0);
		if(debug) cout << "I made kp!" << endl << nao << " is the maximum for iu" << endl;
		for(int iu=0; iu < nao; iu++){
			for(int iv=0; iv <= iu; iv++){
				int iuv=(iu*(iu+1)/2)+iv;
				if(debug) cout << "iu: " << iu << " iv: " << iv << " iuv: " << iuv << " kp(iu): " << iu*(iu+1)/2<< endl;
				for(int m=0; m < wave.get_nmo(); m++){
					if(debug) cout << "DM before: " << wave.get_DM(iuv) << endl;
					if(!wave.set_DM(iuv,wave.get_DM(iuv)+wave.get_MO_occ(m)*CMO[iu+(m*nao)]*CMO[iv+(m*nao)])){
						cout << "Something went wrong while writing the DM! iuv=" << iuv << endl;
						cout << "accessing CMO elements: " << iu+(m*nao) << " and " << iv+(m*nao) << " while CMO has size: " << CMO.size() << endl;
						cout << "CMO[" << iu+(m*nao) << "] = " << CMO[iu*(m*nao)] << "; CMO[" << iv+(m*nao) << "] = " << CMO[iv+(m*nao)] << " DM: " << wave.get_DM(iuv) << endl;
						cout << "DM_size: " << wave.get_DM_size() << endl;
						Enter();
						return false;
					}
					else if (debug) cout << "DM after: " << wave.get_DM(iuv) << endl;
				}
			}
		}
		if(debug){
			cout << "DM is:" << endl;
			for(int p=0; p< wave.get_DM_size() ; p++){
				string temp;
				for( int i=0; i< 5; i++){
					stringstream stream;
					stream << scientific << setw(14) << setprecision(7) << wave.get_DM(i+p) << " ";
					temp += stream.str();
				}
				p+=4;
				cout << temp << endl;
			}
			cout << wave.get_DM_size() << " Elements in DM" << endl;
			Enter();
		}
		// open fchk for writing
		string temp_fchk=fchk_name;
		temp_fchk.append("fchk");
		if(exists(temp_fchk) && !force_overwrite){
			cout << "The fchk already exists, do you want me to overwrite it?" << endl;
			if(!yesno()) return false;
		}
		ofstream fchk(temp_fchk.c_str());
		if(!fchk.is_open()){
			cout << "ERROR while opening .fchk ifile!" << endl;
			Enter();
			return false;
		}
		stringstream st_s;
		string s;
		st_s.str("");
		fchk << "TITLE\n";
		fchk.flush();

		fchk << "SP        RHF                                                         Gen\n";
		fchk.flush();
		s="Number of atoms                            I";
		st_s << setw(17) << wave.get_ncen();
		s+=st_s.str();
		st_s.str("");
		s+="\nInfo1-9                                    I   N=           9\n";
		s+="          53          51           0           0           0         111\n           1           1           2\n";  //I have NO CLUE what the fuck this means...
		fchk << s;
		fchk.flush();
		s="Charge                                     I";
		st_s << setw(17) << wave.get_charge();
		s+=st_s.str();
		st_s.str("");
		s+="\nMultiplicity                               I";
		st_s << setw(17) << wave.get_multi() << endl;
		s+=st_s.str();
		st_s.str("");
		int elcount=0;
		elcount-=wave.get_charge();
		for(int i=0; i< wave.get_ncen(); i++) elcount += wave.get_atom_charge(i);
		s+="Number of electrons                        I";
		st_s << setw(17) << elcount;
		s+=st_s.str();
		st_s.str("");
		s+="\nNumber of alpha electrons                  I";
		st_s << setw(17) << (int)ceil((float)elcount/2);
		s+=st_s.str();
		st_s.str("");
		s+="\nNumber of beta electrons                   I";
		st_s << setw(17) << (int)floor((float)elcount/2);
		s+=st_s.str();
		st_s.str("");
		s+="\nNumber of basis functions                  I";
		st_s << setw(17) << nao;
		s+=st_s.str();
		st_s.str("");
		s+="\nNumber of independent functions            I";
		st_s << setw(17) << nao << endl;
		s+=st_s.str();
		st_s.str("");
		s+="Number of point charges in /Mol/           I                0\nNumber of translation vectors              I                0\n";

		fchk << s;
		fchk.flush();

		s="Atomic numbers                             I   N=";
		st_s << setw(12) << wave.get_ncen() << endl;
		for(int i=0; i< wave.get_ncen(); i++){
			st_s << setw(12) << wave.get_atom_charge(i);
			if(((i+1)%6==0&&i!=0)||i==wave.get_ncen()-1) st_s << endl;
		}
		s+=st_s.str();
		st_s.str("");
		s+="Nuclear charges                            R   N=";
		st_s << setw(12) << wave.get_ncen() << endl;
		for(int i=0; i< wave.get_ncen(); i++){
			st_s << uppercase << scientific << setw(16) << setprecision(8) << (double) wave.get_atom_charge(i);
			if(((i+1)%5==0&&i!=0)||i==wave.get_ncen()-1) st_s << endl;
		}
		s+=st_s.str();
		st_s.str("");
		s+="Current cartesian coordinates              R   N=";
		st_s << setw(12) << wave.get_ncen()*3 << endl;
		s+=st_s.str();
		st_s.str("");
		unsigned int runs=0;
		for(int i=0; i< wave.get_ncen(); i++){
			for (int j=0; j<3; j++){
				st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_coordinate(i,j,debug);
				runs++;
				if(runs%5==0||(i==wave.get_ncen()-1&&j==2)) st_s << endl;
			}
		}
		s+=st_s.str();
		st_s.str("");
		s+="Integer atomic weights                     I   N=";
		st_s << setw(12) << wave.get_ncen() << endl;
		s+=st_s.str();
		st_s.str("");
		for(int i=0; i<wave.get_ncen(); i++){
			st_s << setw(12) << wave.get_atom_integer_mass((unsigned int) i);
			if(((i+1)%6==0&&i!=0)||i==wave.get_ncen()-1) st_s << endl;
		}
		s+=st_s.str();
		st_s.str("");
		s+="Real atomic weights                        R   N=";
		st_s << setw(12) << wave.get_ncen() << endl;
		s+=st_s.str();
		st_s.str("");
		for(int i=0; i<wave.get_ncen(); i++){
			st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_real_mass(i);
			if(((i+1)%5==0&&i!=0)||i==wave.get_ncen()-1) st_s << endl;
		}
		s+=st_s.str();
		st_s.str("");

		fchk << s;
		fchk.flush();

		unsigned int n_contracted_shells=0;
		for(int a=0; a<wave.get_ncen(); a++)
			for(int sh=0; sh<wave.get_atom_shell_count(a); sh++)
				n_contracted_shells++;
		s="Number of contracted shells                I";
		st_s << setw(17) << n_contracted_shells;
		s+=st_s.str();
		st_s.str("");
		unsigned int n_prim_shells=0;
		for(int a=0; a<wave.get_ncen(); a++)
			for(int sh=0; sh<wave.get_atom_shell_count(a); sh++)
				for(int p=0; p<wave.get_atom_shell_primitives(a,sh); p++)
					n_prim_shells++;
		s+="\nNumber of primitive shells                 I";
		st_s << setw(17) << n_prim_shells;
		s+=st_s.str();
		st_s.str("");
		s+="\nPure/Cartesian d shells                    I";
		int d_count=0;
		int f_count=0;
		int max_contraction=0;
		for(int a=0; a<wave.get_ncen(); a++)
			for(int s=0; s<wave.get_atom_shell_count(a); s++){
				if(wave.get_atom_shell_primitives(a,s) > max_contraction) max_contraction = wave.get_atom_shell_primitives(a,s);
				if(wave.get_shell_type(a,s)==3) d_count++;
				else if(wave.get_shell_type(a,s)==4) f_count++;
			}
		st_s << setw(17) << 1;
		s+=st_s.str();
		st_s.str("");
		s+="\nPure/Cartesian f shells                    I";
		st_s << setw(17) << 1;
		s+=st_s.str();
		st_s.str("");
		s+="\nHighest angular momentum                   I";
		st_s << setw(17) << 1+(d_count>0)+(f_count>0);
		s+=st_s.str();
		st_s.str("");
		s+="\nLargest degree of contraction              I";
		st_s << setw(17) << max_contraction;
		s+=st_s.str();
		st_s.str("");
		s+="\nShell types                                I   N=";
		st_s << setw(12) << n_contracted_shells << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for(int a=0; a<wave.get_ncen(); a++)
			for(int sh=0; sh<wave.get_atom_shell_count(a); sh++){
				st_s << setw(12) << wave.get_shell_type(a,sh)-1;
				runs++;
				if((runs%6==0&&runs!=0)||(a==wave.get_ncen()-1&&sh==wave.get_atom_shell_count(a)-1)) st_s << endl;
			}
		s+=st_s.str();
		st_s.str("");
		s+="Number of primitives per shell             I   N=";
		st_s << setw(12) << n_contracted_shells << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for(int a=0; a<wave.get_ncen(); a++)
			for(int sh=0; sh<wave.get_atom_shell_count(a); sh++){
				st_s << setw(12) << wave.get_atom_shell_primitives(a,sh);
				runs++;
				if((runs%6==0&&runs!=0)||(a==wave.get_ncen()-1&&sh==wave.get_atom_shell_count(a)-1)) st_s << endl;
			}
		s+=st_s.str();
		st_s.str("");
		s+="Shell to atom map                          I   N=";
		st_s << setw(12) << n_contracted_shells << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for(int a=0; a<wave.get_ncen(); a++)
			for(int sh=0; sh<wave.get_atom_shell_count(a); sh++){
				st_s << setw(12) << wave.get_shell_center(a,sh);
				runs++;
				if((runs%6==0&&runs!=0)||(a==wave.get_ncen()-1&&sh==wave.get_atom_shell_count(a)-1)) st_s << endl;
			}
		s+=st_s.str();
		st_s.str("");

		fchk << s;
		fchk.flush();

		s="Primitive exponents                        R   N=";
		st_s << setw(12) << n_prim_shells << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for(int a=0; a<wave.get_ncen(); a++){
			int p_run=0;
			for(int sh=0; sh<wave.get_atom_shell_count(a); sh++)
				for(int p=0; p<wave.get_atom_shell_primitives(a,sh); p++){
					st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_basis_set_exponent(a,p_run);
					runs++;
					p_run++;
					if((runs%5==0)||(a==wave.get_ncen()-1&&sh==wave.get_atom_shell_count(a)-1&&p==wave.get_atom_shell_primitives(a,sh)-1)) st_s << endl;
				}
		}
		s+=st_s.str();
		st_s.str("");
		s+="Contraction coefficients                   R   N=";
		st_s << setw(12) << n_prim_shells << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for(int a=0; a<wave.get_ncen(); a++){
			int p_run=0;
			for(int sh=0; sh<wave.get_atom_shell_count(a); sh++)
				for(int p=0; p<wave.get_atom_shell_primitives(a,sh); p++){
					st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_basis_set_coefficient(a,p_run);
					runs++;
					p_run++;
					if((runs%5==0&&runs!=0)||(a==wave.get_ncen()-1&&sh==wave.get_atom_shell_count(a)-1&&p==wave.get_atom_shell_primitives(a,sh)-1)) st_s << endl;
				}
		}
		s+=st_s.str();
		st_s.str("");
		s+="Coordinates of each shell                  R   N=";
		st_s << setw(12) << 3*n_contracted_shells << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for (int a=0; a < wave.get_ncen(); a++)
			for(int sh=0; sh<wave.get_atom_shell_count(a); sh++)
				for(int i=0; i<3; i++){
					st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_coordinate(a,i,false);
					runs++;
					if((runs%5==0&&runs!=0)||(a==wave.get_ncen()-1&&sh==wave.get_atom_shell_count(a)-1&&i==2)) st_s << endl;
				}
		s+=st_s.str();
		st_s.str("");
		s+="Constraint structure                       R   N=";
		st_s << setw(12) << 3*wave.get_ncen() << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for (int a=0; a < wave.get_ncen(); a++)
			for(int i=0; i<3; i++){
				st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_coordinate(a,i,false);
				runs++;
				if((runs%5==0&&runs!=0)||(a==wave.get_ncen()-1&&i==2)) st_s << endl;
			}
		s+=st_s.str();
		st_s.str("");
		s+="Virial Ratio                               R";
		st_s << uppercase << scientific << setw(27) << setprecision(15) << 2.0 << endl;
		s+=st_s.str();
		st_s.str("");
		s+="SCF Energy                                 R";
		double energy=0;
		for(int m=0; m<wave.get_nmo(); m++) if(wave.get_MO_occ(m)>0) energy+=wave.get_MO_energy(m);
		st_s << uppercase << scientific << setw(27) << setprecision(15) << energy << endl;
		s+=st_s.str();
		st_s.str("");
		s+="Total Energy                               R";
		st_s << uppercase << scientific << setw(27) << setprecision(15) << energy << endl;
		s+=st_s.str();
		st_s.str("");
		s+="RMS Density                                R";
		st_s << uppercase << scientific << setw(27) << setprecision(15) << 2.0*pow(10,-9) << endl;
		s+=st_s.str();
		st_s.str("");

		fchk << s;
		fchk.flush();

		s="Alpha Orbital Energies                     R   N=";
		st_s << setw(12) << nao << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for(int m=0; m<nao; m++){
			if(m<wave.get_nmo()) st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_MO_energy(m);
			else st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_MO_energy(wave.get_nmo()-1)+m;
			runs++;
			if((runs%5==0&&runs!=0)||m==nao-1) st_s << endl;
		}
		s+=st_s.str();
		st_s.str("");
		s+="Alpha MO coefficients                      R   N=";
		st_s << setw(12) << nao*nao << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for(int i=0; i< nao*nao; i++){
			if(i<CMO.size()) st_s << uppercase << scientific << setw(16) << setprecision(8) << CMO[i];
			else if(i<nao*nao) st_s << uppercase << scientific << setw(16) << setprecision(8) << 0.0;
			runs++;
			if((runs%5==0&&runs!=0)||i==nao*nao-1) st_s << endl;
		}
		s+=st_s.str();
		st_s.str("");
		s+="Total SCF Density                          R   N=";
		st_s << setw(12) << wave.get_DM_size() << endl;
		s+=st_s.str();
		st_s.str("");
		runs=0;
		for(int i=0; i< wave.get_DM_size(); i++){
			st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_DM(i);
			runs++;
			if((runs%5==0&&runs!=0)||i==wave.get_DM_size()-1) st_s << endl;
		}
		s+=st_s.str();
		st_s.str("");
		fchk << s;
		fchk.flush();
		fchk.close();
	}
	return true;
};

bool free_fchk(const string& fchk_name, const string& basis_set_path, WFN& wave, bool& debug, bool force_overwrite, ofstram &file) {
	if (wave.get_nr_basis_set_loaded() == 0) {
		file << "No basis set loaded, will load a complete basis set now!" << endl;
		if (debug) Enter();
		if (!read_basis_set_vanilla(basis_set_path, wave, debug, false)) {
			file << "Problem during reading of the basis set!" << endl;
			return "WRONG";
		}
	}
	else if (wave.get_nr_basis_set_loaded() < wave.get_ncen()) {
		file << "Not all atoms have a basis set loaded!" << endl
			<< "Do you want to laod the missing atoms?" << flush;
		if (!yesno()) {
			file << "Do you want to load a new basis set?" << flush;
			if (!yesno()) {
				file << "Okay, aborting then!!!" << endl;
				if (debug) Enter();
				return "WRONG";
			}
			file << "deleting old one..." << endl;
			if (debug) Enter();
			if (!delete_basis_set_vanilla(basis_set_path, wave, debug)) {
				file << "ERROR while deleting a basis set!" << endl;
				return "WRONG";
			}
			else if (!read_basis_set_vanilla(basis_set_path, wave, debug, true)) {
				file << "ERROR during reading of the new basis set!" << endl;
				return "WRONG";
			}
		}
		else {
			file << "okay, loading the missing atoms.." << endl;
			if (debug) Enter();
			if (!read_basis_set_missing(basis_set_path, wave, debug)) {
				file << "ERROR during reading of missing basis set!" << endl;
				return "WRONG";
			}
		}
	}
	else if (wave.get_nr_basis_set_loaded() == wave.get_ncen()) {
		file << "There already is a basis set loaded!" << endl
			<< "Do you want me to delete the old one and load a new one?" << flush;
		if (debug) Enter();
		if (!yesno()) {
			file << "Okay, continuing with the old one..." << endl;
			if (debug) Enter();
		}
		else {
			file << "Deleting the old basis set!" << endl;
			if (debug) Enter();
			if (!delete_basis_set_vanilla(basis_set_path, wave, debug)) {
				file << "ERROR during deleting of the basis set!";
				Enter();
				return "WRONG";
			}
			file << "Going to load a new one now!" << endl;
			if (!read_basis_set_vanilla(basis_set_path, wave, debug, true)) {
				file << "Problem during reading of the basis set!" << endl;
				return "WRONG";
			}
		}
	}
	else {
		file << "# of loaded > # atoms" << endl
			<< "Sorry, this should not happen... aborting!!!" << endl;
		return "WRONG";
	}
	//wave.set_modified();
	if (debug) debug_fchk = true;
	vector<double> CMO;
	int nao = 0;
	int nshell = 0;
	int naotr = 0;
	if (debug) {
		file << "Origin: " << wave.get_origin() << endl;
		Enter();
	}
	if (wave.get_origin() == 2 || wave.get_origin() == 4) {
		double pi = 3.14159265358979;
		//---------------normalize basis set---------------------------------
		if (debug) file << "starting to normalize the basis set" << endl;
		vector<double> norm_const;
		//-----------debug output---------------------------------------------------------
		if (debug) {
			file << "exemplary output before norm_const of the first atom with all it's properties: " << endl;
			wave.print_atom_long(0);
			file << "ended normalizing the basis set, now for the MO_coeffs" << endl;
			file << "Status report:" << endl;
			file << "size of norm_const: " << norm_const.size() << endl;
			file << "WFN MO counter: " << wave.get_nmo() << endl;
			file << "Number of atoms: " << wave.get_ncen() << endl;
			file << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
			file << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
			Enter();
		}
		//-----------------------check ordering and order accordingly----------------------
		wave.sort_wfn(wave.check_order(debug), debug);
		//-------------------normalize the basis set shell wise into a copy vector---------
		vector <vector <double> > basis_coefficients;
		basis_coefficients.resize(wave.get_ncen());
		for (int a = 0; a < wave.get_ncen(); a++)
			for (int p = 0; p < wave.get_atom_primitive_count(a); p++) {
				double temp = wave.get_atom_basis_set_exponent(a, p);
				switch (wave.get_atom_primitive_type(a, p)) {
				case 1:
					temp = 2 * temp / pi;
					temp = pow(temp, 0.75);
					temp = temp * wave.get_atom_basis_set_coefficient(a, p);
					basis_coefficients[a].push_back(temp);
					break;
				case 2:
					temp = 128 * pow(temp, 5);
					temp = temp / pow(pi, 3);
					temp = pow(temp, 0.25);
					temp = wave.get_atom_basis_set_coefficient(a, p) * temp;
					basis_coefficients[a].push_back(temp);
					break;
				case 3:
					temp = 2048 * pow(temp, 7);
					temp = temp / (9 * pow(pi, 3));
					temp = pow(temp, 0.25);
					temp = wave.get_atom_basis_set_coefficient(a, p) * temp;
					basis_coefficients[a].push_back(temp);
					break;
				case 4:
					temp = 32768 * pow(temp, 9);
					temp = temp / (225 * pow(pi, 3));
					temp = pow(temp, 0.25);
					temp = wave.get_atom_basis_set_coefficient(a, p) * temp;
					basis_coefficients[a].push_back(temp);
					break;
				case -1:
					file << "Sorry, the type reading went wrong somwhere, look where it may have gone crazy..." << endl;
					Enter();
					break;
				}
			}
		for (int a = 0; a < wave.get_ncen(); a++) {
			double aiaj = 0.0;
			double factor = 0.0;
			for (int s = 0; s < wave.get_atom_shell_count(a); s++) {
				int type_temp = wave.get_shell_type(a, s);
				if (type_temp == -1) {
					file << "ERROR in type assignement!!" << endl;
					Enter();
				}
				if (debug) {
					file << "Shell: " << s << " of atom: " << a << " Shell type: " << type_temp << endl;
					int testcount = 0;
					file << "start: " << wave.get_shell_start(a, s, false) << flush;
					file << " stop: " << wave.get_shell_end(a, s, false) << flush << endl;
					file << "factor: ";
				}
				switch (type_temp) {
				case 1:
					factor = 0;
					for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
						for (int j = wave.get_shell_start(a, s, false); j <= wave.get_shell_end(a, s, false); j++) {
							double aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
							double term = (pi / aiaj);
							term = pow(term, 1.5);
							factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
						}
					}
					if (factor == 0) return false;
					factor = pow(factor, -0.5);
					if (debug) file << factor << endl;
					for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
						if (debug) {
							file << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i) << endl;
							file << "Contraction coefficient after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
						}
						basis_coefficients[a][i] *= factor;
						norm_const.push_back(basis_coefficients[a][i]);
					}
					break;
				case 2:
					factor = 0;
					for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
						for (int j = wave.get_shell_start(a, s, false); j <= wave.get_shell_end(a, s, false); j++) {
							double aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
							double term = 4 * pow(aiaj, 5);
							term = pow(pi, 3) / term;
							term = pow(term, 0.5);
							factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
						}
					}
					if (factor == 0) return false;
					factor = pow(factor, -0.5);
					if (debug) file << factor << endl;
					for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
						if (debug) {
							file << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i) << endl;
							file << "Contraction coefficient after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
						}
						basis_coefficients[a][i] *= factor;
						for (int k = 0; k < 3; k++) norm_const.push_back(basis_coefficients[a][i]);
					}
					break;
				case 3:
					factor = 0;
					for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
						for (int j = wave.get_shell_start(a, s, false); j <= wave.get_shell_end(a, s, false); j++) {
							double aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
							double term = 16 * pow(aiaj, 7);
							term = pow(pi, 3) / term;
							term = pow(term, 0.5);
							factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
						}
					}
					if (factor == 0) return false;
					factor = (pow(factor, -0.5)) / sqrt(3);
					if (debug) file << factor << endl;
					for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
						if (debug) {
							file << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i) << endl;
							file << "Contraction coefficient after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
						}
						basis_coefficients[a][i] *= factor;
						for (int k = 0; k < 3; k++) norm_const.push_back(basis_coefficients[a][i]);
						for (int k = 0; k < 3; k++) norm_const.push_back(sqrt(3) * basis_coefficients[a][i]);
					}
					break;
				case 4:
					factor = 0;
					for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
						for (int j = wave.get_shell_start(a, s, false); j <= wave.get_shell_end(a, s, false); j++) {
							double aiaj = wave.get_atom_basis_set_exponent(a, i) + wave.get_atom_basis_set_exponent(a, j);
							double term = 64 * pow((aiaj), 9);
							term = pow(pi, 3) / term;
							term = pow(term, 0.5);
							factor += basis_coefficients[a][i] * basis_coefficients[a][j] * term;
						}
					}
					if (factor == 0) return false;
					factor = pow(factor, -0.5) / sqrt(15);
					if (debug) file << factor << endl;
					for (int i = wave.get_shell_start(a, s, false); i <= wave.get_shell_end(a, s, false); i++) {
						if (debug) {
							file << "Contraction coefficient before: " << wave.get_atom_basis_set_coefficient(a, i) << endl;
							file << "Contraction coefficient after:  " << factor * wave.get_atom_basis_set_coefficient(a, i) << endl;
						}
						basis_coefficients[a][i] *= factor;
						for (int l = 0; l < 3; l++) 	norm_const.push_back(basis_coefficients[a][i]);
						for (int l = 0; l < 6; l++) 	norm_const.push_back(sqrt(5) * basis_coefficients[a][i]);
						norm_const.push_back(sqrt(15) * basis_coefficients[a][i]);
					}
					break;
				}
				if (debug)	file << "This shell has: " << wave.get_shell_end(a, s, false) - wave.get_shell_start(a, s, false) + 1 << " primitives" << endl;
			}
		}
		//-----------debug output---------------------------------------------------------
		if (debug) {
			file << "exemplary output of the first atom with all it's properties: " << endl;
			wave.print_atom_long(0);
			file << "ended normalizing the basis set, now for the MO_coeffs" << endl;
			file << "Status report:" << endl;
			file << "size of norm_const: " << norm_const.size() << endl;
			file << "WFN MO counter: " << wave.get_nmo() << endl;
			file << "Number of atoms: " << wave.get_ncen() << endl;
			file << "Primitive count of zero MO: " << wave.get_MO_primitive_count(0) << endl;
			file << "Primitive count of first MO: " << wave.get_MO_primitive_count(1) << endl;
			Enter();
			Enter();
		}
		//---------------------To not mix up anything start normalizing WFN_matrix now--------------------------
		int run = 0;
		vector <vector <double> > changed_coefs;
		changed_coefs.resize(wave.get_nmo());
		ofstream norm_cprim;
		if (debug) norm_cprim.open("norm_prim.debug", ofstream::out);
		for (int m = 0; m < wave.get_nmo(); m++) {
			if (debug) norm_cprim << m << ". MO:" << endl;
			for (int p = 0; p < wave.get_MO_primitive_count(m); p++) {
				if (debug) file << p << ". primitive; " << m << ". MO " << "norm nonst: " << norm_const[p] << endl;
				changed_coefs[m].push_back(wave.get_MO_coef(m, p, debug) / norm_const[p]);
				if (debug) {
					file << " temp after normalization: " << changed_coefs[m][p] << endl;
					norm_cprim << " " << changed_coefs[m][p] << endl;
				}
				run++;
			}
		}
		if (debug) {
			norm_cprim.flush();
			norm_cprim.close();
			file << "See norm_cprim.debug for the CPRIM vectors" << endl;
			file << "Total count in CPRIM: " << run << endl;
			Enter();
			Enter();
		}
		//--------------Build CMO of alessandro from the first elements of each shell-------------
		for (int m = 0; m < wave.get_nmo(); m++) {
			int run_2 = 0;
			for (int a = 0; a < wave.get_ncen(); a++) {
				for (int s = 0; s < wave.get_atom_shell_count(a); s++) {
					if (debug) file << "Going to load the " << wave.get_shell_start_in_primitives(a, s) << ". value" << endl;
					switch (wave.get_shell_type(a, s)) {
					case 1:
						CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s)]);
						if (m == 0) nao++;
						if (debug && wave.get_atom_shell_primitives(a, s) != 1)
							file << "Pushing back 1 coefficient for S shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives! Shell start is: " << wave.get_shell_start(a, s, false) << endl;
						break;
					case 2:
						for (int i = 0; i < 3; i++) CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i]);
						if (debug && wave.get_atom_shell_primitives(a, s) != 1)
							file << "Pushing back 3 coefficients for P shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
						if (m == 0) nao += 3;
						break;
					case 3:
						for (int i = 0; i < 6; i++) CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i]);
						if (debug && wave.get_atom_shell_primitives(a, s) != 1)
							file << "Pushing back 6 coefficient for D shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
						if (m == 0) nao += 6;
						break;
					case 4:
						//this hardcoded piece is due to the order of f-type functions in the fchk
						for (int i = 0; i < 3; i++) 	CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i]);
						CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + 6]);
						for (int i = 0; i < 2; i++) 	CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i + 3]);
						for (int i = 0; i < 2; i++) 	CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + i + 7]);
						CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + 5]);
						CMO.push_back(changed_coefs[m][wave.get_shell_start_in_primitives(a, s) + 9]);
						if (debug && wave.get_atom_shell_primitives(a, s) != 1)
							file << "Pushing back 10 coefficient for F shell, this shell has " << wave.get_atom_shell_primitives(a, s) << " primitives!" << endl;
						if (m == 0) nao += 10;
						break;
					}
					run_2++;
				}
				if (debug) file << "finished with atom!" << endl;
			}
			if (debug) file << "finished with MO!" << endl;
			if (nshell != run_2) nshell = run_2;
		}
		if (debug) {
			Enter();
			Enter();
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
			file << CMO.size() << " Elements in CMO" << endl;
			file << norm_const.size() << " = nprim" << endl;
			file << nao << " = nao" << endl;
			file << nshell << " = nshell" << endl;
			Enter();
		}
		//------------------ make the DM -----------------------------
		naotr = nao * (nao + 1) / 2;
		vector<double> kp;
		for (int i = 0; i < naotr; i++) wave.push_back_DM(0.0);
		if (debug) file << "I made kp!" << endl << nao << " is the maximum for iu" << endl;
		for (int iu = 0; iu < nao; iu++) {
			for (int iv = 0; iv <= iu; iv++) {
				int iuv = (iu * (iu + 1) / 2) + iv;
				if (debug) file << "iu: " << iu << " iv: " << iv << " iuv: " << iuv << " kp(iu): " << iu * (iu + 1) / 2 << endl;
				for (int m = 0; m < wave.get_nmo(); m++) {
					if (debug) file << "DM before: " << wave.get_DM(iuv) << endl;
					if (!wave.set_DM(iuv, wave.get_DM(iuv) + wave.get_MO_occ(m) * CMO[iu + (m * nao)] * CMO[iv + (m * nao)])) {
						file << "Something went wrong while writing the DM! iuv=" << iuv << endl;
						file << "accessing CMO elements: " << iu + (m * nao) << " and " << iv + (m * nao) << " while CMO has size: " << CMO.size() << endl;
						file << "CMO[" << iu + (m * nao) << "] = " << CMO[iu * (m * nao)] << "; CMO[" << iv + (m * nao) << "] = " << CMO[iv + (m * nao)] << " DM: " << wave.get_DM(iuv) << endl;
						file << "DM_size: " << wave.get_DM_size() << endl;
						Enter();
						return false;
					}
					else if (debug) file << "DM after: " << wave.get_DM(iuv) << endl;
				}
			}
		}
		if (debug) {
			file << "DM is:" << endl;
			for (int p = 0; p < wave.get_DM_size(); p++) {
				string temp;
				for (int i = 0; i < 5; i++) {
					stringstream stream;
					stream << scientific << setw(14) << setprecision(7) << wave.get_DM(i + p) << " ";
					temp += stream.str();
				}
				p += 4;
				file << temp << endl;
			}
			file << wave.get_DM_size() << " Elements in DM" << endl;
			Enter();
		}
		// open fchk for writing
		string temp_fchk = fchk_name;
		temp_fchk.append("fchk");
		if (exists(temp_fchk) && !force_overwrite) {
			file << "The fchk already exists, do you want me to overwrite it?" << endl;
			if (!yesno()) return false;
		}
		ofstream fchk(temp_fchk.c_str());
		if (!fchk.is_open()) {
			file << "ERROR while opening .fchk ifile!" << endl;
			Enter();
			return false;
		}
		stringstream st_s;
		string s;
		st_s.str("");
		fchk << "TITLE\n";
		fchk.flush();

		fchk << "SP        RHF                                                         Gen\n";
		fchk.flush();
		s = "Number of atoms                            I";
		st_s << setw(17) << wave.get_ncen();
		s += st_s.str();
		st_s.str("");
		s += "\nInfo1-9                                    I   N=           9\n";
		s += "          53          51           0           0           0         111\n           1           1           2\n";  //I have NO CLUE what the fuck this means...
		fchk << s;
		fchk.flush();
		s = "Charge                                     I";
		st_s << setw(17) << wave.get_charge();
		s += st_s.str();
		st_s.str("");
		s += "\nMultiplicity                               I";
		st_s << setw(17) << wave.get_multi() << endl;
		s += st_s.str();
		st_s.str("");
		int elcount = 0;
		elcount -= wave.get_charge();
		for (int i = 0; i < wave.get_ncen(); i++) elcount += wave.get_atom_charge(i);
		s += "Number of electrons                        I";
		st_s << setw(17) << elcount;
		s += st_s.str();
		st_s.str("");
		s += "\nNumber of alpha electrons                  I";
		st_s << setw(17) << (int)ceil((float)elcount / 2);
		s += st_s.str();
		st_s.str("");
		s += "\nNumber of beta electrons                   I";
		st_s << setw(17) << (int)floor((float)elcount / 2);
		s += st_s.str();
		st_s.str("");
		s += "\nNumber of basis functions                  I";
		st_s << setw(17) << nao;
		s += st_s.str();
		st_s.str("");
		s += "\nNumber of independent functions            I";
		st_s << setw(17) << nao << endl;
		s += st_s.str();
		st_s.str("");
		s += "Number of point charges in /Mol/           I                0\nNumber of translation vectors              I                0\n";

		fchk << s;
		fchk.flush();

		s = "Atomic numbers                             I   N=";
		st_s << setw(12) << wave.get_ncen() << endl;
		for (int i = 0; i < wave.get_ncen(); i++) {
			st_s << setw(12) << wave.get_atom_charge(i);
			if (((i + 1) % 6 == 0 && i != 0) || i == wave.get_ncen() - 1) st_s << endl;
		}
		s += st_s.str();
		st_s.str("");
		s += "Nuclear charges                            R   N=";
		st_s << setw(12) << wave.get_ncen() << endl;
		for (int i = 0; i < wave.get_ncen(); i++) {
			st_s << uppercase << scientific << setw(16) << setprecision(8) << (double)wave.get_atom_charge(i);
			if (((i + 1) % 5 == 0 && i != 0) || i == wave.get_ncen() - 1) st_s << endl;
		}
		s += st_s.str();
		st_s.str("");
		s += "Current cartesian coordinates              R   N=";
		st_s << setw(12) << wave.get_ncen() * 3 << endl;
		s += st_s.str();
		st_s.str("");
		unsigned int runs = 0;
		for (int i = 0; i < wave.get_ncen(); i++) {
			for (int j = 0; j < 3; j++) {
				st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_coordinate(i, j, debug);
				runs++;
				if (runs % 5 == 0 || (i == wave.get_ncen() - 1 && j == 2)) st_s << endl;
			}
		}
		s += st_s.str();
		st_s.str("");
		s += "Integer atomic weights                     I   N=";
		st_s << setw(12) << wave.get_ncen() << endl;
		s += st_s.str();
		st_s.str("");
		for (int i = 0; i < wave.get_ncen(); i++) {
			st_s << setw(12) << wave.get_atom_integer_mass((unsigned int)i);
			if (((i + 1) % 6 == 0 && i != 0) || i == wave.get_ncen() - 1) st_s << endl;
		}
		s += st_s.str();
		st_s.str("");
		s += "Real atomic weights                        R   N=";
		st_s << setw(12) << wave.get_ncen() << endl;
		s += st_s.str();
		st_s.str("");
		for (int i = 0; i < wave.get_ncen(); i++) {
			st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_real_mass(i);
			if (((i + 1) % 5 == 0 && i != 0) || i == wave.get_ncen() - 1) st_s << endl;
		}
		s += st_s.str();
		st_s.str("");

		fchk << s;
		fchk.flush();

		unsigned int n_contracted_shells = 0;
		for (int a = 0; a < wave.get_ncen(); a++)
			for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
				n_contracted_shells++;
		s = "Number of contracted shells                I";
		st_s << setw(17) << n_contracted_shells;
		s += st_s.str();
		st_s.str("");
		unsigned int n_prim_shells = 0;
		for (int a = 0; a < wave.get_ncen(); a++)
			for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
				for (int p = 0; p < wave.get_atom_shell_primitives(a, sh); p++)
					n_prim_shells++;
		s += "\nNumber of primitive shells                 I";
		st_s << setw(17) << n_prim_shells;
		s += st_s.str();
		st_s.str("");
		s += "\nPure/Cartesian d shells                    I";
		int d_count = 0;
		int f_count = 0;
		int max_contraction = 0;
		for (int a = 0; a < wave.get_ncen(); a++)
			for (int s = 0; s < wave.get_atom_shell_count(a); s++) {
				if (wave.get_atom_shell_primitives(a, s) > max_contraction) max_contraction = wave.get_atom_shell_primitives(a, s);
				if (wave.get_shell_type(a, s) == 3) d_count++;
				else if (wave.get_shell_type(a, s) == 4) f_count++;
			}
		st_s << setw(17) << 1;
		s += st_s.str();
		st_s.str("");
		s += "\nPure/Cartesian f shells                    I";
		st_s << setw(17) << 1;
		s += st_s.str();
		st_s.str("");
		s += "\nHighest angular momentum                   I";
		st_s << setw(17) << 1 + (d_count > 0) + (f_count > 0);
		s += st_s.str();
		st_s.str("");
		s += "\nLargest degree of contraction              I";
		st_s << setw(17) << max_contraction;
		s += st_s.str();
		st_s.str("");
		s += "\nShell types                                I   N=";
		st_s << setw(12) << n_contracted_shells << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int a = 0; a < wave.get_ncen(); a++)
			for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++) {
				st_s << setw(12) << wave.get_shell_type(a, sh) - 1;
				runs++;
				if ((runs % 6 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1)) st_s << endl;
			}
		s += st_s.str();
		st_s.str("");
		s += "Number of primitives per shell             I   N=";
		st_s << setw(12) << n_contracted_shells << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int a = 0; a < wave.get_ncen(); a++)
			for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++) {
				st_s << setw(12) << wave.get_atom_shell_primitives(a, sh);
				runs++;
				if ((runs % 6 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1)) st_s << endl;
			}
		s += st_s.str();
		st_s.str("");
		s += "Shell to atom map                          I   N=";
		st_s << setw(12) << n_contracted_shells << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int a = 0; a < wave.get_ncen(); a++)
			for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++) {
				st_s << setw(12) << wave.get_shell_center(a, sh);
				runs++;
				if ((runs % 6 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1)) st_s << endl;
			}
		s += st_s.str();
		st_s.str("");

		fchk << s;
		fchk.flush();

		s = "Primitive exponents                        R   N=";
		st_s << setw(12) << n_prim_shells << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int a = 0; a < wave.get_ncen(); a++) {
			int p_run = 0;
			for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
				for (int p = 0; p < wave.get_atom_shell_primitives(a, sh); p++) {
					st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_basis_set_exponent(a, p_run);
					runs++;
					p_run++;
					if ((runs % 5 == 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1 && p == wave.get_atom_shell_primitives(a, sh) - 1)) st_s << endl;
				}
		}
		s += st_s.str();
		st_s.str("");
		s += "Contraction coefficients                   R   N=";
		st_s << setw(12) << n_prim_shells << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int a = 0; a < wave.get_ncen(); a++) {
			int p_run = 0;
			for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
				for (int p = 0; p < wave.get_atom_shell_primitives(a, sh); p++) {
					st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_basis_set_coefficient(a, p_run);
					runs++;
					p_run++;
					if ((runs % 5 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1 && p == wave.get_atom_shell_primitives(a, sh) - 1)) st_s << endl;
				}
		}
		s += st_s.str();
		st_s.str("");
		s += "Coordinates of each shell                  R   N=";
		st_s << setw(12) << 3 * n_contracted_shells << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int a = 0; a < wave.get_ncen(); a++)
			for (int sh = 0; sh < wave.get_atom_shell_count(a); sh++)
				for (int i = 0; i < 3; i++) {
					st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_coordinate(a, i, false);
					runs++;
					if ((runs % 5 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && sh == wave.get_atom_shell_count(a) - 1 && i == 2)) st_s << endl;
				}
		s += st_s.str();
		st_s.str("");
		s += "Constraint structure                       R   N=";
		st_s << setw(12) << 3 * wave.get_ncen() << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int a = 0; a < wave.get_ncen(); a++)
			for (int i = 0; i < 3; i++) {
				st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_atom_coordinate(a, i, false);
				runs++;
				if ((runs % 5 == 0 && runs != 0) || (a == wave.get_ncen() - 1 && i == 2)) st_s << endl;
			}
		s += st_s.str();
		st_s.str("");
		s += "Virial Ratio                               R";
		st_s << uppercase << scientific << setw(27) << setprecision(15) << 2.0 << endl;
		s += st_s.str();
		st_s.str("");
		s += "SCF Energy                                 R";
		double energy = 0;
		for (int m = 0; m < wave.get_nmo(); m++) if (wave.get_MO_occ(m) > 0) energy += wave.get_MO_energy(m);
		st_s << uppercase << scientific << setw(27) << setprecision(15) << energy << endl;
		s += st_s.str();
		st_s.str("");
		s += "Total Energy                               R";
		st_s << uppercase << scientific << setw(27) << setprecision(15) << energy << endl;
		s += st_s.str();
		st_s.str("");
		s += "RMS Density                                R";
		st_s << uppercase << scientific << setw(27) << setprecision(15) << 2.0 * pow(10, -9) << endl;
		s += st_s.str();
		st_s.str("");

		fchk << s;
		fchk.flush();

		s = "Alpha Orbital Energies                     R   N=";
		st_s << setw(12) << nao << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int m = 0; m < nao; m++) {
			if (m < wave.get_nmo()) st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_MO_energy(m);
			else st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_MO_energy(wave.get_nmo() - 1) + m;
			runs++;
			if ((runs % 5 == 0 && runs != 0) || m == nao - 1) st_s << endl;
		}
		s += st_s.str();
		st_s.str("");
		s += "Alpha MO coefficients                      R   N=";
		st_s << setw(12) << nao * nao << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int i = 0; i < nao * nao; i++) {
			if (i < CMO.size()) st_s << uppercase << scientific << setw(16) << setprecision(8) << CMO[i];
			else if (i < nao * nao) st_s << uppercase << scientific << setw(16) << setprecision(8) << 0.0;
			runs++;
			if ((runs % 5 == 0 && runs != 0) || i == nao * nao - 1) st_s << endl;
		}
		s += st_s.str();
		st_s.str("");
		s += "Total SCF Density                          R   N=";
		st_s << setw(12) << wave.get_DM_size() << endl;
		s += st_s.str();
		st_s.str("");
		runs = 0;
		for (int i = 0; i < wave.get_DM_size(); i++) {
			st_s << uppercase << scientific << setw(16) << setprecision(8) << wave.get_DM(i);
			runs++;
			if ((runs % 5 == 0 && runs != 0) || i == wave.get_DM_size() - 1) st_s << endl;
		}
		s += st_s.str();
		st_s.str("");
		fchk << s;
		fchk.flush();
		fchk.close();
	}
	return true;
};
