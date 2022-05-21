#pragma once

#include <iostream>
#include <vector>
#include <string>
#include "convenience.h"

//-----------------Definition of atoms and basis sets--------------------

struct basis_set_entry{
	double coefficient;
	double exponent;
	unsigned int type; //1=S; 2=P; 3=D; 4=F; 5=G
	unsigned int shell;
	basis_set_entry operator=(const basis_set_entry &rhs);
	basis_set_entry();
	basis_set_entry(double g_coefficient, double g_exponent, unsigned int g_type, unsigned int g_shell);
};

inline basis_set_entry::basis_set_entry(){
	coefficient=0.0;
	exponent=0.0;
	type=0;
	shell=0;
};

inline basis_set_entry::basis_set_entry(double g_coefficient, double g_exponent, unsigned int g_type, unsigned int g_shell){
	coefficient=g_coefficient;
	exponent=g_exponent;
	type=g_type;
	shell=g_shell;
};

inline basis_set_entry basis_set_entry::operator=(const basis_set_entry &rhs){
	coefficient = rhs.coefficient;
	exponent = rhs.exponent;
	type = rhs.type;
	return *this;
};

struct atom {
	std::string label;
	int nr, charge, ECP_electrons;
	double x, y, z;
	std::vector<double> frac_coords;
	atom();
	atom(const std::string &l, const int &n, const double &c1, const double &c2, const double &c3, const int &ch);
	atom(const std::string& l, const int& n, const double& c1, const double& c2, const double& c3, const int& ch, const int& ECP_els);
	atom operator=(const atom &rhs);
	void print_values();
	bool push_back_basis_set(const double &coefficient, const double &exponent, const int &type, const int &shell);
	void print_values_long();
	bool get_basis_set_loaded();
	bool is_anharm();
	void assign_ADPs(std::vector<double> &second, std::vector<double> &third, std::vector<double> &fourth);
	void assign_ADPs(std::vector<double>& second);
	void assign_ADPs(double& Uiso);
	std::vector<basis_set_entry> basis_set;
	std::vector<unsigned int> shellcount;
	//The Order is:
	//[0] = second order (U11, U22, U33, U12, U13, U23)
	//[1] = third order  (C111, C112, C113, C122, C123, C133, C222, C223, C233, C333)
	//[2] = fourth order (D1111, D1112, D1113, D1122, D1123, D1133, D1222, D1223, D1233, D1333, D2222, D2223, D2233, D2333, D3333)
	std::vector<std::vector<double>> ADPs;
};						  

inline atom::atom() {
	label = '?';
	nr = 0;
	x = 0.0;
	y = 0.0;
	z = 0.0;
	charge = 0;
	ECP_electrons = 0;
	frac_coords.resize(3);
	frac_coords = { 0,0,0 };
};

inline atom::atom (const std::string &l, const int &n, const double &c1, const double &c2, const double &c3, const int &ch){
	nr = n;
	label = l;
	x = c1;
	y = c2;
	z = c3;
	charge = ch;
	ECP_electrons = 0;
	frac_coords.resize(3);
	frac_coords = { 0,0,0 };
};

inline atom::atom(const std::string& l, const int& n, const double& c1, const double& c2, const double& c3, const int& ch, const int& ECP_els) {
	nr = n;
	label = l;
	x = c1;
	y = c2;
	z = c3;
	charge = ch;
	ECP_electrons = ECP_els;
	frac_coords.resize(3);
	frac_coords = { 0,0,0 };
};

inline atom atom::operator= (const atom &rhs){
	label = rhs.label;
	nr = rhs.nr;
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	charge = rhs.charge;
	basis_set = rhs.basis_set;
	shellcount = rhs.shellcount;
	ECP_electrons = rhs.ECP_electrons;
	frac_coords = rhs.frac_coords;
	return *this;
};

inline void atom::print_values(){
	std::cout << "nr: " << nr << " label: " << label << " x: " << x << " y: " << y << " z: " << z << " charge: " << charge << std::endl;
};

inline bool atom::is_anharm() {
	return ADPs.size() > 1;
};

inline void atom::print_values_long(){
	std::cout << "nr: " << nr << " label: " << label << " x: " << x << " y: " << y << " z: " << z << " charge: " << charge << std::endl;
	if(basis_set.size()>0){
		for (unsigned int i=0; i<basis_set.size(); i++){
			std::cout << "basis set entry " << i << ": expo: " << basis_set[i].exponent << " coef: " << basis_set[i].coefficient << " type: " << basis_set[i].type << " shell: " << basis_set[i].shell << std::endl;
		}
	}
	if(shellcount.size()>0){
		for (unsigned int i=0; i<shellcount.size(); i++){
			std::cout << "shellcount " << i << ": value: " << shellcount[i] << std::endl;
		}
	}
};

inline bool atom::push_back_basis_set(const double &exponent, const double &coefficient, const int &type, const int &shell){
	if(shell == shellcount.size())
		shellcount.push_back(1);
	else
		shellcount[shell]++;
	if( type >= 0 && type < 6 && shell>=0){
		basis_set.push_back(basis_set_entry(coefficient, exponent, type, shell));
		return true;
	}
	else{
		if (type >= 6) err_checkc(false, "h and higher types are not yet supported!");
		std::cout << "This is not a valid basis set entry!" << std::endl;
		std::cout << "Exponent: " << exponent << " coefficient: "<< coefficient << " type: " << type << " shell: " << shell << std::endl;
		return false;
	}
};

inline bool atom::get_basis_set_loaded(){
	if(basis_set.size()>0) return true;
	else return false;
};

inline void atom::assign_ADPs(double& Uiso) {
	ADPs.resize(1);
	ADPs[0].resize(6);
	ADPs[0][0] = ADPs[0][1] = ADPs[0][2] = Uiso;
};

inline void atom::assign_ADPs(std::vector<double>& second) {
	if (second.size() != 6) {
		std::cout << "Wrong size of second order ADP!" << std::endl;
		return;
	}
	else {
		ADPs.resize(1);
		ADPs[0].resize(6);
		ADPs[0] = second;
	}
};

inline void atom::assign_ADPs(std::vector<double> &second, std::vector<double> &third, std::vector<double> &fourth) {
	if (second.size() != 6) {
		std::cout << "Wrong size of second order ADP!" << std::endl;
		return;
	}
	if (third.size() == 0 && second.size() != 0) {
		ADPs.resize(1);
		ADPs[0].resize(6);
		ADPs[0] = second;
	}
	else if (third.size() != 0 && second.size() != 0) {
		if (third.size() != 10) {
			std::cout << "Wrong size of third order ADP!" << std::endl;
			return;
		}
		if (fourth.size() != 15) {
			std::cout << "Wrong size of fourth order ADP!" << std::endl;
			return;
		}
		ADPs.resize(3);
		ADPs[0].resize(6);
		ADPs[1].resize(10);
		ADPs[2].resize(15);
		ADPs[0] = second;
		ADPs[1] = third;
		ADPs[2] = fourth;
	}
};