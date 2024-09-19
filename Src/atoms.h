#pragma once
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
	bool operator==(const basis_set_entry& other) const;
	primitive p;
};

class atom {
public:
	std::string label;
    std::string ID;
	int nr, charge, ECP_electrons;
	double x, y, z;
	vec frac_coords;
	atom();
	atom(const std::string &l, const std::string& id, const int &n, const double &c1, const double &c2, const double &c3, const int &ch);
	atom(const std::string& l, const std::string& id, const int& n, const double& c1, const double& c2, const double& c3, const int& ch, const int& ECP_els);
	atom operator=(const atom &rhs);
	void print_values() const;
	bool push_back_basis_set(const double &coefficient, const double &exponent, const int &type, const int &shell);
	void print_values_long() const;
	bool get_basis_set_loaded() const;
	bool is_anharm() const;
	void assign_ADPs(vec &second, vec &third, vec &fourth);
	void assign_ADPs(vec& second);
	void assign_ADPs(double& Uiso);
    void assign_ID(const std::string& id);
    void set_ID(const std::string& id);
    std::string get_ID() const;
	std::vector<basis_set_entry> basis_set;
	int basis_set_id;
	std::vector<unsigned int> shellcount;
	//The Order is:
	//[0] = second order (U11, U22, U33, U12, U13, U23)
	//[1] = third order  (C111, C112, C113, C122, C123, C133, C222, C223, C233, C333)
	//[2] = fourth order (D1111, D1112, D1113, D1122, D1123, D1133, D1222, D1223, D1233, D1333, D2222, D2223, D2233, D2333, D3333)
	vec2 ADPs;

	bool operator==(const atom& other) const;
};