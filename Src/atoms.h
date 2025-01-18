#pragma once
#include "convenience.h"

//-----------------Definition of atoms and basis sets--------------------

class basis_set_entry{
private:
	double coefficient;
	double exponent;
	unsigned int type; //1=S; 2=P; 3=D; 4=F; 5=G
	unsigned int shell;
	primitive p;
public:
	// Default assignment operator
	basis_set_entry& operator=(const basis_set_entry& rhs) = default;
	basis_set_entry();
	basis_set_entry(double g_coefficient, double g_exponent, unsigned int g_type, unsigned int g_shell);
    // Getter functions
	
    double get_coefficient() const { return coefficient; };
    double get_exponent() const { return exponent; };
    unsigned int get_type() const { return type; };
    unsigned int get_shell() const { return shell; };
    primitive get_primitive() const { return p; };
    // Setter functions
	 
    void set_coefficient(const double& value) { coefficient = value; };
    void set_exponent(const double& value) { exponent = value; };
    void set_type(const unsigned int& value) { type = value; };
    void set_shell(const unsigned int& value) { shell = value; };
    void set_primitive(const primitive& value) { p = value; };

	// Equality operator
	bool operator==(const basis_set_entry& other) const {
		return coefficient == other.get_coefficient() &&
			exponent == other.get_exponent() &&
			type == other.get_type() &&
			shell == other.get_shell() &&
			p == other.get_primitive();
	}
};

class atom {
private:
	std::string label;
	std::string ID;
	int nr, charge, ECP_electrons;
	double x, y, z;
	std::array<double,3> frac_coords;
	std::vector<basis_set_entry> basis_set;
	int basis_set_id;
	std::vector<unsigned int> shellcount;
	//The Order is:
	//[0] = second order (U11, U22, U33, U12, U13, U23)
	//[1] = third order  (C111, C112, C113, C122, C123, C133, C222, C223, C233, C333)
	//[2] = fourth order (D1111, D1112, D1113, D1122, D1123, D1133, D1222, D1223, D1233, D1333, D2222, D2223, D2233, D2333, D3333)
	vec2 ADPs;
public:
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
    basis_set_entry get_basis_set_entry(const int& nr) const { return basis_set[nr]; };
    double get_coordinate(const unsigned int& axis) const;
	int get_charge() const { return charge; };
    int get_ECP_electrons() const { return ECP_electrons; };
    void set_charge(const int& ch) { charge = ch; };
    void set_ECP_electrons(const int& ECP_els) { ECP_electrons = ECP_els; };
    void set_coordinate(const unsigned int& axis, const double& value);
    void set_frac_coords(const std::array<double, 3>& frac);
    std::string get_label() const { return label; };
    int get_nr() const { return nr; };
    void set_label(const std::string& l) { label = l; };
    void set_nr(const int& n) { nr = n; };
    void clear_basis_set() { basis_set.clear(); };
    int get_basis_set_size() const { return (int)basis_set.size(); };
    double get_basis_set_exponent(const int& nr) const { return basis_set[nr].get_exponent(); };
    double get_basis_set_coefficient(const int& nr) const { return basis_set[nr].get_coefficient(); };
    void set_basis_set_exponent(const int& nr, const double& value) { basis_set[nr].set_exponent(value); };
    void set_basis_set_coefficient(const int& nr, const double& value) { basis_set[nr].set_coefficient(value); };
    void set_basis_set_type(const int& nr, const int& value) { basis_set[nr].set_type(value); };
    void set_basis_set_shell(const int& nr, const int& value) { basis_set[nr].set_shell(value); };
    std::vector<basis_set_entry> get_basis_set() const { return basis_set; };
    void erase_basis_set(const unsigned int& nr) { basis_set.erase(basis_set.begin() + nr); };
    int get_basis_set_type(const int& nr) const { return basis_set[nr].get_type(); };
    int get_basis_set_shell(const int& nr) const { return basis_set[nr].get_shell(); };
    void set_basis_set_id(const int& id) { basis_set_id = id; };
    int get_basis_set_id() const { return basis_set_id; };
    void set_shellcount(const std::vector<unsigned int>& sc) { shellcount = sc; };
    std::vector<unsigned int> get_shellcount() const { return shellcount; };
    unsigned int get_shellcount_size() const { return (unsigned int)shellcount.size(); };
    void set_shellcount(const unsigned int& nr, const unsigned int& value) { shellcount[nr] = value; };
    unsigned int get_shellcount(const unsigned int& nr) const { return shellcount[nr]; };
    void clear_shellcount() { shellcount.clear(); };
    void set_ADPs(const vec2& adps) { ADPs = adps; };
    vec2 get_ADPs() const { return ADPs; };


	bool operator==(const atom& other) const;
};