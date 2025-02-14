#include "pch.h"
#include "convenience.h"
#include "atoms.h"

//-----------------Definition of atoms and basis sets--------------------

basis_set_entry::basis_set_entry() {
	coefficient = 0.0;
	exponent = 0.0;
	type = 0;
	shell = 0;
};

basis_set_entry::basis_set_entry(double g_coefficient, double g_exponent, unsigned int g_type, unsigned int g_shell) {
	coefficient = g_coefficient;
	exponent = g_exponent;
	type = g_type;
	shell = g_shell;
	p = primitive(0, type, exponent, coefficient);
};

atom::atom() {
	label = '?';
	ID = "0000000000000";
	nr = 0;
	x = 0.0;
	y = 0.0;
	z = 0.0;
	charge = 0;
	basis_set_id = 0;
	ECP_electrons = 0;
	frac_coords = { 0.,0.,0. };
};

atom::atom(const std::string& l, 
	const std::string& id, 
	const int& n, 
	const double& c1, 
	const double& c2, 
	const double& c3, 
	const int& ch) : nr(n), label(l), ID(id), x(c1), y(c2), z(c3), charge(ch), ECP_electrons(0), basis_set_id(0), frac_coords({ 0,0,0 }) 
{};

atom::atom(const std::string& l, 
	const std::string& id, 
	const int& n, 
	const double& c1, 
	const double& c2, 
	const double& c3, 
	const int& ch, 
	const int& ECP_els) : nr(n), label(l), ID(id), x(c1), y(c2), z(c3), charge(ch), ECP_electrons(ECP_els), basis_set_id(0), frac_coords({ 0,0,0 })
{};

atom atom::operator= (const atom& rhs) {
	label = rhs.label;
	ID = rhs.ID;
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

void atom::print_values() const {
	std::cout << "nr: " << nr << " label: " << label << " x: " << x << " y: " << y << " z: " << z << " charge: " << charge << std::endl;
};

bool atom::is_anharm() const {
	return ADPs.size() > 1;
};

void atom::print_values_long() const {
	std::cout << "nr: " << nr << " label: " << label << " x: " << x << " y: " << y << " z: " << z << " charge: " << charge << std::endl;
	if (basis_set.size() > 0) {
		for (unsigned int i = 0; i < basis_set.size(); i++) {
			std::cout << "basis set entry " << i << ": expo: " << basis_set[i].get_exponent() << " coef: " << basis_set[i].get_coefficient() << " type: " << basis_set[i].get_type() << " shell: " << basis_set[i].get_shell() << std::endl;
		}
	}
	if (shellcount.size() > 0) {
		for (unsigned int i = 0; i < shellcount.size(); i++) {
			std::cout << "shellcount " << i << ": value: " << shellcount[i] << std::endl;
		}
	}
};

bool atom::push_back_basis_set(const double& exponent, const double& coefficient, const int& type, const int& shell) {
	if (shell == shellcount.size())
		shellcount.push_back((unsigned int)1);
	else
		shellcount[shell]++;
	if (type >= 0 && shell >= 0) {
		basis_set.push_back(basis_set_entry(coefficient, exponent, type, shell));
		return true;
	}
	else {
		if (type >= 10) err_checkf(false, "h and higher types are not yet supported!", std::cout);
		std::cout << "This is not a valid basis set entry!" << std::endl;
		std::cout << "Exponent: " << exponent << " coefficient: " << coefficient << " type: " << type << " shell: " << shell << std::endl;
		return false;
	}
};

bool atom::get_basis_set_loaded() const {
	if (basis_set.size() > 0) return true;
	else return false;
};

void atom::assign_ADPs(double& Uiso) {
	ADPs.resize(1);
	ADPs[0].resize(6);
	ADPs[0][0] = ADPs[0][1] = ADPs[0][2] = Uiso;
};

void atom::assign_ADPs(vec& second) {
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

void atom::assign_ADPs(vec& second, vec& third, vec& fourth) {
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

void atom::assign_ID(const std::string& id) {
	ID = id;
};

void atom::set_ID(const std::string& id) {
	ID = id;
};

std::string atom::get_ID() const {
	return ID;
};

double atom::get_coordinate(const unsigned int& axis) const {
    if (axis == 0) return x;
    else if (axis == 1) return y;
    else if (axis == 2) return z;
    else return 0.0;
};

void atom::set_coordinate(const unsigned int& axis, const double& value) {
    if (axis == 0) x = value;
    else if (axis == 1) y = value;
    else if (axis == 2) z = value;
};

void atom::set_frac_coords(const std::array<double, 3>& frac) {
    frac_coords = frac;
};

bool atom::operator==(const atom& other) const {
	if (this == &other) {
		return true;
	}

	return label == other.label &&
		nr == other.nr &&
		charge == other.charge &&
		ECP_electrons == other.ECP_electrons &&
		x == other.x &&
		y == other.y &&
		z == other.z &&
		frac_coords == other.frac_coords &&
		basis_set == other.basis_set &&
		shellcount == other.shellcount &&
		ADPs == other.ADPs;
}