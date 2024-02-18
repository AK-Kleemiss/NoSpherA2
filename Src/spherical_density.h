#pragma once

#include <vector>
#include <iostream>

inline void not_implemented_SA(const std::string& file, const int& line, const std::string& function, const std::string& error_mesasge, std::ostream& log_file)
{
	log_file << function << " at: " << file << ":" << line << " " << error_mesasge << std::endl;
	log_file.flush();
	exit(-1);
};
#define err_not_impl_SA() not_implemented_SA(__FILE__, __LINE__, __func__, "Virtual_function", std::cout);

class Spherical_Atom {
protected:
	int atomic_number;
	int ECP_mode;
	const int first_ex();
	virtual const int previous_element_coef();
	const int* nex, * ns, * np, * nd, * nf, * occ, * n;
	const double* z, * c;
	int charge;
	int _prev_coef, _offset, _first_ex;
	virtual void calc_orbs(int& nr_ex,
		int& nr_coef,
		const double& dist,
		const int& offset,
		const int* n_vector,
		const int lower_m,
		const int upper_m,
		double* Orb) {
		err_not_impl_SA();
	};
	virtual double calc_type(
		int& nr_ex,
		int& nr_coef,
		const double& k_vector,
		const int& offset,
		const int* n_vector,
		const int lower_m,
		const int upper_m,
		const int& max,
		const int& min) {
		err_not_impl_SA();
		return -1;
	};
public:
	Spherical_Atom(const int g_atom_number, const int ECP_m = 1) : atomic_number(g_atom_number), 
		_offset((atomic_number-1) * 19), _first_ex(0), _prev_coef(0), c(NULL), n(NULL), nd(NULL), ns(NULL), np(NULL), nf(NULL), nex(NULL), occ(NULL), z(NULL) {
		ECP_mode = ECP_m;
		charge = 0;
	};
	Spherical_Atom() : _first_ex(0), _offset(0), c(NULL), n(NULL), nd(NULL), ns(NULL), np(NULL), nf(NULL), nex(NULL), occ(NULL), z(NULL) {
		ECP_mode = 1;
		atomic_number = 1;
		charge = 0;
	};
	virtual const double get_radial_density(double& dist) {
		err_not_impl_SA();
		return -1;
	};
	virtual const double get_form_factor(const double& k_vector) {
		err_not_impl_SA();
		return -1;
	};
	virtual const double get_core_form_factor(const double& k_vector, const int& core_els) {
		err_not_impl_SA();
		return -1;
	};
	virtual const double get_custom_form_factor(
		const double& k_vector,
		const int& max_s,
		const int& max_p,
		const int& max_d,
		const int& max_f,
		const int& min_s,
		const int& min_p,
		const int& min_d,
		const int& min_f) {
		err_not_impl_SA();
		return -1;
	};
	const int get_atomic_number() const { return atomic_number; };
	const int get_charge() const { return charge; };
};

class Thakkar : public Spherical_Atom {
protected:
	void calc_orbs(int& nr_ex,
		int& nr_coef,
		const double& dist,
		const int& offset,
		const int* n_vector,
		const int lower_m,
		const int upper_m,
		double* Orb) override;
	void calc_custom_orbs(int& nr_ex,
		int& nr_coef,
		const double& dist,
		const int& offset,
		const int* n_vector,
		const int lower_m,
		const int upper_m,
		const int& max,
		const int& min,
		double* Orb);
	double calc_type(
		int& nr_ex,
		int& nr_coef,
		const double& k_vector,
		const int& offset,
		const int* n_vector,
		const int lower_m,
		const int upper_m,
		const int& max,
		const int& min) override;
public:
	Thakkar(const int g_atom_number, const int ECP_mode = 1);
	Thakkar();
	const double get_radial_density(double& dist) override;
	const double get_radial_custom_density(
		const double& dist,
		const int& max_s,
		const int& max_p,
		const int& max_d,
		const int& max_f,
		const int& min_s,
		const int& min_p,
		const int& min_d,
		const int& min_f);
	const double get_form_factor(const double& k_vector) override;
	const double get_core_form_factor(const double& k_vector, const int& core_els) override;
	const double get_core_density(const double& dist, const int& core_els);
	const double get_custom_form_factor(
		const double& k_vector,
		const int& max_s,
		const int& max_p,
		const int& max_d,
		const int& max_f,
		const int& min_s,
		const int& min_p,
		const int& min_d,
		const int& min_f) override;
};

class Thakkar_Anion : public Thakkar {
public:
	Thakkar_Anion(const int g_atom_number);
};

class Thakkar_Cation : public Thakkar {
public:
	Thakkar_Cation(const int g_atom_number);
};

class Gaussian_Atom : public Spherical_Atom {
protected:
	const int* ng, * nh;
	int first_atomic_number;
	const int previous_element_coef() override;
	void calc_orbs(int& nr_ex,
		int& nr_coef,
		const double& dist,
		const int& offset,
		const int* n_vector,
		const int lower_m,
		const int upper_m,
		double* Orb) override;
	double calc_type(
		int& nr_ex,
		int& nr_coef,
		const double& k_vector,
		const int& offset,
		const int* n_vector,
		const int lower_m,
		const int upper_m,
		const int& max,
		const int& min) override;
public:
	Gaussian_Atom(const int g_atom_number, std::string& basis);
	Gaussian_Atom() = default;
	const double get_radial_density(double& dist) override;
	const double get_form_factor(const double& k_vector) override;
	const double get_core_form_factor(const double& k_vector, const int& core_els) override;
	const double get_custom_form_factor(
		const double& k_vector,
		const int& max_s,
		const int& max_p,
		const int& max_d,
		const int& max_f,
		const int& min_s,
		const int& min_p,
		const int& min_d,
		const int& min_f) override;
	const double get_custom_form_factor(
		const double& k_vector,
		const int& max_s,
		const int& max_p,
		const int& max_d,
		const int& max_f,
		const int& max_g,
		const int& max_h,
		const int& min_s,
		const int& min_p,
		const int& min_d,
		const int& min_f,
		const int& min_g,
		const int& min_h);
};