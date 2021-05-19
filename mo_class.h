#pragma once
#ifndef __MO_CLASS_H__
#define __MO_CLASS_H__

#include <vector>
#include <string>


class MO {
	private:
		int nr;
		double occ;
		double ener;
		std::vector<double> coefficients;
	public:
		MO();
		MO (int number,double occupation,double energy);
		bool push_back_coef(double val, const int nex);
		void push_back_coef(double val);
		bool erase_coef(int nr, const int nex);
		double get_coefficient(const int nr, const bool debug) const;
		double get_coefficient_f(const int nr) const;
		double* get_coefficient_ptr();
		void change_coefficient(int nr);
		bool change_coefficient(int nr, double value, bool debug);
		void set_nr(const int inr) { nr=inr; };
		void set_occ(const int iocc) { occ=iocc; };
		double get_occ() const { return occ; };
		void set_ener(const double iener) {ener=iener;};
		int get_primitive_count() const {return coefficients.size();};
		std::string hdr();
		double get_energy() const;
		double* get_ptr_coefficients() { return &coefficients[0]; };
		void assign_coefs(std::vector<double>& values) { coefficients.resize(values.size()); coefficients = values; }
};

#endif
