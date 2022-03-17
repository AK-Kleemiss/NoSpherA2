#pragma once

#include <vector>
#include <string>

class WFN;

class cube
{
	public:
		cube();
		cube(int x, int y, int z, int g_na = 0, bool grow_values = false);
		cube(std::string filepath, bool read, WFN &wave, bool &return_value, bool expert = false);
		cube(std::string filepath, bool read, bool expert = false);
		cube(int g_na, const std::vector <int> &g_size, const std::vector <double> &g_origin, const std::vector < std::vector<double> > &g_vectors, const std::vector<std::vector<std::vector<double> > > &g_values);
		cube(const cube &given);
		std::string path;
		int get_size(int direction) const;
		bool get_loaded() const{ return loaded; };
		cube operator + (cube &right) const;
		cube operator - (cube &right) const;
		cube operator * (cube &right) const;
		cube operator / (cube &right) const;
		bool operator += (cube &right);
		bool operator -= (cube &right);
		bool operator *= (cube &right);
		bool operator /= (cube &right);
		void operator = (cube &right);
		bool mask (cube &right);
		bool mask (cube &right, double thresh);
		bool negative_mask (cube &right);
		double rrs(cube &right);
		double sum();
		double diff_sum();
		std::vector<double> double_sum();
		double get_value(int x, int y, int z) const;
		bool set_value(int x, int y, int z, double value);
		bool read_file(bool full, bool header, WFN &wave, bool expert = false);
		bool read_file(bool full, bool header, bool expert = false);
		bool write_file(WFN &wave, bool force = false, bool absolute = false);
		bool write_file(WFN &wave, std::string &given_path, bool debug);
		bool write_xdgraph(WFN &wave, std::string &given_path, bool debug);
		bool fractal_dimension(const double stepsize);
		double get_vector(int i, int j) const;
		bool set_vector(int i, int j, double value);
		double get_origin(unsigned int i) const;
		bool set_origin(unsigned int i, double value);
		int get_na() const {return na;};
		void set_na(int g_na) {na=g_na;};
		std::string super_cube(WFN &wave);
		void set_comment1 (std::string input) {comment1=input;};
		void set_comment2 (std::string input) {comment2=input;};
		void set_zero();
		std::string get_comment1 () const {return comment1;};
		std::string get_comment2 () const {return comment2;};
	private:
		int na;
		bool loaded;
		std::string comment1;
		std::string comment2;
		std::vector <int> size;
		std::vector <double> origin;
		std::vector < std::vector <double> > vectors;
		std::vector < std::vector < std::vector <double> > > values;
};

#include "wfn_class.h"