#pragma once
#ifndef __WFN_CLASS_H__
#define __WFN_CLASS_H__


class MO;
#include "atoms.h"
#include "cube.h"

#include <vector>
#include <string>
#include <fstream>

class WFN {
	private:
		int ncen;
		int nfunc;
		int nmo;
		int nex;
		int charge;
		unsigned int multi;
		int origin; //1=CRYSTAL; 2=WFN; 3=CUBE; 4=FFN; 5=FCHK
		std::string basis_set_name;
		std::string comment;
		std::string path;
		std::string method;
		
		std::vector<MO> MOs;
		std::vector<int> centers;
		std::vector<int> types;
		std::vector<double> exponents;
		std::vector<double> DensityMatrix;
		std::vector<double> SpinDensityMatrix;

		bool erase_center(int nr);
		bool erase_type(int nr);
		bool erase_exponent(int nr);
		bool push_back_center(int cent);
		bool push_back_type(int type);
		bool push_back_exponent(double e);
		bool modified;
		bool d_f_switch;
		bool distance_switch;
	public:
		WFN() ;
		WFN(int given_origin);
		std::vector<cube> cub;
		std::vector<atom> atoms;
		
		//--------------------MO handling--------------------------------------------
		void change_MO_coef(int nr);
		bool change_MO_coef(int nr_mo, int nr_primitive, double value, bool debug);
		double get_MO_coef(int nr_mo, int nr_primtive, bool debug) const; 
		int get_MO_primitive_count(int nr_mo) const; 
		bool push_back_MO(int nr, double occ, double ener);
		bool push_back_MO_coef(int nr, double value, int nr2);
		double get_MO_energy(int mo) const;
		int get_center(int nr) const {return centers[nr];};
		int get_type(int nr) const {return types[nr];};
		double get_MO_occ(int nr);
		
		//--------------------in and output----------------------------------------
		void change_basis_set_name(std::string name) {basis_set_name=name;};
		bool add_primitive(int cent, int type, double e, double * values);
		bool add_exp(int cent, int type, double e);
		bool read_wfn(std::string fileName, bool debug);
		bool read_wfn(std::string fileName, bool debug, std::ofstream &file);
		bool write_wfn(const std::string &fileName, bool debug, bool occupied);
		bool set_path(std::string given_path) { path=given_path; return true; };
		void print_primitive(int nr);
		void assign_charge(int i_charge) { charge=i_charge; };
		void assign_multi(int i_multi) { multi=i_multi; };
		int get_charge() const { return charge; };
		int get_multi() const { return multi; };
		int get_nex() const { return nex; };
		int get_ncen() const { return ncen; };
		int get_nmo() const { return nmo; };
		int get_nmo(bool only_occ);
		int get_origin() const { return origin; };
		double get_exponent(int nr) const {return exponents[nr];};
		unsigned int get_nr_electrons(bool &debug);
		/*int get_type(int nr) { if(nr>types.size()||nr<0) return -1; else return types[nr]; };
		int get_center(int nr) { if(nr>centers.size()||nr<0) return -1; else return centers[nr]; };					NOT NEEDED AT THIS POINT!
		double get_exponent(int nr) { if(nr>exponents.size()||nr<0) return -1; else return exponents[nr]; }; */
		std::string get_centers(bool bohr);
		std::string get_basis_set_name(){ return basis_set_name; };
		void set_basis_set_name( std::string input ){basis_set_name=input;};
		std::string get_path() { return path; };
		std::string hdr(bool occupied);
		void set_method(std::string input) { method = input; };
		std::string get_method() { return method; };
		bool erase_atom(int nr);
		void list_primitives();
		void list_centers();
		bool remove_center(int nr);
		bool remove_primitive(int nr);
		void change_center(int nr);
		void change_type(int nr);
		void change_exponent(int nr);
		/*bool change_center(int nr, int value);
		bool change_type(int nr, int value);					NOT NEEDED AT THIS POINT!
		bool change_exponent(int nr, double value);*/
		void set_modified(){ modified=true; };
		bool get_modified() { return modified; };
		void set_d_f_switch(bool in) { d_f_switch=in; };
		bool get_d_f_switch() { return d_f_switch; };
		int check_order(bool debug);
		bool sort_wfn(int order, bool debug);
		void set_dist_switch(){ distance_switch = true; };
		void operator=(const WFN &right);
		int calculate_charge();
		int calculate_charge(std::ofstream &file);
		bool guess_multiplicity(bool expert = false);
		bool guess_multiplicity(std::ofstream& file,bool expert = false);
		//-------------------atom handling--------------------------------------------------------------
		double get_atom_coordinate(int nr, int axis, bool debug);
		std::string get_atom_label(int nr);
		int get_nr_basis_set_loaded();
		bool get_atom_basis_set_loaded(int nr);
		double get_atom_basis_set_exponent(int nr_atom, int nr_prim);
		double get_atom_basis_set_coefficient(int nr_atom, int nr_prim);
		bool change_atom_basis_set_exponent (int nr_atom, int nr_prim, double value);
		bool change_atom_basis_set_coefficient (int nr_atom, int nr_prim, double value);
		int get_atom_primitive_count(int nr);
		int get_atom_primitive_type(int nr_atom, int nr_prim)
			{if(nr_atom<atoms.size()&&nr_atom>=0&&nr_prim<atoms[nr_atom].basis_set.size()&&nr_prim>=0)
					return atoms[nr_atom].basis_set[nr_prim].type; else return -1;};
		bool erase_atom_primitive(int nr, int nr_prim);
		int get_atom_shell_count(int nr);
		int get_atom_shell_primitives (int nr_atom, int nr_shell);
		int get_shell_type(int nr_atom, int nr_shell);
		int get_shell_center(int nr_atom, int nr_shell);
		int get_basis_set_shell(int nr_atom, int nr_prim);
		int get_shell_start(int nr_atom, int nr_shell, bool debug);
		int get_shell_start_in_primitives (int nr_atom, int nr_shell);
		int get_shell_end(int nr_atom, int nr_shell, bool debug);
		bool push_back_atom(std::string label, double x, double y, double z, int charge);
		//atom get_atom(int nr) { if(nr >= 0 && nr < ncen) return atoms[nr]; else return atom(); };
		bool push_back_atom_basis_set (int nr, double exponent, double coefficient, int type, int shell){
			if(nr<=ncen&&nr>=0) return atoms[nr].push_back_basis_set(exponent, coefficient, type, shell);
			else return false;
		};
		void print_atom_long(int nr){if(nr<=ncen&&nr>=0) atoms[nr].print_values_long();};
		int get_atom_charge(int nr);
		unsigned int get_atom_integer_mass(unsigned int atomnr);
		double get_atom_real_mass(int atomnr);
		//----------DM Handling--------------------------------
		void push_back_DM(double value);
		bool set_DM(int nr, double value);
		double get_DM(int nr);
		int get_DM_size() { return DensityMatrix.size(); };
		//----------S_DM Handling--------------------------------
		void push_back_SDM(double value);
		bool set_SDM(int nr, double value);
		double get_SDM(int nr);
		int get_SDM_size() { return SpinDensityMatrix.size(); };
		//-----------Cube handling-------------------------
		bool push_back_cube(std::string filepath, bool full, bool expert = false);
		void push_back_cube(cube given) {cub.push_back(given);};
		void pop_back_cube();
		
		//-----------Pointer to members---------------------
		int * get_ptr_types() {return &types[0];};
		int * get_ptr_centers() {return &centers[0];};
		double * get_ptr_exponents() {return &exponents[0];};
		double * get_ptr_mo_coefficients(int mo);

};


#include "mo_class.h"


#endif
