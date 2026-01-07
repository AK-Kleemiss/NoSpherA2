#pragma once
#include "convenience.h"
#include "wfn_class.h"

//Temporary basis set format
struct LibCintBasis {
    vec coefficients;
    vec exponents;
    int env_idx;
    ivec shellcount;
    ivec shelltypes;
};


class Int_Params {
private:
    unsigned long long ncen = 0;           // Number of centers (atoms)
    int wfn_origin = 0;
    std::vector<atom> atoms;

    std::map<int, LibCintBasis> basis_sets;  //Maps atomic number -> (vec(coefs), vec(expon))

    ivec _atm;              // Flattened ATM array
    ivec _bas;              // Flattened BAS array
    vec _env;

    unsigned long long nao = 0;
    unsigned long long nbas = 0;
    void calc_integration_parameters();
    void collect_basis_data();
    //   void collect_basis_data_from_gbw();
       //void collect_basis_data_internal();

    void populate_atm();
    void populate_bas();
    void populate_env();

    static double gaussian_int(int n, double exp) {
        double n1 = (n + 1) * 0.5;
        return std::tgamma(n1) / (2.0 * pow(exp, n1));
    }
public:
    Int_Params();
    Int_Params(const WFN& wavy);
    Int_Params(const Int_Params& first, const Int_Params& second);
    ~Int_Params() = default;// Combine two Int_Params objects

    // Getters
    int* get_ptr_atm() { return _atm.data(); };
    int* get_ptr_bas() { return _bas.data(); };
    double* get_ptr_env() { return _env.data(); };

    ivec get_atm() const { return _atm; };
    ivec get_bas() const { return _bas; };
    vec get_env() const { return _env; };

    unsigned long long get_nbas() const { return nbas; };  // Number of basis segments
    unsigned long long get_nao() const { return nao; };    // Number of atomic orbitals  
    unsigned long long get_natoms() const { return ncen; };

    static vec normalize_gto(vec coef, const vec& exp, const int l);

    void print_data(std::string name); //FOR DEBUG PURPOSES
    std::vector<atom> get_atoms() const {
        return atoms;
    }
};

