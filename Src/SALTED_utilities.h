#pragma once

#include "convenience.h"


#if defined(_WIN32) || defined(__RASCALINE__)
#include "rascaline.hpp"
#include "metatensor.h"
#endif

namespace SALTED_Utils
{
    std::vector<cvec2> complex_to_real_transformation(std::vector<int> sizes);
    std::vector<std::string> filter_species(const std::vector<std::string>& atomic_symbols, const std::vector<std::string>& species);
    void set_lmax_nmax(std::unordered_map<std::string, int>& lmax, std::unordered_map<std::string, int>& nmax, std::array<std::vector<primitive>, 35> basis_set, std::vector<std::string> species);
    int get_lmax_max(std::unordered_map<std::string, int>& lmax);
}


// ALL below RASCALINE

#if defined(_WIN32) || defined(__RASCALINE__)
class Rascaline_Descriptors
{
public:
    Rascaline_Descriptors(std::string filepath, int nrad, int nang, double atomic_gaussian_width,
                          double rcut, int n_atoms, std::vector<std::string> neighspe, std::vector<std::string> species,
                          double center_atom_weight = 1.0, double spline_accuracy = 1e-6, double cutoff_width = 0.1);
    std::string filepath;
    int nrad;
    int nang;
    double atomic_gaussian_width;
    double rcut;
    int n_atoms;
    std::vector<std::string> neighspe;
    std::vector<std::string> species;
    double center_atom_weight = 1.0;
    double spline_accuracy = 1e-6;
    double cutoff_width = 0.1;
      
    cvec4 calculate_expansion_coeffs();
private:
    int nspe;
  

    struct RadialBasis
    {
        std::string type;
        double spline_accuracy;
    };

    struct CutoffFunction
    {
        std::string type;
        double width;
    };

    struct HyperParametersDensity
    {
        double cutoff;
        int max_radial;
        int max_angular;
        double atomic_gaussian_width;
        double center_atom_weight;
        RadialBasis radial_basis;
        CutoffFunction cutoff_function;
    };

    //Helper functions
    std::string to_json(const HyperParametersDensity& params);
    std::string gen_parameters();

    //Used to generate metatensor::TensorMap and save the buffer location into the descriptor_buffer
    metatensor::TensorMap get_feats_projs();
    //Reads the descriptor buffer and fills the expansion coefficients vector
    cvec4 get_expansion_coeffs(std::vector<uint8_t> descriptor_buffer);
};
#endif