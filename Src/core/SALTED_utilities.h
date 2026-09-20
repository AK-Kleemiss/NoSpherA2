#pragma once

#include "convenience.h"
#include "cube.h"
#include "featomic.hpp"
#include "constants.h"
#include "metatensor.h"

//Predefine the SALTEDConfig struct, since it is used in the SALTED_Utils namespace
struct SALTEDConfig;

// Stores one contiguous slab per angular momentum. equicomb fixes l1/l2 for
// substantial stretches of work, so this avoids the thousands of tiny vectors
// in the former atom/channel/l/m representation.
class SALTEDDescriptors {
public:
    SALTEDDescriptors() = default;

    SALTEDDescriptors(const int natoms, const int nchannels, const int lmax)
        : _nchannels(nchannels), _offsets(static_cast<size_t>(lmax + 1)) {
        size_t offset = 0;
        for (int l = 0; l <= lmax; ++l) {
            _offsets[static_cast<size_t>(l)] = offset;
            offset += static_cast<size_t>(natoms) * nchannels * (2 * static_cast<size_t>(l) + 1);
        }
        _values.assign(offset, cdouble{ 0.0, 0.0 });
    }

    cdouble* block(const int atom, const int channel, const int l) noexcept {
        return _values.data() + _offsets[static_cast<size_t>(l)]
            + (static_cast<size_t>(atom) * _nchannels + channel) * (2 * static_cast<size_t>(l) + 1);
    }

    const cdouble* block(const int atom, const int channel, const int l) const noexcept {
        return _values.data() + _offsets[static_cast<size_t>(l)]
            + (static_cast<size_t>(atom) * _nchannels + channel) * (2 * static_cast<size_t>(l) + 1);
    }

    std::vector<cdouble>& values() noexcept { return _values; }

    //Flat view for the GPU path: block(atom, channel, l) is offsets[l] +
    //(atom * nchannels + channel) * (2l+1), which the device reproduces
    const std::vector<cdouble>& values() const noexcept { return _values; }
    const std::vector<size_t>& offsets() const noexcept { return _offsets; }
    int nchannels() const noexcept { return _nchannels; }

    void clear() noexcept {
        _nchannels = 0;
        _offsets.clear();
        _values.clear();
    }

    void shrink_to_fit() {
        _offsets.shrink_to_fit();
        _values.shrink_to_fit();
    }

private:
    int _nchannels = 0;
    std::vector<size_t> _offsets;
    std::vector<cdouble> _values;
};

namespace SALTED_Utils
{
    std::vector<cvec2> complex_to_real_transformation(std::vector<int> sizes);
    //Removes the atoms the model cannot predict from wavy; returns which of the input atoms were removed
    std::vector<char> filter_input(WFN& wavy, options& opt, const SALTEDConfig& config);
    void set_lmax_nmax(std::unordered_map<std::string, int> &lmax, std::unordered_map<std::string, int> &nmax, const std::array<std::vector<primitive>, 118> &basis_set, std::vector<std::string> species);
    int get_lmax_max(std::unordered_map<std::string, int> &lmax);

    inline featomic::SimpleSystem gen_featomic_system(const WFN& wfn)
    {
        featomic::SimpleSystem featomic_system;
        for (const atom& a : *wfn.get_atoms_ptr())
        {
            d3 xyz = { constants::bohr2ang(a.get_coordinate(0)),
                                          constants::bohr2ang(a.get_coordinate(1)),
                                          constants::bohr2ang(a.get_coordinate(2)) };
            featomic_system.add_atom(a.get_charge(), xyz);
        }
        return featomic_system;
    }

    struct FeatomicHyperParameters
    {
        struct RadialBasis
        {
            std::string type = "Gto";
            double spline_accuracy = 1e-6;
        };

        struct CutoffFunction
        {
            std::string type = "ShiftedCosine";
            double width = 0.1;
        };

        double cutoff_radius = 2.0;
        int max_radial = 3;
        int max_angular = 3;
        double atomic_gaussian_width = 0.5;
        double center_atom_weight = 1.0;

        std::vector<std::string> species = { "H", "C", "N", "O" };
        std::vector<std::string> neighspe = { "H", "C", "N", "O" };

        RadialBasis radial_basis{};
        CutoffFunction cutoff_function{};

        std::string to_json() const;
    };

    SALTEDDescriptors calculate_SALTED_descriptors(const featomic::SimpleSystem& featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters);
    metatensor::TensorMap calculate_SOAP_Powerspectrum(featomic::SimpleSystem featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters);
}

//Flat copy of an aux basis for evaluating the fitted density on many points, see aux_density.h
struct aux_density_table
{
    int n_at = 0, n_sh = 0, n_pr = 0, n_coef = 0;

    vec cx, cy, cz, r2_max, pr_exp, pr_norm;
    // alpha^(l + 3/2), needed for Fourier-Bessel transform
    vec pr_exp_l32;
    ivec sh_start, sh_atom, sh_l, pr_start, coef_off;
    // nuclear charges less ECP electrons, for the potential
    ivec Z;
    // distinct (exponent, l) pairs and the slot of each primitive: the radial Fourier factor depends on nothing else
    vec uniq_exp, uniq_exp_l32;
    ivec uniq_l, pr_uniq;

    // coefficient index -> shell / m
    ivec coef_shell;
    ivec coef_m;

    explicit aux_density_table(const std::vector<atom>& atoms);
    double operator()(const double x, const double y, const double z, const double* coefs) const;
    double operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz) const;
    double operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz, double& lap) const;
    //rho with its gradient and Hessian (row-major 3x3)
    double operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz, double* H) const;
    //Electrostatic potential of nuclei and fitted density, see aux_density::esp_at
    double esp(const double x, const double y, const double z, const double* coefs) const;
    double lap(const double x, const double y, const double z, const double* coefs) const;
    double eli(const double x, const double y, const double z, const double* coefs) const;

    // Convenience function for one atom
    cdouble fourier_atom(double kx, double ky, double kz, const double* coefs, int atom_idx) const;

    // Integral:
    //
    // ∫_0^inf R_l(r) r^(l+2) dr
    //
    // where
    //
    // R_l(r) = sum_p pr_norm[p] exp(-alpha_p r²)
    //
    double shell_radial_moment(int shell) const;

    // Full 3D integral of an s-shell including Y_00.
    // Only valid for l = 0.
    double shell_population_integral(int shell) const;
};
//The fitted density of an aux_density_table (RI fit or SALTED prediction, the same path) on np points:
//rho, with the gradient when gx/gy/gz are given, the Laplacian with lap and the row-major 3x3 Hessian
//with hess (9 per point, needs the gradient arrays); GPU kernel when one is present and enabled, else OpenMP
void calc_aux_density(const aux_density_table& t, const vec& coefficients, const int np, const double* x, const double* y, const double* z, double* rho, double* gx = nullptr, double* gy = nullptr, double* gz = nullptr, double* lap = nullptr, double* hess = nullptr);

vec calc_atomic_density(const std::vector<atom> &atoms, const vec &coefs);

// Scale the l=0 coefficients so the predicted density integrates to the exact
// electron count. Returns the applied factor (1.0 if nothing was done).
double apply_charge_constraint(const std::vector<atom> &atoms, vec &coefs,
                               int net_charge, bool spherical_fill_used,
                               int n_filled, double filled_eeq_charge,
                               double applied_fill_charge, std::ostream &file);

cube calc_cube_ML(const vec& data, WFN &dummy, const int& atom_nr = -1);
void calc_cube_ML(const vec& data, WFN& dummy, cube& cube_data, const int& atom_nr = -1);

void create_SALTED_training_data(const WFN& orbital, const WFN& aux, const options& opts);
