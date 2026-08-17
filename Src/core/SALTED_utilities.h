#pragma once

#include "convenience.h"
#include "cube.h"
#include "featomic.hpp"
#include "constants.h"
#include "metatensor.h"

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
    std::vector<std::string> filter_species(const std::vector<std::string> &atomic_symbols, const std::vector<std::string> &species);
    void set_lmax_nmax(std::unordered_map<std::string, int> &lmax, std::unordered_map<std::string, int> &nmax, const std::array<std::vector<primitive>, 118> &basis_set, std::vector<std::string> species);
    int get_lmax_max(std::unordered_map<std::string, int> &lmax);

    inline featomic::SimpleSystem gen_featomic_system(const std::filesystem::path& filepath)
    {
        featomic::SimpleSystem featomic_system;
        WFN wfn = WFN(filepath);
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

//Calc density from RI fit coefficients
const double calc_density_ML(const double& x,
                            const double& y,
                            const double& z,
                            const vec& coefficients,
                            const std::vector<atom>& atoms);
//Perform the calculation for only one atom
const double calc_density_ML(const double &x,
                             const double &y,
                             const double &z,
                             const vec &coefficients,
                             const std::vector<atom> &atoms,
                             const int &atom_nr);

vec calc_atomic_density(const std::vector<atom> &atoms, const vec &coefs);

cube calc_cube_ML(const vec& data, WFN &dummy, const int& atom_nr = -1);
void calc_cube_ML(const vec& data, WFN& dummy, cube& cube_data, const int& atom_nr = -1);

void create_SALTED_training_data(const WFN& orbital, const WFN& aux);
