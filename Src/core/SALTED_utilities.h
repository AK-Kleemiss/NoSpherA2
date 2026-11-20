#pragma once

#include "convenience.h"
#include "cube.h"
#include "aux_density_table.h"
#include "featomic.hpp"
#include "constants.h"
#include "metatensor.h"

struct SALTEDConfig;

// One contiguous slab per angular momentum, since equicomb holds l1/l2 fixed over long stretches.
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

	cvec& values() noexcept { return _values; }

	//Flat view for the GPU path: block(atom, channel, l) is offsets[l] +
	//(atom * nchannels + channel) * (2l+1), which the device reproduces
	const cvec& values() const noexcept { return _values; }
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
	cvec _values;
};

namespace SALTED_Utils
{
	cvec3 complex_to_real_transformation(ivec sizes);
	//Removes the atoms the model cannot predict from wavy; returns which of the input atoms were removed
	std::vector<char> filter_input(WFN& wavy, options& opt, const SALTEDConfig& config);
	void set_lmax_nmax(std::unordered_map<std::string, int> &lmax, std::unordered_map<std::string, int> &nmax, const BasisSet& basis_set, std::vector<std::string> species);
	int get_lmax_max(std::unordered_map<std::string, int> &lmax);

	// skip[a] != 0 leaves atom a out, z[a] replaces its nuclear charge as the featomic type
	inline featomic::SimpleSystem gen_featomic_system(const WFN& wfn, const std::vector<char>& skip = {}, const ivec& z = {})
	{
		featomic::SimpleSystem featomic_system;
		const auto& atoms = *wfn.get_atoms_ptr();
		for (size_t i = 0; i < atoms.size(); i++)
		{
			if (i < skip.size() && skip[i]) continue;
			const atom& a = atoms[i];
			d3 xyz = { constants::bohr2ang(a.get_coordinate(0)),
										  constants::bohr2ang(a.get_coordinate(1)),
										  constants::bohr2ang(a.get_coordinate(2)) };
			featomic_system.add_atom(z.empty() ? a.get_charge() : z[i], xyz);
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

	// centres: the featomic atoms that get a descriptor row, all when empty
	SALTEDDescriptors calculate_SALTED_descriptors(const featomic::SimpleSystem& featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters, const std::vector<int32_t>& centres = {});
	metatensor::TensorMap calculate_SOAP_Powerspectrum(featomic::SimpleSystem featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters);
}

cube calc_cube_ML(const vec& data, WFN &dummy, const int& atom_nr = -1);
void calc_cube_ML(const vec& data, WFN& dummy, cube& cube_data, const int& atom_nr = -1);

void create_SALTED_training_data(const WFN& orbital, const WFN& aux, const options& opts);
