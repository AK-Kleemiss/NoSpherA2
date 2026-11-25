#pragma once

#include "convenience.h"
#include "cube.h"
#include "aux_density_table.h"
#include "constants.h"

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
}

cube calc_cube_ML(const vec& data, WFN &dummy, const int& atom_nr = -1);
void calc_cube_ML(const vec& data, WFN& dummy, cube& cube_data, const int& atom_nr = -1);

void create_SALTED_training_data(const WFN& orbital, const WFN& aux, const options& opts);
//geometry_aid.h still takes the featomic declarations from here
#include "SALTED_descriptors.h"
