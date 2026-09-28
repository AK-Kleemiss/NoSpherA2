#pragma once
#include "convenience.h"
#include "cell.h"
#include "extinction.h"
#include "scattering_factors.h"
#include "wfn_class.h"

class structure_factors {

public:
	structure_factors() = default;
	structure_factors(const options& opt_in);

//private:

	// Store cristallographic quality criteria
	struct quality_criteria {
		double GooF1;
		double GooF2;
		double weighted_GooF1;
		double weighted_GooF2;
		double R1;
	};

	// Store scattering data for each reflection
	struct scatter_data {
		vec F_obs;
		vec F_obs2;
		vec sigma_obs;
		vec sigma_obs2;
		vec abs_F_obs;
		cvec F_calc;
		vec F_calc2;
		double scale;
		ivec hkl_mask;
		hkl_list hkl;
		hkl_list hkl_enlarged;
		cvec anom_correction;
	};

	// Store information about the model (e.g. number of atoms, reflections)
	struct model_data {
		int ncen;
		int nr;
		int nr_enlarged;
		int n_params = 177;
	};


	// Stores the Debye-Waller factors
	cvec2 DW_facts;
	// Stores the rotational phase factors
	cvec2 phase_facts;
	// Stores the translational phase factors
	cvec2 translation_phase_facts;
	// Used for conveying options to the structure_factors class
	const options* opt;
	// Store information about the model
	model_data model_data;
	// Store scattering data for each reflection
	scatter_data scatter_data;
	// Store cristallographic quality criteria
	quality_criteria quality_criteria;
	// Store the unit cell
	cell unit_cell;
	// Store the atoms of the asymmetric unit (or the grown unit)
	std::vector<asym_atom> asym_atoms;
	// Store the full rotational symmetry in case of a grown structure, for rotating ADPs
	ivec3 original_rotations;
	// Dummy wavefunction used for reading in atoms and basis sets
	WFN dummy_wave;
	// List used for grid generation
	ivec asym_atom_list;
	// Store the basis set name
	std::string basis_set_name = "sto-3g";
	// Store the k points
	vec2 k_pt;
	// The XCW log file
	std::ofstream XCW_log;
	// Store the number of reflections in the fit set
	int nr_fit;
	// Store the I/sigma cutoff 
	double i_sigma_cutoff = 2.0;
	// Store the inverse scale factor for GooF calculations
	double inv_scale;
	// Extincition model
	extinction::model ext_model = extinction::model::none;
	// Store the wavelength
	double wavelength;
	// Stuff regarding extinction refinement
	bool extinction_aniso = false;
	bool extinction_refine = true;
	double extinction_start = 1e-4;
	vec ext_p_;
	vec ext_c_, ext_cos2t_, ext_a_;
	vec ext_y_, ext_sqrt_y_, ext_g_, ext_m_, ext_dyc_;
	// Ordered vector of hkl
	std::vector<i3> hkl_ordered_;


	int n_params() const {
		return model_data.n_params + static_cast<int>(extinction_refine ? ext_p_.size() : 0);
	}
	void setup_extinction(const std::filesystem::path& cif);
	void ensure_hkl_ordered();
	double read_cif_wavelength(const std::filesystem::path& cif);



	// closing class
};