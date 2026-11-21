#pragma once
#include "convenience.h"
#include "constants.h"
#include "wfn_class.h"
#include "featomic.hpp"
#include "metatensor.h"

class SALTEDDescriptors;

namespace SALTED_Utils
{
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
