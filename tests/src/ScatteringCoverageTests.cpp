#include "pch.h"
#include <gtest/gtest.h>
#include "core/convenience.h"
#include "core/constants.h"
#include "core/scattering_factors.h"
#include "core/tsc_block.h"
#include "core/cell.h"
#include "core/wfn_class.h"
#include "core/atoms.h"
#include "core/cube.h"
#include "core/spherical_density.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <cmath>
#include <chrono>
#include <complex>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

//convert_to_ED lives at namespace scope in scattering_factors.cpp but is not declared in the header
void convert_to_ED(const ivec& asym_atom_list, const WFN& wave, cvec2& sf, const vec& stl);
void convert_to_ED(const ivec& asym_atom_list, const WFN& wave, cvec2& sf, const cell& unit_cell, const std::vector<i3>& hkl_vector);
void convert_to_ED(const ivec& asym_atom_list, const WFN& wave, cvec2& sf, const cell& unit_cell, const hkl_list& hkl);

namespace
{
	namespace fs = std::filesystem;

	//cubic P1 cell of this edge, so stl = |hkl| / (2 kA) and k = 2 pi |hkl| / kA (bohr)
	constexpr double kA = 4.0;

	fs::path tmp(const std::string& name)
	{
		return fs::temp_directory_path() / ("nosphera2_sfcov_" + name);
	}

	//a scratch directory that is the working directory for the test's lifetime, removed afterwards
	struct scoped_cwd
	{
		fs::path old;
		fs::path dir;
		explicit scoped_cwd(const std::string& name) : old(fs::current_path()), dir(tmp(name))
		{
			fs::remove_all(dir);
			fs::create_directories(dir);
			fs::current_path(dir);
		}
		~scoped_cwd()
		{
			fs::current_path(old);
			fs::remove_all(dir);
		}
	};

	void write_p1_cif(const fs::path& cif, const std::string& atom_rows = "", const std::string& extra_headers = "")
	{
		std::ofstream out(cif);
		out << "data_test\n_cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n"
			<< "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n_cell_volume 64.0\n"
			<< "loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n";
		if (!atom_rows.empty())
			out << "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
				<< extra_headers << atom_rows;
	}

	//number printed right after a log tag
	long parse_after(const std::string& log, const std::string& tag)
	{
		const std::size_t at = log.find(tag);
		if (at == std::string::npos)
			return -1;
		return std::stol(log.substr(at + tag.size()));
	}

	int first_nonzero(const d3& v)
	{
		for (int i = 0; i < 3; i++)
			if (v[i] != 0.0)
				return v[i] < 0.0 ? -1 : 1;
		return 0;
	}

	//the reflections generate_fractional_hkl must end up with on the cubic kA cell for half-integer
	//steps: the sphere d >= 0.9 dmin is |hkl|^2 <= (kA / (0.9 dmin))^2, the generator only walks l >= 0
	//(the twin law {1,0,0,0,1,0,0,0,-1} closes that to the full sphere) and of each Friedel pair present
	//the lexicographically smaller one, the one whose first non-zero index is negative, is kept
	struct fractional_sphere
	{
		hkl_list_d candidates;
		hkl_list_d expected;
		fractional_sphere(const double dmin, const bool full_sphere)
		{
			const double r2max = std::pow(kA / (0.9 * dmin), 2);
			for (int H = -12; H <= 12; H++)
				for (int K = -12; K <= 12; K++)
					for (int L = -12; L <= 12; L++)
					{
						if (0.25 * (H * H + K * K + L * L) > r2max)
							continue;
						if (!full_sphere && L < 0)
							continue;
						candidates.insert(d3{ 0.5 * H, 0.5 * K, 0.5 * L });
					}
			for (const d3& v : candidates)
			{
				if (first_nonzero(v) == 0)
					continue;
				const d3 minus{ -v[0], -v[1], -v[2] };
				if (candidates.find(minus) == candidates.end() || first_nonzero(v) < 0)
					expected.insert(v);
			}
		}
	};

	double stl_cubic(const i3& hkl)
	{
		return std::sqrt(double(hkl[0] * hkl[0] + hkl[1] * hkl[1] + hkl[2] * hkl[2])) / (2.0 * kA);
	}

	//what the IAM path must put in a row: the Thakkar factor of the element, or its Mott-Bethe form
	cdouble iam_expected(const i3& hkl, const int Z, const bool ED)
	{
		const double stl = stl_cubic(hkl);
		const double f = Thakkar(Z).get_form_factor(constants::bohr2ang(constants::FOUR_PI * stl));
		if (ED)
			return cdouble(constants::ED_fact * (Z - f) / (stl * stl), 0.0);
		return cdouble(f, 0.0);
	}

	//C1 at (0.25, 0.25, 0.25) and O1 at (0.5, 0.5, 0.5) in the CIF, only the carbon in the wavefunction
	struct iam_setup
	{
		fs::path cif;
		std::vector<WFN> wavy;
		options opt;
		svec known;
		explicit iam_setup(const std::string& name) : cif(tmp(name + ".cif"))
		{
			write_p1_cif(cif, "C1 C 0.25 0.25 0.25\nO1 O 0.5 0.5 0.5\n");
			wavy.emplace_back(e_origin::NOT_YET_DEFINED);
			const double f = constants::ang2bohr(kA);
			wavy[0].push_back_atom("C1", 0.25 * f, 0.25 * f, 0.25 * f, 6);
			opt.cif = cif;
			opt.iam_switch = true;
			opt.groups[0] = ivec{ 0 };
			opt.m_hkl_list = { i3{ 1, 0, 0 }, i3{ 0, 2, 0 }, i3{ 1, 1, 1 } };
			opt.tsc_block_size = 0;
			opt.no_date = true;
			opt.use_gpu = false;
		}
		~iam_setup()
		{
			fs::remove(cif);
		}
		itsc_block run(std::ostream& log, vec2* kpts = nullptr, salted_part_prep* prep = nullptr)
		{
			return calculate_scattering_factors<itsc_block, std::vector<WFN>&>(opt, wavy, log, known, 0, kpts, prep);
		}
	};

	void expect_iam_rows(const itsc_block& block, const std::size_t scatterer, const int Z, const bool ED)
	{
		const cvec& row = block.get_sf_for_scatterer(scatterer);
		for (std::size_t r = 0; r < block.reflection_size(); r++)
		{
			const cdouble e = iam_expected(block.get_indices(r), Z, ED);
			EXPECT_NEAR(row[r].real(), e.real(), 1e-10) << "reflection " << r;
			EXPECT_NEAR(row[r].imag(), e.imag(), 1e-12) << "reflection " << r;
		}
	}

	//a normalised s Gaussian rho = pi^-3/2 exp(-r^2) on one hydrogen, F(k) = exp(-k^2 / 4)
	void write_gaussian_h_wfn(const fs::path& path)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_MO(0, 1.0, -0.5);
		w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
		double coef = std::pow(constants::PI, -0.75);
		w.add_primitive(1, 1, 0.5, &coef);
		ASSERT_TRUE(w.write_wfn(path, false, false));
	}

	vec2 read_tsc_rows(const fs::path& path)
	{
		std::ifstream in(path);
		std::string line;
		while (std::getline(in, line) && line != "DATA:")
		{
		}
		vec2 rows;
		while (std::getline(in, line))
		{
			if (line.empty())
				continue;
			for (char& c : line)
				if (c == ',')
					c = ' ';
			std::istringstream ss(line);
			vec v(5);
			ss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4];
			rows.push_back(v);
		}
		return rows;
	}
}

//generate_fractional_hkl: a half-integer sphere, walked for l >= 0 and Friedel-reduced, origin removed
TEST(ScatteringCoverageHklTests, FractionalHklSphereHalfSteps)
{
	cell c(kA, kA, kA, 90.0, 90.0, 90.0);
	const double dmin = 2.0;
	const fractional_sphere ref(dmin, false);
	hkl_list_d hkl;
	std::ostringstream log;
	generate_fractional_hkl(dmin, hkl, {}, c, log, d3{ 0.5, 0.5, 0.5 }, false);
	EXPECT_EQ(hkl, ref.expected);
	EXPECT_EQ(parse_after(log.str(), "Nr of reflections generated: "), (long)ref.candidates.size());
	EXPECT_EQ(parse_after(log.str(), "Nr of reflections to be used: "), (long)ref.expected.size());
	EXPECT_EQ(parse_after(log.str(), "Number of symmetry operations: "), 1);
	EXPECT_NE(log.str().find("Generating hkl indices up to d=:"), std::string::npos);
	EXPECT_EQ(hkl.count(d3{ 0.0, 0.0, 0.0 }), 0u);
}

//the twin law z -> -z closes the l >= 0 walk to the full sphere, the debug branch reports each stage
TEST(ScatteringCoverageHklTests, FractionalHklTwinLawAndDebugLog)
{
	cell c(kA, kA, kA, 90.0, 90.0, 90.0);
	const double dmin = 3.0;
	const fractional_sphere half(dmin, false);
	const fractional_sphere full(dmin, true);
	hkl_list_d hkl;
	std::ostringstream log;
	const vec2 twin{ vec{ 1, 0, 0, 0, 1, 0, 0, 0, -1 } };
	generate_fractional_hkl(dmin, hkl, twin, c, log, d3{ 0.5, 0.5, 0.5 }, true);
	EXPECT_EQ(hkl, full.expected);
	EXPECT_EQ(hkl.size(), (full.candidates.size() - 1) / 2);
	EXPECT_EQ(parse_after(log.str(), "Number of reflections before twin: "), (long)half.candidates.size());
	EXPECT_EQ(parse_after(log.str(), "Number of reflections after twin: "), (long)full.candidates.size());
	EXPECT_EQ(parse_after(log.str(), "Number of reflections after sym gen: "), (long)full.candidates.size());
	EXPECT_NE(log.str().find("Read 1 symmetry elements!"), std::string::npos);
	EXPECT_EQ(log.str().find("Number of symmetry operations: "), std::string::npos);
}

//read_hkl_full: text and one-character lines are skipped, F and sigma(F) follow from F^2 and sigma(F^2)
TEST(ScatteringCoverageHklTests, ReadHklFullDerivesFAndSigma)
{
	const fs::path path = tmp("full.hkl");
	{
		std::ofstream out(path);
		out << "CELL 0.71073 4 4 4 90 90 90\n"
			<< "x\n"
			<< "   1   0   0  100.00   10.00\n"
			<< "   0   0   1  400.00   20.00\n";
	}
	cell c(kA, kA, kA, 90.0, 90.0, 90.0);
	hkl_list hkl;
	std::vector<scattering_data> obs;
	std::ostringstream log;
	const hkl_list enlarged = read_hkl_full(path, hkl, {}, c, log, obs, false);
	fs::remove(path);
	EXPECT_EQ(hkl, (hkl_list{ i3{ 0, 0, 1 }, i3{ 1, 0, 0 } }));
	EXPECT_EQ(enlarged, hkl);
	//obs follows the (h,k,l) order of the set, not the file: (0,0,1) comes first
	ASSERT_EQ(obs.size(), 2u);
	EXPECT_DOUBLE_EQ(obs[1].F_obs2, 100.0);
	EXPECT_DOUBLE_EQ(obs[1].F_obs, 10.0);
	EXPECT_DOUBLE_EQ(obs[1].abs_F_obs, 10.0);
	EXPECT_DOUBLE_EQ(obs[1].sigma_obs2, 10.0);
	EXPECT_DOUBLE_EQ(obs[1].sigma_obs, 0.5);
	EXPECT_DOUBLE_EQ(obs[0].F_obs, 20.0);
	EXPECT_DOUBLE_EQ(obs[0].sigma_obs, 0.5);
	EXPECT_EQ(parse_after(log.str(), "Nr of reflections read from file: "), 2);
	EXPECT_EQ(parse_after(log.str(), "Number of symmetry operations: "), 1);
	EXPECT_EQ(parse_after(log.str(), "Nr of reflections to be used: "), 2);
}

//a row too short for h k l, a row without an F^2 value and a file holding only 0 0 0 each stop the run
TEST(ScatteringCoverageHklTests, ReadHklFullRejectsBrokenRows)
{
	const fs::path short_row = tmp("short.hkl");
	const fs::path no_dot = tmp("nodot.hkl");
	const fs::path only_origin = tmp("origin.hkl");
	{
		std::ofstream(short_row) << "   1   0\n";
		std::ofstream(no_dot) << "   1   0   0    12345   6789\n";
		std::ofstream(only_origin) << "   0   0   0  100.00   10.00\n";
	}
	cell c(kA, kA, kA, 90.0, 90.0, 90.0);
	hkl_list hkl;
	std::vector<scattering_data> obs;
	std::ostringstream log;
	EXPECT_EXIT({ read_hkl_full(short_row, hkl, {}, c, log, obs, false); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT({ read_hkl_full(no_dot, hkl, {}, c, log, obs, false); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT({ read_hkl_full(only_origin, hkl, {}, c, log, obs, true); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	fs::remove(short_row);
	fs::remove(no_dot);
	fs::remove(only_origin);
}

//a disorder group that is not a number is reported and read as 0
TEST(ScatteringCoverageCifTests, ReadAtomsUnreadableDisorderGroupIsZero)
{
	const fs::path cif = tmp("group_a.cif");
	write_p1_cif(cif, "C1 C 0.25 0.35 0.45 A\n", "_atom_site_disorder_group\n");
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	WFN w(e_origin::NOT_YET_DEFINED);
	const double f = constants::ang2bohr(kA);
	w.push_back_atom("C1", 0.25 * f, 0.35 * f, 0.45 * f, 6);
	std::ifstream in(cif);
	ivec types, asym_to_type, asym;
	bvec needs_grid(1, false);
	std::ostringstream log;
	const svec labels = read_atoms_from_CIF(in, {}, c, w, {}, types, asym_to_type, asym, needs_grid, log, false);
	in.close();
	fs::remove(cif);
	EXPECT_NE(log.str().find("Could not read disorder group for atom C1; treating it as 0."), std::string::npos);
	EXPECT_EQ(labels, (svec{ "C1" }));
	EXPECT_EQ(w.get_id_for_atom(0).to_hex_string(), atomID(0.25, 0.35, 0.45, 0, 6).to_hex_string());
}

//a row whose label carries no element and whose type is another element is dropped, with and without debug
TEST(ScatteringCoverageCifTests, ReadAtomsRejectsRowByLabelAndType)
{
	const fs::path cif = tmp("q1n.cif");
	write_p1_cif(cif, "Q1 N 0.25 0.35 0.45\n");
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	const double f = constants::ang2bohr(kA);
	for (const bool debug : { true, false })
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C1", 0.25 * f, 0.35 * f, 0.45 * f, 6);
		std::ifstream in(cif);
		ivec types, asym_to_type, asym;
		bvec needs_grid(1, false);
		std::ostringstream log;
		const svec labels = read_atoms_from_CIF(in, {}, c, w, {}, types, asym_to_type, asym, needs_grid, log, debug, true);
		in.close();
		EXPECT_TRUE(labels.empty()) << debug;
		EXPECT_TRUE(asym.empty()) << debug;
		EXPECT_EQ(types, (ivec{ 6 })) << debug;
		if (debug)
		{
			EXPECT_NE(log.str().find("Element symbol not found in label, this is a problem!\n checking type... ALSO FAILED! WILL IGNORE ATOM!"), std::string::npos);
			EXPECT_NE(log.str().find("I did not find this atom! Tolerances were: "), std::string::npos);
		}
		else
			EXPECT_NE(log.str().find("Atom q1 was not matching by element determined by label reduction or type field, skipping!"), std::string::npos);
	}
	//without permission for an empty asymmetric unit the same CIF ends the run
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C1", 0.25 * f, 0.35 * f, 0.45 * f, 6);
		std::ifstream in(cif);
		ivec types, asym_to_type, asym;
		bvec needs_grid(1, false);
		std::ostringstream log;
		EXPECT_EXIT({ read_atoms_from_CIF(in, {}, c, w, {}, types, asym_to_type, asym, needs_grid, log, false, false); },
			::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	}
	fs::remove(cif);
}

//hydrogen rows: a label without h is accepted when the type says H, a label with t or d is accepted
//unconditionally, and a carbon type on a hydrogen site is rejected
TEST(ScatteringCoverageCifTests, ReadAtomsHydrogenLabelRules)
{
	struct probe
	{
		std::string row;
		bool accepted;
		std::string message;
	};
	const std::vector<probe> probes{
		{ "Q1 H 0.25 0.35 0.45\n", true, " will check type..." },
		{ "T1 H 0.25 0.35 0.45\n", true, "" },
		{ "Q1 C 0.25 0.35 0.45\n", false, " will check type... ALSO FAILED! WILL IGNORE ATOM!" } };
	const double f = constants::ang2bohr(kA);
	for (const probe& p : probes)
	{
		const fs::path cif = tmp("hrule.cif");
		write_p1_cif(cif, p.row);
		std::ostringstream sink;
		cell c(cif, sink, false, true);
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("H1", 0.25 * f, 0.35 * f, 0.45 * f, 1);
		std::ifstream in(cif);
		ivec types, asym_to_type, asym;
		bvec needs_grid(1, false);
		std::ostringstream log;
		const svec labels = read_atoms_from_CIF(in, {}, c, w, {}, types, asym_to_type, asym, needs_grid, log, true, true);
		in.close();
		fs::remove(cif);
		if (p.accepted)
		{
			EXPECT_EQ(labels, (svec{ p.row.substr(0, 2) })) << p.row;
			EXPECT_EQ(asym, (ivec{ 0 })) << p.row;
		}
		else
		{
			EXPECT_TRUE(labels.empty()) << p.row;
		}
		if (!p.message.empty())
			EXPECT_NE(log.str().find(p.message), std::string::npos) << p.row;
		else
			EXPECT_EQ(log.str().find("check type"), std::string::npos) << p.row;
		EXPECT_EQ(types, (ivec{ 1 })) << p.row;
	}
}

//convert_to_ED is Mott-Bethe on a complex f: (Z - Re f) and -Im f, both over stl^2
TEST(ScatteringCoverageEdTests, ConvertToEdIsMottBethe)
{
	WFN w(e_origin::NOT_YET_DEFINED);
	w.push_back_atom("C1", 0.0, 0.0, 0.0, 6);
	w.push_back_atom("O1", 1.0, 0.0, 0.0, 8);
	const ivec asym{ 0, 1 };
	const double E = constants::ED_fact;
	{
		cvec2 sf{ cvec{ cdouble(2.0, 1.0), cdouble(1.0, 0.0) }, cvec{ cdouble(5.0, -2.0), cdouble(0.0, 0.0) } };
		convert_to_ED(asym, w, sf, vec{ 0.25, 0.5 });
		EXPECT_NEAR(sf[0][0].real(), E * 4.0 / 0.0625, 1e-12);
		EXPECT_NEAR(sf[0][0].imag(), -E * 1.0 / 0.0625, 1e-12);
		EXPECT_NEAR(sf[0][1].real(), E * 5.0 / 0.25, 1e-12);
		EXPECT_NEAR(sf[0][1].imag(), 0.0, 1e-12);
		EXPECT_NEAR(sf[1][0].real(), E * 3.0 / 0.0625, 1e-12);
		EXPECT_NEAR(sf[1][0].imag(), E * 2.0 / 0.0625, 1e-12);
		EXPECT_NEAR(sf[1][1].real(), E * 8.0 / 0.25, 1e-12);
	}
	//the two cell overloads derive stl from the indices: (2,0,0) is stl 0.25, (0,0,4) is 0.5 on the 4 A cell
	cell c(kA, kA, kA, 90.0, 90.0, 90.0);
	{
		cvec2 sf{ cvec{ cdouble(2.0, 1.0), cdouble(1.0, 0.0) } };
		convert_to_ED(ivec{ 0 }, w, sf, c, std::vector<i3>{ i3{ 2, 0, 0 }, i3{ 0, 0, 4 } });
		EXPECT_NEAR(sf[0][0].real(), E * 4.0 / 0.0625, 1e-9);
		EXPECT_NEAR(sf[0][0].imag(), -E * 1.0 / 0.0625, 1e-9);
		EXPECT_NEAR(sf[0][1].real(), E * 5.0 / 0.25, 1e-9);
	}
	{
		//a set orders (0,0,4) before (2,0,0)
		cvec2 sf{ cvec{ cdouble(1.0, 0.0), cdouble(2.0, 1.0) } };
		convert_to_ED(ivec{ 0 }, w, sf, c, hkl_list{ i3{ 2, 0, 0 }, i3{ 0, 0, 4 } });
		EXPECT_NEAR(sf[0][0].real(), E * 5.0 / 0.25, 1e-9);
		EXPECT_NEAR(sf[0][1].real(), E * 4.0 / 0.0625, 1e-9);
		EXPECT_NEAR(sf[0][1].imag(), -E * 1.0 / 0.0625, 1e-9);
	}
}

//the IAM path of calculate_scattering_factors: one Thakkar row per asymmetric atom, keyed by atomID
TEST(ScatteringCoverageIamTests, IamRowsAreThakkarFactors)
{
	iam_setup s("iam");
	std::ostringstream log;
	const itsc_block block = s.run(log);
	EXPECT_NE(log.str().find("Number of protons: 6"), std::string::npos);
	ASSERT_EQ(block.scatterer_size(), 1u);
	ASSERT_EQ(block.reflection_size(), 3u);
	EXPECT_EQ(std::get<atomID>(block.get_scatterer(0)).to_hex_string(), atomID(0.25, 0.25, 0.25, 0, 6).to_hex_string());
	EXPECT_EQ(block.get_index_vector(), (ivec2{ { 0, 1, 1 }, { 2, 0, 1 }, { 0, 0, 1 } }));
	expect_iam_rows(block, 0, 6, false);
	EXPECT_FALSE(s.opt.tsc_written_by_stream);
}

//electron diffraction in the IAM path is Mott-Bethe with the tabulated charge, and -label puts labels on the rows
TEST(ScatteringCoverageIamTests, IamElectronDiffractionAndLabels)
{
	iam_setup s("iam_ed");
	s.opt.electron_diffraction = true;
	s.opt.label_tsc_output = true;
	std::ostringstream log;
	const itsc_block block = s.run(log);
	ASSERT_EQ(block.scatterer_size(), 1u);
	EXPECT_EQ(block.get_scatterers_string(), (svec{ "C1" }));
	expect_iam_rows(block, 0, 6, true);
}

//with prep_out the IAM path stops after the bookkeeping and hands back the sorted reflections, stl and k
TEST(ScatteringCoverageIamTests, IamPrepOutHandsBackReflectionTables)
{
	iam_setup s("iam_prep");
	salted_part_prep prep;
	std::ostringstream log;
	const itsc_block block = s.run(log, nullptr, &prep);
	EXPECT_TRUE(block.is_empty());
	EXPECT_EQ(prep.hkl_v, (std::vector<i3>{ i3{ 0, 2, 0 }, i3{ 1, 0, 0 }, i3{ 1, 1, 1 } }));
	EXPECT_EQ(prep.labels, (svec{ "C1" }));
	EXPECT_EQ(prep.asym_atom_list, (ivec{ 0 }));
	EXPECT_EQ(prep.atom_type_list, (ivec{ 6 }));
	EXPECT_EQ(prep.asym_atom_to_type_list, (ivec{ 0 }));
	ASSERT_EQ(prep.stl_of_reflection.size(), 3u);
	ASSERT_EQ(prep.k_of_reflection.size(), 3u);
	for (std::size_t r = 0; r < 3; r++)
	{
		const double stl = stl_cubic(prep.hkl_v[r]);
		EXPECT_NEAR(prep.stl_of_reflection[r], stl, 1e-12);
		EXPECT_NEAR(prep.k_of_reflection[r], constants::bohr2ang(constants::FOUR_PI * stl), 1e-12);
	}
	EXPECT_FALSE(s.opt.tsc_written_by_stream);
}

//an empty k-point vector is filled with 2 pi h / a (bohr^-1) and a filled one is taken as is
TEST(ScatteringCoverageIamTests, IamKpointsOutAndIn)
{
	iam_setup s("iam_kpts");
	vec2 kp;
	std::ostringstream log;
	const itsc_block first = s.run(log, &kp);
	ASSERT_EQ(kp.size(), 3u);
	ASSERT_EQ(kp[0].size(), 3u);
	const double g = constants::TWO_PI / constants::ang2bohr(kA);
	for (std::size_t r = 0; r < 3; r++)
	{
		const i3 hkl = first.get_indices(r);
		for (int x = 0; x < 3; x++)
			EXPECT_NEAR(kp[x][r], g * hkl[x], 1e-9) << r << " " << x;
	}
	EXPECT_NE(log.str().find("Number of k-points to evaluate: 3"), std::string::npos);
	std::ostringstream log2;
	const itsc_block second = s.run(log2, &kp);
	EXPECT_EQ(log2.str().find("Generating k-points"), std::string::npos);
	ASSERT_EQ(second.reflection_size(), 3u);
	expect_iam_rows(second, 0, 6, false);
}

//with a block size the IAM rows are streamed into experimental.tscb in the working directory
TEST(ScatteringCoverageIamTests, IamStreamsBlocksToTscb)
{
	scoped_cwd cwd("iam_stream");
	iam_setup s("iam_stream");
	s.opt.tsc_block_size = 2;
	std::ostringstream log;
	const itsc_block block = s.run(log);
	EXPECT_NE(log.str().find("Streaming tsc in blocks of 2 reflections"), std::string::npos);
	EXPECT_TRUE(s.opt.tsc_written_by_stream);
	ASSERT_TRUE(fs::exists(cwd.dir / "experimental.tscb"));
	const itsc_block back(cwd.dir / "experimental.tscb");
	ASSERT_EQ(back.scatterer_size(), 1u);
	ASSERT_EQ(back.reflection_size(), 3u);
	EXPECT_EQ(std::get<atomID>(back.get_scatterer(0)).to_hex_string(), atomID(0.25, 0.25, 0.25, 0, 6).to_hex_string());
	expect_iam_rows(back, 0, 6, false);
	//the streamed table is on disk, the returned block only names the scatterer: no rows, so the
	//size accessors report 0 for the mismatched table
	EXPECT_EQ(block.get_scatterers().size(), 1u);
	EXPECT_EQ(block.scatterer_size(), 0u);
	EXPECT_EQ(block.reflection_size(), 0u);
}

//the stream writer is binary only: a text table (no tscb, or -old_tsc) with a block size is kept
//in memory and returned with the ED rows, and nothing is written under experimental.tsc*
TEST(ScatteringCoverageIamTests, IamTextTableIsNotStreamed)
{
	scoped_cwd cwd("iam_stream_txt");
	for (const bool old_tsc : { false, true })
	{
		iam_setup s("iam_stream_txt");
		s.opt.tsc_block_size = 1000;
		s.opt.binary_tsc = old_tsc;
		s.opt.old_tsc = old_tsc;
		s.opt.electron_diffraction = true;
		std::ostringstream log;
		const itsc_block block = s.run(log);
		EXPECT_EQ(log.str().find("Streaming tsc in blocks"), std::string::npos) << old_tsc;
		EXPECT_NE(log.str().find("Text tsc requested: the table is kept in memory"), std::string::npos) << old_tsc;
		EXPECT_FALSE(s.opt.tsc_written_by_stream) << old_tsc;
		EXPECT_FALSE(fs::exists(cwd.dir / "experimental.tsc")) << old_tsc;
		EXPECT_FALSE(fs::exists(cwd.dir / "experimental.tscb")) << old_tsc;
		ASSERT_EQ(block.scatterer_size(), 1u) << old_tsc;
		ASSERT_EQ(block.reflection_size(), 3u) << old_tsc;
		expect_iam_rows(block, 0, 6, true);
	}
}

//-cif_based -mtc reads the part's own CIF from the combined list instead of -cif
TEST(ScatteringCoverageIamTests, IamCifBasedCombinedUsesPartCif)
{
	iam_setup s("iam_cifbased");
	s.opt.cif_based_combined_tsc_calc = true;
	s.opt.combined_tsc_calc_cifs = { s.cif };
	s.opt.cif = "";
	std::ostringstream log;
	const itsc_block block = s.run(log);
	ASSERT_EQ(block.scatterer_size(), 1u);
	expect_iam_rows(block, 0, 6, false);
}

//an already known atom leaves nothing to compute: an empty block, allowed only with allow_empty_asym
TEST(ScatteringCoverageIamTests, IamKnownAtomsGiveEmptyBlock)
{
	iam_setup s("iam_known");
	s.opt.allow_empty_asym = true;
	s.known = { "C1" };
	salted_part_prep prep;
	prep.labels = { "stale" };
	std::ostringstream log;
	const itsc_block block = s.run(log, nullptr, &prep);
	EXPECT_TRUE(block.is_empty());
	EXPECT_TRUE(prep.labels.empty());
	EXPECT_TRUE(prep.hkl_v.empty());
}

//needs_Thakkar_fill: the atoms the wavefunction lacks come from -wfn as Thakkar rows on the same reflections
TEST(ScatteringCoverageIamTests, ThakkarFillAppendsMissingAtomsFromXyz)
{
	iam_setup s("iam_fill");
	const fs::path xyz = tmp("iam_fill.xyz");
	std::ofstream(xyz) << "2\nfill\nC 1.0 1.0 1.0\nO 2.0 2.0 2.0\n";
	s.opt.wfn = xyz;
	s.opt.needs_Thakkar_fill = true;
	s.opt.label_tsc_output = true;
	std::ostringstream log;
	const itsc_block block = s.run(log);
	fs::remove(xyz);
	EXPECT_NE(log.str().find("Performing the remaining calculation of spherical atoms..."), std::string::npos);
	EXPECT_NE(log.str().find("Number of protons: 14"), std::string::npos);
	ASSERT_EQ(block.scatterer_size(), 2u);
	ASSERT_EQ(block.reflection_size(), 3u);
	EXPECT_EQ(block.get_scatterers_string(), (svec{ "C1", "O1" }));
	expect_iam_rows(block, 0, 6, false);
	expect_iam_rows(block, 1, 8, false);
	EXPECT_FALSE(s.opt.needs_Thakkar_fill);
	EXPECT_TRUE(s.opt.iam_switch);
	EXPECT_FALSE(s.opt.allow_empty_asym);
	EXPECT_FALSE(s.opt.spherical_fill);
	EXPECT_EQ(s.opt.m_hkl_list, (hkl_list{ i3{ 1, 0, 0 }, i3{ 0, 2, 0 }, i3{ 1, 1, 1 } }));
}

//stream_mtc_salted declines without a block size, without -SALTED, with fewer than two parts or under -IAM
TEST(ScatteringCoverageStreamTests, StreamMtcSaltedGates)
{
	scoped_cwd cwd("mtc_gates");
	std::vector<WFN> none;
	std::ostringstream log;
	options opt;
	opt.SALTED = true;
	opt.combined_tsc_calc_files = { "a.xyz", "b.xyz" };
	opt.iam_switch = true;
	EXPECT_FALSE(stream_mtc_salted(opt, none, log, nullptr));
	opt.iam_switch = false;
	opt.tsc_block_size = 0;
	EXPECT_FALSE(stream_mtc_salted(opt, none, log, nullptr));
	opt.tsc_block_size = 1000;
	opt.SALTED = false;
	EXPECT_FALSE(stream_mtc_salted(opt, none, log, nullptr));
	opt.SALTED = true;
	opt.combined_tsc_calc_files = { "a.xyz" };
	EXPECT_FALSE(stream_mtc_salted(opt, none, log, nullptr));
	EXPECT_TRUE(log.str().empty());
	EXPECT_TRUE(fs::is_empty(cwd.dir));
	EXPECT_FALSE(opt.tsc_written_by_stream);
}

//calc_SF on the CPU: sum of rho exp(i k.d) over the points, the debug and timing lines around it
TEST(ScatteringCoverageTransformTests, CalcSfSumsPhasesAndReportsPreparationTime)
{
	const double pi = constants::PI;
	const vec2 k_pt{ vec{ pi, 2.0 * pi }, vec{ 0.0, 0.0 }, vec{ 0.0, 0.0 } };
	const vec2 dx{ vec{ 0.5 }, vec{ 0.5, -0.5 } };
	const vec2 dy{ vec{ 0.0 }, vec{ 0.0, 0.0 } };
	const vec2 dz{ vec{ 0.0 }, vec{ 0.0, 0.0 } };
	const vec2 dens{ vec{ 1.0 }, vec{ 0.5, 0.5 } };
	cvec2 sf;
	std::ostringstream log;
	_time_point start = get_time() - std::chrono::seconds(2);
	_time_point end1;
	calc_SF(3, k_pt, dx, dy, dz, dens, sf, log, start, end1, true, false, false, false);
	ASSERT_EQ(sf.size(), 2u);
	ASSERT_EQ(sf[0].size(), 2u);
	//atom 0: one point at x = 0.5 -> exp(i pi/2) = i and exp(i pi) = -1
	EXPECT_NEAR(sf[0][0].real(), 0.0, 1e-12);
	EXPECT_NEAR(sf[0][0].imag(), 1.0, 1e-12);
	EXPECT_NEAR(sf[0][1].real(), -1.0, 1e-12);
	EXPECT_NEAR(sf[0][1].imag(), 0.0, 1e-12);
	//atom 1: half weight at +-0.5 -> cos(pi/2) = 0 and cos(pi) = -1, no imaginary part
	EXPECT_NEAR(sf[1][0].real(), 0.0, 1e-12);
	EXPECT_NEAR(sf[1][0].imag(), 0.0, 1e-12);
	EXPECT_NEAR(sf[1][1].real(), -1.0, 1e-12);
	EXPECT_NEAR(sf[1][1].imag(), 0.0, 1e-12);
	EXPECT_NE(log.str().find("Initialized FFs"), std::string::npos);
	EXPECT_NE(log.str().find("asym atom list size: 2 total grid size: 3"), std::string::npos);
	EXPECT_NE(log.str().find("Time to prepare: 2 s"), std::string::npos);
	EXPECT_GE(end1, start);
	//under a second the preparation time is reported in milliseconds
	std::ostringstream log2;
	cvec2 sf2;
	_time_point now = get_time();
	calc_SF(3, k_pt, dx, dy, dz, dens, sf2, log2, now, end1, false, false, false, false);
	EXPECT_NE(log2.str().find(" ms"), std::string::npos);
	EXPECT_EQ(log2.str().find("Initialized FFs"), std::string::npos);
	EXPECT_NEAR(sf2[0][0].imag(), 1.0, 1e-12);
}

//calc_sfac_diffuse: a normalised Gaussian on one hydrogen at the origin gives F = exp(-k^2/4), real,
//on the half-integer reflection set, written as a non-integer tsc into the working directory
TEST(ScatteringCoverageDiffuseTests, GaussianDiffuseTableMatchesAnalyticTransform)
{
	scoped_cwd cwd("diffuse");
	const fs::path wfn = cwd.dir / "gauss.wfn";
	const fs::path cif = cwd.dir / "diffuse.cif";
	write_gaussian_h_wfn(wfn);
	write_p1_cif(cif, "H1 H 0.0 0.0 0.0\n");
	options opt;
	opt.wfn = wfn;
	opt.cif = cif;
	opt.dmin = 4.0;
	opt.sfac_diffuse = d3{ 0.5, 0.5, 0.5 };
	opt.partition_type = PartitionType::Becke;
	opt.use_gpu = false;
	opt.no_date = true;
	std::ostringstream log;
	calc_sfac_diffuse(opt, log);
	ASSERT_TRUE(fs::exists(cwd.dir / "experimental.tsc"));
	const fractional_sphere ref(opt.dmin, false);
	EXPECT_EQ(parse_after(log.str(), "Nr of reflections to be used: "), (long)ref.expected.size());
	const vec2 rows = read_tsc_rows(cwd.dir / "experimental.tsc");
	ASSERT_EQ(rows.size(), ref.expected.size());
	hkl_list_d seen;
	const double g = constants::TWO_PI / constants::ang2bohr(kA);
	for (const vec& row : rows)
	{
		seen.insert(d3{ row[0], row[1], row[2] });
		const double k2 = g * g * (row[0] * row[0] + row[1] * row[1] + row[2] * row[2]);
		EXPECT_NEAR(row[3], std::exp(-k2 / 4.0), 2e-3) << row[0] << " " << row[1] << " " << row[2];
		EXPECT_NEAR(row[4], 0.0, 2e-3) << row[0] << " " << row[1] << " " << row[2];
	}
	EXPECT_EQ(seen, ref.expected);
	std::ifstream in(cwd.dir / "experimental.tsc");
	std::string head;
	std::getline(in, head);
	EXPECT_EQ(head, "TITLE: diffuse");
	std::getline(in, head);
	std::getline(in, head);
	EXPECT_EQ(head, "SCATTERERS: H1");
}

//the cube path with the reflections generated from -dmin, Mott-Bethe conversion, labels and the timing table
TEST(ScatteringCoverageCubeTests, CubeElectronDiffractionWithGeneratedReflections)
{
	const fs::path cif = tmp("cube_ed.cif");
	write_p1_cif(cif, "H1 H 0.0 0.0 0.0\n");
	WFN w(e_origin::NOT_YET_DEFINED);
	w.push_back_MO(0, 1.0, -1.0);
	w.push_back_atom("H1", 0.0, 0.0, 0.0, 1);
	double coef = 1.0;
	w.add_primitive(1, 1, 1.0, &coef);
	const int n = 101;
	const double h = 0.2;
	const double origin = -10.0;
	cube density({ n, n, n }, 1, true);
	for (int k = 0; k < 3; k++)
	{
		density.set_origin(k, origin);
		density.set_vector(k, k, h);
	}
	const double norm = std::pow(constants::PI, -1.5);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
			{
				const double x = origin + i * h, y = origin + j * h, z = origin + k * h;
				density.set_value(i, j, k, norm * std::exp(-(x * x + y * y + z * z)));
			}
	options opt;
	opt.cif = cif;
	opt.groups[0] = ivec{ 0 };
	opt.partition_type = PartitionType::Becke;
	opt.use_gpu = false;
	opt.no_date = false;
	opt.electron_diffraction = true;
	opt.label_tsc_output = true;
	//ED callers get the sphere at dmin/2 - 0.001: 2.499 A on the 4 A cell is |hkl|^2 <= 2, nine Friedel-unique reflections
	opt.dmin = 5.0;
	std::ostringstream log;
	const itsc_block sf = calculate_scattering_factors_from_cube(opt, w, density, log);
	fs::remove(cif);
	EXPECT_NE(log.str().find("Time Breakdown!"), std::string::npos);
	ASSERT_EQ(sf.scatterer_size(), 1u);
	EXPECT_EQ(sf.get_scatterers_string(), (svec{ "H1" }));
	ASSERT_EQ(sf.reflection_size(), 9u);
	EXPECT_EQ(opt.m_hkl_list.size(), 9u);
	const double g = constants::TWO_PI / constants::ang2bohr(kA);
	const cvec& row = sf.get_sf_for_scatterer(0);
	for (std::size_t r = 0; r < sf.reflection_size(); r++)
	{
		const i3 hkl = sf.get_indices(r);
		const int n2 = hkl[0] * hkl[0] + hkl[1] * hkl[1] + hkl[2] * hkl[2];
		EXPECT_GE(n2, 1) << hkl[0] << hkl[1] << hkl[2];
		EXPECT_LE(n2, 2) << hkl[0] << hkl[1] << hkl[2];
		const double stl2 = n2 / (4.0 * kA * kA);
		const double f = std::exp(-g * g * n2 / 4.0);
		const double expected = constants::ED_fact * (1.0 - f) / stl2;
		EXPECT_NEAR(row[r].real(), expected, 0.03 * constants::ED_fact / stl2) << hkl[0] << hkl[1] << hkl[2];
		EXPECT_NEAR(row[r].imag(), 0.0, 0.03 * constants::ED_fact / stl2) << hkl[0] << hkl[1] << hkl[2];
	}
}

//a NaN structure factor ends the run when the block is built
TEST(ScatteringCoverageBlockTests, NanStructureFactorIsFatal)
{
	const cvec2 sf{ cvec{ cdouble(std::nan(""), 0.0) } };
	const hkl_list hkl{ i3{ 1, 0, 0 } };
	EXPECT_EXIT({ itsc_block b(sf, svec{ "C1" }, hkl); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//the member prologue is byte for byte the head of the file write_tscb_file produces
TEST(ScatteringCoverageBlockTests, MemberPrologueMatchesFileHead)
{
	scoped_cwd cwd("prologue");
	const cvec2 sf{ cvec{ cdouble(1.0, 2.0), cdouble(3.0, 4.0) }, cvec{ cdouble(5.0, 6.0), cdouble(7.0, 8.0) } };
	const hkl_list hkl{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 } };
	const itsc_block b(sf, svec{ "C1", "O1" }, hkl, "head");
	b.write_tscb_file("x.cif", "full.tscb");
	std::ostringstream prologue(std::ios::binary);
	b.write_tscb_prologue(prologue, 2);
	const std::string head = prologue.str();
	std::ifstream in(cwd.dir / "full.tscb", std::ios::binary);
	std::stringstream ss;
	ss << in.rdbuf();
	const std::string file = ss.str();
	ASSERT_GT(head.size(), 8u);
	ASSERT_GE(file.size(), head.size());
	EXPECT_EQ(file.substr(0, head.size()), head);
	//int32 header length, the header, then the scatterers and the int32 reflection count
	EXPECT_EQ(head.substr(4, 4), "head");
	EXPECT_EQ(head.substr(head.size() - 4), std::string("\x02\x00\x00\x00", 4));
	//what follows the prologue in the file is the reflection data: 2 x (3 int32 + 2 cdouble)
	EXPECT_EQ(file.size() - head.size(), 2u * (3 * sizeof(std::int32_t) + 2 * sizeof(cdouble)));
}

//the non-integer writer on an integer block: fixed/setprecision do not act on int, the indices stay integers
//(the double-indexed block with three decimals is covered by ScatteringFactorBlockTests.NonIntegerIndicesWriter)
TEST(ScatteringCoverageBlockTests, NonIntegerWriterOnIntegerBlock)
{
	scoped_cwd cwd("nonint");
	const cvec2 sf{ cvec{ cdouble(1.0, -0.5) } };
	const itsc_block b(sf, svec{ "C1" }, hkl_list{ i3{ 1, 0, -2 } });
	b.write_tsc_file_non_integer("nonint.cif", "nonint.tsc");
	std::ifstream in(cwd.dir / "nonint.tsc");
	std::stringstream ss;
	ss << in.rdbuf();
	const std::string text = ss.str();
	EXPECT_NE(text.find("TITLE: nonint\nSYMM: expanded\nSCATTERERS: C1\nDATA:\n"), std::string::npos);
	EXPECT_NE(text.find("1 0 -2 1.00000000e+00,-5.00000000e-01 \n"), std::string::npos);
}

//appending a block with another number of reflections is fatal, one with other indices too
TEST(ScatteringCoverageBlockTests, AppendRejectsIncompatibleReflections)
{
	const cvec2 one{ cvec{ cdouble(1.0, 0.0) } };
	const cvec2 two{ cvec{ cdouble(1.0, 0.0), cdouble(2.0, 0.0) } };
	const itsc_block a(one, svec{ "C1" }, hkl_list{ i3{ 1, 0, 0 } });
	const itsc_block b(two, svec{ "O1" }, hkl_list{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 } });
	const itsc_block c(one, svec{ "N1" }, hkl_list{ i3{ 0, 0, 1 } });
	std::ostringstream log;
	EXPECT_EXIT({ itsc_block x = a; x.append(b, log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT({ itsc_block x = a; x.append(c, log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//validate_dimensions: the index constructor takes the arrays as given, the writers reject the broken ones
TEST(ScatteringCoverageBlockTests, ValidateDimensionsThrowsOnBrokenTables)
{
	scoped_cwd cwd("dims");
	const cvec2 one{ cvec{ cdouble(1.0, 0.0) } };
	{
		const itsc_block two_dims(one, svec{ "C1" }, ivec2{ { 1 }, { 0 } });
		EXPECT_THROW(two_dims.write_tsc_file("a.cif", "a.tsc"), std::runtime_error);
	}
	{
		const itsc_block ragged(one, svec{ "C1" }, ivec2{ { 1 }, { 0 }, { 0, 1 } });
		EXPECT_THROW(ragged.write_tscb_file("a.cif", "a.tscb"), std::runtime_error);
	}
	{
		const itsc_block short_row(one, svec{ "C1" }, ivec2{ { 1, 2 }, { 0, 0 }, { 0, 0 } });
		EXPECT_THROW(short_row.write_tsc_file("a.cif", "a.tsc"), std::runtime_error);
	}
	{
		//a duplicate label is dropped at construction, so the row count no longer matches the scatterers
		//and scatterer_size() reports 0 for the mismatch
		const cvec2 twice{ cvec{ cdouble(1.0, 0.0) }, cvec{ cdouble(2.0, 0.0) } };
		const itsc_block duplicate(twice, svec{ "C1", "C1" }, hkl_list{ i3{ 1, 0, 0 } });
		EXPECT_EQ(duplicate.get_scatterers().size(), 1u);
		EXPECT_EQ(duplicate.scatterer_size(), 0u);
		EXPECT_THROW(duplicate.write_tsc_file("a.cif", "a.tsc"), std::runtime_error);
	}
	EXPECT_FALSE(fs::exists(cwd.dir / "a.tsc") && fs::file_size(cwd.dir / "a.tsc") > 0);
}
