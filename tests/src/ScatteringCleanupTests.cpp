// Tests written for the split of scattering_factors.cpp: functions that had no gtest of their own.
#include "pch.h"
#include <gtest/gtest.h>

#include "core/base/convenience.h"
#include "core/tsc/scattering_factors.h"
#include "core/wfn/atoms.h"
#include "core/xcw/cell.h"

#include <filesystem>
#include <fstream>

namespace
{
	namespace fs = std::filesystem;

	fs::path tmp(const std::string& name)
	{
		return fs::temp_directory_path() / ("nosphera2_sfclean_" + name);
	}

	void write_cif(const fs::path& cif, const std::string& atom_rows)
	{
		std::ofstream out(cif);
		out << "data_test\n_cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n"
			<< "loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n"
			<< "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_disorder_group\n"
			<< atom_rows;
	}
}

//the keys -merge holds a merged .tshb against: quoted and wrapped rows are one row, an unreadable disorder group is part 0,
//a type that is no element is left out and an atomID seen in an earlier CIF is not repeated
TEST(ScatteringCleanupTests, CifScattererKeysReadsEveryRowOnce)
{
	const fs::path a = tmp("keys_a.cif"), b = tmp("keys_b.cif");
	write_cif(a, "C1 C 0.25 0.35 0.45 .\n'O1' O 0.5 0.6\n 0.7 1\nQ1 Q 0.1 0.1 0.1 .\nN1 N 0.2 0.2 0.2 x\n\nloop_\n_atom_site_aniso_label\nC1\n");
	write_cif(b, "C1b C 0.25 0.35 0.45 ?\nH1 H 0.9 0.8 0.7 2\n");
	const auto keys = cif_scatterer_keys({ a, b });
	fs::remove(a);
	fs::remove(b);
	ASSERT_EQ(keys.size(), 4u);
	const svec labels = { "C1", "O1", "N1", "H1" };
	const atomID ids[4] = { atomID(0.25, 0.35, 0.45, 0, 6), atomID(0.5, 0.6, 0.7, 1, 8), atomID(0.2, 0.2, 0.2, 0, 7), atomID(0.9, 0.8, 0.7, 2, 1) };
	for (int i = 0; i < 4; i++)
	{
		EXPECT_EQ(keys[i].first, labels[i]);
		EXPECT_EQ(keys[i].second.to_hex_string(), ids[i].to_hex_string()) << labels[i];
	}
}

//read_mask: the fab rows follow the sorted hkl list whatever the file order, a header line is skipped and a reflection the
//file does not list stays 0
TEST(ScatteringCleanupTests, ReadMaskFollowsTheHklOrder)
{
	const fs::path cif = tmp("mask_p1.cif"), fab = tmp("mask_p1.fab");
	{
		std::ofstream c(cif);
		c << "data_test\n_cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n"
			<< "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
			<< "loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n";
		std::ofstream out(fab);
		out << "   h   k   l  real  imag\n   1   0   0   1.5   -0.5\n   0   0   1   2.5   0.25\n";
	}
	std::ostringstream log;
	cell c(cif, log, false, true);
	const hkl_list hkl = { i3{ 1, 0, 0 }, i3{ 0, 0, 1 }, i3{ 0, 1, 0 } };
	vec2 mask;
	read_mask(fab, hkl, hkl, c, log, mask, 3);
	fs::remove(cif);
	fs::remove(fab);
	ASSERT_EQ(mask.size(), 3u);
	EXPECT_EQ(mask[0], (vec{ 2.5, 0.25 }));
	EXPECT_EQ(mask[1], (vec{ 0.0, 0.0 }));
	EXPECT_EQ(mask[2], (vec{ 1.5, -0.5 }));
}

//read_mask: a fab file lists one reflection of each set of equivalents, F(h) = exp(2 pi i h.t) F(h R); under the 2_1 screw
//-x, y+1/2, -z the image (-1 1 0) of (1 1 0) carries the phase exp(i pi) = -1
TEST(ScatteringCleanupTests, ReadMaskTakesTheScrewPhase)
{
	const fs::path cif = tmp("mask_p21.cif"), fab = tmp("mask_p21.fab");
	{
		std::ofstream out(cif);
		out << "data_test\n_cell_length_a 4.0\n_cell_length_b 5.0\n_cell_length_c 6.0\n"
			<< "_cell_angle_alpha 90\n_cell_angle_beta 100\n_cell_angle_gamma 90\n"
			<< "loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n'-x, y+1/2, -z'\n";
		std::ofstream f(fab);
		f << "  -1   1   0   2.0   0.0\n";
	}
	std::ostringstream log;
	cell c(cif, log, false, true);
	const hkl_list hkl = { i3{ 1, 1, 0 } };
	const hkl_list enlarged = { i3{ 1, 1, 0 }, i3{ -1, 1, 0 }, i3{ 1, -1, 0 }, i3{ -1, -1, 0 } };
	vec2 mask;
	read_mask(fab, hkl, enlarged, c, log, mask, 1);
	fs::remove(cif);
	fs::remove(fab);
	ASSERT_EQ(mask.size(), 1u);
	EXPECT_NEAR(mask[0][0], -2.0, 1e-12);
	EXPECT_NEAR(mask[0][1], 0.0, 1e-12);
}
