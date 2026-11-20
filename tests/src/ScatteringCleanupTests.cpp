// Tests written for the split of scattering_factors.cpp: functions that had no gtest of their own.
#include "pch.h"
#include <gtest/gtest.h>

#include "core/convenience.h"
#include "core/scattering_factors.h"
#include "core/atoms.h"

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
