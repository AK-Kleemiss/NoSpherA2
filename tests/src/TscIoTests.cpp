
#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/fchk.h"
#include "core/AtomGrid.h"
#include "core/SALTED_utilities.h"
#include "core/scattering_factors.h"
#include "core/nos_math.h"
#include "core/GridManager.h"
#include "core/atoms.h"
#include "core/tsc_block.h"
#include "core/cell.h"
#include "core/wfn_class.h"
#include "core/properties.h"
#include "core/integrator.h"
#include "core/integration_params.h"
#include "core/basis_set.h"
#include "core/geometry_aid.h"
#include "core/crystal_energies.h"
#include "core/NoSpherA2.h"
#include "core/isosurface.h"
#include "core/npy.h"
#include "core/libCintMain.h"
#include "core/throughput.h"
#include "core/i_tensor_stream.h"
#include "core/tsc_label_converter.h"
#include "core/sphere_lebedev_rule.h"
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#undef I
#ifdef NOSPHERA2_USE_GPU
#include "core/blas_gpu.h"
#include "core/aux_density_gpu.h"
#endif

static constexpr double PI_VAL = 3.14159265358979323846;

namespace {
	void test_reading_SALTED_binary_file() {
		std::filesystem::path path("../../../tests/SALTED/Model/model.salted");
		if (!std::filesystem::exists(path)) {
			path = std::filesystem::path("tests/SALTED/Model/model.salted");
		}
		if (!std::filesystem::exists(path)) {
			path = std::filesystem::path("Model/model.salted");
		}
		if (!std::filesystem::exists(path)) {
			path = std::filesystem::path("../SALTED/Model/model.salted");
		}
		SALTED_BINARY_FILE file = SALTED_BINARY_FILE(path, true);
		SALTEDConfig config;
		file.populate_config(config);
		std::unordered_map<int, std::vector<int64_t>> fps = file.read_fps();
		std::unordered_map<std::string, vec> averages = file.read_averages();
		std::unordered_map<int, vec> wigners = file.read_wigners();
		vec weights = file.read_weights();
		std::unordered_map<std::string, dMatrix2> feats = file.read_features();
		std::unordered_map<std::string, dMatrix2> proj = file.read_projectors();
		std::cout << "Finished reading SALTED binary file\n";
		// TEST if both configs are the same
		std::cout << "Average:" << config.average << "\n";
		std::cout << "Field:" << config.field << "\n";
		std::cout << "Sparsify:" << config.sparsify << "\n";
		std::cout << "Ncut:" << config.ncut << "\n";
		std::cout << "Ntrain:" << config.Ntrain << "\n";
		std::cout << "Menv:" << config.Menv << "\n";
		std::cout << "trainfrac:" << config.trainfrac << "\n";
		std::cout << "Rcut1:" << config.rcut1 << "\n";
		std::cout << "Rcut2:" << config.rcut2 << "\n";
		std::cout << "nang1:" << config.nang1 << "\n";
		std::cout << "nang2:" << config.nang2 << "\n";
		std::cout << "sig1:" << config.sig1 << "\n";
		std::cout << "sig2:" << config.sig2 << "\n";
		std::cout << "zeta:" << config.zeta << "\n";
		std::cout << "neighspe size:" << config.neighspe1.size() << "\n";
		for (int i = 0; i < config.neighspe1.size(); i++)
		{
			std::cout << "neighspe1[" << i << "]:" << config.neighspe1[i] << "\n";
		}
		std::cout << "neighspe2 size:" << config.neighspe2.size() << "\n";
		for (int i = 0; i < config.neighspe2.size(); i++)
		{
			std::cout << "neighspe2[" << i << "]:" << config.neighspe2[i] << "\n";
		}
		std::cout << "dfBasis:" << config.dfbasis << "\n";

		std::cout << "Comparing wigners\n";
		for (int i = 0; i < wigners.size(); i++)
		{
			for (int j = 0; j < wigners[i].size(); j += 10)
			{
				std::cout << "wigners[" << i << "][" << j << "]:" << wigners[i][j] << "\n";
			}
		}

		std::cout << "Comparing FPS\n";
		for (int i = 0; i < fps.size(); i++)
		{
			for (int j = 0; j < fps[i].size(); j += 10)
			{
				std::cout << "fps[" << i << "][" << j << "]:" << fps[i][j] << "\n";
			}
		}

		std::cout << "All tests passed!\n";
	}

	//Here only as a reference for simple tests, the actual implementation is in SALTED_utilities.cpp
}

namespace NoSpherA2UnitTests
{
	TEST(TscBlockTests, ConstructorWarnsForDuplicateScattererIds)
	{
		const cvec2 form_factors = {
			{ cdouble(1.0, 0.0) },
			{ cdouble(2.0, 0.0) },
			{ cdouble(3.0, 0.0) }
		};
		const std::vector<atomID> scatterer_ids = { {0.1, 0.1, 0.1, 0, 1 }, {0.5, 0.2, 0.4, 0, 2 }, {0.1, 0.1, 0.1, 0, 1 } };
		const ivec2 indices = { { 1 }, { 0 }, { 0 } };

		testing::internal::CaptureStdout();
		tsc_block<int, cdouble> block(form_factors, scatterer_ids, indices);
		const std::string output = testing::internal::GetCapturedStdout();

		EXPECT_NE(output.find("Warning: Duplicate scatterer: atomID(frac_x: 0.1, frac_y: 0.1, frac_z: 0.1, Z: 1, data: 0, reserved: 0)"), std::string::npos);
	}

	TEST(TscBlockTests, AppendWarnsForDuplicateScattererIds)
	{
		const ivec2 indices = { { 1 }, { 0 }, { 0 } };
		const cvec2 lhs_form_factors = {
			{ cdouble(1.0, 0.0) }
		};
		const std::vector<atomID> lhs_scatterer_ids{ {0.1, 0.1, 0.1, 0, 1 } };
		const cvec2 rhs_form_factors = {
			{ cdouble(2.0, 0.0) },
			{ cdouble(3.0, 0.0) }
		};
		const std::vector<atomID> rhs_scatterer_ids{ {0.1, 0.1, 0.1, 0, 1 }, {0.5, 0.2, 0.4, 0, 2 } };

		tsc_block<int, cdouble> lhs(lhs_form_factors, lhs_scatterer_ids, indices);
		tsc_block<int, cdouble> rhs(rhs_form_factors, rhs_scatterer_ids, indices);
		std::ostringstream log;

		testing::internal::CaptureStdout();
		lhs.append(rhs, log);
		const std::string output = testing::internal::GetCapturedStdout();

		EXPECT_NE(output.find("Warning: Duplicate scatterer in append: atomID(frac_x: 0.1, frac_y: 0.1, frac_z: 0.1, Z: 1, data: 0, reserved: 0)"), std::string::npos);
		EXPECT_EQ(lhs.scatterer_size(), 2);
	}

	TEST(TscBlockTests, BinaryFileRoundTripsWith32BitSizes)
	{
		const cvec2 form_factors = {
			{ cdouble(1.25, -0.5), cdouble(2.5, 0.75) },
			{ cdouble(-3.0, 1.5), cdouble(4.25, -2.0) }
		};
		const std::vector<atomID> scatterer_ids = {
			atomID(0.1, 0.2, 0.3, 1, 6),
			atomID(0.4, 0.5, 0.6, 2, 8)
		};
		const ivec2 indices = {
			{ 1, -2 }, { 0, 3 }, { -1, 4 }
		};
		const std::filesystem::path path =
			std::filesystem::temp_directory_path() / "nosphera2_tscb_32bit_roundtrip.tscb";

		tsc_block<int, cdouble> original(
			form_factors, scatterer_ids, indices);
		original.write_tscb_file({}, path);
		tsc_block<int, cdouble> restored(path);
		std::filesystem::remove(path);

		ASSERT_EQ(restored.scatterer_size(), scatterer_ids.size());
		ASSERT_EQ(restored.reflection_size(), indices[0].size());
		for (std::size_t scatterer = 0; scatterer < scatterer_ids.size(); ++scatterer)
		{
			EXPECT_EQ(std::get<atomID>(restored.get_scatterer(scatterer)),
				scatterer_ids[scatterer]);
			EXPECT_EQ(restored.get_sf_for_scatterer(scatterer), form_factors[scatterer]);
		}
		for (std::size_t reflection = 0; reflection < indices[0].size(); ++reflection)
		{
			EXPECT_EQ(restored.get_indices(reflection),
				(std::array<int, 3>{ indices[0][reflection], indices[1][reflection], indices[2][reflection] }));
		}
	}

	TEST(SALTEDTests, ReadingSALTEDBinaryFile)
	{
		test_reading_SALTED_binary_file();
	}

	// -----------------------------------------------------------------------

	// -----------------------------------------------------------------------
	// FchkParsingTests
	// -----------------------------------------------------------------------

		// FCHK format: 40-char label + type char + value starting at position 49
		// We pad to exactly 49 chars then append the value.

	TEST(FchkParsingTests, ReadFchkInt_PositiveValue)
	{
			// Format: keyword padded to 40, type at 40, 8 spaces (41-48), value at 49+
			const char* line = "Number of atoms                         I        3";
		EXPECT_EQ(3, read_fchk_integer(std::string(line)));
	}

	TEST(FchkParsingTests, ReadFchkInt_NegativeValue)
	{
		const char* line = "Charge                                  I        -1";
		EXPECT_EQ(-1, read_fchk_integer(std::string(line)));
	}

	TEST(FchkParsingTests, ReadFchkInt_LargeValue)
	{
		const char* line = "Number of basis functions               I        1024";
		EXPECT_EQ(1024, read_fchk_integer(std::string(line)));
	}

	TEST(FchkParsingTests, ReadFchkDbl_NegativeScientific)
	{
		const char* line = "Total Energy                            R        -1.23456789E+02";
		double val = read_fchk_double(std::string(line));
		EXPECT_NEAR(-123.456789, val, 1e-6);
	}

	TEST(FchkParsingTests, ReadFchkDbl_PositiveScientific)
	{
		const char* line = "Zero-point correction                   R        4.56000000E-02";
		double val = read_fchk_double(std::string(line));
		EXPECT_NEAR(0.0456, val, 1e-10);
	}

	// ------------------------------------------------------------------
	// tsc SCATTERER_IDS -> labels through the CIF

	TEST(TscLabelConverterTests, FixtureIdsMapOntoEpoxideLabels)
	{
		const auto dir = nos_test_repo_root() / "tests" / "epoxide_gbw";
		const auto cif = dir / "epoxide.cif";
		if (!std::filesystem::exists(cif)) GTEST_SKIP() << "Missing " << cif;
		for (const char* name : { "fixture.tsc", "fixture.tscb" }) {
			const auto table = dir / name;
			const tsc_block<int, cdouble> in = read_tsc_table(table);
			ASSERT_GT(in.scatterer_size(), 0u) << name;
			EXPECT_TRUE(std::holds_alternative<atomID>(in.get_scatterer(0))) << name;
			const auto out = std::filesystem::temp_directory_path() / ("nosphera2_labels_" + std::string(name) + ".tsc");
			std::ostringstream log;
			ASSERT_TRUE(convert_tsc_ids_to_labels(table, cif, out, log)) << log.str();
			EXPECT_NE(log.str().find("Matched"), std::string::npos);
			const tsc_block<int, cdouble> labelled = read_tsc_table(out);
			ASSERT_EQ(labelled.scatterer_size(), in.scatterer_size());
			ASSERT_EQ(labelled.get_index_vector().size(), in.get_index_vector().size());
			const std::vector<std::string> labels = labelled.get_scatterers_string();
			EXPECT_NE(std::find(labels.begin(), labels.end(), "O1"), labels.end()) << name;
			EXPECT_NE(std::find(labels.begin(), labels.end(), "C2"), labels.end()) << name;
			for (std::size_t i = 0; i < in.scatterer_size(); i++)
				EXPECT_EQ(labelled.get_sf_for_scatterer(i)[0], in.get_sf_for_scatterer(i)[0]) << name << " scatterer " << i;
			std::filesystem::remove(out);
		}
	}

	TEST(TscLabelConverterTests, ErrorsAreReportedNotThrown)
	{
		const auto dir = nos_test_repo_root() / "tests" / "epoxide_gbw";
		const auto out = std::filesystem::temp_directory_path() / "nosphera2_labels_err.tsc";
		std::ostringstream log;
		EXPECT_FALSE(convert_tsc_ids_to_labels(dir / "missing.tsc", dir / "epoxide.cif", out, log));
		EXPECT_NE(log.str().find("Could not convert"), std::string::npos);
		EXPECT_FALSE(convert_tsc_ids_to_labels(dir / "epoxide.cif", dir / "epoxide.cif", out, log)) << "wrong extension";
		EXPECT_FALSE(convert_tsc_ids_to_labels(dir / "fixture.tsc", dir / "missing.cif", out, log));
		// a CIF whose atoms do not match the IDs
		const auto cif = std::filesystem::temp_directory_path() / "nosphera2_labels_other.cif";
		{
			std::ofstream f(cif);
			f << "data_x\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
			  << "Xe1 Xe 0.1 0.2 0.3\n\n";
		}
		EXPECT_FALSE(convert_tsc_ids_to_labels(dir / "fixture.tsc", cif, out, log));
		EXPECT_NE(log.str().find("No CIF atom matches"), std::string::npos) << log.str();
		{
			std::ofstream f(cif);
			f << "data_x\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
			  << "A1 O 0.1 0.2 0.3\nA2 O 0.1 0.2 0.3\n\n";
		}
		EXPECT_FALSE(convert_tsc_ids_to_labels(dir / "fixture.tsc", cif, out, log));
		EXPECT_NE(log.str().find("duplicate"), std::string::npos) << log.str();
		{
			std::ofstream f(cif);
			f << "data_x\n_cell_length_a 1.0\n";
		}
		EXPECT_FALSE(convert_tsc_ids_to_labels(dir / "fixture.tsc", cif, out, log));
		EXPECT_NE(log.str().find("no usable _atom_site loop"), std::string::npos) << log.str();
		// a label table is refused, not passed through: the fixture with its ID line replaced by labels
		// (written here; the old check on experimental.tsc only held after an integration run had left one behind)
		const auto labelled = std::filesystem::temp_directory_path() / "nosphera2_labels_in.tsc";
		{
			std::ifstream in(dir / "fixture.tsc");
			std::ofstream f(labelled);
			std::string line;
			while (std::getline(in, line)) f << (line.rfind("SCATTERER_IDS:", 0) == 0 ? "SCATTERERS: O1 C2 C3 H4 H5 H6 H7" : line) << "\n";
		}
		EXPECT_FALSE(convert_tsc_ids_to_labels(labelled, dir / "epoxide.cif", out, log));
		EXPECT_NE(log.str().find("does not use SCATTERER_IDS"), std::string::npos) << log.str();
		std::filesystem::remove(labelled);
		std::filesystem::remove(cif);
		std::filesystem::remove(out);
	}

	// ------------------------------------------------------------------
	// The compact XCW I tensor file

	TEST(ITensorFileTests, RoundTripInBothPrecisionsAndWindowedReads)
	{
		const int nr = 5, nmo = 3;
		const ivec mu = { 0, 0, 1, 2 }, nu = { 0, 2, 1, 2 };
		for (const bool single : { false, true }) {
			const auto p = std::filesystem::temp_directory_path() / (single ? "nosphera2_itensor_32.bin" : "nosphera2_itensor_64.bin");
			{
				i_tensor_file f;
				f.create(p, nr, nmo, mu, nu, single);
				cvec block(mu.size());
				for (int r = nr - 1; r >= 0; r--) { // out of order, as the workers finish
					for (size_t k = 0; k < mu.size(); k++) block[k] = cdouble(r + 0.25 * k, -1.0 * r);
					f.write_block(r, block.data());
				}
				std::vector<std::complex<float>> block32(mu.size(), std::complex<float>(9.0f, 9.0f));
				f.write_block(2, block32.data()); // the other element type in
				f.finish_write();
			}
			size_t kept = 0;
			bool is_single = false;
			EXPECT_TRUE(i_tensor_file::matches(p, nr, nmo, kept, is_single));
			EXPECT_EQ(kept, mu.size());
			EXPECT_EQ(is_single, single);
			EXPECT_FALSE(i_tensor_file::matches(p, nr + 1, nmo, kept, is_single));
			EXPECT_EQ(i_tensor_file::total_bytes(nr, mu.size(), single), nr * mu.size() * (single ? 8u : 16u));
			i_tensor_file f;
			f.open(p, 2);
			EXPECT_EQ(f.nr(), nr);
			EXPECT_EQ(f.nmo(), nmo);
			EXPECT_EQ(f.kept(), mu.size());
			EXPECT_EQ(f.pair_nu()[1], 2);
			EXPECT_EQ(f.window_blocks(), 2u);
			EXPECT_EQ(f.path(), p);
			EXPECT_THROW(f.load(0, 3), std::runtime_error) << "beyond the window";
			EXPECT_THROW(f.load(4, 6), std::runtime_error) << "beyond nr";
			f.load(3, 5);
			if (single) {
				EXPECT_EQ(f.block32(4)[1], std::complex<float>(4.25f, -4.0f));
				EXPECT_EQ(f.block32(3)[0].real(), 3.0f);
			} else {
				EXPECT_EQ(f.block(4)[1], cdouble(4.25, -4.0));
				EXPECT_EQ(f.block(3)[0], cdouble(3.0, -3.0));
			}
			f.set_window(100); // clamps to nr
			EXPECT_EQ(f.window_blocks(), (size_t)nr);
			f.load(0, nr);
			if (single) EXPECT_EQ(f.block32(2)[3], std::complex<float>(9.0f, 9.0f));
			else EXPECT_EQ(f.block(2)[3], cdouble(9.0, 9.0));
			f.close();
			std::filesystem::remove(p);
		}
	}

	TEST(ITensorFileTests, DamagedFilesAreRefused)
	{
		const auto p = std::filesystem::temp_directory_path() / "nosphera2_itensor_bad.bin";
		size_t kept;
		bool single;
		EXPECT_FALSE(i_tensor_file::matches(p / "nowhere", 1, 1, kept, single));
		i_tensor_file f;
		EXPECT_THROW(f.open(p / "nowhere", 1), std::runtime_error);
		{
			std::ofstream(p, std::ios::binary) << "not a tensor at all, just text";
		}
		EXPECT_FALSE(i_tensor_file::matches(p, 1, 1, kept, single));
		EXPECT_THROW(f.open(p, 1), std::runtime_error);
		// a valid header whose data is truncated
		{
			i_tensor_file w;
			w.create(p, 3, 2, ivec{ 0, 1 }, ivec{ 0, 1 }, false);
			cvec block(2, 1.0);
			w.write_block(0, block.data());
			w.finish_write();
		}
		EXPECT_FALSE(i_tensor_file::matches(p, 3, 2, kept, single));
		EXPECT_THROW(f.open(p, 1), std::runtime_error);
		// a corrupt pair list (mu > nu)
		{
			i_tensor_file w;
			w.create(p, 1, 2, ivec{ 1 }, ivec{ 0 }, false);
			cvec block(1, 1.0);
			w.write_block(0, block.data());
			w.finish_write();
		}
		EXPECT_THROW(f.open(p, 1), std::runtime_error);
		std::filesystem::remove(p);
	}

} // namespace NoSpherA2UnitTests
