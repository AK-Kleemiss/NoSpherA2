#include "pch.h"
#include <gtest/gtest.h>

#include "core/convenience.h"
#include "core/NoSpherA2.h"
#include "core/cell.h"
#include "core/constants.h"
#include "core/gaussian_atom.h"
#include "core/scattering_factors.h"
#include "core/spherical_density.h"
#include "core/spherical_harmonic.h"
#include "core/tshb_block.h"
#include "core/GridManager.h"
#include "core/AtomGrid.h"
#include "core/wfn_class.h"
#include "core/basis_set.h"
#include "core/SALTED_utilities.h"
#include "core/npy.h"

#include <complex>
#include <fstream>
#include <limits>
#include <map>
#include <random>

//tshb_block (Src/core/tshb_block.h): the per-scatterer densities a .tshb carries, keyed and merged as tsc_block keys and
//merges rows, and the joint file of a -mtc run against the table the same run makes without -tshb
namespace
{
	namespace fs = std::filesystem;

	tshb_record record(const std::string& label, const int Z, const double c, const int part = 0)
	{
		tshb_record r;
		r.label = label;
		r.part = part;
		r.Z = Z;
		r.shells.push_back({ 0, 1.5, { c } });
		r.shells.push_back({ 1, 0.7, { 0.1 * c, -0.2 * c, 0.3 * c } });
		r.slaters.push_back({ 2.0, 1.0, 3.0 });
		r.gaussians.push_back({ 0.5, 4.0 });
		return r;
	}

	//the console a block prints its warnings on
	struct cout_capture
	{
		std::ostringstream text;
		std::streambuf* previous;
		cout_capture() : previous(std::cout.rdbuf(text.rdbuf())) {}
		~cout_capture() { std::cout.rdbuf(previous); }
	};

	fs::path scratch(const std::string& name)
	{
		const fs::path dir = fs::temp_directory_path() / ("nos_tshb_" + name + "_" +
			std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
		fs::create_directories(dir);
		return dir;
	}

	//run_app() in its own directory: the log and every file it writes land there
	struct app_dir
	{
		fs::path dir, previous;
		explicit app_dir(const std::string& name) : dir(scratch(name)), previous(fs::current_path()) { fs::current_path(dir); }
		~app_dir()
		{
			std::error_code ec;
			fs::current_path(previous, ec);
			fs::remove_all(dir, ec);
		}
		int run(std::vector<std::string> args) const
		{
			args.insert(args.begin(), "NoSpherA2");
			std::vector<char*> argv;
			for (std::string& a : args) argv.push_back(a.data());
			return run_app(static_cast<int>(argv.size()), argv.data());
		}
		std::string log() const
		{
			std::stringstream s;
			s << std::ifstream(dir / "NoSpherA2.log").rdbuf();
			return s.str();
		}
	};

	std::string bytes(const fs::path& p)
	{
		std::stringstream s;
		s << std::ifstream(p, std::ios::binary).rdbuf();
		return s.str();
	}

	std::string fixture(const std::string& rel)
	{
		return (nos_test_repo_root() / "tests" / rel).string();
	}

	//f of one record at k (A^-1, cartesian), independent of aux_density_table: the shells as the closed form
	//i^l pi^(3/2) (k/2)^l e^(-k^2/4a) a^-(l+3/2) Y_l(k^; c), a Slater term c R^n e^(-zR) as 4 pi c (n+1)! Im[(z - ik)^-(n+2)] / k,
	//an l = 0 Gaussian as c (pi/a)^(3/2) e^(-k^2/4a)
	cdouble record_form_factor(const tshb_record& r, const double* k)
	{
		const double H = std::sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
		const double d[3] = { H > 0 ? k[0] / H : 0.0, H > 0 ? k[1] / H : 0.0, H > 0 ? k[2] / H : 1.0 };
		cdouble f = 0.0;
		for (const tshb_record::shell& s : r.shells)
		{
			if (H == 0.0 && s.l > 0)
				continue;
			const double radial = std::pow(0.5 * H, s.l) * std::exp(-H * H / (4.0 * s.alpha)) / std::pow(s.alpha, s.l + 1.5);
			const double v = constants::PI3_2 * radial * constants::spherical_harmonic(s.l, d[0], d[1], d[2], s.c.data());
			const cdouble il[4] = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };
			f += il[s.l % 4] * v;
		}
		for (const auto& t : r.slaters)
		{
			const int n = static_cast<int>(std::lround(t[1]));
			double fact = 1.0;
			for (int i = 2; i <= n + 1; i++) fact *= i;
			if (H == 0.0)
				f += 4.0 * constants::PI * t[0] * fact * (n + 2) / std::pow(t[2], n + 3);
			else
				f += 4.0 * constants::PI * t[0] * fact * std::imag(std::pow(cdouble(t[2], -H), -(n + 2))) / H;
		}
		for (const auto& g : r.gaussians)
			f += g[0] * std::pow(constants::PI / g[1], 1.5) * std::exp(-H * H / (4.0 * g[1]));
		return f;
	}

	//every scatterer of the table once in the .tshb and the reverse; f(h) of each record against its row, max |dF|; the
	//largest |F| of the table comes back in f_max so a comparison of two empty tables cannot pass
	double compare_with_table(const tshb_block& tshb, const itsc_block& table, const cell& unit_cell, double& f_max, double* rel_max = nullptr)
	{
		std::map<std::string, size_t> row;
		const svec keys = table.get_scatterers_string();
		for (size_t i = 0; i < keys.size(); i++)
			EXPECT_TRUE(row.emplace(keys[i], i).second) << "table key twice: " << keys[i];
		const svec tkeys = tshb.get_scatterers_string();
		EXPECT_EQ(tkeys.size(), keys.size());
		const double scale = unit_cell.get_hkl_scale();
		double worst = 0.0;
		f_max = 0.0;
		for (size_t a = 0; a < tkeys.size(); a++)
		{
			const auto it = row.find(tkeys[a]);
			EXPECT_NE(it, row.end()) << tshb.get_record(a).label << " is in the .tshb, not in the table";
			if (it == row.end())
				continue;
			EXPECT_TRUE(tshb.get_record(a).f_terms.empty());
			const cvec& sf = table.get_sf_for_scatterer(it->second);
			for (size_t r = 0; r < table.reflection_size(); r++)
			{
				const auto h = table.get_indices(r);
				double k[3] = { 0, 0, 0 };
				for (int x = 0; x < 3; x++)
					for (int j = 0; j < 3; j++)
						k[x] += unit_cell.get_rcm(x, j) * h[j] / scale;
				for (double& v : k) v = constants::ang2bohr(v);
				const cdouble ft = record_form_factor(tshb.get_record(a), k);
				const double d = std::abs(ft - sf[r]);
				//a NaN would vanish in std::max: the first non-finite dF names its side and record, the rest stay inf
				if (!std::isfinite(d) && std::isfinite(worst)) {
					const tshb_record& rec = tshb.get_record(a);
					std::ostringstream s;
					for (const auto& sh : rec.shells) s << " (l " << sh.l << " a " << sh.alpha << ")";
					ADD_FAILURE() << rec.label << " h " << h[0] << " " << h[1] << " " << h[2] << ": .tshb f " << ft << ", table f " << sf[r]
						<< "; shells" << s.str() << ", " << rec.slaters.size() << " Slater, " << rec.gaussians.size() << " Gaussian terms";
				}
				worst = std::isfinite(d) ? std::max(worst, d) : std::numeric_limits<double>::infinity();
				f_max = std::max(f_max, std::abs(sf[r]));
				if (rel_max && std::abs(sf[r]) > 0.0)
					*rel_max = std::max(*rel_max, d / std::abs(sf[r]));
			}
		}
		return worst;
	}
}

TEST(TshbBlockTests, AppendKeepsTheFirstRecordOfAKeyAndWarns)
{
	tshb_block a({ record("C1", 6, 1.0), record("O1", 8, 2.0) }, svec{ "C1", "O1" });
	const tshb_block b({ record("O1", 8, 9.0), record("N1", 7, 3.0) }, svec{ "O1", "N1" });
	cout_capture out;
	a.append(b, std::cout);
	ASSERT_EQ(a.scatterer_size(), 3u);
	EXPECT_EQ(a.get_scatterers_string(), (svec{ "C1", "O1", "N1" }));
	EXPECT_EQ(a.get_record(1).shells[0].c[0], 2.0) << "a later O1 overwrote the first";
	EXPECT_NE(out.text.str().find("Warning: Duplicate scatterer in append: O1"), std::string::npos);
	EXPECT_EQ(b.scatterer_size(), 2u);
}

TEST(TshbBlockTests, PlusEqualsAndTheMoveOverloadMergeAlike)
{
	const tshb_block first({ record("C1", 6, 1.0) }, svec{ "C1" });
	tshb_block second({ record("C1", 6, 5.0), record("H1", 1, 0.5) }, svec{ "C1", "H1" });
	tshb_block by_copy;
	by_copy += first;
	by_copy += second;
	tshb_block by_move;
	by_move += tshb_block(first);
	cout_capture out;
	by_move += std::move(second);
	EXPECT_NE(out.text.str().find("Warning: Duplicate scatterer in append: C1"), std::string::npos);
	ASSERT_EQ(by_move.scatterer_size(), 2u);
	for (size_t i = 0; i < 2; i++)
		EXPECT_EQ(by_move.get_record(i), by_copy.get_record(i));
	EXPECT_EQ(by_move.get_record(0).shells[0].c[0], 1.0);
}

TEST(TshbBlockTests, ConstructorDropsADuplicateKeyLikeTscBlock)
{
	cout_capture out;
	const tshb_block a({ record("C1", 6, 1.0), record("C1", 6, 2.0) }, svec{ "C1", "C1" });
	EXPECT_EQ(a.scatterer_size(), 1u);
	EXPECT_EQ(a.get_record(0).shells[0].c[0], 1.0);
	EXPECT_NE(out.text.str().find("Warning: Duplicate scatterer: C1"), std::string::npos);
}

TEST(TshbBlockTests, MixedLabelAndIdKeysAreRefused)
{
	const ScattererLabels keys{ std::string("C1"), atomID(0.1, 0.2, 0.3, 0, 6) };
	EXPECT_THROW(tshb_block({ record("C1", 6, 1.0), record("C2", 6, 1.0) }, keys), std::runtime_error);
}

//-merge marks the result complete only when every input says so and every scatterer of the run is in it
TEST(TshbBlockTests, MergeIsCompleteOnlyWhenEveryInputIsAndTheRunIsCovered)
{
	const fs::path dir = scratch("coverage");
	const atomID c1(0.1, 0.1, 0.1, 0, 6), o1(0.2, 0.2, 0.2, 1, 8), n1(0.3, 0.3, 0.3, 2, 7), h1(0.4, 0.4, 0.4, 0, 1);
	tshb_block({ record("C1", 6, 1.0), record("O1", 8, 1.0, 1) }, std::vector<atomID>{ c1, o1 }).write(dir / "a.tshb");
	tshb_block({ record("N1", 7, 1.0, 2) }, std::vector<atomID>{ n1 }).write(dir / "b.tshb");
	tshb_block partial({ record("H1", 1, 1.0) }, std::vector<atomID>{ h1 });
	partial.set_complete(false);
	partial.write(dir / "h.tshb");
	EXPECT_FALSE(tshb_block(dir / "h.tshb").is_complete());
	const std::vector<std::pair<std::string, atomID>> run3{ { "C1", c1 }, { "O1", o1 }, { "N1", n1 } };
	const std::vector<std::pair<std::string, atomID>> run4{ { "C1", c1 }, { "O1", o1 }, { "N1", n1 }, { "H1", h1 } };
	auto merged = [&](const pathvec& files, const std::vector<std::pair<std::string, atomID>>* run)
	{
		cout_capture quiet;
		EXPECT_TRUE(merge_tshbs(files, dir / "m.tshb", run));
		const std::string text = bytes(dir / "m.tshb");
		const bool complete = tshb_block(dir / "m.tshb").is_complete();
		EXPECT_NE(text.find(complete ? "COVERAGE: complete" : "COVERAGE: partial"), std::string::npos);
		return complete;
	};
	EXPECT_TRUE(merged({ dir / "a.tshb", dir / "b.tshb" }, &run3));
	EXPECT_FALSE(merged({ dir / "a.tshb", dir / "b.tshb" }, &run4)) << "H1 is in no input";
	EXPECT_FALSE(merged({ dir / "a.tshb", dir / "b.tshb", dir / "h.tshb" }, &run4)) << "h.tshb is partial";
	EXPECT_FALSE(merged({ dir / "a.tshb", dir / "b.tshb" }, nullptr)) << "without the run's atoms nothing says the inputs are all of it";
	//label keys are held against the labels
	tshb_block({ record("C1", 6, 1.0) }, svec{ "C1" }).write(dir / "l.tshb");
	const std::vector<std::pair<std::string, atomID>> run_c{ { "C1", c1 } };
	EXPECT_TRUE(merged({ dir / "l.tshb" }, &run_c));
	fs::remove_all(dir);
}

//a truncated file is an error at every length, and random byte flips give a block or std::runtime_error: never a crash,
//another exception or an allocation sized by a corrupt count
TEST(TshbBlockTests, CorruptFilesAreErrorsNotCrashes)
{
	const fs::path dir = scratch("fuzz");
	tshb_record big = record("S1", 16, 1.0, 2);
	big.shells.push_back({ 2, 0.3, { 1, 2, 3, 4, 5 } });
	const std::vector<atomID> ids{ atomID(0.1, 0.2, 0.3, 1, 6), atomID(0.4, 0.5, 0.6, 2, 16) };
	tshb_block({ record("C1", 6, 1.0, 1), big }, ids).write(dir / "good.tshb");
	const std::string good = bytes(dir / "good.tshb");
	const fs::path bad = dir / "bad.tshb";
	auto write = [&bad](const std::string& b) { std::ofstream(bad, std::ios::binary | std::ios::trunc).write(b.data(), b.size()); };
	ASSERT_EQ(tshb_block(dir / "good.tshb").scatterer_size(), 2u);
	for (size_t len = 0; len < good.size(); len++)
	{
		write(good.substr(0, len));
		EXPECT_THROW(tshb_block{ bad }, std::runtime_error) << "truncated at " << len << " of " << good.size();
	}
	std::mt19937 gen(20261026);
	std::uniform_int_distribution<size_t> at(0, good.size() - 1);
	std::uniform_int_distribution<int> byte(0, 255), flips(1, 4);
	cout_capture quiet;
	int parsed = 0, refused = 0;
	for (int trial = 0; trial < 3000; trial++)
	{
		std::string b = good;
		for (int f = flips(gen); f > 0; f--)
			b[at(gen)] = static_cast<char>(byte(gen));
		write(b);
		try
		{
			const tshb_block back(bad);
			EXPECT_LE(back.scatterer_size(), 2u);
			parsed++;
		}
		catch (const std::runtime_error&)
		{
			refused++;
		}
		catch (const std::exception& e)
		{
			ADD_FAILURE() << "trial " << trial << ": " << typeid(e).name() << " " << e.what();
		}
	}
	EXPECT_EQ(parsed + refused, 3000);
	EXPECT_GT(refused, 0);
	fs::remove_all(dir);
}

TEST(TshbBlockTests, Version4RoundTripKeepsIdsPartsAndEd)
{
	const fs::path dir = scratch("v4");
	const std::vector<atomID> ids{ atomID(0.1, 0.2, 0.3, 1, 6), atomID(0.4, 0.5, 0.6, 2, 8), atomID(0.7, 0.8, 0.9, 0, 1) };
	const tshb_block a({ record("C1", 6, 1.0, 1), record("O1", 8, 2.0, 2), record("H1", 1, 3.0, 0) }, ids, true);
	a.write(dir / "a.tshb");
	const tshb_block b = tshb_block::read(dir / "a.tshb");
	ASSERT_EQ(b.scatterer_size(), 3u);
	EXPECT_EQ(b.get_scatterers(), a.get_scatterers());
	EXPECT_TRUE(b.get_ED());
	for (size_t i = 0; i < 3; i++)
		EXPECT_EQ(b.get_record(i), a.get_record(i));
	b.write(dir / "b.tshb");
	EXPECT_EQ(bytes(dir / "a.tshb"), bytes(dir / "b.tshb"));
	//the header the data implies: ED, complete, the parts present
	const std::string text = bytes(dir / "a.tshb");
	EXPECT_NE(text.find("ED: true\nCOVERAGE: complete\nPARTS: 0 1 2"), std::string::npos);
	const tshb_block labels({ record("C1", 6, 1.0) }, svec{ "C1" });
	labels.write(dir / "l.tshb");
	EXPECT_EQ(tshb_block::read(dir / "l.tshb").get_scatterers_string(), (svec{ "C1" }));
	fs::remove_all(dir);
}

//what Gaussian_Molecule::write_tshb wrote before version 4, byte by byte: read and keyed by label
TEST(TshbBlockTests, Version3IsReadKeyedByLabel)
{
	const fs::path dir = scratch("v3");
	{
		std::ofstream f(dir / "v3.tshb", std::ios::binary);
		auto u32 = [&f](const uint32_t v) { f.write(reinterpret_cast<const char*>(&v), 4); };
		auto str = [&](const std::string& s) { u32(static_cast<uint32_t>(s.size())); f.write(s.data(), s.size()); };
		auto f64 = [&f](const double v) { f.write(reinterpret_cast<const char*>(&v), 8); };
		f.write("TSHB", 4);
		u32(3);
		str("FRAME: cartesian\nUNITS: angstrom\nSOURCE: NoSpherA2\nED: true");
		u32(1);
		str("S1");
		u32(1);
		const int32_t l = 1;
		f.write(reinterpret_cast<const char*>(&l), 4);
		f64(2.5);
		f64(0.1); f64(0.2); f64(0.3);
		u32(1); f64(4.0); f64(9.0);
		u32(1); f64(1.0); f64(2.0); f64(3.0);
		u32(1); f64(10.0); f64(0.125);
		const int32_t Z = 16;
		f.write(reinterpret_cast<const char*>(&Z), 4);
	}
	const tshb_block b(dir / "v3.tshb");
	ASSERT_EQ(b.scatterer_size(), 1u);
	EXPECT_EQ(b.get_scatterers_string(), (svec{ "S1" }));
	EXPECT_TRUE(b.get_ED());
	const tshb_record& r = b.get_record(0);
	EXPECT_EQ(r.Z, 16);
	EXPECT_EQ(r.part, 0);
	ASSERT_EQ(r.shells.size(), 1u);
	EXPECT_EQ(r.shells[0].l, 1);
	EXPECT_EQ(r.shells[0].c, (vec{ 0.1, 0.2, 0.3 }));
	EXPECT_EQ(r.gaussians, (std::vector<std::array<double, 2>>{ { 4.0, 9.0 } }));
	EXPECT_EQ(r.slaters, (std::vector<std::array<double, 3>>{ { 1.0, 2.0, 3.0 } }));
	EXPECT_EQ(r.f_terms, (std::vector<std::array<double, 2>>{ { 10.0, 0.125 } }));
	fs::remove_all(dir);
}

//a Thakkar-filled atom as Slater records, with the fitted delta_k (Mg 1.6, Fe 2.5) or GTO (Cs, Au, U) differences as
//Gaussians, gives what the fill writes into the table
TEST(TshbBlockTests, FillRecordsReproduceTheSphericalFormFactor)
{
	const double b = constants::bohr2ang(1.0);
	for (const auto& [Z, q] : std::vector<std::pair<int, double>>{ { 1, 0.0 }, { 6, 0.0 }, { 8, -0.45 }, { 11, 0.7 }, { 17, -0.6 }, { 26, 0.9 }, { 7, 0.3 },
		{ 12, 1.6 }, { 26, 2.5 }, { 55, 0.8 }, { 79, 0.0 }, { 92, 2.5 } })
	{
		const HE_Spherical_Atom atom(Z, q);
		std::vector<std::array<double, 3>> slaters;
		std::vector<std::array<double, 2>> gaussians;
		atom.terms(slaters, gaussians);
		EXPECT_EQ(gaussians.empty(), !atom.uses_gto_tables() && !atom.uses_delta_series()) << Z << " " << q;
		tshb_record r;
		for (const auto& t : slaters) r.slaters.push_back(tshb_slater(t));
		for (const auto& g : gaussians) r.gaussians.push_back({ g[0] / (b * b * b), g[1] / (b * b) });
		double worst = 0.0;
		for (const double stl : { 0.0, 0.05, 0.2, 0.5, 1.0, 1.5 })
		{
			const double k_bohr = constants::bohr2ang(constants::FOUR_PI * stl); //the fill's k
			const double k[3] = { 0.0, 0.0, constants::ang2bohr(k_bohr) };
			worst = std::max(worst, std::abs(record_form_factor(r, k) - atom.get_form_factor(k_bohr)));
		}
		EXPECT_LT(worst, 1e-10) << "Z " << Z << " q " << q;
	}
}

//the joint .tshb of a -mtc run against the table of the same run without -tshb: every atom once, f(h) to round-off
static void joint_against_table(const std::string& name, const std::vector<std::string>& args, const std::string& cif)
{
	std::vector<std::string> with = args;
	with.push_back("-tshb");
	app_dir run(name);
	ASSERT_EQ(run.run(with), 0);
	EXPECT_NE(run.log().find("Wrote experimental.tshb:"), std::string::npos) << run.log();
	EXPECT_NE(run.log().find("Writing tsc file...  ... done!"), std::string::npos);
	EXPECT_FALSE(fs::exists(run.dir / "experimental.tscb"));
	const tshb_block tshb(run.dir / "experimental.tshb");
	fs::rename(run.dir / "NoSpherA2.log", run.dir / "tshb.log");
	ASSERT_EQ(run.run(args), 0);
	EXPECT_FALSE(fs::exists(run.dir / "experimental.tshb.tmp"));
	const itsc_block table(run.dir / "experimental.tscb");
	const cell unit_cell(cif, std::cout, false);
	double f_max = 0.0;
	const double worst = compare_with_table(tshb, table, unit_cell, f_max);
	std::cout << name << ": " << tshb.scatterer_size() << " atoms, " << table.reflection_size() << " reflections, max |F| " << f_max << ", max |dF| " << worst << std::endl;
	EXPECT_GT(f_max, 1.0);
	EXPECT_LT(worst, 1e-8);
}

TEST(TshbJointTests, HybridModeMatchesItsTable)
{
	const std::string cif = fixture("Hybrid/ZP2.cif");
	joint_against_table("hybrid", { "-cif", cif, "-dmin", "0.9", "-acc", "1", "-ri_fit", "combo_basis_fit",
		"-mtc", fixture("Hybrid/ZP2_part1.gbw"), "0.1", fixture("Hybrid/ZP2_part2.gbw"), "0.2",
		"-mtc_mult", "1", "1", "-mtc_charge", "0", "0", "-mtc_ECP", "0", "0", "-no_date" }, cif);
}

//An ECP gbw run without -ECP: the reader sets the wavefunction's mode to 1, and the streamed RI blocks must add that
//core, not the -ECP option's default 0. The same run with -ECP 1 is the reference. The non-streamed table is the
//ECP_auto_core TOML test; this covers the streamed call, which only writes the binary .tscb
TEST(EcpCoreTests, StreamedRiTableTakesTheReadersEcpMode)
{
	const std::vector<std::string> args = { "-wfn", fixture("TFVC/Rh.gbw"), "-cif", fixture("TFVC/Rh.cif"),
		"-hkl", fixture("TFVC/Rh_ECP_subset.hkl"), "-acc", "1", "-ri_fit", "def2-universal-jkfit", "-tsc_block", "30",
		"-no_date" };
	app_dir run("ecp_stream");
	ASSERT_EQ(run.run(args), 0);
	EXPECT_NE(run.log().find("ECP block removes 28 core electrons"), std::string::npos) << run.log();
	EXPECT_NE(run.log().find("Streaming tsc in blocks of 30 reflections"), std::string::npos) << run.log();
	const itsc_block detected(run.dir / "experimental.tscb");
	ASSERT_TRUE(fs::remove(run.dir / "experimental.tscb"));
	std::vector<std::string> explicit_mode = args;
	explicit_mode.insert(explicit_mode.end(), { "-ECP", "1" });
	ASSERT_EQ(run.run(explicit_mode), 0);
	const itsc_block reference(run.dir / "experimental.tscb");
	ASSERT_EQ(detected.scatterer_size(), reference.scatterer_size());
	ASSERT_EQ(detected.reflection_size(), reference.reflection_size());
	double f_max = 0.0, worst = 0.0;
	for (size_t a = 0; a < reference.scatterer_size(); a++)
		for (size_t r = 0; r < reference.reflection_size(); r++)
		{
			f_max = std::max(f_max, std::abs(reference.get_sf_for_scatterer(a)[r]));
			worst = std::max(worst, std::abs(detected.get_sf_for_scatterer(a)[r] - reference.get_sf_for_scatterer(a)[r]));
		}
	EXPECT_GT(f_max, 10.0);
	EXPECT_LT(worst, 1e-8) << "max |df| between the detected and the explicit ECP mode";
}

//The def2 ECP valence density plus the def2 core against an all-electron X2C wavefunction of the same HgI2 (Hg 60, I 28
//core electrons), Hirshfeld, s = 0.18-0.71 A^-1: they agree to 0.062 e. The mode-0 core was 11 e off on Hg, 1.2 e on I
TEST(EcpCoreTests, Def2CoreMatchesAllElectronX2c)
{
	const std::vector<std::string> common = { "-cif", fixture("ECP_HgI2/HgI2.cif"), "-hkl", fixture("ECP_HgI2/HgI2_subset.hkl"),
		"-acc", "1", "-no_date" };
	const std::string x2c = fixture("ECP_HgI2/HgI2_x2c.gbw"), def2_gbw = fixture("ECP_HgI2/HgI2_def2.gbw");
	app_dir run("ecp_x2c");
	std::vector<std::string> args = common;
	args.insert(args.end(), { "-wfn", x2c });
	ASSERT_EQ(run.run(args), 0);
	EXPECT_EQ(run.log().find("ECP block removes"), std::string::npos) << run.log();
	const itsc_block all_electron(run.dir / "experimental.tscb");
	ASSERT_TRUE(fs::remove(run.dir / "experimental.tscb"));
	args = common;
	args.insert(args.end(), { "-wfn", def2_gbw });
	ASSERT_EQ(run.run(args), 0);
	EXPECT_NE(run.log().find("ECP block removes 116 core electrons"), std::string::npos) << run.log();
	const itsc_block def2(run.dir / "experimental.tscb");
	ASSERT_EQ(def2.get_scatterers_string(), all_electron.get_scatterers_string());
	ASSERT_EQ(def2.reflection_size(), all_electron.reflection_size());
	ASSERT_EQ(def2.scatterer_size(), 3u);
	for (size_t a = 0; a < def2.scatterer_size(); a++)
	{
		double f_max = 0.0, worst = 0.0;
		for (size_t r = 0; r < def2.reflection_size(); r++)
		{
			f_max = std::max(f_max, std::abs(all_electron.get_sf_for_scatterer(a)[r]));
			worst = std::max(worst, std::abs(def2.get_sf_for_scatterer(a)[r] - all_electron.get_sf_for_scatterer(a)[r]));
		}
		EXPECT_GT(f_max, 30.0) << def2.get_scatterers_string()[a];
		EXPECT_LT(worst, 0.15) << def2.get_scatterers_string()[a] << ": max |df| def2 + core against all-electron";
	}
}

//-coef with -ECP 2: the coefficients fit the valence density, so the table, the charge table and the .tshb all add
//cysteine's 22 xTB core electrons, in mode 2; the table without -ECP is the one without core. Below Rb modes 1 and 2
//add the same core, so an iodine atom with zero coefficients, whose table is its core alone, tells them apart
TEST(EcpCoreTests, CoefficientFileTakesTheEcpCoreEverywhere)
{
	//the core add_ECP_contribution adds in mode m for the xTB core counts -ECP 2 puts on the atoms
	auto core = [](const int Z, const int m, const double stl) {
		const int n = constants::ECP_core_electrons(constants::ECP_electrons_xTB, Z);
		const double k = constants::FOUR_PI * constants::bohr2ang(stl);
		if (n == 0) return 0.0;
		if (m == 0) return n * std::exp(-stl / 4.0);
		return Thakkar(Z, m).get_core_form_factor(k, n) + Spherical_Gaussian_Density(Z, m).get_form_factor(k);
	};
	const std::string cif = fixture("SALTED/test_cysteine.cif");
	const std::vector<std::string> base = { "-wfn", fixture("SALTED/test_cysteine.xyz"), "-cif", cif, "-dmin", "0.73",
		"-coef", fixture("SALTED/test_cysteine.npy"), "cc-pvqz-jkfit", "-no_gpu_salted", "-no_date" };
	auto electrons = [](const std::string& log) {
		const size_t p = log.rfind("Total number of analytical Electrons: ");
		EXPECT_NE(p, std::string::npos) << log;
		return p == std::string::npos ? 0.0 : std::stod(log.substr(p + 38));
	};
	app_dir run("coef_ecp");
	ASSERT_EQ(run.run(base), 0);
	const double n0 = electrons(run.log());
	const itsc_block none(run.dir / "experimental.tscb");
	ASSERT_TRUE(fs::remove(run.dir / "experimental.tscb"));
	std::vector<std::string> ecp = base;
	ecp.insert(ecp.end(), { "-ECP", "2" });
	ASSERT_EQ(run.run(ecp), 0);
	const double n2 = electrons(run.log());
	const itsc_block table(run.dir / "experimental.tscb");
	ecp.push_back("-tshb");
	ASSERT_EQ(run.run(ecp), 0);
	const tshb_block tshb(run.dir / "experimental.tshb");
	EXPECT_NEAR(n2 - n0, 22.0, 1e-3);
	const cell unit_cell(cif, std::cout, false);
	ASSERT_EQ(table.get_scatterers_string(), none.get_scatterers_string());
	ASSERT_EQ(table.reflection_size(), none.reflection_size());
	//every reflection, so also the highest s, where the mode-0 core lies 3.7 e above the mode-2 one
	size_t high = 0;
	double worst = 0.0, sum[4] = { 0, 0, 0, 0 };
	for (size_t r = 0; r < table.reflection_size(); r++)
		if (unit_cell.get_stl_of_hkl(table.get_indices(r)) > unit_cell.get_stl_of_hkl(table.get_indices(high))) high = r;
	for (size_t a = 0; a < table.scatterer_size(); a++)
		for (size_t r = 0; r < table.reflection_size(); r++)
		{
			const double s = unit_cell.get_stl_of_hkl(table.get_indices(r));
			const cdouble d = table.get_sf_for_scatterer(a)[r] - none.get_sf_for_scatterer(a)[r];
			worst = std::max(worst, std::abs(d - core(std::get<atomID>(table.get_scatterer(a)).Z(), 2, s)));
			if (r != high) continue;
			sum[3] += d.real();
			for (int m = 0; m < 3; m++) sum[m] += core(std::get<atomID>(table.get_scatterer(a)).Z(), m, s);
		}
	std::cout << "core at s = " << unit_cell.get_stl_of_hkl(table.get_indices(high)) << ": table " << sum[3] << ", mode 0 "
		<< sum[0] << ", mode 1 " << sum[1] << ", mode 2 " << sum[2] << ", max |df| against mode 2 " << worst << std::endl;
	EXPECT_LT(worst, 1e-6) << "max |df| of the table's core against the mode-2 core";
	double f_max = 0.0;
	EXPECT_LT(compare_with_table(tshb, table, unit_cell, f_max), 1e-8);
	EXPECT_GT(f_max, 1.0);
	npy::npy_data<double> zero;
	WFN iodine(e_origin::NOT_YET_DEFINED);
	iodine.push_back_atom("I", 0.0, 0.0, 0.0, 53);
	std::vector<std::shared_ptr<BasisSet>> aux{ BasisSetLibrary::get_basis_set("def2-universal-jkfit") };
	zero.data.assign(aux_density_table(generate_aux_wfn(iodine, aux).get_atoms()).n_coef, 0.0);
	zero.shape = { static_cast<unsigned long>(zero.data.size()) };
	npy::write_npy<double>("I.npy", zero);
	std::ofstream("I.xyz") << "1\n\nI 0.0 0.0 0.0\n";
	std::ofstream("I.cif") << "data_I\nloop_\n _space_group_symop_id\n _space_group_symop_operation_xyz\n 1 +X,+Y,+Z\n"
		"_cell_length_a 10.0\n_cell_length_b 10.0\n_cell_length_c 10.0\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
		"_cell_volume 1000.0\nloop_\n _atom_site_label\n _atom_site_type_symbol\n _atom_site_fract_x\n _atom_site_fract_y\n"
		" _atom_site_fract_z\n _atom_site_U_iso_or_equiv\nI1 I 0.0 0.0 0.0 0.02\n";
	ASSERT_EQ(run.run({ "-wfn", "I.xyz", "-cif", "I.cif", "-dmin", "0.73", "-coef", "I.npy", "def2-universal-jkfit", "-ECP", "2",
		"-no_gpu_salted", "-no_date" }), 0);
	const itsc_block i_table(run.dir / "experimental.tscb");
	const cell i_cell(run.dir / "I.cif", std::cout, false);
	ASSERT_EQ(i_table.scatterer_size(), 1u);
	high = 0;
	worst = 0.0;
	for (size_t r = 0; r < i_table.reflection_size(); r++)
	{
		const double s = i_cell.get_stl_of_hkl(i_table.get_indices(r));
		worst = std::max(worst, std::abs(i_table.get_sf_for_scatterer(0)[r] - core(53, 2, s)));
		if (s > i_cell.get_stl_of_hkl(i_table.get_indices(high))) high = r;
	}
	const double s = i_cell.get_stl_of_hkl(i_table.get_indices(high));
	std::cout << "iodine at s = " << s << ": table " << i_table.get_sf_for_scatterer(0)[high].real() << ", mode 0 " << core(53, 0, s)
		<< ", mode 1 " << core(53, 1, s) << ", mode 2 " << core(53, 2, s) << ", max |df| against mode 2 " << worst << std::endl;
	EXPECT_LT(worst, 1e-6) << "max |df| of iodine's table against its mode-2 core";
}

//-cmtc: four fragments with their own CIFs; C_3, N_3 and O_3 are the same PART 0 atoms in fragments 2 and 3
static std::vector<std::string> rubredoxin_fragment(const int f, const std::string& groups)
{
	const std::string r = fixture("rubredoxin_cmtc/residues/") + std::to_string(f);
	return { r + ".gbw", r + ".cif", groups };
}

static std::vector<std::string> rubredoxin(const std::vector<int>& fragments)
{
	const char* groups[5] = { "", "0", "0.1", "0.2", "0" };
	std::vector<std::string> a{ "-hkl", fixture("rubredoxin_cmtc/1yk4_h.hkl"), "-acc", "1", "-ri_fit", "combo_basis_fit", "-no_date", "-cmtc" };
	for (const int f : fragments)
		for (const std::string& x : rubredoxin_fragment(f, groups[f])) a.push_back(x);
	for (const char* key : { "-mtc_mult", "-mtc_charge", "-mtc_ECP" })
	{
		a.push_back(key);
		for (size_t i = 0; i < fragments.size(); i++) a.push_back(key == std::string("-mtc_mult") ? "1" : "0");
	}
	return a;
}

TEST(TshbJointTests, RubredoxinFragmentsMatchTheirTableWithSharedAtomsOnce)
{
	joint_against_table("rubredoxin", rubredoxin({ 1, 2, 3, 4 }), fixture("rubredoxin_cmtc/residues/1.cif"));
}

//per-fragment .tshb files merged with -merge are the joint file of the -cmtc run, byte for byte: first wins in both
TEST(TshbJointTests, MergeOfTheFragmentFilesIsTheJointFile)
{
	std::vector<std::vector<std::string>> runs;
	runs.push_back(rubredoxin({ 1, 2, 3, 4 }));
	for (const int f : { 1, 2, 3, 4 })
		runs.push_back(rubredoxin({ f }));
	std::vector<std::string> merge_cif{ "-merge", "f1.tshb", "f2.tshb", "f3.tshb", "f4.tshb", "-cif" };
	for (const int f : { 1, 2, 3, 4 })
		merge_cif.push_back(fixture("rubredoxin_cmtc/residues/") + std::to_string(f) + ".cif");
	app_dir run("merge");
	for (size_t i = 0; i < runs.size(); i++)
	{
		runs[i].push_back("-tshb");
		ASSERT_EQ(run.run(runs[i]), 0) << i;
		fs::rename(run.dir / "experimental.tshb", run.dir / (i == 0 ? std::string("joint.tshb") : "f" + std::to_string(i) + ".tshb"));
	}
	//without the CIFs nothing says the four fragments are the whole structure
	ASSERT_EQ(run.run({ "-merge", "f1.tshb", "f2.tshb", "f3.tshb", "f4.tshb" }), 0);
	EXPECT_FALSE(tshb_block(run.dir / "combined.tshb").is_complete());
	ASSERT_EQ(run.run(merge_cif), 0);
	ASSERT_TRUE(fs::exists(run.dir / "combined.tshb"));
	EXPECT_TRUE(tshb_block(run.dir / "combined.tshb").is_complete());
	EXPECT_EQ(bytes(run.dir / "combined.tshb"), bytes(run.dir / "joint.tshb"));
	size_t sum = 0;
	for (const int f : { 1, 2, 3, 4 })
		sum += tshb_block(run.dir / ("f" + std::to_string(f) + ".tshb")).scatterer_size();
	const tshb_block joint(run.dir / "joint.tshb");
	EXPECT_EQ(joint.scatterer_size() + 3, sum) << "C_3, N_3 and O_3 are in two fragments and once in the joint file";
}

//SALTED with atoms the model lacks (Co, F in Co110): those come in as Slater records in their EEQ charge, the others as
//predicted shells, and all of it matches the table the run streams without -tshb
TEST(TshbJointTests, SaltedThakkarFillMatchesTheStreamedTable)
{
	const std::string model = fixture("SALTED/Model"), cif = fixture("SALTED/Co110.cif"), xyz = fixture("SALTED/Co110.xyz");
	app_dir run("salted_fill");
	const std::vector<std::string> args{ "-SALTED", model, "-cif", cif, "-wfn", xyz, "-dmin", "0.8", "-no_date" };
	std::vector<std::string> with = args;
	with.push_back("-tshb");
	ASSERT_EQ(run.run(with), 0);
	const std::string log = run.log();
	EXPECT_NE(log.find("Spherical fill: 7 atom(s) go into the .tshb as Slater records"), std::string::npos) << log;
	EXPECT_NE(log.find("Wrote experimental.tshb: 49 atom(s), 7 of them spherical"), std::string::npos) << log;
	EXPECT_FALSE(fs::exists(run.dir / "experimental.tscb"));
	const tshb_block tshb(run.dir / "experimental.tshb");
	ASSERT_EQ(run.run(args), 0);
	const itsc_block table(run.dir / "experimental.tscb");
	double f_max = 0.0;
	const double worst = compare_with_table(tshb, table, cell(cif, std::cout, false), f_max);
	std::cout << "Co110: " << tshb.scatterer_size() << " atoms, " << table.reflection_size() << " reflections, max |F| " << f_max << ", max |dF| " << worst << std::endl;
	EXPECT_GT(f_max, 1.0);
	EXPECT_LT(worst, 1e-8);
}

//the same structure with Rh in place of Co: Rh's fitted GTO difference goes into its record as Gaussians, and the .tshb
//matches the table of the run without -tshb
TEST(TshbJointTests, AGtoFillAtomGoesIntoTheTshbAsGaussians)
{
	const std::string model = fixture("SALTED/Model"), cif_in = fixture("SALTED/Co110.cif"), xyz_in = fixture("SALTED/Co110.xyz");
	app_dir run("salted_gto");
	auto rewrite = [](const std::string& from, const fs::path& to, const std::string& a, const std::string& b)
	{
		std::stringstream in;
		in << std::ifstream(from).rdbuf();
		std::string text = in.str();
		const size_t at = text.find(a);
		ASSERT_NE(at, std::string::npos);
		text.replace(at, a.size(), b);
		std::ofstream(to) << text;
	};
	rewrite(cif_in, run.dir / "rh.cif", " Co Co ", " Rh Rh ");
	rewrite(xyz_in, run.dir / "rh.xyz", "\nCo ", "\nRh ");
	const std::vector<std::string> args{ "-SALTED", model, "-cif", "rh.cif", "-wfn", "rh.xyz", "-dmin", "0.8", "-no_date" };
	std::vector<std::string> with = args;
	with.push_back("-tshb");
	ASSERT_EQ(run.run(with), 0);
	EXPECT_NE(run.log().find("Wrote experimental.tshb: 49 atom(s), 7 of them spherical"), std::string::npos) << run.log();
	EXPECT_FALSE(fs::exists(run.dir / "experimental.tscb"));
	const tshb_block tshb(run.dir / "experimental.tshb");
	size_t n_gauss = 0;
	for (size_t i = 0; i < tshb.scatterer_size(); i++)
		n_gauss += tshb.get_record(i).Z == 45 ? tshb.get_record(i).gaussians.size() : 0;
	EXPECT_GT(n_gauss, 0u);
	ASSERT_EQ(run.run(args), 0);
	const itsc_block table(run.dir / "experimental.tscb");
	double f_max = 0.0;
	const double worst = compare_with_table(tshb, table, cell((run.dir / "rh.cif").string(), std::cout, false), f_max);
	std::cout << "Rh: max |F| " << f_max << ", max |dF| " << worst << std::endl;
	EXPECT_LT(worst, 1e-8);
}

//a run the .tshb cannot stand in for writes the table, no .tshb, and says why
TEST(TshbJointTests, AGridPartitionKeepsTheTableAndSaysWhy)
{
	const std::string wfn = fixture("epoxide_gbw/epoxide.gbw"), cif = fixture("epoxide_gbw/epoxide.cif");
	app_dir run("grid");
	ASSERT_EQ(run.run({ "-wfn", wfn, "-cif", cif, "-dmin", "1.5", "-acc", "1", "-tshb", "-no_date" }), 0);
	EXPECT_NE(run.log().find("-tshb: the table is written instead, the partition is no fitted density"), std::string::npos) << run.log();
	EXPECT_TRUE(fs::exists(run.dir / "experimental.tscb"));
	EXPECT_FALSE(fs::exists(run.dir / "experimental.tshb"));
}

//disorder under SALTED: PART 1 and PART 2 predicted from their own wavefunctions into one .tshb, each part's atoms once and
//against the table of the same run; the two single-part files merged with the CIF are the joint file byte for byte
TEST(TshbJointTests, SaltedDisorderedPartsMatchTheirTableAndMergeToTheJointFile)
{
	const std::string cif = fs::absolute(fixture("Hybrid/ZP2.cif")).string(), model = fixture("SALTED/Model");
	auto parts = [&](const std::vector<int>& p)
	{
		std::vector<std::string> a{ "-SALTED", model, "-cif", cif, "-dmin", "0.9", "-no_date", "-mtc" };
		for (const int i : p)
		{
			a.push_back(fixture("Hybrid/ZP2_part") + std::to_string(i) + ".gbw");
			a.push_back("0." + std::to_string(i));
		}
		for (const std::string f : { "-mtc_mult", "-mtc_charge", "-mtc_ECP" })
		{
			a.push_back(f);
			for (size_t i = 0; i < p.size(); i++) a.push_back(f == "-mtc_mult" ? "1" : "0");
		}
		return a;
	};
	joint_against_table("salted_disorder", parts({ 1, 2 }), cif);
	//the fixture paths are relative to the starting directory, so every argument list is made before app_dir moves away
	std::vector<std::vector<std::string>> runs{ parts({ 1, 2 }), parts({ 1 }), parts({ 2 }) };
	const std::string names[3] = { "joint.tshb", "p1.tshb", "p2.tshb" };
	app_dir run("salted_disorder_merge");
	for (int i = 0; i < 3; i++)
	{
		std::vector<std::string>& a = runs[i];
		a.push_back("-tshb");
		ASSERT_EQ(run.run(a), 0) << i;
		fs::rename(run.dir / "experimental.tshb", run.dir / names[i]);
	}
	const tshb_block joint(run.dir / "joint.tshb"), p1(run.dir / "p1.tshb"), p2(run.dir / "p2.tshb");
	EXPECT_EQ(joint.scatterer_size(), 50u);
	EXPECT_EQ(p1.scatterer_size() + p2.scatterer_size(), 50u) << "no PART 0 atoms: nothing shared";
	ASSERT_EQ(run.run({ "-merge", "p1.tshb", "p2.tshb", "-cif", cif }), 0);
	EXPECT_TRUE(tshb_block(run.dir / "combined.tshb").is_complete());
	EXPECT_EQ(bytes(run.dir / "combined.tshb"), bytes(run.dir / "joint.tshb"));
}

//a part without orbitals (an .xyz) in a -ri_fit -tshb run: its atoms go in as the fill's spherical records, neutral and at
//full occupancy, and match the table of the run without -tshb, which takes IAM for that part; the file is complete, and
//the per-part files merged with the CIF are the joint file byte for byte
TEST(TshbJointTests, AnOrbitalLessPartGoesInAsSphericalRecords)
{
	const std::string cif = fs::absolute(fixture("Hybrid/ZP2.cif")).string();
	auto parts = [&](const std::vector<std::string>& files)
	{
		std::vector<std::string> a{ "-cif", cif, "-dmin", "0.9", "-acc", "1", "-ri_fit", "combo_basis_fit", "-no_date", "-mtc" };
		for (const std::string& f : files)
		{
			a.push_back(fs::absolute(fixture("Hybrid/" + f)).string());
			a.push_back(f.find("part1") != std::string::npos ? "0.1" : "0.2");
		}
		for (const std::string f : { "-mtc_mult", "-mtc_charge", "-mtc_ECP" })
		{
			a.push_back(f);
			for (size_t i = 0; i < files.size(); i++) a.push_back(f == "-mtc_mult" ? "1" : "0");
		}
		return a;
	};
	std::vector<std::vector<std::string>> runs{ parts({ "ZP2_part1.gbw", "ZP2_part2.xyz" }), parts({ "ZP2_part1.gbw" }), parts({ "ZP2_part2.xyz" }) };
	const std::string names[3] = { "joint.tshb", "p1.tshb", "p2.tshb" };
	app_dir run("orbitalless");
	ASSERT_EQ(run.run(runs[0]), 0);
	const itsc_block table(run.dir / "experimental.tscb");
	for (int i = 0; i < 3; i++)
	{
		runs[i].push_back("-tshb");
		ASSERT_EQ(run.run(runs[i]), 0) << i;
		EXPECT_EQ(run.log().find("-tshb: the table is written instead"), std::string::npos) << run.log();
		fs::rename(run.dir / "experimental.tshb", run.dir / names[i]);
	}
	EXPECT_NE(run.log().find("Spherical fill: 25 atom(s) go into the .tshb as Slater records"), std::string::npos) << run.log();
	const tshb_block joint = tshb_block::read(run.dir / "joint.tshb");
	EXPECT_EQ(joint.scatterer_size(), 50u);
	EXPECT_TRUE(joint.is_complete());
	EXPECT_NE(bytes(run.dir / "joint.tshb").find("COVERAGE: complete\nPARTS: 1 2"), std::string::npos);
	size_t spherical = 0;
	for (size_t i = 0; i < joint.scatterer_size(); i++)
		if (joint.get_record(i).shells.empty())
		{
			spherical++;
			EXPECT_EQ(joint.get_record(i).part, 2);
			EXPECT_GT(joint.get_record(i).Z, 0);
		}
	EXPECT_EQ(spherical, 25u);
	double f_max = 0.0, rel = 0.0;
	const double worst = compare_with_table(joint, table, cell(cif, std::cout, false), f_max, &rel);
	std::cout << "orbital-less part: " << joint.scatterer_size() << " atoms, " << table.reflection_size() << " reflections, max |F| " << f_max
		<< ", max |dF| " << worst << ", max |dF|/|F| " << rel << std::endl;
	EXPECT_LT(worst, 1e-8);
	ASSERT_EQ(run.run({ "-merge", "p1.tshb", "p2.tshb", "-cif", cif }), 0);
	EXPECT_TRUE(tshb_block(run.dir / "combined.tshb").is_complete());
	EXPECT_EQ(bytes(run.dir / "combined.tshb"), bytes(run.dir / "joint.tshb"));
}

//XCW with -tshb: every lambda step writes the RI fit of its converged density (AutoAux without -ri_fit) as NA2_<lambda>.tshb
//instead of the table, complete and keyed as the table's scatterers
TEST(TshbJointTests, XcwStepsWriteTheirTshb)
{
	const std::string cif = fs::absolute(fixture("P1_test/P1_test.cif")).string(), hkl = fs::absolute(fixture("P1_test/P1_subset.hkl")).string();
	app_dir run("xcw");
	std::ofstream(run.dir / "settings.txt") << "normal normal_conv params 177 basis_set sto-3g max_iter 100 F charge 0 mult 1 rhf start 0 step_size 0.01 end 0.01";
	ASSERT_EQ(run.run({ "-do_XCW", "-no_gpu", "-cif", cif, "-hkl", hkl, "-XCW_settings", "settings.txt", "-tshb", "-no_date" }), 0);
	for (const std::string step : { "NA2_0000000", "NA2_0010000" })
	{
		ASSERT_TRUE(fs::exists(run.dir / (step + ".tshb"))) << step;
		EXPECT_FALSE(fs::exists(run.dir / (step + ".tscb"))) << step;
		const tshb_block b = tshb_block::read(run.dir / (step + ".tshb"));
		EXPECT_EQ(b.scatterer_size(), 23u);
		EXPECT_TRUE(b.is_complete());
		EXPECT_NE(bytes(run.dir / (step + ".tshb")).find("COVERAGE: complete"), std::string::npos);
		for (size_t i = 0; i < b.scatterer_size(); i++)
			EXPECT_FALSE(b.get_record(i).shells.empty()) << b.get_record(i).label;
	}
}

//XCW -xcw_partition ri -tshb: each step's .tshb is the RI partition's own fit (XR.c mapped onto NoSpherA2's aux table shell
//by shell), so its F is the RI route's table to round-off (the separate density fit, NOSPHERA2_XCW_TSHB_FIT, is 1e-7 off). The
//table's scatterers are matched to their atoms by id: by closest f (NOSPHERA2_XCW_TSCB_BY_F) the RI table gave two atoms of
//one element each other's f, 5e-2 of max |F|
TEST(XcwPartitionTests, RiTshbIsTheRiTable)
{
	const std::string cif = fs::absolute(fixture("P1_test/P1_test.cif")).string(), hkl = fs::absolute(fixture("P1_test/P1_subset.hkl")).string();
	app_dir run("xcw_ri_tshb");
	std::ofstream(run.dir / "settings.txt") << "normal normal_conv params 177 basis_set sto-3g max_iter 100 F charge 0 mult 1 rhf start 0 step_size 0.01 end 0.01";
	auto worst = [&](std::vector<std::string> extra) {
		std::vector<std::string> a{ "-do_XCW", "-no_gpu", "-cif", cif, "-hkl", hkl, "-XCW_settings", "settings.txt", "-xcw_partition", "ri", "-tshb",
			"-tune", "NOSPHERA2_XCW_TSHB_TABLE", "-no_date" };
		a.insert(a.end(), extra.begin(), extra.end());
		EXPECT_EQ(run.run(a), 0);
		double w = 0.0;
		for (const std::string step : { "NA2_0000000", "NA2_0010000" }) {
			double f_max = 0.0, rel = 0.0;
			const double d = compare_with_table(tshb_block::read(run.dir / (step + ".tshb")), itsc_block(run.dir / (step + ".tscb")), cell(cif, std::cout, false), f_max, &rel);
			EXPECT_GT(f_max, 0.0) << step << ".tscb is all zero, max |dF| " << d;
			w = std::max(w, d / f_max);
		}
		return w;
	};
	const double own = worst({});
	std::stringstream scf;
	scf << std::ifstream(run.dir / "SCF.log").rdbuf();
	EXPECT_NE(scf.str().find("the RI partition's fit"), std::string::npos);
	std::cout << "RI .tshb against the RI table: max |dF| / max |F| " << own << std::endl;
	EXPECT_LT(own, 1e-10);
	EXPECT_GT(worst({ "-tune", "NOSPHERA2_XCW_TSHB_FIT" }), 1e-12);
	EXPECT_GT(worst({ "-tune", "NOSPHERA2_XCW_TSCB_BY_F" }), 1e-3);
}

//XCW's RI partition (-xcw_partition ri, AutoAux): the perturbation is the exact transpose of F_calc (NOSPHERA2_XCW_ADJOINT_CHECK,
//random D and weights), unrestrained to round-off and with Hirshfeld moment restraints, whose targets go through the grid, to its
//single-precision blocks; restrained to MBIS the targets come from the converged step and the scan still runs to its end
TEST(XcwPartitionTests, RiPartitionPerturbationIsTheTransposeOfFcalc)
{
	const std::string cif = fs::absolute(fixture("P1_test/P1_test.cif")).string(), hkl = fs::absolute(fixture("P1_test/P1_subset.hkl")).string();
	app_dir run("xcw_ri");
	std::ofstream(run.dir / "settings.txt") << "normal normal_conv params 177 basis_set sto-3g max_iter 100 F charge 0 mult 1 rhf start 0 step_size 0.01 end 0.01";
	auto check = [&](std::vector<std::string> extra, const double tol)
	{
		std::vector<std::string> a{ "-do_XCW", "-no_gpu", "-cif", cif, "-hkl", hkl, "-XCW_settings", "settings.txt", "-xcw_partition", "ri",
			"-tune", "NOSPHERA2_XCW_ADJOINT_CHECK=1", "-no_date" };
		a.insert(a.end(), extra.begin(), extra.end());
		std::stringstream out;
		std::streambuf* old = std::cout.rdbuf(out.rdbuf());
		const int rc = run.run(a);
		std::cout.rdbuf(old);
		EXPECT_EQ(rc, 0);
		const std::string log = out.str() + run.log();
		EXPECT_NE(log.find("XCW partition: RI fit"), std::string::npos) << log;
		const size_t at = log.find("XCW adjoint check (ri): ");
		EXPECT_NE(at, std::string::npos);
		if (at != std::string::npos)
			EXPECT_LT(std::stod(log.substr(at + 24)), tol) << log.substr(at, 60);
		return log;
	};
	const std::string plain = check({ "-tune", "NOSPHERA2_XCW_RI_PACK_CHECK" }, 1e-10);
	const size_t pk = plain.find("XCW RI packed store check: max |dg| / max |g| ");
	ASSERT_NE(pk, std::string::npos);
	EXPECT_LT(std::stod(plain.substr(pk + 46)), 1e-12);
	check({ "-tune", "NOSPHERA2_XCW_RI_FULL" }, 1e-10);
	check({ "-multipole_moments", "hirshfeld", "2" }, 1e-5);
	const std::string mbis = check({ "-multipole_moments", "mbis", "1" }, 1e-5);
	EXPECT_NE(mbis.find("MBIS shares:"), std::string::npos);
	EXPECT_TRUE(fs::exists(run.dir / "NA2_0010000.tscb"));
	std::stringstream scf;
	scf << std::ifstream(run.dir / "SCF.log").rdbuf();
	EXPECT_NE(scf.str().find("shares refit to the converged density"), std::string::npos);
	for (const std::string s : { "mulliken", "nuclear", "sanderson" })
		check({ "-multipole_moments", s }, 1e-10);
}

//TFVC, MBIS and EMBIS shares on XCW's grid follow the density (xcw_weights.cpp): on the P1 subset R1(gt) at lambda 0 is within
//10 % of Hirshfeld's 0.0253 (P1_test_XCW.good) and falls at 0.05, the perturbation is the transpose of F_calc at the held shares,
//and at the converged lambda 0 density the charges are the stand-alone ones on GridManager's grid of the written wfn: TFVC as
//it gives them, MBIS and EMBIS fitted there to 1e-8 (the stand-alone 1e-3 stop leaves them up to 5e-2 e off)
static void density_following_shares(const std::string& scheme)
{
	const std::string cif = fs::absolute(fixture("P1_test/P1_test.cif")).string(), hkl = fs::absolute(fixture("P1_test/P1_subset.hkl")).string();
	app_dir run("xcw_" + scheme);
	std::ofstream(run.dir / "settings.txt") << "normal normal_conv params 177 basis_set sto-3g max_iter 100 F charge 0 mult 1 rhf start 0 step_size 0.05 end 0.05";
	auto xcw = [&](std::vector<std::string> extra) {
		std::vector<std::string> a{ "-do_XCW", "-no_gpu", "-cif", cif, "-hkl", hkl, "-XCW_settings", "settings.txt", "-" + scheme, "-no_date" };
		a.insert(a.end(), extra.begin(), extra.end());
		std::stringstream out;
		std::streambuf* old = std::cout.rdbuf(out.rdbuf());
		const int rc = run.run(a);
		std::cout.rdbuf(old);
		EXPECT_EQ(rc, 0);
		return out.str() + run.log();
	};
	//the round-10 path: plain fixed-point fits to 1e-6, the Becke and TFVC pass, refits between solves only
	const std::string old_log = xcw({ "-tune", "NOSPHERA2_XCW_TFVC_OLD", "-tune", "NOSPHERA2_XCW_WEIGHT_MIX=0", "-tune", "NOSPHERA2_XCW_WEIGHT_DTOL=0" });
	const std::string log = xcw({ "-tune", "NOSPHERA2_XCW_ADJOINT_CHECK=1" });
	auto row = [](const std::string& text, const std::string& lambda, const int col) {
		const size_t at = text.find("\n" + lambda);
		EXPECT_NE(at, std::string::npos) << lambda;
		std::istringstream l(text.substr(at + 1));
		double v[4] = { 0, 0, 0, 1 };
		for (int i = 0; i < 4; i++) l >> v[i];
		return v[col];
	};
	auto r1 = [&](const std::string& lambda) { return row(log, lambda, 3); };
	//the converged result is the same fixed point within the SCF and refit tolerances: GooF to 1e-4, so its square, the
	//criterion (about 20 here), to 1e-3, R1 to 1e-5
	for (const std::string l : { " 0.00000\t", " 0.05000\t" }) {
		EXPECT_NEAR(row(log, l, 1), row(old_log, l, 1), 1e-3) << l;
		EXPECT_NEAR(row(log, l, 3), row(old_log, l, 3), 1e-5) << l;
	}
	const double r0 = r1(" 0.00000\t"), r5 = r1(" 0.05000\t");
	EXPECT_LT(std::abs(r0 - 0.0253), 0.1 * 0.0253) << r0;
	EXPECT_LT(r5, r0);
	const size_t at = log.find("XCW adjoint check (grid): ");
	ASSERT_NE(at, std::string::npos);
	EXPECT_LT(std::stod(log.substr(at + 26)), 1e-5);
	std::stringstream scf;
	scf << std::ifstream(run.dir / "SCF.log").rdbuf();
	const size_t c = scf.str().find("XCW grid charges at lambda 0.00000:");
	ASSERT_NE(c, std::string::npos);
	std::istringstream line(scf.str().substr(c + 35, scf.str().find('\n', c) - c - 35));
	vec q;
	std::string label;
	double v;
	while (line >> label >> v) q.push_back(v);
	WFN w(run.dir / "NA2_0000000.wfn");
	ASSERT_EQ(static_cast<int>(q.size()), w.get_ncen());
	GridConfiguration config;
	config.partition_type = scheme == "tfvc" ? PartitionType::TFVC : scheme == "mbis" ? PartitionType::MBIS : PartitionType::EMBIS;
	GridManager gm(config);
	ivec atoms(w.get_ncen());
	for (int a = 0; a < w.get_ncen(); a++) atoms[a] = a;
	const cell uc(cif, std::cout, false);
	gm.setup3DGridsForMolecule(w, atoms, {}, uc);
	vec ref(w.get_ncen());
	if (scheme == "tfvc") {
		const PartitionResults p = gm.calculatePartitionedCharges(w, uc);
		for (int a = 0; a < w.get_ncen(); a++) ref[a] = w.get_atom_charge(a) - p.atom_charges[PartitionResults::CHARGE_ORDER::S_TFVC][a];
	}
	else {
		GridData& G = gm.getGridData();
		std::ostream none(nullptr);
		const auto m = make_MBIS_vectors(w, G.helper_grids, G.helper_num_points_per_atom, false, none, {}, 1e-8);
		//Anderson mixing reaches the same fixed point in fewer sweeps
		int plain = 0, mixed = 0;
		const auto m0 = make_MBIS_vectors(w, G.helper_grids, G.helper_num_points_per_atom, false, none, {}, 1e-8, &plain);
		const auto m6 = make_MBIS_vectors(w, G.helper_grids, G.helper_num_points_per_atom, false, none, {}, 1e-8, &mixed, 6);
		for (int a = 0; a < w.get_ncen(); a++)
			for (int i = 0; i < static_cast<int>(m0[a].second.size()); i++) {
				EXPECT_NEAR(m0[a].second[i], m6[a].second[i], 1e-6);
				EXPECT_NEAR(m0[a].first[i], m6[a].first[i], 1e-6);
			}
		EXPECT_LT(mixed, plain);
		const auto e = scheme == "embis" ? make_EMBIS_tensors(w, G.helper_grids, G.helper_num_points_per_atom, false, m, none, {}, 1e-8) : std::vector<std::pair<vec2, vec>>();
		for (int a = 0; a < w.get_ncen(); a++) {
			ref[a] = w.get_atom_charge(a);
			const vec& pop = scheme == "embis" ? e[a].second : m[a].second;
			for (int i = 0; i < static_cast<int>(pop.size()); i++) ref[a] -= pop[i];
		}
	}
	double worst = 0;
	for (int a = 0; a < w.get_ncen(); a++) worst = std::max(worst, std::abs(q[a] - ref[a]));
	std::cout << scheme << ": R1(gt) " << r0 << " -> " << r5 << ", max |q_XCW - q_stand-alone| " << worst << std::endl;
	EXPECT_LT(worst, 2e-4);
}

TEST(XcwPartitionTests, TfvcSharesFollowTheDensity) { density_following_shares("tfvc"); }

//RI restrained to TFVC on the P1 subset through lambda 0.03: N8-H8's line density loses its minimum on the way, and with
//make_chi's fallback chi = 1 the shares flipped between solves until the 20-solve cap (NOSPHERA2_XCW_CHI_FALLBACK); with
//chi at the inflection every step settles without a second solve
TEST(XcwPartitionTests, RiTfvcRefitsSettle)
{
	const std::string cif = fs::absolute(fixture("P1_test/P1_test.cif")).string(), hkl = fs::absolute(fixture("P1_test/P1_subset.hkl")).string();
	app_dir run("xcw_ri_tfvc");
	std::ofstream(run.dir / "settings.txt") << "normal normal_conv params 177 basis_set sto-3g max_iter 100 F charge 0 mult 1 rhf start 0 step_size 0.01 end 0.03";
	std::stringstream out;
	std::streambuf* old = std::cout.rdbuf(out.rdbuf());
	const int rc = run.run({ "-do_XCW", "-no_gpu", "-cif", cif, "-hkl", hkl, "-XCW_settings", "settings.txt", "-xcw_partition", "ri", "-multipole_moments", "tfvc", "2", "-no_date" });
	std::cout.rdbuf(old);
	ASSERT_EQ(rc, 0);
	const std::string log = out.str() + run.log();
	int steps = 0;
	for (size_t at = log.find("TFVC shares:"); at != std::string::npos; at = log.find("TFVC shares:", at + 1)) {
		const size_t e = log.find('\n', at);
		const std::string line = log.substr(at, e - at);
		if (line.find("more solves") == std::string::npos) continue;
		steps++;
		EXPECT_NE(line.find(", 0 more solves"), std::string::npos) << line;
	}
	EXPECT_EQ(steps, 3);
}

//get_tfvc_weight is get_integration_weights' TFVC share, the same operations: equal to round-off on random points and chi
TEST(XcwPartitionTests, TfvcShareKernelIsTheStandAloneOne)
{
	std::mt19937 gen(5);
	std::uniform_real_distribution<double> u(-4.0, 4.0), c(0.5, 2.0);
	const int n = 12;
	vec x(n), y(n), z(n), chi(n * n), dab(n * n), pa(n), pb(n), dist(n);
	ivec Z(n);
	for (int a = 0; a < n; a++) { x[a] = u(gen); y[a] = u(gen); z[a] = u(gen); Z[a] = 1 + a % 8; }
	for (int a = 0; a < n; a++)
		for (int b = 0; b < n; b++) {
			dab[a * n + b] = std::hypot(x[b] - x[a], y[b] - y[a], z[b] - z[a]);
			if (b > a) { chi[a * n + b] = c(gen); chi[b * n + a] = 1.0 / chi[a * n + b]; }
		}
	double worst = 0;
	for (int p = 0; p < 2000; p++) {
		const int g = p % n;
		const double px = x[g] + u(gen), py = y[g] + u(gen), pz = z[g] + u(gen);
		const double ref = get_integration_weights(n, Z.data(), x.data(), y.data(), z.data(), g, px, py, pz, pb, pa, chi, dab.data())[1];
		worst = std::max(worst, std::abs(ref - get_tfvc_weight(n, x.data(), y.data(), z.data(), g, px, py, pz, pa, dist, chi, dab.data())));
	}
	EXPECT_LT(worst, 1e-14);
}
TEST(XcwPartitionTests, MbisSharesFollowTheDensity) { density_following_shares("mbis"); }
TEST(XcwPartitionTests, EmbisSharesFollowTheDensity) { density_following_shares("embis"); }
