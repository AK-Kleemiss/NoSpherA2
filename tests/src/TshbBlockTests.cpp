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

#include <complex>
#include <fstream>
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
	double compare_with_table(const tshb_block& tshb, const itsc_block& table, const cell& unit_cell, double& f_max)
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
				worst = std::max(worst, std::abs(record_form_factor(tshb.get_record(a), k) - sf[r]));
				f_max = std::max(f_max, std::abs(sf[r]));
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

//a Thakkar-filled atom as Slater records gives what the fill writes into the table, in every charge state the Slater
//tables carry; the tabulated delta_k and GTO routes refuse
TEST(TshbBlockTests, FillSlaterRecordsReproduceTheSphericalFormFactor)
{
	for (const auto& [Z, q] : std::vector<std::pair<int, double>>{ { 1, 0.0 }, { 6, 0.0 }, { 8, -0.45 }, { 11, 0.7 }, { 17, -0.6 }, { 26, 0.9 }, { 7, 0.3 } })
	{
		const HE_Spherical_Atom atom(Z, q);
		std::vector<std::array<double, 3>> terms;
		ASSERT_TRUE(atom.slater_terms(terms)) << Z << " " << q;
		tshb_record r;
		for (const auto& t : terms) r.slaters.push_back(tshb_slater(t));
		double worst = 0.0;
		for (const double stl : { 0.0, 0.05, 0.2, 0.5, 1.0, 1.5 })
		{
			const double k_bohr = constants::bohr2ang(constants::FOUR_PI * stl); //the fill's k
			const double k[3] = { 0.0, 0.0, constants::ang2bohr(k_bohr) };
			worst = std::max(worst, std::abs(record_form_factor(r, k) - atom.get_form_factor(k_bohr)));
		}
		EXPECT_LT(worst, 1e-10) << "Z " << Z << " q " << q;
	}
	std::vector<std::array<double, 3>> terms;
	EXPECT_FALSE(HE_Spherical_Atom(55, 0.8).slater_terms(terms)) << "Cs takes the GTO tables";
	EXPECT_TRUE(terms.empty());
	EXPECT_FALSE(HE_Spherical_Atom(12, 1.6).slater_terms(terms)) << "Mg2+ walks the delta_k series";
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
	app_dir run("merge");
	for (size_t i = 0; i < runs.size(); i++)
	{
		runs[i].push_back("-tshb");
		ASSERT_EQ(run.run(runs[i]), 0) << i;
		fs::rename(run.dir / "experimental.tshb", run.dir / (i == 0 ? std::string("joint.tshb") : "f" + std::to_string(i) + ".tshb"));
	}
	ASSERT_EQ(run.run({ "-merge", "f1.tshb", "f2.tshb", "f3.tshb", "f4.tshb" }), 0);
	ASSERT_TRUE(fs::exists(run.dir / "combined.tshb"));
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

//the same structure with Rh in place of Co: Rh takes the tabulated GTO ion data, so the run writes the table instead
TEST(TshbJointTests, AGtoFillAtomKeepsTheTableAndSaysWhy)
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
	ASSERT_EQ(run.run({ "-SALTED", model, "-cif", "rh.cif", "-wfn", "rh.xyz", "-dmin", "0.8", "-tshb", "-no_date" }), 0);
	const std::string log = run.log();
	EXPECT_NE(log.find("-tshb: the table is written instead, Rh is filled from the tabulated GTO ion data"), std::string::npos) << log;
	EXPECT_TRUE(fs::exists(run.dir / "experimental.tscb"));
	EXPECT_FALSE(fs::exists(run.dir / "experimental.tshb"));
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
