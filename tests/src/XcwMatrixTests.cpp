#include "pch.h"
#include <gtest/gtest.h>

#include "core/NoSpherA2.h"

#include <cmath>
#include <fstream>
#include <regex>
#include <sstream>

//The XCW permutation matrix (tests/xcw_matrix): NH3 with 22 parameters, three lambdas each, one row per pairwise
//combination of weight scheme, extinction, anharmonic ADPs, grown cluster, F/F^2 and target, plus all-default (D) and
//all-on (ON). Every row runs in process on one thread and is compared against matrix.good (written by
//run_cases.py --cases matrix --golden) on final lambda, R1, wR2, GooF, scale and extinction - never on time
namespace
{
	namespace fs = std::filesystem;

	struct matrix_row
	{
		std::string id, cif, settings, ext;
		bool grown = false;
		double lambda = 0, r1 = 0, wr2 = 0, goof = 0, scale = 0;
	};

	//ctest lists the row by its id, not as a byte dump
	void PrintTo(const matrix_row& r, std::ostream* os) { *os << r.id; }

	fs::path matrix_dir() { return nos_test_repo_root() / "tests" / "xcw_matrix"; }

	std::vector<matrix_row> matrix_rows()
	{
		std::vector<matrix_row> rows;
		std::ifstream in(matrix_dir() / "matrix.good");
		for (std::string line; std::getline(in, line);)
		{
			if (line.empty() || line[0] == '#') continue;
			std::vector<std::string> f;
			std::stringstream s(line);
			for (std::string t; std::getline(s, t, '\t');) f.push_back(t);
			if (f.size() != 10) continue;
			matrix_row r;
			r.id = f[0], r.cif = f[1], r.grown = f[2] == "1", r.ext = f[8], r.settings = f[9];
			r.lambda = std::stod(f[3]), r.r1 = std::stod(f[4]), r.wr2 = std::stod(f[5]), r.goof = std::stod(f[6]), r.scale = std::stod(f[7]);
			rows.push_back(r);
		}
		return rows;
	}

	std::string text(const fs::path& p)
	{
		std::stringstream s;
		s << std::ifstream(p).rdbuf();
		return s.str();
	}

	void expect_rel(const double got, const double want, const double tol, const char* what)
	{
		EXPECT_LE(std::abs(got - want), tol * std::max(std::abs(want), 1e-12)) << what << " " << got << " vs " << want;
	}
}

class XcwMatrixTests : public ::testing::TestWithParam<matrix_row> {};

TEST_P(XcwMatrixTests, MatchesGolden)
{
	const matrix_row& row = GetParam();
	const fs::path dir = fs::temp_directory_path() / ("nos_xcw_matrix_" + row.id + "_" +
		std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
	fs::create_directories(dir);
	for (const std::string f : { std::string("nh3.hkl"), row.cif, std::string("nh3_grown.xyz") })
		fs::copy_file(matrix_dir() / f, dir / f);
	std::ofstream(dir / "s.txt") << row.settings;
	const fs::path previous = fs::current_path();
	fs::current_path(dir);
	std::vector<std::string> args = { "NoSpherA2", "-do_XCW", "-no_gpu", "-no_date", "-print_F", "-cpus", "1", "-cif", row.cif,
		"-hkl", "nh3.hkl", "-XCW_settings", "s.txt" };
	if (row.grown) args.insert(args.end(), { "-xyz", "nh3_grown.xyz" });
	std::vector<char*> argv;
	for (std::string& a : args) argv.push_back(a.data());
	const int rc = run_app(static_cast<int>(argv.size()), argv.data());
	fs::current_path(previous);
	const std::string log = text(dir / "NoSpherA2.log"), scf = text(dir / "SCF.log");
	std::vector<fs::path> fcalc;
	for (const auto& e : fs::directory_iterator(dir))
		if (e.path().filename().string().rfind("NA2_", 0) == 0 && e.path().filename().string().find("_Fcalc.txt") != std::string::npos)
			fcalc.push_back(e.path());
	std::sort(fcalc.begin(), fcalc.end());
	const std::string first_fcalc_line = fcalc.empty() ? "" : text(fcalc.back()).substr(0, text(fcalc.back()).find('\n'));
	std::error_code ec;
	fs::remove_all(dir, ec);
	ASSERT_EQ(rc, 0);

	std::smatch m;
	ASSERT_TRUE(std::regex_search(log, m, std::regex("XCW final model: ([^\\n]*)"))) << "no final-model line";
	std::map<std::string, std::string> kv;
	std::stringstream s(m[1].str());
	for (std::string k, v; s >> k >> v;) kv[k] = v;
	EXPECT_NEAR(std::stod(kv.at("lambda")), row.lambda, 1e-6);
	EXPECT_NEAR(std::stod(kv.at("R1")), row.r1, 2e-4);
	EXPECT_NEAR(std::stod(kv.at("wR2")), row.wr2, 2e-4);
	expect_rel(std::stod(kv.at("GooF")), row.goof, 2e-3, "GooF");

	ASSERT_NE(first_fcalc_line.find("scale ="), std::string::npos) << "no Fcalc file";
	expect_rel(std::stod(first_fcalc_line.substr(first_fcalc_line.find("scale =") + 7)), row.scale, 2e-3, "scale");

	//the last lambda's extinction line in SCF.log, or none at all for "-"
	const std::regex ext_line("lambda [\\d.]+: extinction[^\\n]*\\s([-\\d.eE+]+)[ \\t]*(\\n|$)");
	std::string ext;
	for (std::sregex_iterator i(scf.begin(), scf.end(), ext_line), end; i != end; ++i) ext = (*i)[1].str();
	if (row.ext == "-")
		EXPECT_TRUE(ext.empty()) << "unexpected extinction " << ext;
	else
	{
		ASSERT_FALSE(ext.empty()) << "no extinction line";
		expect_rel(std::stod(ext), std::stod(row.ext), 2e-3, "extinction");
	}
}

INSTANTIATE_TEST_SUITE_P(Rows, XcwMatrixTests, ::testing::ValuesIn(matrix_rows()),
	[](const ::testing::TestParamInfo<matrix_row>& info) {
		std::string n = "m" + info.param.id;
		for (char& c : n)
			if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';
		return n;
	});
