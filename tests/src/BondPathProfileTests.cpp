#include "pch.h"
#include "core/base/convenience.h"
#include "core/wfn/wfn_class.h"
#include "core/analysis/bondwise_analysis.h"
#include "core/analysis/topology.h"
#include "core/analysis/b2c.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
	std::vector<std::string> read_lines(const std::filesystem::path& p)
	{
		std::ifstream in(p);
		std::vector<std::string> lines;
		for (std::string l; std::getline(in, l);) lines.push_back(l);
		return lines;
	}
}

//One file per path named after its nuclei, a second path between the same pair numbered _2, an
//unassigned path (a < 0) skipped; s is the arc length and rho, the Laplacian are the density's own
TEST(BondPathProfile, FilesHoldTheDensityAlongEachPath)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wavy(f);
	ASSERT_EQ(wavy.get_ncen(), 2);
	const d3 o = wavy.get_atom_pos(0), h = wavy.get_atom_pos(1);
	topology::bond_path bp;
	for (int k = 0; k <= 4; k++)
		bp.points.push_back(d3{ o[0] + k * (h[0] - o[0]) / 4, o[1] + k * (h[1] - o[1]) / 4, o[2] + k * (h[2] - o[2]) / 4 });
	bp.bcp_index = 3;
	bp.a = 0;
	bp.b = 1;
	bp.bpl = bp.gpl = array_length(o, h);
	topology::bond_path off = bp, lost = bp;
	off.gap_b = 0.25;
	lost.a = -1;
	const std::filesystem::path dir = std::filesystem::temp_directory_path() / "bond_path_profile_test";
	std::filesystem::create_directories(dir);
	const std::string stem = (dir / "oh").string();
	std::ostringstream log;
	write_bond_path_profiles(wavy, { bp, lost, off }, stem, &log);
	const std::filesystem::path first = stem + "_bond_path_0_1.dat", second = stem + "_bond_path_0_1_2.dat";
	ASSERT_TRUE(std::filesystem::exists(first));
	ASSERT_TRUE(std::filesystem::exists(second));
	EXPECT_FALSE(std::filesystem::exists(stem + "_bond_path_-1_1.dat"));
	const std::vector<std::string> lines = read_lines(first);
	ASSERT_EQ(lines.size(), bp.points.size() + 1);
	EXPECT_EQ(lines[0], "# s x y z rho laplacian, bohr and a.u.; the bond critical point is point 4");
	double s = 0.0;
	for (size_t k = 0; k < bp.points.size(); k++) {
		if (k) s += array_length(bp.points[k], bp.points[k - 1]);
		d3 grad{ 0.0, 0.0, 0.0 };
		double lap = 0.0;
		const double rho = calculate_density(wavy, bp.points[k], grad, lap);
		std::ostringstream want;
		want << std::scientific << std::setprecision(6) << std::setw(14) << s;
		for (int c = 0; c < 3; c++) want << std::setw(14) << bp.points[k][c];
		want << std::setw(14) << rho << std::setw(14) << lap;
		EXPECT_EQ(lines[k + 1], want.str()) << "point " << k;
	}
	EXPECT_EQ(read_lines(second), lines) << "the second path is the same line, only its log entry differs";
	const std::string text = log.str();
	EXPECT_EQ(std::count(text.begin(), text.end(), '\n'), 2) << text;
	EXPECT_NE(text.find("Bond path 0 (1) - 1 (8): BPL "), std::string::npos) << text;
	EXPECT_EQ(text.find("end maxima"), text.rfind("end maxima")) << "only the second path has an end off its nucleus";
	EXPECT_NE(text.find(", end maxima 0.0000 / 0.2500 bohr off the nuclei"), std::string::npos) << text;
	std::filesystem::remove_all(dir);
}

//digest_options starts every run from basin_flags_reset, so a flag set by one in-process run must not reach the next
TEST(BasinFlags, ResetRestoresTheDefaults)
{
	beta_spheres_set_enabled(false);
	basin_adaptive_step_set_enabled(true);
	basin_step_scale_set(0.5);
	basin_timing_set_enabled(true);
	basin_gpu_set_mode(0);
	basin_flags_reset();
	EXPECT_TRUE(beta_spheres_enabled());
	EXPECT_FALSE(basin_adaptive_step_enabled());
	EXPECT_EQ(basin_step_scale(), 1.0);
	EXPECT_FALSE(basin_timing_enabled());
}
