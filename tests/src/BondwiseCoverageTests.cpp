//Coverage for the uncovered branches of Src/core/bondwise_analysis.cpp: the do_bonds box
//geometry in every size mode, the property cubes of a debug run, the Roby table row of a
//homonuclear bond, the wavefunction mode of run_QTAIM_ELI_mask and the debug listing of
//ELI_analysis. Every expectation is derived by hand from the input.example spec and from the
//closed forms of Gaussian densities; BondwiseTests.cpp keeps the parser and golden cases.
#include "pch.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "core/convenience.h"
#include "core/wfn_class.h"
#include "core/cube.h"
#include "core/atoms.h"
#include "core/constants.h"
#include "core/bondwise_analysis.h"

namespace
{
	struct CoutCapture
	{
		std::ostringstream buffer;
		std::streambuf* old;
		CoutCapture() : old(std::cout.rdbuf(buffer.rdbuf())) {}
		~CoutCapture() { std::cout.rdbuf(old); }
		std::string str() const { return buffer.str(); }
	};

	//a scratch directory named after the test; the cwd goes back and the directory is removed when the test passed
	struct Scratch
	{
		std::filesystem::path dir;
		std::filesystem::path old_cwd;
		explicit Scratch(const std::string& name) : old_cwd(std::filesystem::current_path())
		{
			dir = std::filesystem::temp_directory_path() / ("NoSpherA2_BondwiseCoverage_" + name);
			std::filesystem::remove_all(dir);
			std::filesystem::create_directories(dir);
		}
		~Scratch()
		{
			std::error_code ec;
			std::filesystem::current_path(old_cwd, ec);
			if (!::testing::Test::HasFailure())
				std::filesystem::remove_all(dir, ec);
		}
	};

	//C at the origin, N two bohr up the z axis, O one bohr off it; one s Gaussian rho = A exp(-2 r^2) on C.
	//do_bonds axes for atoms 1 2 3: x = (-1,0,0), y = (0,1,0), z = (0,0,-1)
	WFN three_atom_wfn()
	{
		WFN wavy(e_origin::NOT_YET_DEFINED);
		wavy.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		wavy.push_back_atom("N", 0.0, 0.0, 2.0, 7);
		wavy.push_back_atom("O", 1.0, 0.0, 2.0, 8);
		wavy.push_back_MO(1, 2.0, -0.5);
		std::vector<primitive> prims;
		prims.emplace_back(1, 1, 1.0, 1.0);
		wavy.push_back_spherical_shell(0, 0, vec2{ vec{ 1.0 } }, prims, 0, 1);
		wavy.set_exp_cutoff();
		return wavy;
	}

	//the same three atoms with two doubly occupied s orbitals on C, exp(-r^2 / 2) and exp(-r^2 / 4), so
	//rho = 2 exp(-r^2) + 2 exp(-r^2 / 2). Two orbitals because the ELI-D of one alone diverges (its pair
	//volume rho tau - |grad rho|^2 / 4 is identically zero); exponents below 0.92 because the primitive
	//cutoff set from density_accuracy drops exp(-a r^2) beyond a r^2 = 9.9 and the 5^3 box of the
	//property test reaches r^2 = 6
	WFN two_shell_wfn()
	{
		WFN wavy(e_origin::NOT_YET_DEFINED);
		wavy.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		wavy.push_back_atom("N", 0.0, 0.0, 2.0, 7);
		wavy.push_back_atom("O", 1.0, 0.0, 2.0, 8);
		wavy.push_back_MO(1, 2.0, -0.5);
		wavy.push_back_MO(2, 2.0, -0.4);
		std::vector<primitive> prims;
		prims.emplace_back(1, 1, 0.5, 1.0);
		prims.emplace_back(1, 1, 0.25, 1.0);
		for (int mo = 0; mo < 2; mo++)
			for (int p = 0; p < 2; p++)
				wavy.push_back_spherical_shell(mo, 0, vec2{ vec{ mo == p ? 1.0 : 0.0 } }, prims, p, 1);
		wavy.set_exp_cutoff();
		return wavy;
	}

	//H at (-1,0,0) and (+1,0,0) bohr, one s Gaussian (exponent 1) each, the bonding and the antibonding
	//combination both doubly occupied: rho = 4 (exp(-2 r_a^2) + exp(-2 r_b^2)), two atomic Gaussians, so
	//the density has exactly three critical points (two nuclear maxima, one saddle at the origin) and a
	//finite ELI-D everywhere (a single orbital has no Pauli kinetic energy and ELI-D diverges)
	WFN h2_wfn()
	{
		WFN wavy(e_origin::NOT_YET_DEFINED);
		wavy.push_back_atom("H", -1.0, 0.0, 0.0, 1);
		wavy.push_back_atom("H", 1.0, 0.0, 0.0, 1);
		wavy.push_back_MO(1, 2.0, -0.5);
		wavy.push_back_MO(2, 2.0, -0.4);
		std::vector<primitive> prims;
		prims.emplace_back(1, 1, 1.0, 1.0);
		prims.emplace_back(2, 1, 1.0, 1.0);
		for (int mo = 0; mo < 2; mo++)
		{
			wavy.push_back_spherical_shell(mo, 0, vec2{ vec{ 1.0 } }, prims, 0, 1);
			wavy.push_back_spherical_shell(mo, 0, vec2{ vec{ mo == 0 ? 1.0 : -1.0 } }, prims, 1, 1);
		}
		wavy.set_exp_cutoff();
		return wavy;
	}

	//the H2 model with the antibonding orbital built from exponent-2 Gaussians h instead: MO1 = g_a + g_b,
	//MO2 = h_a - h_b, rho = 2 (g_a + g_b)^2 + 2 (h_a - h_b)^2, still mirror symmetric so the QTAIM basins split
	//at x = 0. The primitive cutoff set from density_accuracy drops a Gaussian beyond a r^2 = 9.9; in h2_wfn
	//both orbitals collapse onto g_a beyond 3.15 bohr from the far hydrogen and the ELI-D pair volume
	//rho tau - |grad rho|^2 / 4 vanishes there, here it stays 16 r_a^2 g_a^2 h_a^2 > 0, so the ELI-D is finite
	//and positive on every evaluated point
	WFN h2_two_shell_wfn()
	{
		WFN wavy(e_origin::NOT_YET_DEFINED);
		wavy.push_back_atom("H", -1.0, 0.0, 0.0, 1);
		wavy.push_back_atom("H", 1.0, 0.0, 0.0, 1);
		wavy.push_back_MO(1, 2.0, -0.5);
		wavy.push_back_MO(2, 2.0, -0.4);
		std::vector<primitive> prims;
		prims.emplace_back(1, 1, 1.0, 1.0);
		prims.emplace_back(2, 1, 1.0, 1.0);
		prims.emplace_back(1, 1, 2.0, 1.0);
		prims.emplace_back(2, 1, 2.0, 1.0);
		const double coef[2][4] = { { 1.0, 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0, -1.0 } };
		for (int mo = 0; mo < 2; mo++)
			for (int p = 0; p < 4; p++)
				wavy.push_back_spherical_shell(mo, 0, vec2{ vec{ coef[mo][p] } }, prims, p, 1);
		wavy.set_exp_cutoff();
		return wavy;
	}

	cube read_cube(const std::filesystem::path& p)
	{
		WFN dummy;
		std::ostringstream log;
		return cube(p, true, dummy, log);
	}

	double axis_length(const cube& c, int axis)
	{
		double s = 0.0;
		for (int j = 0; j < 3; j++)
			s += c.get_vector(axis, j) * c.get_vector(axis, j);
		return std::sqrt(s);
	}

	//the centre of the grid points: origin + sum_i ((np_i - 1) / 2) v_i
	d3 box_centre(const cube& c)
	{
		d3 r{};
		for (int i = 0; i < 3; i++)
		{
			r[i] = c.get_origin(i);
			for (int j = 0; j < 3; j++)
				r[i] += 0.5 * (c.get_size(j) - 1) * c.get_vector(j, i);
		}
		return r;
	}

	void expect_box(const cube& c, const std::array<int, 3>& np, const d3& centre, double tol = 1e-5)
	{
		for (int i = 0; i < 3; i++)
		{
			EXPECT_EQ(c.get_size(i), np[i]) << "axis " << i;
			EXPECT_NEAR(box_centre(c)[i], centre[i], tol) << "centre component " << i;
		}
	}

	//the three axes point along -x, +y and -z of the three-atom wavefunction
	void expect_axes(const cube& c)
	{
		const double n0 = axis_length(c, 0), n1 = axis_length(c, 1), n2 = axis_length(c, 2);
		EXPECT_NEAR(c.get_vector(0, 0) / n0, -1.0, 1e-9);
		EXPECT_NEAR(c.get_vector(0, 1), 0.0, 1e-9);
		EXPECT_NEAR(c.get_vector(0, 2), 0.0, 1e-9);
		EXPECT_NEAR(c.get_vector(1, 1) / n1, 1.0, 1e-9);
		EXPECT_NEAR(c.get_vector(1, 0), 0.0, 1e-9);
		EXPECT_NEAR(c.get_vector(1, 2), 0.0, 1e-9);
		EXPECT_NEAR(c.get_vector(2, 2) / n2, -1.0, 1e-9);
		EXPECT_NEAR(c.get_vector(2, 0), 0.0, 1e-9);
		EXPECT_NEAR(c.get_vector(2, 1), 0.0, 1e-9);
	}

	//runs do_bonds on atoms 1 2 3 of wavy with the cube files under the scratch dir; returns the file prefix
	std::string run_bond(WFN& wavy, const Scratch& s, int mode, bool mleng, bool mres, std::array<double, 3> res, bool cub,
		std::array<double, 3> box, bool debug, bool bohr, int run, bool rho, bool rdg, bool eli, bool lap, std::string* out = nullptr)
	{
		wavy.set_path(s.dir / "three");
		std::string captured;
		bond result;
		{
			CoutCapture cap;
			result = do_bonds(wavy, mode, mleng, mres, res.data(), cub, box.data(), 1, 2, 3, debug, bohr, run, rho, rdg, eli, lap);
			captured = cap.str();
		}
		if (out)
			*out = captured;
		EXPECT_TRUE(result.success) << captured;
		const std::string prefix = (s.dir / "three").generic_string() + "_C1_N2_O3_" + std::to_string(run);
		EXPECT_EQ(result.filename, prefix);
		return prefix;
	}

	//the numbers that follow key on the first line containing it; empty when the line is absent
	vec numbers_after(const std::string& text, const std::string& key)
	{
		vec out;
		const auto at = text.find(key);
		if (at == std::string::npos)
			return out;
		const auto eol = text.find('\n', at);
		std::istringstream in(text.substr(at + key.size(), eol == std::string::npos ? std::string::npos : eol - at - key.size()));
		double v;
		while (in >> v)
			out.push_back(v);
		return out;
	}

	struct CriticalPoint
	{
		std::string type;
		std::string owner;
		bool converged = false;
		d3 pos{};
		double rho = 0.0;
		double lap = 0.0;
		double ellipticity = 0.0;
		d3 eig{};
		std::array<d3, 3> vecs{};
		bool has_vecs = false;
	};

	std::string trim(const std::string& s)
	{
		const auto b = s.find_first_not_of(" \t\r\n");
		if (b == std::string::npos)
			return {};
		const auto e = s.find_last_not_of(" \t\r\n");
		return s.substr(b, e - b + 1);
	}

	//the "Density Critical Points" block ELI_analysis prints, up to the QTAIM tables
	std::vector<CriticalPoint> parse_critical_points(const std::string& out)
	{
		std::vector<CriticalPoint> cps;
		std::istringstream in(out);
		std::string line;
		bool inside = false;
		auto fill3 = [](const std::string& l, const std::string& key, d3& target) {
			const vec v = numbers_after(l, key);
			if (v.size() == 3)
				for (int k = 0; k < 3; k++) target[k] = v[k];
		};
		while (std::getline(in, line))
		{
			if (line.rfind("Density Critical Points", 0) == 0) { inside = true; continue; }
			if (!inside)
				continue;
			if (line.rfind("QTAIM Analysis", 0) == 0)
				break;
			if (line.rfind("  CP ", 0) == 0)
			{
				CriticalPoint cp;
				const auto lb = line.find('['), rb = line.find(']');
				if (lb != std::string::npos && rb != std::string::npos)
					cp.type = trim(line.substr(lb + 1, rb - lb - 1));
				cp.converged = line.find("NOT converged") == std::string::npos;
				const auto c = line.find("converged");
				if (c != std::string::npos)
					cp.owner = trim(line.substr(c + 9));
				cps.push_back(cp);
				continue;
			}
			if (cps.empty())
				continue;
			CriticalPoint& cp = cps.back();
			if (line.rfind("    Position  :", 0) == 0) fill3(line, "Position  :", cp.pos);
			else if (line.rfind("    Rho       :", 0) == 0) { const vec v = numbers_after(line, "Rho       :"); if (!v.empty()) cp.rho = v[0]; }
			else if (line.rfind("    HessRho_EigVals:", 0) == 0) fill3(line, "HessRho_EigVals:", cp.eig);
			else if (line.rfind("    HessRho_EigVecs v1:", 0) == 0) { fill3(line, "v1:", cp.vecs[0]); cp.has_vecs = true; }
			else if (line.rfind("                    v2:", 0) == 0) fill3(line, "v2:", cp.vecs[1]);
			else if (line.rfind("                    v3:", 0) == 0) fill3(line, "v3:", cp.vecs[2]);
			else if (line.rfind("    DelSqRho  :", 0) == 0) { const vec v = numbers_after(line, "DelSqRho  :"); if (!v.empty()) cp.lap = v[0]; }
			else if (line.rfind("    Bond Ellipticity:", 0) == 0) { const vec v = numbers_after(line, "Bond Ellipticity:"); if (!v.empty()) cp.ellipticity = v[0]; }
		}
		return cps;
	}

	//label -> (electrons, charge) of the "QTAIM Analysis (atomic quadrature grids)" table
	std::map<std::string, std::pair<double, double>> parse_qtaim_table(const std::string& out)
	{
		std::map<std::string, std::pair<double, double>> rows;
		const auto at = out.find("QTAIM Analysis (atomic quadrature grids):");
		if (at == std::string::npos)
			return rows;
		std::istringstream in(out.substr(at));
		std::string line;
		std::getline(in, line); //title
		std::getline(in, line); //header
		while (std::getline(in, line))
		{
			if (line.rfind("  total in basins:", 0) == 0)
				break;
			//the label column holds spaces - "NNA near H0", "O0 LP", "H2-Hg0 bond" - so it cannot be read
			//with a single >>: that reads "NNA", then fails on "near" as a double and drops the row without
			//saying so, which is how an NNA row would disappear from a table a test is asserting on.
			std::istringstream row(line);
			std::vector<std::string> tok;
			for (std::string t; row >> t; )
				tok.push_back(t);
			if (tok.size() < 4)
				continue;
			//the label runs from token 1 to the first token that is a number: electrons and charge are the
			//first two numeric columns and four more follow them, so counting from the end would depend on
			//how many columns the table happens to print
			size_t first_num = 1;
			auto is_number = [](const std::string& t) {
				try { size_t used = 0; (void)std::stod(t, &used); return used == t.size(); }
				catch (const std::exception&) { return false; }
			};
			while (first_num < tok.size() && !is_number(tok[first_num]))
				first_num++;
			if (first_num < 2 || first_num + 1 >= tok.size() || !is_number(tok[0]))
				continue; //not a data row
			std::string label = tok[1];
			for (size_t t = 2; t < first_num; t++)
				label += " " + tok[t];
			rows[label] = { std::stod(tok[first_num]), std::stod(tok[first_num + 1]) };
		}
		return rows;
	}

	std::filesystem::path epoxide_fixture()
	{
		const auto p = nos_test_repo_root() / "tests" / "epoxide_gbw" / "epoxide.gbw";
		return std::filesystem::exists(p) ? p : std::filesystem::path{};
	}
}

//mode 2 with grid-point counts and box lengths in every axis (mode_res, mode_leng): np = res, |v| = box / np,
//origin = midpoint(atom1, atom2) - sum_i ((np_i - 1) |v_i| / 2) axis_i = (0,0,1) - (-0.75, 1.25, -1.75)
TEST(BondwiseCoverageBoxTests, ModeTwoPointCountsAndBoxLengths)
{
	Scratch s("ModeTwoPointCountsAndBoxLengths");
	WFN wavy = three_atom_wfn();
	const std::string prefix = run_bond(wavy, s, 2, true, true, { 4.0, 6.0, 8.0 }, false, { 2.0, 3.0, 4.0 }, false, false, 7, true, false, false, false);
	const std::filesystem::path rho_file(prefix + "_rho.cube");
	ASSERT_TRUE(std::filesystem::exists(rho_file));
	const cube c = read_cube(rho_file);
	expect_box(c, { 4, 6, 8 }, { 0.0, 0.0, 1.0 });
	expect_axes(c);
	EXPECT_NEAR(c.get_origin(0), 0.75, 1e-5);
	EXPECT_NEAR(c.get_origin(1), -1.25, 1e-5);
	EXPECT_NEAR(c.get_origin(2), 2.75, 1e-5);
	for (int i = 0; i < 3; i++)
		EXPECT_NEAR(axis_length(c, i), 0.5, 1e-5) << "axis " << i;
	EXPECT_FALSE(std::filesystem::exists(prefix + "_rdg.cube"));
	EXPECT_FALSE(std::filesystem::exists(prefix + "_signed_rho.cube"));
}

//mode 3 with box lengths of zero: every axis takes the default ceil(15 |a3 - a1|) / 10 = ceil(15 sqrt 5) / 10 = 3.4
//and the box is centred on the midpoint of atoms 1 and 3
TEST(BondwiseCoverageBoxTests, ModeThreeDefaultLengthFromAtomOneToThree)
{
	Scratch s("ModeThreeDefaultLengthFromAtomOneToThree");
	WFN wavy = three_atom_wfn();
	const std::string prefix = run_bond(wavy, s, 3, true, true, { 5.0, 5.0, 5.0 }, false, { 0.0, 0.0, 0.0 }, false, false, 1, true, false, false, false);
	const cube c = read_cube(prefix + "_rho.cube");
	expect_box(c, { 5, 5, 5 }, { 0.5, 0.0, 1.0 });
	expect_axes(c);
	for (int i = 0; i < 3; i++)
		EXPECT_NEAR(c.get_size(i) * axis_length(c, i), 3.4, 1e-5) << "axis " << i;
}

//mode 4 with spacings and multipliers (!mode_res, !mode_leng): the box sits on the centroid (1/3, 0, 4/3) and
//the lengths are multiplier x mean centroid distance h = 1.0208445 (x: 2h, z: 1h) or the default
//ceil(30 h) / 10 = 3.1 where the multiplier is zero (y); np = round(length / spacing) + 1 = (5, 7, 3)
TEST(BondwiseCoverageBoxTests, ModeFourSpacingFromMultipliers)
{
	Scratch s("ModeFourSpacingFromMultipliers");
	WFN wavy = three_atom_wfn();
	const std::string prefix = run_bond(wavy, s, 4, false, false, { 0.5, 0.5, 0.5 }, false, { 2.0, 0.0, 1.0 }, false, false, 1, true, false, false, false);
	const cube c = read_cube(prefix + "_rho.cube");
	expect_box(c, { 5, 7, 3 }, { 1.0 / 3.0, 0.0, 4.0 / 3.0 });
	expect_axes(c);
}

//mode 2 with spacings and box lengths (!mode_res, mode_leng): x and z fall back to ceil(15 |a2 - a1|) / 10 = 3.0,
//y takes the given 2.5, so np = round(length / spacing) + 1 = (7, 6, 4) around the bond midpoint
TEST(BondwiseCoverageBoxTests, ModeTwoSpacingWithBoxLengths)
{
	Scratch s("ModeTwoSpacingWithBoxLengths");
	WFN wavy = three_atom_wfn();
	const std::string prefix = run_bond(wavy, s, 2, true, false, { 0.5, 0.5, 1.0 }, false, { 0.0, 2.5, 0.0 }, false, false, 1, true, false, false, false);
	const cube c = read_cube(prefix + "_rho.cube");
	expect_box(c, { 7, 6, 4 }, { 0.0, 0.0, 1.0 });
	expect_axes(c);
}

//mode 3 with point counts and multipliers (mode_res, !mode_leng): lengths h, 2h and the default 3.4 for the
//zero multiplier with h = |a3 - a1| = sqrt 5, np = res = (3, 4, 5), centred between atoms 1 and 3
TEST(BondwiseCoverageBoxTests, ModeThreePointCountsWithMultipliers)
{
	Scratch s("ModeThreePointCountsWithMultipliers");
	WFN wavy = three_atom_wfn();
	const std::string prefix = run_bond(wavy, s, 3, false, true, { 3.0, 4.0, 5.0 }, false, { 1.0, 2.0, 0.0 }, false, false, 1, true, false, false, false);
	const cube c = read_cube(prefix + "_rho.cube");
	expect_box(c, { 3, 4, 5 }, { 0.5, 0.0, 1.0 });
	expect_axes(c);
	const double h = std::sqrt(5.0);
	EXPECT_NEAR(c.get_size(0) * axis_length(c, 0), h, 1e-5);
	EXPECT_NEAR(c.get_size(1) * axis_length(c, 1), 2.0 * h, 1e-5);
	EXPECT_NEAR(c.get_size(2) * axis_length(c, 2), 3.4, 1e-5);
}

//"cube selection 1 = all selections in x-direction will be applied in the y and z direction": with spacings and
//multipliers the three axes must share the x length 3.0 and the spacing, i.e. |v_1| = |v_2| = |v_0|
TEST(BondwiseCoverageBoxTests, CubicSpacingAxesShareOneLength)
{
	Scratch s("CubicSpacingAxesShareOneLength");
	WFN wavy = three_atom_wfn();
	const std::string prefix = run_bond(wavy, s, 2, false, false, { 0.5, 0.5, 0.5 }, true, { 0.0, 0.0, 0.0 }, false, false, 1, true, false, false, false);
	const cube c = read_cube(prefix + "_rho.cube");
	expect_box(c, { 7, 7, 7 }, { 0.0, 0.0, 1.0 });
	for (int i = 1; i < 3; i++)
		EXPECT_NEAR(axis_length(c, i), axis_length(c, 0), 1e-9) << "axis " << i;
}

//"orientation_selection 1 = atom1 centered": a 2 x 2 x 2 box of 4 points per axis (steps of 0.5) around atom 1 at
//the origin has its origin at a1 - 0.75 (x + y + z) = (0.75, -0.75, 0.75)
TEST(BondwiseCoverageBoxTests, ModeOneCentresOnAtomOne)
{
	Scratch s("ModeOneCentresOnAtomOne");
	WFN wavy = three_atom_wfn();
	const std::string prefix = run_bond(wavy, s, 1, true, true, { 4.0, 4.0, 4.0 }, false, { 2.0, 2.0, 2.0 }, false, false, 1, true, false, false, false);
	const cube c = read_cube(prefix + "_rho.cube");
	expect_box(c, { 4, 4, 4 }, { 0.0, 0.0, 0.0 });
	EXPECT_NEAR(c.get_origin(0), 0.75, 1e-5);
	EXPECT_NEAR(c.get_origin(1), -0.75, 1e-5);
	EXPECT_NEAR(c.get_origin(2), 0.75, 1e-5);
}

//"resolution selection 0 = res will contain distance between gridpoints": the written cube vectors must have the
//requested lengths 0.5, 0.5 and 1.0
TEST(BondwiseCoverageBoxTests, SpacingModeUsesTheRequestedSpacing)
{
	Scratch s("SpacingModeUsesTheRequestedSpacing");
	WFN wavy = three_atom_wfn();
	const std::string prefix = run_bond(wavy, s, 2, true, false, { 0.5, 0.5, 1.0 }, false, { 0.0, 2.5, 0.0 }, false, false, 1, true, false, false, false);
	const cube c = read_cube(prefix + "_rho.cube");
	EXPECT_NEAR(axis_length(c, 0), 0.5, 1e-5);
	EXPECT_NEAR(axis_length(c, 1), 0.5, 1e-5);
	EXPECT_NEAR(axis_length(c, 2), 1.0, 1e-5);
}

//with bohr = true the box lengths are Angstrom ("1 = box will contain length in angstrom"): a 2 A box of 4 points
//has vectors of 0.5 A = 0.94486 bohr and stays centred on the bond midpoint (0, 0, 1) bohr
TEST(BondwiseCoverageBoxTests, BohrBoxStaysCentredOnTheBond)
{
	Scratch s("BohrBoxStaysCentredOnTheBond");
	WFN wavy = three_atom_wfn();
	const std::string prefix = run_bond(wavy, s, 2, true, true, { 4.0, 4.0, 4.0 }, false, { 2.0, 2.0, 2.0 }, false, true, 1, true, false, false, false);
	const cube c = read_cube(prefix + "_rho.cube");
	expect_box(c, { 4, 4, 4 }, { 0.0, 0.0, 1.0 });
	for (int i = 0; i < 3; i++)
		EXPECT_NEAR(axis_length(c, i), constants::ang2bohr(0.5), 1e-5) << "axis " << i;
}

//a debug run with rho, rdg, eli and lap on the two-shell density rho = 2 exp(-r^2) + 2 exp(-r^2 / 2): the debug lines
//of do_bonds and compute_dens appear, the signed density is -rho everywhere (rho falls with r, so the tangential
//Hessian eigenvalue rho'(r) / r is negative and it is the middle one; at r = 0 all three equal rho''(0) < 0),
//lap = 2 (4 r^2 - 6) exp(-r^2) + 2 (r^2 - 3) exp(-r^2 / 2), |grad rho| = r (4 exp(-r^2) + 2 exp(-r^2 / 2)),
//rdg = |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3)), and the orbital ELI-D 1/2 rho (48 / g)^(3/8) with the pair volume
//g = rho tau - |grad rho|^2 / 4 = 4 |phi_1 grad phi_2 - phi_2 grad phi_1|^2 = r^2 exp(-3 r^2 / 2); at the grid point
//on the nucleus g = 0 and the code writes 0 for the diverging ELI-D. The primitives are never cut off on this box,
//so rho > 0 everywhere and the rdg = 101 marker of an empty point does not occur
TEST(BondwiseCoveragePropertyTests, DebugRunWritesEveryPropertyCube)
{
	Scratch s("DebugRunWritesEveryPropertyCube");
	WFN wavy = two_shell_wfn();
	std::string out;
	//5 points over 2.5 bohr: an odd count puts one grid point on the nucleus at (0, 0, 0)
	const std::string prefix = run_bond(wavy, s, 2, true, true, { 5.0, 5.0, 5.0 }, false, { 2.5, 2.5, 2.5 }, true, false, 1, true, true, true, true, &out);
	EXPECT_NE(out.find("The Atoms found corresponding to your selection are:"), std::string::npos);
	EXPECT_NE(out.find("mode_leng=true; using boxsize"), std::string::npos);
	EXPECT_NE(out.find("gvector before: "), std::string::npos);
	for (const char* suffix : { "_rho.cube", "_signed_rho.cube", "_rdg.cube", "_eli.cube", "_lap.cube" })
		EXPECT_TRUE(std::filesystem::exists(prefix + suffix)) << suffix;
	const cube rho = read_cube(prefix + "_rho.cube");
	const cube signed_rho = read_cube(prefix + "_signed_rho.cube");
	const cube rdg = read_cube(prefix + "_rdg.cube");
	const cube eli = read_cube(prefix + "_eli.cube");
	const cube lap = read_cube(prefix + "_lap.cube");
	expect_box(rho, { 5, 5, 5 }, { 0.0, 0.0, 1.0 });
	const double kf = std::cbrt(3.0 * constants::PI * constants::PI);
	int checked = 0;
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			for (int k = 0; k < 5; k++)
			{
				const double r = rho.get_value(i, j, k);
				EXPECT_DOUBLE_EQ(signed_rho.get_value(i, j, k), -r);
				const d3 p = rho.get_pos(i, j, k);
				const double r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
				const double ga = std::exp(-r2), gb = std::exp(-0.5 * r2);
				const double expected_rho = 2.0 * ga + 2.0 * gb;
				EXPECT_NEAR(r, expected_rho, 1e-5 * expected_rho);
				const double l = lap.get_value(i, j, k);
				EXPECT_NEAR(l, 2.0 * (4.0 * r2 - 6.0) * ga + 2.0 * (r2 - 3.0) * gb,
					1e-5 * (std::abs(l) + 2.0 * (4.0 * r2 + 6.0) * ga + 2.0 * (r2 + 3.0) * gb) + 1e-12);
				const double expected_rdg = std::sqrt(r2) * (4.0 * ga + 2.0 * gb) / (2.0 * kf * std::pow(expected_rho, 4.0 / 3.0));
				EXPECT_NEAR(rdg.get_value(i, j, k), expected_rdg, 2e-5 * expected_rdg + 1e-6);
				const double e = eli.get_value(i, j, k);
				EXPECT_TRUE(std::isfinite(e));
				if (r2 < 1e-6)
				{
					EXPECT_DOUBLE_EQ(e, 0.0) << "ELI-D on the nucleus";
					continue;
				}
				const double expected_eli = 0.5 * expected_rho * std::pow(48.0 / (r2 * std::exp(-1.5 * r2)), 3.0 / 8.0);
				EXPECT_NEAR(e, expected_eli, 2e-5 * expected_eli);
				checked++;
			}
	EXPECT_EQ(checked, 5 * 5 * 5 - 1);
}

//the epoxide C-C bond has equal nuclear charges, so the row is written by the el_a <= el_b branch: the higher index
//first ("4 -   1"), and the ionic index with its sign flipped. The printed numbers keep their definitions:
//s_AB = n_A + n_B - n_AB, Tot = sqrt(Cov^2 + Ion^2), Pyth = 100 Cov^2 / Tot^2, Arak = 200 |asin(Cov / Tot)| / pi,
//and the near-C2v molecule makes the bond almost homopolar: |Ion| and |n_A - n_B| small. The Arakai identity is
//checked as Tot sin(pi Arak / 200) = Cov: with Cov / Tot within 1e-4 of one the asin of the three-decimal
//table values is undetermined to half a percent, the sine is not
TEST(BondwiseCoverageRobyTests, EpoxideCarbonCarbonRowKeepsTheIdentities)
{
	const auto p = epoxide_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
	std::string out;
	{
		CoutCapture cap;
		WFN wavy(p);
		Roby_information roby(wavy);
		out = cap.str();
	}
	EXPECT_NE(out.find("Roby-Gould Bond Indices (RGBI) Analysis"), std::string::npos);
	const std::string key = "  C -  C";
	const auto at = out.find(key);
	ASSERT_NE(at, std::string::npos) << out;
	EXPECT_EQ(out.find(key, at + 1), std::string::npos) << "one C-C bond expected";
	const auto bol = out.rfind('\n', at) + 1;
	EXPECT_EQ(out.substr(bol, at - bol), "   4 -   1  ");
	const vec n = numbers_after(out, key);
	ASSERT_EQ(n.size(), 9u);
	const double n_A = n[0], n_B = n[1], n_AB = n[2], s_AB = n[3], cov = n[4], ion = n[5], tot = n[6], pyth = n[7], arak = n[8];
	EXPECT_NEAR(s_AB, n_A + n_B - n_AB, 2e-3);
	EXPECT_NEAR(tot, std::hypot(cov, ion), 2e-3);
	ASSERT_GT(tot, 0.1);
	EXPECT_NEAR(pyth, 100.0 * cov * cov / (tot * tot), 0.3);
	EXPECT_NEAR(tot * std::sin(arak * constants::PI / 200.0), cov, 2e-3);
	EXPECT_GT(cov, 0.0);
	EXPECT_LT(std::abs(ion), 0.05);
	EXPECT_LT(std::abs(n_A - n_B), 0.05);
}

//wavefunction mode: the two-shell H2 model is written as a .wfn, read back and gridded with radius 1.1 A and 0.5 A steps,
//which is 8 x 6 x 6 points from (-1 - r, -r, -r) with r = ang2bohr(1.1) = 2.0787 bohr and steps (2 + 2r) / 8 and
//2r / 6 - the counts are even because readxyzMinMax_fromWFN forces them to be, which is what puts x index 4 exactly
//on the H-H midplane, the mirror plane of this molecule. Only points strictly inside r of a nucleus are evaluated, so
//the y = z = -r planes and the x = -1 - r plane stay empty and the basin of H0 covers the x columns 1..4 (the last of
//them being that midplane, which the ascent hands to H0) and the y, z rows 1..5: a 4 x 5 x 5 cube from
//(-1 - r + (2 + 2r) / 8, -2r/3, -2r/3). Of its 100 voxels eight carry the background: the four corners at x = -2.3090
//lie 2.3568 bohr from H0, and the four at x = 0 lie 2.2002 bohr from both nuclei, all outside r; the other 92 hold
//the finite, positive ELI-D (the x = -2.31 column lies 3.85 bohr from H1, which is why the model needs its second shell)
TEST(BondwiseCoverageMaskTests, WfnModeMasksTheFirstHydrogenBasin)
{
	Scratch s("WfnModeMasksTheFirstHydrogenBasin");
	const auto wfn_path = s.dir / "h2.wfn";
	{
		WFN model = h2_two_shell_wfn();
		ASSERT_TRUE(model.write_wfn(wfn_path, false, false));
	}
	options opt;
	opt.properties.radius = 1.1;
	opt.properties.resolution = 0.5;
	std::ostringstream log;
	{
		CoutCapture cap;
		run_QTAIM_ELI_mask(wfn_path, {}, { 0 }, -1.0, opt, log);
	}
	EXPECT_NE(log.str().find("Loading wavefunction: "), std::string::npos) << log.str();
	EXPECT_NE(log.str().find("Calculating density and ELI grid (8 x 6 x 6)"), std::string::npos) << log.str();
	const auto masked = s.dir / "eli_qtaim_masked.cube";
	ASSERT_TRUE(std::filesystem::exists(masked));
	const cube c = read_cube(masked);
	EXPECT_NE(c.get_comment1().find("QTAIM-masked ELI"), std::string::npos);
	EXPECT_NE(c.get_comment2().find("Selected atoms: 0"), std::string::npos);
	const double r = constants::ang2bohr(1.1);
	const double step_x = (2.0 + 2.0 * r) / 8.0, step_yz = 2.0 * r / 6.0;
	ASSERT_EQ(c.get_size(0), 4);
	ASSERT_EQ(c.get_size(1), 5);
	ASSERT_EQ(c.get_size(2), 5);
	//1e-6 and not 0 because the cube header carries the origin and the step in six decimals (the residual is 5.3e-08)
	EXPECT_NEAR(c.get_origin(0) + 3.0 * step_x, 0.0, 1e-6) << "the last x column is the mirror plane itself";
	EXPECT_NEAR(c.get_origin(0), -1.0 - r + step_x, 1e-5);
	EXPECT_NEAR(c.get_origin(1), -r + step_yz, 1e-5);
	EXPECT_NEAR(c.get_origin(2), -r + step_yz, 1e-5);
	EXPECT_NEAR(c.get_vector(0, 0), step_x, 1e-5);
	EXPECT_NEAR(c.get_vector(1, 1), step_yz, 1e-5);
	EXPECT_NEAR(c.get_vector(2, 2), step_yz, 1e-5);
	int kept = 0, background = 0;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 5; j++)
			for (int k = 0; k < 5; k++)
			{
				const double v = c.get_value(i, j, k);
				const bool corner = (i == 0 || i == 3) && (j == 0 || j == 4) && (k == 0 || k == 4);
				if (v == -1.0)
				{
					EXPECT_TRUE(corner) << "background at " << i << " " << j << " " << k;
					background++;
					continue;
				}
				EXPECT_FALSE(corner) << "corner " << j << " " << k << " lies outside the evaluation radius";
				EXPECT_TRUE(std::isfinite(v));
				EXPECT_GT(v, 0.0);
				kept++;
			}
	EXPECT_EQ(kept, 92);
	EXPECT_EQ(background, 8);
}

//debug run of ELI_analysis on the two-Gaussian H2 model: rho = 4 (exp(-2 r_a^2) + exp(-2 r_b^2)) has the two nuclear
//maxima (at x = +-0.99933, the tail of the other Gaussian pulls them in by 2 exp(-8)) and one saddle at the origin
//with Hessian eigenvalues (-4, -4, 12) rho, Laplacian 4 rho and ellipticity 0; the debug listing adds the eigenvectors
//and the x axis carries the positive curvature. The QTAIM basins split at the midplane, so both hydrogens hold the
//same electron count and charge = Z - electrons; -debug also drops rho.cube and eli.cube into the working directory
TEST(BondwiseCoverageEliTests, DebugListsTheCriticalPointsOfTwoGaussians)
{
	Scratch s("DebugListsTheCriticalPointsOfTwoGaussians");
	std::filesystem::current_path(s.dir);
	WFN wavy = h2_wfn();
	wavy.set_path(s.dir / "h2.wfn");
	options opt;
	opt.debug = true;
	opt.properties.radius = 1.6;
	opt.properties.resolution = 0.25;
	std::string out;
	{
		CoutCapture cap;
		ELI_analysis(wavy, opt);
		out = cap.str();
	}
	EXPECT_TRUE(std::filesystem::exists(s.dir / "rho.cube"));
	EXPECT_TRUE(std::filesystem::exists(s.dir / "eli.cube"));
	EXPECT_NE(out.find("Calcualting grid of size 18 x 14 x 14"), std::string::npos) << out;
	EXPECT_NE(out.find("Density Critical Points (3 found):"), std::string::npos);
	const std::vector<CriticalPoint> cps = parse_critical_points(out);
	ASSERT_EQ(cps.size(), 3u) << out;
	int attractors = 0, bonds = 0;
	for (const CriticalPoint& cp : cps)
	{
		EXPECT_TRUE(cp.has_vecs) << "debug prints the eigenvectors";
		EXPECT_NEAR(cp.pos[1], 0.0, 2e-3);
		EXPECT_NEAR(cp.pos[2], 0.0, 2e-3);
		ASSERT_GT(cp.rho, 0.0);
		if (cp.type == "attractor")
		{
			attractors++;
			EXPECT_NEAR(std::abs(cp.pos[0]), 0.99933, 2e-3);
			EXPECT_EQ(cp.owner, cp.pos[0] < 0.0 ? "H0" : "H1");
			continue;
		}
		ASSERT_EQ(cp.type, "bond") << out;
		bonds++;
		EXPECT_TRUE(cp.converged);
		EXPECT_NE(cp.owner.find("bond"), std::string::npos);
		EXPECT_NE(cp.owner.find("H0"), std::string::npos);
		EXPECT_NE(cp.owner.find("H1"), std::string::npos);
		EXPECT_NEAR(cp.pos[0], 0.0, 1e-3);
		EXPECT_NEAR(cp.eig[0] / cp.rho, -4.0, 5e-3);
		EXPECT_NEAR(cp.eig[1] / cp.rho, -4.0, 5e-3);
		EXPECT_NEAR(cp.eig[2] / cp.rho, 12.0, 5e-3);
		EXPECT_NEAR(cp.lap / cp.rho, 4.0, 5e-3);
		EXPECT_NEAR(cp.ellipticity, 0.0, 1e-3);
		EXPECT_NEAR(std::abs(cp.vecs[2][0]), 1.0, 1e-3) << "the positive curvature runs along the bond";
		EXPECT_NEAR(cp.vecs[0][0] * cp.vecs[0][0] + cp.vecs[0][1] * cp.vecs[0][1] + cp.vecs[0][2] * cp.vecs[0][2], 1.0, 1e-3);
	}
	EXPECT_EQ(attractors, 2);
	EXPECT_EQ(bonds, 1);
	const auto rows = parse_qtaim_table(out);
	//H2 has two basins, so the table has two rows. A parser that dropped a row whose label carries a space
	//("NNA near H0", "H0-H1 bond") would still satisfy the two lookups below while hiding a third basin.
	ASSERT_EQ(rows.size(), 2u) << out;
	ASSERT_TRUE(rows.count("H0") && rows.count("H1")) << out;
	const double e0 = rows.at("H0").first, e1 = rows.at("H1").first;
	EXPECT_GT(e0, 0.0);
	EXPECT_NEAR(e0, e1, 5e-3 * (e0 + e1));
	EXPECT_NEAR(rows.at("H0").second, 1.0 - e0, 2e-4);
	EXPECT_NEAR(rows.at("H1").second, 1.0 - e1, 2e-4);
}

//The H2 table above has two single-word labels, so the row-count assertion in it cannot by itself show
//that the parser stopped dropping rows. This does, on a table written out here: two of its three rows
//carry a label with a space, which is what the real tables print for NNA, lone-pair and bond basins.
TEST(BondwiseCoverageEliTests, QtaimTableParserKeepsLabelsThatCarrySpaces)
{
	const std::string out =
		"QTAIM Analysis (atomic quadrature grids):\n"
		"  basin  label               electrons     charge      volume         maximum        x          y          z\n"
		"      1  Hg0                   79.4496     0.5504  12885.2903     367149.1320      0.000      0.000      0.000\n"
		"      2  NNA near H0            0.1234    -0.1234     42.0000          1.2300      0.500      0.000      0.000\n"
		"      3  H2-Hg0 bond            0.9000     0.1000     99.0000          0.4500     -0.500      0.000      0.000\n"
		"  total in basins:     80.4730   outside every basin:     0.0100\n";
	const auto rows = parse_qtaim_table(out);
	ASSERT_EQ(rows.size(), 3u);
	EXPECT_NEAR(rows.at("Hg0").first, 79.4496, 1e-9);
	EXPECT_NEAR(rows.at("NNA near H0").first, 0.1234, 1e-9);
	EXPECT_NEAR(rows.at("NNA near H0").second, -0.1234, 1e-9);
	EXPECT_NEAR(rows.at("H2-Hg0 bond").first, 0.9000, 1e-9);
	EXPECT_NEAR(rows.at("H2-Hg0 bond").second, 0.1000, 1e-9);
	//and the same table read the way it used to be read, kept here so this test cannot quietly become a
	//tautology: a single >> for the label keeps one row of the three and says nothing about the two it lost
	size_t naive = 0;
	std::istringstream in(out);
	for (std::string line; std::getline(in, line); )
	{
		std::istringstream row(line);
		int idx = 0; std::string label; double electrons = 0.0, charge = 0.0;
		if (row >> idx >> label >> electrons >> charge)
			naive++;
	}
	EXPECT_EQ(naive, 1u) << "the old parser kept only the row whose label is a single word";
}

namespace {

	//every "total in basins: X   outside every basin: Y" line, in the order printed
	std::vector<std::pair<double, double>> parse_basin_totals(const std::string& out)
	{
		std::vector<std::pair<double, double>> totals;
		const std::string key = "  total in basins:", mid = "outside every basin:";
		for (size_t at = out.find(key); at != std::string::npos; at = out.find(key, at + 1))
		{
			const size_t m = out.find(mid, at);
			if (m == std::string::npos)
				break;
			totals.emplace_back(std::strtod(out.c_str() + at + key.size(), nullptr),
								std::strtod(out.c_str() + m + mid.size(), nullptr));
		}
		return totals;
	}

}  // namespace

//The two basin tables of an ECP wavefunction close on two different electron counts, and both are
//right: HgH2 with the def2 ECP on Hg has 22 electrons in its orbitals, and the QTAIM arm fills the
//60-electron core with a Thakkar density before integrating, so its table closes on 82 while the
//ELI-D arm - which only ever sees the orbitals - closes on 22. Comparing the wrong pair is how a
//correct ECP run gets reported as a 60-electron error, and it is the question a colleague asks first.
//
//Conservation is what is asserted rather than any particular population, because the cube resolution
//only decides where the boundaries fall: every quadrature point is charged either to a basin or to
//"outside every basin", so the two must add to the electron count however coarse the grid is. That
//is what makes this checkable at a resolution the suite can afford.
TEST(BondwiseCoverageEliTests, TheEcpBasinTablesCloseOnTheirOwnElectronCounts)
{
	const std::filesystem::path p = nos_test_repo_root() / "tests" / "ELI_heavy" / "hgh2_ecp.gbw";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/ELI_heavy/hgh2_ecp.gbw not found";
	Scratch s("TheEcpBasinTablesCloseOnTheirOwnElectronCounts");
	std::filesystem::current_path(s.dir);

	WFN wavy(p);
	ASSERT_TRUE(wavy.get_has_ECPs());
	const double explicit_electrons = wavy.count_nr_electrons();
	const double cores = static_cast<double>(wavy.get_nr_ECP_electrons());
	ASSERT_NEAR(explicit_electrons, 22.0, 1e-9);
	ASSERT_NEAR(cores, 60.0, 1e-9);

	options opt;
	opt.properties.radius = 2.0;
	opt.properties.resolution = 0.4;
	std::string out;
	{
		CoutCapture cap;
		ELI_analysis(wavy, opt);
		out = cap.str();
	}

	const auto totals = parse_basin_totals(out);
	ASSERT_EQ(totals.size(), 2u) << "expected a QTAIM table and an ELI-D table\n" << out;
	const double qtaim = totals[0].first + totals[0].second;
	const double eli = totals[1].first + totals[1].second;
	EXPECT_NEAR(qtaim, explicit_electrons + cores, 1e-2 * (explicit_electrons + cores))
		<< "the QTAIM arm integrates the Thakkar cores too, so it must close on "
		<< explicit_electrons + cores << " electrons\n" << out;
	EXPECT_NEAR(eli, explicit_electrons, 1e-2 * explicit_electrons)
		<< "the ELI-D arm sees only the orbitals, so it must close on " << explicit_electrons
		<< " electrons\n" << out;
	EXPECT_NEAR(qtaim - eli, cores, 2e-2 * cores)
		<< "the difference between the two arms is the ECP core the one fills and the other cannot see";
}

// A free atom is the system ELI-D is calibrated against, and until af99fa2a the analysis could not
// run on one: the ELI label branch looked for the two atoms nearest each maximum and aborted the
// process when there was only one. Six of the twenty-two reader inputs in the robustness matrix
// failed there and nowhere else, five of them the only ELI coverage their reader has. This is the
// cheapest of the six, and it checks the two things a lone atom must satisfy: the basins plus
// whatever left them close on the electron count the file itself declares, and no basin is a bond,
// because there is nothing to bond to. If the abort comes back this test kills the whole binary
// rather than failing, which is the nature of err_checkf and worth knowing when it happens.
TEST(BondwiseCoverageEliTests, AFreeAtomIsAnalysedAndItsBasinsCloseOnItsOwnElectrons)
{
	const std::filesystem::path p = nos_test_repo_root() / "tests" / "molden_file" / "f_ref.wfx";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/molden_file/f_ref.wfx not found";
	Scratch s("AFreeAtomIsAnalysedAndItsBasinsCloseOnItsOwnElectrons");
	std::filesystem::current_path(s.dir);

	WFN wavy(p);
	ASSERT_EQ(wavy.get_ncen(), 1) << "the fixture is meant to be a single fluoride ion";
	const double electrons = wavy.count_nr_electrons();
	ASSERT_NEAR(electrons, 10.0, 1e-9) << "f_ref.wfx declares 10 electrons (F-)";

	options opt;
	opt.properties.radius = 2.5;
	opt.properties.resolution = 0.3;
	std::string out;
	{
		CoutCapture cap;
		ELI_analysis(wavy, opt);
		out = cap.str();
	}

	const auto totals = parse_basin_totals(out);
	ASSERT_EQ(totals.size(), 2u) << "expected a QTAIM table and an ELI-D table\n" << out;
	EXPECT_NEAR(totals[0].first + totals[0].second, electrons, 1e-2 * electrons)
		<< "QTAIM on one atom: everything is in that atom's basin\n" << out;
	EXPECT_NEAR(totals[1].first + totals[1].second, electrons, 1e-2 * electrons)
		<< "ELI-D on one atom: core shells plus valence, and nothing may go missing\n" << out;
	EXPECT_EQ(out.find(" bond"), std::string::npos)
		<< "a lone atom has no bond basin, so no label may carry one\n" << out;
	//the label is the atom's own label with its index appended, so this fixture prints "F10":
	//its single atom is named F1. Read the name from the wavefunction rather than spelling it out.
	const std::string own = wavy.get_atom_label(0) + "0";
	EXPECT_NE(out.find(own), std::string::npos) << "every basin belongs to atom " << own << "\n" << out;
}

namespace
{
	//"  basin  label  <first number>  ..." rows of the table that follows `after`, up to the
	//"total in basins" line that closes it. The first number is the population in every one of
	//the four tables, whatever the remaining columns are.
	std::map<std::string, double> parse_basin_rows(const std::string& out, const std::string& after)
	{
		std::map<std::string, double> rows;
		const size_t start = out.find(after);
		if (start == std::string::npos)
			return rows;
		const size_t stop = out.find("total in basins", start);
		std::istringstream in(out.substr(start, stop == std::string::npos ? std::string::npos : stop - start));
		std::string line;
		while (std::getline(in, line))
		{
			std::istringstream row(line);
			int index = 0;
			std::string label;
			double value = 0;
			if (row >> index >> label >> value)
				rows[label] = value;   //a two-word label ("H2-Li4 bond") keeps its first word, which is enough
		}
		return rows;
	}
}  // namespace

//An invariant that needs no reference calculation: nh3li.gbw is NH3 with a lithium on the
//three-fold axis, so H1, H2 and H3 are related by symmetry and must carry the same population.
//Anything else is the code's own error, and it is measurable without AIMAll or DGrid.
//
//The QTAIM arm holds the orbit to 7e-4 electrons and is asserted here. The ELI-D arm on the very
//same runs does not, and the resolution scan says why it is not asserted: the spread over the three
//equivalent hydrogens is 0.2077 e at resolution 0.4, 0.7420 e at 0.3 (one hydrogen at 2.6756
//against two at 1.93), 0.1546 e at 0.2 and 0.0661 e at 0.15, while the conserved total stays at
//12.67-12.68 throughout and the basin count grows 6, 6, 8, 9 with spurious 0.001-electron "bond"
//basins. A defect the conserved total cannot see, in the voxel watershed that decides the
//boundaries rather than in the field or the quadrature - the same quadrature integrates the QTAIM
//basins of this molecule to 7e-4. Left as measured, not asserted, because no resolution the suite
//can afford makes the ELI numbers symmetric.
TEST(BondwiseCoverageEliTests, QtaimHoldsTheThreefoldOrbitOfNH3Li)
{
	const std::filesystem::path p = nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
	Scratch s("QtaimHoldsTheThreefoldOrbitOfNH3Li");
	std::filesystem::current_path(s.dir);

	WFN wavy(p);
	ASSERT_EQ(wavy.get_ncen(), 5);
	options opt;
	opt.properties.radius = 2.5;
	opt.properties.resolution = 0.4;   //the QTAIM basins follow the analytic gradient, so this only sizes the ELI arm
	std::string out;
	{
		CoutCapture cap;
		ELI_analysis(wavy, opt);
		out = cap.str();
	}

	const auto rows = parse_basin_rows(out, "QTAIM Analysis");
	ASSERT_EQ(rows.count("H1") + rows.count("H2") + rows.count("H3"), 3u)
		<< "expected one QTAIM basin per hydrogen\n" << out;
	const double h1 = rows.at("H1"), h2 = rows.at("H2"), h3 = rows.at("H3");
	EXPECT_NEAR(h1, h2, 2e-3) << "H1 and H2 are related by the three-fold axis\n" << out;
	EXPECT_NEAR(h1, h3, 2e-3) << "H1 and H3 are related by the three-fold axis\n" << out;
	EXPECT_NEAR(h2, h3, 2e-3) << "H2 and H3 are related by the three-fold axis\n" << out;

	//and the whole table closes on the molecule's own electron count, 13 for NH3Li
	const auto totals = parse_basin_totals(out);
	ASSERT_GE(totals.size(), 1u) << out;
	const double electrons = wavy.count_nr_electrons();
	EXPECT_NEAR(totals[0].first + totals[0].second, electrons, 1e-2 * electrons) << out;
}
