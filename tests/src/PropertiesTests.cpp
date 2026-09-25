#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/atoms.h"
#include "core/cell.h"
#include "core/cube.h"
#include "core/wfn_class.h"
#include "core/spherical_density.h"
#include "core/GridManager.h"
#include "core/properties.h"
#include "core/b2c.h"

namespace {

// h2-like model: two protons at (-R, 0, 0) and (R, 0, 0), one unnormalised s primitive
// exp(-a r^2) on each, a doubly occupied bonding MO and an antibonding MO of occupation n1
// (empty by default), both normalised through the analytic gaussian overlaps. With n1 = 2 the
// Fermi-hole curvature rho tau - |grad rho|^2 / 4 = n0 n1 |psi0 grad psi1 - psi1 grad psi0|^2
// (Lagrange identity over the occupied pair) is a product with no cancellation, which gives
// ELF and ELI-D in closed form; with a single orbital it is 0 and ELI-D is rounding noise
struct H2Model {
	double a = 1.0;
	double R = 1.0;
	double n1 = 0.0;
	double cp = 0.0;
	double cm = 0.0;
	WFN wavy{ e_origin::NOT_YET_DEFINED };

	H2Model(double alpha, double half_distance, bool with_virtual = true, double antibonding_occ = 0.0) : a(alpha), R(half_distance), n1(with_virtual ? antibonding_occ : 0.0)
	{
		const double SAA = std::pow(constants::PI / (2.0 * a), 1.5);
		const double SAB = SAA * std::exp(-a * (2.0 * R) * (2.0 * R) / 2.0);
		cp = 1.0 / std::sqrt(2.0 * (SAA + SAB));
		cm = 1.0 / std::sqrt(2.0 * (SAA - SAB));
		wavy.push_back_atom("H", -R, 0.0, 0.0, 1);
		wavy.push_back_atom("H", R, 0.0, 0.0, 1);
		wavy.push_back_MO(1, 2.0, -0.5);
		if (with_virtual)
			wavy.push_back_MO(2, antibonding_occ, 0.3);
		double on_A[2] = { cp, cm };
		double on_B[2] = { cp, -cm };
		wavy.add_primitive(1, 1, a, on_A);
		wavy.add_primitive(2, 1, a, on_B);
		constants::exp_cutoff = -23.5; // read_* calls set_exp_cutoff() on this global; keep each model order independent
	}
	d3 A() const { return { -R, 0.0, 0.0 }; }
	d3 B() const { return { R, 0.0, 0.0 }; }
	double gA(const d3 &p) const { return std::exp(-a * ((p[0] + R) * (p[0] + R) + p[1] * p[1] + p[2] * p[2])); }
	double gB(const d3 &p) const { return std::exp(-a * ((p[0] - R) * (p[0] - R) + p[1] * p[1] + p[2] * p[2])); }
	// grad g = -2a (p - c) g, lap g = (4 a^2 |p - c|^2 - 6a) g
	d3 grad_g(const d3 &p, const d3 &c, double g) const { return { -2.0 * a * (p[0] - c[0]) * g, -2.0 * a * (p[1] - c[1]) * g, -2.0 * a * (p[2] - c[2]) * g }; }
	double lap_g(const d3 &p, const d3 &c, double g) const { return (4.0 * a * a * ((p[0] - c[0]) * (p[0] - c[0]) + (p[1] - c[1]) * (p[1] - c[1]) + (p[2] - c[2]) * (p[2] - c[2])) - 6.0 * a) * g; }
	double psi0(const d3 &p) const { return cp * (gA(p) + gB(p)); }
	double psi1(const d3 &p) const { return cm * (gA(p) - gB(p)); }
	d3 grad_psi0(const d3 &p) const
	{
		const d3 dA = grad_g(p, A(), gA(p)), dB = grad_g(p, B(), gB(p));
		return { cp * (dA[0] + dB[0]), cp * (dA[1] + dB[1]), cp * (dA[2] + dB[2]) };
	}
	d3 grad_psi1(const d3 &p) const
	{
		const d3 dA = grad_g(p, A(), gA(p)), dB = grad_g(p, B(), gB(p));
		return { cm * (dA[0] - dB[0]), cm * (dA[1] - dB[1]), cm * (dA[2] - dB[2]) };
	}
	double lap_psi0(const d3 &p) const { return cp * (lap_g(p, A(), gA(p)) + lap_g(p, B(), gB(p))); }
	double lap_psi1(const d3 &p) const { return cm * (lap_g(p, A(), gA(p)) - lap_g(p, B(), gB(p))); }
	double rho(const d3 &p) const { return 2.0 * psi0(p) * psi0(p) + n1 * psi1(p) * psi1(p); }
	// grad rho = sum_i 2 n_i psi_i grad psi_i
	d3 grad_rho(const d3 &p) const
	{
		const d3 g0 = grad_psi0(p), g1 = grad_psi1(p);
		const double f0 = 4.0 * psi0(p), f1 = 2.0 * n1 * psi1(p);
		return { f0 * g0[0] + f1 * g1[0], f0 * g0[1] + f1 * g1[1], f0 * g0[2] + f1 * g1[2] };
	}
	// lap rho = sum_i 2 n_i (|grad psi_i|^2 + psi_i lap psi_i)
	double lap_rho(const d3 &p) const
	{
		const d3 g0 = grad_psi0(p), g1 = grad_psi1(p);
		return 4.0 * (g0[0] * g0[0] + g0[1] * g0[1] + g0[2] * g0[2] + psi0(p) * lap_psi0(p)) + 2.0 * n1 * (g1[0] * g1[0] + g1[1] * g1[1] + g1[2] * g1[2] + psi1(p) * lap_psi1(p));
	}
	// rho tau - |grad rho|^2 / 4 with tau = sum_i n_i |grad psi_i|^2: 2 n1 |psi0 grad psi1 - psi1 grad psi0|^2
	double curvature(const d3 &p) const
	{
		const d3 g0 = grad_psi0(p), g1 = grad_psi1(p);
		const double s0 = psi0(p), s1 = psi1(p);
		double w = 0.0;
		for (int k = 0; k < 3; k++)
			w += (s0 * g1[k] - s1 * g0[k]) * (s0 * g1[k] - s1 * g0[k]);
		return 2.0 * n1 * w;
	}
	// ELI-D for a closed shell: rho_s (12 / g_s)^(3/8) with rho_s = rho / 2 and g_s = curvature / 4
	double eli(const d3 &p) const { return 0.5 * rho(p) * std::pow(48.0 / curvature(p), 3.0 / 8.0); }
	// ELF: 1 / (1 + chi^2), chi = D / D_h with the Pauli kinetic energy density D = curvature / (2 rho) and
	// D_h = rho^(5/3) / constants::ctelf. The code's ctelf = 10 / (3 (6 pi^2)^(2/3)) = 0.2194 is 2^(2/3) below the
	// Becke-Edgecombe closed-shell 10 / (3 (3 pi^2)^(2/3)) = 0.3483; the test takes the code's constant so the
	// functional form is what it checks, the constant itself is Florian's call
	double elf(const d3 &p) const
	{
		const double chi = constants::ctelf * std::pow(rho(p), -5.0 / 3.0) * curvature(p) / (2.0 * rho(p));
		return 1.0 / (1.0 + chi * chi);
	}
};

// n^3 voxels of spacing h centred on the origin, axis aligned
cube make_grid(int n, double h)
{
	cube grid({ n, n, n }, 0, true);
	for (int k = 0; k < 3; k++)
	{
		grid.set_origin(k, -0.5 * (n - 1) * h);
		grid.set_vector(k, k, h);
	}
	grid.calc_dv(); // set_vector does not refresh dv
	return grid;
}

// the 16 property cubes Calc_Prop indexes by cube_type, on the geometry of grid; Rho and the
// listed types carry storage, every other one is an unloaded placeholder
std::vector<cube> property_cubes(const cube &grid, std::initializer_list<int> loaded)
{
	std::vector<cube> cubes;
	for (int t = 0; t < 16; t++)
	{
		bool on = t == cube_type::Rho;
		for (const int l : loaded)
			if (l == t) on = true;
		cubes.emplace_back(grid.get_sizes(), 2, on);
		for (int k = 0; k < 3; k++)
		{
			cubes.back().set_origin(k, grid.get_origin(k));
			cubes.back().set_vector(k, k, grid.get_vector(k, k));
		}
		cubes.back().calc_dv();
	}
	return cubes;
}

double dist(const d3 &p, const d3 &q)
{
	return std::sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]));
}

bool inside(const WFN &w, const d3 &p, double radius_ang)
{
	for (int a = 0; a < w.get_ncen(); a++)
		if (dist(p, w.get_atom_pos(a)) < constants::ang2bohr(radius_ang))
			return true;
	return false;
}

std::filesystem::path temp_dir(const std::string &name)
{
	const std::filesystem::path dir = std::filesystem::temp_directory_path() / ("nosphera2_properties_" + name);
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir);
	return dir;
}

constexpr int N = 17;
constexpr double H = 0.25;
constexpr double RADIUS = 1.0; // angstrom = 1.89 bohr: the corners of the 2 bohr box lie outside

} // namespace

// Calc_Rho evaluates the closed-form two-gaussian density and zeroes every voxel farther than
// the radius from both nuclei; the grid corner and the nucleus pin both branches
TEST(PropertiesGridTests, CalcRhoMatchesAnalyticDensityAndMasksOutsideRadius)
{
	const H2Model m(1.0, 1.0);
	cube rho = make_grid(N, H);
	std::ostringstream log;
	Calc_Rho(rho, m.wavy, RADIUS, log, false);
	int masked = 0;
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = rho.get_pos(x, y, z);
				if (inside(m.wavy, p, RADIUS))
					EXPECT_NEAR(rho.get_value(x, y, z), m.rho(p), 1e-12);
				else
				{
					EXPECT_EQ(rho.get_value(x, y, z), 0.0);
					masked++;
				}
			}
	EXPECT_GT(masked, 0);
	EXPECT_EQ(rho.get_value(0, 0, 0), 0.0);
	const double e4 = std::exp(-4.0);
	EXPECT_NEAR(rho.get_value(4, 8, 8), 2.0 * m.cp * m.cp * (1.0 + e4) * (1.0 + e4), 1e-12);
	EXPECT_TRUE(rho.get_loaded());
}

// Calc_Rho_spherical_harmonics has no radius mask: every voxel carries the analytic density
TEST(PropertiesGridTests, CalcRhoSphericalHarmonicsFillsEveryVoxel)
{
	const H2Model m(1.0, 1.0);
	cube rho = make_grid(N, H);
	std::ostringstream log;
	Calc_Rho_spherical_harmonics(rho, m.wavy, log);
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
				EXPECT_NEAR(rho.get_value(x, y, z), m.rho(rho.get_pos(x, y, z)), 1e-12);
	EXPECT_GT(rho.get_value(0, 0, 0), 0.0);
}

// Calc_Cube with a constant function counts the voxels inside the mask; with wrap every voxel
// receives its 27 periodic images, so a huge radius makes each entry exactly 27, and x^2 + y
// summed over the images pins the image period L = n h per axis: 27 x^2 + 18 L^2 + 27 y
TEST(PropertiesGridTests, CalcCubeMasksByRadiusAndWrapSumsTwentySevenImages)
{
	const H2Model m(1.0, 1.0);
	cube one = make_grid(N, H);
	std::ostringstream log;
	const auto unity = [](const d3 &) { return 1.0; };
	Calc_Cube(one, m.wavy, unity, RADIUS, log, false);
	int expected = 0;
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
				expected += inside(m.wavy, one.get_pos(x, y, z), RADIUS) ? 1 : 0;
	EXPECT_NEAR(one.sum(), expected * one.get_dv(), 1e-12); // sum() already carries dv
	EXPECT_GT(expected, 100);
	EXPECT_LT(expected, N * N * N);

	cube wrapped = make_grid(5, H);
	Calc_Cube(wrapped, m.wavy, unity, 100.0, log, true);
	for (int x = 0; x < 5; x++)
		for (int y = 0; y < 5; y++)
			for (int z = 0; z < 5; z++)
				EXPECT_NEAR(wrapped.get_value(x, y, z), 27.0, 1e-12);

	const auto quad = [](const d3 &p) { return p[0] * p[0] + p[1]; };
	Calc_Cube(wrapped, m.wavy, quad, 100.0, log, true);
	const double L = 5 * H;
	for (int x = 0; x < 5; x++)
		for (int y = 0; y < 5; y++)
			for (int z = 0; z < 5; z++)
			{
				const d3 p = wrapped.get_pos(x, y, z);
				EXPECT_NEAR(wrapped.get_value(x, y, z), 27.0 * p[0] * p[0] + 18.0 * L * L + 27.0 * p[1], 1e-10);
			}
}

// Calc_Eli on the two-orbital model matches the closed-form ELI-D inside the radius and is zero
// outside. The evaluator forms rho tau - |grad rho|^2 / 4 by cancellation, which at the box edge
// (gB^2 / gA^2 ~ e^-16) loses ~7 digits, hence the 1e-6 relative tolerance; a wrong constant,
// exponent or sign moves the value by O(1)
TEST(PropertiesGridTests, CalcEliMatchesAnalyticEliDAndMask)
{
	const H2Model m(1.0, 1.0, true, 2.0);
	cube eli = make_grid(N, H);
	std::ostringstream log;
	Calc_Eli(eli, m.wavy, RADIUS, log, false);
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = eli.get_pos(x, y, z);
				const double expected = inside(m.wavy, p, RADIUS) ? m.eli(p) : 0.0;
				EXPECT_NEAR(eli.get_value(x, y, z), expected, 1e-6 * expected);
			}
	// bond midpoint: psi1 = 0, grad psi0 = 0, so the curvature is 4 psi0^2 |grad psi1|^2 with
	// grad psi1 = cm (-2)(-2R, 0, 0) e^-1 and psi0 = 2 cp e^-1 for a = R = 1
	const double s0 = 2.0 * m.cp * std::exp(-1.0), g1 = 4.0 * m.cm * std::exp(-1.0);
	EXPECT_NEAR(eli.get_value(8, 8, 8), 0.5 * 2.0 * s0 * s0 * std::pow(48.0 / (4.0 * s0 * s0 * g1 * g1), 3.0 / 8.0), 1e-10);
	EXPECT_EQ(eli.get_value(0, 0, 0), 0.0);
}

// Calc_RhoEli fills both cubes in one pass from computeRhoELI and must give the closed-form
// density and ELI-D; the cubes are pre-filled with junk to prove both are reset first
TEST(PropertiesGridTests, CalcRhoEliWithoutFieldMatchesAnalyticValues)
{
	const H2Model m(1.0, 1.0, true, 2.0);
	cube rho = make_grid(N, H), eli = make_grid(N, H);
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				rho.set_value(x, y, z, 7.0);
				eli.set_value(x, y, z, 7.0);
			}
	Calc_RhoEli(rho, eli, m.wavy, RADIUS);
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = rho.get_pos(x, y, z);
				const bool in = inside(m.wavy, p, RADIUS);
				EXPECT_NEAR(rho.get_value(x, y, z), in ? m.rho(p) : 0.0, 1e-12);
				const double expected = in ? m.eli(p) : 0.0;
				EXPECT_NEAR(eli.get_value(x, y, z), expected, 1e-6 * expected);
			}
}

// with a density_field rho comes from the field while ELI-D still comes from the orbitals;
// a wavefunction without orbitals leaves ELI-D at zero instead of dividing by nothing
TEST(PropertiesGridTests, CalcRhoEliWithDensityFieldTakesRhoFromFieldAndEliFromOrbitals)
{
	const H2Model m(1.0, 1.0, true, 2.0);
	density_field field;
	field.rho = [](const d3 &p) { return 1.0 + p[0] + 2.0 * p[1]; };
	field.grad = [](const d3 &, d3 &g) { g = { 1.0, 2.0, 0.0 }; };

	cube rho = make_grid(N, H), eli = make_grid(N, H);
	Calc_RhoEli(rho, eli, m.wavy, RADIUS, &field);
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = rho.get_pos(x, y, z);
				if (!inside(m.wavy, p, RADIUS))
				{
					EXPECT_EQ(rho.get_value(x, y, z), 0.0);
					EXPECT_EQ(eli.get_value(x, y, z), 0.0);
					continue;
				}
				EXPECT_NEAR(rho.get_value(x, y, z), 1.0 + p[0] + 2.0 * p[1], 1e-12);
				EXPECT_NEAR(eli.get_value(x, y, z), m.eli(p), 1e-6 * m.eli(p));
			}

	WFN bare(e_origin::NOT_YET_DEFINED);
	bare.push_back_atom("H", -1.0, 0.0, 0.0, 1);
	bare.push_back_atom("H", 1.0, 0.0, 0.0, 1);
	ASSERT_EQ(bare.get_nmo(), 0);
	cube rho2 = make_grid(N, H), eli2 = make_grid(N, H);
	Calc_RhoEli(rho2, eli2, bare, RADIUS, &field);
	EXPECT_NEAR(rho2.get_value(8, 8, 8), 1.0, 1e-12);
	EXPECT_EQ(eli2.max_value(), 0.0);
	EXPECT_NEAR(eli2.sum(), 0.0, 0.0);
}

// every loaded property cube of a Calc_Prop run on the two-orbital model against its closed form:
// Elf and Eli relative to the value (the evaluators form the Fermi-hole curvature by
// cancellation, see CalcEliMatchesAnalyticEliDAndMask), Lap absolute, all zero outside the radius
void expect_prop_cubes_analytic(const std::vector<cube> &cubes, const H2Model &m)
{
	const cube &grid = cubes[cube_type::Rho];
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = grid.get_pos(x, y, z);
				const bool in = inside(m.wavy, p, RADIUS);
				if (cubes[cube_type::Elf].get_loaded())
					EXPECT_NEAR(cubes[cube_type::Elf].get_value(x, y, z), in ? m.elf(p) : 0.0, 1e-8);
				if (cubes[cube_type::Eli].get_loaded())
					EXPECT_NEAR(cubes[cube_type::Eli].get_value(x, y, z), in ? m.eli(p) : 0.0, in ? 1e-6 * m.eli(p) : 0.0);
				if (cubes[cube_type::Lap].get_loaded())
					EXPECT_NEAR(cubes[cube_type::Lap].get_value(x, y, z), in ? m.lap_rho(p) : 0.0, 1e-10);
			}
}

// Calc_Prop with only Elf requested takes the computeELF branch and leaves Rho untouched; the
// two-orbital ELF is neither 0 nor 1 across the box, so the Pauli kinetic term is really exercised
TEST(PropertiesGridTests, CalcPropElfOnlyMatchesAnalyticElf)
{
	const H2Model m(1.0, 1.0, true, 2.0);
	const cube grid = make_grid(N, H);
	std::vector<cube> cubes = property_cubes(grid, { cube_type::Elf });
	cubes[cube_type::Rho].set_value(8, 8, 8, 3.0);
	std::ostringstream log;
	Calc_Prop(cubes, m.wavy, RADIUS, log, true, false);
	expect_prop_cubes_analytic(cubes, m);
	EXPECT_LT(cubes[cube_type::Elf].get_value(8, 8, 8), 0.9);
	EXPECT_GT(cubes[cube_type::Elf].get_value(8, 8, 8), 0.1);
	EXPECT_EQ(cubes[cube_type::Rho].get_value(8, 8, 8), 3.0);
	EXPECT_FALSE(cubes[cube_type::RDG].get_loaded());
}

// Calc_Prop with only Eli requested takes the computeELI branch
TEST(PropertiesGridTests, CalcPropEliOnlyMatchesAnalyticEliD)
{
	const H2Model m(1.0, 1.0, true, 2.0);
	const cube grid = make_grid(N, H);
	std::vector<cube> cubes = property_cubes(grid, { cube_type::Eli });
	std::ostringstream log;
	Calc_Prop(cubes, m.wavy, RADIUS, log, true, false);
	expect_prop_cubes_analytic(cubes, m);
}

// Calc_Prop with only Lap requested takes the computeLap branch; the one-orbital model keeps the
// hand-derived midpoint value as an anchor for the model's own lap_rho
TEST(PropertiesGridTests, CalcPropLapOnlyMatchesAnalyticLaplacian)
{
	const H2Model m(1.0, 1.0);
	const cube grid = make_grid(N, H);
	std::vector<cube> cubes = property_cubes(grid, { cube_type::Lap });
	std::ostringstream log;
	Calc_Prop(cubes, m.wavy, RADIUS, log, true, false);
	expect_prop_cubes_analytic(cubes, m);
	// rho = 2 cp^2 (gA + gB)^2, grad(gA + gB) = 0 at the origin, lap(gA + gB) = -4 e^-1 for a = R = 1
	const double analytic = -32.0 * m.cp * m.cp * std::exp(-2.0);
	EXPECT_NEAR(cubes[cube_type::Lap].get_value(8, 8, 8), analytic, 1e-10);
	EXPECT_NEAR(m.lap_rho({ 0.0, 0.0, 0.0 }), analytic, 1e-14);
}

// Eli and Lap together take the computeLapELI branch
TEST(PropertiesGridTests, CalcPropEliAndLapMatchesAnalyticValues)
{
	const H2Model m(1.0, 1.0, true, 2.0);
	const cube grid = make_grid(N, H);
	std::vector<cube> cubes = property_cubes(grid, { cube_type::Eli, cube_type::Lap });
	std::ostringstream log;
	Calc_Prop(cubes, m.wavy, RADIUS, log, true, false);
	expect_prop_cubes_analytic(cubes, m);
}

// Elf, Eli and Lap together take the computeLapELIELF branch
TEST(PropertiesGridTests, CalcPropElfEliLapMatchesAnalyticValues)
{
	const H2Model m(1.0, 1.0, true, 2.0);
	const cube grid = make_grid(N, H);
	std::vector<cube> cubes = property_cubes(grid, { cube_type::Elf, cube_type::Eli, cube_type::Lap });
	std::ostringstream log;
	Calc_Prop(cubes, m.wavy, RADIUS, log, true, false);
	expect_prop_cubes_analytic(cubes, m);
}

// Elf and Eli without Lap take the computeELIELF branch
TEST(PropertiesGridTests, CalcPropElfAndEliMatchesAnalyticValues)
{
	const H2Model m(1.0, 1.0, true, 2.0);
	const cube grid = make_grid(N, H);
	std::vector<cube> cubes = property_cubes(grid, { cube_type::Elf, cube_type::Eli });
	std::ostringstream log;
	Calc_Prop(cubes, m.wavy, RADIUS, log, true, false);
	expect_prop_cubes_analytic(cubes, m);
}

// the RDG run replaces Rho by sign(lambda_2) rho: negative on the nucleus where every hessian
// eigenvalue is negative, positive on the perpendicular bisector far from the bond where two are
// positive; RDG is the reduced gradient inside the radius and the 101 mask value outside
TEST(PropertiesGridTests, CalcPropRdgSignsRhoByMiddleEigenvalueAndMasksWith101)
{
	const H2Model m(1.0, 1.0);
	const cube grid = make_grid(N, H);
	std::vector<cube> cubes = property_cubes(grid, { cube_type::RDG, cube_type::Lap, cube_type::Eli });
	std::ostringstream log;
	Calc_Prop(cubes, m.wavy, RADIUS, log, true, false);
	const cube &rho = cubes[cube_type::Rho];
	const cube &rdg = cubes[cube_type::RDG];
	EXPECT_NEAR(rho.get_value(4, 8, 8), -m.rho(grid.get_pos(4, 8, 8)), 1e-12);
	EXPECT_NEAR(rho.get_value(8, 8, 8), -m.rho(grid.get_pos(8, 8, 8)), 1e-12);
	EXPECT_NEAR(rho.get_value(8, 14, 8), m.rho(grid.get_pos(8, 14, 8)), 1e-12);
	EXPECT_EQ(rdg.get_value(0, 0, 0), 101.0);
	EXPECT_EQ(rho.get_value(0, 0, 0), 0.0);
	const double alpha = 1.0 / (2.0 * std::cbrt(3.0 * constants::PI * constants::PI));
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = grid.get_pos(x, y, z);
				if (!inside(m.wavy, p, RADIUS))
				{
					EXPECT_EQ(rdg.get_value(x, y, z), 101.0);
					continue;
				}
				const d3 g = m.grad_rho(p);
				const double expected = alpha * std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) / std::pow(m.rho(p), 4.0 / 3.0);
				EXPECT_NEAR(rdg.get_value(x, y, z), expected, 1e-8 * std::max(1.0, expected));
				EXPECT_NEAR(std::abs(rho.get_value(x, y, z)), m.rho(p), 1e-12);
			}
}

// Calc_ESP on a single gaussian atom: Z/r from the nucleus minus the potential of a spherical
// gaussian charge cloud, N erf(sqrt(2a) r)/r, with N its electron count
TEST(PropertiesGridTests, CalcEspMatchesGaussianChargeCloudPotential)
{
	const double a = 0.8;
	WFN w(e_origin::NOT_YET_DEFINED);
	w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
	w.push_back_MO(1, 2.0, -0.9);
	const double c = 1.0 / std::sqrt(std::pow(constants::PI / (2.0 * a), 1.5));
	double coef[1] = { c };
	w.add_primitive(1, 1, a, coef);
	const double Nel = 2.0 * c * c * std::pow(constants::PI / (2.0 * a), 1.5);
	ASSERT_NEAR(Nel, 2.0, 1e-12);

	// even n keeps every voxel off the nucleus
	cube esp = make_grid(12, 0.3);
	std::ostringstream log;
	Calc_ESP(esp, w, 2.0, true, log, false);
	EXPECT_TRUE(log.str().empty());
	for (int x = 0; x < 12; x++)
		for (int y = 0; y < 12; y++)
			for (int z = 0; z < 12; z++)
			{
				const d3 p = esp.get_pos(x, y, z);
				const double r = dist(p, { 0.0, 0.0, 0.0 });
				if (r >= constants::ang2bohr(2.0))
				{
					EXPECT_EQ(esp.get_value(x, y, z), 0.0);
					continue;
				}
				const double expected = 2.0 / r - Nel * std::erf(std::sqrt(2.0 * a) * r) / r;
				EXPECT_NEAR(esp.get_value(x, y, z), expected, 1e-9);
			}
	EXPECT_GT(esp.get_value(6, 6, 6), 0.0);
}

// Calc_MO of the antibonding orbital reproduces cm (gA - gB), including its node on the
// bisector plane and the sign flip across it; the mask keeps the far corners at zero
TEST(PropertiesGridTests, CalcMoMatchesAnalyticAntibondingOrbital)
{
	const H2Model m(1.0, 1.0);
	cube mo = make_grid(N, H);
	std::ostringstream log;
	Calc_MO(mo, 1, m.wavy, RADIUS, log, false);
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = mo.get_pos(x, y, z);
				const double expected = inside(m.wavy, p, RADIUS) ? m.psi1(p) : 0.0;
				EXPECT_NEAR(mo.get_value(x, y, z), expected, 1e-12);
			}
	EXPECT_NEAR(mo.get_value(8, 8, 8), 0.0, 1e-15);
	EXPECT_GT(mo.get_value(4, 8, 8), 0.0);
	EXPECT_NEAR(mo.get_value(12, 8, 8), -mo.get_value(4, 8, 8), 1e-12);
	EXPECT_EQ(mo.get_value(0, 0, 0), 0.0);
}

// the four Fukui cubes are |psi_lumo|^2, |psi_homo|^2, their mean and their difference; a
// cube set without the Fukui_plus storage makes Calc_Fukui return without touching anything
TEST(PropertiesFukuiTests, CalcFukuiCubesAreFrontierOrbitalDensities)
{
	const H2Model m(1.0, 1.0);
	const cube grid = make_grid(N, H);
	std::vector<cube> cubes = property_cubes(grid, { cube_type::Fukui_plus, cube_type::Fukui_minus, cube_type::Fukui_zero, cube_type::Dual_Descriptor });
	std::ostringstream log;
	Calc_Fukui(cubes, m.wavy, 0, 1, RADIUS, log, false);
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = grid.get_pos(x, y, z);
				const bool in = inside(m.wavy, p, RADIUS);
				const double fp = in ? m.psi1(p) * m.psi1(p) : 0.0;
				const double fm = in ? m.psi0(p) * m.psi0(p) : 0.0;
				EXPECT_NEAR(cubes[cube_type::Fukui_plus].get_value(x, y, z), fp, 1e-12);
				EXPECT_NEAR(cubes[cube_type::Fukui_minus].get_value(x, y, z), fm, 1e-12);
				EXPECT_NEAR(cubes[cube_type::Fukui_zero].get_value(x, y, z), 0.5 * (fp + fm), 1e-12);
				EXPECT_NEAR(cubes[cube_type::Dual_Descriptor].get_value(x, y, z), fp - fm, 1e-12);
			}
	EXPECT_NEAR(cubes[cube_type::Fukui_plus].get_value(8, 8, 8), 0.0, 1e-15);
	EXPECT_LT(cubes[cube_type::Dual_Descriptor].get_value(8, 8, 8), 0.0);

	std::vector<cube> none = property_cubes(grid, { cube_type::Eli });
	none[cube_type::Eli].set_value(8, 8, 8, 5.0);
	Calc_Fukui(none, m.wavy, 0, 1, RADIUS, log, false);
	EXPECT_EQ(none[cube_type::Eli].get_value(8, 8, 8), 5.0);
	EXPECT_FALSE(none[cube_type::Fukui_plus].get_loaded());
}

// condensed Fukui functions of normalised orbitals integrate to one under every partition and
// split evenly between two equivalent protons; the printout carries the table and its sum row
// and an invalid result prints nothing at all
TEST(PropertiesFukuiTests, CondensedFukuiSumsToOneAndSplitsEvenlyOnH2)
{
	const H2Model m(1.0, 1.0);
	std::ostringstream grid_log;
	const CondensedFukuiResults r = Calc_Condensed_Fukui(m.wavy, 0, 1, cell(), 1, grid_log);
	ASSERT_TRUE(r.valid);
	ASSERT_EQ(r.labels.size(), 2u);
	ASSERT_EQ(r.f_plus.size(), 5u);
	ASSERT_EQ(r.f_minus.size(), 5u);
	for (int s = 0; s < 5; s++)
	{
		ASSERT_EQ(r.f_plus[s].size(), 2u);
		EXPECT_NEAR(r.f_plus[s][0] + r.f_plus[s][1], 1.0, 2e-2) << "partition " << s;
		EXPECT_NEAR(r.f_minus[s][0] + r.f_minus[s][1], 1.0, 2e-2) << "partition " << s;
		EXPECT_NEAR(r.f_plus[s][0], r.f_plus[s][1], 1e-6) << "partition " << s;
		EXPECT_NEAR(r.f_minus[s][0], r.f_minus[s][1], 1e-6) << "partition " << s;
	}

	// a LUMO sitting on atom B only: f+ lands on B, f- stays even, so a homo/lumo swap shows
	WFN skew(e_origin::NOT_YET_DEFINED);
	skew.push_back_atom("H", -1.0, 0.0, 0.0, 1);
	skew.push_back_atom("H", 1.0, 0.0, 0.0, 1);
	skew.push_back_MO(1, 2.0, -0.5);
	skew.push_back_MO(2, 0.0, 0.3);
	const double cB = 1.0 / std::sqrt(std::pow(constants::PI / 2.0, 1.5));
	double on_A[2] = { m.cp, 0.0 };
	double on_B[2] = { m.cp, cB };
	skew.add_primitive(1, 1, 1.0, on_A);
	skew.add_primitive(2, 1, 1.0, on_B);
	const CondensedFukuiResults k = Calc_Condensed_Fukui(skew, 0, 1, cell(), 1, grid_log);
	ASSERT_TRUE(k.valid);
	for (int s = 0; s < 5; s++)
	{
		EXPECT_GT(k.f_plus[s][1], 0.7) << "partition " << s;
		EXPECT_LT(k.f_plus[s][0], 0.3) << "partition " << s;
		EXPECT_NEAR(k.f_minus[s][0], k.f_minus[s][1], 1e-6) << "partition " << s;
	}

	std::ostringstream table;
	print_condensed_fukui(r, table);
	const std::string t = table.str();
	EXPECT_NE(t.find("Condensed (atom-summed) Fukui functions"), std::string::npos);
	EXPECT_NE(t.find("Hirshfeld"), std::string::npos);
	EXPECT_NE(t.find("EMBIS"), std::string::npos);
	EXPECT_NE(t.find("dual descriptor (f+ - f-)"), std::string::npos);
	EXPECT_NE(t.find("sum"), std::string::npos);
	EXPECT_NE(t.find("H"), std::string::npos);

	std::ostringstream empty;
	print_condensed_fukui(CondensedFukuiResults(), empty);
	EXPECT_TRUE(empty.str().empty());
	const CondensedFukuiResults no_lumo = Calc_Condensed_Fukui(m.wavy, 0, -1, cell(), 1, grid_log);
	EXPECT_FALSE(no_lumo.valid);
}

// fukui_analysis on a written wfn reports the frontier pair, the gap and writes the
// _fukui.dat summary beside the wavefunction with one block per partition
TEST(PropertiesIoTests, FukuiAnalysisWritesSummaryBesideWavefunction)
{
	const std::filesystem::path dir = temp_dir("fukui_io");
	const std::filesystem::path wfn_path = dir / "h2_frontier.wfn";
	const H2Model m(1.0, 1.0);
	ASSERT_TRUE(m.wavy.write_wfn(wfn_path, false, false));

	const bool hidden = constants::hide_timings;
	constants::hide_timings = true;
	options opt;
	opt.wfn = wfn_path;
	opt.no_date = true;
	opt.accuracy = 1;
	std::ostringstream log;
	fukui_analysis(opt, log);
	constants::hide_timings = hidden;

	const std::string s = log.str();
	EXPECT_NE(s.find("Read 2 atoms and 2 molecular orbitals (1 occupied)"), std::string::npos);
	EXPECT_NE(s.find("HOMO = MO 0"), std::string::npos);
	EXPECT_NE(s.find("LUMO = MO 1"), std::string::npos);
	EXPECT_NE(s.find("HOMO-LUMO gap: 0.800000 Hartree"), std::string::npos);
	EXPECT_NE(s.find("Summary (Hirshfeld partition):"), std::string::npos);
	EXPECT_NE(s.find("Wrote h2_frontier_fukui.dat"), std::string::npos);
	EXPECT_EQ(s.find("Unrestricted wavefunction."), std::string::npos);

	const std::filesystem::path dat = dir / "h2_frontier_fukui.dat";
	ASSERT_TRUE(std::filesystem::exists(dat));
	std::ifstream in(dat);
	std::stringstream buf;
	buf << in.rdbuf();
	in.close();
	const std::string d = buf.str();
	EXPECT_NE(d.find("HOMO_index 0"), std::string::npos);
	EXPECT_NE(d.find("LUMO_index 1"), std::string::npos);
	EXPECT_NE(d.find("LUMO_energy 0.300000"), std::string::npos);
	EXPECT_NE(d.find("unrestricted 0"), std::string::npos);
	for (const char *name : { "hirshfeld", "becke", "tfvc", "mbis", "embis" })
		EXPECT_NE(d.find(std::string("partition ") + name), std::string::npos) << name;
	// per atom: label f+ f- df, and f+ of the two protons agree
	std::istringstream lines(d.substr(d.find("partition hirshfeld")));
	std::string line, label;
	double fp[2] = { -1.0, -1.0 }, fm[2] = { -1.0, -1.0 }, df[2] = { 0.0, 0.0 };
	std::getline(lines, line);
	for (int a = 0; a < 2; a++)
	{
		std::getline(lines, line);
		std::istringstream iss(line);
		iss >> label >> fp[a] >> fm[a] >> df[a];
		EXPECT_EQ(label, "H");
	}
	EXPECT_NEAR(fp[0] + fp[1], 1.0, 2e-2);
	EXPECT_NEAR(fm[0] + fm[1], 1.0, 2e-2);
	EXPECT_NEAR(df[0], fp[0] - fm[0], 2e-6);
	std::filesystem::remove_all(dir);
}

// a wavefunction holding only occupied orbitals has no LUMO: fukui_analysis says so and
// writes nothing
TEST(PropertiesIoTests, FukuiAnalysisWithoutVirtualsExplainsAndWritesNothing)
{
	const std::filesystem::path dir = temp_dir("fukui_novirt");
	const std::filesystem::path wfn_path = dir / "h2_occ.wfn";
	const H2Model m(1.0, 1.0);
	ASSERT_TRUE(m.wavy.write_wfn(wfn_path, false, true));

	options opt;
	opt.wfn = wfn_path;
	opt.no_date = true;
	opt.accuracy = 1;
	std::ostringstream log;
	fukui_analysis(opt, log);
	const std::string s = log.str();
	EXPECT_NE(s.find("Read 2 atoms and 1 molecular orbitals (1 occupied)"), std::string::npos);
	EXPECT_NE(s.find("ERROR: no HOMO/LUMO pair found (HOMO index 0, LUMO index -1)"), std::string::npos);
	EXPECT_NE(s.find("stores no virtual orbitals"), std::string::npos);
	EXPECT_EQ(s.find("Wrote"), std::string::npos);
	EXPECT_FALSE(std::filesystem::exists(dir / "h2_occ_fukui.dat"));
	std::filesystem::remove_all(dir);
}

// do_combine_mo writes the sum and difference cube of one orbital from each fragment into the
// working directory plus a vmd script; half the sum of the two files read back is fragment 1's
// orbital and half their difference fragment 2's, both on the union box of the two fragments
// (the earlier "values mismatch" was Calc_MO's default wrap = true summing the molecular MO over the images of the box)
TEST(PropertiesIoTests, DoCombineMoWritesSumAndDifferenceCubes)
{
	const std::filesystem::path dir = temp_dir("combine_mo");
	const double a1 = 0.7, a2 = 1.1;
	WFN w1(e_origin::NOT_YET_DEFINED), w2(e_origin::NOT_YET_DEFINED);
	w1.push_back_atom("H", 0.0, 0.0, 0.0, 1);
	w1.push_back_MO(1, 2.0, -0.5);
	double c1[1] = { 0.9 };
	w1.add_primitive(1, 1, a1, c1);
	w2.push_back_atom("H", 0.5, 0.0, 0.0, 1);
	w2.push_back_MO(1, 2.0, -0.4);
	double c2[1] = { 0.6 };
	w2.add_primitive(1, 1, a2, c2);
	ASSERT_TRUE(w1.write_wfn(dir / "fragA.wfn", false, false));
	ASSERT_TRUE(w2.write_wfn(dir / "fragB.wfn", false, false));

	const std::filesystem::path old_cwd = std::filesystem::current_path();
	std::filesystem::current_path(dir);
	options opt;
	opt.combine_mo = { dir / "fragA.wfn", dir / "fragB.wfn" };
	opt.cmo1 = { 1 };
	opt.cmo2 = { 1 };
	opt.properties.radius = 0.5;
	const bool hidden = constants::hide_timings;
	constants::hide_timings = true;
	do_combine_mo(opt);
	constants::hide_timings = hidden;
	std::filesystem::current_path(old_cwd);

	const std::filesystem::path plus = dir / "fragA_1+fragB_1.cube";
	const std::filesystem::path minus = dir / "fragA_1-fragB_1.cube";
	ASSERT_TRUE(std::filesystem::exists(plus));
	ASSERT_TRUE(std::filesystem::exists(minus));
	ASSERT_TRUE(std::filesystem::exists(dir / "read_files.vmd"));
	std::ifstream vmd(dir / "read_files.vmd");
	std::stringstream vbuf;
	vbuf << vmd.rdbuf();
	vmd.close();
	EXPECT_NE(vbuf.str().find("mol new {fragA_1+fragB_1.cube}"), std::string::npos);
	EXPECT_NE(vbuf.str().find("mol addfile {fragA_1-fragB_1.cube}"), std::string::npos);

	WFN dummy_p(e_origin::NOT_YET_DEFINED), dummy_m(e_origin::NOT_YET_DEFINED);
	std::ostringstream log;
	cube sum(plus, true, dummy_p, log);
	cube diff(minus, true, dummy_m, log);
	ASSERT_TRUE(sum.get_loaded());
	ASSERT_TRUE(diff.get_loaded());
	ASSERT_EQ(sum.get_sizes(), diff.get_sizes());
	// union box: -0.5 A around x in [0, 0.5] bohr and around 0 in y, z; 0.1 A steps rounded up
	const double pad = constants::ang2bohr(0.5);
	EXPECT_NEAR(sum.get_origin(0), -pad, 1e-5);
	EXPECT_NEAR(sum.get_origin(1), -pad, 1e-5);
	// ceil of a step count that is an exact integer in principle may land one voxel either way
	EXPECT_NEAR(sum.get_size(0), constants::bohr2ang(0.5 + 2.0 * pad) / 0.1, 1.0);
	EXPECT_NEAR(sum.get_size(1), constants::bohr2ang(2.0 * pad) / 0.1, 1.0);
	double max_abs = 0.0;
	for (int x = 0; x < sum.get_size(0); x++)
		for (int y = 0; y < sum.get_size(1); y++)
			for (int z = 0; z < sum.get_size(2); z++)
			{
				const d3 p = sum.get_pos(x, y, z);
				const double mo1 = 0.9 * std::exp(-a1 * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]));
				const double mo2 = 0.6 * std::exp(-a2 * ((p[0] - 0.5) * (p[0] - 0.5) + p[1] * p[1] + p[2] * p[2]));
				const double s = sum.get_value(x, y, z), d = diff.get_value(x, y, z);
				// the cube file carries six significant digits
				EXPECT_NEAR(0.5 * (s + d), mo1, 2e-5);
				EXPECT_NEAR(0.5 * (s - d), mo2, 2e-5);
				max_abs = std::max(max_abs, std::abs(s));
			}
	EXPECT_GT(max_abs, 1.0);
	std::filesystem::remove_all(dir);
}

// the two-cube Calc_Static_Def is rho minus the sum of the free-atom Thakkar densities inside
// the radius; the vector version fed a spherical cube of exactly those sums must agree to
// rounding, and Calc_Spherical_Dens' interpolated table must agree with the direct sums
TEST(PropertiesPromolTests, StaticDeformationDensityMatchesThakkarSumsBothWays)
{
	const H2Model m(1.0, 1.0);
	const cube grid = make_grid(N, H);
	std::ostringstream log;
	cube rho = make_grid(N, H);
	Calc_Rho(rho, m.wavy, RADIUS, log, false);

	cube def = make_grid(N, H);
	Calc_Static_Def(def, rho, m.wavy, RADIUS, log, false);

	const Thakkar hydrogen(1);
	std::vector<cube> cubes = property_cubes(grid, { cube_type::DEF, cube_type::spherical_density });
	cubes[cube_type::Rho] = rho;
	cube &spher = cubes[cube_type::spherical_density];
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = grid.get_pos(x, y, z);
				spher.set_value(x, y, z, hydrogen.get_radial_density(dist(p, m.wavy.get_atom_pos(0))) + hydrogen.get_radial_density(dist(p, m.wavy.get_atom_pos(1))));
			}
	Calc_Static_Def(cubes, m.wavy, RADIUS, log, false);

	int checked = 0;
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = grid.get_pos(x, y, z);
				if (!inside(m.wavy, p, RADIUS))
				{
					EXPECT_EQ(def.get_value(x, y, z), 0.0);
					EXPECT_EQ(cubes[cube_type::DEF].get_value(x, y, z), 0.0);
					continue;
				}
				const double expected = m.rho(p) - spher.get_value(x, y, z);
				EXPECT_NEAR(def.get_value(x, y, z), expected, 1e-12);
				EXPECT_NEAR(cubes[cube_type::DEF].get_value(x, y, z), expected, 1e-12);
				checked++;
			}
	EXPECT_GT(checked, 100);
	// the free H atom holds one electron more than the gaussian bond puts near a proton's shell,
	// so the deformation density is negative right on the nucleus
	EXPECT_LT(def.get_value(4, 8, 8), 0.0);
}

// Calc_Spherical_Dens on the table interpolator reproduces the exact Thakkar sums
TEST(PropertiesPromolTests, SphericalDensInterpolatesThakkarSums)
{
	const H2Model m(1.0, 1.0);
	const cube grid = make_grid(N, H);
	std::ostringstream log;
	const Thakkar hydrogen(1);
	cube interp = make_grid(N, H);
	Calc_Spherical_Dens(interp, m.wavy, RADIUS, log, false);
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = grid.get_pos(x, y, z);
				const double exact = inside(m.wavy, p, RADIUS)
					? hydrogen.get_radial_density(dist(p, m.wavy.get_atom_pos(0))) + hydrogen.get_radial_density(dist(p, m.wavy.get_atom_pos(1)))
					: 0.0;
				EXPECT_NEAR(interp.get_value(x, y, z), exact, 2e-3 * std::max(1e-3, exact));
			}
}

// the three-cube Hirshfeld deformation density equals the vector version when the spherical
// cube holds the exact Thakkar sums, and both follow rho_A / sum rho * rho - rho_A around the
// chosen atom only; Calc_Hirshfeld_atom over both atoms partitions rho exactly
TEST(PropertiesPromolTests, HirshfeldDeformationAndAtomPartitionFollowStockholderWeights)
{
	const H2Model m(1.0, 1.0);
	const cube grid = make_grid(N, H);
	std::ostringstream log;
	const Thakkar hydrogen(1);

	std::vector<cube> cubes = property_cubes(grid, { cube_type::HDEF, cube_type::Hirsh, cube_type::spherical_density });
	Calc_Rho(cubes[cube_type::Rho], m.wavy, 100.0, log, false);
	cube &spher = cubes[cube_type::spherical_density];
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = grid.get_pos(x, y, z);
				spher.set_value(x, y, z, hydrogen.get_radial_density(dist(p, m.wavy.get_atom_pos(0))) + hydrogen.get_radial_density(dist(p, m.wavy.get_atom_pos(1))));
			}
	Calc_Hirshfeld(cubes, m.wavy, RADIUS, 0, log, false);
	cube hdef3 = make_grid(N, H);
	Calc_Hirshfeld(hdef3, cubes[cube_type::Rho], spher, m.wavy, RADIUS, 0, log, false);

	const d3 A = m.wavy.get_atom_pos(0);
	int checked = 0;
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = grid.get_pos(x, y, z);
				const double v = cubes[cube_type::HDEF].get_value(x, y, z);
				EXPECT_NEAR(hdef3.get_value(x, y, z), v, 1e-12);
				if (dist(p, A) >= constants::ang2bohr(RADIUS))
				{
					EXPECT_EQ(v, 0.0);
					continue;
				}
				const double rhoA = hydrogen.get_radial_density(dist(p, A));
				EXPECT_NEAR(v, rhoA / spher.get_value(x, y, z) * m.rho(p) - rhoA, 1e-12);
				checked++;
			}
	EXPECT_GT(checked, 100);
	// the far nucleus is more than the radius from atom 0 and stays zero
	EXPECT_EQ(cubes[cube_type::HDEF].get_value(12, 8, 8), 0.0);
	EXPECT_NE(cubes[cube_type::HDEF].get_value(4, 8, 8), 0.0);

	Calc_Hirshfeld_atom(cubes, m.wavy, 100.0, 0, log, false);
	const cube atom0 = cubes[cube_type::Hirsh];
	Calc_Hirshfeld_atom(cubes, m.wavy, 100.0, 1, log, false);
	const cube &atom1 = cubes[cube_type::Hirsh];
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = grid.get_pos(x, y, z);
				EXPECT_NEAR(atom0.get_value(x, y, z) + atom1.get_value(x, y, z), m.rho(p), 1e-12);
				EXPECT_GT(atom0.get_value(x, y, z), 0.0);
			}
	EXPECT_NEAR(atom0.get_value(4, 8, 8), atom1.get_value(12, 8, 8), 1e-12);
	EXPECT_GT(atom0.get_value(4, 8, 8), atom1.get_value(4, 8, 8));
}

// the seed finder wants at least three voxels per axis, a floor above the maximum leaves no
// seed, and the default floor keeps the single stationary voxel of a gaussian centred on a
// grid point with a vanishing central-difference gradient
TEST(PropertiesCriticalPointTests, SeedFinderRejectsTinyCubesAndFindsGaussianMaximum)
{
	cube tiny({ 2, 5, 5 }, 0, true);
	EXPECT_TRUE(find_cube_critical_point_seeds(&tiny, false).empty());

	cube g = make_grid(11, 0.5);
	for (int x = 0; x < 11; x++)
		for (int y = 0; y < 11; y++)
			for (int z = 0; z < 11; z++)
			{
				const d3 p = g.get_pos(x, y, z);
				g.set_value(x, y, z, std::exp(-(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])));
			}
	const std::vector<critical_point_seed> seeds = find_cube_critical_point_seeds(&g, false);
	ASSERT_EQ(seeds.size(), 1u);
	EXPECT_EQ(seeds[0].grid_index, (i3{ 5, 5, 5 }));
	EXPECT_NEAR(seeds[0].position[0], 0.0, 1e-12);
	EXPECT_NEAR(seeds[0].value, 1.0, 1e-12);
	EXPECT_NEAR(seeds[0].gradient_norm, 0.0, 1e-12);
	EXPECT_FALSE(seeds[0].is_nuclear_seed);

	EXPECT_TRUE(find_cube_critical_point_seeds(&g, false, 2.0).empty());
}

// with the gaussian shifted off the grid the bracketing voxel has a finite gradient: the default
// cutoff (none) keeps it and an explicit epsilon below that gradient drops it
TEST(PropertiesCriticalPointTests, SeedFinderExplicitGradientCutoffDropsOffGridMaximum)
{
	cube g = make_grid(11, 0.5);
	for (int x = 0; x < 11; x++)
		for (int y = 0; y < 11; y++)
			for (int z = 0; z < 11; z++)
			{
				const d3 p = g.get_pos(x, y, z);
				g.set_value(x, y, z, std::exp(-((p[0] - 0.1) * (p[0] - 0.1) + p[1] * p[1] + p[2] * p[2])));
			}
	const std::vector<critical_point_seed> loose = find_cube_critical_point_seeds(&g, false);
	ASSERT_EQ(loose.size(), 1u);
	EXPECT_EQ(loose[0].grid_index, (i3{ 5, 5, 5 }));
	EXPECT_GT(loose[0].gradient_norm, 1e-3);
	const std::vector<critical_point_seed> strict = find_cube_critical_point_seeds(&g, false, -1.0, 1e-30);
	EXPECT_TRUE(strict.empty());
	const std::vector<critical_point_seed> wide = find_cube_critical_point_seeds(&g, false, -1.0, 1.0);
	EXPECT_EQ(wide.size(), 1u);
}

// the analytic refinement on the closer, tighter H2 (a bonded pair by the covalent radii, so
// bond seeds are injected) finds the two nuclear attractors and the (3,-1) bond point at the
// midpoint with the closed-form density, laplacian and degenerate curvatures
TEST(PropertiesCriticalPointTests, AnalyzeCriticalPointsFindsAttractorsAndBondPoint)
{
	const H2Model m(2.0, 0.7);
	cube rho = make_grid(21, 0.2);
	std::ostringstream log;
	Calc_Rho(rho, m.wavy, 100.0, log, false);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, m.wavy, false);

	int attractors = 0, bonds = 0;
	for (const critical_point &cp : cps)
	{
		if (cp.type == "attractor")
		{
			attractors++;
			EXPECT_EQ(cp.negative_eigenvalues, 3);
			EXPECT_TRUE(std::abs(cp.position[0]) > 0.6 && std::abs(cp.position[0]) < 0.8) << cp.position[0];
			EXPECT_NEAR(cp.position[1], 0.0, 1e-4);
			EXPECT_NEAR(cp.position[2], 0.0, 1e-4);
			EXPECT_NEAR(cp.density, m.rho(cp.position), 1e-10);
		}
		else if (cp.type == "bond")
		{
			bonds++;
			EXPECT_TRUE(cp.converged);
			EXPECT_EQ(cp.negative_eigenvalues, 2);
			EXPECT_EQ(cp.positive_eigenvalues, 1);
			EXPECT_EQ(cp.zero_eigenvalues, 0);
			EXPECT_NEAR(cp.position[0], 0.0, 1e-5);
			EXPECT_NEAR(cp.position[1], 0.0, 1e-5);
			EXPECT_NEAR(cp.position[2], 0.0, 1e-5);
			EXPECT_NEAR(cp.density, m.rho({ 0.0, 0.0, 0.0 }), 1e-8);
			// rho = 2 cp^2 f^2 with f = gA + gB, grad f = 0 at the midpoint:
			// hessian = 4 cp^2 f H(f), H(f) = 2 g diag(4 a^2 R^2 - 2 a, -2 a, -2 a), g = e^{-a R^2}
			const double g = std::exp(-m.a * m.R * m.R);
			const double pref = 4.0 * m.cp * m.cp * 2.0 * g * 2.0 * g;
			const double ax = pref * (4.0 * m.a * m.a * m.R * m.R - 2.0 * m.a);
			const double tr = pref * (-2.0 * m.a);
			EXPECT_NEAR(cp.laplacian, ax + 2.0 * tr, 1e-7);
			EXPECT_NEAR(cp.hessian_eigenvalues[0] + cp.hessian_eigenvalues[1] + cp.hessian_eigenvalues[2], ax + 2.0 * tr, 1e-7);
			EXPECT_NEAR(cp.ellipticity, 0.0, 1e-6);
		}
		EXPECT_NEAR(cp.gradient_norm, 0.0, 1e-6);
	}
	EXPECT_EQ(attractors, 2);
	EXPECT_EQ(bonds, 1);
	EXPECT_EQ(cps.size(), 3u);
	// a floor above the bond density keeps the nuclear seeds (they skip the floor) and drops
	// the bond point and the bond seeds sitting below it
	const std::vector<critical_point> floored = analyze_cube_critical_points(&rho, m.wavy, false, 1.5 * m.rho({ 0.0, 0.0, 0.0 }));
	int floored_attractors = 0;
	for (const critical_point &cp : floored)
	{
		EXPECT_NE(cp.type, "bond");
		if (cp.type == "attractor")
			floored_attractors++;
	}
	EXPECT_EQ(floored_attractors, 2);
}

// near-grid ascent on the H2 density gives two mirror-image basins whose maxima sit on the
// nuclei; integrate_values_in_basins (with its debug line) splits the cube sum in half and the
// QTAIM labels name the nearest atom
TEST(PropertiesBasinTests, NearGridAscentSplitsH2IntoTwoMirrorBasins)
{
	const H2Model m(1.0, 1.0);
	cube rho = make_grid(N, H);
	std::ostringstream log;
	Calc_Rho(rho, m.wavy, 100.0, log, false);
	const std::vector<atom> atoms = m.wavy.get_atoms();

	std::pair<cubei, std::vector<d4>> result = topological_cube_analysis(&rho, atoms, false, false, 0.0, 0.0);
	const cubei &basins = result.first;
	const std::vector<d4> &maxima = result.second;
	ASSERT_EQ(maxima.size(), 2u);
	EXPECT_EQ(basins.max_value(), 2);
	for (const d4 &mx : maxima)
	{
		EXPECT_NEAR(std::abs(mx[0]), 1.0, 1e-12);
		EXPECT_NEAR(mx[1], 0.0, 1e-12);
		EXPECT_NEAR(mx[2], 0.0, 1e-12);
		EXPECT_NEAR(mx[3], m.rho({ mx[0], 0.0, 0.0 }), 1e-12);
	}
	const int left = basins.get_value(4, 8, 8), right = basins.get_value(12, 8, 8);
	EXPECT_NE(left, right);
	EXPECT_NE(left, 0);
	EXPECT_NE(right, 0);
	int mirrored = 0, total = 0;
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const int b = basins.get_value(x, y, z);
				EXPECT_NE(b, 0);
				const int mirror = basins.get_value(N - 1 - x, y, z);
				if (x != 8 && (b == left ? mirror == right : mirror == left))
					mirrored++;
				total++;
			}
	EXPECT_EQ(mirrored, total - N * N);

	svec labels = assign_labels_to_basins(maxima, atoms, false, 0);
	ASSERT_EQ(labels.size(), 2u);
	const bool first_left = maxima[0][0] < 0.0;
	EXPECT_EQ(labels[0], first_left ? "H0" : "H1");
	EXPECT_EQ(labels[1], first_left ? "H1" : "H0");

	const vec eds = integrate_values_in_basins(&rho, &basins, labels, true);
	ASSERT_EQ(eds.size(), 2u);
	// the x = 0 plane (index 8) is a gradient tie and goes to one basin by tie-break; off the plane
	// the two basins are mirror images
	double plane[2] = { 0.0, 0.0 };
	for (int y = 0; y < N; y++)
		for (int z = 0; z < N; z++)
			plane[basins.get_value(8, y, z) - 1] += rho.get_value(8, y, z) * rho.get_dv();
	EXPECT_NEAR(eds[0] - plane[0], eds[1] - plane[1], 1e-12);
	EXPECT_NEAR(eds[0] + eds[1], rho.sum(), 1e-12);
	// the 2 bohr box holds most of the two electrons
	EXPECT_GT(eds[0] + eds[1], 1.8);
}

// seeds own their voxel and their position is reported as the maximum even when the grid
// maximum sits one voxel away; an off-nucleus maximum away from any atom is an NNA in the labels
TEST(PropertiesBasinTests, SeedsOwnTheirBasinsAndOffNucleusMaximaAreLabelledNna)
{
	const H2Model m(1.0, 1.0);
	cube rho = make_grid(N, H);
	std::ostringstream log;
	Calc_Rho(rho, m.wavy, 100.0, log, false);
	const std::vector<atom> atoms = m.wavy.get_atoms();

	const std::vector<d3> seeds{ d3{ -1.0, 0.15, 0.0 }, d3{ 1.0, -0.15, 0.0 } };
	std::pair<cubei, std::vector<d4>> result = topological_cube_analysis(&rho, atoms, false, false, 0.0, 0.0, -1.0, 5e-3, &seeds);
	ASSERT_EQ(result.second.size(), 2u);
	EXPECT_NEAR(result.second[0][0], -1.0, 1e-12);
	EXPECT_NEAR(result.second[0][1], 0.15, 1e-12);
	EXPECT_NEAR(result.second[1][0], 1.0, 1e-12);
	EXPECT_NEAR(result.second[1][1], -0.15, 1e-12);
	// the seed's value is that of the voxel it rounds to, one step off the nucleus
	EXPECT_NEAR(result.second[0][3], rho.get_value(4, 9, 8), 1e-12);
	EXPECT_EQ(result.first.get_value(4, 8, 8), 1);
	EXPECT_EQ(result.first.get_value(4, 9, 8), 1);
	EXPECT_EQ(result.first.get_value(12, 8, 8), 2);
	EXPECT_EQ(result.first.max_value(), 2);

	// a seed outside the grid is ignored
	const std::vector<d3> outside{ d3{ 5.0, 0.0, 0.0 } };
	std::pair<cubei, std::vector<d4>> ignored = topological_cube_analysis(&rho, atoms, false, false, 0.0, 0.0, -1.0, 5e-3, &outside);
	EXPECT_EQ(ignored.second.size(), 2u);

	const std::vector<d4> maxima{ d4{ -1.0, 0.0, 0.0, 1.0 }, d4{ 0.0, 1.2, 0.0, 0.1 } };
	const svec labels = assign_labels_to_basins(maxima, atoms, false, 0);
	ASSERT_EQ(labels.size(), 2u);
	EXPECT_EQ(labels[0], "H0");
	EXPECT_EQ(labels[1], "NNA near H0");
}

// with the analytic gradient and no seeds every trajectory is unresolved and the highest voxel
// it crossed becomes a maximum; a trajectory ending within a voxel and a half of an existing
// one joins it, so H2 still comes out as two basins on the nuclei
TEST(PropertiesBasinTests, GradientTrajectoriesWithoutSeedsCreateMaximaAtNuclei)
{
	const H2Model m(1.0, 1.0);
	cube rho = make_grid(N, H);
	std::ostringstream log;
	Calc_Rho(rho, m.wavy, 100.0, log, false);
	const std::vector<atom> atoms = m.wavy.get_atoms();

	std::pair<cubei, std::vector<d4>> result = topological_cube_analysis(&rho, atoms, false, false, 0.0, 0.0, -1.0, 5e-3, nullptr, &m.wavy);
	ASSERT_EQ(result.second.size(), 2u);
	for (const d4 &mx : result.second)
	{
		EXPECT_NEAR(std::abs(mx[0]), 1.0, 1e-12);
		EXPECT_NEAR(mx[1], 0.0, 1e-12);
		EXPECT_NEAR(mx[2], 0.0, 1e-12);
	}
	const cubei &basins = result.first;
	EXPECT_NE(basins.get_value(4, 8, 8), basins.get_value(12, 8, 8));
	EXPECT_EQ(basins.max_value(), 2);
	int count[3] = { 0, 0, 0 };
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
				count[basins.get_value(x, y, z)]++;
	EXPECT_EQ(count[0], 0);
	// mirror symmetry up to the bisector plane, which either side may take
	EXPECT_NEAR(count[1], count[2], N * N);
	EXPECT_GT(count[1], N * N * 7);
	EXPECT_GT(count[2], N * N * 7);
}

// an assignment radius keeps only voxels near an atom (basin 0 elsewhere) and a value floor
// removes the faint tail; a one-voxel bump in the tail is a basin of its own with the merge
// switched off and folds into the nuclear basin at a persistence of 5e-3
TEST(PropertiesBasinTests, AssignmentRadiusFloorAndPersistenceMerge)
{
	const H2Model m(1.0, 1.0);
	cube rho = make_grid(N, H);
	std::ostringstream log;
	Calc_Rho(rho, m.wavy, 100.0, log, false);
	const std::vector<atom> atoms = m.wavy.get_atoms();

	std::pair<cubei, std::vector<d4>> limited = topological_cube_analysis(&rho, atoms, false, false, 0.0, 0.0, 0.5);
	EXPECT_EQ(limited.second.size(), 2u);
	EXPECT_EQ(limited.first.get_value(0, 0, 0), 0);
	EXPECT_EQ(limited.first.get_value(8, 8, 8), 0);
	EXPECT_NE(limited.first.get_value(4, 8, 8), 0);
	const double r_bohr = constants::ang2bohr(0.5);
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++)
			for (int z = 0; z < N; z++)
			{
				const d3 p = rho.get_pos(x, y, z);
				const bool close_to_atom = dist(p, atoms[0].get_pos()) <= r_bohr || dist(p, atoms[1].get_pos()) <= r_bohr;
				EXPECT_EQ(limited.first.get_value(x, y, z) != 0, close_to_atom);
			}

	const double floor = rho.get_value(8, 8, 8);
	std::pair<cubei, std::vector<d4>> floored = topological_cube_analysis(&rho, atoms, false, false, floor, 0.0);
	EXPECT_EQ(floored.second.size(), 2u);
	EXPECT_EQ(floored.first.get_value(8, 8, 8), 0);
	EXPECT_EQ(floored.first.get_value(0, 0, 0), 0);
	EXPECT_NE(floored.first.get_value(4, 8, 8), 0);

	// bump at (2.0, 0, 0) beyond the right nucleus, a tenth of a percent over its highest neighbour
	cube bumped = rho;
	double best = 0.0;
	for (int dx = -1; dx <= 1; dx++)
		for (int dy = -1; dy <= 1; dy++)
			for (int dz = -1; dz <= 1; dz++)
				if (dx || dy || dz)
					best = std::max(best, rho.get_value(16 + dx < N ? 16 + dx : 16, 8 + dy, 8 + dz));
	bumped.set_value(16, 8, 8, best * 1.001);
	std::pair<cubei, std::vector<d4>> kept = topological_cube_analysis(&bumped, atoms, false, false, 0.0, 0.0, -1.0, 0.0);
	ASSERT_EQ(kept.second.size(), 3u);
	bool bump_found = false;
	for (const d4 &mx : kept.second)
		if (std::abs(mx[0] - 2.0) < 1e-12 && std::abs(mx[1]) < 1e-12)
			bump_found = true;
	EXPECT_TRUE(bump_found);
	EXPECT_EQ(kept.first.max_value(), 3);
	std::pair<cubei, std::vector<d4>> merged = topological_cube_analysis(&bumped, atoms, false, false, 0.0, 0.0, -1.0, 5e-3);
	EXPECT_EQ(merged.second.size(), 2u);
	EXPECT_EQ(merged.first.max_value(), 2);
	EXPECT_EQ(merged.first.get_value(16, 8, 8), merged.first.get_value(12, 8, 8));
	//This bump's persistence is (1.001 - 1) / 1.001 = 9.99e-4 of its height, which is the scale a
	//shard of a flat valence shell sits at - and 5e-3 is wide enough to eat it. Inside a shell every
	//saddle is about as deep as the one down to the core, so single linkage then chains the shards
	//INTO the core basin and the core reports whole electrons too many (Cl2's chlorine 14.8951 e
	//against the 10 its closed shells hold). The ELI-D call site therefore passes 3e-4 now, and the
	//length-based unify_shell_basins folds the shell instead; this asserts the boundary the constant
	//has to stay on the right side of, because nothing else in the suite would notice it moving back.
	std::pair<cubei, std::vector<d4>> shipped = topological_cube_analysis(&bumped, atoms, false, false, 0.0, 0.0, -1.0, 3e-4);
	EXPECT_EQ(shipped.second.size(), 3u) << "at 3e-4 a bump 9.99e-4 above its saddle survives";
	EXPECT_NE(shipped.first.get_value(16, 8, 8), shipped.first.get_value(12, 8, 8));
}

// ELI labels: the proton's basin by its nucleus, a maximum inside the core shell of a heavier
// atom is a core, one much closer to one heavy atom than the other a lone pair, otherwise a bond
TEST(PropertiesBasinTests, EliLabelsCoreLonePairBondAndProton)
{
	std::vector<atom> atoms;
	atoms.emplace_back("C", atomID(), 1, 0.0, 0.0, 0.0, 6);
	atoms.emplace_back("O", atomID(), 2, 2.2, 0.0, 0.0, 8);
	atoms.emplace_back("H", atomID(), 3, -2.0, 0.0, 0.0, 1);
	const std::vector<d4> maxima{
		d4{ 0.1, 0.0, 0.0, 50.0 },   // 0.1 bohr from C: core (radius 0.25)
		d4{ 1.1, 0.0, 0.0, 2.0 },    // midpoint of C-O: bond
		d4{ 3.0, 0.0, 0.0, 1.5 },    // 0.8 bohr past O, 3.0 from C: ratio 0.07, lone pair
		d4{ -1.6, 0.0, 0.0, 0.9 },   // 0.4 bohr from H: the proton's basin
		d4{ 2.1, 0.0, 0.0, 40.0 }    // 0.1 bohr from O: core (radius 0.55)
	};
	const svec labels = assign_labels_to_basins(maxima, atoms, false, 1);
	ASSERT_EQ(labels.size(), 5u);
	EXPECT_EQ(labels[0], "C0 core");
	EXPECT_EQ(labels[1], "C0-O1 bond");
	EXPECT_EQ(labels[2], "O1 LP");
	EXPECT_EQ(labels[3], "H2");
	EXPECT_EQ(labels[4], "O1 core");
}

// core_shell_radius steps by period; unify_core_basins folds every maximum inside an atom's
// core radius into one basin per atom keeping the highest, renumbers the cube and reports the
// number merged
TEST(PropertiesBasinTests, UnifyCoreBasinsMergesMaximaInsideTheCoreRadius)
{
	EXPECT_EQ(core_shell_radius(1), 0.0);
	EXPECT_EQ(core_shell_radius(2), 0.0);
	EXPECT_EQ(core_shell_radius(6), 0.25);
	EXPECT_EQ(core_shell_radius(10), 0.25);
	//Na-Ar shares the 1.0 bohr band with K-Kr: 0.55 was measured to sit INSIDE the L shell at the
	//electropositive end of the row (Na's L-shell ELI-D maximum is 0.740 bohr out, Al's 0.582), so
	//about 7 e of a 10 e core stayed unfolded. The gap between the furthest maximum that must fold
	//in (0.740) and the nearest that must not (1.472) is 0.732 bohr wide.
	EXPECT_EQ(core_shell_radius(11), 1.0);
	EXPECT_EQ(core_shell_radius(17), 1.0);
	EXPECT_EQ(core_shell_radius(26), 1.0);
	EXPECT_EQ(core_shell_radius(53), 1.4);
	EXPECT_EQ(core_shell_radius(82), 1.8);

	std::vector<atom> atoms;
	atoms.emplace_back("C", atomID(), 1, 0.0, 0.0, 0.0, 6);
	atoms.emplace_back("H", atomID(), 2, 2.0, 0.0, 0.0, 1);
	std::vector<d4> maxima{
		d4{ 0.1, 0.0, 0.0, 30.0 },   // inside the C core
		d4{ 2.0, 0.0, 0.0, 0.4 },    // on H, whose core radius is 0 - never merged
		d4{ -0.1, 0.0, 0.0, 31.0 },  // inside the C core too, the higher one
		d4{ 1.0, 0.0, 0.0, 0.3 }     // bond region, outside every core
	};
	cubei basins({ 4, 1, 1 }, 0, true);
	for (int b = 0; b < 4; b++)
		basins.set_value(b, 0, 0, b + 1);
	const int merged = unify_core_basins(basins, maxima, atoms);
	EXPECT_EQ(merged, 1);
	ASSERT_EQ(maxima.size(), 3u);
	EXPECT_NEAR(maxima[0][3], 31.0, 0.0);
	EXPECT_NEAR(maxima[0][0], -0.1, 0.0);
	EXPECT_NEAR(maxima[1][3], 0.4, 0.0);
	EXPECT_NEAR(maxima[2][3], 0.3, 0.0);
	EXPECT_EQ(basins.get_value(0, 0, 0), 1);
	EXPECT_EQ(basins.get_value(2, 0, 0), 1);
	EXPECT_EQ(basins.get_value(1, 0, 0), 2);
	EXPECT_EQ(basins.get_value(3, 0, 0), 3);
	EXPECT_EQ(basins.max_value(), 3);

	std::vector<d4> apart{ d4{ 0.1, 0.0, 0.0, 30.0 }, d4{ 1.0, 0.0, 0.0, 0.3 } };
	cubei two({ 2, 1, 1 }, 0, true);
	two.set_value(0, 0, 0, 1);
	two.set_value(1, 0, 0, 2);
	EXPECT_EQ(unify_core_basins(two, apart, atoms), 0);
	EXPECT_EQ(apart.size(), 2u);
	EXPECT_EQ(two.get_value(1, 0, 0), 2);
}

// Outside the cores the same sphere of maxima appears with nothing to fold it: a spherically symmetric
// ELI-D shell sampled on a cubic grid is handed out one basin per voxel, and Co2 - two atoms - kept 845
// basins that way, 824 of them under 0.01 e. The persistence merge cannot fix it: swept from 5e-3 to
// 2e-1, the first threshold that dented Co2 at all (3e-2, 845 -> 687) already took one of OH's two REAL
// oxygen lone pairs, because both a shell's grid saddles and the saddle between two genuine lone pairs
// are shallow. Only a length separates them, and these are the measured geometries.
TEST(PropertiesBasinTests, UnifyShellBasinsFoldsAShatteredShellAndKeepsTwoRealLonePairs)
{
	std::vector<d4> maxima;
	// Co2's shell: one sphere at 5.35 bohr with neighbours 0.378 bohr apart (two voxels at 0.1 A) and
	// values alternating 2.2671 / 2.3067 as measured - 1.7 % apart, so no value test keeps them together
	const int n_shell = 89;  // 2 pi 5.35 / 0.378
	for (int i = 0; i < n_shell; i++) {
		const double a = constants::TWO_PI * i / n_shell;
		maxima.push_back(d4{ 5.35 * std::cos(a), 5.35 * std::sin(a), 0.0, i % 2 ? 2.3067 : 2.2671 });
	}
	// OH's two REAL oxygen lone pairs: 1.890 bohr apart, degenerate to 0.2 %
	maxima.push_back(d4{ 2.268, -0.945, -0.378, 1.6704 });
	maxima.push_back(d4{ 2.268, 0.945, -0.378, 1.6671 });
	// ZP2's duplicated F1 lone pair: 0.84 bohr apart, degenerate to 0.2 %, holding 1.4276 and 1.2907 e
	// where one lone pair holds about 2.7. The prediction in the other direction - these MUST merge.
	maxima.push_back(d4{ 0.0, 0.0, 12.0, 1.6460 });
	maxima.push_back(d4{ 0.0, 0.84, 12.0, 1.6427 });

	const std::vector<d4> before = maxima;
	const int nb = static_cast<int>(maxima.size());
	cubei basins({ nb, 1, 1 }, 0, true);
	for (int b = 0; b < nb; b++)
		basins.set_value(b, 0, 0, b + 1);
	ivec map;
	const int merged = unify_shell_basins(basins, maxima, &map);
	EXPECT_EQ(merged, n_shell);  // 88 of the shell, plus one of ZP2's pair
	ASSERT_EQ(maxima.size(), 4u);

	// the shell keeps its highest maximum, and it is still on the sphere
	EXPECT_NEAR(maxima[0][3], 2.3067, 0.0);
	EXPECT_NEAR(std::sqrt(maxima[0][0] * maxima[0][0] + maxima[0][1] * maxima[0][1]), 5.35, 1e-9);
	// both oxygen lone pairs survive, separately
	EXPECT_NEAR(maxima[1][3], 1.6704, 0.0);
	EXPECT_NEAR(maxima[2][3], 1.6671, 0.0);
	EXPECT_NEAR(maxima[1][1], -0.945, 0.0);
	EXPECT_NEAR(maxima[2][1], 0.945, 0.0);
	// ZP2's duplicate is one basin now, keeping the higher maximum
	EXPECT_NEAR(maxima[3][3], 1.6460, 0.0);

	// the cube and the basin map agree: the whole shell is basin 1, the lone pairs are 2 and 3
	ASSERT_EQ(map.size(), static_cast<size_t>(nb) + 1);
	for (int b = 0; b < n_shell; b++) {
		EXPECT_EQ(basins.get_value(b, 0, 0), 1) << "shell voxel " << b;
		EXPECT_EQ(map[b + 1], 1) << "shell maximum " << b;
	}
	EXPECT_EQ(map[n_shell + 1], 2);
	EXPECT_EQ(map[n_shell + 2], 3);
	EXPECT_EQ(map[n_shell + 3], 4);
	EXPECT_EQ(map[n_shell + 4], 4);
	EXPECT_EQ(basins.get_value(n_shell + 3, 0, 0), 4);
	EXPECT_EQ(basins.max_value(), 4);

	// a distance of zero is the off switch and must change nothing at all
	std::vector<d4> untouched = before;
	cubei same({ nb, 1, 1 }, 0, true);
	for (int b = 0; b < nb; b++)
		same.set_value(b, 0, 0, b + 1);
	EXPECT_EQ(unify_shell_basins(same, untouched, nullptr, 0.0, 0.05), 0);
	EXPECT_EQ(untouched.size(), before.size());
	EXPECT_EQ(same.max_value(), nb);

	// and the margin: the default 1.2 bohr sits between ZP2's 0.84 and OH's 1.890, so a cutoff past
	// 1.890 eats a real lone pair. That is the failure this test exists to catch.
	std::vector<d4> too_far = before;
	cubei wide({ nb, 1, 1 }, 0, true);
	for (int b = 0; b < nb; b++)
		wide.set_value(b, 0, 0, b + 1);
	EXPECT_EQ(unify_shell_basins(wide, too_far, nullptr, 2.0, 0.05), n_shell + 1);
	ASSERT_EQ(too_far.size(), 3u);
}

// the legacy interactive b2c() with its selection read from a redirected cin writes the log and
// the selected-basins cube next to the input cube
TEST(PropertiesBasinTests, LegacyB2cWritesLogAndSelectedBasinCube)
{
	const std::filesystem::path dir = temp_dir("legacy_b2c");
	const H2Model m(1.0, 1.0);
	cube rho = make_grid(9, 0.5);
	std::ostringstream log;
	Calc_Rho(rho, m.wavy, 100.0, log, false);
	rho.set_path(dir / "h2.cube");

	std::istringstream answers("1\n0\n");
	std::streambuf *old_cin = std::cin.rdbuf(answers.rdbuf());
	const bool ok = b2c(&rho, m.wavy.get_atoms(), false, false);
	std::cin.rdbuf(old_cin);
	EXPECT_TRUE(ok);
	EXPECT_TRUE(std::filesystem::exists(dir / "h2.b2c_log"));
	EXPECT_TRUE(std::filesystem::exists(dir / "h2_1_basins.cube"));
	std::filesystem::remove_all(dir);
}

// -def on its own goes through properties_calculation: the spherical density must be computed
// for it too, otherwise the "deformation" cube is a copy of rho (it was, until 20 Sep 2026)
TEST(PropertiesDriverTests, StaticDeformationAloneSubtractsTheSphericalAtoms)
{
	const auto input = nos_test_repo_root() / "tests" / "epoxide_gbw" / "epoxide.gbw";
	if (!std::filesystem::exists(input)) GTEST_SKIP() << "Missing " << input;
	const std::filesystem::path dir = temp_dir("def_alone");
	std::filesystem::copy_file(input, dir / "epoxide.gbw");
	options opt;
	opt.wfn = (dir / "epoxide.gbw").string();
	opt.properties.def = true;
	opt.properties.resolution = 0.5;
	opt.properties.radius = 1.5;
	properties_calculation(opt);
	std::ostringstream log;
	WFN wave(dir / "epoxide.gbw", false);
	const cube rho(dir / "epoxide_rho.cube", true, wave, log);
	const cube def(dir / "epoxide_def.cube", true, wave, log);
	const double rho_max = rho.max_value();
	const double def_max = std::max(def.max_value(), -def.min_value());
	EXPECT_GT(rho_max, 1.0);
	// the Thakkar atoms take the nuclear peaks out; what is left is the bonding rearrangement
	EXPECT_GT(def_max, 0.02);
	EXPECT_LT(def_max, 0.1 * rho_max);
	std::filesystem::remove_all(dir);
}
