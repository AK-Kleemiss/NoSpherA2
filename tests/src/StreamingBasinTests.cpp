#include "pch.h"
#include "core/b2c.h"
#include "core/wfn_class.h"
#include "core/atoms.h"
#include "core/constants.h"
#include "core/cube.h"
#include "core/nos_math.h"
#include "core/properties.h"

#include <cmath>
#include <filesystem>
#include <sstream>

//The basin analysis used to take its attractors from a cube: wherever a trajectory started on a
//voxel ended nowhere, a maximum was declared there. That invents non-nuclear attractors out of
//the sampling - epoxide at 0.05 A grew one of 0.027 e that no critical point of the analytic
//density lies under - and it is the reason a finer grid used to produce more debris rather than
//less. The streaming path asks the field instead, and these tests pin the question it asks.
namespace
{
	//Two spherical bumps a fixed distance apart: maxima at the centres, one saddle exactly
	//between them, and nothing else anywhere. Everything the criterion has to separate, with
	//the answers known in closed form
	scalar_field two_bumps(const d3 &c1, const d3 &c2, const double a)
	{
		return [c1, c2, a](const d3 &p, d3 &g) {
			const d3 d1{ p[0] - c1[0], p[1] - c1[1], p[2] - c1[2] };
			const d3 d2{ p[0] - c2[0], p[1] - c2[1], p[2] - c2[2] };
			const double e1 = std::exp(-a * (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]));
			const double e2 = std::exp(-a * (d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]));
			for (int k = 0; k < 3; k++) g[k] = -2.0 * a * (d1[k] * e1 + d2[k] * e2);
			return e1 + e2;
		};
	}

	WFN load(const std::filesystem::path &p) { return WFN(p); }

	//A coarse cube, only so the critical-point search has seeds to start from. The point of the
	//streaming path is that this resolution no longer decides the answer
	cube seed_cube(const WFN &wavy, const double spacing, const double radius)
	{
		const std::vector<atom> atoms = wavy.get_atoms();
		const double pad = constants::ang2bohr(radius);
		d3 lo{ 1e30, 1e30, 1e30 }, hi{ -1e30, -1e30, -1e30 };
		for (const atom &a : atoms) {
			const d3 p = a.get_pos();
			for (int k = 0; k < 3; k++) {
				lo[k] = std::min(lo[k], p[k] - pad);
				hi[k] = std::max(hi[k], p[k] + pad);
			}
		}
		std::array<int, 3> n{};
		for (int k = 0; k < 3; k++) n[k] = static_cast<int>((hi[k] - lo[k]) / spacing) + 1;
		cube rho(n, 0, true);
		for (int k = 0; k < 3; k++) {
			rho.set_origin(k, lo[k]);
			rho.set_vector(k, k, spacing);
		}
		rho.calc_dv();
		std::ostringstream log;
		Calc_Rho(rho, wavy, radius, log, false);
		return rho;
	}
}

//The criterion itself, where the answers are known: a bump's centre is a maximum, the saddle
//between two bumps is not, and a point on a flank belongs to the bump it is on
TEST(StreamingBasins, MaximumCriterionSeparatesMaximaFromSaddles)
{
	//Six bohr apart, so each bump pulls the other's maximum off its centre by ~1e-7 and the
	//centres are the maxima to the precision asked for below
	const d3 c1{ 0.0, 0.0, 0.0 }, c2{ 6.0, 0.0, 0.0 };
	const scalar_field f = two_bumps(c1, c2, 0.5);

	d3 p{ 0.2, 0.1, -0.1 };
	ASSERT_TRUE(converge_to_maximum(f, p)) << "a point beside a bump has to find its top";
	for (int k = 0; k < 3; k++) EXPECT_NEAR(p[k], c1[k], 1e-3);

	//The saddle at the midpoint satisfies grad = 0, which is exactly what a gradient-only test
	//would accept; the Hessian has one positive eigenvalue along the axis and rejects it
	d3 saddle{ 3.0, 0.0, 0.0 };
	EXPECT_FALSE(converge_to_maximum(f, saddle)) << "a saddle is not an attractor";

	//not named "far": windows.h still defines that as an empty macro
	d3 flank{ 5.6, 0.3, 0.0 };
	ASSERT_TRUE(converge_to_maximum(f, flank));
	for (int k = 0; k < 3; k++) EXPECT_NEAR(flank[k], c2[k], 1e-3);
}

//Every nucleus is an attractor of the density and nothing else in OH- is. The old path at 0.1 A
//also found only these two; the difference is that this one cannot find more however the cube
//is chosen, because there is no cube in the answer
TEST(StreamingBasins, HydroxideHasExactlyTwoAttractors)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy = load(wfn);
	ASSERT_EQ(wavy.get_ncen(), 2);
	const cube rho = seed_cube(wavy, 0.25, 3.0);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
	const std::vector<d4> maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
	ASSERT_EQ(maxima.size(), 2u) << "a non-nuclear attractor appeared in a molecule that has none";
	for (int a = 0; a < wavy.get_ncen(); a++) {
		const d3 p = wavy.get_atom_pos(a);
		double best = 1e30;
		for (const d4 &m : maxima) best = std::min(best, array_length(p, d3{ m[0], m[1], m[2] }));
		EXPECT_LT(best, 1e-8) << "atom " << a << " is not among the attractors";
	}
}

//The local virial theorem is pointwise: it holds at every point of every density however that
//density is evaluated, so it costs nothing and catches what nothing else looks at. V used to be
//computed as K - G, which is L - the right magnitude at a bond critical point and the wrong
//sign - and it was printed that way in every critical-point block of every job.
TEST(StreamingBasins, CriticalPointEnergyDensitiesObeyTheLocalVirialTheorem)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy = load(wfn);
	const cube rho = seed_cube(wavy, 0.25, 3.0);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
	int checked = 0;
	for (const critical_point &cp : cps) {
		//The energy densities need ELF in (0,1) to be defined at all; where it is not, they are
		//left NaN on purpose and there is nothing to check
		if (!std::isfinite(cp.virial_field)) continue;
		const double scale = std::max(1.0, std::abs(cp.kinetic_hamiltonian) + std::abs(cp.kinetic_lagrangian));
		EXPECT_NEAR(cp.lagrangian_density, -0.25 * cp.laplacian, 1e-10 * std::max(1.0, std::abs(cp.laplacian)));
		EXPECT_NEAR(cp.kinetic_hamiltonian - cp.kinetic_lagrangian, cp.lagrangian_density, 1e-10 * scale);
		EXPECT_NEAR(cp.virial_field, -(cp.kinetic_hamiltonian + cp.kinetic_lagrangian), 1e-10 * scale);
		//The potential energy density is negative everywhere. This alone fails on the old code
		//at any bond critical point, where the Laplacian is negative and L is therefore positive
		EXPECT_LT(cp.virial_field, 0.0) << "the virial field came out positive at a " << cp.type;
		EXPECT_GE(cp.kinetic_lagrangian, 0.0) << "G is positive definite";
		checked++;
	}
	EXPECT_GT(checked, 0) << "no critical point carried energy densities, so nothing was checked";
}

//The debris case, asked directly. At 0.05 A the cube path grew a non-nuclear attractor at this
//position on epoxide, 0.027 e over 2.0 bohr^3, while the critical-point search at the same
//resolution found seven attractors and all seven were nuclei. Handed the same position as a
//candidate, the criterion has to throw it out and keep every nucleus.
TEST(StreamingBasins, EpoxideGridDebrisIsNotAnAttractor)
{
	const std::filesystem::path gbw = nos_test_repo_root() / "tests" / "epoxide_gbw" / "epoxide.gbw";
	if (!std::filesystem::exists(gbw)) GTEST_SKIP() << "fixture missing: " << gbw.string();
	const WFN wavy = load(gbw);
	ASSERT_EQ(wavy.get_ncen(), 7);

	critical_point debris{};
	debris.position = d3{ -0.595, 14.378, 3.275 };
	debris.type = "attractor";
	debris.converged = true;

	const std::vector<d4> maxima = streaming_density_attractors(wavy, { debris }, nullptr, nullptr, false);
	EXPECT_EQ(maxima.size(), static_cast<size_t>(wavy.get_ncen()))
		<< "the grid's debris maximum survived a test of the analytic field";
	for (const d4 &m : maxima)
		EXPECT_GT(array_length(d3{ m[0], m[1], m[2] }, debris.position), 0.1)
			<< "an attractor was kept at the debris position";
}

//End to end, with no cube in the integration at all: the populations of a real wavefunction have
//to add up to its electrons. Ten electrons in OH-, and a quadrature that loses them makes every
//charge and every delocalization index below it meaningless
TEST(StreamingBasins, HydroxideStreamsToTheRightElectronCount)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy = load(wfn);
	const cube rho = seed_cube(wavy, 0.25, 3.0);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
	const std::vector<d4> maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
	ASSERT_EQ(maxima.size(), 2u);

	vec volumes;
	double outside = 0.0;
	basin_overlaps ovl;
	const vec pop = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, 3, false, volumes, outside, nullptr, nullptr, 1, nullptr, &ovl);
	ASSERT_EQ(pop.size(), 2u);
	double total = outside;
	for (const double p : pop) total += p;
	EXPECT_NEAR(total, 10.0, 0.05) << "the streaming quadrature lost electrons";
	EXPECT_LT(outside, 0.01) << "electrons ended up in no basin at all";
	//The hydrogen of a hydroxide keeps well under an electron and the oxygen carries the rest;
	//a boundary put in the wrong place shows up here long before the total does
	EXPECT_NEAR(pop[0], 0.61, 0.05);
	EXPECT_NEAR(pop[1], 9.37, 0.05);

	//sum_A S^A = I is not something the code chooses: it follows from the basins tiling space,
	//so it measures the integration and nothing else
	const delocalization_result r = delocalization_indices(wavy, ovl);
	EXPECT_LT(r.identity_error, 0.02) << "the basins do not add up to the whole of space";
	for (size_t b = 0; b < r.lambda.size(); b++) {
		double half = 0.0;
		for (size_t p = 0; p < r.pairs.size(); p++)
			if (r.pairs[p][0] == static_cast<int>(b) || r.pairs[p][1] == static_cast<int>(b)) half += 0.5 * r.di[p];
		EXPECT_NEAR(r.lambda[b] + half - r.population[b], 0.0, 0.02) << "basin " << b + 1 << " breaks the sum rule";
	}
}

//The partition has to be exact even where the rule is not. Whatever the quadrature's own error
//on the total, the basins only decide who gets each cell's weight, so sum(populations) + outside
//is the plain weighted sum over the grid and cannot depend on how many basins there are. Putting
//a spurious attractor at the bond midpoint is the hardest version of that: it carves a third
//basin out of the middle of the molecule, so most cells near the bond become straddling cells and
//go through the bisection. If a split cell ever handed out more than its own weight - both sides
//rounded up, a boundary counted by the climb and again by the split - this total would rise with
//the extra basin. It must not move at all.
TEST(StreamingBasins, BasinPartitionConservesTheQuadratureWeight)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy = load(wfn);
	const cube rho = seed_cube(wavy, 0.25, 3.0);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
	std::vector<d4> maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
	ASSERT_EQ(maxima.size(), 2u);

	auto total_of = [&](const std::vector<d4> &m) {
		vec volumes;
		double outside = 0.0;
		const vec pop = integrate_basins_on_atomic_grids(nullptr, nullptr, m, wavy, 3, false, volumes, outside, nullptr, nullptr, 1, nullptr, nullptr);
		double t = outside;
		for (const double p : pop) t += p;
		return t;
	};

	const double two = total_of(maxima);
	const d3 a = wavy.get_atom_pos(0), b = wavy.get_atom_pos(1);
	maxima.push_back(d4{ 0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), 1.0 });
	const double three = total_of(maxima);

	//Relative, because the absolute size of the total is the rule's business and not this test's
	EXPECT_NEAR(three, two, 1e-9 * std::max(1.0, std::abs(two)))
		<< "adding a basin changed the integrated total by " << three - two << " electrons, so a cell's weight is not being conserved across the split";
}
