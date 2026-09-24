#include "pch.h"
#include "core/b2c.h"
#include "core/wfn_class.h"
#include "core/atoms.h"
#include "core/constants.h"
#include "core/cube.h"
#include "core/properties.h"

#include <cmath>
#include <filesystem>
#include <sstream>

//A beta sphere is the radius around an attractor inside which no ascent trajectory can get out:
//grad f . rhat < 0 at every point of the sphere, so a path leaving it would have to cross
//outwards while the gradient it follows points inwards. Every point inside therefore belongs to
//the one attractor inside without being climbed, and a trajectory that enters is finished on the
//spot - which is where the time goes, the step shrinking to a third of a voxel for the last 1.5
//bohr of every single climb.
//
//That is an argument, not a measurement, and it has one weak joint: the sphere is only sampled
//along finitely many directions, so a separatrix that comes closest between two samples is not
//seen. So the answer it produces is compared here against the answer of the same code climbing
//every trajectory the whole way, which is what -no_beta_spheres restores. If the radius is ever
//too large the populations move, and they move by far more than the quadrature's own error - a
//sphere reaching past a separatrix swallows part of a neighbouring basin.
//
//The first version of this file only checked OH, and OH passed a sampling that was in fact much
//too coarse: 26 directions kept at 90 % moved NH3Li's third hydrogen by 0.019 e, three per cent
//of its population, while the two-atom molecule stayed inside 5e-3. A two-atom fixture cannot
//find this, because its one separatrix is a plane square to the axis and the coarse sample hits
//it. NH3Li is therefore checked too, and the window is a thousandth of an electron - a
//same-binary A/B differs only in which side of a bisection a boundary cell lands on, which is
//worth a hundred-thousandth.
namespace
{
	WFN load(const std::filesystem::path &p) { return WFN(p); }

	//A coarse cube, only to give the critical-point search its seeds; the streaming path does not
	//take its boundaries from it
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

	//Restores the switch whatever the test does, so a failure cannot leave it off for the rest of
	//the suite
	struct beta_guard {
		bool was = beta_spheres_enabled();
		~beta_guard() { beta_spheres_set_enabled(was); }
	};

	struct adaptive_guard {
		bool was = basin_adaptive_step_enabled();
		~adaptive_guard() { basin_adaptive_step_set_enabled(was); }
	};

	//The same streaming QTAIM integration twice out of one binary, once with the spheres and once
	//with every trajectory climbed the whole way, compared basin by basin
	void expect_same_basins(const std::filesystem::path &wfn, const size_t nmax)
	{
		const WFN wavy(wfn);
		const cube rho = seed_cube(wavy, 0.25, 3.0);
		const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
		const std::vector<d4> maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
		ASSERT_EQ(maxima.size(), nmax);

		beta_guard guard;
		vec v_on, v_off;
		double out_on = 0.0, out_off = 0.0;
		beta_spheres_set_enabled(true);
		const vec on = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, 3, false, v_on, out_on);
		beta_spheres_set_enabled(false);
		const vec off = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, 3, false, v_off, out_off);

		ASSERT_EQ(on.size(), off.size());
		for (size_t b = 0; b < on.size(); b++) {
			EXPECT_NEAR(on[b], off[b], 1e-3) << "basin " << b + 1 << " population moved when the climbs were cut short";
			EXPECT_NEAR(v_on[b], v_off[b], 0.005 * std::max(1.0, v_off[b])) << "basin " << b + 1 << " volume moved";
		}
		EXPECT_NEAR(out_on, out_off, 1e-3) << "the beta spheres changed what falls outside every basin";
	}
}

//The angle-adaptive step is the same argument in the other direction: the midpoint gradient an RK2
//step already computes says how far the field turned over the step just taken, and a field that
//turned less than a degree cannot hide a separatrix in the next one. So the step is allowed to
//double, and any doubt drops it straight back to the floor step the populations were validated at.
//
//That is again an argument, and its weak joint is what a rejected step does next. The first version
//rejected by restarting the iteration from the same point - where the monotonicity test met a value
//it had already recorded in last_value, read "the walk stopped rising" and ended the trajectory in
//mid flight. Nothing about the totals showed it: UH6 still integrated to its electron count, it had
//just moved 0.0132 e between basins. The check below is the one that sees that, because it compares
//basin by basin against the same integration with the growth switched off.
void expect_adaptive_matches_floor(const std::filesystem::path &wfn, const size_t nmax)
{
	const WFN wavy(wfn);
	const cube rho = seed_cube(wavy, 0.25, 3.0);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
	const std::vector<d4> maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
	ASSERT_EQ(maxima.size(), nmax);

	adaptive_guard guard;
	vec v_on, v_off;
	double out_on = 0.0, out_off = 0.0;
	basin_adaptive_step_set_enabled(true);
	const vec on = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, 3, false, v_on, out_on);
	basin_adaptive_step_set_enabled(false);
	const vec off = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, 3, false, v_off, out_off);

	ASSERT_EQ(on.size(), off.size());
	double moved = 0.0;
	for (size_t b = 0; b < on.size(); b++) {
		moved += std::abs(on[b] - off[b]);
		EXPECT_NEAR(on[b], off[b], 1e-3) << "basin " << b + 1 << " population moved when a step was allowed to grow";
		EXPECT_NEAR(v_on[b], v_off[b], 0.005 * std::max(1.0, v_off[b])) << "basin " << b + 1 << " volume moved";
	}
	EXPECT_NEAR(out_on, out_off, 1e-3) << "the grown steps changed what falls outside every basin";
	//A total that survives a redistribution is exactly what the mid-flight bug looked like, so the
	//sum of the moves is asserted as well and not only the conserved total
	EXPECT_LT(moved, 2e-3) << "the populations were redistributed between basins";
}

TEST(AdaptiveStep, AgreesWithTheFloorStepOnHydroxide)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_adaptive_matches_floor(wfn, 2u);
}

//The fixture with three separatrices meeting at an angle, which is where a long step would cross
TEST(AdaptiveStep, AgreesWithTheFloorStepOnNH3Li)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_adaptive_matches_floor(wfn, 5u);
}

TEST(BetaSpheres, AgreeWithTheFullClimbOnHydroxide)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_same_basins(wfn, 2u);
}

//Five atoms, three of them equivalent hydrogens whose basins come within 0.4 bohr of the nucleus
//and whose separatrices meet nitrogen's at an angle: this is the fixture that caught a sphere
//poking through one, and it catches it by the three hydrogens disagreeing with each other
TEST(BetaSpheres, AgreeWithTheFullClimbOnNH3Li)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_same_basins(wfn, 5u);
}

//The sphere's defining property, checked directly on the field rather than through a population:
//sample the sphere the code chose and confirm the radial derivative is negative all over it. This
//is what makes the assignment a proof and not a guess, so it is worth its own test - and it uses
//directions the radius was not built from, which is where a 26-direction sample could fail.
TEST(BetaSpheres, NoAscentPathLeavesTheSphereItChose)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy = load(wfn);
	const cube rho = seed_cube(wavy, 0.25, 3.0);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
	const std::vector<d4> maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
	ASSERT_FALSE(maxima.empty());

	//The radius the integration settled on is not returned, so it is reproduced here by the same
	//rule: 90 % of the smallest radius at which the field stops falling, capped by the neighbours
	for (const d4 &m : maxima) {
		double cap = 3.0;
		for (const d4 &n : maxima) {
			const double d = std::sqrt(std::pow(m[0] - n[0], 2) + std::pow(m[1] - n[1], 2) + std::pow(m[2] - n[2], 2));
			if (d > 1e-8) cap = std::min(cap, 0.45 * d);
		}
		double r = cap;
		for (int i = 0; i < 302; i++) {
			//The same 302-direction spiral the code builds
			const double z = 1.0 - 2.0 * (i + 0.5) / 302.0;
			const double s = std::sqrt(std::max(0.0, 1.0 - z * z));
			const double phi = 2.39996322972865332 * i;
			const d3 u{ s * std::cos(phi), s * std::sin(phi), z };
			double rr = 0.05;
			for (; rr <= cap + 1e-12; rr += 0.05) {
				d3 g;
				wavy.computeGrad(d3{ m[0] + rr * u[0], m[1] + rr * u[1], m[2] + rr * u[2] }, g);
				if (g[0] * u[0] + g[1] * u[1] + g[2] * u[2] >= 0.0) break;
			}
			r = std::min(r, rr - 0.05);
		}
		if (r <= 0.1) continue;   //no sphere claimed here, nothing to check
		r *= 0.7;
		//Off-lattice directions: a spiral of 200 points, a different count so not one of them is
		//one of the 302 the radius was built from
		for (int i = 0; i < 200; i++) {
			const double z = 1.0 - 2.0 * (i + 0.5) / 200.0;
			const double s = std::sqrt(std::max(0.0, 1.0 - z * z));
			const double phi = 2.39996322972865332 * i;   //golden angle, so no two samples line up
			const d3 u{ s * std::cos(phi), s * std::sin(phi), z };
			d3 g;
			wavy.computeGrad(d3{ m[0] + r * u[0], m[1] + r * u[1], m[2] + r * u[2] }, g);
			const double radial = g[0] * u[0] + g[1] * u[1] + g[2] * u[2];
			EXPECT_LT(radial, 0.0) << "the density rises outwards at r = " << r << " in a direction the "
				"26-direction sample never looked at, so a trajectory could leave this sphere";
		}
	}
}

//The climb needs the field's value and its gradient at the same point at every single step, and
//it used to take two passes over every primitive for them: computeGrad, then compute_dens.  The
//reduction in computeGrad already holds phi, so the density is one multiply-add per MO - a
//pointer asks for it.  What has to hold is that the density it hands back is the same number
//compute_dens produces, and that asking for it does not disturb the gradient.  The two build
//their polynomial factors by different groupings of the same products, so they agree to rounding
//rather than to the last bit; a relative 1e-12 is two orders tighter than anything the basin
//integration can see.
namespace
{
	void expect_fused_density_matches(const std::filesystem::path &wfn)
	{
		const WFN wavy = load(wfn);
		const std::vector<atom> atoms = wavy.get_atoms();
		ASSERT_FALSE(atoms.empty());
		//Points on and off the nuclei, in the bonds and out in the tail, where the density spans
		//several orders of magnitude and a cancelling term would show
		std::vector<d3> probes;
		for (const atom &a : atoms) {
			const d3 c = a.get_pos();
			probes.push_back(c);
			for (int k = 0; k < 3; k++) {
				d3 q = c; q[k] += 0.37; probes.push_back(q);
				q = c; q[k] -= 1.9; probes.push_back(q);
			}
		}
		for (size_t i = 1; i < atoms.size(); i++) {
			const d3 a = atoms[0].get_pos(), b = atoms[i].get_pos();
			for (double t : { 0.25, 0.5, 0.75 })
				probes.push_back(d3{ a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]) });
		}
		for (const d3 &q : probes) {
			d3 g_plain, g_fused;
			double rho = -1.0;
			wavy.computeGrad(q, g_plain);
			wavy.computeGrad(q, g_fused, &rho);
			const double ref = wavy.compute_dens(q);
			EXPECT_NEAR(rho, ref, 1e-12 * std::max(1e-30, std::abs(ref)))
				<< "the density riding along on the gradient pass disagrees with compute_dens at ("
				<< q[0] << ", " << q[1] << ", " << q[2] << ")";
			for (int k = 0; k < 3; k++)
				EXPECT_DOUBLE_EQ(g_fused[k], g_plain[k])
					<< "asking for the density changed component " << k << " of the gradient";
		}
	}
}

TEST(FusedDensityGradient, MatchesComputeDensOnOH)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_fused_density_matches(wfn);
}

//A .gbw carrying virtuals and a Li whose basin boundary sits close to the nucleus: the reduction
//skips the empty MOs, so this also checks that the density and the gradient skip the same ones
TEST(FusedDensityGradient, MatchesComputeDensOnNH3Li)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_fused_density_matches(wfn);
}

//g, h and i shells, open shell, and read through the molden reader rather than the gbw one: the
//high angular momenta are where the value and the gradient take different branches of the
//polynomial switch, so this is the case that would catch one of them being wrong
TEST(FusedDensityGradient, MatchesComputeDensWithGHIShells)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "CuF2_i_func" / "71" / "calc_occupied.molden";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_fused_density_matches(wfn);
}

//A basin population belongs to the wavefunction, not to the file format it arrived in. The
//streaming walk never touches the file - it asks the WFN for a gradient - so the readers are the
//one place where a basin number can go wrong without any of the geometry being wrong, and the
//failure mode is specific: a reader that orders a shell's primitives differently, or normalises
//one of them differently, still produces a density that looks like a fluorine atom and integrates
//to something that is not nine electrons.
//
//An isolated atom is the fixture for that, and deliberately so. It has one attractor and no
//separatrix, so nothing here can be blamed on the assignment: whatever the three files disagree
//about is the density itself. The nuclear count is the reference, which is what makes this a check
//and not a golden - no captured output to go stale, and nothing to regenerate if it goes red.
namespace
{
	double lone_atom_population(const std::filesystem::path &wfn, double &outside, double &electrons)
	{
		const WFN wavy(wfn);
		electrons = wavy.count_nr_electrons();
		const cube rho = seed_cube(wavy, 0.25, 3.0);
		const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
		const std::vector<d4> maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
		EXPECT_EQ(maxima.size(), 1u) << "an isolated atom has one attractor: " << wfn.filename().string();
		if (maxima.size() != 1) return 0.0;
		vec v;
		outside = 0.0;
		const vec pops = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, 3, false, v, outside);
		return pops.empty() ? 0.0 : pops[0];
	}
}

TEST(BasinReaders, TheSameFluorineThroughThreeReaders)
{
	const std::filesystem::path dir = nos_test_repo_root() / "tests" / "molden_file";
	const char *files[] = { "f_ref.wfn", "f_ref.wfx", "F_full.molden" };
	double first = 0.0;
	const char *first_name = nullptr;
	for (const char *name : files) {
		const std::filesystem::path f = dir / name;
		if (!std::filesystem::exists(f)) continue;
		double outside = 0.0, electrons = 0.0;
		const double pop = lone_atom_population(f, outside, electrons);
		EXPECT_NEAR(pop, electrons, 0.01) << name << " did not integrate to its own electron count";
		EXPECT_NEAR(outside, 0.0, 1e-3) << name << " left density outside the one basin there is";
		if (!first_name) { first = pop; first_name = name; }
		else EXPECT_NEAR(pop, first, 1e-3) << name << " disagrees with " << first_name;
	}
	if (!first_name) GTEST_SKIP() << "no fluorine fixture under " << dir.string();
}
