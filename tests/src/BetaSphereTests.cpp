#include "pch.h"
#include "core/tuning.h"
#include "core/b2c.h"
#include "core/wfn_class.h"
#include "core/atoms.h"
#include "core/constants.h"
#include "core/cube.h"
#include "core/fiber_loop.h"
#include "core/properties.h"
#ifdef NOSPHERA2_USE_GPU
#include "core/aux_density_gpu.h"
#endif

#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <sstream>
#include <string>

//Inside a beta sphere grad f . rhat < 0 everywhere on its surface, so no ascent path leaves it and a
//trajectory that enters is assigned without climbing. The sphere is sampled along finitely many
//directions and can miss a separatrix, so the populations are compared against the full climb
//(-no_beta_spheres). NH3Li is needed: a two-atom fixture has one planar separatrix that a coarse
//sample still hits. Same-binary A/B differs only by bisection side, far inside the 1e-3 e window.
namespace
{
	WFN load(const std::filesystem::path &p) { return WFN(p); }

	//seeds for the critical-point search only; the streaming path takes no boundaries from it
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

	//restores the switch even when the test fails, so it cannot stay off for the rest of the suite
	struct beta_guard {
		bool was = beta_spheres_enabled();
		~beta_guard() { beta_spheres_set_enabled(was); }
	};

	struct adaptive_guard {
		bool was = basin_adaptive_step_enabled();
		~adaptive_guard() { basin_adaptive_step_set_enabled(was); }
	};

	//one binary, spheres on and off, compared basin by basin
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

//The midpoint gradient of the RK2 step bounds how far the field turned, so a step may double and any
//doubt drops it to the validated floor step. A rejected step that restarts from the same point ends the
//trajectory on the monotonicity test; totals hide that, only a basin-by-basin comparison sees it.
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
	//a redistribution conserves the total, so the sum of the moves is asserted as well
	EXPECT_LT(moved, 2e-3) << "the populations were redistributed between basins";
}

TEST(AdaptiveStep, AgreesWithTheFloorStepOnHydroxide)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_adaptive_matches_floor(wfn, 2u);
}

//A fallback is the fate of a proposal and cannot outnumber proposals; misreading mult, which is raised
//at the end of a step, reverts floor steps and moves no population. Only ELI-D stresses it: its broad
//flat maxima (stall_reach) stop trajectories rising far from an attractor; QTAIM never falls back.
TEST(AdaptiveStep, NeverFallsBackMoreOftenThanItProposesOnELID)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "ELI_heavy" / "uh6.gbw";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy(wfn);
	const cube rho = seed_cube(wavy, 0.3, 2.5);

	//ELI-D's own maxima, since the walk climbs ELI-D
	cube eli(rho.get_sizes(), 0, true);
	eli.set_origin(0, rho.get_origin(0)); eli.set_origin(1, rho.get_origin(1)); eli.set_origin(2, rho.get_origin(2));
	for (int k = 0; k < 3; k++) eli.set_vector(k, k, rho.get_vector(k, k));
	eli.calc_dv();
	std::ostringstream log;
	Calc_Eli(eli, wavy, 3.0, log, false);
	const std::vector<d4> maxima = topological_cube_analysis(&eli, wavy.get_atoms(), false, false, 0.0, 1e-12, -1.0, 5e-3, nullptr, &wavy).second;
	ASSERT_GT(maxima.size(), 3u) << "no ELI-D attractors, so nothing is climbed";

	adaptive_guard guard;
	basin_adaptive_step_set_enabled(true);
	basin_adaptive_step_counters_reset();
	vec volumes;
	double outside = 0.0;
	integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, 3, true, volumes, outside);

	long long steps = 0, proposed = 0, turned = 0, fell = 0;
	basin_adaptive_step_counters(steps, proposed, turned, fell);
	basin_adaptive_step_counters_reset();
	ASSERT_GT(steps, 0) << "the growth was enabled and no step was counted, so this asserts nothing";
	ASSERT_GT(proposed, 0) << "no longer step was ever proposed, so the counters cannot be compared";
	EXPECT_LE(fell, proposed)
		<< fell << " steps were reverted as grown against " << proposed
		<< " grown proposals, so the fallback is blaming steps that ran at the floor";
	EXPECT_LE(turned + fell, proposed) << "more proposals were rejected than were ever made";
}

//sets a -tune knob for the scope and restores it, set or unset; a leaked knob changes every later test
struct env_guard {
	std::string name;
	std::string old;
	bool had = false;
	env_guard(const char *n, const char *v) : name(n)
	{
		if (const char *e = tuning(n)) { old = e; had = true; }
		set_tuning(name, v);
	}
	~env_guard() { set_tuning(name, had ? old.c_str() : nullptr); }
};

//checks the NOS_ADP_* overrides reach the walk
TEST(AdaptiveStep, KnobsComeFromTheEnvironment)
{
	//declared first, destroyed last: the env_guards have cleared the variables, so re-reading restores the defaults
	struct knob_restore {
		bool was = basin_adaptive_step_enabled();
		~knob_restore() { basin_adaptive_step_set_enabled(true); basin_adaptive_step_set_enabled(was); }
	} restore;

	double cap = 0.0, grow = 0.0, keep = 0.0, reach = 0.0;
	basin_adaptive_step_set_enabled(true);
	basin_adaptive_step_knobs(cap, grow, keep, reach);
	const double shipped_grow = grow, shipped_cap = cap;
	EXPECT_GT(shipped_grow, 0.9) << "the shipped cosine gate is a tight one; this test assumes it";

	{
		env_guard g("NOS_ADP_GROW", "0.99");
		basin_adaptive_step_set_enabled(true);
		basin_adaptive_step_knobs(cap, grow, keep, reach);
		EXPECT_DOUBLE_EQ(grow, 0.99) << "NOS_ADP_GROW was set and the walk would still use " << grow;
		EXPECT_DOUBLE_EQ(cap, shipped_cap) << "setting one knob moved another";
	}
	basin_adaptive_step_set_enabled(true);
	basin_adaptive_step_knobs(cap, grow, keep, reach);
	EXPECT_DOUBLE_EQ(grow, shipped_grow) << "clearing the variable did not restore the validated default";

	//Junk must not be parsed into a zero: a cosine gate of 0 would grow every step in the suite
	{
		env_guard g("NOS_ADP_GROW", "not-a-number");
		basin_adaptive_step_set_enabled(true);
		basin_adaptive_step_knobs(cap, grow, keep, reach);
		EXPECT_DOUBLE_EQ(grow, shipped_grow) << "unparseable knob was not ignored";
	}
	{
		env_guard g("NOS_ADP_CAP", "-4");
		basin_adaptive_step_set_enabled(true);
		basin_adaptive_step_knobs(cap, grow, keep, reach);
		EXPECT_DOUBLE_EQ(cap, shipped_cap) << "a negative cap was accepted";
	}
}

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

//three equivalent hydrogens whose separatrices meet nitrogen's at an angle: a sphere poking through one
//shows as the hydrogens disagreeing
TEST(BetaSpheres, AgreeWithTheFullClimbOnNH3Li)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_same_basins(wfn, 5u);
}

//radial derivative negative all over the chosen sphere, sampled in directions the radius was not built from
TEST(BetaSpheres, NoAscentPathLeavesTheSphereItChose)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy = load(wfn);
	const cube rho = seed_cube(wavy, 0.25, 3.0);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
	const std::vector<d4> maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
	ASSERT_FALSE(maxima.empty());

	//the radius is not returned, so it is rebuilt by the same rule: 90 % of the smallest radius at which the
	//field stops falling, capped by the neighbours
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
		r *= basin_beta_margin();   //whatever margin the integration is drawing them at
		//200-point spiral: none of its directions is one of the 302
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

//The margin must reach the code from NOS_BETA_MARGIN, and only values in (0, 1] are accepted: above 1
//the sphere is wider than the radius measured safe, which is a wrong population, not a looser setting.
TEST(BetaSpheres, MarginComesFromTheEnvironmentAndStaysAMargin)
{
	const double shipped = basin_beta_margin();
	EXPECT_GT(shipped, 0.0);
	EXPECT_LE(shipped, 1.0) << "the shipped margin is already outside the range this test enforces";

	{
		env_guard g("NOS_BETA_MARGIN", "0.85");
		EXPECT_DOUBLE_EQ(basin_beta_margin(), 0.85) << "NOS_BETA_MARGIN was set and the spheres would "
			"still be drawn at " << basin_beta_margin();
	}
	EXPECT_DOUBLE_EQ(basin_beta_margin(), shipped) << "clearing the variable did not restore the default";

	for (const char *bad : { "1.4", "0", "-0.8", "not-a-number", "" }) {
		env_guard g("NOS_BETA_MARGIN", bad);
		EXPECT_DOUBLE_EQ(basin_beta_margin(), shipped) << "NOS_BETA_MARGIN=" << bad << " was accepted";
	}
}

//The computeGrad reduction already holds phi, so the density is one multiply-add per MO. It must equal
//compute_dens to rounding (different product groupings) and must leave the gradient untouched.
namespace
{
	void expect_fused_density_matches(const std::filesystem::path &wfn)
	{
		const WFN wavy = load(wfn);
		const std::vector<atom> atoms = wavy.get_atoms();
		ASSERT_FALSE(atoms.empty());
		//on and off the nuclei, in the bonds and in the tail, where a cancelling term would show
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

//virtuals in a .gbw: the density and the gradient must skip the same empty MOs
TEST(FusedDensityGradient, MatchesComputeDensOnNH3Li)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_fused_density_matches(wfn);
}

//g, h and i shells, open shell, molden reader: high l is where value and gradient take different
//branches of the polynomial switch
TEST(FusedDensityGradient, MatchesComputeDensWithGHIShells)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "CuF2_i_func" / "71" / "calc_occupied.molden";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	expect_fused_density_matches(wfn);
}

//A reader that orders or normalises a shell's primitives differently still gives a plausible density
//with the wrong integral. A lone atom has one attractor and no separatrix, so any disagreement is the
//density itself; the nuclear count is the reference.
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
		EXPECT_GT(outside, 0.0) << name << " has no density beyond the isosurface";
		EXPECT_LT(outside, 0.01) << name << " lost density inside the isosurface";
		EXPECT_NEAR(pop + outside, electrons, 0.01) << name << " did not conserve the quadrature";
		if (!first_name) { first = pop; first_name = name; }
		else EXPECT_NEAR(pop, first, 1e-3) << name << " disagrees with " << first_name;
	}
	if (!first_name) GTEST_SKIP() << "no fluorine fixture under " << dir.string();
}

//open shell: a per-MO occupancy taken from the wrong place integrates to another spin state's count
TEST(BasinReaders, TheOpenShellFluorineAtThreeSpinStates)
{
	const std::filesystem::path dir = nos_test_repo_root() / "tests" / "molden_file";
	const char *files[] = { "F_open.molden", "F_s1.molden", "F_s32.molden" };
	int seen = 0;
	for (const char *name : files) {
		const std::filesystem::path f = dir / name;
		if (!std::filesystem::exists(f)) continue;
		double outside = 0.0, electrons = 0.0;
		const double pop = lone_atom_population(f, outside, electrons);
		EXPECT_GT(electrons, 0.0) << name << " reported no electrons at all";
		EXPECT_NEAR(pop, electrons, 0.01) << name << " integrated to " << pop << " and not to its own "
			<< electrons << " electrons";
		EXPECT_GT(outside, 0.0) << name << " has no density beyond the isosurface";
		EXPECT_LT(outside, 0.01) << name << " lost density inside the isosurface";
		EXPECT_NEAR(pop + outside, electrons, 0.01) << name << " did not conserve the quadrature";
		seen++;
	}
	if (!seen) GTEST_SKIP() << "no open-shell fluorine fixture under " << dir.string();
}

//A lone atom has no second-nearest atom, so the ELI-D labelling must not demand one, and its non-core
//basin is its valence shell, never a bond.
TEST(BasinLabels, ALoneAtomHasNoBondBasin)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "molden_file" / "F_full.molden";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy(wfn);
	ASSERT_EQ(wavy.get_ncen(), 1) << "this fixture is meant to be one atom";
	const std::vector<atom> atoms = wavy.get_atoms();
	const double x = wavy.get_atom_coordinate(0, 0);
	const double y = wavy.get_atom_coordinate(0, 1);
	const double z = wavy.get_atom_coordinate(0, 2);
	//one maximum on the nucleus (core) and one 1.4 bohr out along z (valence shell)
	const std::vector<d4> maxima = { { x, y, z, 1.0 }, { x, y, z + 1.4, 0.1 } };

	const svec eli = assign_labels_to_basins(maxima, atoms, false, 1);
	ASSERT_EQ(eli.size(), maxima.size());
	for (size_t i = 0; i < eli.size(); i++) {
		EXPECT_NE(eli[i].find(atoms[0].get_label()), std::string::npos)
			<< "ELI label " << i << " (\"" << eli[i] << "\") does not name the only atom there is";
		EXPECT_EQ(eli[i].find("bond"), std::string::npos)
			<< "ELI label " << i << " (\"" << eli[i] << "\") calls a lone atom's basin a bond";
	}
	EXPECT_NE(eli[0].find("core"), std::string::npos) << "the maximum on the nucleus is the core basin: \"" << eli[0] << "\"";

	//the QTAIM branch, so a fix cannot trade one branch for the other
	const svec qtaim = assign_labels_to_basins(maxima, atoms, false, 0);
	ASSERT_EQ(qtaim.size(), maxima.size());
	EXPECT_NE(qtaim[0].find(atoms[0].get_label()), std::string::npos) << "QTAIM label 0: \"" << qtaim[0] << "\"";
	EXPECT_NE(qtaim[1].find("NNA"), std::string::npos)
		<< "a maximum 1.4 bohr off the only nucleus is a non-nuclear attractor: \"" << qtaim[1] << "\"";
}

//constants::density_accuracy is load-bearing, not slack: ELI-D g = rho tau - |grad rho|^2/4 is a
//difference of large terms and amplifies a truncation invisible in rho. Only the knob and its range
//are checked.
TEST(ExpCutoff, ScreeningAccuracyIsSettableAndRangeChecked)
{
	const WFN wavy = load(nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn");
	wavy.set_exp_cutoff();
	const double shipped = constants::exp_cutoff;
	ASSERT_LT(shipped, 0.0);
	{
		const env_guard g("NOS_DENSITY_ACCURACY", "1e-2");
		wavy.set_exp_cutoff();
		//a looser accuracy is a cutoff nearer zero: fewer primitives survive it
		EXPECT_GT(constants::exp_cutoff, shipped);
	}
	//out of range or unreadable keeps the shipped accuracy: >= 1 is not an accuracy, <= 0 breaks the log
	for (const char *bad : { "1.0", "2", "0", "-1e-3", "not-a-number", "" })
	{
		const env_guard g("NOS_DENSITY_ACCURACY", bad);
		wavy.set_exp_cutoff();
		EXPECT_DOUBLE_EQ(constants::exp_cutoff, shipped) << "accepted " << bad;
	}
	//the guard cleared it, so no later test screens differently
	wavy.set_exp_cutoff();
	EXPECT_DOUBLE_EQ(constants::exp_cutoff, shipped);
}

//The basin boundary error is set by the angular order of the atomic grid, worst on ionic centres; the
//delocalization residual cannot see it, since sum_A S^A is the identity for any partition. Refining
//must keep converging: |N(3) - N(5)| > |N(4) - N(5)|. config.accuracy is max(accuracy, 3), so level 2
//does not exist here.
TEST(QuadratureAccuracy, RefiningTheGridKeepsConverging)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy(wfn);
	const cube rho = seed_cube(wavy, 0.25, 3.0);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
	const std::vector<d4> maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
	ASSERT_EQ(maxima.size(), 2u);

	//by nucleus, not by index
	size_t ox = maxima.size();
	for (int a = 0; a < wavy.get_ncen(); a++) {
		if (wavy.get_atom_charge(a) != 8) continue;
		for (size_t b = 0; b < maxima.size(); b++)
			if (std::pow(maxima[b][0] - wavy.get_atom_coordinate(a, 0), 2) +
				std::pow(maxima[b][1] - wavy.get_atom_coordinate(a, 1), 2) +
				std::pow(maxima[b][2] - wavy.get_atom_coordinate(a, 2), 2) < 0.01) ox = b;
	}
	ASSERT_LT(ox, maxima.size()) << "no basin sits on the oxygen";

	vec n[3], v[3];
	double out[3] = { 0.0, 0.0, 0.0 };
	for (int lvl = 3; lvl <= 5; lvl++)
		n[lvl - 3] = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, lvl, false, v[lvl - 3], out[lvl - 3]);

	const double e3 = std::abs(n[0][ox] - n[2][ox]);
	const double e4 = std::abs(n[1][ox] - n[2][ox]);
	std::cout << "  oxygen basin: level 3 " << n[0][ox] << ", level 4 " << n[1][ox]
		<< ", level 5 " << n[2][ox] << "  (|3-5| " << e3 << ", |4-5| " << e4 << ")" << std::endl;
	//strict: two levels that agree exactly mean one was not the level asked for
	EXPECT_GT(e3, e4) << "refining the grid stopped converging on the oxygen basin";

	//Electrons from the occupations, not the nuclear charges: OH.wfn is the anion. The anion's diffuse tail
	//reaches past the atomic grids' radial range at every level, hence the 0.02 e window.
	double total[3] = { 0.0, 0.0, 0.0 };
	for (int lvl = 0; lvl < 3; lvl++) {
		total[lvl] = out[lvl];
		for (const double b : n[lvl]) total[lvl] += b;
		EXPECT_NEAR(total[lvl], wavy.count_nr_electrons(), 0.02) << "level " << lvl + 3;
	}
	//the tail is set by the grid's extent, so refining must not move the total
	EXPECT_NEAR(total[0], total[2], 1e-3) << "the integrated total moved with the accuracy level";
	EXPECT_NEAR(total[1], total[2], 1e-3) << "the integrated total moved with the accuracy level";
}

//A stalled trajectory must still be assigned, at any distance: a homopolar bond critical point lies
//more than a bohr from either nucleus. Removing Co2's non-nuclear attractor makes its charge stall
//between the cobalts. The stall goes by the trajectory's start, not where it stopped: at a saddle
//between equivalent atoms the nearest attractor is a floating-point coin flip.
TEST(BasinStalls, AStalledTrajectoryIsStillAssigned)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "molden_file" / "Co2.molden";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	const WFN wavy(wfn);
	const cube rho = seed_cube(wavy, 0.1, 3.0);
	const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
	const std::vector<d4> all = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);

	//nuclear attractors by nucleus; the count is not asserted because it depends on the seed resolution
	std::vector<d4> nuclei;
	for (const d4 &m : all) {
		for (int a = 0; a < wavy.get_ncen(); a++)
			if (std::pow(m[0] - wavy.get_atom_coordinate(a, 0), 2) +
				std::pow(m[1] - wavy.get_atom_coordinate(a, 1), 2) +
				std::pow(m[2] - wavy.get_atom_coordinate(a, 2), 2) < 0.01) { nuclei.push_back(m); break; }
	}
	for (const d4 &m : all)
		std::cout << "  attractor at " << m[0] << " " << m[1] << " " << m[2] << " rho " << m[3] << std::endl;
	ASSERT_EQ(nuclei.size(), static_cast<size_t>(wavy.get_ncen())) << "an attractor is missing from a nucleus";
	ASSERT_GT(all.size(), nuclei.size()) << "Co2 no longer resolves any non-nuclear attractor";

	vec n_all, n_nuc, v_all, v_nuc;
	double out_all = 0.0, out_nuc = 0.0;
	n_all = integrate_basins_on_atomic_grids(nullptr, nullptr, all, wavy, 3, false, v_all, out_all);
	n_nuc = integrate_basins_on_atomic_grids(nullptr, nullptr, nuclei, wavy, 3, false, v_nuc, out_nuc);
	ASSERT_EQ(n_all.size(), all.size());
	ASSERT_EQ(n_nuc.size(), nuclei.size());

	//A stall with a gradient still on it is a floor step too long for the local curvature: the step must
	//halve below its floor. The count is asserted, as populations move by less than any tolerance here.
	EXPECT_EQ(basin_stalls_on_a_slope(), 0)
		<< "a density trajectory stopped rising with a gradient still on it: that is a step too long for "
		   "the local curvature, not a critical point, and the step must halve rather than give up";

	double total_all = out_all, total_nuc = out_nuc, non_nuclear = 0.0;
	for (const double b : n_all) total_all += b;
	for (const double b : n_nuc) total_nuc += b;
	for (size_t b = nuclei.size(); b < n_all.size(); b++) non_nuclear += n_all[b];
	std::cout << "  with them:  " << n_all[0] << " + " << n_all[1] << ", non-nuclear " << non_nuclear
		<< ", outside " << out_all << "\n  without:    " << n_nuc[0] << " + " << n_nuc[1]
		<< ", outside " << out_nuc << ", total " << total_nuc << " vs " << total_all << std::endl;
	//the dropped attractors have to hold something, or removing them proves nothing
	ASSERT_GT(non_nuclear, 1.0) << "the non-nuclear attractors hold nothing to redistribute";

	//the NNA's electrons stay in the molecule wherever the partition puts them
	EXPECT_NEAR(total_nuc, total_all, 5e-3) << "the NNA's electrons went missing when it was not an attractor";
	EXPECT_NEAR(out_nuc, out_all, 5e-3) << "removing a non-nuclear attractor lost density inside the isosurface";
	//and land in the cobalts, roughly half of 4.15 e each
	EXPECT_GT(n_nuc[0], n_all[0] + 0.5) << "the first cobalt did not gain the non-nuclear charge";
	EXPECT_GT(n_nuc[1], n_all[1] + 0.5) << "the second cobalt did not gain the non-nuclear charge";
	//a homonuclear diatomic cannot prefer one end; a stall-point tie-break fails this
	EXPECT_NEAR(n_nuc[0], n_nuc[1], 1e-3) << "the two cobalts split the stalled charge unevenly";
}

//-basin_gpu: the same streaming integration with the field from the device. The fields agree to 1e-10
//(BasinFieldGpuTests), so only a probe starting on a separatrix can change basin; the window is 1e-4 e per basin.
static void expect_gpu_matches_host(const std::filesystem::path &wfn, const bool eli_field)
{
#ifdef NOSPHERA2_USE_GPU
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	if (!aux_density_gpu_available()) GTEST_SKIP() << "no device";
	struct guard { ~guard() { basin_gpu_set_mode(-1); aux_density_gpu_set_enabled(false); } } restore;
	aux_density_gpu_set_enabled(true);
	const WFN wavy(wfn);
	{
		//the device has to take the field at all, or the second pass is the first one again
		const double p[3]{ 0.1, 0.2, 0.3 };
		double g[3];
		ASSERT_TRUE(wavy.field_grad_gpu(false, 1, p, nullptr, g)) << "the device declined the field";
	}
	const cube rho = seed_cube(wavy, 0.25, 3.0);
	std::vector<d4> maxima;
	if (!eli_field) {
		const std::vector<critical_point> cps = analyze_cube_critical_points(&rho, wavy, false, std::max(1e-8, rho.max_value() * 1e-6));
		maxima = streaming_density_attractors(wavy, cps, nullptr, nullptr, false);
	}
	else {
		cube eli(rho.get_sizes(), 0, true);
		for (int k = 0; k < 3; k++) { eli.set_origin(k, rho.get_origin(k)); eli.set_vector(k, k, rho.get_vector(k, k)); }
		eli.calc_dv();
		std::ostringstream log;
		Calc_Eli(eli, wavy, 3.0, log, false);
		maxima = topological_cube_analysis(&eli, wavy.get_atoms(), false, false, 0.0, 1e-12, -1.0, 5e-3, nullptr, &wavy).second;
	}
	ASSERT_FALSE(maxima.empty());
	vec v_host, v_dev;
	double o_host = 0.0, o_dev = 0.0;
	//QTAIM also takes the overlap matrices, whose orbitals the fibers share per thread
	//(not s_host: winsock defines that)
	basin_overlaps ov_host, ov_dev;
	basin_gpu_set_mode(0);
	const vec host = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, 3, eli_field, v_host, o_host, nullptr, nullptr, 1, nullptr, eli_field ? nullptr : &ov_host);
	basin_gpu_set_mode(1);
	const vec dev = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wavy, 3, eli_field, v_dev, o_dev, nullptr, nullptr, 1, nullptr, eli_field ? nullptr : &ov_dev);
	ASSERT_EQ(host.size(), dev.size());
	ASSERT_EQ(ov_host.S.size(), ov_dev.S.size());
	for (size_t b = 0; b < ov_host.S.size(); b++)
		for (size_t t = 0; t < ov_host.S[b].size(); t++)
			EXPECT_NEAR(ov_dev.S[b][t], ov_host.S[b][t], 1e-4) << "basin " << b + 1 << " overlap " << t;
	for (size_t b = 0; b < host.size(); b++) {
		std::cout << "  basin " << b + 1 << ": host " << host[b] << ", device " << dev[b] << ", diff " << dev[b] - host[b] << std::endl;
		EXPECT_NEAR(dev[b], host[b], 1e-4) << "basin " << b + 1 << " population";
		EXPECT_NEAR(v_dev[b], v_host[b], 1e-3 * std::max(1.0, v_host[b])) << "basin " << b + 1 << " volume";
	}
	EXPECT_NEAR(o_dev, o_host, 1e-4) << "outside every basin";
#else
	(void)wfn; (void)eli_field;
	GTEST_SKIP() << "built without a GPU backend";
#endif
}

//Every body runs once and keeps its own stack across parks; a round ends once all live fibers have parked
TEST(FiberLoop, BodiesKeepTheirStacksAcrossParks)
{
	const int n = 1000, parks = 5;
	std::atomic<int> next{ 0 };
	int rounds = 0, wrong = 0;
	std::vector<int> ran(n, 0);
	fiber_loop(16, next, n, [&](const int i) {
		double mine[32];
		for (int k = 0; k < 32; k++) mine[k] = i * 32.0 + k;
		for (int y = 0; y < parks; y++) { fiber_yield(); mine[y] += 0.25; }
		for (int k = 0; k < 32; k++) wrong += mine[k] != i * 32.0 + k + (k < parks ? 0.25 : 0.0);
		ran[i]++;
	}, [&]() { rounds++; });
	EXPECT_EQ(wrong, 0);
	EXPECT_EQ(std::count(ran.begin(), ran.end(), 1), n);
	EXPECT_GE(rounds, n / 16 * parks);
}

TEST(BasinGpu, QTAIMMatchesHostOnNH3Li)
{
	expect_gpu_matches_host(nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw", false);
}

TEST(BasinGpu, ELIDMatchesHostOnNH3Li)
{
	expect_gpu_matches_host(nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw", true);
}

//The ELI-D maxima search with its seed climbs on the device finds the host's maxima
TEST(BasinGpu, ELIDMaximaMatchHostOnNH3Li)
{
#ifdef NOSPHERA2_USE_GPU
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	if (!aux_density_gpu_available()) GTEST_SKIP() << "no device";
	struct guard { ~guard() { basin_gpu_set_mode(-1); aux_density_gpu_set_enabled(false); } } restore;
	aux_density_gpu_set_enabled(true);
	const WFN wavy(wfn);
	basin_gpu_set_mode(0);
	const std::vector<d4> host = analytic_eli_maxima(wavy);
	basin_gpu_set_mode(1);
	const std::vector<d4> dev = analytic_eli_maxima(wavy);
	ASSERT_EQ(dev.size(), host.size());
	for (const d4 &h : host) {
		double best = 1e30, f = 0.0;
		for (const d4 &d : dev) {
			const double r2 = std::pow(d[0] - h[0], 2) + std::pow(d[1] - h[1], 2) + std::pow(d[2] - h[2], 2);
			if (r2 < best) { best = r2; f = d[3]; }
		}
		EXPECT_LT(std::sqrt(best), 1e-3) << "host maximum " << h[0] << " " << h[1] << " " << h[2];
		EXPECT_NEAR(f, h[3], 1e-6);
	}
#else
	GTEST_SKIP() << "built without a GPU backend";
#endif
}
