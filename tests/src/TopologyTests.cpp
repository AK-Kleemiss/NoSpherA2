#include "pch.h"

#include "core/topology.h"

//topology.h: ring and cage critical points from a Newton-Raphson search on grad rho = 0 with the
//analytic Hessian, and Poincare-Hopf as the test of whether the set is complete.
//
//Two kinds of reference, and deliberately no golden file: a golden file records what this code
//already does, so it cannot catch a systematic error in it.
//  1. Exactly known matrices, for the signature logic.
//  2. A sum of spherical Gaussians, whose topology is fixed by symmetry rather than by our output.
//     Six equal Gaussians on an octahedron must give 6 attractors, 12 bond, 8 ring and 1 cage
//     point - the octahedron's vertices, edges, faces and centre - so 6 - 12 + 8 - 1 = 1.  The
//     cage point is at the origin, the one point O_h leaves fixed, where the Hessian must be
//     isotropic; the ring points lie on the eight <111> rays and the bond points on the twelve
//     <110> rays.  None of that is read off a previous run.
//     Two equal Gaussians at +-1/sqrt(2a) give an exactly degenerate critical point at the origin:
//     H_xx = 2 (4 a^2 s^2 - 2a) exp(-a s^2) vanishes at s^2 = 1/(2a) while H_yy = H_zz = -4a
//     exp(-a s^2) stays negative, so the signature is (-, -, 0) analytically.
namespace
{
	//A sum of equal spherical Gaussians, rho = sum_i exp(-a |p - R_i|^2), with its analytic gradient
	//and Hessian.  Satisfies density_source.h's PointDensityWithHessian, so the same search that
	//runs on a WFN runs on this.
	struct gaussian_sum
	{
		std::vector<d3> centres;
		double a = 1.0;

		double hessian(const d3& p, d3& grad, double* H) const
		{
			double rho = 0.0;
			grad = { 0.0, 0.0, 0.0 };
			std::fill(H, H + 9, 0.0);
			for (const d3& c : centres) {
				const d3 d{ p[0] - c[0], p[1] - c[1], p[2] - c[2] };
				const double w = std::exp(-a * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
				rho += w;
				for (int i = 0; i < 3; i++) {
					grad[i] += -2.0 * a * d[i] * w;
					for (int j = 0; j < 3; j++)
						H[3 * i + j] += (4.0 * a * a * d[i] * d[j] - (i == j ? 2.0 * a : 0.0)) * w;
				}
			}
			return rho;
		}
	};

	//The six octahedron vertices at distance d along the axes, as nuclei of charge 1
	std::vector<topology::nucleus> octahedron_nuclei(const double d)
	{
		std::vector<topology::nucleus> n;
		for (int k = 0; k < 3; k++)
			for (const int sgn : { +1, -1 }) {
				d3 p{ 0.0, 0.0, 0.0 };
				p[k] = sgn * d;
				n.push_back({ p, 1 });
			}
		return n;
	}

	gaussian_sum octahedron_source(const std::vector<topology::nucleus>& n, const double a)
	{
		gaussian_sum s;
		s.a = a;
		for (const topology::nucleus& x : n) s.centres.push_back(x.pos);
		return s;
	}

	//Hydrogen's covalent radius is 0.23 A, so the default bond_scale of 1.3 calls nothing further
	//apart than 1.1 bohr bonded and the test densities would get no bond seeds at all.  These are
	//not real molecules; widening the seeding criterion is the honest way to say so.  It only adds
	//seeds - which pair of centres is bonded is read back off the bond critical points that the
	//search actually finds, not off this number.
	topology::options gaussian_options()
	{
		topology::options opt;
		opt.bond_scale = 10.0;
		return opt;
	}
}

//---------------------------------------------------------------------------------------------
// The signature logic, against matrices whose eigenvalues are written down rather than computed
//---------------------------------------------------------------------------------------------
TEST(Topology, ClassifyHessianSignatures)
{
	d3 w;
	int neg = 0, pos = 0, zero = 0;
	//(3,-3) attractor: all three curvatures negative
	const double a[9]{ -1.0, 0, 0, 0, -2.0, 0, 0, 0, -3.0 };
	EXPECT_EQ(topology::classify_hessian(a, w, neg, pos, zero), topology::cp_kind::attractor);
	EXPECT_EQ(neg, 3);
	EXPECT_EQ(zero, 0);
	//eigenvalues come back ascending, which is what the ellipticity lambda1/lambda2 - 1 assumes
	EXPECT_NEAR(w[0], -3.0, 1E-12);
	EXPECT_NEAR(w[2], -1.0, 1E-12);
	//(3,-1) bond
	const double b[9]{ -1.0, 0, 0, 0, -2.0, 0, 0, 0, 0.5 };
	EXPECT_EQ(topology::classify_hessian(b, w, neg, pos, zero), topology::cp_kind::bond);
	EXPECT_EQ(neg, 2);
	EXPECT_EQ(pos, 1);
	//(3,+1) ring
	const double c[9]{ -1.0, 0, 0, 0, 2.0, 0, 0, 0, 0.5 };
	EXPECT_EQ(topology::classify_hessian(c, w, neg, pos, zero), topology::cp_kind::ring);
	EXPECT_EQ(neg, 1);
	EXPECT_EQ(pos, 2);
	//(3,+3) cage
	const double e[9]{ 1.0, 0, 0, 0, 2.0, 0, 0, 0, 0.5 };
	EXPECT_EQ(topology::classify_hessian(e, w, neg, pos, zero), topology::cp_kind::cage);
	EXPECT_EQ(pos, 3);
	//A rotation of the ring matrix: same eigenvalues, no longer diagonal, so the eigen solver and
	//not the diagonal is being read
	const double s2 = std::sqrt(0.5);
	//R = rotation by 45 deg about z applied to diag(-1, 2, 0.5): R diag R^T
	const double f[9]{
		0.5 * (-1.0 + 2.0), 0.5 * (-1.0 - 2.0), 0.0,
		0.5 * (-1.0 - 2.0), 0.5 * (-1.0 + 2.0), 0.0,
		0.0, 0.0, 0.5 };
	EXPECT_EQ(topology::classify_hessian(f, w, neg, pos, zero), topology::cp_kind::ring);
	EXPECT_NEAR(w[0], -1.0, 1E-10);
	EXPECT_NEAR(w[1], 0.5, 1E-10);
	EXPECT_NEAR(w[2], 2.0, 1E-10);
	EXPECT_GT(s2, 0.0); //keep the constant referenced
}

//A zero eigenvalue must be reported as degenerate, not rounded into one of the four signatures.
//Binning it would let the Poincare-Hopf sum come out right for the wrong reason.
TEST(Topology, DegenerateHessianIsReportedNotBinned)
{
	d3 w;
	int neg = 0, pos = 0, zero = 0;
	//(-, -, 0): the exact signature of two coalescing maxima. Would be an attractor if the zero
	//were rounded down, a bond point if it were rounded up
	const double h[9]{ -2.0, 0, 0, 0, -2.0, 0, 0, 0, 0.0 };
	EXPECT_EQ(topology::classify_hessian(h, w, neg, pos, zero), topology::cp_kind::degenerate);
	EXPECT_EQ(zero, 1);
	EXPECT_EQ(neg, 2);
	//A curvature below the relative tolerance is zero: 2E-7 against max|lambda| = 2 and the
	//default 1E-6 relative tolerance
	const double h2[9]{ -2.0, 0, 0, 0, -2.0, 0, 0, 0, 2E-7 };
	EXPECT_EQ(topology::classify_hessian(h2, w, neg, pos, zero), topology::cp_kind::degenerate);
	//and one above it is not
	const double h3[9]{ -2.0, 0, 0, 0, -2.0, 0, 0, 0, 2E-5 };
	EXPECT_EQ(topology::classify_hessian(h3, w, neg, pos, zero), topology::cp_kind::bond);
	//A single degenerate point makes the whole analysis incomplete even when the sum adds to 1
	topology::result r;
	topology::cp attractor;
	attractor.kind = topology::cp_kind::attractor;
	topology::cp degenerate;
	degenerate.kind = topology::cp_kind::degenerate;
	r.points = { attractor, degenerate };
	topology::tally(r);
	EXPECT_EQ(r.sum, 1);
	EXPECT_TRUE(r.balanced);
	EXPECT_FALSE(r.complete) << "a degenerate point must void the completeness verdict";
	EXPECT_NE(r.diagnosis.find("degenerate"), std::string::npos) << r.diagnosis;
}

//The degenerate point of a real density, placed analytically: two equal Gaussians separated by
//2/sqrt(2a) have H_xx = 0 at their midpoint while H_yy = H_zz = -4a exp(-a s^2).
TEST(Topology, AnalyticallyDegenerateDensityCriticalPoint)
{
	const double a = 1.0;
	const double s = 1.0 / std::sqrt(2.0 * a);
	gaussian_sum src;
	src.a = a;
	src.centres = { d3{ s, 0.0, 0.0 }, d3{ -s, 0.0, 0.0 } };
	const std::vector<topology::nucleus> nuc{ { { s, 0.0, 0.0 }, 1 }, { { -s, 0.0, 0.0 }, 1 } };

	const topology::cp point = topology::describe_cp(src, d3{ 0.0, 0.0, 0.0 }, nuc, gaussian_options());
	EXPECT_NEAR(point.gradient_norm, 0.0, 1E-14) << "the midpoint of two equal Gaussians is stationary by symmetry";
	//the two transverse curvatures, in closed form
	const double transverse = -4.0 * a * std::exp(-a * s * s);
	EXPECT_NEAR(point.eigenvalues[0], transverse, 1E-10);
	EXPECT_NEAR(point.eigenvalues[1], transverse, 1E-10);
	EXPECT_NEAR(point.eigenvalues[2], 0.0, 1E-12) << "the longitudinal curvature vanishes at s = 1/sqrt(2a)";
	EXPECT_EQ(point.kind, topology::cp_kind::degenerate);
	EXPECT_EQ(point.zero, 1);

	//and the full search on that density must refuse to certify it
	const topology::result r = topology::analyze_topology(src, nuc, gaussian_options());
	EXPECT_GT(r.n_degenerate, 0);
	EXPECT_FALSE(r.complete);
	EXPECT_NE(r.diagnosis.find("degenerate"), std::string::npos) << r.diagnosis;
}

//---------------------------------------------------------------------------------------------
// The cage branch, on a density whose topology symmetry fixes
//---------------------------------------------------------------------------------------------
TEST(Topology, OctahedronHasEightRingsAndOneCage)
{
	const double d = 2.0, a = 1.0;
	const std::vector<topology::nucleus> nuc = octahedron_nuclei(d);
	const gaussian_sum src = octahedron_source(nuc, a);
	const topology::result r = topology::analyze_topology(src, nuc, gaussian_options());

	EXPECT_EQ(r.n_attractor, 6) << "one maximum per vertex";
	EXPECT_EQ(r.n_bond, 12) << "one bond point per octahedron edge";
	EXPECT_EQ(r.n_ring, 8) << "one ring point per triangular face";
	EXPECT_EQ(r.n_cage, 1) << "one cage point, at the centre";
	EXPECT_EQ(r.n_degenerate, 0);
	EXPECT_EQ(r.sum, 1) << "6 - 12 + 8 - 1";
	EXPECT_TRUE(r.complete) << r.diagnosis;
	EXPECT_FALSE(r.escalated) << "the topology-driven seeding should find all of these on its own";
	EXPECT_EQ(r.n_nna, 0) << "every maximum of this density sits on a centre";
	//the bond graph's cycle rank, 12 - 6 + 1 = 7, equals n_ring - n_cage = 8 - 1
	EXPECT_EQ(r.required_ring_minus_cage, 7);
	EXPECT_EQ(r.found_ring_minus_cage, 7);

	//O_h leaves only the origin fixed, so the cage point is there and its Hessian is isotropic
	const topology::cp* cage = nullptr;
	for (const topology::cp& p : r.points) if (p.kind == topology::cp_kind::cage) cage = &p;
	ASSERT_NE(cage, nullptr);
	for (int k = 0; k < 3; k++) EXPECT_NEAR(cage->position[k], 0.0, 1E-6) << "cage point off the centre of symmetry";
	EXPECT_NEAR(cage->eigenvalues[0], cage->eigenvalues[2], 1E-8 * std::abs(cage->eigenvalues[2]) + 1E-10)
		<< "the Hessian at an O_h fixed point must be isotropic";
	EXPECT_GT(cage->eigenvalues[0], 0.0);

	//the eight ring points lie on the <111> rays, the twelve bond points on the <110> rays
	for (const topology::cp& p : r.points) {
		if (p.kind == topology::cp_kind::ring) {
			EXPECT_NEAR(std::abs(p.position[0]), std::abs(p.position[1]), 1E-6);
			EXPECT_NEAR(std::abs(p.position[1]), std::abs(p.position[2]), 1E-6);
			EXPECT_GT(std::abs(p.position[0]), 1E-3);
		}
		if (p.kind == topology::cp_kind::bond) {
			//one component zero, the other two equal in magnitude
			d3 m{ std::abs(p.position[0]), std::abs(p.position[1]), std::abs(p.position[2]) };
			std::sort(m.begin(), m.end());
			EXPECT_NEAR(m[0], 0.0, 1E-6);
			EXPECT_NEAR(m[1], m[2], 1E-6);
			EXPECT_GT(m[1], 1E-3);
		}
	}
}

//The topology must not depend on how diffuse the Gaussians are, only on the arrangement.  This is
//the check that caught the real bug: with an absolute gradient tolerance alone, seeds in the
//density tail satisfy |grad rho| <= 1E-7 without moving and are reported as critical points - on
//this density that invented 128 ring points and a Poincare-Hopf sum of 121.
TEST(Topology, OctahedronTopologyIsIndependentOfTheExponent)
{
	for (const auto& [a, d] : std::vector<std::pair<double, double>>{ { 1.0, 2.0 }, { 0.8, 2.0 }, { 1.0, 2.5 }, { 0.6, 2.0 } }) {
		SCOPED_TRACE("a = " + std::to_string(a) + ", d = " + std::to_string(d));
		const std::vector<topology::nucleus> nuc = octahedron_nuclei(d);
		const topology::result r = topology::analyze_topology(octahedron_source(nuc, a), nuc, gaussian_options());
		EXPECT_EQ(r.n_attractor, 6);
		EXPECT_EQ(r.n_bond, 12);
		EXPECT_EQ(r.n_ring, 8);
		EXPECT_EQ(r.n_cage, 1);
		EXPECT_EQ(r.sum, 1);
		EXPECT_TRUE(r.complete) << r.diagnosis;
	}
}

//Newton must reach a known saddle from a start nowhere near it: the cage point of the octahedron
//is at the origin, and a start 1.3 bohr away and off every symmetry element still lands on it.
TEST(Topology, NewtonReachesTheKnownSaddleFromAPoorStart)
{
	const double d = 2.0;
	const std::vector<topology::nucleus> nuc = octahedron_nuclei(d);
	const gaussian_sum src = octahedron_source(nuc, 1.0);
	const topology::options opt = gaussian_options();

	//The ring point on the (1,1,1) ray, located by the search, is the reference for the poor start
	const topology::result r = topology::analyze_topology(src, nuc, opt);
	d3 ring111{ 0.0, 0.0, 0.0 };
	bool have_ring = false;
	for (const topology::cp& p : r.points)
		if (p.kind == topology::cp_kind::ring && p.position[0] > 0 && p.position[1] > 0 && p.position[2] > 0) {
			ring111 = p.position;
			have_ring = true;
			break;
		}
	ASSERT_TRUE(have_ring);

	//A start far from the origin, on no symmetry element, must still find the cage point exactly
	d3 p{ 0.9, -0.7, 0.55 };
	int iterations = 0;
	ASSERT_TRUE(topology::newton_to_critical_point(src, p, opt, iterations));
	topology::cp point = topology::describe_cp(src, p, nuc, opt);
	//Newton converges on the nearest stationary point of any index; whichever of the two it picks,
	//it must be one of them exactly and not somewhere in between
	const bool at_cage = array_length(p, d3{ 0.0, 0.0, 0.0 }) < 1E-6;
	const bool at_ring = std::abs(std::abs(p[0]) - ring111[0]) < 1E-6 &&
		std::abs(std::abs(p[1]) - ring111[0]) < 1E-6 &&
		std::abs(std::abs(p[2]) - ring111[0]) < 1E-6;
	EXPECT_TRUE(at_cage || at_ring) << p[0] << " " << p[1] << " " << p[2];
	EXPECT_LT(point.gradient_norm, opt.gradient_tolerance);
	EXPECT_GT(iterations, 1) << "a poor start should take real iterations, not be accepted where it stands";

	//And a start deep in the tail must be rejected, not reported as a critical point
	d3 tail{ 40.0, 40.0, 40.0 };
	int tail_iterations = 0;
	EXPECT_FALSE(topology::newton_to_critical_point(src, tail, opt, tail_iterations))
		<< "a seed in the density tail has a tiny absolute gradient and is not a critical point";
}

//---------------------------------------------------------------------------------------------
// A real wavefunction with a ring
//---------------------------------------------------------------------------------------------
//Benzene would be the textbook case (12 - 12 + 1 - 0 = 1) but no benzene wavefunction ships with
//the tests.  tests/cytidine_tonto/cyt.wfn does the job and is checked against a count nobody typed
//in: for a connected molecular graph the number of bonds is fixed by Euler, E = V - 1 + (rings), so
//cytidine's 30 nuclei and 2 rings (pyrimidine + ribose) force 31 bond points and
//30 - 31 + 2 - 0 = 1.  Both the bond count and the sum are therefore predictions, not a golden file.
//
//epoxide.molden, the obvious small ring, cannot be used - see NotEveryWavefunctionHasNuclearMaxima
//below.  It is a valence-only pseudopotential wavefunction, so rho has no maximum at C or O at all.
TEST(Topology, PoincareHopfHoldsOnACyclicMolecule)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "cytidine_tonto" / "cyt.wfn";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wavy(f);
	const std::vector<topology::nucleus> nuc = topology::nuclei_of(wavy);
	ASSERT_EQ(nuc.size(), 30u) << "expected cytidine, C9H13N3O5";

	const topology::result r = topology::analyze_topology(wavy, nuc, topology::options{});
	SCOPED_TRACE(r.diagnosis);
	EXPECT_EQ(r.n_attractor, 30) << "one nuclear attractor per atom, none left over";
	EXPECT_EQ(r.n_ring, 2) << "the pyrimidine and the ribose ring";
	EXPECT_EQ(r.n_cage, 0) << "two separate rings enclose no cage";
	//E = V - 1 + rings for a connected graph, so the bond count is not a free parameter
	EXPECT_EQ(r.n_bond, r.n_attractor - 1 + r.n_ring);
	EXPECT_EQ(r.n_degenerate, 0);
	EXPECT_EQ(r.sum, 1) << "30 - 31 + 2 - 0";
	EXPECT_TRUE(r.complete) << r.diagnosis;
	EXPECT_FALSE(r.escalated) << "topology-driven seeding alone should close the sum on a molecule";
	EXPECT_EQ(r.graph_vertices, 30);
	EXPECT_EQ(r.graph_edges, 31);
	EXPECT_EQ(r.graph_components, 1) << "cytidine is one covalent fragment";
	EXPECT_EQ(r.required_ring_minus_cage, 2);

	//Both ring points must lie inside their ring, well away from every bond point
	for (const topology::cp& p : r.points) {
		if (p.kind != topology::cp_kind::ring) continue;
		EXPECT_GT(p.nearest_nucleus_distance, 1.0) << "a ring point sits in the middle of its ring";
		EXPECT_GT(p.density, 0.0);
		EXPECT_GT(p.laplacian, 0.0) << "two positive curvatures dominate at a ring point";
	}
	//No debris: every attractor belongs to a nucleus.  The tolerance is not the same for the two
	//kinds of centre - a heavy-atom maximum sits on the nucleus to machine precision, while a
	//hydrogen maximum in a finite Gaussian basis is displaced a tenth of a bohr towards its bond
	//partner, because the basis cannot reproduce the nuclear cusp of so light a nucleus.
	EXPECT_EQ(r.n_nna, 0);
	for (const topology::cp& p : r.points) {
		if (p.kind != topology::cp_kind::attractor) continue;
		ASSERT_GE(p.nearest_nucleus, 0);
		const double limit = nuc[p.nearest_nucleus].Z > 1 ? 0.01 : 0.2;
		EXPECT_LT(p.nearest_nucleus_distance, limit) << "attractor " << p.nearest_nucleus << " off its nucleus";
	}
}

//Not every file that reads as a wavefunction has a nuclear maximum, and the search must not pretend
//it does.  tests/molden_file/epoxide.molden carries 9 occupied orbitals holding 18 electrons for
//C2H4O, which has 24 - the heavy atoms wear pseudopotentials (largest oxygen s exponent 69, against
//the ~10^4 an all-electron oxygen needs, and the contraction coefficients start negative).  A
//valence-only density has a local *minimum* at each heavy nucleus, so the true topology of this
//density has (3,+3) points on C and O and maxima out in the valence shell.  That is what the search
//reports, and Poincare-Hopf correctly refuses to close: the value of the check is exactly that it
//says so instead of printing a plausible-looking table.
TEST(Topology, NotEveryWavefunctionHasNuclearMaxima)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "epoxide.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wavy(f);
	const std::vector<topology::nucleus> nuc = topology::nuclei_of(wavy);
	ASSERT_EQ(nuc.size(), 7u) << "expected C2H4O";
	//18 electrons for a molecule of 24: the cores are not in this file
	EXPECT_LT(wavy.count_nr_electrons(), 24) << "if this ever becomes an all-electron file, use it above";

	const topology::result r = topology::analyze_topology(wavy, nuc, topology::options{});
	SCOPED_TRACE(r.diagnosis);
	//rho at a bare pseudopotential nucleus is a minimum, not a maximum
	int cage_on_a_heavy_nucleus = 0;
	for (const topology::cp& p : r.points)
		if (p.kind == topology::cp_kind::cage && p.nearest_nucleus >= 0 &&
			p.nearest_nucleus_distance < 0.1 && nuc[p.nearest_nucleus].Z > 1) cage_on_a_heavy_nucleus++;
	EXPECT_EQ(cage_on_a_heavy_nucleus, 3) << "one density minimum on each of O, C and C";
	EXPECT_GT(r.n_nna, 0) << "the valence-shell maxima are attractors that belong to no nucleus";
	EXPECT_NE(r.sum, 1) << "the molecular Poincare-Hopf form cannot hold for a coreless density";
	EXPECT_FALSE(r.complete);
	EXPECT_FALSE(r.diagnosis.empty()) << "an incomplete result has to say what is missing";
}

//The diagnosis has to name the class that is short, not just print a wrong sum.  Drop the ring
//points from a complete result and the report must say ring seeding is the gap.
TEST(Topology, MissingRingPointIsDiagnosedAsARingGap)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "cytidine_tonto" / "cyt.wfn";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wavy(f);
	const std::vector<topology::nucleus> nuc = topology::nuclei_of(wavy);
	topology::result r = topology::analyze_topology(wavy, nuc, topology::options{});
	ASSERT_TRUE(r.complete) << r.diagnosis;

	r.points.erase(std::remove_if(r.points.begin(), r.points.end(),
		[](const topology::cp& p) { return p.kind == topology::cp_kind::ring; }), r.points.end());
	topology::tally(r);
	EXPECT_EQ(r.sum, -1) << "30 - 31 + 0 - 0";
	EXPECT_FALSE(r.balanced);
	EXPECT_FALSE(r.complete);
	EXPECT_NE(r.diagnosis.find("Ring seeding is the likely gap"), std::string::npos) << r.diagnosis;
	//and the graph invariant it is compared against is still there to be quoted
	EXPECT_EQ(r.required_ring_minus_cage, 2);
	EXPECT_EQ(r.found_ring_minus_cage, 0);
}

//A non-nuclear attractor is a (3,-3) point away from every nucleus, and it has to be reported with
//its position, its rho and its distance to the nearest nucleus: an integration that misses it
//assigns its basin's charge to a neighbour.  Constructed rather than hunted for - a Gaussian
//placed between two nuclei with no nucleus of its own is exactly that debris.
TEST(Topology, NonNuclearAttractorIsReportedWithItsDistance)
{
	gaussian_sum src;
	src.a = 1.5;
	//two "nuclei" far apart and a third Gaussian in the middle that no nucleus accounts for
	src.centres = { d3{ -3.0, 0.0, 0.0 }, d3{ 3.0, 0.0, 0.0 }, d3{ 0.0, 0.0, 0.0 } };
	const std::vector<topology::nucleus> nuc{ { { -3.0, 0.0, 0.0 }, 1 }, { { 3.0, 0.0, 0.0 }, 1 } };
	topology::options opt = gaussian_options();

	const topology::result r = topology::analyze_topology(src, nuc, opt);
	ASSERT_EQ(r.n_nna, 1) << "the unclaimed maximum at the origin must be flagged";
	const topology::cp* nna = nullptr;
	for (const topology::cp& p : r.points) if (p.is_nna) nna = &p;
	ASSERT_NE(nna, nullptr);
	EXPECT_EQ(nna->kind, topology::cp_kind::attractor);
	for (int k = 0; k < 3; k++) EXPECT_NEAR(nna->position[k], 0.0, 1E-6);
	EXPECT_NEAR(nna->nearest_nucleus_distance, 3.0, 1E-5);
	EXPECT_GT(nna->density, 0.0);
	//An NNA splits one bond path into two, so it adds one attractor and one bond point and leaves
	//the Poincare-Hopf sum at 1: the relation cannot be used to detect one
	EXPECT_EQ(r.n_attractor, 3);
	EXPECT_EQ(r.n_bond, 2);
	EXPECT_EQ(r.sum, 1);
	EXPECT_TRUE(r.complete) << r.diagnosis;
}

//The report must print the sum, the counts and, when it does not close, the diagnosis - the point
//of the exercise is that a user sees the gap rather than a plausible-looking number.
TEST(Topology, ReportNamesTheAssumedFormAndTheDeficit)
{
	topology::result r;
	topology::cp attractor;
	attractor.kind = topology::cp_kind::attractor;
	r.points.assign(2, attractor);
	r.graph_vertices = 2;
	r.graph_edges = 0;
	r.graph_components = 2;
	topology::tally(r);
	EXPECT_EQ(r.sum, 2);
	EXPECT_FALSE(r.balanced);

	std::ostringstream out;
	topology::report_topology(r, { { { 0.0, 0.0, 0.0 }, 1 }, { { 10.0, 0.0, 0.0 }, 1 } }, out);
	const std::string text = out.str();
	//which Poincare-Hopf form is being assumed, stated in the output
	EXPECT_NE(text.find("molecular form"), std::string::npos);
	EXPECT_NE(text.find("INCOMPLETE"), std::string::npos);
	EXPECT_NE(text.find("deficit"), std::string::npos);
	//two covalent fragments is a named possibility, not silently folded into a wrong sum
	EXPECT_NE(text.find("fragment"), std::string::npos) << text;
	//the thresholds that did the rejecting are printed too
	EXPECT_NE(text.find("Accepted when"), std::string::npos);
}
