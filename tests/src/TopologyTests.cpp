#include "pch.h"

#include "core/topology.h"

//References whose topology symmetry fixes, not a golden file: six equal Gaussians on an octahedron
//give 6 attractors, 12 bond, 8 ring and 1 cage point (vertices, edges, faces, centre), 6 - 12 + 8 - 1 = 1.
namespace
{
	//rho = sum_i exp(-a |p - R_i|^2) with analytic gradient and Hessian, a PointDensityWithHessian
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

	//bond_scale 1.3 seeds no H pair beyond 1.1 bohr; widened so these densities get bond seeds at all.
	//It only adds seeds, bonding is read off the bond points found.
	topology::options gaussian_options()
	{
		topology::options opt;
		opt.bond_scale = 10.0;
		return opt;
	}
}

TEST(Topology, ClassifyHessianSignatures)
{
	d3 w;
	int neg = 0, pos = 0, zero = 0;
	const double a[9]{ -1.0, 0, 0, 0, -2.0, 0, 0, 0, -3.0 };
	EXPECT_EQ(topology::classify_hessian(a, w, neg, pos, zero), topology::cp_kind::attractor);
	EXPECT_EQ(neg, 3);
	EXPECT_EQ(zero, 0);
	//eigenvalues ascending, as the ellipticity lambda1/lambda2 - 1 assumes
	EXPECT_NEAR(w[0], -3.0, 1E-12);
	EXPECT_NEAR(w[2], -1.0, 1E-12);
	const double b[9]{ -1.0, 0, 0, 0, -2.0, 0, 0, 0, 0.5 };
	EXPECT_EQ(topology::classify_hessian(b, w, neg, pos, zero), topology::cp_kind::bond);
	EXPECT_EQ(neg, 2);
	EXPECT_EQ(pos, 1);
	const double c[9]{ -1.0, 0, 0, 0, 2.0, 0, 0, 0, 0.5 };
	EXPECT_EQ(topology::classify_hessian(c, w, neg, pos, zero), topology::cp_kind::ring);
	EXPECT_EQ(neg, 1);
	EXPECT_EQ(pos, 2);
	const double e[9]{ 1.0, 0, 0, 0, 2.0, 0, 0, 0, 0.5 };
	EXPECT_EQ(topology::classify_hessian(e, w, neg, pos, zero), topology::cp_kind::cage);
	EXPECT_EQ(pos, 3);
	//not diagonal, so the eigensolver is read rather than the diagonal
	const double s2 = std::sqrt(0.5);
	//R diag(-1, 2, 0.5) R^T, R a 45 deg rotation about z
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

//A zero eigenvalue must not be binned, or Poincare-Hopf could close for the wrong reason
TEST(Topology, DegenerateHessianIsReportedNotBinned)
{
	d3 w;
	int neg = 0, pos = 0, zero = 0;
	//(-, -, 0), two coalescing maxima: an attractor if rounded down, a bond point if rounded up
	const double h[9]{ -2.0, 0, 0, 0, -2.0, 0, 0, 0, 0.0 };
	EXPECT_EQ(topology::classify_hessian(h, w, neg, pos, zero), topology::cp_kind::degenerate);
	EXPECT_EQ(zero, 1);
	EXPECT_EQ(neg, 2);
	//2E-7 against max|lambda| = 2 is below the default 1E-6 relative tolerance
	const double h2[9]{ -2.0, 0, 0, 0, -2.0, 0, 0, 0, 2E-7 };
	EXPECT_EQ(topology::classify_hessian(h2, w, neg, pos, zero), topology::cp_kind::degenerate);
	const double h3[9]{ -2.0, 0, 0, 0, -2.0, 0, 0, 0, 2E-5 };
	EXPECT_EQ(topology::classify_hessian(h3, w, neg, pos, zero), topology::cp_kind::bond);
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

//Two equal Gaussians 2/sqrt(2a) apart: H_xx = 0 at the midpoint, H_yy = H_zz = -4a exp(-a s^2)
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
	const double transverse = -4.0 * a * std::exp(-a * s * s);
	EXPECT_NEAR(point.eigenvalues[0], transverse, 1E-10);
	EXPECT_NEAR(point.eigenvalues[1], transverse, 1E-10);
	EXPECT_NEAR(point.eigenvalues[2], 0.0, 1E-12) << "the longitudinal curvature vanishes at s = 1/sqrt(2a)";
	EXPECT_EQ(point.kind, topology::cp_kind::degenerate);
	EXPECT_EQ(point.zero, 1);

	const topology::result r = topology::analyze_topology(src, nuc, gaussian_options());
	EXPECT_GT(r.n_degenerate, 0);
	EXPECT_FALSE(r.complete);
	EXPECT_NE(r.diagnosis.find("degenerate"), std::string::npos) << r.diagnosis;
}

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
	//cycle rank of the bond graph, 12 - 6 + 1 = 7 = n_ring - n_cage
	EXPECT_EQ(r.required_ring_minus_cage, 7);
	EXPECT_EQ(r.found_ring_minus_cage, 7);

	//O_h fixes only the origin, so the cage point is there with an isotropic Hessian
	const topology::cp* cage = nullptr;
	for (const topology::cp& p : r.points) if (p.kind == topology::cp_kind::cage) cage = &p;
	ASSERT_NE(cage, nullptr);
	for (int k = 0; k < 3; k++) EXPECT_NEAR(cage->position[k], 0.0, 1E-6) << "cage point off the centre of symmetry";
	EXPECT_NEAR(cage->eigenvalues[0], cage->eigenvalues[2], 1E-8 * std::abs(cage->eigenvalues[2]) + 1E-10)
		<< "the Hessian at an O_h fixed point must be isotropic";
	EXPECT_GT(cage->eigenvalues[0], 0.0);

	//ring points on the <111> rays, bond points on the <110> rays
	for (const topology::cp& p : r.points) {
		if (p.kind == topology::cp_kind::ring) {
			EXPECT_NEAR(std::abs(p.position[0]), std::abs(p.position[1]), 1E-6);
			EXPECT_NEAR(std::abs(p.position[1]), std::abs(p.position[2]), 1E-6);
			EXPECT_GT(std::abs(p.position[0]), 1E-3);
		}
		if (p.kind == topology::cp_kind::bond) {
			d3 m{ std::abs(p.position[0]), std::abs(p.position[1]), std::abs(p.position[2]) };
			std::sort(m.begin(), m.end());
			EXPECT_NEAR(m[0], 0.0, 1E-6);
			EXPECT_NEAR(m[1], m[2], 1E-6);
			EXPECT_GT(m[1], 1E-3);
		}
	}
}

//With an absolute gradient tolerance alone, tail seeds pass |grad rho| <= 1E-7 without moving,
//so a change of exponent would invent critical points
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

TEST(Topology, NewtonReachesTheKnownSaddleFromAPoorStart)
{
	const double d = 2.0;
	const std::vector<topology::nucleus> nuc = octahedron_nuclei(d);
	const gaussian_sum src = octahedron_source(nuc, 1.0);
	const topology::options opt = gaussian_options();

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

	//a start off every symmetry element
	d3 p{ 0.9, -0.7, 0.55 };
	int iterations = 0;
	ASSERT_TRUE(topology::newton_to_critical_point(src, p, opt, iterations));
	topology::cp point = topology::describe_cp(src, p, nuc, opt);
	//Newton reaches the nearest stationary point of any index, so either one, exactly
	const bool at_cage = array_length(p, d3{ 0.0, 0.0, 0.0 }) < 1E-6;
	const bool at_ring = std::abs(std::abs(p[0]) - ring111[0]) < 1E-6 &&
		std::abs(std::abs(p[1]) - ring111[0]) < 1E-6 &&
		std::abs(std::abs(p[2]) - ring111[0]) < 1E-6;
	EXPECT_TRUE(at_cage || at_ring) << p[0] << " " << p[1] << " " << p[2];
	EXPECT_LT(point.gradient_norm, opt.gradient_tolerance);
	EXPECT_GT(iterations, 1) << "a poor start should take real iterations, not be accepted where it stands";

	d3 tail{ 40.0, 40.0, 40.0 };
	int tail_iterations = 0;
	EXPECT_FALSE(topology::newton_to_critical_point(src, tail, opt, tail_iterations))
		<< "a seed in the density tail has a tiny absolute gradient and is not a critical point";
}

//---------------------------------------------------------------------------------------------
// A real wavefunction with a ring
//---------------------------------------------------------------------------------------------
//For a connected graph E = V - 1 + rings, so cytidine's 30 nuclei and 2 rings force 31 bond points
//and 30 - 31 + 2 - 0 = 1; both are predictions, not a golden file.
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
	EXPECT_EQ(r.n_bond, r.n_attractor - 1 + r.n_ring);
	EXPECT_EQ(r.n_degenerate, 0);
	EXPECT_EQ(r.sum, 1) << "30 - 31 + 2 - 0";
	EXPECT_TRUE(r.complete) << r.diagnosis;
	EXPECT_FALSE(r.escalated) << "topology-driven seeding alone should close the sum on a molecule";
	EXPECT_EQ(r.graph_vertices, 30);
	EXPECT_EQ(r.graph_edges, 31);
	EXPECT_EQ(r.graph_components, 1) << "cytidine is one covalent fragment";
	EXPECT_EQ(r.required_ring_minus_cage, 2);

	for (const topology::cp& p : r.points) {
		if (p.kind != topology::cp_kind::ring) continue;
		EXPECT_GT(p.nearest_nucleus_distance, 1.0) << "a ring point sits in the middle of its ring";
		EXPECT_GT(p.density, 0.0);
		EXPECT_GT(p.laplacian, 0.0) << "two positive curvatures dominate at a ring point";
	}
	//A heavy-atom maximum sits on its nucleus to machine precision; a hydrogen one is displaced towards
	//its bond partner, as a finite Gaussian basis cannot reproduce so light a cusp.
	EXPECT_EQ(r.n_nna, 0);
	for (const topology::cp& p : r.points) {
		if (p.kind != topology::cp_kind::attractor) continue;
		ASSERT_GE(p.nearest_nucleus, 0);
		const double limit = nuc[p.nearest_nucleus].Z > 1 ? 0.01 : 0.2;
		EXPECT_LT(p.nearest_nucleus_distance, limit) << "attractor " << p.nearest_nucleus << " off its nucleus";
	}
}

//Seeds are searched in parallel and merged serially in seed order, so the result is bitwise thread-independent
TEST(Topology, ThreadCountDoesNotChangeThePointSet)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "cytidine_tonto" / "cyt.wfn";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wavy(f);
	const std::vector<topology::nucleus> nuc = topology::nuclei_of(wavy);
	const int threads = omp_get_max_threads();
	omp_set_num_threads(1);
	const topology::result serial = topology::analyze_topology(wavy, nuc, topology::options{});
	omp_set_num_threads(std::max(threads, 4));
	const topology::result parallel = topology::analyze_topology(wavy, nuc, topology::options{});
	omp_set_num_threads(threads);
	ASSERT_EQ(serial.points.size(), parallel.points.size());
	for (size_t i = 0; i < serial.points.size(); i++) {
		EXPECT_EQ(serial.points[i].kind, parallel.points[i].kind) << "point " << i;
		EXPECT_EQ(serial.points[i].from, parallel.points[i].from) << "point " << i;
		for (int k = 0; k < 3; k++) EXPECT_EQ(serial.points[i].position[k], parallel.points[i].position[k]) << "point " << i;
	}
	EXPECT_EQ(serial.sum, parallel.sum);
}

//rho has a cusp maximum at every nucleus of an all-electron density, whatever Poincare-Hopf adds up to.
//The heavy nuclei need a gradient tolerance scaled to rho: at Hessians of 1E8-1E9 an absolute 1E-7 is
//below the floating-point floor of the analytic gradient.
TEST(Topology, EveryNucleusOfAnAllElectronDensityIsAnAttractor)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "Fe_gbw" / "Fe.gbw";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wavy(f);
	const std::vector<topology::nucleus> nuc = topology::nuclei_of(wavy);
	ASSERT_EQ(nuc.size(), 21u) << "expected the iron thiolate, Fe(SCH3)4";
	double charge = 0.0;
	for (const topology::nucleus& n : nuc) charge += n.Z;
	//an ion is fine, a missing core is not: the smallest def2 ECP core is 10 electrons
	ASSERT_LT(charge - wavy.count_nr_electrons(), 5.0)
		<< "this fixture has to be all-electron for the cusp argument to apply";

	const topology::result r = topology::analyze_topology(wavy, nuc, topology::options{});
	SCOPED_TRACE(r.diagnosis);
	std::vector<int> attractors_of(nuc.size(), 0);
	for (const topology::cp& p : r.points)
		if (p.kind == topology::cp_kind::attractor && !p.is_nna && p.nearest_nucleus >= 0)
			attractors_of[p.nearest_nucleus]++;
	for (size_t a = 0; a < nuc.size(); a++)
		EXPECT_EQ(attractors_of[a], 1) << "nucleus " << a + 1 << " (Z=" << nuc[a].Z
			<< ") has no maximum of its own";
	//rho at the Fe and S nuclei is three orders above a carbon's
	int heavy = 0;
	for (const topology::cp& p : r.points)
		if (p.kind == topology::cp_kind::attractor && p.nearest_nucleus >= 0 && nuc[p.nearest_nucleus].Z > 15) {
			heavy++;
			EXPECT_GT(p.density, 1E3) << "Z=" << nuc[p.nearest_nucleus].Z;
			EXPECT_LT(p.nearest_nucleus_distance, 1E-3) << "a heavy maximum sits on its nucleus";
		}
	EXPECT_EQ(heavy, 5) << "one iron and four sulfurs";
}

//epoxide.molden is a pseudopotential wavefunction (18 of 24 electrons): rho has a minimum, a (3,+3)
//point, at C and O and maxima in the valence shell, so Poincare-Hopf must refuse to close.
TEST(Topology, NotEveryWavefunctionHasNuclearMaxima)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "epoxide.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wavy(f);
	const std::vector<topology::nucleus> nuc = topology::nuclei_of(wavy);
	ASSERT_EQ(nuc.size(), 7u) << "expected C2H4O";
	//18 electrons for 24: the cores are not in this file
	EXPECT_LT(wavy.count_nr_electrons(), 24) << "if this ever becomes an all-electron file, use it above";

	const topology::result r = topology::analyze_topology(wavy, nuc, topology::options{});
	SCOPED_TRACE(r.diagnosis);
	//rho at a bare pseudopotential nucleus is a minimum
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
	EXPECT_EQ(r.required_ring_minus_cage, 2);
	EXPECT_EQ(r.found_ring_minus_cage, 0);
}

//A Gaussian between two nuclei with no nucleus of its own is a non-nuclear attractor; missing it in an
//integration gives its basin's charge to a neighbour.
TEST(Topology, NonNuclearAttractorIsReportedWithItsDistance)
{
	gaussian_sum src;
	src.a = 1.5;
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
	//An NNA splits a bond path in two, adding one attractor and one bond point, so the sum cannot detect it
	EXPECT_EQ(r.n_attractor, 3);
	EXPECT_EQ(r.n_bond, 2);
	EXPECT_EQ(r.sum, 1);
	EXPECT_TRUE(r.complete) << r.diagnosis;
}

//One bonded pair with one attractor found: sum 0 against 1, a real deficit
TEST(Topology, ReportNamesTheAssumedFormAndTheDeficit)
{
	topology::result r;
	topology::cp p;
	p.kind = topology::cp_kind::attractor;
	r.points.push_back(p);
	p.kind = topology::cp_kind::bond;
	r.points.push_back(p);
	r.graph_vertices = 2;
	r.graph_edges = 1;
	r.graph_components = 1;
	r.covalent_fragments = 1;
	topology::tally(r);
	EXPECT_EQ(r.sum, 0);
	EXPECT_EQ(r.target, 1);
	EXPECT_FALSE(r.balanced);

	std::ostringstream out;
	topology::report_topology(r, { { { 0.0, 0.0, 0.0 }, 1 }, { { 1.4, 0.0, 0.0 }, 1 } }, out);
	const std::string text = out.str();
	EXPECT_NE(text.find("molecular form"), std::string::npos);
	EXPECT_NE(text.find("INCOMPLETE"), std::string::npos);
	EXPECT_NE(text.find("deficit"), std::string::npos) << text;
	EXPECT_NE(text.find("attractor or ring point"), std::string::npos) << text;
	EXPECT_NE(text.find("Accepted when"), std::string::npos);
}

//water.wfx is centrosymmetric about its Mn for 45 of 48 nuclei (the fifth water has no image), and every
//nucleus carries a cusp maximum, so a nucleus with a maximum has an image with one.  Two accepted points
//closer than merge_distance but of different signature are two points: merging them would let a bond
//seed replace a nuclear maximum on one side only.
TEST(Topology, ACentrosymmetricDensityHasACentrosymmetricCriticalPointSet)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "grown" / "water.wfx";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wavy(f);
	const std::vector<topology::nucleus> nuc = topology::nuclei_of(wavy);
	ASSERT_EQ(nuc.size(), 48u) << "expected the grown manganese complex with five waters";
	//centre from the fixture's own Mn coordinates: a literal truncated to seven decimals pairs nothing
	ASSERT_EQ(nuc[0].Z, 25) << "nucleus 1 is the inversion centre of this fixture";
	const d3 centre = nuc[0].pos;
	std::vector<int> partner(nuc.size(), -1);
	for (size_t a = 0; a < nuc.size(); a++) {
		const d3 im{ 2 * centre[0] - nuc[a].pos[0], 2 * centre[1] - nuc[a].pos[1], 2 * centre[2] - nuc[a].pos[2] };
		for (size_t b = 0; b < nuc.size(); b++)
			if (nuc[b].Z == nuc[a].Z && array_length(im, nuc[b].pos) < 1E-8) { partner[a] = (int)b; break; }
	}
	const size_t paired = (size_t)std::count_if(partner.begin(), partner.end(), [](const int p) { return p >= 0; });
	ASSERT_EQ(paired, 45u) << "the inversion centre or the fixture changed: 45 of the 48 nuclei pair "
		"under it, the three that do not being the unpaired fifth water";

	const topology::options opt{};
	const topology::result r = topology::analyze_topology(wavy, nuc, opt);
	SCOPED_TRACE(r.diagnosis);
	std::vector<int> attractors_of(nuc.size(), 0);
	std::vector<double> rho_of(nuc.size(), 0.0), doff_of(nuc.size(), 0.0);
	for (const topology::cp& p : r.points)
		if (p.kind == topology::cp_kind::attractor && !p.is_nna && p.nearest_nucleus >= 0) {
			attractors_of[p.nearest_nucleus]++;
			//the attractor nearest the nucleus decides, not a second one in the same basin
			if (rho_of[p.nearest_nucleus] == 0.0 || p.nearest_nucleus_distance < doff_of[p.nearest_nucleus]) {
				rho_of[p.nearest_nucleus] = p.density;
				doff_of[p.nearest_nucleus] = p.nearest_nucleus_distance;
			}
		}
	double worst_rho_rel = 0.0;
	std::string worst_pair = "none";
	for (size_t a = 0; a < nuc.size(); a++) {
		const int b = partner[a];
		if (b < 0) continue;
		EXPECT_EQ(attractors_of[a], attractors_of[(size_t)b])
			<< "nucleus " << a + 1 << " and its inversion image " << b + 1 << " (both Z=" << nuc[a].Z
			<< ") disagree on whether they own a maximum";
		if (attractors_of[a] > 0 && attractors_of[(size_t)b] > 0 && rho_of[a] > 0.0) {
			const double rel = std::abs(rho_of[a] - rho_of[(size_t)b]) / rho_of[a];
			if (rel > worst_rho_rel) {
				worst_rho_rel = rel;
				std::ostringstream o;
				o << nuc[a].Z << ": " << a + 1 << " rho " << rho_of[a] << " at " << doff_of[a]
					<< " bohr off its nucleus, image " << b + 1 << " rho " << rho_of[(size_t)b] << " at "
					<< doff_of[(size_t)b] << " bohr";
				worst_pair = o.str();
			}
		}
	}
	std::cout << "worst inversion pair, relative drho " << std::scientific << worst_rho_rel
		<< "   Z=" << worst_pair << "\n";

	//rho at each nucleus against its image: the geometry pairs them to 1E-13 bohr, so any difference
	//belongs to the wavefunction, an SCF on a cluster that is not itself centrosymmetric.
	double worst_nucleus_rel = 0.0;
	std::string worst_nucleus_pair;
	for (size_t a = 0; a < nuc.size(); a++) {
		const int b = partner[a];
		if (b < 0 || (size_t)b <= a) continue;
		const double rho_a = topology::describe_cp(wavy, nuc[a].pos, nuc, opt).density;
		const double rho_b = topology::describe_cp(wavy, nuc[(size_t)b].pos, nuc, opt).density;
		const double rel = std::abs(rho_a - rho_b) / rho_a;
		if (rel > worst_nucleus_rel) {
			worst_nucleus_rel = rel;
			std::ostringstream o;
			o << a + 1 << "/" << b + 1 << " Z=" << nuc[a].Z << " rho " << rho_a << " vs " << rho_b;
			worst_nucleus_pair = o.str();
		}
	}
	std::cout << "rho at a nucleus against rho at its image, 22 pairs, no search: worst relative drho "
		<< worst_nucleus_rel << "   " << worst_nucleus_pair << "\n";
	ASSERT_GT(worst_nucleus_rel, 0.0) << "no pair was compared, so the bound below is vacuous";

	//The search may not add asymmetry beyond the density's own at those nuclei; the factor 3 is headroom.
	EXPECT_LT(worst_rho_rel, 3.0 * worst_nucleus_rel) << "rho at a nuclear maximum differs from its "
		"image's by more than the density does at those nuclei, worst maximum pair Z=" << worst_pair
		<< ", worst nucleus pair " << worst_nucleus_pair;
	for (size_t a = 0; a < nuc.size(); a++)
		if (partner[a] >= 0)
			EXPECT_EQ(attractors_of[a], 1) << "nucleus " << a + 1 << " (Z=" << nuc[a].Z << ")";
}

//Poincare-Hopf is necessary, not sufficient: a spurious bond and ring point cancel in the sum.  HgH2's
//bond graph is a path with cycle rank 0, so n_ring - n_cage must be 0.  Counts entered by hand via tally().
namespace
{
	//fragments defaults to C; a closed-shell contact gives fragments > C, a missing bridging bond point < C
	topology::result graph_case(int n_attractor, int n_bond, int n_ring, int V, int E, int C, int fragments = -1)
	{
		topology::result r;
		topology::cp p;
		for (int k = 0; k < n_attractor; k++) { p.kind = topology::cp_kind::attractor; r.points.push_back(p); }
		for (int k = 0; k < n_bond; k++) { p.kind = topology::cp_kind::bond; r.points.push_back(p); }
		for (int k = 0; k < n_ring; k++) { p.kind = topology::cp_kind::ring; r.points.push_back(p); }
		r.graph_vertices = V;
		r.graph_edges = E;
		r.graph_components = C;
		r.covalent_fragments = fragments >= 0 ? fragments : std::max(C, 1);
		//as bond_graph_cycles() computes it, which needs point positions these have none of
		r.required_ring_minus_cage = E - V + C;
		topology::tally(r);
		return r;
	}
}

TEST(Topology, ABalancedSumThatContradictsTheBondGraphIsNotComplete)
{
	//HgH2: 3 nuclei, 2 bonds, one component, cycle rank 0, yet two ring points
	topology::result r = graph_case(3, 4, 2, 3, 2, 1);
	ASSERT_EQ(r.sum, 1);
	ASSERT_TRUE(r.balanced) << "the arm is vacuous unless the alternating sum still closes";
	ASSERT_EQ(r.required_ring_minus_cage, 0);
	ASSERT_EQ(r.found_ring_minus_cage, 2);
	EXPECT_FALSE(r.graph_consistent);
	EXPECT_FALSE(r.complete) << "a path graph has no ring, so 2 ring points with a balanced sum is "
		"two spurious bond points and two spurious ring points cancelling: " << r.diagnosis;
	EXPECT_NE(r.diagnosis.find("cycle rank"), std::string::npos) << r.diagnosis;

	//-topology's exit code is the flag, so the printed verdict must agree
	std::ostringstream out;
	topology::report_topology(r, { { { 0.0, 0.0, 0.0 }, 80 }, { { 3.1, 0.0, 0.0 }, 1 }, { { -3.1, 0.0, 0.0 }, 1 } }, out);
	EXPECT_NE(out.str().find("INCOMPLETE"), std::string::npos) << out.str();
}

TEST(Topology, ANucleusWithoutAnAttractorIsNotComplete)
{
	//Au2Br2's shape scaled down: the sum closes while a nucleus has no attractor, and every nucleus with
	//core electrons has a cusp, so n_attractor must be n_nuclei + n_NNA
	topology::result r = graph_case(3, 2, 0, 4, 3, 1);
	ASSERT_EQ(r.sum, 1);
	ASSERT_TRUE(r.balanced);
	ASSERT_EQ(r.required_ring_minus_cage, 0);
	ASSERT_EQ(r.found_ring_minus_cage, 0) << "the ring arm must not be what fails here";
	EXPECT_FALSE(r.graph_consistent);
	EXPECT_FALSE(r.complete) << r.diagnosis;
	EXPECT_NE(r.diagnosis.find("nuclear seeding"), std::string::npos) << r.diagnosis;
}

//An ECP core is named from rho at the nucleus, not from an ECP flag the file may not record.
TEST(Topology, AnEcpCoreIsNamedInsteadOfBlamedOnTheSeeding)
{
	topology::result r = graph_case(3, 2, 0, 4, 3, 1);
	ASSERT_FALSE(r.complete) << "the arm is vacuous unless the set is refused";
	ASSERT_NE(r.diagnosis.find("nuclear seeding"), std::string::npos) << r.diagnosis;
	EXPECT_EQ(r.diagnosis.find("pseudopotential"), std::string::npos)
		<< "an all-electron shortfall must NOT be excused as an ECP: " << r.diagnosis;

	topology::result ecp = r;
	ecp.coreless_nuclei = { 0 };
	topology::tally(ecp);
	EXPECT_FALSE(ecp.complete) << "naming the cause must not turn the refusal into a pass";
	EXPECT_NE(ecp.diagnosis.find("pseudopotential"), std::string::npos) << ecp.diagnosis;
	EXPECT_NE(ecp.diagnosis.find("atom number 1"), std::string::npos)
		<< "the note has to say WHICH nucleus, 1-based as the table prints them: " << ecp.diagnosis;

	//epoxide's counts with two coreless heavy atoms: complete, so no note on every ECP run
	topology::result fine = graph_case(7, 7, 1, 7, 7, 1);
	fine.coreless_nuclei = { 0, 3 };
	topology::tally(fine);
	ASSERT_TRUE(fine.complete) << fine.diagnosis;
	EXPECT_TRUE(fine.diagnosis.empty()) << fine.diagnosis;
}

TEST(Topology, AGenuineRingStaysComplete)
{
	//epoxide: 7 nuclei, 7 bonds in one component, cycle rank 1, one ring point
	topology::result r = graph_case(7, 7, 1, 7, 7, 1);
	ASSERT_EQ(r.sum, 1);
	ASSERT_EQ(r.required_ring_minus_cage, 1);
	EXPECT_TRUE(r.graph_consistent);
	EXPECT_TRUE(r.complete) << r.diagnosis;

	//no bond graph at all (one atom, or no nuclei) is not refused either
	topology::result lone = graph_case(1, 0, 0, 0, 0, 0);
	EXPECT_TRUE(lone.graph_consistent);
	EXPECT_TRUE(lone.complete) << lone.diagnosis;
}

//C isolated molecules with rho -> 0 between them sum to C; TFVC/water.gbw is water plus a distant He.
TEST(Topology, SeparatedFragmentsSumToTheirNumber)
{
	//water + He: 4 nuclei, 2 bond points, 2 components, 2 covalent fragments
	topology::result r = graph_case(4, 2, 0, 4, 2, 2);
	EXPECT_EQ(r.target, 2) << "two separated fragments have index sum 2";
	EXPECT_EQ(r.sum, 2);
	EXPECT_TRUE(r.balanced);
	EXPECT_TRUE(r.graph_consistent);
	EXPECT_TRUE(r.complete) << r.diagnosis;

	//one fragment still expects 1, so the target is not simply relaxed
	topology::result one = graph_case(4, 3, 0, 4, 3, 1);
	EXPECT_EQ(one.target, 1);
	EXPECT_TRUE(one.complete) << one.diagnosis;

	//an explicit target survives: the periodic Morse sum is 0 however many fragments
	topology::result r2;
	topology::cp p;
	for (int k = 0; k < 4; k++) { p.kind = topology::cp_kind::attractor; r2.points.push_back(p); }
	for (int k = 0; k < 2; k++) { p.kind = topology::cp_kind::bond; r2.points.push_back(p); }
	r2.graph_vertices = 4; r2.graph_edges = 2; r2.graph_components = 2; r2.covalent_fragments = 2;
	topology::options morse;
	morse.poincare_hopf_target = 0;
	topology::tally(r2, morse);
	EXPECT_EQ(r2.target, 0) << "an explicit Poincare-Hopf target must not be overwritten";
}

//The target comes from the found bond paths, since an H-bonded dimer sums to 1.  Dropping a bridging bond
//point raises the sum and splits a component together, so sum == C cannot see it; bond paths may be
//less disconnected than the covalent graph, never more.
TEST(Topology, FoundBondPathsMayNotBeMoreDisconnectedThanTheGeometry)
{
	//one covalent molecule of 4 nuclei, 2 bond points found, paths in 2 pieces
	topology::result missing = graph_case(4, 2, 0, 4, 2, 2, 1);
	EXPECT_EQ(missing.target, 2);
	EXPECT_TRUE(missing.balanced) << "the arm is vacuous unless the sum still closes against 2";
	EXPECT_FALSE(missing.graph_consistent);
	EXPECT_FALSE(missing.complete) << missing.diagnosis;
	EXPECT_NE(missing.diagnosis.find("bridging bond point"), std::string::npos) << missing.diagnosis;

	//a dimer held by a closed-shell contact: 2 covalent fragments, 1 connected set of paths, sum 1
	topology::result hbond = graph_case(6, 5, 0, 6, 5, 1, 2);
	EXPECT_EQ(hbond.target, 1);
	EXPECT_TRUE(hbond.graph_consistent);
	EXPECT_TRUE(hbond.complete) << hbond.diagnosis;
}

//He sits 6.3 A from the nearest H, beyond bond_scale 1.3 for any covalent radius in the table
TEST(Topology, CovalentFragmentCountIsGeometryOnly)
{
	const std::vector<topology::nucleus> water_and_he{
		{ { 0.0, 0.0, 0.2318 }, 8 }, { { 0.0, 1.39815, -0.89997 }, 1 }, { { 0.0, -1.39815, -0.89997 }, 1 },
		{ { 0.0, 13.22808, 0.0 }, 2 } };   //tests/TFVC/water.gbw's geometry, in bohr
	topology::options opt;
	EXPECT_EQ(topology::covalent_fragment_count(water_and_he, opt), 2);

	const std::vector<topology::nucleus> water(water_and_he.begin(), water_and_he.begin() + 3);
	EXPECT_EQ(topology::covalent_fragment_count(water, opt), 1);
	EXPECT_EQ(topology::covalent_fragment_count({}, opt), 1) << "an empty system is not a deficit";
	EXPECT_EQ(topology::covalent_fragment_count({ water_and_he[0] }, opt), 1);
	EXPECT_EQ(topology::covalent_fragment_count({ water_and_he[1], water_and_he[3] }, opt), 2);

	//bond_scale is the knob, not a hard-coded distance
	opt.bond_scale = 10.0;
	EXPECT_EQ(topology::covalent_fragment_count(water_and_he, opt), 1);
}
