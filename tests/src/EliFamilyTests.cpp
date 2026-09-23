#include "pch.h"

#include "core/eli_family.h"

//eli_family.h: Kohout's ELI-D for alpha-alpha, beta-beta and triplet-coupled pairs, and ELI-q.
//Two kinds of check: the analytic uniform-electron-gas limit, which fixes the normalisation
//constants without any reference program, and per-member regression rows taken off DGrid 5.2's
//own grids (the reference implementation) for a closed-shell and an open-shell wavefunction.
//The DGrid conventions these rows encode are written out in eli_family.h.
namespace
{
	//One voxel of a DGrid FIELD-DATA grid: position in bohr and DGrid's value for each member.
	struct dg_row
	{
		d3 p;
		double rho_a, elid_aa, elid_bb, eliq_aa, elid_t;
	};

	void check_against_dgrid(const char* molden, const std::vector<dg_row>& rows, const double tol)
	{
		const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / molden;
		if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
		WFN wave(f);
		const double factor = eli_family::triplet_density_factor(wave);
		for (const dg_row& r : rows) {
			eli_family::SpinFields s;
			eli_family::spin_fields(wave, r.p, s);
			const double ga = eli_family::g_same_spin(s, 0), gb = eli_family::g_same_spin(s, 1);
			SCOPED_TRACE(std::to_string(r.p[0]) + " " + std::to_string(r.p[1]) + " " + std::to_string(r.p[2]));
			EXPECT_NEAR(s.rho[0], r.rho_a, tol * std::abs(r.rho_a));
			EXPECT_NEAR(eli_family::eli_d(s.rho[0], ga), r.elid_aa, tol * std::abs(r.elid_aa));
			EXPECT_NEAR(eli_family::eli_d(s.rho[1], gb), r.elid_bb, tol * std::abs(r.elid_bb));
			EXPECT_NEAR(eli_family::eli_q(s.rho[0], ga), r.eliq_aa, tol * std::abs(r.eliq_aa));
			EXPECT_NEAR(eli_family::eli_d(factor * (s.rho[0] + s.rho[1]), eli_family::g_triplet(s)),
				r.elid_t, tol * std::abs(r.elid_t));
		}
	}
}

//The uniform electron gas: grad rho = 0 and T_s = (3/5)(6 pi^2)^(2/3) rho_s^(5/3), so
//g^s = (3/5)(6 pi^2)^(2/3) rho_s^(8/3) and ELI-D is the density-independent constant
//[12/((3/5)(6 pi^2)^(2/3))]^(3/8).  This pins every normalisation constant in the family.
TEST(EliFamily, UniformElectronGasLimit)
{
	const double c = eli_family::ueg_g_coefficient();
	EXPECT_NEAR(c, 0.6 * std::pow(6.0 * constants::PI2, 2.0 / 3.0), 1E-12);
	//[12/((3/5)(6 pi^2)^(2/3))]^(3/8) = 1.108596490563... - worked out independently of the header.
	EXPECT_NEAR(eli_family::ueg_eli_d(), 1.10859649056, 1E-10);
	for (const double rho_s : { 0.01, 0.1, 1.0, 12.5 }) {
		eli_family::SpinFields f;
		f.rho[0] = f.rho[1] = rho_s;
		f.T[0] = f.T[1] = c * std::pow(rho_s, 5.0 / 3.0);
		const double g = eli_family::g_same_spin(f, 0);
		EXPECT_NEAR(g, c * std::pow(rho_s, 8.0 / 3.0), 1E-10 * g);
		//ELI-D is density independent in the UEG, ELI-q is its -8/3 power
		EXPECT_NEAR(eli_family::eli_d(rho_s, g), eli_family::ueg_eli_d(), 1E-10);
		EXPECT_NEAR(eli_family::eli_q(rho_s, g), eli_family::ueg_eli_q(), 1E-10 * eli_family::ueg_eli_q());
		EXPECT_NEAR(eli_family::eli_q(rho_s, g), std::pow(eli_family::eli_d(rho_s, g), -8.0 / 3.0), 1E-10);
		//no spin polarisation and no density gradient: g^(t) = 3 g^alpha exactly
		EXPECT_NEAR(eli_family::g_triplet(f), 3.0 * g, 1E-10 * g);
	}
}

//ELI-q is an exact monotone rescaling of ELI-D for a single determinant, Y_q = Y_D^(-8/3);
//there is no extra information in it, only a different restriction scheme.  Checked on the
//real field rather than on the UEG, where every point has the same value.
TEST(EliFamily, EliQIsEliDToTheMinusEightThirds)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F_full.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	for (double x = -1.5; x <= 1.5; x += 0.5) {
		eli_family::SpinFields s;
		eli_family::spin_fields(wave, d3{ x, 0.3, -0.2 }, s);
		const double g = eli_family::g_same_spin(s, 0);
		const double d = eli_family::eli_d(s.rho[0], g);
		ASSERT_GT(d, 0.0);
		EXPECT_NEAR(eli_family::eli_q(s.rho[0], g), std::pow(d, -8.0 / 3.0), 1E-10 * std::pow(d, -8.0 / 3.0));
	}
}

//A restricted wavefunction has to come out spin-symmetric, and its ELI-D(alpha-alpha) has to be
//what WFN::computeELI already returns - the existing ELI-D is Kohout's Y_D^alpha.
TEST(EliFamily, RestrictedMatchesExistingEliD)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F2.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	for (double z = 0.0; z <= 2.0; z += 0.5) {
		const d3 p{ 0.4, -0.3, z };
		eli_family::SpinFields s;
		eli_family::spin_fields(wave, p, s);
		EXPECT_NEAR(s.rho[0], s.rho[1], 1E-12);
		EXPECT_NEAR(s.T[0], s.T[1], 1E-12);
		const double mine = eli_family::eli_d(s.rho[0], eli_family::g_same_spin(s, 0));
		EXPECT_NEAR(mine, wave.computeELI(p), 1E-8 * std::max(1.0, mine));
	}
}

//rho^(t)/rho = 1 - N_beta/(2(N-1)), DGrid 5.2's convention (Kohout's Eq. 5 gives a different
//normalisation; the header says so).  F2 has N = 14, N_beta = 7 -> 19/26.
TEST(EliFamily, TripletDensityFactor)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F2.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	EXPECT_NEAR(eli_family::triplet_density_factor(wave), 19.0 / 26.0, 1E-12);
	const std::filesystem::path g = nos_test_repo_root() / "tests" / "molden_file" / "F_open.molden";
	if (!std::filesystem::exists(g)) GTEST_SKIP() << "missing fixture " << g.string();
	WFN open(g);
	//N = 9, N_beta = 4 -> 1 - 4/16 = 3/4
	EXPECT_NEAR(eli_family::triplet_density_factor(open), 0.75, 1E-12);
}

//ELIA has no field function; the status string has to say why, so a user does not think it is
//merely unimplemented.
TEST(EliFamily, EliaReportsWhyItIsNotComputable)
{
	const std::string s = eli_family::elia_status();
	EXPECT_NE(s.find("2-matrix"), std::string::npos);
	EXPECT_NE(s.find("on-top"), std::string::npos);
}

//----------------------------------------------------------------------------------------------
//Regression rows straight off DGrid 5.2's own grids.  How to regrow them:
//  1. DGrid converts the molden first, and it picks the convention from the header - "program=
//     ORCA" -> <name>.molden.orca, contraction coefficients already carry the primitive
//     normalisation; "program= MOLDEN" -> <name>.molden.md, they multiply normalised primitives
//     and the contracted shell is renormalised.  NoSpherA2's molden reader assumes the ORCA
//     convention unconditionally, so ALL FOUR fixtures below are ORCA-convention files.  Do not
//     add a MOLDEN-convention one (F2.molden and epoxide.molden are two) until the reader tells
//     them apart - on F2.molden the two codes differ by 695x at the nucleus, which is a reader
//     question and not an ELI one.
//  2. compute / using wfn_1 / <property> / save field_1 / mesh=0.3 rho=0.0001, once per property:
//     "rho alpha", "ELI-D alpha-alpha", "ELI-D beta-beta", "ELI-q alpha-alpha", "ELI-D triplet-pair".
//  3. The FIELD-DATA header's origin and I/J/K columns are FULL SPANS, so step = vec/(n-1), and
//     the value list runs with I fastest, then J, then K.  (Confirmed by landing the global maxima
//     on the nuclei to the voxel.)  The points below are those voxel centres, in bohr.
//
//Measured over 400 voxels with rho > 1e-3, NoSpherA2 against DGrid, relative, max over the set:
//                    rho(median)  ELI-D aa   ELI-D bb   ELI-q aa   ELI-D triplet
//  F_open              6.2e-8      4.3e-5     5.9e-5     1.1e-4       4.7e-5
//  F_s1                7.5e-8      4.0e-5     2.2e-4     1.1e-4       5.0e-5
//  F_s32               8.9e-8      3.6e-5     2.1e-3     9.7e-5       4.2e-5
//  F_full              6.7e-8      4.3e-5    (skipped)   1.1e-4       4.3e-5
//The median ~6e-8 is DGrid's 9-significant-digit ASCII output; the tail is NoSpherA2's exp_cutoff
//primitive skipping, worst in the depleted beta channel of the S = 3/2 case.  The tolerances below
//are those worst cases with headroom - a regression fence, not a precision claim.
namespace
{
	//F_full: DGrid 5.2 on F_full.molden.orca, mesh 0.3, grid [25, 25, 25]; rho agrees to 1.7e-05 here.
	//Restricted, so DGrid builds no beta orbital set and its elid_r_b_bb is an artefact (its beta tau
	//is 0, making its g_beta = -1/4|grad rho_b|^2 negative).  The beta column below is DGrid's alpha
	//field, which is what beta must equal for a closed shell - that is the statement being pinned.
	const std::vector<dg_row> f_full_rows = {
		//{x, y, z}, rho_alpha, ELI-D(aa), ELI-D(bb), ELI-q(aa), ELI-D(triplet)
		{ {1.21000000, 0.91000000, 1.51000000}, 1.00107620e-02, 1.13291071e+00, 1.13291071e+00, 7.16932351e-01, 1.08386715e+00 },
		{ {1.21000000, 1.51000000, 0.01000000}, 1.69256091e-02, 1.21240590e+00, 1.21240590e+00, 5.98326938e-01, 1.15992101e+00 },
		{ {0.61000000, -1.49000000, -0.59000000}, 3.05759474e-02, 1.31125800e+00, 1.31125800e+00, 4.85472828e-01, 1.25449382e+00 },
		{ {1.21000000, 0.01000000, -0.59000000}, 8.88643700e-02, 1.48095384e+00, 1.48095384e+00, 3.50931706e-01, 1.41684354e+00 },
		{ {0.01000000, 0.01000000, 0.01000000}, 1.63065130e+02, 9.12533809e+00, 9.12533809e+00, 2.75002155e-03, 8.73030339e+00 },
	};

	//F_open: DGrid 5.2 on F_open.molden.orca, mesh 0.3, grid [25, 25, 24]; rho agrees to 3.3e-05 here.
	const std::vector<dg_row> f_open_rows = {
		//{x, y, z}, rho_alpha, ELI-D(aa), ELI-D(bb), ELI-q(aa), ELI-D(triplet)
		{ {1.81000000, 0.01000000, 1.05000000}, 1.12642273e-02, 1.15018569e+00, 1.23496210e+00, 6.88576362e-01, 1.15597674e+00 },
		{ {1.21000000, -0.29000000, -1.35000000}, 2.19825823e-02, 1.25504461e+00, 1.05686869e+00, 5.45643393e-01, 1.13180525e+00 },
		{ {-1.49000000, 0.01000000, -0.75000000}, 3.48282921e-02, 1.33371853e+00, 1.50072447e+00, 4.63975987e-01, 1.37260816e+00 },
		{ {-0.59000000, 0.31000000, 1.05000000}, 1.19785207e-01, 1.51611960e+00, 1.09743088e+00, 3.29643185e-01, 1.29753837e+00 },
		{ {0.01000000, 0.01000000, -0.15000000}, 1.58359749e+01, 2.17522564e+00, 2.81194455e+00, 1.25889207e-01, 2.39074421e+00 },
	};

	//F_s1: DGrid 5.2 on F_s1.molden.orca, mesh 0.3, grid [24, 25, 24]; rho agrees to 1.0e-05 here.
	const std::vector<dg_row> f_s1_rows = {
		//{x, y, z}, rho_alpha, ELI-D(aa), ELI-D(bb), ELI-q(aa), ELI-D(triplet)
		{ {1.35000000, 0.31000000, 1.35000000}, 1.69621244e-02, 1.21274567e+00, 6.15420583e-01, 5.97880031e-01, 9.97031252e-01 },
		{ {-1.05000000, 0.31000000, -1.35000000}, 2.86606937e-02, 1.30011234e+00, 7.24630856e-01, 4.96650624e-01, 1.09063940e+00 },
		{ {-1.05000000, 0.91000000, -0.75000000}, 4.49159808e-02, 1.37709351e+00, 1.66646709e+00, 4.26020849e-01, 1.37746539e+00 },
		{ {1.05000000, -0.59000000, 0.15000000}, 1.30524867e-01, 1.52459937e+00, 1.57006885e+00, 3.24776578e-01, 1.47913254e+00 },
		{ {-0.15000000, 0.01000000, 0.15000000}, 5.94428704e+00, 1.32990696e+00, 1.84285316e+00, 4.67530534e-01, 1.53370550e+00 },
	};

	//F_s32: DGrid 5.2 on F_s32.molden.orca, mesh 0.3, grid [24, 24, 24]; rho agrees to 1.4e-05 here.
	const std::vector<dg_row> f_s32_rows = {
		//{x, y, z}, rho_alpha, ELI-D(aa), ELI-D(bb), ELI-q(aa), ELI-D(triplet)
		{ {-1.35000000, -0.15000000, -1.35000000}, 1.78317647e-02, 1.22073774e+00, 8.27634687e+00, 5.87498831e-01, 1.15110125e+00 },
		{ {-1.05000000, -0.15000000, -1.35000000}, 3.03984416e-02, 1.31025136e+00, 8.34532181e+00, 4.86468077e-01, 1.26242499e+00 },
		{ {-0.45000000, 0.45000000, 1.35000000}, 5.77760063e-02, 1.41823925e+00, 9.51689524e+00, 3.93853459e-01, 1.39873752e+00 },
		{ {-0.45000000, -1.05000000, -0.15000000}, 1.55327685e-01, 1.53911620e+00, 1.37842696e+01, 3.16671932e-01, 1.56318908e+00 },
		{ {-0.15000000, -0.15000000, -0.15000000}, 3.09115361e+00, 9.77527412e-01, 1.38395867e+00, 1.06248501e+00, 1.16173152e+00 },
	};
}

//Closed shell, N = 10.
TEST(EliFamilyDGrid, RestrictedFminus) { check_against_dgrid("F_full.molden", f_full_rows, 1E-3); }

//Doublet, N = 9.  aa, bb and the triplet are three genuinely different fields here.
TEST(EliFamilyDGrid, DoubletF) { check_against_dgrid("F_open.molden", f_open_rows, 1E-3); }

//Triplet cation, N = 8, S = 1.
TEST(EliFamilyDGrid, TripletFplus) { check_against_dgrid("F_s1.molden", f_s1_rows, 1E-3); }

//Quartet dication, N = 7, S = 3/2.  Its beta channel is depleted enough that NoSpherA2's exp_cutoff
//primitive skipping shows at the 2e-3 level, hence the looser fence.
TEST(EliFamilyDGrid, QuartetFpp) { check_against_dgrid("F_s32.molden", f_s32_rows, 5E-3); }

//The member list is read off the orbitals, never off a stated multiplicity.  A restricted file gets
//ELI-D(aa) and ELI-q(aa) and nothing else - ELI-D(bb) is the same field and ELI-D(triplet) a constant
//multiple of it - while an unrestricted one gets all six.
TEST(EliFamily, VariantSelectionFollowsTheOrbitals)
{
	const std::filesystem::path root = nos_test_repo_root() / "tests" / "molden_file";
	if (!std::filesystem::exists(root / "F_full.molden") || !std::filesystem::exists(root / "F_open.molden"))
		GTEST_SKIP() << "missing fixtures in " << root.string();
	WFN closed(root / "F_full.molden");
	std::string warn;
	std::vector<eli_family::Member> v = eli_family::eli_variants_for(closed, &warn);
	EXPECT_EQ(v, (std::vector<eli_family::Member>{ eli_family::Member::eli_d_aa, eli_family::Member::eli_q_aa }));

	WFN open(root / "F_open.molden");
	v = eli_family::eli_variants_for(open, &warn);
	EXPECT_EQ(v.size(), 6u);
	EXPECT_NE(std::find(v.begin(), v.end(), eli_family::Member::eli_d_triplet), v.end());
	EXPECT_NE(std::find(v.begin(), v.end(), eli_family::Member::elia_singlet), v.end());

	//Only ELI-D has maxima of its own.  ELI-q is ELI-D^(-8/3), so an ascent on it runs into ELI-D's
	//minima and shatters into tail basins - DGrid's own ELID_core ascent returns ~60 (H2O) / ~100
	//(F doublet) of them, holding 0.0000-0.0005 e each.
	EXPECT_TRUE(eli_family::basins_independent(eli_family::Member::eli_d_aa));
	EXPECT_TRUE(eli_family::basins_independent(eli_family::Member::eli_d_triplet));
	EXPECT_FALSE(eli_family::basins_independent(eli_family::Member::eli_q_aa));
	EXPECT_FALSE(eli_family::basins_independent(eli_family::Member::elia_singlet));
	EXPECT_STREQ(eli_family::member_name(eli_family::Member::eli_d_triplet), "eli_d_triplet");
	EXPECT_STREQ(eli_family::member_column(eli_family::Member::eli_d_triplet), "ELI-D_t");
}

//A multiplicity label that disagrees with the occupations must not change the member list, only
//raise a warning.  F_full is restricted, so its ELI family is the closed-shell one whatever the
//label claims; inventing a spin split from a label is how a spin heuristic doubles every index
//while the populations still look right.
TEST(EliFamily, MultiplicityLabelIsOnlyACrossCheck)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F_full.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	wave.set_multi(3);   //lie about it
	std::string warn;
	const std::vector<eli_family::Member> v = eli_family::eli_variants_for(wave, &warn);
	EXPECT_EQ(v, (std::vector<eli_family::Member>{ eli_family::Member::eli_d_aa, eli_family::Member::eli_q_aa }));
	EXPECT_NE(warn.find("occupations"), std::string::npos);
	EXPECT_NE(warn.find("restricted"), std::string::npos);
}

//ELIA / the singlet ELI-q of [III] Eq. 53 collapses to 1 - zeta^2 for any single determinant, so it
//is identically 1 for a closed shell and carries no pair information at all.  Pinned so nobody
//mistakes it for a second opinion on the electron pair structure.
TEST(EliFamily, SingletEliQIsOneMinusZetaSquared)
{
	eli_family::SpinFields f{};
	f.rho[0] = f.rho[1] = 0.37;
	EXPECT_NEAR(eli_family::elia_singlet_eli_q(f), 1.0, 1E-14);
	f.rho[0] = 0.8049; f.rho[1] = 0.2214;
	const double rho = f.rho[0] + f.rho[1], zeta = (f.rho[0] - f.rho[1]) / rho;
	EXPECT_NEAR(eli_family::elia_singlet_eli_q(f), 1.0 - zeta * zeta, 1E-14);
	f.rho[0] = f.rho[1] = 0.0;
	EXPECT_EQ(eli_family::elia_singlet_eli_q(f), 0.0);   //no density, no value

	const std::filesystem::path g = nos_test_repo_root() / "tests" / "molden_file" / "F_full.molden";
	if (!std::filesystem::exists(g)) GTEST_SKIP() << "missing fixture " << g.string();
	WFN closed(g);
	for (double z = 0.2; z <= 2.0; z += 0.6) {
		eli_family::SpinFields s;
		eli_family::spin_fields(closed, d3{ 0.3, -0.2, z }, s);
		EXPECT_NEAR(eli_family::elia_singlet_eli_q(s), 1.0, 1E-12);
	}
}
