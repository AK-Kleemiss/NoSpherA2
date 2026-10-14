#include "pch.h"

#include "core/eli_family.h"
#include "core/b2c.h"

//Kohout's ELI family: the UEG limit fixes the normalisation constants analytically, regression rows
//come from DGrid 5.2 grids (the conventions they encode are in eli_family.h).
namespace
{
	//One DGrid FIELD-DATA voxel: position in bohr and DGrid's value per member.
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

//UEG: grad rho = 0, T_s = (3/5)(6 pi^2)^(2/3) rho_s^(5/3), so g^s = (3/5)(6 pi^2)^(2/3) rho_s^(8/3)
//and ELI-D = [12/((3/5)(6 pi^2)^(2/3))]^(3/8), density independent.
TEST(EliFamily, UniformElectronGasLimit)
{
	const double c = eli_family::ueg_g_coefficient();
	EXPECT_NEAR(c, 0.6 * std::pow(6.0 * constants::PI2, 2.0 / 3.0), 1E-12);
	//[12/((3/5)(6 pi^2)^(2/3))]^(3/8) = 1.108596490563..., worked out independently of the header
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

//Y_q = Y_D^(-8/3) for a single determinant; checked on a real field, the UEG value is constant.
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

//WFN::computeELI is Kohout's Y_D^alpha.
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

//rho^(t)/rho = 1 - N_beta/(2(N-1)), DGrid 5.2's convention, not Kohout's Eq. 5. F2: N = 14, N_beta = 7.
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

TEST(EliFamily, EliaReportsWhyItIsNotComputable)
{
	const std::string s = eli_family::elia_status();
	EXPECT_NE(s.find("2-matrix"), std::string::npos);
	EXPECT_NE(s.find("on-top"), std::string::npos);
}

//Rows off DGrid 5.2 FIELD-DATA grids: compute / using wfn_1 / <property> / save field_1 / mesh=0.3 rho=0.0001.
//DGrid guesses the molden coefficient convention from [Title]; these fixtures are orca_2mkl output with the
//ORCA title, so both codes read them bare (MoldenConventionTests decides it for a new fixture). FIELD-DATA
//spans are full, step = vec/(n-1), I fastest. The tolerances fence exp_cutoff primitive skipping.
namespace
{
	//F_full.molden.orca, grid [25, 25, 25]. Restricted, so DGrid's own bb field is an artefact (beta tau 0);
	//the bb column is its alpha field, which beta must equal for a closed shell.
	const std::vector<dg_row> f_full_rows = {
		//{x, y, z}, rho_alpha, ELI-D(aa), ELI-D(bb), ELI-q(aa), ELI-D(triplet)
		{ {1.21000000, 0.91000000, 1.51000000}, 1.00107620e-02, 1.13291071e+00, 1.13291071e+00, 7.16932351e-01, 1.08386715e+00 },
		{ {1.21000000, 1.51000000, 0.01000000}, 1.69256091e-02, 1.21240590e+00, 1.21240590e+00, 5.98326938e-01, 1.15992101e+00 },
		{ {0.61000000, -1.49000000, -0.59000000}, 3.05759474e-02, 1.31125800e+00, 1.31125800e+00, 4.85472828e-01, 1.25449382e+00 },
		{ {1.21000000, 0.01000000, -0.59000000}, 8.88643700e-02, 1.48095384e+00, 1.48095384e+00, 3.50931706e-01, 1.41684354e+00 },
		{ {0.01000000, 0.01000000, 0.01000000}, 1.63065130e+02, 9.12533809e+00, 9.12533809e+00, 2.75002155e-03, 8.73030339e+00 },
	};

	//F_open.molden.orca, grid [25, 25, 24]
	const std::vector<dg_row> f_open_rows = {
		//{x, y, z}, rho_alpha, ELI-D(aa), ELI-D(bb), ELI-q(aa), ELI-D(triplet)
		{ {1.81000000, 0.01000000, 1.05000000}, 1.12642273e-02, 1.15018569e+00, 1.23496210e+00, 6.88576362e-01, 1.15597674e+00 },
		{ {1.21000000, -0.29000000, -1.35000000}, 2.19825823e-02, 1.25504461e+00, 1.05686869e+00, 5.45643393e-01, 1.13180525e+00 },
		{ {-1.49000000, 0.01000000, -0.75000000}, 3.48282921e-02, 1.33371853e+00, 1.50072447e+00, 4.63975987e-01, 1.37260816e+00 },
		{ {-0.59000000, 0.31000000, 1.05000000}, 1.19785207e-01, 1.51611960e+00, 1.09743088e+00, 3.29643185e-01, 1.29753837e+00 },
		{ {0.01000000, 0.01000000, -0.15000000}, 1.58359749e+01, 2.17522564e+00, 2.81194455e+00, 1.25889207e-01, 2.39074421e+00 },
	};

	//F_s1.molden.orca, grid [24, 25, 24]
	const std::vector<dg_row> f_s1_rows = {
		//{x, y, z}, rho_alpha, ELI-D(aa), ELI-D(bb), ELI-q(aa), ELI-D(triplet)
		{ {1.35000000, 0.31000000, 1.35000000}, 1.69621244e-02, 1.21274567e+00, 6.15420583e-01, 5.97880031e-01, 9.97031252e-01 },
		{ {-1.05000000, 0.31000000, -1.35000000}, 2.86606937e-02, 1.30011234e+00, 7.24630856e-01, 4.96650624e-01, 1.09063940e+00 },
		{ {-1.05000000, 0.91000000, -0.75000000}, 4.49159808e-02, 1.37709351e+00, 1.66646709e+00, 4.26020849e-01, 1.37746539e+00 },
		{ {1.05000000, -0.59000000, 0.15000000}, 1.30524867e-01, 1.52459937e+00, 1.57006885e+00, 3.24776578e-01, 1.47913254e+00 },
		{ {-0.15000000, 0.01000000, 0.15000000}, 5.94428704e+00, 1.32990696e+00, 1.84285316e+00, 4.67530534e-01, 1.53370550e+00 },
	};

	//F_s32.molden.orca, grid [24, 24, 24]
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

//Doublet, N = 9: aa, bb and triplet are three different fields.
TEST(EliFamilyDGrid, DoubletF) { check_against_dgrid("F_open.molden", f_open_rows, 1E-3); }

//Triplet cation, N = 8, S = 1.
TEST(EliFamilyDGrid, TripletFplus) { check_against_dgrid("F_s1.molden", f_s1_rows, 1E-3); }

//Quartet dication, N = 7, S = 3/2: exp_cutoff skipping shows in the depleted beta channel.
TEST(EliFamilyDGrid, QuartetFpp) { check_against_dgrid("F_s32.molden", f_s32_rows, 5E-3); }

//Members follow the orbitals: restricted gets ELI-D(aa) and ELI-q(aa) only (bb is the same field, the
//triplet a constant multiple), unrestricted all six.
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

	//ELI-q = ELI-D^(-8/3) has no maxima of its own: an ascent on it shatters into tail basins.
	EXPECT_TRUE(eli_family::basins_independent(eli_family::Member::eli_d_aa));
	EXPECT_TRUE(eli_family::basins_independent(eli_family::Member::eli_d_triplet));
	EXPECT_FALSE(eli_family::basins_independent(eli_family::Member::eli_q_aa));
	EXPECT_FALSE(eli_family::basins_independent(eli_family::Member::elia_singlet));
	EXPECT_STREQ(eli_family::member_name(eli_family::Member::eli_d_triplet), "eli_d_triplet");
	EXPECT_STREQ(eli_family::member_column(eli_family::Member::eli_d_triplet), "ELI-D_t");
}

//One electron: the bb members, the triplet (rho^(t) = 0 below two electrons) and the singlet ELI-q
//(1 - zeta^2 = 0) are identically zero and must not be listed as members.
TEST(EliFamily, AnEmptySpinChannelIsNotAMember)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "ptb_H_file" / "H.gbw";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	ASSERT_NEAR(wave.count_nr_electrons(), 1.0, 1e-9) << "this fixture is the one-electron case";
	std::string warn;
	const std::vector<eli_family::Member> v = eli_family::eli_variants_for(wave, &warn);
	EXPECT_EQ(v, (std::vector<eli_family::Member>{ eli_family::Member::eli_d_aa,
		eli_family::Member::eli_q_aa }));
	EXPECT_NE(warn.find("holds no electrons"), std::string::npos) << warn;
	//rho^(t) = 0 below two electrons
	EXPECT_DOUBLE_EQ(eli_family::triplet_density_factor(wave), 0.0);
	//the bound is two electrons, not "unrestricted": the F doublet keeps all six
	const std::filesystem::path open = nos_test_repo_root() / "tests" / "molden_file" / "F_open.molden";
	if (!std::filesystem::exists(open)) GTEST_SKIP() << "missing fixture " << open.string();
	WFN doublet(open);
	EXPECT_EQ(eli_family::eli_variants_for(doublet, &warn).size(), 6u);
	EXPECT_GT(eli_family::triplet_density_factor(doublet), 0.0);
}

//A multiplicity label that disagrees with the occupations only warns; a spin split taken from a label
//doubles every index while the populations still look right.
TEST(EliFamily, MultiplicityLabelIsOnlyACrossCheck)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F_full.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	wave.set_multi(3);
	std::string warn;
	const std::vector<eli_family::Member> v = eli_family::eli_variants_for(wave, &warn);
	EXPECT_EQ(v, (std::vector<eli_family::Member>{ eli_family::Member::eli_d_aa, eli_family::Member::eli_q_aa }));
	EXPECT_NE(warn.find("occupations"), std::string::npos);
	EXPECT_NE(warn.find("restricted"), std::string::npos);
}

//The singlet ELI-q of [III] Eq. 53 is 1 - zeta^2 for any single determinant: 1 for a closed shell.
TEST(EliFamily, SingletEliQIsOneMinusZetaSquared)
{
	eli_family::SpinFields f{};
	f.rho[0] = f.rho[1] = 0.37;
	EXPECT_NEAR(eli_family::elia_singlet_eli_q(f), 1.0, 1E-14);
	f.rho[0] = 0.8049; f.rho[1] = 0.2214;
	const double rho = f.rho[0] + f.rho[1], zeta = (f.rho[0] - f.rho[1]) / rho;
	EXPECT_NEAR(eli_family::elia_singlet_eli_q(f), 1.0 - zeta * zeta, 1E-14);
	f.rho[0] = f.rho[1] = 0.0;
	EXPECT_EQ(eli_family::elia_singlet_eli_q(f), 0.0);

	const std::filesystem::path g = nos_test_repo_root() / "tests" / "molden_file" / "F_full.molden";
	if (!std::filesystem::exists(g)) GTEST_SKIP() << "missing fixture " << g.string();
	WFN closed(g);
	for (double z = 0.2; z <= 2.0; z += 0.6) {
		eli_family::SpinFields s;
		eli_family::spin_fields(closed, d3{ 0.3, -0.2, z }, s);
		EXPECT_NEAR(eli_family::elia_singlet_eli_q(s), 1.0, 1E-12);
	}
}

//computeELISpinGrad: value = the eli_family member, gradient = its central difference, aux[3] = rho_s * Y_q.
namespace
{
	const d3 spin_points[] = { { 0.3, -0.2, 0.5 }, { 1.1, 0.4, -0.7 }, { 0.05, 0.1, 0.15 }, { 2.0, 0.5, 1.0 }, { -0.6, 0.9, 0.2 } };
}

TEST(EliSpin, PointValuesMatchTheFamilyMembers)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F_open.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	const double tf = eli_family::triplet_density_factor(wave);
	for (const d3 &p : spin_points) {
		eli_family::SpinFields s;
		eli_family::spin_fields(wave, p, s);
		const double ref[3] = { eli_family::eli_d(s.rho[0], eli_family::g_same_spin(s, 0)), eli_family::eli_d(s.rho[1], eli_family::g_same_spin(s, 1)),
			eli_family::eli_d(tf * (s.rho[0] + s.rho[1]), eli_family::g_triplet(s)) };
		for (int field = 0; field < 3; field++) {
			double y, aux[4];
			d3 g;
			wave.computeELISpinGrad(p, field, tf, y, g, aux);
			SCOPED_TRACE("field " + std::to_string(field) + " at " + std::to_string(p[0]) + " " + std::to_string(p[1]) + " " + std::to_string(p[2]));
			EXPECT_NEAR(y, ref[field], 1E-9 * std::max(1.0, ref[field]));
			EXPECT_NEAR(aux[1], s.rho[0], 1E-10 * std::max(1.0, s.rho[0]));
			EXPECT_NEAR(aux[2], s.rho[1], 1E-10 * std::max(1.0, s.rho[1]));
			EXPECT_NEAR(aux[0], s.rho[0] + s.rho[1], 1E-10 * std::max(1.0, aux[0]));
			if (field < 2) {
				const double q = eli_family::eli_q(s.rho[field], eli_family::g_same_spin(s, field));
				EXPECT_NEAR(aux[3] / s.rho[field], q, 1E-9 * std::max(1.0, q));
			}
			else
				EXPECT_EQ(aux[3], 0.0);
		}
	}
}

TEST(EliSpin, AnalyticGradientMatchesFiniteDifference)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F_open.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	const double tf = eli_family::triplet_density_factor(wave), h = 1E-4;
	for (const d3 &p : spin_points)
		for (int field = 0; field < 3; field++) {
			double y;
			d3 g;
			wave.computeELISpinGrad(p, field, tf, y, g);
			for (int k = 0; k < 3; k++) {
				d3 a = p, b = p, ga, gb;
				a[k] += h;
				b[k] -= h;
				double ya, yb;
				wave.computeELISpinGrad(a, field, tf, ya, ga);
				wave.computeELISpinGrad(b, field, tf, yb, gb);
				const double fd = (ya - yb) / (2 * h);
				EXPECT_NEAR(g[k], fd, 1E-6 * std::max(1.0, std::abs(fd))) << "field " << field << " axis " << k << " at " << p[0] << " " << p[1] << " " << p[2];
			}
		}
}

//Restricted: both channels carry occ/2, so ELI-D(alpha-alpha) = ELI-D(beta-beta) = computeELI
TEST(EliSpin, RestrictedChannelsAreTheExistingEliD)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F_full.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	const double tf = eli_family::triplet_density_factor(wave);
	for (const d3 &p : spin_points) {
		const double ref = wave.computeELI(p);
		for (int field = 0; field < 2; field++) {
			double y;
			d3 g;
			wave.computeELISpinGrad(p, field, tf, y, g);
			EXPECT_NEAR(y, ref, 1E-8 * std::max(1.0, ref)) << "field " << field;
		}
	}
}

//The F doublet's basins hold every alpha and beta electron up to what the quadrature leaves outside.
TEST(EliSpin, BasinSpinPopulationsAddUpToTheElectronCounts)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F_open.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	const double tf = eli_family::triplet_density_factor(wave);
	for (int field = 0; field < 3; field++) {
		const eli_spin_field eval = [&wave, field, tf](const d3 &p, double &y, d3 &g, double *aux) { wave.computeELISpinGrad(p, field, tf, y, g, aux); };
		const std::vector<d4> maxima = analytic_eli_maxima(wave, false, &eval);
		ASSERT_FALSE(maxima.empty());
		vec volumes;
		vec2 spin;
		double outside = 0.0;
		const vec pop = integrate_basins_on_atomic_grids(nullptr, nullptr, maxima, wave, 2, true, volumes, outside, nullptr, nullptr, 1, nullptr, nullptr, nullptr, &eval, &spin);
		ASSERT_EQ(spin.size(), pop.size());
		double na = 0.0, nb = 0.0, n = 0.0;
		for (size_t b = 0; b < pop.size(); b++) {
			na += spin[b][0];
			nb += spin[b][1];
			n += pop[b];
			EXPECT_NEAR(spin[b][0] + spin[b][1], pop[b], 1E-9 * std::max(1.0, pop[b])) << "field " << field << " basin " << b;
		}
		EXPECT_NEAR(na, 5.0, 5E-3) << "field " << field;
		EXPECT_NEAR(nb, 4.0, 5E-3) << "field " << field;
		EXPECT_NEAR(n + outside, 9.0, 2E-3) << "field " << field;
	}
}

//NH3Li, UKS: points near each nucleus, off-axis and on bonds, where several centres contribute.
TEST(EliSpin, MoleculeValuesAndGradientsOnNH3Li)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	ASSERT_GT(wave.get_MO_op_count(1), 0);
	const double tf = eli_family::triplet_density_factor(wave), h = 1E-5;
	std::vector<d3> pts;
	const d3 p0 = wave.get_atom_pos(0);
	for (int a = 0; a < wave.get_ncen(); a++) {
		const d3 pa = wave.get_atom_pos(a);
		pts.push_back({ pa[0] + 0.05, pa[1] + 0.03, pa[2] - 0.04 });
		pts.push_back({ pa[0] + 0.4, pa[1] - 0.3, pa[2] + 0.5 });
		if (a) pts.push_back({ 0.5 * (pa[0] + p0[0]) + 0.02, 0.5 * (pa[1] + p0[1]), 0.5 * (pa[2] + p0[2]) - 0.03 });
	}
	for (const d3 &p : pts) {
		eli_family::SpinFields s;
		eli_family::spin_fields(wave, p, s);
		const double ref[3] = { eli_family::eli_d(s.rho[0], eli_family::g_same_spin(s, 0)), eli_family::eli_d(s.rho[1], eli_family::g_same_spin(s, 1)),
			eli_family::eli_d(tf * (s.rho[0] + s.rho[1]), eli_family::g_triplet(s)) };
		for (int field = 0; field < 3; field++) {
			SCOPED_TRACE("field " + std::to_string(field) + " at " + std::to_string(p[0]) + " " + std::to_string(p[1]) + " " + std::to_string(p[2]));
			double y;
			d3 g;
			wave.computeELISpinGrad(p, field, tf, y, g);
			EXPECT_NEAR(y, ref[field], 1E-9 * std::max(1.0, ref[field]));
			for (int k = 0; k < 3; k++) {
				d3 a = p, b = p, ga, gb;
				a[k] += h;
				b[k] -= h;
				double ya, yb;
				wave.computeELISpinGrad(a, field, tf, ya, ga);
				wave.computeELISpinGrad(b, field, tf, yb, gb);
				const double fd = (ya - yb) / (2 * h);
				EXPECT_NEAR(g[k], fd, 1E-5 * std::max(1.0, std::abs(fd))) << "axis " << k;
			}
		}
	}
}

//ROKS: a singly occupied MO is alpha, a doubly occupied one carries one electron of each spin. Rb doublet:
//N_alpha 19, N_beta 18, and the alpha excess is the 5s density.
TEST(EliSpin, RestrictedOpenShellIsSplitByOccupation)
{
	double n[2];
	eli_family::mo_spin_occupations(2.0, 0, eli_family::SpinSplit::restricted_open, n);
	EXPECT_EQ(n[0], 1.0); EXPECT_EQ(n[1], 1.0);
	eli_family::mo_spin_occupations(1.0, 0, eli_family::SpinSplit::restricted_open, n);
	EXPECT_EQ(n[0], 1.0); EXPECT_EQ(n[1], 0.0);
	eli_family::mo_spin_occupations(1.0, 0, eli_family::SpinSplit::halves, n);
	EXPECT_EQ(n[0], 0.5); EXPECT_EQ(n[1], 0.5);
	eli_family::mo_spin_occupations(1.0, 1, eli_family::SpinSplit::unrestricted, n);
	EXPECT_EQ(n[0], 0.0); EXPECT_EQ(n[1], 1.0);

	const std::filesystem::path f = nos_test_repo_root() / "tests" / "ECP_SF" / "Rb.gbw";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	ASSERT_EQ(eli_family::spin_split(wave), eli_family::SpinSplit::restricted_open);
	std::string warn;
	EXPECT_EQ(eli_family::eli_variants_for(wave, &warn).size(), 6u) << warn;
	EXPECT_NEAR(eli_family::triplet_density_factor(wave), 1.0 - 18.0 / (2.0 * 36.0), 1E-12);
	const double tf = eli_family::triplet_density_factor(wave);
	for (const d3 &p : { d3{ 0.0, 0.0, 4.5 }, d3{ 1.0, -2.0, 3.0 }, d3{ 0.2, 0.1, -0.3 } }) {
		eli_family::SpinFields s;
		eli_family::spin_fields(wave, p, s);
		double y, aux[4];
		d3 g;
		wave.computeELISpinGrad(p, 0, tf, y, g, aux);
		EXPECT_NEAR(aux[1], s.rho[0], 1E-10 * std::max(1.0, s.rho[0]));
		EXPECT_NEAR(aux[2], s.rho[1], 1E-10 * std::max(1.0, s.rho[1]));
		EXPECT_NEAR(aux[1] + aux[2], wave.compute_dens(p), 1E-9 * std::max(1.0, aux[0]));
		EXPECT_GT(aux[1], aux[2]);
		EXPECT_NEAR(y, eli_family::eli_d(s.rho[0], eli_family::g_same_spin(s, 0)), 1E-9 * std::max(1.0, y));
	}
}

//Rb alpha-alpha: core, n = 4 shell, and a degenerate sphere of maxima just inside the rho = 1e-4
//isosurface; the merges make each one basin.
static int rb_alpha_alpha_basins(const WFN &wave, vec &pop, vec2 &spin, double &outside)
{
	const double tf = eli_family::triplet_density_factor(wave);
	const eli_spin_field eval = [&wave, tf](const d3 &p, double &y, d3 &g, double *aux) { wave.computeELISpinGrad(p, 0, tf, y, g, aux); };
	const std::vector<d4> all = analytic_eli_maxima(wave, false, &eval);
	std::vector<d4> maxima = all;
	cubei none;
	ivec core_map, shell_map, edge_map;
	unify_core_basins(none, maxima, *wave.get_atoms_ptr(), &core_map);
	unify_shell_basins(none, maxima, &shell_map, 1.2, 0.05, wave.get_atoms_ptr());
	for (size_t b = 1; b < core_map.size(); b++) core_map[b] = shell_map[core_map[b]];
	if (unify_boundary_basins(maxima, wave, eval, &edge_map) > 0)
		for (size_t b = 1; b < core_map.size(); b++) core_map[b] = edge_map[core_map[b]];
	vec volumes;
	pop = integrate_basins_on_atomic_grids(nullptr, nullptr, all, wave, 2, true, volumes, outside, nullptr, nullptr, 1, nullptr, nullptr, &core_map, &eval, &spin);
	return static_cast<int>(maxima.size());
}

TEST(EliSpin, RestrictedOpenShellSphereIsOneBasin)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "ECP_SF" / "Rb.gbw";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	vec pop;
	vec2 spin;
	double outside = 0.0;
	EXPECT_EQ(rb_alpha_alpha_basins(wave, pop, spin, outside), 3);
	ASSERT_EQ(spin.size(), pop.size());
	double n = 0.0;
	for (size_t b = 0; b < pop.size(); b++) {
		n += pop[b];
		EXPECT_NEAR(spin[b][0] + spin[b][1], pop[b], 1E-9 * std::max(1.0, pop[b])) << "basin " << b;
	}
	EXPECT_NEAR(n + outside, wave.count_nr_electrons(), 5E-3);
	//outside is the tail beyond the isosurface
	EXPECT_LT(outside, 0.15);
}

//Broken symmetry: with N_alpha = N_beta only the orbitals tell a spin-polarised determinant from a
//restricted one written out twice.
TEST(EliFamily, BrokenSymmetryIsDecidedByTheOrbitals)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F_full.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN wave(f);
	const d3 p{ 0.3, -0.2, 0.5 };
	const double ref = wave.computeELI(p);
	const int n0 = wave.get_nmo();
	for (int i = 0; i < n0; i++) {
		MO a = wave.get_MO(i);
		a.set_occ(0.5 * wave.get_MO_occ(i));
		a.set_op(0);
		wave.push_back_MO(a);
	}
	for (int i = 0; i < n0; i++) {
		MO b = wave.get_MO(n0 + i);
		b.set_op(1);
		wave.push_back_MO(b);
	}
	for (int i = 0; i < n0; i++) wave.delete_MO(0);
	ASSERT_EQ(wave.get_nmo(), 2 * n0);
	ASSERT_EQ(eli_family::spin_split(wave), eli_family::SpinSplit::unrestricted);
	EXPECT_FALSE(eli_family::alpha_beta_orbitals_differ(wave));
	std::string warn;
	EXPECT_EQ(eli_family::eli_variants_for(wave, &warn).size(), 2u) << warn;
	double y;
	d3 g;
	wave.computeELISpinGrad(p, 0, eli_family::triplet_density_factor(wave), y, g);
	EXPECT_NEAR(y, ref, 1E-8 * std::max(1.0, ref));

	//the beta HOMO's largest coefficient, 2 % off: same electron counts, different orbitals
	int homo = -1, big = 0;
	for (int m = 0; m < wave.get_nmo(); m++)
		if (wave.get_MO_op(m) == 1 && wave.get_MO_occ(m) > 0.0) homo = m;
	ASSERT_GE(homo, 0);
	for (int j = 1; j < wave.get_nex(); j++)
		if (std::abs(wave.get_MO_coef(homo, j)) > std::abs(wave.get_MO_coef(homo, big))) big = j;
	wave.set_MO_coef(homo, big, 1.02 * wave.get_MO_coef(homo, big));
	EXPECT_TRUE(eli_family::alpha_beta_orbitals_differ(wave));
	EXPECT_EQ(eli_family::eli_variants_for(wave, &warn).size(), 6u) << warn;
	//a sign flip of a whole orbital is the same orbital
	wave.set_MO_coef(homo, big, wave.get_MO_coef(homo, big) / 1.02);
	for (int j = 0; j < wave.get_nex(); j++) wave.set_MO_coef(homo, j, -wave.get_MO_coef(homo, j));
	EXPECT_FALSE(eli_family::alpha_beta_orbitals_differ(wave));
	//two equally occupied beta orbitals mixed by 45 degrees: every orbital differs, rho_beta does not
	int lo = -1;
	for (int m = 0; m < homo; m++)
		if (wave.get_MO_op(m) == 1 && wave.get_MO_occ(m) == wave.get_MO_occ(homo)) lo = m;
	ASSERT_GE(lo, 0);
	for (int j = 0; j < wave.get_nex(); j++) {
		const double a = wave.get_MO_coef(lo, j), b = wave.get_MO_coef(homo, j);
		wave.set_MO_coef(lo, j, (a + b) / std::sqrt(2.0));
		wave.set_MO_coef(homo, j, (a - b) / std::sqrt(2.0));
	}
	EXPECT_FALSE(eli_family::alpha_beta_orbitals_differ(wave));
	EXPECT_EQ(eli_family::eli_variants_for(wave, &warn).size(), 2u) << warn;
}
