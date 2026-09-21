// Fast tests for the hkl generation, the tsc/tscb block and merge code and the
// Thakkar-style form factor helpers. Every input is written under the system temp
// directory and removed again; nothing touches the repository tree.
#include "pch.h"
#include <gtest/gtest.h>

#include "core/convenience.h"
#include "core/constants.h"
#include "core/scattering_factors.h"
#include "core/tsc_block.h"
#include "core/cell.h"
#include "core/wfn_class.h"
#include "core/atoms.h"
#include "core/cube.h"

#include <filesystem>
#include <fstream>
#include <sstream>

namespace
{
	namespace fs = std::filesystem;
	using block = tsc_block<int, cdouble>;

	fs::path tmp(const std::string& name)
	{
		return fs::temp_directory_path() / ("nosphera2_sfac_" + name);
	}

	//a fresh directory that becomes the working directory for functions writing into cwd
	//(kpts.dat, combined.tsc); restores the old cwd and removes the directory afterwards
	struct scoped_cwd
	{
		fs::path old;
		fs::path dir;
		explicit scoped_cwd(const std::string& name) : old(fs::current_path()), dir(tmp(name))
		{
			fs::remove_all(dir);
			fs::create_directories(dir);
			fs::current_path(dir);
		}
		~scoped_cwd()
		{
			fs::current_path(old);
			fs::remove_all(dir);
		}
	};

	std::string slurp(const fs::path& p)
	{
		std::ifstream in(p);
		std::stringstream ss;
		ss << in.rdbuf();
		return ss.str();
	}

	void write_p63_cif(const fs::path& cif)
	{
		std::ofstream out(cif);
		out << "data_test\n_cell_length_a 8.0\n_cell_length_b 8.0\n_cell_length_c 6.0\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 120\n_cell_volume 332.554\n"
			<< "loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n'-y, x-y, z'\n'-x+y, -x, z'\n'-x, -y, z+1/2'\n'y, -x+y, z+1/2'\n'x-y, x, z+1/2'\n";
	}

	//cubic P1 cell, a = 4 A, with an optional atom loop given as raw rows
	constexpr double kA = 4.0;
	void write_p1_cif(const fs::path& cif, const std::string& atom_rows = "", const std::string& extra_headers = "")
	{
		std::ofstream out(cif);
		out << "data_test\n_cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n_cell_volume 64.0\n"
			<< "loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n";
		if (!atom_rows.empty())
			out << "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n" << extra_headers << atom_rows;
	}

	std::array<i3, 6> p63_images(const i3& v)
	{
		const int h = v[0], k = v[1], l = v[2];
		return { i3{ h, k, l }, i3{ k, -h - k, l }, i3{ -h - k, h, l }, i3{ -h, -k, l }, i3{ -k, h + k, l }, i3{ h + k, -h, l } };
	}

	bool friedel_half(const i3& v)
	{
		return v[2] > 0 || (v[2] == 0 && (v[1] > 0 || (v[1] == 0 && v[0] > 0)));
	}

	i3 canonical(const i3& v)
	{
		return friedel_half(v) ? v : i3{ -v[0], -v[1], -v[2] };
	}

	hkl_list canonical_set(const hkl_list& in)
	{
		hkl_list out;
		for (const i3& v : in)
			out.insert(canonical(v));
		return out;
	}

	//independent count of the Friedel half of the resolution sphere: d*^2 from a textbook
	//formula, not from the cell's reciprocal metric, and the box wide enough for either cell
	template <typename F>
	hkl_list brute_sphere_of(const F& dstar2, const double longest_axis, const double dmin)
	{
		const double d_keep = dmin * (1.0 - 1e-3);
		const double s_max = 1.0 / (d_keep * d_keep);
		const int ext = int(longest_axis / d_keep) + 1;
		hkl_list out;
		for (int h = -ext; h <= ext; h++)
			for (int k = -ext; k <= ext; k++)
				for (int l = -ext; l <= ext; l++)
				{
					const i3 v{ h, k, l };
					if (dstar2(h, k, l) <= s_max && friedel_half(v))
						out.insert(v);
				}
		return out;
	}

	//cubic cell of edge a
	hkl_list brute_sphere(const double a, const double dmin)
	{
		return brute_sphere_of([a](int h, int k, int l) { return (h * h + k * k + l * l) / (a * a); }, a, dmin);
	}

	//hexagonal cell: d*^2 = 4 (h^2 + hk + k^2) / (3 a^2) + l^2 / c^2
	hkl_list brute_sphere_hex(const double a, const double c, const double dmin)
	{
		return brute_sphere_of([a, c](int h, int k, int l) { return 4.0 * (h * h + h * k + k * k) / (3.0 * a * a) + l * l / (c * c); }, std::max(a, c), dmin);
	}

	//spherical Bessel functions j_0..j_3 in closed form, for x > 0
	double sph_bessel(const int l, const double x)
	{
		const double s = std::sin(x), c = std::cos(x);
		switch (l)
		{
		case 0: return s / x;
		case 1: return s / (x * x) - c / x;
		case 2: return (3.0 / (x * x) - 1.0) * s / x - 3.0 * c / (x * x);
		default: return (15.0 / (x * x * x) - 6.0 / x) * s / x - (15.0 / (x * x) - 1.0) * c / x;
		}
	}

	//Simpson's rule over [0, R] with n (even) intervals; f(0) is passed in because the
	//integrands below vanish there while their closed forms divide by x
	template <typename F>
	double simpson(const F& f, const double f0, const double R, const int n)
	{
		const double dr = R / n;
		double sum = f0 + f(R);
		for (int i = 1; i < n; i++)
			sum += (i % 2 ? 4.0 : 2.0) * f(i * dr);
		return sum * dr / 3.0;
	}

	options make_options(std::vector<std::string> args)
	{
		args.insert(args.begin(), "NoSpherA2");
		std::vector<char*> argv;
		for (auto& a : args) argv.push_back(a.data());
		int argc = (int)args.size();
		options opt(argc, argv.data(), std::cout);
		opt.digest_options();
		return opt;
	}

	block two_scatterer_block(const std::string& header = "")
	{
		const hkl_list hkl{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 } };
		cvec2 sf{ { cdouble(1.0, 0.5), cdouble(2.0, -0.5) }, { cdouble(3.0, 1.0), cdouble(4.0, -1.0) } };
		return block(sf, svec{ "C1", "O1" }, hkl, header);
	}

	void write_le_int(std::ofstream& out, const std::int32_t v)
	{
		out.write(reinterpret_cast<const char*>(&v), sizeof(v));
	}
}

//the accuracy ladder shared with the XCW I tensor screening
TEST(ScatteringFactorTests, CutoffLadder)
{
	EXPECT_DOUBLE_EQ(cutoff(0), 1e-10);
	EXPECT_DOUBLE_EQ(cutoff(2), 1e-10);
	EXPECT_DOUBLE_EQ(cutoff(3), 1e-14);
	EXPECT_DOUBLE_EQ(cutoff(4), 1e-30);
	EXPECT_DOUBLE_EQ(cutoff(9), 1e-30);
}

//the Fourier-Bessel transform of a Gaussian radial part against the closed form
//H^l exp(-H^2 / 4b) / (2^l b^(l + 3/2)) for l = 0..3, and the closed form itself against
//a quadrature of its definition, 4 pi int r^(l+2) exp(-b r^2) j_l(H r) dr = pi^(3/2) FBI
TEST(ScatteringFactorTests, FourierBesselIntegralClosedForm)
{
	for (int l = 0; l < 4; l++)
	{
		const double b = 0.7 + 0.3 * l;
		const primitive p(1, l, b, 1.0);
		for (const double H : { 0.0, 0.5, 2.0 })
		{
			const double expected = std::pow(H, l) * std::exp(-H * H / (4.0 * b)) / (std::pow(2.0, l) * std::pow(b, l + 1.5));
			EXPECT_NEAR(fourier_bessel_integral(p, H, l), expected, 1e-12 * std::max(1.0, expected)) << "l=" << l << " H=" << H;
		}
		const double H = 2.0;
		const double quad = 4.0 * constants::PI * simpson([&](double r) { return std::pow(r, l + 2) * std::exp(-b * r * r) * sph_bessel(l, H * r); }, 0.0, 10.0, 4000);
		EXPECT_NEAR(quad, constants::PI3_2 * fourier_bessel_integral(p, H, l), 1e-8 * quad) << "l=" << l;
	}
}

//sfac_bessel carries the i^l phase: real, +imag, -real, -imag for l = 0..3, with the
//magnitude PI^(3/2) * FBI * normalised coefficient * Y_l along the k direction
TEST(ScatteringFactorTests, SfacBesselPhaseCycle)
{
	const double H = 1.3;
	const double k[5] = { 0.0, 0.0, 1.0, H, 0.0 };
	for (int l = 0; l < 4; l++)
	{
		const primitive p(1, l, 0.9, 0.8);
		vec coefs(2 * l + 1, 0.0);
		coefs[l] = 1.0;
		const double magnitude = constants::PI3_2 * fourier_bessel_integral(p, H, l) * p.get_normalized_coefficient() * constants::spherical_harmonic(l, k, coefs.data());
		const cdouble got = sfac_bessel(p, k, coefs.data());
		const cdouble expected = (l % 2 == 0) ? cdouble(magnitude, 0.0) : cdouble(0.0, magnitude);
		const double sign = (l >= 2) ? -1.0 : 1.0;
		EXPECT_NEAR(got.real(), sign * expected.real(), 1e-12) << "l=" << l;
		EXPECT_NEAR(got.imag(), sign * expected.imag(), 1e-12) << "l=" << l;
	}
	//l = 1 along z with the m = 0 coefficient is the analytic sqrt(3/4pi) * z
	const primitive p1(1, 1, 0.9, 0.8);
	const double c1[3] = { 0.0, 1.0, 0.0 };
	const double analytic = constants::PI3_2 * H * std::exp(-H * H / 3.6) / (2.0 * std::pow(0.9, 2.5)) * p1.get_normalized_coefficient() * std::sqrt(3.0 / (4.0 * constants::PI));
	EXPECT_NEAR(sfac_bessel(p1, k, c1).imag(), analytic, 1e-12);
	EXPECT_NEAR(sfac_bessel(p1, k, c1).real(), 0.0, 1e-15);
}

//primitive normalisation matches the textbook Gaussian norm for an s function:
//(2a/pi)^(3/4), so the derived getters agree with what sfac_bessel multiplies with
TEST(ScatteringFactorTests, PrimitiveNormalisation)
{
	const double a = 1.7;
	const primitive s(1, 0, a, 0.5);
	//the class holds the radial norm only: the angular 1/sqrt(4 pi) lives in the spherical harmonic,
	//so the s value is (2a/pi)^(3/4) * 2 sqrt(pi) = (128 a^3 / pi)^(1/4)
	EXPECT_NEAR(s.get_normalized_coefficient(), 0.5 * std::pow(2.0 * a / constants::PI, 0.75) * 2.0 * std::sqrt(constants::PI), 1e-12);
	EXPECT_NEAR(s.get_exp_l_plus_3_2(), std::pow(a, 1.5), 1e-12);
	EXPECT_EQ(s.get_type(), 0);
	EXPECT_DOUBLE_EQ(s.get_exp(), a);
	//for every l the radial norm N must give int (N r^l exp(-a r^2))^2 r^2 dr = 1 by quadrature
	for (int l = 0; l < 4; l++)
	{
		const primitive p(1, l, a, 1.0);
		const double N = p.get_normalized_coefficient();
		const double quad = simpson([&](double r) { return N * N * std::pow(r, 2 * l + 2) * std::exp(-2.0 * a * r * r); }, 0.0, 8.0, 4000);
		EXPECT_NEAR(quad, 1.0, 1e-9) << "l=" << l;
	}
}

//k points are the reciprocal cell times hkl: for a cubic cell |k(h00)| = 2 pi h / a in bohr^-1,
//and saving them then reading kpts.dat back returns the same vectors and hkl
TEST(ScatteringFactorTests, MakeKPtsAndKptsRoundTrip)
{
	const fs::path cif = tmp("kpts.cif");
	write_p1_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	hkl_list hkl{ i3{ 1, 0, 0 }, i3{ 0, 2, 0 }, i3{ 1, 1, 1 } };
	scoped_cwd cwd("kpts_dir");
	vec2 k;
	std::ostringstream log;
	make_k_pts(false, true, c, hkl, k, log, true);
	EXPECT_NE(log.str().find("K_point_vector is here! size: 3"), std::string::npos);
	ASSERT_EQ(k.size(), 3u);
	ASSERT_EQ(k[0].size(), 3u);
	const double g = constants::TWO_PI / constants::ang2bohr(kA);
	//set order: (0,2,0) < (1,0,0) < (1,1,1)
	EXPECT_NEAR(std::hypot(k[0][0], k[1][0], k[2][0]), 2.0 * g, 1e-10);
	EXPECT_NEAR(std::hypot(k[0][1], k[1][1], k[2][1]), g, 1e-10);
	EXPECT_NEAR(std::hypot(k[0][2], k[1][2], k[2][2]), std::sqrt(3.0) * g, 1e-10);
	EXPECT_NEAR(k[1][0], 2.0 * g, 1e-10);
	ASSERT_TRUE(fs::exists("kpts.dat"));
	vec2 k2;
	hkl_list hkl2;
	std::ostringstream log2;
	make_k_pts(true, false, c, hkl2, k2, log2, false);
	EXPECT_NE(log2.str().find("expecting 3 k points"), std::string::npos);
	EXPECT_EQ(hkl2, hkl);
	ASSERT_EQ(k2.size(), 3u);
	for (int i = 0; i < 3; i++)
		for (int n = 0; n < 3; n++)
			EXPECT_DOUBLE_EQ(k2[i][n], k[i][n]);
}

//a hexagonal cell exercises the off-diagonal reciprocal matrix: |k(100)| = 2 pi a* with
//a* = 2 / (sqrt(3) a), |k(001)| = 2 pi / c, k(100).k(010) = |k(100)|^2 cos 60, k(100) _|_ k(001)
TEST(ScatteringFactorTests, MakeKPtsHexagonalMetric)
{
	const fs::path cif = tmp("kpts_hex.cif");
	write_p63_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	//set order: (0,0,1), (0,1,0), (1,0,0)
	hkl_list hkl{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 }, i3{ 0, 0, 1 } };
	vec2 k;
	std::ostringstream log;
	make_k_pts(false, false, c, hkl, k, log, false);
	ASSERT_EQ(k.size(), 3u);
	ASSERT_EQ(k[0].size(), 3u);
	const double a_star = constants::TWO_PI * 2.0 / (std::sqrt(3.0) * constants::ang2bohr(8.0));
	const double c_star = constants::TWO_PI / constants::ang2bohr(6.0);
	auto dot = [&k](int i, int j) { return k[0][i] * k[0][j] + k[1][i] * k[1][j] + k[2][i] * k[2][j]; };
	EXPECT_NEAR(std::sqrt(dot(2, 2)), a_star, 1e-10);
	EXPECT_NEAR(std::sqrt(dot(1, 1)), a_star, 1e-10);
	EXPECT_NEAR(std::sqrt(dot(0, 0)), c_star, 1e-10);
	EXPECT_NEAR(dot(2, 1), 0.5 * a_star * a_star, 1e-10);
	EXPECT_NEAR(dot(2, 0), 0.0, 1e-10);
	EXPECT_NEAR(dot(1, 0), 0.0, 1e-10);
}

//the sphere generator of a cubic P1 cell agrees with a brute-force enumeration of the
//Friedel half for two resolutions, and the debug print reports the index extreme
TEST(ScatteringFactorTests, GenerateHklSphereMatchesBruteForce)
{
	const fs::path cif = tmp("sphere.cif");
	write_p1_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	for (const double dmin : { 2.0, 1.1 })
	{
		hkl_list hkl;
		std::ostringstream log;
		generate_hkl(dmin, hkl, {}, c, log, true);
		EXPECT_EQ(hkl, brute_sphere(kA, dmin)) << "dmin=" << dmin;
		EXPECT_NE(log.str().find("extreme: " + std::to_string(int(kA / (dmin * (1.0 - 1e-3)) + 1e-4))), std::string::npos);
		EXPECT_NE(log.str().find("Read 1 symmetry elements!"), std::string::npos);
	}
	//the hexagonal cell has a non-zero hk cross term in the reciprocal metric, which the
	//cubic case never touches; the Friedel half of the sphere is still the whole answer
	const fs::path hex_cif = tmp("sphere_hex.cif");
	write_p63_cif(hex_cif);
	cell hex(hex_cif, sink, false, true);
	fs::remove(hex_cif);
	hkl_list hkl;
	std::ostringstream log;
	generate_hkl(1.5, hkl, {}, hex, log, false);
	EXPECT_EQ(hkl, brute_sphere_hex(8.0, 6.0, 1.5));
	EXPECT_GT(hkl.size(), 100u);
}

//a merohedral twin law (k,h,l) of a P6_3 cell maps the sphere onto itself, so the twin
//branch must expand and reduce back to the untwinned half list modulo Friedel choice
TEST(ScatteringFactorTests, GenerateHklTwinLawInSphereLeavesListUnchanged)
{
	const fs::path cif = tmp("twin.cif");
	write_p63_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	hkl_list plain, twinned;
	std::ostringstream log;
	generate_hkl(1.5, plain, {}, c, log, false);
	const vec2 twin{ { 0, 1, 0, 1, 0, 0, 0, 0, 1 } };
	std::ostringstream dlog;
	generate_hkl(1.5, twinned, twin, c, dlog, true);
	EXPECT_EQ(canonical_set(twinned), canonical_set(plain));
	EXPECT_EQ(twinned.size(), plain.size());
	EXPECT_EQ(twinned.count(i3{ 0, 0, 0 }), 0u);
	const std::string s = dlog.str();
	EXPECT_NE(s.find("Number of reflections before twin: " + std::to_string(plain.size())), std::string::npos);
	EXPECT_NE(s.find("Number of reflections after twin: "), std::string::npos);
	EXPECT_NE(s.find("Number of reflections after sym gen: "), std::string::npos);
	EXPECT_NE(s.find("Read 6 symmetry elements!"), std::string::npos);
}

//a twin law that is not a symmetry operation of the cell, here (h+2k,k,l), has images
//outside the sphere: they must be kept and the list grows beyond the sphere. The matrix is
//deliberately asymmetric so a transposed read of the nine numbers, (h,2h+k,l), is caught
TEST(ScatteringFactorTests, GenerateHklTwinLawOutsideSphereGrowsList)
{
	const fs::path cif = tmp("twin2.cif");
	write_p1_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	hkl_list twinned;
	std::ostringstream log;
	const vec2 twin{ { 1, 2, 0, 0, 1, 0, 0, 0, 1 } };
	generate_hkl(2.0, twinned, twin, c, log, false);
	const hkl_list plain = brute_sphere(kA, 2.0);
	EXPECT_GT(twinned.size(), plain.size());
	for (const i3& v : plain)
		EXPECT_EQ(canonical_set(twinned).count(v), 1u);
	//(0,1,0) twinned is (2,1,0), outside the 2 A sphere of a 4 A cell; the transpose would
	//send (1,0,0) to (1,2,0) instead, equally outside, so exactly one of the two may appear
	EXPECT_EQ(twinned.count(i3{ 2, 1, 0 }) + twinned.count(i3{ -2, -1, 0 }), 1u);
	EXPECT_EQ(twinned.count(i3{ 1, 2, 0 }) + twinned.count(i3{ -1, -2, 0 }), 0u);
	for (const i3& v : twinned)
		EXPECT_EQ(twinned.count(i3{ -v[0], -v[1], -v[2] }), 0u) << "Friedel pair kept twice";
}

//the box-only generator for P1: h,k in [-1,1], l in [0,1], Friedel-reduced, 000 removed = 13,
//and the ED variant doubles each axis and warns = 62
TEST(ScatteringFactorTests, GenerateHklBoxCountsAndEdWarning)
{
	const fs::path cif = tmp("box.cif");
	write_p1_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	const ivec2 box{ { -1, 1 }, { -1, 1 }, { -1, 1 } };
	hkl_list plain, ed;
	std::ostringstream log, elog;
	generate_hkl(box, plain, {}, c, log, true, false);
	EXPECT_EQ(plain.size(), 13u);
	EXPECT_EQ(plain.count(i3{ 0, 0, 0 }), 0u);
	for (const i3& v : plain)
	{
		EXPECT_LE(std::abs(v[0]), 1);
		EXPECT_EQ(plain.count(i3{ -v[0], -v[1], -v[2] }), 0u);
	}
	EXPECT_NE(log.str().find("Generating hkl between [-1, 1] ; [-1, 1] ; [-1, 1]"), std::string::npos);
	EXPECT_NE(log.str().find("Number of reflections after sym gen: 18"), std::string::npos);
	EXPECT_EQ(log.str().find("Warning: an ED table"), std::string::npos);
	generate_hkl(box, ed, {}, c, elog, false, true);
	EXPECT_EQ(ed.size(), 62u);
	EXPECT_NE(elog.str().find("Warning: an ED table generated from an index box alone"), std::string::npos);
	EXPECT_NE(elog.str().find("[-2, 2] ; [-2, 2] ; [-2, 2]"), std::string::npos);
}

//box generation in P6_3 expands the box by the point group: the result is the orbit of
//the box rows with l >= 0, modulo Friedel, and a twin law that is a group operation adds nothing
TEST(ScatteringFactorTests, GenerateHklBoxExpandsByPointGroup)
{
	const fs::path cif = tmp("box63.cif");
	write_p63_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	const ivec2 box{ { 0, 2 }, { 0, 1 }, { 0, 1 } };
	hkl_list got, twinned;
	std::ostringstream log;
	generate_hkl(box, got, {}, c, log, false, false);
	generate_hkl(box, twinned, vec2{ { -1, 0, 0, 0, -1, 0, 0, 0, 1 } }, c, log, false, false);
	hkl_list expected;
	for (int h = -2; h <= 2; h++)
		for (int k = -1; k <= 1; k++)
			for (int l = 0; l <= 1; l++)
				for (const i3& g : p63_images({ h, k, l }))
					expected.insert(canonical(g));
	expected.erase(i3{ 0, 0, 0 });
	EXPECT_EQ(canonical_set(got), expected);
	EXPECT_EQ(got.size(), expected.size());
	EXPECT_EQ(canonical_set(twinned), expected);
}

//generate_hkl_from_options routes: dmin alone is the sphere, dmin plus box the box orbit in
//the sphere, box alone the box generator, ED the sphere at dmin/2, none the -hkl file
TEST(ScatteringFactorTests, GenerateHklFromOptionsRoutes)
{
	const fs::path cif = tmp("routes.cif");
	write_p1_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	std::ostringstream log;
	hkl_list sphere;
	generate_hkl_from_options(make_options({ "-dmin", "2.0" }), sphere, c, log);
	EXPECT_EQ(sphere, brute_sphere(kA, 2.0));
	EXPECT_EQ(sphere.size(), 16u);
	hkl_list boxed;
	generate_hkl_from_options(make_options({ "-dmin", "2.0", "-hkl_min_max", "0", "1", "0", "1", "0", "1" }), boxed, c, log);
	EXPECT_EQ(boxed.size(), 7u);
	for (const i3& v : boxed)
		for (int i = 0; i < 3; i++)
			EXPECT_TRUE(v[i] == 0 || v[i] == 1);
	EXPECT_NE(log.str().find("Keeping only the symmetry images of the index box"), std::string::npos);
	hkl_list box_only;
	generate_hkl_from_options(make_options({ "-hkl_min_max", "0", "1", "0", "1", "0", "1" }), box_only, c, log);
	EXPECT_EQ(box_only.size(), 13u);
	hkl_list ed;
	generate_hkl_from_options(make_options({ "-dmin", "2.0", "-ED" }), ed, c, log);
	EXPECT_EQ(ed, brute_sphere(kA, 0.999));
	const fs::path hklfile = tmp("routes.hkl");
	{
		std::ofstream out(hklfile);
		out << "   1   0   0   10.00    1.00\n   0   1   0   10.00    1.00\n";
	}
	hkl_list read;
	generate_hkl_from_options(make_options({ "-hkl", hklfile.string() }), read, c, log);
	fs::remove(hklfile);
	EXPECT_EQ(read, (hkl_list{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 } }));
}

//read_hkl skips text and short lines, drops 000, applies the twin law and reports counts
TEST(ScatteringFactorTests, ReadHklSkipsTextDropsOriginAppliesTwin)
{
	const fs::path cif = tmp("readhkl.cif");
	write_p63_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	const fs::path hklfile = tmp("read.hkl");
	{
		std::ofstream out(hklfile);
		out << "CELL 0.71073 8 8 6 90 90 120\n\n   1   0   0   10.00    1.00\n   0   1   0   20.00    2.00\n   0   0   0    0.00    0.00\n   2  -1   1    5.00    0.50\n";
	}
	hkl_list plain;
	std::ostringstream log;
	read_hkl(hklfile, plain, {}, c, log, true);
	EXPECT_EQ(plain, (hkl_list{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 }, i3{ 2, -1, 1 } }));
	EXPECT_NE(log.str().find("popping back 0 0 0"), std::string::npos);
	EXPECT_NE(log.str().find("Nr of reflections read from file: 3"), std::string::npos);
	EXPECT_NE(log.str().find("Read 6 symmetry elements!"), std::string::npos);
	//the twin loop inserts into the set it walks, so the law has to be an involution here;
	//(h,h-k,l) is one and its transpose (h+k,-k,l) is a different one, which tells the two apart
	hkl_list twinned;
	std::ostringstream tlog;
	read_hkl(hklfile, twinned, vec2{ { 1, 0, 0, 1, -1, 0, 0, 0, 1 } }, c, tlog, false);
	fs::remove(hklfile);
	EXPECT_EQ(twinned, (hkl_list{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 }, i3{ 2, -1, 1 }, i3{ 1, 1, 0 }, i3{ 0, -1, 0 }, i3{ 2, 3, 1 } }));
	EXPECT_NE(tlog.str().find("Number of symmetry operations: 6"), std::string::npos);
	EXPECT_NE(tlog.str().find("Nr of reflections to be used: 6"), std::string::npos);
}

//read_hkl_full parses F^2 and sigma, keeps the sign of a negative F^2 in F, derives
//sigma(F) = sigma(F^2) / 2|F|, and returns the symmetry-expanded list
TEST(ScatteringFactorTests, ReadHklFullObservationsAndSymmetryExpansion)
{
	const fs::path cif = tmp("readfull.cif");
	write_p63_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	const fs::path hklfile = tmp("full.hkl");
	{
		std::ofstream out(hklfile);
		out << "   1   0   0  100.00    2.00\n   0   0   1 -100.00    4.00\n   0   0   0    0.00    0.00\n";
	}
	hkl_list hkl;
	std::vector<scattering_data> obs;
	std::ostringstream log;
	const hkl_list enlarged = read_hkl_full(hklfile, hkl, {}, c, log, obs, true);
	fs::remove(hklfile);
	ASSERT_EQ(obs.size(), 3u);
	EXPECT_DOUBLE_EQ(obs[0].F_obs2, 100.0);
	EXPECT_DOUBLE_EQ(obs[0].abs_F_obs, 10.0);
	EXPECT_DOUBLE_EQ(obs[0].F_obs, 10.0);
	EXPECT_DOUBLE_EQ(obs[0].sigma_obs2, 2.0);
	EXPECT_DOUBLE_EQ(obs[0].sigma_obs, 0.1);
	EXPECT_DOUBLE_EQ(obs[1].F_obs2, 100.0);
	EXPECT_DOUBLE_EQ(obs[1].F_obs, -10.0);
	EXPECT_DOUBLE_EQ(obs[1].abs_F_obs, 10.0);
	EXPECT_DOUBLE_EQ(obs[1].sigma_obs, 0.2);
	EXPECT_EQ(hkl, (hkl_list{ i3{ 1, 0, 0 }, i3{ 0, 0, 1 } }));
	//the sixfold axis turns (1,0,0) into the six in-plane vectors; (0,0,1) is invariant
	hkl_list expected{ i3{ 0, 0, 1 } };
	for (const i3& g : p63_images({ 1, 0, 0 }))
		expected.insert(g);
	EXPECT_EQ(enlarged, expected);
	EXPECT_NE(log.str().find("popping back 0 0 0"), std::string::npos);
}

//an F^2 that fills its eight-character field (3I4,2F8.2) must not bleed into l
TEST(ScatteringFactorTests, ReadHklFullStrongReflectionKeepsL)
{
	const fs::path cif = tmp("strong.cif");
	write_p1_cif(cif);
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	fs::remove(cif);
	const fs::path hklfile = tmp("strong.hkl");
	{
		std::ofstream out(hklfile);
		out << "   1   0   012345.67   10.00\n";
	}
	hkl_list hkl;
	std::vector<scattering_data> obs;
	std::ostringstream log;
	read_hkl_full(hklfile, hkl, {}, c, log, obs, false);
	fs::remove(hklfile);
	EXPECT_EQ(hkl, (hkl_list{ i3{ 1, 0, 0 } }));
	ASSERT_EQ(obs.size(), 1u);
	EXPECT_NEAR(obs[0].F_obs2, 12345.67, 1e-2);
}

//the WFN CIF reader matches rows to wavefunction atoms by position, stores fractional
//coordinates and the CIF-derived atomID, and builds the type lists
TEST(ScatteringFactorCifTests, ReadAtomsMatchesWfnAtoms)
{
	const fs::path cif = tmp("atoms.cif");
	//distinct coordinates per axis so a swapped axis or an atom matched to the wrong row shows
	write_p1_cif(cif, "C1 C 0.25 0.35 0.45\nO1 O 0.5 0.6 0.7\n");
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	WFN w(e_origin::NOT_YET_DEFINED);
	const double f = constants::ang2bohr(kA);
	w.push_back_atom("C1", 0.25 * f, 0.35 * f, 0.45 * f, 6);
	w.push_back_atom("O1", 0.5 * f, 0.6 * f, 0.7 * f, 8);
	std::ifstream in(cif);
	ivec types, asym_to_type, asym;
	bvec needs_grid(2, false);
	std::ostringstream log;
	const svec labels = read_atoms_from_CIF(in, {}, c, w, {}, types, asym_to_type, asym, needs_grid, log, true);
	in.close();
	fs::remove(cif);
	EXPECT_EQ(labels, (svec{ "C1", "O1" }));
	EXPECT_EQ(asym, (ivec{ 0, 1 }));
	EXPECT_EQ(types, (ivec{ 6, 8 }));
	EXPECT_EQ(asym_to_type, (ivec{ 0, 1 }));
	EXPECT_TRUE(needs_grid[0] && needs_grid[1]);
	EXPECT_DOUBLE_EQ(w.get_atom(0).get_frac_coordinate(0), 0.25);
	EXPECT_DOUBLE_EQ(w.get_atom(0).get_frac_coordinate(1), 0.35);
	EXPECT_DOUBLE_EQ(w.get_atom(0).get_frac_coordinate(2), 0.45);
	EXPECT_DOUBLE_EQ(w.get_atom(1).get_frac_coordinate(0), 0.5);
	EXPECT_DOUBLE_EQ(w.get_atom(1).get_frac_coordinate(1), 0.6);
	EXPECT_DOUBLE_EQ(w.get_atom(1).get_frac_coordinate(2), 0.7);
	EXPECT_EQ(w.get_id_for_atom(1).to_hex_string(), atomID(0.5, 0.6, 0.7, 0, 8).to_hex_string());
	EXPECT_NE(w.get_id_for_atom(1).to_hex_string(), atomID(0.5, 0.7, 0.6, 0, 8).to_hex_string());
	EXPECT_NE(log.str().find("ASYM:"), std::string::npos);
	EXPECT_NE(log.str().find("There are 2 types of atoms"), std::string::npos);
}

//a loop_ directly after the atom rows (no blank line) ends the atom loop instead of being read as a row
TEST(ScatteringFactorCifTests, ReadAtomsStopsAtNextLoop)
{
	const fs::path cif = tmp("nextloop.cif");
	write_p1_cif(cif, "C1 C 0.25 0.35 0.45\nloop_\n_atom_site_aniso_label\n_atom_site_aniso_U_11\nC1 0.01\n");
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	WFN w(e_origin::NOT_YET_DEFINED);
	const double f = constants::ang2bohr(kA);
	w.push_back_atom("C1", 0.25 * f, 0.35 * f, 0.45 * f, 6);
	std::ifstream in(cif);
	ivec types, asym_to_type, asym;
	bvec needs_grid(1, false);
	std::ostringstream log;
	const svec labels = read_atoms_from_CIF(in, {}, c, w, {}, types, asym_to_type, asym, needs_grid, log, false);
	in.close();
	fs::remove(cif);
	EXPECT_EQ(labels, (svec{ "C1" }));
	EXPECT_EQ(types, (ivec{ 6 }));
}

//disorder groups: only rows of the requested parts are matched, so the same site labelled
//O1A (part 1) and O1B (part 2) resolves to whichever part was asked for
TEST(ScatteringFactorCifTests, ReadAtomsFiltersDisorderGroups)
{
	const fs::path cif = tmp("parts.cif");
	write_p1_cif(cif, "C1 C 0.25 0.25 0.25 .\nO1A O 0.5 0.5 0.5 1\nO1B O 0.5 0.5 0.5 2\n", "_atom_site_disorder_group\n");
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	const double f = constants::ang2bohr(kA);
	for (const int part : { 1, 2 })
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C1", 0.25 * f, 0.25 * f, 0.25 * f, 6);
		w.push_back_atom("O1", 0.5 * f, 0.5 * f, 0.5 * f, 8);
		std::ifstream in(cif);
		ivec types, asym_to_type, asym;
		bvec needs_grid(2, false);
		std::ostringstream log;
		const svec labels = read_atoms_from_CIF(in, ivec{ 0, part }, c, w, {}, types, asym_to_type, asym, needs_grid, log, true);
		EXPECT_EQ(labels, (svec{ "C1", part == 1 ? "O1A" : "O1B" }));
		EXPECT_EQ(w.get_id_for_atom(1).to_hex_string(), atomID(0.5, 0.5, 0.5, part, 8).to_hex_string());
		EXPECT_NE(log.str().find("Wrong part!"), std::string::npos);
	}
	fs::remove(cif);
}

//known atoms are skipped by label or by hex ID, quoted labels and rows wrapped over two
//lines are read as one row, and a hydrogen labelled D is accepted
TEST(ScatteringFactorCifTests, ReadAtomsKnownAtomsQuotedAndWrappedRows)
{
	const fs::path cif = tmp("known.cif");
	write_p1_cif(cif, "'C1' C 0.25 0.25\n 0.25\nO1 O 0.5 0.5 0.5\nD1 H 0.75 0.75 0.75\n");
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	const double f = constants::ang2bohr(kA);
	const std::string o_id = atomID(0.5, 0.5, 0.5, 0, 8).to_hex_string();
	for (const svec& known : { svec{ "C1" }, svec{ o_id } })
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C1", 0.25 * f, 0.25 * f, 0.25 * f, 6);
		w.push_back_atom("O1", 0.5 * f, 0.5 * f, 0.5 * f, 8);
		w.push_back_atom("H1", 0.75 * f, 0.75 * f, 0.75 * f, 1);
		std::ifstream in(cif);
		ivec types, asym_to_type, asym;
		bvec needs_grid(3, false);
		std::ostringstream log;
		const svec labels = read_atoms_from_CIF(in, {}, c, w, known, types, asym_to_type, asym, needs_grid, log, true);
		EXPECT_EQ(labels, (svec{ known[0] == "C1" ? "O1" : "C1", "D1" }));
		EXPECT_NE(log.str().find("I already know this one!"), std::string::npos);
		//the skipped atom still gets a type so its spherical density can be built
		EXPECT_EQ(types.size(), 3u);
	}
	fs::remove(cif);
}

//a row whose wavefunction atom is missing is reported with its tolerances and left out
TEST(ScatteringFactorCifTests, ReadAtomsReportsUnmatchedRow)
{
	const fs::path cif = tmp("unmatched.cif");
	write_p1_cif(cif, "C1 C 0.25 0.25 0.25\nN1 N 0.1 0.9 0.4\n");
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	const double f = constants::ang2bohr(kA);
	WFN w(e_origin::NOT_YET_DEFINED);
	w.push_back_atom("C1", 0.25 * f, 0.25 * f, 0.25 * f, 6);
	std::ifstream in(cif);
	ivec types, asym_to_type, asym;
	bvec needs_grid(1, false);
	std::ostringstream log;
	const svec labels = read_atoms_from_CIF(in, {}, c, w, {}, types, asym_to_type, asym, needs_grid, log, true);
	in.close();
	fs::remove(cif);
	EXPECT_EQ(labels, (svec{ "C1" }));
	EXPECT_NE(log.str().find("I did not find this atom! Tolerances were:"), std::string::npos);
	EXPECT_EQ(types, (ivec{ 6 }));
}

//a label without the element symbol is rescued by the type field
TEST(ScatteringFactorCifTests, TypeSymbolRescuesLabelWithoutElement)
{
	const fs::path cif = tmp("q.cif");
	write_p1_cif(cif, "Q1 C 0.25 0.25 0.25\n");
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	const double f = constants::ang2bohr(kA);
	WFN w(e_origin::NOT_YET_DEFINED);
	w.push_back_atom("C1", 0.25 * f, 0.25 * f, 0.25 * f, 6);
	std::ifstream in(cif);
	ivec types, asym_to_type, asym;
	bvec needs_grid(1, false);
	std::ostringstream log;
	const svec labels = read_atoms_from_CIF(in, {}, c, w, {}, types, asym_to_type, asym, needs_grid, log, false, true);
	in.close();
	fs::remove(cif);
	EXPECT_EQ(labels, (svec{ "Q1" }));
}

//the wavefunction-free reader for XCW: types are Z (get_Z_from_label is 0-based, plus one), positions in bohr from the cell,
//a short row is skipped, and the debug trace names the columns
TEST(ScatteringFactorCifTests, ReadAtomsWithoutWfnParsesTypesAndPositions)
{
	const fs::path cif = tmp("xcw.cif");
	write_p1_cif(cif, "C1 C 0.25 0.25 0.25\nN9 N 0.1\nO1 O 0.5 0.6 0.7\n_extra_item 1\n");
	std::ostringstream sink;
	cell c(cif, sink, false, true);
	std::ifstream in(cif);
	int ncen = -1;
	bvec needs_grid;
	std::vector<asym_atom> atoms;
	testing::internal::CaptureStdout();
	read_atoms_from_CIF(in, c, ncen, needs_grid, atoms, true);
	const std::string out = testing::internal::GetCapturedStdout();
	in.close();
	fs::remove(cif);
	ASSERT_EQ(ncen, 2);
	ASSERT_EQ(atoms.size(), 2u);
	EXPECT_EQ(needs_grid, (bvec{ true, true }));
	EXPECT_EQ(atoms[0].label, "C1");
	EXPECT_EQ(atoms[0].type, 6);
	EXPECT_EQ(atoms[1].label, "O1");
	EXPECT_EQ(atoms[1].type, 8);
	EXPECT_EQ(atoms[1].frac_pos, (d3{ 0.5, 0.6, 0.7 }));
	for (int i = 0; i < 3; i++)
		EXPECT_NEAR(atoms[1].pos[i], atoms[1].frac_pos[i] * constants::ang2bohr(kA), 1e-9) << "axis " << i;
	EXPECT_NE(out.find("Found atom_site loop"), std::string::npos);
	EXPECT_NE(out.find("Skipping malformed atom line:\nN9 N 0.1"), std::string::npos);
	EXPECT_NE(out.find("Total atoms parsed: 2"), std::string::npos);
}

//a bad stream is refused before anything is parsed
TEST(ScatteringFactorCifTests, ReadAtomsWithoutWfnRejectsBadStream)
{
	const fs::path cif = tmp("nonexistent_xcw.cif");
	fs::remove(cif);
	std::ifstream in(cif);
	cell c;
	int ncen = 0;
	bvec needs_grid;
	std::vector<asym_atom> atoms;
	EXPECT_THROW(read_atoms_from_CIF(in, c, ncen, needs_grid, atoms, false), std::runtime_error);
}

//a block built from an hkl_list orders reflections like the set and the logged accessors
//return the same rows as the plain ones
TEST(ScatteringFactorBlockTests, HklListCtorAndAccessors)
{
	const block b = two_scatterer_block();
	EXPECT_EQ(b.reflection_size(), 2u);
	EXPECT_EQ(b.scatterer_size(), 2u);
	EXPECT_FALSE(b.is_empty());
	EXPECT_EQ(b.get_indices(0), (std::array<int, 3>{ 0, 1, 0 }));
	EXPECT_EQ(b.get_index(0, 1), 1);
	std::ostringstream log;
	//the rows are stored in the set's reflection order, so sf[1] is O1 as given
	EXPECT_EQ(b.get_sf_for_scatterer(1), (cvec{ cdouble(3.0, 1.0), cdouble(4.0, -1.0) }));
	EXPECT_EQ(b.get_sf_for_scatterer(1, log), b.get_sf_for_scatterer(1));
	EXPECT_EQ(std::get<std::string>(b.get_scatterer(1, log)), "O1");
	EXPECT_EQ(b.get_scatterers_as<std::string>(), (svec{ "C1", "O1" }));
	EXPECT_EQ(b.get_scatterers_string(), (svec{ "C1", "O1" }));
	EXPECT_EQ(b.get_index_vector()[2], (ivec{ 0, 0 }));
	EXPECT_FALSE(b.get_AD());
}

//an empty block reports zero sizes, and mismatched sf rows give a zero reflection count
TEST(ScatteringFactorBlockTests, EmptyAndInconsistentSizes)
{
	const block empty;
	EXPECT_TRUE(empty.is_empty());
	EXPECT_EQ(empty.reflection_size(), 0u);
	EXPECT_EQ(empty.scatterer_size(), 0u);
	const hkl_list hkl{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 } };
	cvec2 sf{ { cdouble(1.0, 0.0), cdouble(2.0, 0.0) }, { cdouble(3.0, 0.0) } };
	const block ragged(sf, svec{ "C1", "O1" }, hkl);
	EXPECT_EQ(ragged.reflection_size(), 0u);
	EXPECT_THROW(ragged.write_tsc_file("x.cif", tmp("ragged.tsc")), std::runtime_error);
	fs::remove(tmp("ragged.tsc"));
	cvec2 one{ { cdouble(1.0, 0.0), cdouble(2.0, 0.0) } };
	const block short_labels(one, svec{ "C1", "O1" }, hkl);
	EXPECT_EQ(short_labels.scatterer_size(), 0u);
	EXPECT_THROW(short_labels.write_tscb_file("x.cif", tmp("short.tscb")), std::runtime_error);
	fs::remove(tmp("short.tscb"));
}

//atomID scatterers are written as SCATTERER_IDS hex tokens, AD as a header line, and the
//text writer formats one scientific complex pair per scatterer
TEST(ScatteringFactorBlockTests, TextWriterIdsAndAd)
{
	const hkl_list hkl{ i3{ 1, 0, 0 } };
	cvec2 sf{ { cdouble(1.5, -0.25) } };
	const atomID id(0.1, 0.2, 0.3, 0, 7);
	block b(sf, std::vector<atomID>{ id }, hkl);
	b.set_AD(true);
	EXPECT_TRUE(b.get_AD());
	EXPECT_EQ(b.get_scatterers_string(), (svec{ id.to_hex_string() }));
	const fs::path out = tmp("ids.tsc");
	b.write_tsc_file("some_dir/mycif.cif", out);
	const std::string text = slurp(out);
	fs::remove(out);
	EXPECT_NE(text.find("TITLE: mycif\nSYMM: expanded\nAD: TRUE\nSCATTERER_IDS: " + id.to_hex_string() + "\nDATA:\n"), std::string::npos);
	EXPECT_NE(text.find("1 0 0 1.50000000e+00,-2.50000000e-01 "), std::string::npos);
}

//the non-integer writer on a double-indexed block prints indices with three decimals
TEST(ScatteringFactorBlockTests, NonIntegerIndicesWriter)
{
	const hkl_list_d hkl{ d3{ 0.5, 0.0, 1.25 } };
	cvec2 sf{ { cdouble(2.0, 0.0) } };
	const tsc_block<double, cdouble> b(sf, svec{ "C1" }, hkl);
	EXPECT_DOUBLE_EQ(b.get_index(2, 0), 1.25);
	const fs::path out = tmp("frac.tsc");
	b.write_tsc_file_non_integer("frac.cif", out);
	const std::string text = slurp(out);
	fs::remove(out);
	EXPECT_NE(text.find("SCATTERERS: C1\n"), std::string::npos);
	EXPECT_NE(text.find("0.500 0.000 1.250 2.00000000e+00,0.00000000e+00"), std::string::npos);
}

//the CIF block names the partition and lists the reflections with real and imaginary parts
TEST(ScatteringFactorBlockTests, CifBlockNamesPartition)
{
	const block b = two_scatterer_block();
	options opt;
	opt.partition_type = PartitionType::Becke;
	const std::string text = b.get_tsc_cif_block(opt);
	EXPECT_NE(text.find("_aspheric_ffs_partitioning.name     'Becke'"), std::string::npos);
	EXPECT_NE(text.find("0 1 0 '[1 3 ]' '[0.5 1 ]'"), std::string::npos);
	EXPECT_NE(text.find("1 0 0 '[2 4 ]' '[-0.5 -1 ]'"), std::string::npos);
	opt.partition_type = PartitionType::EMBIS;
	EXPECT_NE(b.get_tsc_cif_block(opt).find("'EMBIS'"), std::string::npos);
	opt.partition_type = PartitionType::RI;
	EXPECT_NE(b.get_tsc_cif_block(opt).find("'RI-Fit'"), std::string::npos);
}

//append: into an empty block it adopts the argument, otherwise it adds new scatterers,
//warns about duplicates by name and refuses a different reflection list
TEST(ScatteringFactorBlockTests, AppendAdoptsAddsAndWarns)
{
	block target;
	target.append(two_scatterer_block(), std::cout);
	EXPECT_EQ(target.scatterer_size(), 2u);
	const hkl_list hkl{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 } };
	cvec2 sf{ { cdouble(5.0, 0.0), cdouble(6.0, 0.0) }, { cdouble(7.0, 0.0), cdouble(8.0, 0.0) } };
	block more(sf, svec{ "O1", "N1" }, hkl);
	testing::internal::CaptureStdout();
	target.append(std::move(more), std::cout);
	const std::string out = testing::internal::GetCapturedStdout();
	EXPECT_NE(out.find("Warning: Duplicate scatterer in append: O1"), std::string::npos);
	EXPECT_EQ(target.get_scatterers_as<std::string>(), (svec{ "C1", "O1", "N1" }));
	EXPECT_EQ(target.get_sf_for_scatterer(2)[1], cdouble(8.0, 0.0));
	EXPECT_EQ(target.get_sf_for_scatterer(1)[0], cdouble(3.0, 1.0));
}

//mixing string labels and IDs in one block is caught when the block is written
TEST(ScatteringFactorBlockTests, MixedLabelTypesRejected)
{
	block target = two_scatterer_block();
	const hkl_list hkl{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 } };
	cvec2 sf{ { cdouble(5.0, 0.0), cdouble(6.0, 0.0) } };
	block ids(sf, std::vector<atomID>{ atomID(0.1, 0.1, 0.1, 0, 2) }, hkl);
	target.append(ids, std::cout);
	EXPECT_EQ(target.scatterer_size(), 3u);
	EXPECT_THROW(target.write_tscb_file("x.cif", tmp("mixed.tscb")), std::runtime_error);
	EXPECT_THROW(target.write_tsc_file("x.cif", tmp("mixed.tsc")), std::runtime_error);
	fs::remove(tmp("mixed.tscb"));
	fs::remove(tmp("mixed.tsc"));
}

//a tscb round trip with string labels restores labels, indices and values, and one
//trailing byte after the last reflection is refused as an incompatible layout
TEST(ScatteringFactorBlockTests, TscbRoundTripAndTrailingBytes)
{
	const block b = two_scatterer_block("TITLE test\nSCATTERER_IDS\nSYMM: expanded");
	const fs::path out = tmp("round.tscb");
	b.write_tscb_file("x.cif", out);
	const block back(out);
	EXPECT_EQ(back.get_scatterers_as<std::string>(), (svec{ "C1", "O1" }));
	EXPECT_EQ(back.get_index_vector(), b.get_index_vector());
	for (std::size_t s = 0; s < 2; s++)
		EXPECT_EQ(back.get_sf_for_scatterer(s), b.get_sf_for_scatterer(s));
	{
		std::ofstream app(out, std::ios::binary | std::ios::app);
		app.put('x');
	}
	EXPECT_THROW(block trailing(out), std::runtime_error);
	fs::remove(out);
}

//hand-written tscb prefixes with negative sizes are refused before anything is allocated
TEST(ScatteringFactorBlockTests, TscbNegativeSizesRejected)
{
	const fs::path out = tmp("bad.tscb");
	{
		std::ofstream o(out, std::ios::binary);
		write_le_int(o, -1);
	}
	EXPECT_THROW(block b(out), std::runtime_error);
	{
		std::ofstream o(out, std::ios::binary);
		write_le_int(o, 0);
		write_le_int(o, -5);
	}
	EXPECT_THROW(block b(out), std::runtime_error);
	{
		std::ofstream o(out, std::ios::binary);
		write_le_int(o, 0);
		write_le_int(o, 2);
		o << "C1";
		write_le_int(o, -3);
	}
	EXPECT_THROW(block b(out), std::runtime_error);
	fs::remove(out);
	EXPECT_THROW(block b(tmp("missing.tscb")), std::runtime_error);
}

//checked merge: a tscb input is converted beside itself, the first file fixes the
//reflection list (000 dropped), a Friedel mate in the second file is conjugated, an
//unmatched reflection is skipped and a new scatterer is appended with zeros elsewhere
TEST(ScatteringFactorMergeTests, CheckedMergeConjugatesFriedelMates)
{
	scoped_cwd cwd("merge_checked");
	const hkl_list hkl_a{ i3{ 1, 0, 0 }, i3{ 0, 1, 0 }, i3{ 0, 0, 0 } };
	//the block orders reflections like the set: (0,0,0), (0,1,0), (1,0,0)
	cvec2 sf_a{ { cdouble(9.0, 9.0), cdouble(2.0, 2.0), cdouble(1.0, 1.0) }, { cdouble(9.0, 9.0), cdouble(4.0, 4.0), cdouble(3.0, 3.0) } };
	block a(sf_a, svec{ "C1", "O1" }, hkl_a, "TITLE a\nSYMM: expanded");
	a.write_tscb_file("a.cif", "a.tscb");
	{
		std::ofstream out("b.tsc");
		out << "TITLE: b\nSYMM: expanded\nSCATTERERS: O1 N1\nDATA:\n-1 0 0 5.0,0.5 6.0,0.6\n2 2 2 7.0,0.7 8.0,0.8\n";
	}
	testing::internal::CaptureStdout();
	const bool ok = merge_tscs("merge", pathvec{ "a.tscb", "b.tsc" }, false);
	const std::string out = testing::internal::GetCapturedStdout();
	ASSERT_TRUE(ok);
	EXPECT_NE(out.find("Converting to:"), std::string::npos);
	EXPECT_TRUE(fs::exists("a.tsc"));
	EXPECT_NE(out.find("Read 2 atoms, 1 are new."), std::string::npos);
	ASSERT_TRUE(fs::exists("combined.tscb"));
	const block combined("combined.tscb");
	EXPECT_EQ(combined.get_scatterers_as<std::string>(), (svec{ "C1", "O1", "N1" }));
	ASSERT_EQ(combined.reflection_size(), 2u);
	//set order of the first file: (0,1,0) then (1,0,0); 000 was dropped
	EXPECT_EQ(combined.get_indices(0), (std::array<int, 3>{ 0, 1, 0 }));
	EXPECT_EQ(combined.get_indices(1), (std::array<int, 3>{ 1, 0, 0 }));
	EXPECT_EQ(combined.get_sf_for_scatterer(0)[1], cdouble(1.0, 1.0));
	EXPECT_EQ(combined.get_sf_for_scatterer(1)[1], cdouble(5.0, -0.5));
	EXPECT_EQ(combined.get_sf_for_scatterer(2)[1], cdouble(6.0, -0.6));
	EXPECT_EQ(combined.get_sf_for_scatterer(1)[0], cdouble(4.0, 4.0));
	EXPECT_EQ(combined.get_sf_for_scatterer(2)[0], cdouble(0.0, 0.0));
}

//the text output path keeps the first file's header and lists all scatterers
TEST(ScatteringFactorMergeTests, CheckedMergeWritesTextTable)
{
	scoped_cwd cwd("merge_text");
	{
		std::ofstream out("a.tsc");
		out << "TITLE: a\nSYMM: expanded\nAD: TRUE\nSCATTERERS: C1\nDATA:\n1 0 0 1.0,0.0\n1 0 0 9.0,9.0\n";
		std::ofstream out2("b.tsc");
		out2 << "TITLE: b\nSCATTERERS: N1\nDATA:\n1 0 0 2.0,0.0\n";
	}
	ASSERT_TRUE(merge_tscs("merge", pathvec{ "a.tsc", "b.tsc" }, true));
	const std::string text = slurp("combined.tsc");
	EXPECT_NE(text.find("TITLE: combined\n"), std::string::npos);
	EXPECT_NE(text.find("SCATTERERS: C1 N1\nDATA:\n"), std::string::npos);
	//the duplicate row of the first file is ignored, the first value wins
	EXPECT_NE(text.find("1 0 0 1.00000000e+00,0.00000000e+00 2.00000000e+00,0.00000000e+00 "), std::string::npos);
	EXPECT_EQ(text.find("9.00000000e+00"), std::string::npos);
}

//unchecked merge trusts row order: equal counts merge by position, more or fewer rows fail
TEST(ScatteringFactorMergeTests, UncheckedMergeRequiresEqualRowCounts)
{
	scoped_cwd cwd("merge_unchecked");
	{
		std::ofstream out("a.tsc");
		out << "TITLE: a\nSCATTERERS: C1\nDATA:\n1 0 0 1.0,0.0\n0 1 0 2.0,0.0\n";
		std::ofstream out2("b.tsc");
		out2 << "TITLE: b\nSCATTERERS: N1\nDATA:\n0 1 0 3.0,0.0\n1 0 0 4.0,0.0\n";
		std::ofstream out3("c.tsc");
		out3 << "TITLE: c\nSCATTERERS: N1\nDATA:\n1 0 0 4.0,0.0\n";
		std::ofstream out4("d.tsc");
		out4 << "TITLE: d\nSCATTERERS: N1\nDATA:\n1 0 0 4.0,0.0\n0 1 0 4.0,0.0\n2 0 0 4.0,0.0\n";
	}
	ASSERT_TRUE(merge_tscs_without_checks("merge", pathvec{ "a.tsc", "b.tsc" }, false));
	const block combined("combined.tscb");
	ASSERT_EQ(combined.reflection_size(), 2u);
	EXPECT_EQ(combined.get_indices(0), (std::array<int, 3>{ 1, 0, 0 }));
	EXPECT_EQ(combined.get_sf_for_scatterer(1)[0], cdouble(3.0, 0.0));
	testing::internal::CaptureStderr();
	EXPECT_FALSE(merge_tscs_without_checks("merge", pathvec{ "a.tsc", "c.tsc" }, false));
	EXPECT_FALSE(merge_tscs_without_checks("merge", pathvec{ "a.tsc", "d.tsc" }, false));
	const std::string err = testing::internal::GetCapturedStderr();
	EXPECT_NE(err.find("Could not merge TSC files: A merged TSC file contains a different number of reflections"), std::string::npos);
	EXPECT_NE(err.find("Could not merge TSC files: A merged TSC file contains more reflections than the first file"), std::string::npos);
}

//every malformed input is reported on stderr and answered with false instead of a throw
TEST(ScatteringFactorMergeTests, MalformedInputsReturnFalse)
{
	scoped_cwd cwd("merge_bad");
	{
		std::ofstream("nodata.tsc") << "TITLE: x\nSCATTERERS: C1\n1 0 0 1.0,0.0\n";
		std::ofstream("noscat.tsc") << "TITLE: x\nDATA:\n1 0 0 1.0,0.0\n";
		std::ofstream("badvalue.tsc") << "TITLE: x\nSCATTERERS: C1\nDATA:\n1 0 0 1.0\n";
		std::ofstream("badcount.tsc") << "TITLE: x\nSCATTERERS: C1 O1\nDATA:\n1 0 0 1.0,0.0\n";
	}
	EXPECT_FALSE(merge_tscs("", pathvec{ "nodata.tsc" }, false));
	EXPECT_FALSE(merge_tscs("merge", pathvec{}, false));
	testing::internal::CaptureStderr();
	EXPECT_FALSE(merge_tscs("merge", pathvec{ "nodata.tsc" }, false));
	EXPECT_FALSE(merge_tscs("merge", pathvec{ "noscat.tsc" }, false));
	EXPECT_FALSE(merge_tscs("merge", pathvec{ "badvalue.tsc" }, false));
	EXPECT_FALSE(merge_tscs("merge", pathvec{ "badcount.tsc" }, false));
	EXPECT_FALSE(merge_tscs("merge", pathvec{ "absent.tsc" }, false));
	const std::string err = testing::internal::GetCapturedStderr();
	EXPECT_NE(err.find("TSC file contains no DATA section"), std::string::npos);
	EXPECT_NE(err.find("TSC file contains no scatterers"), std::string::npos);
	EXPECT_NE(err.find("Invalid TSC complex value: 1.0"), std::string::npos);
	EXPECT_NE(err.find("Scatterer and form-factor counts differ in TSC row"), std::string::npos);
	EXPECT_NE(err.find("Failed to open TSC file: absent.tsc"), std::string::npos);
	EXPECT_FALSE(fs::exists("combined.tscb"));
}

//form factors from a cube: a normalised Gaussian density exp(-|r - d|^2)/pi^(3/2) centred
//d = (0.5, 0, 0) bohr off the single H atom of a cubic P1 cell. The tsc holds the atom-relative
//transform, sum rho exp(+i k.(r - r_atom)), so F(hkl) = exp(-k^2/4) exp(i k.d) with
//k = 2 pi hkl / a: the offset gives every h != 0 reflection an imaginary part whose sign
//fixes the phase convention and whose axis fixes the k-vector ordering. The Becke electron
//count of the lone atom is one
TEST(ScatteringFactorCubeTests, GaussianCubeGivesAnalyticFormFactors)
{
	const fs::path cif = tmp("cube.cif");
	write_p1_cif(cif, "H1 H 0.0 0.0 0.0\n");
	WFN w(e_origin::NOT_YET_DEFINED);
	w.push_back_MO(0, 1.0, -1.0);
	w.push_back_atom("H1", 0.0, 0.0, 0.0, 1);
	double coef = 1.0;
	w.add_primitive(1, 1, 1.0, &coef);
	const int n = 101;
	const double h = 0.2;
	const double origin = -10.0;
	const double dx = 0.5;
	cube density({ n, n, n }, 1, true);
	for (int k = 0; k < 3; k++)
	{
		density.set_origin(k, origin);
		density.set_vector(k, k, h);
	}
	const double norm = std::pow(constants::PI, -1.5);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
			{
				const double x = origin + i * h - dx, y = origin + j * h, z = origin + k * h;
				density.set_value(i, j, k, norm * std::exp(-(x * x + y * y + z * z)));
			}
	options opt = make_options({ "-cif", cif.string(), "-group", "0", "-becke", "-no_gpu", "-no_date" });
	opt.m_hkl_list = { i3{ 0, 0, 0 }, i3{ 1, 0, 0 }, i3{ 0, 1, 0 }, i3{ 2, 0, 0 }, i3{ 1, 1, 0 } };
	std::ostringstream log;
	const itsc_block sf = calculate_scattering_factors_from_cube(opt, w, density, log);
	fs::remove(cif);
	EXPECT_NE(log.str().find("Table of Charges in electrons"), std::string::npos);
	const std::string tag = "Total number of partitioned electrons from cube: ";
	const std::size_t at = log.str().find(tag);
	ASSERT_NE(at, std::string::npos);
	EXPECT_NEAR(std::stod(log.str().substr(at + tag.size())), 1.0, 0.03);
	ASSERT_EQ(sf.scatterer_size(), 1u);
	ASSERT_EQ(sf.reflection_size(), 5u);
	EXPECT_EQ(std::get<atomID>(sf.get_scatterer(0)).to_hex_string(), atomID(0.0, 0.0, 0.0, 0, 1).to_hex_string());
	const double g = constants::TWO_PI / constants::ang2bohr(kA);
	const cvec& row = sf.get_sf_for_scatterer(0);
	for (std::size_t r = 0; r < 5; r++)
	{
		const std::array<int, 3> hkl = sf.get_indices(r);
		const double k2 = g * g * (hkl[0] * hkl[0] + hkl[1] * hkl[1] + hkl[2] * hkl[2]);
		const cdouble expected = std::exp(-k2 / 4.0) * std::polar(1.0, g * hkl[0] * dx);
		EXPECT_NEAR(row[r].real(), expected.real(), 0.03) << hkl[0] << hkl[1] << hkl[2];
		EXPECT_NEAR(row[r].imag(), expected.imag(), 0.03) << hkl[0] << hkl[1] << hkl[2];
	}
}

//get_interpolated_value is exact up to the last node (origin + (size - 1) * step) and 0 beyond it;
//the old expectation of 1.0 at x = 1.5 lay outside a 2-node grid and was extrapolation
TEST(ScatteringFactorCubeTests, CubeInterpolationLastCell)
{
	cube c({ 2, 2, 2 }, 0, true);
	for (int k = 0; k < 3; k++)
	{
		c.set_origin(k, 0.0);
		c.set_vector(k, k, 1.0);
	}
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
				c.set_value(i, j, k, 1.0);
	EXPECT_NEAR(c.get_interpolated_value(0.5, 0.5, 0.5), 1.0, 1e-12);
	EXPECT_NEAR(c.get_interpolated_value(1.0, 0.5, 0.5), 1.0, 1e-12);
	EXPECT_EQ(c.get_interpolated_value(1.5, 0.5, 0.5), 0.0);
}
