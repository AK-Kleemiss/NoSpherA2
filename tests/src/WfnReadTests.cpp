#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/atoms.h"
#include "core/cell.h"
#include "core/wfn_class.h"
#include "core/mo_class.h"
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#undef I

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

namespace
{
	static std::filesystem::path tmp_path(const std::string& name)
	{
		return std::filesystem::temp_directory_path() / ("nosphera2_wfnread_" + name);
	}

	static void write_text(const std::filesystem::path& p, const std::string& text)
	{
		std::ofstream f(p);
		f << text;
	}

	//two H atoms, three MOs (the third restarts the energies, which the wfn/wfx readers take as a second spin),
	//three primitives whose coefficients per MO are known so every reader can be checked value by value
	static WFN make_h2_wfn()
	{
		WFN w(e_origin::wfn);
		w.push_back_atom("H", 0.0, 0.0, -0.7, 1);
		w.push_back_atom("H", 0.0, 0.0, 0.7, 1);
		w.push_back_MO(1, 2.0, -0.6);
		w.push_back_MO(2, 0.0, 0.3);
		w.push_back_MO(3, 1.0, -0.5, 1);
		double c0[3] = { 0.5, 0.7, 0.4 };
		double c1[3] = { 0.5, -0.7, 0.4 };
		double c2[3] = { 0.1, 0.2, -0.3 };
		w.add_primitive(1, 1, 1.2, c0);
		w.add_primitive(2, 1, 1.2, c1);
		w.add_primitive(1, 2, 0.9, c2);
		return w;
	}

	//analytic overlap of cartesian primitives: per axis the binomial expansion of (x - A)^la (x - B)^lb about
	//P = (aA + bB)/p, sum_ij C(la,i) C(lb,j) PA^(la-i) PB^(lb-j) (i+j-1)!! / (2p)^((i+j)/2) sqrt(pi/p) for even
	//i + j, times exp(-ab/p |AB|^2); any l, so the def2-SVP-like epoxide gbw with d functions is covered
	static double axis_overlap(const int la, const int lb, const double PA, const double PB, const double p)
	{
		auto binom = [](const int n, const int k) { return std::tgamma(n + 1.0) / (std::tgamma(k + 1.0) * std::tgamma(n - k + 1.0)); };
		double sum = 0.0;
		for (int i = 0; i <= la; i++)
			for (int j = 0; j <= lb; j++)
			{
				if ((i + j) % 2) continue;
				double dfact = 1.0;
				for (int k = i + j - 1; k > 1; k -= 2) dfact *= k;
				sum += binom(la, i) * binom(lb, j) * std::pow(PA, la - i) * std::pow(PB, lb - j) * dfact / std::pow(2.0 * p, 0.5 * (i + j));
			}
		return sum * std::sqrt(constants::PI / p);
	}

	static vec2 sp_primitive_overlap(const WFN& w)
	{
		const int nex = w.get_nex();
		vec2 S(nex, vec(nex));
		for (int a = 0; a < nex; a++)
			for (int b = 0; b <= a; b++)
			{
				int la[3], lb[3];
				constants::type2vector(w.get_type(a), la);
				constants::type2vector(w.get_type(b), lb);
				const double al = w.get_exponent(a), be = w.get_exponent(b), p = al + be;
				double s = 1.0, AB2 = 0.0;
				for (int k = 0; k < 3; k++)
				{
					const double A = w.get_atom_coordinate(w.get_center(a) - 1, k), B = w.get_atom_coordinate(w.get_center(b) - 1, k), P = (al * A + be * B) / p;
					AB2 += (A - B) * (A - B);
					s *= axis_overlap(la[k], lb[k], P - A, P - B, p);
				}
				S[a][b] = S[b][a] = s * std::exp(-al * be / p * AB2);
			}
		return S;
	}

	//C S C^T = 1 over every MO: a wrong contraction norm, a dropped primitive or a shifted coefficient column all break it
	static void expect_orthonormal(const WFN& w, const double tol)
	{
		const int nex = w.get_nex(), nmo = w.get_nmo();
		const vec2 S = sp_primitive_overlap(w);
		vec2 SC(nmo, vec(nex, 0.0));
		for (int i = 0; i < nmo; i++)
			for (int a = 0; a < nex; a++)
				for (int b = 0; b < nex; b++)
					SC[i][a] += S[a][b] * w.get_MO_coef(i, b);
		for (int i = 0; i < nmo; i++)
			for (int j = 0; j <= i; j++)
			{
				double o = 0.0;
				for (int a = 0; a < nex; a++)
					o += w.get_MO_coef(i, a) * SC[j][a];
				EXPECT_NEAR(o, i == j ? 1.0 : 0.0, tol) << "MOs " << i << " " << j;
			}
	}

	static void expect_same_density(const WFN& a, const WFN& b, const double rtol)
	{
		ASSERT_EQ(a.get_ncen(), b.get_ncen());
		for (int i = 0; i < a.get_ncen(); i++)
			for (const double r : { 0.3, 1.1, 2.5 })
			{
				const d3 p{ a.get_atom_coordinate(i, 0) + r, a.get_atom_coordinate(i, 1) + 0.5 * r, a.get_atom_coordinate(i, 2) - 0.7 * r };
				const double da = a.compute_dens(p);
				EXPECT_NEAR(da, b.compute_dens(p), rtol * std::max(1.0, std::abs(da))) << "atom " << i << " r " << r;
			}
	}

	//------------------------------------------------------------ mutators ------------------------------------------------------------

	//push_back_atom with a charge below one still counts the slot, erase_atom shifts the rest down and
	//get_centers prints either unit; the bohr/angstrom strings must differ by exactly the conversion factor
	TEST(WfnReadTests, AtomPushEraseAndCentreListing)
	{
		WFN w(e_origin::wfn);
		EXPECT_TRUE(w.push_back_atom("O", 1.0, 2.0, 3.0, 8));
		EXPECT_FALSE(w.push_back_atom("X", 0.0, 0.0, 0.0, 0));
		EXPECT_EQ(w.get_ncen(), 2);
		EXPECT_TRUE(w.push_back_atom(atom("H", atomID(), 3, 0.5, 0.0, 0.0, 1)));
		EXPECT_EQ(w.get_ncen(), 3);
		EXPECT_TRUE(w.erase_atom(1));
		EXPECT_EQ(w.get_ncen(), 2);
		EXPECT_EQ(w.get_atom_label(1), "H");
		EXPECT_NEAR(w.get_atom_pos(1)[0], 0.5, 1e-12);
		EXPECT_NEAR(w.get_atom_coordinate(0, 2), 3.0, 1e-12);
		EXPECT_EQ(w.get_atom_charge(0), 8);
		EXPECT_EQ(w.get_nr_electrons(), 9u);
		const std::string bohr = w.get_centers(true), ang = w.get_centers(false);
		EXPECT_EQ(bohr, "O 1.000000 2.000000 3.000000\nH 0.500000 0.000000 0.000000\n");
		std::istringstream is(ang);
		std::string label;
		double x, y, z;
		is >> label >> x >> y >> z;
		EXPECT_EQ(label, "O");
		EXPECT_NEAR(x, constants::bohr2ang(1.0), 1e-6);
		EXPECT_NEAR(z, constants::bohr2ang(3.0), 1e-6);
		//stdout listings, no return value: exercised for crash freedom only
		w.list_centers();
		w.print_atom_long(0);
	}

	//push_back_MO by number, by MO object and with an operator; delete_MO and clear_MOs keep nmo in step and
	//get_MO_op_count splits the spins; hdr(true) counts only occupied MOs
	TEST(WfnReadTests, MoPushDeleteClearAndHeader)
	{
		WFN w(e_origin::wfn);
		w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
		EXPECT_TRUE(w.push_back_MO(1, 2.0, -0.9));
		EXPECT_TRUE(w.push_back_MO(MO(2, 0.0, 0.4)));
		EXPECT_TRUE(w.push_back_MO(3, 1.0, -0.8, 1));
		EXPECT_EQ(w.get_nmo(), 3);
		EXPECT_EQ(w.get_nmo(true), 2);
		EXPECT_EQ(w.get_MO_op_count(0), 2);
		EXPECT_EQ(w.get_MO_op_count(1), 1);
		EXPECT_EQ(w.get_MO_op(2), 1);
		EXPECT_NEAR(w.get_MO(1).get_energy(), 0.4, 1e-12);
		EXPECT_NEAR(w.get_MO_occ(2), 1.0, 1e-12);
		EXPECT_EQ(w.hdr(false), "GTO                   3 MOL ORBITALS      0 PRIMITIVES        1 NUCLEI\n");
		EXPECT_EQ(w.hdr(true), "GTO                   2 MOL ORBITALS      0 PRIMITIVES        1 NUCLEI\n");
		w.delete_MO(1);
		EXPECT_EQ(w.get_nmo(), 2);
		EXPECT_NEAR(w.get_MO_energy(1), -0.8, 1e-12);
		w.clear_MOs();
		EXPECT_EQ(w.get_nmo(), 0);
		EXPECT_EQ(w.get_MO_op_count(0), 0);
	}

	//add_primitive appends one coefficient per MO, add_exp on a missing centre fails, remove_primitive drops the
	//slot from centres, types, exponents and every MO so the remaining coefficients shift down by one
	TEST(WfnReadTests, PrimitiveAddRemoveShiftsCoefficients)
	{
		WFN w = make_h2_wfn();
		EXPECT_EQ(w.get_nex(), 3);
		EXPECT_EQ(w.get_MO_primitive_count(2), 3);
		EXPECT_EQ(w.get_center(2), 1);
		EXPECT_EQ(w.get_type(2), 2);
		EXPECT_NEAR(w.get_exponent(2), 0.9, 1e-12);
		EXPECT_NEAR(w.get_MO_coef(1, 1), -0.7, 1e-12);
		EXPECT_NEAR(w.get_MO_coef_f(2, 2), -0.3, 1e-12);
		EXPECT_NEAR(w.get_MO_coef_ptr(0)[2], 0.1, 1e-12);
		EXPECT_TRUE(w.set_MO_coef(0, 2, 0.25));
		EXPECT_NEAR(w.get_MO_coef(0, 2), 0.25, 1e-12);
		EXPECT_NEAR(w.get_maximum_MO_coefficient(true), 0.5, 1e-12);
		EXPECT_NEAR(w.get_maximum_MO_coefficient(false), 0.7, 1e-12);
		EXPECT_TRUE(w.remove_primitive(1));
		EXPECT_EQ(w.get_nex(), 2);
		EXPECT_EQ(w.get_types().size(), 2u);
		EXPECT_EQ(w.get_exponents().size(), 2u);
		EXPECT_EQ(w.get_centers().size(), 2u);
		EXPECT_EQ(w.get_MO_primitive_count(0), 2);
		EXPECT_EQ(w.get_center(0), 2);
		EXPECT_NEAR(w.get_MO_coef(1, 0), -0.7, 1e-12);
		w.list_primitives();
	}

	//add_exp on a centre that does not exist fails without leaving a trace: nex stays in step with the
	//primitive arrays so the next remove_primitive still shifts the MO coefficients
	TEST(WfnReadTests, AddExpOnMissingCentreLeavesNexUnchanged)
	{
		WFN w = make_h2_wfn();
		EXPECT_FALSE(w.add_exp(5, 1, 2.0));
		EXPECT_EQ(w.get_nex(), 3);
		EXPECT_TRUE(w.remove_primitive(1));
		EXPECT_EQ(w.get_nex(), 2);
		EXPECT_NEAR(w.get_MO_coef(1, 0), -0.7, 1e-12);
	}

	//a cartesian p shell pushed with an order and scale: coefficient (cart, prim) = shell[order[cart]][prim] * scale[cart],
	//the exponents, centres and types are appended only for MO 0 and the density along x is the analytic
	//occ * (cs * exp(-as x^2) + cx * x * exp(-ap x^2))^2 of the pushed primitives
	TEST(WfnReadTests, CartesianShellPushOrderScaleAndDensity)
	{
		WFN w(e_origin::wfn);
		w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
		w.push_back_MO(1, 1.0, -0.5);
		w.push_back_MO(2, 0.0, 0.2);
		const std::vector<primitive> prims{ primitive(1, 1, 1.3, 1.0), primitive(1, 2, 0.8, 1.0) };
		const vec2 s_shell{ { 0.6 } };
		const vec2 p_shell{ { 0.1 }, { 0.2 }, { 0.3 } };
		const int order[3] = { 2, 0, 1 };
		const double scale[3] = { 1.0, 2.0, 3.0 };
		for (int mo = 0; mo < 2; mo++)
		{
			w.push_back_cartesian_shell(mo, 0, s_shell, prims, 0, 1);
			w.push_back_cartesian_shell(mo, 1, p_shell, prims, 1, 1, order, scale);
		}
		EXPECT_EQ(w.get_nex(), 4);
		EXPECT_EQ(w.get_MO_primitive_count(1), 4);
		EXPECT_EQ(w.get_type(0), 1);
		EXPECT_EQ(w.get_type(1), 2);
		EXPECT_EQ(w.get_type(3), 4);
		EXPECT_NEAR(w.get_exponent(0), 1.3, 1e-12);
		EXPECT_NEAR(w.get_exponent(2), 0.8, 1e-12);
		EXPECT_EQ(w.get_center(3), 1);
		EXPECT_NEAR(w.get_MO_coef(0, 1), 0.3 * 1.0, 1e-12);
		EXPECT_NEAR(w.get_MO_coef(0, 2), 0.1 * 2.0, 1e-12);
		EXPECT_NEAR(w.get_MO_coef(0, 3), 0.2 * 3.0, 1e-12);
		EXPECT_NEAR(w.get_MO_coef(1, 3), 0.6, 1e-12);
		const double x = 0.7;
		const double expected = std::pow(0.6 * std::exp(-1.3 * x * x) + 0.3 * x * std::exp(-0.8 * x * x), 2);
		EXPECT_NEAR(w.compute_dens(d3{ x, 0.0, 0.0 }), expected, 1e-10);
		EXPECT_NEAR(w.compute_dens(d3{ 0.0, 0.0, 0.0 }), 0.36, 1e-10);
	}

	//push_back_spherical_shell maps the m components onto cartesians through sph2cart: for l = 0 the identity, for
	//l = 1 the molden m order 0, +1, -1 lands on z, x, y (constants.cpp p_mat), so 0.3/0.5/0.7 become pz/px/py with
	//the type run 2, 3, 4; the density along +y then carries the py coefficient alone
	TEST(WfnReadTests, SphericalShellPushIsPermutationForP)
	{
		WFN w(e_origin::wfn);
		w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
		w.push_back_MO(1, 2.0, -0.5);
		const std::vector<primitive> prims{ primitive(1, 1, 1.1, 1.0), primitive(1, 2, 0.7, 1.0) };
		w.push_back_spherical_shell(0, 0, vec2{ { 0.9 } }, prims, 0, 1);
		w.push_back_spherical_shell(0, 1, vec2{ { 0.3 }, { 0.5 }, { 0.7 } }, prims, 1, 1);
		ASSERT_EQ(w.get_nex(), 4);
		EXPECT_NEAR(w.get_MO_coef(0, 0), 0.9, 1e-12);
		EXPECT_NEAR(w.get_MO_coef(0, 1), 0.5, 1e-12);
		EXPECT_NEAR(w.get_MO_coef(0, 2), 0.7, 1e-12);
		EXPECT_NEAR(w.get_MO_coef(0, 3), 0.3, 1e-12);
		for (int i = 1; i < 4; i++)
		{
			EXPECT_EQ(w.get_type(i), i + 1);
			EXPECT_NEAR(w.get_exponent(i), 0.7, 1e-12);
		}
		w.set_exp_cutoff();
		const double y = 0.6;
		const double expected = 2.0 * std::pow(0.9 * std::exp(-1.1 * y * y) + 0.7 * y * std::exp(-0.7 * y * y), 2);
		EXPECT_NEAR(w.compute_dens(d3{ 0.0, y, 0.0 }), expected, 1e-10);
	}

	//extract_xyz converts only when the stored and requested units differ: a hand-built wavefunction carries no
	//bohr flag so "bohr" grows by ang2bohr and "angstrom" is a copy; a read xyz is flagged bohr and shrinks by bohr2ang
	TEST(WfnReadTests, ExtractXyzConvertsUnitsBothWays)
	{
		WFN w = make_h2_wfn();
		ASSERT_FALSE(w.get_isBohr());
		const std::vector<asym_atom> grown = w.extract_xyz("bohr");
		ASSERT_EQ(grown.size(), 2u);
		EXPECT_EQ(grown[1].label, "H");
		EXPECT_NEAR(grown[1].pos[2], constants::ang2bohr(0.7), 1e-12);
		const std::vector<asym_atom> same = w.extract_xyz("angstrom");
		EXPECT_NEAR(same[1].pos[2], 0.7, 1e-12);
		const std::filesystem::path xyz = tmp_path("units.xyz");
		write_text(xyz, "1\ncomment\nO 0.0 0.0 1.0\n");
		WFN a(e_origin::xyz);
		std::ostringstream log;
		ASSERT_TRUE(a.read_xyz(xyz, log));
		std::filesystem::remove(xyz);
		EXPECT_NEAR(a.get_atom_coordinate(0, 2), constants::ang2bohr(1.0), 1e-12);
		EXPECT_NEAR(a.extract_xyz("angstrom")[0].pos[2], 1.0, 1e-12);
		EXPECT_NEAR(a.extract_xyz("bohr")[0].pos[2], constants::ang2bohr(1.0), 1e-12);
	}

	//the per-atom basis set container: entries land on the right atom, the shell counter grows with new shell
	//indices, clear_atom_basis_set empties one atom and the loaded count follows; plus the small flag setters
	TEST(WfnReadTests, AtomBasisSetBookkeepingAndFlags)
	{
		WFN w = make_h2_wfn();
		EXPECT_EQ(w.get_nr_basis_set_loaded(), 0);
		EXPECT_TRUE(w.push_back_atom_basis_set(0, 3.0, 0.5, 1, 0));
		EXPECT_TRUE(w.push_back_atom_basis_set(0, 0.5, 0.6, 1, 0));
		EXPECT_TRUE(w.push_back_atom_basis_set(0, 0.8, 1.0, 2, 1));
		EXPECT_FALSE(w.push_back_atom_basis_set(7, 0.8, 1.0, 2, 0));
		EXPECT_EQ(w.get_atom_basis_set_size(0), 3);
		EXPECT_EQ(w.get_atom_basis_set_size(1), 0);
		EXPECT_EQ(w.get_nr_basis_set_loaded(), 1);
		EXPECT_EQ(w.get_atom_shell_count(0), 2);
		EXPECT_EQ(w.get_shell_type(0, 1), 2);
		EXPECT_EQ(w.get_shell_start(0, 1), 2);
		EXPECT_EQ(w.get_shell_end(0, 0), 1);
		w.clear_atom_basis_set(0);
		EXPECT_EQ(w.get_atom_basis_set_size(0), 0);
		EXPECT_EQ(w.get_nr_basis_set_loaded(), 0);
		w.change_basis_set_name("def2-svp");
		EXPECT_EQ(w.get_basis_set_name(), "def2-svp");
		w.set_d_f_switch(true);
		EXPECT_TRUE(w.get_d_f_switch());
		w.set_charge(-1);
		w.set_multi(3);
		EXPECT_EQ(w.get_charge(), -1);
		EXPECT_EQ(w.get_multi(), 3u);
		EXPECT_EQ(w.get_nr_electrons(), 3u);
		EXPECT_EQ(w.get_SDM_size(), 0);
		EXPECT_EQ(w.get_origin(), e_origin::wfn);
	}

	//get_norm_const on single-primitive s, p, d and f shells with unit contraction coefficient: the shell
	//factor is one, so every entry is the primitive norm (2a/pi)^3/4, (128 a^5/pi^3)^1/4, ... with the
	//sqrt(3), sqrt(5), sqrt(15) multipliers of the off-axis cartesians, 1 + 3 + 6 + 10 entries in total
	TEST(WfnReadTests, NormConstSinglePrimitiveShellsGivePrimitiveNorms)
	{
		WFN w(e_origin::wfn);
		w.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		const double a_s = 1.3, a_p = 0.9, a_d = 0.7, a_f = 0.5;
		ASSERT_TRUE(w.push_back_atom_basis_set(0, a_s, 1.0, 1, 0));
		ASSERT_TRUE(w.push_back_atom_basis_set(0, a_p, 1.0, 2, 1));
		ASSERT_TRUE(w.push_back_atom_basis_set(0, a_d, 1.0, 3, 2));
		ASSERT_TRUE(w.push_back_atom_basis_set(0, a_f, 1.0, 4, 3));
		std::ostringstream log;
		const vec n = w.get_norm_const(log);
		ASSERT_EQ(n.size(), 20u);
		const double pi = constants::PI, pi3 = pi * pi * pi;
		const double ns = std::pow(2.0 * a_s / pi, 0.75);
		const double np = std::pow(128.0 * std::pow(a_p, 5) / pi3, 0.25);
		const double nd = std::pow(2048.0 * std::pow(a_d, 7) / (9.0 * pi3), 0.25);
		const double nf = std::pow(32768.0 * std::pow(a_f, 9) / (225.0 * pi3), 0.25);
		EXPECT_NEAR(n[0], ns, 1e-12);
		for (int i = 1; i < 4; i++)
			EXPECT_NEAR(n[i], np, 1e-12);
		for (int i = 4; i < 7; i++)
			EXPECT_NEAR(n[i], nd, 1e-12);
		for (int i = 7; i < 10; i++)
			EXPECT_NEAR(n[i], std::sqrt(3.0) * nd, 1e-12);
		for (int i = 10; i < 13; i++)
			EXPECT_NEAR(n[i], nf, 1e-12);
		for (int i = 13; i < 19; i++)
			EXPECT_NEAR(n[i], std::sqrt(5.0) * nf, 1e-12);
		EXPECT_NEAR(n[19], std::sqrt(15.0) * nf, 1e-12);
	}

	//a two-primitive contracted s shell: the returned coefficients must keep the ratio c_i (2a_i/pi)^3/4 of the
	//input and the contracted function must integrate to one under the analytic gaussian overlap (pi/(a_i+a_j))^3/2
	TEST(WfnReadTests, NormConstContractedSShellIntegratesToOne)
	{
		WFN w(e_origin::wfn);
		w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
		const double a[2] = { 2.0, 0.5 }, c[2] = { 0.6, 0.4 };
		ASSERT_TRUE(w.push_back_atom_basis_set(0, a[0], c[0], 1, 0));
		ASSERT_TRUE(w.push_back_atom_basis_set(0, a[1], c[1], 1, 0));
		std::ostringstream log;
		const vec n = w.get_norm_const(log, true);
		ASSERT_EQ(n.size(), 2u);
		const double r0 = c[0] * std::pow(2.0 * a[0] / constants::PI, 0.75), r1 = c[1] * std::pow(2.0 * a[1] / constants::PI, 0.75);
		EXPECT_NEAR(n[0] / n[1], r0 / r1, 1e-12);
		double overlap = 0.0;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				overlap += n[i] * n[j] * std::pow(constants::PI / (a[i] + a[j]), 1.5);
		EXPECT_NEAR(overlap, 1.0, 1e-12);
	}

	//remove_center drops the atom (1-based) together with its primitives and leaves the others intact
	TEST(WfnReadTests, RemoveCenterDropsItsPrimitives)
	{
		WFN w = make_h2_wfn();
		EXPECT_FALSE(w.remove_center(3));
		EXPECT_TRUE(w.remove_center(2));
		EXPECT_EQ(w.get_ncen(), 1);
		EXPECT_EQ(w.get_nex(), 2);
		EXPECT_EQ(w.get_MO_primitive_count(0), 2);
		EXPECT_EQ(w.get_center(0), 1);
		EXPECT_EQ(w.get_center(1), 1);
		EXPECT_NEAR(w.get_MO_coef(2, 1), -0.3, 1e-12);
	}

	//------------------------------------------------------------ wfn / wfx / xyz ------------------------------------------------------------

	//write_wfn then read_wfn on the in-code wavefunction: atoms, primitives, coefficients, occupations and energies
	//come back, the energy restart of MO 3 marks the file unrestricted and the operators split 2/1
	TEST(WfnReadIoTests, WfnRoundTripKeepsEverything)
	{
		WFN w = make_h2_wfn();
		const std::filesystem::path p = tmp_path("roundtrip.wfn");
		ASSERT_TRUE(w.write_wfn(p, false, false));
		WFN r(e_origin::wfn);
		std::ostringstream log;
		ASSERT_TRUE(r.read_wfn(p, true, log));
		std::filesystem::remove(p);
		EXPECT_NE(log.str().find("e_nmo: 3, e_nex: 3, e_nuc : 2"), std::string::npos);
		EXPECT_EQ(r.get_ncen(), 2);
		EXPECT_EQ(r.get_nex(), 3);
		EXPECT_EQ(r.get_nmo(), 3);
		EXPECT_TRUE(r.get_isBohr());
		EXPECT_TRUE(r.get_is_unrestricted());
		EXPECT_EQ(r.get_MO_op_count(0), 2);
		EXPECT_EQ(r.get_MO_op_count(1), 1);
		EXPECT_EQ(r.get_atom_label(0), "H");
		EXPECT_NEAR(r.get_atom_coordinate(1, 2), 0.7, 1e-8);
		for (int p_i = 0; p_i < 3; p_i++)
		{
			EXPECT_EQ(r.get_center(p_i), w.get_center(p_i));
			EXPECT_EQ(r.get_type(p_i), w.get_type(p_i));
			EXPECT_NEAR(r.get_exponent(p_i), w.get_exponent(p_i), 1e-7);
			for (int m = 0; m < 3; m++)
				EXPECT_NEAR(r.get_MO_coef(m, p_i), w.get_MO_coef(m, p_i), 1e-8);
		}
		for (int m = 0; m < 3; m++)
		{
			EXPECT_NEAR(r.get_MO_occ(m), w.get_MO_occ(m), 1e-8);
			EXPECT_NEAR(r.get_MO_energy(m), w.get_MO_energy(m), 1e-6);
		}
		expect_same_density(w, r, 1e-7);
		//only occupied MOs: the header and the file shrink to two MOs
		ASSERT_TRUE(w.write_wfn(p, false, true));
		WFN o(e_origin::wfn);
		ASSERT_TRUE(o.read_wfn(p, false, log));
		std::filesystem::remove(p);
		EXPECT_EQ(o.get_nmo(), 2);
		EXPECT_NEAR(o.get_MO_occ(1), 1.0, 1e-8);
	}

	//a hand-written Fortran-column wfn: fields that touch, a 'D' exponent and an '=' glued to the orbital energy
	//are all cut by width, not by whitespace
	TEST(WfnReadIoTests, WfnHandWrittenFortranColumns)
	{
		const std::string text =
			"hand written\n"
			"GAUSSIAN 1 MOL ORBITALS 2 PRIMITIVES 1 NUCLEI\n"
			"H      1    (CENTRE  1)   0.00000000  0.00000000 -0.70000000  CHARGE =  1.0\n"
			"CENTRE ASSIGNMENTS    1  1\n"
			"TYPE ASSIGNMENTS      1  2\n"
			"EXPONENTS  0.1200000D+01 0.9000000D+00\n"
			"MO    1    MO 0.0        OCC NO =    2.0000000  ORB. ENERGY  =-0.50000000\n"
			"  0.50000000E+00  0.30000000E+00\n"
			"END DATA\n";
		const std::filesystem::path p = tmp_path("hand.wfn");
		write_text(p, text);
		WFN r(e_origin::wfn);
		std::ostringstream log;
		ASSERT_TRUE(r.read_wfn(p, false, log));
		std::filesystem::remove(p);
		EXPECT_EQ(r.get_comment(), "hand written");
		EXPECT_EQ(r.get_ncen(), 1);
		EXPECT_EQ(r.get_atom_label(0), "H");
		EXPECT_NEAR(r.get_atom_coordinate(0, 2), -0.7, 1e-8);
		EXPECT_EQ(r.get_atom_charge(0), 1);
		ASSERT_EQ(r.get_nex(), 2);
		EXPECT_NEAR(r.get_exponent(0), 1.2, 1e-12);
		EXPECT_NEAR(r.get_exponent(1), 0.9, 1e-12);
		EXPECT_EQ(r.get_type(1), 2);
		ASSERT_EQ(r.get_nmo(), 1);
		EXPECT_NEAR(r.get_MO_occ(0), 2.0, 1e-12);
		EXPECT_NEAR(r.get_MO_energy(0), -0.5, 1e-7);
		EXPECT_NEAR(r.get_MO_coef(0, 0), 0.5, 1e-12);
		EXPECT_NEAR(r.get_MO_coef(0, 1), 0.3, 1e-12);
		EXPECT_FALSE(r.get_is_unrestricted());
	}

	//read_wfn into a wavefunction that already holds atoms refuses and says so instead of appending
	TEST(WfnReadIoTests, ReadWfnRefusesSecondLoad)
	{
		WFN w = make_h2_wfn();
		const std::filesystem::path p = tmp_path("second.wfn");
		ASSERT_TRUE(w.write_wfn(p, false, false));
		std::ostringstream log;
		EXPECT_FALSE(w.read_wfn(p, false, log));
		std::filesystem::remove(p);
		EXPECT_FALSE(log.str().empty());
		EXPECT_EQ(w.get_ncen(), 2);
	}

	//write_wfx then read_wfx: charge and multiplicity travel through the tags, the energy restart marks the
	//second spin and a comment line dropped into the coefficient block is skipped rather than parsed
	TEST(WfnReadIoTests, WfxRoundTripSkipsForeignLinesInCoefficientBlock)
	{
		WFN w = make_h2_wfn();
		w.set_charge(-1);
		w.set_multi(2);
		const std::filesystem::path p = tmp_path("roundtrip.wfx");
		ASSERT_TRUE(w.write_wfx(p, false));
		std::ifstream in(p);
		std::stringstream buf;
		buf << in.rdbuf();
		in.close();
		std::string text = buf.str();
		const size_t at = text.find("<Molecular Orbital Primitive Coefficients>");
		ASSERT_NE(at, std::string::npos);
		text.insert(text.find('\n', at) + 1, "comment line without a tag\n");
		write_text(p, text);
		WFN r(e_origin::wfx);
		std::ostringstream log;
		ASSERT_TRUE(r.read_wfx(p, false, log));
		std::filesystem::remove(p);
		EXPECT_EQ(r.get_charge(), -1);
		EXPECT_EQ(r.get_multi(), 2u);
		EXPECT_EQ(r.get_ncen(), 2);
		EXPECT_EQ(r.get_atom_label(1), "H2");
		EXPECT_EQ(r.get_nex(), 3);
		EXPECT_EQ(r.get_nmo(), 3);
		EXPECT_TRUE(r.get_is_unrestricted());
		EXPECT_EQ(r.get_MO_op_count(1), 1);
		for (int m = 0; m < 3; m++)
		{
			EXPECT_NEAR(r.get_MO_occ(m), w.get_MO_occ(m), 1e-8);
			EXPECT_NEAR(r.get_MO_energy(m), w.get_MO_energy(m), 1e-8);
			for (int p_i = 0; p_i < 3; p_i++)
				EXPECT_NEAR(r.get_MO_coef(m, p_i), w.get_MO_coef(m, p_i), 1e-8);
		}
		expect_same_density(w, r, 1e-8);
	}

	//read_xyz converts angstrom to bohr, tolerates trailing blanks and a further frame's count line, and the
	//debug listing reports the atoms; write_xyz writes it back in angstrom with the NoSpherA2 comment
	TEST(WfnReadIoTests, XyzReadToleratesTrailingLinesAndWritesBack)
	{
		const std::filesystem::path p = tmp_path("frames.xyz");
		write_text(p, "2\nwater fragment\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\n\n2\n");
		WFN w(e_origin::xyz);
		std::ostringstream log;
		ASSERT_TRUE(w.read_xyz(p, log, true));
		std::filesystem::remove(p);
		EXPECT_EQ(w.get_ncen(), 2);
		EXPECT_EQ(w.get_atom_charge(0), 8);
		EXPECT_EQ(w.get_atom_label(1), "H");
		EXPECT_NEAR(w.get_atom_coordinate(1, 2), constants::ang2bohr(0.96), 1e-10);
		EXPECT_TRUE(w.get_isBohr());
		EXPECT_EQ(w.get_origin(), e_origin::xyz);
		EXPECT_NE(log.str().find("e_nuc"), std::string::npos);
		const std::filesystem::path q = tmp_path("written.xyz");
		ASSERT_TRUE(w.write_xyz(q));
		std::ifstream in(q);
		std::string line;
		std::getline(in, line);
		EXPECT_EQ(line, "2");
		std::getline(in, line);
		EXPECT_NE(line.find("NoSpherA2"), std::string::npos);
		std::getline(in, line);
		std::getline(in, line);
		in.close();
		std::istringstream is(line);
		std::string label;
		double x, y, z;
		is >> label >> x >> y >> z;
		EXPECT_EQ(label, "H");
		EXPECT_NEAR(z, 0.96, 1e-7);
		WFN back(q);
		std::filesystem::remove(q);
		EXPECT_NEAR(back.get_atom_coordinate(1, 2), w.get_atom_coordinate(1, 2), 1e-7);
	}

	//the charge/multiplicity constructor dispatches on the extension and keeps the given charge, which the
	//electron count then reflects
	TEST(WfnReadIoTests, ConstructorWithChargeDispatchesXyz)
	{
		const std::filesystem::path p = tmp_path("ctor.xyz");
		write_text(p, "1\nlithium cation\nLi 0.0 0.0 0.0\n");
		WFN w(p, 1, 1);
		std::filesystem::remove(p);
		EXPECT_EQ(w.get_charge(), 1);
		EXPECT_EQ(w.get_multi(), 1u);
		EXPECT_EQ(w.get_ncen(), 1);
		EXPECT_EQ(w.get_atom_charge(0), 3);
		EXPECT_EQ(w.get_nr_electrons(), 2u);
		EXPECT_EQ(w.get_origin(), e_origin::xyz);
	}

	//------------------------------------------------------------ molden ------------------------------------------------------------

	//a cartesian molden with two s shells on one H and a single MO: the coefficient matrix is padded square,
	//the density is the analytic occ * (c1 exp(-a1 r^2) + c2 exp(-a2 r^2))^2 and the [5D][7F][9G] copy reads
	//identically because s shells have no spherical/cartesian difference
	TEST(WfnReadIoTests, MoldenSyntheticSShellsPadAndMatchAnalyticDensity)
	{
		const std::string head =
			"[Molden Format]\n"
			"[Title]\n"
			" synthetic H\n"
			"[Atoms] AU\n"
			"H     1     1   0.000000   0.000000   0.000000\n"
			"[GTO]\n"
			"  1 0\n"
			" s    1 1.00\n"
			"     1.2000000000      1.0000000000\n"
			" s    1 1.00\n"
			"     0.4000000000      1.0000000000\n"
			"\n";
		const std::string mos =
			"[MO]\n"
			" Sym= A\n"
			" Ene= -0.50\n"
			" Spin= Alpha\n"
			" Occup= 1.0\n"
			"   1   0.6\n"
			"   2   0.3\n";
		const std::filesystem::path p = tmp_path("synth.molden"), q = tmp_path("synth_sph.molden");
		write_text(p, head + mos);
		write_text(q, head + "[5D]\n[7F]\n[9G]\n" + mos);
		WFN cart(e_origin::molden), sph(e_origin::molden);
		std::ostringstream log;
		ASSERT_TRUE(cart.read_molden(p, log, true));
		ASSERT_TRUE(sph.read_molden(q, log, false));
		std::filesystem::remove(p);
		std::filesystem::remove(q);
		EXPECT_NE(log.str().find("File is valid"), std::string::npos);
		//the reader appends the title to the constructor's default comment "Test" (wfn_class.cpp:51, 1292)
		EXPECT_NE(cart.get_comment().find("synthetic H"), std::string::npos);
		EXPECT_EQ(cart.get_ncen(), 1);
		EXPECT_EQ(cart.get_nex(), 2);
		EXPECT_EQ(cart.get_nmo(), 1);
		EXPECT_TRUE(cart.get_d_f_switch());
		EXPECT_FALSE(sph.get_d_f_switch());
		EXPECT_FALSE(cart.get_is_unrestricted());
		EXPECT_NEAR(cart.get_MO_coef(0, 0), 0.6, 1e-12);
		EXPECT_NEAR(cart.get_MO_coef(0, 1), 0.3, 1e-12);
		for (const double r : { 0.0, 0.5, 1.3 })
		{
			const double expected = std::pow(0.6 * std::exp(-1.2 * r * r) + 0.3 * std::exp(-0.4 * r * r), 2);
			EXPECT_NEAR(cart.compute_dens(d3{ r, 0.0, 0.0 }), expected, 1e-10);
			EXPECT_NEAR(sph.compute_dens(d3{ 0.0, r, 0.0 }), expected, 1e-10);
		}
	}

	//the F2 molden without spherical flags takes the cartesian shell path: 2 x (4 s + 4 p x 3) = 32 primitives, 8
	//restricted MOs of a valence-only (7 electrons per F) basis and a density that is mirror symmetric about the bond
	//midpoint; the file's contraction coefficients carry the primitive norms (both contractions integrate to one by
	//hand), so the eight MOs must be orthonormal under the analytic primitive overlap, which a consistent sign flip
	//or rescaling of every p column would pass the mirror check with
	TEST(WfnReadIoTests, MoldenCartesianF2IsMirrorSymmetric)
	{
		const std::filesystem::path p = nos_test_repo_root() / "tests" / "molden_file" / "F2.molden";
		if (!std::filesystem::exists(p))
			GTEST_SKIP() << "fixture missing: " << p;
		WFN w(e_origin::molden);
		std::ostringstream log;
		ASSERT_TRUE(w.read_molden(p, log, false));
		EXPECT_EQ(w.get_ncen(), 2);
		EXPECT_EQ(w.get_nex(), 32);
		EXPECT_EQ(w.get_nmo(), 8);
		EXPECT_FALSE(w.get_is_unrestricted());
		EXPECT_TRUE(w.get_d_f_switch());
		EXPECT_NEAR(w.get_atom_coordinate(1, 0), 2.8345891869307, 1e-9);
		const double mid = 0.5 * w.get_atom_coordinate(1, 0);
		for (const double d : { 0.3, 0.9, 1.6 })
			for (const double y : { 0.0, 0.4 })
			{
				const double left = w.compute_dens(d3{ mid - d, y, 0.1 }), right = w.compute_dens(d3{ mid + d, y, 0.1 });
				EXPECT_NEAR(left, right, 1e-6 * std::max(1.0, left));
			}
		double occ_sum = 0.0;
		for (int m = 0; m < w.get_nmo(); m++)
			occ_sum += w.get_MO_occ(m);
		EXPECT_NEAR(occ_sum, 14.0, 1e-12);
		EXPECT_GT(w.compute_dens(d3{ mid, 0.0, 0.0 }), 0.0);
		expect_orthonormal(w, 1e-8);
	}

	//the spherical F molden against the wfn written from it: the same density at every probe point to the
	//eight significant digits the wfn columns carry
	TEST(WfnReadIoTests, MoldenSphericalMatchesWfnGolden)
	{
		const std::filesystem::path m = nos_test_repo_root() / "tests" / "molden_file" / "F_full.molden";
		const std::filesystem::path g = nos_test_repo_root() / "tests" / "molden_file" / "f_ref.wfn";
		if (!std::filesystem::exists(m) || !std::filesystem::exists(g))
			GTEST_SKIP() << "fixture missing";
		WFN a(m), b(g);
		EXPECT_EQ(a.get_origin(), e_origin::molden);
		EXPECT_EQ(b.get_origin(), e_origin::wfn);
		EXPECT_FALSE(a.get_d_f_switch());
		EXPECT_EQ(a.get_ncen(), 1);
		expect_same_density(a, b, 1e-5);
	}

	//a molden with a Spin= Beta block reads as unrestricted with one operator per spin, and the density matrix
	//sums both spins: 1.0 * 0.6^2 + 0.5 * 0.6^2 on the single basis function
	TEST(WfnReadIoTests, MoldenBetaSpinReadsUnrestricted)
	{
		const std::string text =
			"[Molden Format]\n"
			"[Atoms] AU\n"
			"H     1     1   0.000000   0.000000   0.000000\n"
			"[GTO]\n"
			"  1 0\n"
			" s    1 1.00\n"
			"     1.2000000000      1.0000000000\n"
			"\n"
			"[MO]\n"
			" Sym= A\n"
			" Ene= -0.50\n"
			" Spin= Alpha\n"
			" Occup= 1.0\n"
			"   1   0.6\n"
			" Sym= A\n"
			" Ene= -0.40\n"
			" Spin= Beta\n"
			" Occup= 0.5\n"
			"   1   0.6\n";
		const std::filesystem::path p = tmp_path("beta.molden");
		write_text(p, text);
		WFN w(e_origin::molden);
		std::ostringstream log;
		ASSERT_TRUE(w.read_molden(p, log, false));
		std::filesystem::remove(p);
		EXPECT_TRUE(w.get_is_unrestricted());
		EXPECT_EQ(w.get_MO_op_count(0), 1);
		EXPECT_EQ(w.get_MO_op_count(1), 1);
		const dMatrix2 dm = w.get_dm();
		ASSERT_EQ(dm.extent(0), 1u);
		ASSERT_EQ(dm.extent(1), 1u);
		EXPECT_NEAR(dm(0, 0), 1.5 * 0.36, 1e-12);
	}

	//more MOs than basis functions still reads: the density matrix stays nbf x nbf and only the occupied MO counts
	TEST(WfnReadIoTests, MoldenMoreMosThanFunctionsPadsColumns)
	{
		const std::string text =
			"[Molden Format]\n"
			"[Atoms] AU\n"
			"H     1     1   0.000000   0.000000   0.000000\n"
			"[GTO]\n"
			"  1 0\n"
			" s    1 1.00\n"
			"     1.2000000000      1.0000000000\n"
			"\n"
			"[MO]\n"
			" Ene= -0.50\n"
			" Spin= Alpha\n"
			" Occup= 1.0\n"
			"   1   0.6\n"
			" Ene= 0.10\n"
			" Spin= Alpha\n"
			" Occup= 0.0\n"
			"   1   0.8\n";
		const std::filesystem::path p = tmp_path("wide.molden");
		write_text(p, text);
		WFN w(e_origin::molden);
		std::ostringstream log;
		ASSERT_TRUE(w.read_molden(p, log, false));
		std::filesystem::remove(p);
		EXPECT_EQ(w.get_nmo(), 2);
		EXPECT_NEAR(w.compute_dens(d3{ 0.0, 0.0, 0.0 }), 0.36, 1e-10);
		const dMatrix2 dm = w.get_dm();
		ASSERT_EQ(dm.extent(0), 1u);
		ASSERT_EQ(dm.extent(1), 1u);
		EXPECT_NEAR(dm(0, 0), 0.36, 1e-12);
	}

	//------------------------------------------------------------ tonto ------------------------------------------------------------

	//the four ways into read_tonto for the same OH- job (a file named stdout, the energies file, the orbitals file and
	//explicit partner files) all give the same density, which also matches the wfn written from that job
	TEST(WfnReadIoTests, TontoEntryPathsAgreeWithWfnGolden)
	{
		const std::filesystem::path src = nos_test_repo_root() / "tests" / "cytidine_tonto";
		const std::filesystem::path dir = tmp_path("tonto_oh");
		if (!std::filesystem::exists(src / "stdout_OH") || !std::filesystem::exists(src / "OH.wfn"))
			GTEST_SKIP() << "fixture missing";
		std::filesystem::create_directories(dir);
		const auto copy = std::filesystem::copy_options::overwrite_existing;
		std::filesystem::copy_file(src / "stdout_OH", dir / "stdout", copy);
		std::filesystem::copy_file(src / "OH.orbital_energies,restricted", dir / "OH.orbital_energies,restricted", copy);
		std::filesystem::copy_file(src / "OH.molecular_orbitals,restricted", dir / "OH.molecular_orbitals,restricted", copy);
		std::ostringstream log;
		WFN by_stdout(e_origin::tonto), by_energies(e_origin::tonto), by_orbitals(e_origin::tonto), by_args(e_origin::tonto);
		ASSERT_TRUE(by_stdout.read_tonto(dir / "stdout", log, true));
		ASSERT_TRUE(by_energies.read_tonto(dir / "OH.orbital_energies,restricted", log));
		ASSERT_TRUE(by_orbitals.read_tonto(dir / "OH.molecular_orbitals,restricted", log));
		ASSERT_TRUE(by_args.read_tonto(dir / "stdout", log, false, dir / "OH.orbital_energies,restricted", dir / "OH.molecular_orbitals,restricted"));
		std::filesystem::remove_all(dir);
		EXPECT_NE(log.str().find("File is valid"), std::string::npos);
		EXPECT_EQ(by_stdout.get_method(), "rhf");
		EXPECT_EQ(by_stdout.get_ncen(), 2);
		EXPECT_EQ(by_stdout.get_nr_electrons(), 10u);
		EXPECT_EQ(by_stdout.get_nmo(true), 5);
		EXPECT_FALSE(by_stdout.get_is_unrestricted());
		EXPECT_EQ(by_stdout.get_origin(), e_origin::tonto);
		EXPECT_TRUE(by_stdout.get_isBohr());
		double occ_sum = 0.0;
		for (int m = 0; m < by_stdout.get_nmo(); m++)
			occ_sum += by_stdout.get_MO_occ(m);
		EXPECT_NEAR(occ_sum, 10.0, 1e-12);
		EXPECT_EQ(by_energies.get_nex(), by_stdout.get_nex());
		expect_same_density(by_stdout, by_energies, 1e-10);
		expect_same_density(by_stdout, by_orbitals, 1e-10);
		expect_same_density(by_stdout, by_args, 1e-10);
		WFN golden(src / "OH.wfn");
		expect_same_density(by_stdout, golden, 1e-5);
	}

	//the OH radical as ROHF (rhf, multiplicity 2) gets four doubly and one singly occupied MO; as UHF the stdout entry
	//falls back to the alpha files, swaps to beta for the second spin and gives 5 alpha and 4 beta electrons
	TEST(WfnReadIoTests, TontoRohfAndUhfOccupations)
	{
		const std::filesystem::path src = nos_test_repo_root() / "tests" / "cytidine_tonto";
		if (!std::filesystem::exists(src / "stdout_OH_rad") || !std::filesystem::exists(src / "stdout_OH_rad_uhf"))
			GTEST_SKIP() << "fixture missing";
		std::ostringstream log;
		WFN rohf(e_origin::tonto), uhf(e_origin::tonto);
		ASSERT_TRUE(rohf.read_tonto(src / "stdout_OH_rad", log));
		ASSERT_TRUE(uhf.read_tonto(src / "stdout_OH_rad_uhf", log, true));
		EXPECT_EQ(rohf.get_multi(), 2u);
		EXPECT_EQ(rohf.get_nr_electrons(), 9u);
		EXPECT_FALSE(rohf.get_is_unrestricted());
		EXPECT_EQ(rohf.get_nmo(true), 5);
		EXPECT_NEAR(rohf.get_MO_occ(3), 2.0, 1e-12);
		EXPECT_NEAR(rohf.get_MO_occ(4), 1.0, 1e-12);
		EXPECT_NEAR(rohf.get_MO_occ(5), 0.0, 1e-12);
		EXPECT_TRUE(uhf.get_is_unrestricted());
		EXPECT_EQ(uhf.get_method(), "uhf");
		EXPECT_EQ(uhf.get_nmo(), 2 * rohf.get_nmo());
		EXPECT_EQ(uhf.get_MO_op_count(0), rohf.get_nmo());
		EXPECT_EQ(uhf.get_MO_op_count(1), rohf.get_nmo());
		double alpha = 0.0, beta = 0.0;
		for (int m = 0; m < uhf.get_nmo(); m++)
			(uhf.get_MO_op(m) == 0 ? alpha : beta) += uhf.get_MO_occ(m);
		EXPECT_NEAR(alpha, 5.0, 1e-12);
		EXPECT_NEAR(beta, 4.0, 1e-12);
		EXPECT_NE(log.str().find("al/be els:5 4"), std::string::npos);
		//both describe the same radical, so the total densities agree to the ROHF/UHF spin-polarisation difference;
		//5 % of max(1, rho) is a sanity band that catches a swapped or unscaled beta block, not a numeric check,
		//and the fixture set carries no UHF wfn golden (OH_rad_uhf.wfn is empty) to tighten it against
		expect_same_density(rohf, uhf, 5e-2);
	}

	//------------------------------------------------------------ gbw ------------------------------------------------------------

	//the epoxide gbw (C2H4O, the %basis block of epoxide.inp): O and C carry 3s2p1d = 7 + 12 + 6 = 25 primitives,
	//each H 2s1p = 4 + 3 = 7, 103 in all over 62 contracted functions (pure d), 12 doubly occupied MOs for 24 electrons; ORCA stores every contraction
	//normalised, so all 35 canonical MOs are orthonormal under the analytic primitive overlap; the debug read logs
	//every stage and yields the same wavefunction as the quiet read
	TEST(WfnReadIoTests, GbwReadIsOrthonormalAndDebugMatchesQuiet)
	{
		const std::filesystem::path p = nos_test_repo_root() / "tests" / "epoxide_gbw" / "epoxide.gbw";
		if (!std::filesystem::exists(p))
			GTEST_SKIP() << "fixture missing: " << p;
		std::ostringstream log, quiet_log;
		WFN loud(e_origin::gbw), quiet(e_origin::gbw);
		ASSERT_TRUE(loud.read_gbw(p, log, true));
		ASSERT_TRUE(quiet.read_gbw(p, quiet_log, false));
		EXPECT_NE(log.str().find("I read the geometry of 7 atoms successfully"), std::string::npos);
		EXPECT_NE(log.str().find("I read the coefficients successfully"), std::string::npos);
		EXPECT_NE(log.str().find("I read the energies successfully"), std::string::npos);
		EXPECT_EQ(loud.get_ncen(), 7);
		EXPECT_EQ(loud.get_nex(), 103);
		EXPECT_EQ(loud.get_nmo(), 62);
		EXPECT_EQ(loud.get_nmo(true), 12);
		EXPECT_EQ(loud.get_nr_electrons(), 24u);
		EXPECT_FALSE(loud.get_is_unrestricted());
		double occ_sum = 0.0;
		for (int m = 0; m < loud.get_nmo(); m++)
			occ_sum += loud.get_MO_occ(m);
		EXPECT_NEAR(occ_sum, 24.0, 1e-12);
		for (int m = 1; m < loud.get_nmo(); m++)
			EXPECT_LE(loud.get_MO_energy(m - 1), loud.get_MO_energy(m)) << "MO " << m;
		expect_orthonormal(loud, 1e-8);
		EXPECT_EQ(loud.get_nex(), quiet.get_nex());
		EXPECT_EQ(loud.get_nmo(), quiet.get_nmo());
		EXPECT_EQ(loud.get_nr_basis_set_loaded(), 7);
		for (int m = 0; m < loud.get_nmo(); m++)
			EXPECT_NEAR(loud.get_MO_energy(m), quiet.get_MO_energy(m), 1e-12);
		expect_same_density(loud, quiet, 1e-12);
	}

	//------------------------------------------------------------ occ bridge ------------------------------------------------------------

	//H2 in a one-primitive s + p basis through OCC's SCF (deliberately compact exponents, so the orbital energies are
	//positive, but the SCF is well defined)
	static occ::qm::Wavefunction h2_occ_wavefunction()
	{
		spdlog::set_level(spdlog::level::err);
		const std::vector<occ::core::Atom> atoms{ { 1, 0.0, 0.0, -0.7 }, { 1, 0.0, 0.0, 0.7 } };
		std::vector<occ::gto::Shell> shells;
		for (const auto& at : atoms)
			for (int l = 0; l <= 1; l++)
			{
				shells.emplace_back(l, vec{ l == 0 ? 1.2 : 0.9 }, vec2{ { 1.0 } }, std::array<double, 3>{ at.x, at.y, at.z });
				shells.back().kind = occ::gto::Shell::Kind::Spherical;
				shells.back().incorporate_shell_norm();
			}
		occ::gto::AOBasis basis(atoms, shells, "sp");
		basis.set_pure(true);
		occ::qm::HartreeFock hf(basis);
		occ::qm::SCF<occ::qm::HartreeFock> scf(hf, occ::qm::SpinorbitalKind::Restricted);
		scf.set_charge_multiplicity(0, 1);
		scf.compute_initial_guess();
		scf.compute_scf_energy();
		return scf.wavefunction();
	}

	//OCC's wavefunction into a WFN and back through wfn_to_occ_wavefunction: the rebuilt basis has the same size,
	//the electron bookkeeping and the orbital energies survive and the WFN built from the rebuilt wavefunction has
	//the same density as the first one
	TEST(WfnReadIoTests, OccRoundTripKeepsDensity)
	{
		const occ::qm::Wavefunction wf = h2_occ_wavefunction();
		WFN w(wf, false);
		EXPECT_EQ(w.get_origin(), e_origin::OCC);
		EXPECT_EQ(w.get_ncen(), 2);
		EXPECT_EQ(w.get_nex(), 8);
		EXPECT_EQ(w.get_nmo(), 8);
		EXPECT_EQ(w.get_nmo(true), 1);
		EXPECT_NEAR(w.get_MO_occ(0), 2.0, 1e-12);
		occ::qm::Wavefunction back;
		w.wfn_to_occ_wavefunction(back);
		EXPECT_EQ(back.basis.nbf(), wf.basis.nbf());
		EXPECT_EQ(back.basis.shells().size(), wf.basis.shells().size());
		EXPECT_EQ(back.mo.n_alpha, 1);
		EXPECT_EQ(back.mo.n_beta, 1);
		EXPECT_EQ(back.num_electrons, 2);
		ASSERT_EQ(back.mo.energies.size(), wf.mo.energies.size());
		for (int i = 0; i < wf.mo.energies.size(); i++)
			EXPECT_NEAR(back.mo.energies(i), wf.mo.energies(i), 1e-12);
		WFN again(back, false);
		expect_same_density(w, again, 1e-6);
	}

	//the rebuilt OCC wavefunction carries the density matrix OCC converged, in the same normalised AO basis
	//(the rebuilt shells are normalised like OCC's loaders do, with the dominant contraction coefficient positive)
	TEST(WfnReadIoTests, OccRoundTripKeepsDensityMatrix)
	{
		const occ::qm::Wavefunction wf = h2_occ_wavefunction();
		WFN w(wf, false);
		occ::qm::Wavefunction back;
		w.wfn_to_occ_wavefunction(back);
		ASSERT_EQ(back.mo.D.rows(), wf.mo.D.rows());
		for (int i = 0; i < wf.mo.D.rows(); i++)
			for (int j = 0; j < wf.mo.D.cols(); j++)
				EXPECT_NEAR(back.mo.D(i, j), wf.mo.D(i, j), 1e-6) << i << "," << j;
	}

	//------------------------------------------------------------ death tests ------------------------------------------------------------

	//a wfn header without counts, an unknown extension and a tonto name without a recognised extension all end in
	//err_checkf's exit(-1)
	TEST(WfnReadDeathTest, BadHeaderUnknownExtensionAndTontoNameExit)
	{
		const std::filesystem::path bad = tmp_path("bad.wfn"), odd = tmp_path("odd.abc"), tonto = tmp_path("job.energies");
		write_text(bad, "title\nGAUSSIAN MOL ORBITALS PRIMITIVES NUCLEI\n");
		write_text(odd, "nothing\n");
		write_text(tonto, "nothing\n");
		EXPECT_EXIT({ WFN w(e_origin::wfn); std::ostringstream log; w.read_wfn(bad, false, log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w(odd); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w(e_origin::tonto); std::ostringstream log; w.read_tonto(tonto, log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w(e_origin::wfn); std::ostringstream log; w.read_wfx(tmp_path("missing.wfx"), false, log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		std::filesystem::remove(bad);
		std::filesystem::remove(odd);
		std::filesystem::remove(tonto);
	}

	//out-of-range MO, atom and coefficient accessors and get_norm_const without a basis all exit instead of
	//indexing past the vectors
	TEST(WfnReadDeathTest, OutOfRangeAccessorsExit)
	{
		EXPECT_EXIT({ WFN w = make_h2_wfn(); w.get_MO(3); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w = make_h2_wfn(); w.delete_MO(3); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w = make_h2_wfn(); w.erase_atom(2); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w = make_h2_wfn(); w.get_atom_pos(2); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w = make_h2_wfn(); w.get_MO_coef(3, 0); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w = make_h2_wfn(); w.get_MO_energy(3); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w = make_h2_wfn(); std::ostringstream log; w.get_norm_const(log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w = make_h2_wfn(); w.push_back_atom_basis_set(0, 1.0, 1.0, 1, 0); std::ostringstream log; w.get_norm_const(log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	}

	//an xyz with one atom line more than its count and a molden whose first line is not the format tag both exit
	TEST(WfnReadDeathTest, MalformedXyzAndMoldenExit)
	{
		const std::filesystem::path xyz = tmp_path("extra.xyz"), molden = tmp_path("notmolden.molden");
		write_text(xyz, "1\ncomment\nH 0 0 0\nH 0 0 1\n");
		write_text(molden, "[Atoms] AU\n");
		EXPECT_EXIT({ WFN w(e_origin::xyz); std::ostringstream log; w.read_xyz(xyz, log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		EXPECT_EXIT({ WFN w(e_origin::molden); std::ostringstream log; w.read_molden(molden, log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
		std::filesystem::remove(xyz);
		std::filesystem::remove(molden);
	}

	//every reader stops with an error on a file cut at 50 % and at 90 %: the fuzzing of 18 Sep 2026
	//found access violations (molden, gbw, xtb) and endless getline loops (wfx) on exactly these;
	//a regression shows as a crash code or a hanging death test rather than the err_checkf exit
	TEST(WfnReadDeathTest, TruncatedFilesOfEveryFormatExit)
	{
		const std::filesystem::path root = nos_test_repo_root() / "tests";
		const std::filesystem::path fixtures[] = {
			root / "molden_file" / "f_ref.wfn", root / "molden_file" / "f_ref.wfx", root / "molden_file" / "epoxide.molden",
			root / "epoxide_gbw" / "epoxide.gbw", root / "ptb_H_file" / "wfn.xtb" };
		for (const auto& src : fixtures)
		{
			if (!std::filesystem::exists(src)) { ADD_FAILURE() << "Missing " << src; continue; }
			std::ifstream in(src, std::ios::binary);
			const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
			for (const double keep : { 0.5, 0.9 })
			{
				// a molden MO list is open-ended and the 90 % cut of epoxide.molden ends inside a coefficient
				// ("16  0.5" for 0.5656): a shorter valid file, nothing a text reader can refuse
				if (keep > 0.5 && src.extension() == ".molden") continue;
				const std::filesystem::path cut = tmp_path("truncated" + src.extension().string());
				std::ofstream(cut, std::ios::binary) << bytes.substr(0, static_cast<size_t>(keep * bytes.size()));
				EXPECT_EXIT({ WFN w(cut); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*") << src.filename() << " cut at " << keep;
				std::filesystem::remove(cut);
			}
		}
	}
}
