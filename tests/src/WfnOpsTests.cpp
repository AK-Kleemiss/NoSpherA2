#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/atoms.h"
#include "core/wfn_class.h"
#include "core/cube.h"
#include "core/mo_class.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

namespace
{
	//cartesian type index whose exponents are (a, b, c), independent of the numbering convention
	static int type_of(const int a, const int b, const int c)
	{
		int v[3];
		for (int t = 1; t <= 286; t++)
		{
			constants::type2vector(t, v);
			if (v[0] == a && v[1] == b && v[2] == c)
				return t;
		}
		return -1;
	}

	struct prim_spec
	{
		int cent, l[3];
		double e, c[3];
	};

	//two atoms, three MOs (occ 2, 1, 0), primitives from s up to i so every per-axis l branch runs
	static WFN make_wfn()
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
		w.push_back_atom("H", 1.4, 0.3, -0.2, 1);
		w.push_back_MO(1, 2.0, -0.9);
		w.push_back_MO(2, 1.0, -0.3);
		w.push_back_MO(3, 0.0, 0.5);
		prim_spec specs[] = {
			{1, {0, 0, 0}, 1.20, {0.70, -0.20, 0.10}},
			{1, {1, 0, 0}, 0.80, {0.30, 0.40, -0.20}},
			{1, {0, 2, 0}, 0.90, {-0.20, 0.30, 0.50}},
			{1, {4, 0, 0}, 0.70, {0.15, -0.25, 0.20}},
			{1, {0, 5, 0}, 0.75, {0.10, 0.20, -0.30}},
			{1, {0, 0, 6}, 0.80, {-0.08, 0.12, 0.40}},
			{2, {0, 0, 0}, 1.00, {0.40, 0.50, -0.10}},
			{2, {0, 0, 1}, 0.60, {0.25, -0.35, 0.30}},
			{2, {0, 4, 0}, 0.65, {-0.12, 0.18, 0.20}},
			{2, {0, 0, 5}, 0.70, {0.09, 0.14, 0.10}},
			{2, {6, 0, 0}, 0.75, {0.06, -0.10, 0.30}},
			{2, {1, 1, 1}, 0.55, {0.20, 0.10, -0.40}},
		};
		for (prim_spec& p : specs)
			w.add_primitive(p.cent, type_of(p.l[0], p.l[1], p.l[2]), p.e, p.c);
		w.set_exp_cutoff();
		return w;
	}

	//one He nucleus with a single normalised s orbital of exponent a holding occ electrons
	static WFN make_he_ion(const double a, const double occ)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
		w.push_back_MO(1, occ, -0.9);
		double c = std::pow(2 * a / constants::PI, 0.75);
		w.add_primitive(1, 1, a, &c);
		w.set_exp_cutoff();
		return w;
	}

	static double tau_fd(const WFN& w, const d3& P, const double h)
	{
		double tau = 0;
		for (int mo = 0; mo < w.get_nmo(); mo++)
		{
			double g2 = 0;
			for (int k = 0; k < 3; k++)
			{
				d3 Pp = P, Pm = P;
				Pp[k] += h;
				Pm[k] -= h;
				const double g = (w.computeMO(Pp, mo) - w.computeMO(Pm, mo)) / (2 * h);
				g2 += g * g;
			}
			tau += w.get_MO_occ(mo) * g2;
		}
		return tau;
	}

	static const d3 probe{ 0.6, 0.2, -0.1 };

	//one exp per distinct (centre, exponent): make_wfn has 0.80 twice on He, so 11 groups for 12 primitives;
	//the table reproduces exp(-a r^2) of every primitive and zeroes a whole centre once its most diffuse
	//primitive is below exp_cutoff, exactly where the per-primitive test would have skipped it
	TEST(WfnOpsTests, ExpTableMatchesPerPrimitiveExp)
	{
		WFN w = make_wfn();
		w.get_coef_primitive_major();
		EXPECT_EQ(w.get_exp_group_count(), w.get_nex() - 1);
		vec r2{ 0.3, 2.0 };
		vec ex(w.get_exp_group_count(), -1.0);
		w.exp_table([&r2](const int c) { return r2[c]; }, ex.data());
		for (int j = 0; j < w.get_nex(); j++)
			EXPECT_DOUBLE_EQ(ex[w.get_exp_group(j)], std::exp(-w.get_exponent(j) * r2[w.get_center(j) - 1]));
		r2[1] = -constants::exp_cutoff / 0.55 + 1.0; //past the cutoff even for the most diffuse H primitive
		w.exp_table([&r2](const int c) { return r2[c]; }, ex.data());
		for (int j = 0; j < w.get_nex(); j++)
			if (w.get_center(j) == 2)
				EXPECT_EQ(ex[w.get_exp_group(j)], 0.0);
			else
				EXPECT_GT(ex[w.get_exp_group(j)], 0.0);
	}

	//compute_g_cartesian ignores the MO coefficients: it is nmo times the squared sum of all primitives,
	//so the hand-built primitive sum from type2vector must reproduce it
	TEST(WfnOpsTests, GCartesianIsNmoTimesSquaredPrimitiveSum)
	{
		WFN w = make_wfn();
		double chi_sum = 0;
		int v[3];
		for (int j = 0; j < w.get_nex(); j++)
		{
			constants::type2vector(w.get_type(j), v);
			const int iat = w.get_center(j) - 1;
			double r2 = 0, val = 1;
			for (int k = 0; k < 3; k++)
			{
				const double d = probe[k] - w.get_atom_coordinate(iat, k);
				r2 += d * d;
				val *= std::pow(d, v[k]);
			}
			chi_sum += val * std::exp(-w.get_exponent(j) * r2);
		}
		vec2 d(w.get_ncen(), vec(16, 0.0));
		vec phi(w.get_nmo(), 0.0);
		EXPECT_NEAR(w.compute_g_cartesian(probe, d, phi), w.get_nmo() * chi_sum * chi_sum, 1e-12);
		for (int mo = 0; mo < w.get_nmo(); mo++)
			EXPECT_NEAR(phi[mo], chi_sum, 1e-12);
	}

	//spin density is occ-weighted alpha minus beta orbital density; a restricted set gives the total density
	TEST(WfnOpsTests, SpinDensityIsAlphaMinusBeta)
	{
		WFN w = make_wfn();
		EXPECT_NEAR(w.compute_spin_dens(probe), w.compute_dens(probe), 1e-12);
		WFN u(e_origin::NOT_YET_DEFINED);
		u.push_back_atom("He", 0.0, 0.0, 0.0, 2);
		u.push_back_MO(1, 1.0, -0.9, 0);
		u.push_back_MO(2, 1.0, -0.3, 1);
		u.push_back_MO(3, 1.0, -0.2, 1);
		double c1[] = { 0.7, 0.3, -0.2 }, c2[] = { 0.2, -0.5, 0.4 }, c3[] = { -0.1, 0.4, 0.3 };
		u.add_primitive(1, 1, 1.1, c1);
		u.add_primitive(1, type_of(1, 0, 0), 0.8, c2);
		u.add_primitive(1, type_of(0, 0, 2), 0.9, c3);
		u.set_exp_cutoff();
		const double a = u.computeMO(probe, 0), b1 = u.computeMO(probe, 1), b2 = u.computeMO(probe, 2);
		const double expected = a * a - b1 * b1 - b2 * b2;
		EXPECT_NEAR(u.compute_spin_dens(probe), expected, 1e-12);
		vec2 d(u.get_ncen(), vec(16, 0.0));
		vec phi(u.get_nmo(), 0.0);
		EXPECT_NEAR(u.compute_spin_dens(probe, d, phi), expected, 1e-12);
		EXPECT_NEAR(u.compute_dens(probe), a * a + b1 * b1 + b2 * b2, 1e-12);
		EXPECT_EQ(u.get_MO_op_count(0), 1);
		EXPECT_EQ(u.get_MO_op_count(1), 2);
	}

	//computeValues: density equals compute_dens, Hessian and Laplacian equal central differences of it,
	//normGrad is the reduced gradient alpha |grad rho| / rho^(4/3)
	TEST(WfnOpsTests, ComputeValuesMatchesFiniteDifferences)
	{
		WFN w = make_wfn();
		const double h = 1e-3;
		double rho = 0, ng = 0, elf = 0, eli = 0, lap = 0, hess[9]{};
		w.computeValues(probe, rho, ng, hess, elf, eli, lap);
		const double rho_ref = w.compute_dens(probe);
		EXPECT_NEAR(rho, rho_ref, 1e-12);
		double grad2 = 0, trace = 0;
		for (int k = 0; k < 3; k++)
		{
			d3 Pp = probe, Pm = probe;
			Pp[k] += h;
			Pm[k] -= h;
			const double rp = w.compute_dens(Pp), rm = w.compute_dens(Pm);
			const double g = (rp - rm) / (2 * h);
			grad2 += g * g;
			const double hkk = (rp - 2 * rho_ref + rm) / (h * h);
			EXPECT_NEAR(hess[4 * k], hkk, 1e-4 * (1 + std::abs(hkk))) << "diagonal " << k;
			trace += hkk;
			for (int m = k + 1; m < 3; m++)
			{
				d3 Ppp = probe, Ppm = probe, Pmp = probe, Pmm = probe;
				Ppp[k] += h; Ppp[m] += h;
				Ppm[k] += h; Ppm[m] -= h;
				Pmp[k] -= h; Pmp[m] += h;
				Pmm[k] -= h; Pmm[m] -= h;
				const double hkm = (w.compute_dens(Ppp) - w.compute_dens(Ppm) - w.compute_dens(Pmp) + w.compute_dens(Pmm)) / (4 * h * h);
				EXPECT_NEAR(hess[3 * k + m], hkm, 1e-4 * (1 + std::abs(hkm))) << "off-diagonal " << k << m;
				EXPECT_DOUBLE_EQ(hess[3 * k + m], hess[3 * m + k]);
			}
		}
		EXPECT_NEAR(lap, trace, 1e-4 * (1 + std::abs(trace)));
		//reduced gradient s = |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3)), the prefactor written out so a wrong constant fails
		const double alpha = 1 / (2 * std::cbrt(3 * constants::PI * constants::PI));
		EXPECT_NEAR(ng, alpha * std::sqrt(grad2) / std::pow(rho_ref, 4.0 / 3.0), 1e-5 * (1 + ng));
	}

	//every routine that returns ELF or ELI-D must agree, and both must follow their defining formulas
	//built from finite-difference orbital gradients
	TEST(WfnOpsTests, ElfEliAgreeAcrossVariantsAndFormula)
	{
		WFN w = make_wfn();
		const double h = 1e-3;
		double rho = 0, ng = 0, elf = 0, eli = 0, lap = 0, hess[9]{};
		w.computeValues(probe, rho, ng, hess, elf, eli, lap);
		double elf2 = 0, eli2 = 0;
		w.computeELIELF(probe, elf2, eli2);
		double elf3 = 0, eli3 = 0, lap3 = 0;
		w.computeLapELIELF(probe, elf3, eli3, lap3);
		double eli4 = 0, lap4 = 0;
		w.computeLapELI(probe, eli4, lap4);
		double rho5 = 0, eli5 = 0;
		w.computeRhoELI(probe, rho5, eli5);
		double eli6 = 0;
		d3 eli_grad{};
		w.computeELIGrad(probe, eli6, eli_grad);
		EXPECT_NEAR(elf2, elf, 1e-10);
		EXPECT_NEAR(elf3, elf, 1e-10);
		EXPECT_NEAR(w.computeELF(probe), elf, 1e-10);
		EXPECT_NEAR(eli2, eli, 1e-10);
		EXPECT_NEAR(eli3, eli, 1e-10);
		EXPECT_NEAR(eli4, eli, 1e-10);
		EXPECT_NEAR(eli5, eli, 1e-10);
		EXPECT_NEAR(eli6, eli, 1e-10);
		EXPECT_NEAR(w.computeELI(probe), eli, 1e-10);
		EXPECT_NEAR(rho5, rho, 1e-12);
		double grad2 = 0;
		for (int k = 0; k < 3; k++)
		{
			d3 Pp = probe, Pm = probe;
			Pp[k] += h;
			Pm[k] -= h;
			const double g = (w.compute_dens(Pp) - w.compute_dens(Pm)) / (2 * h);
			grad2 += g * g;
		}
		const double tau = tau_fd(w, probe, 1e-4);
		//ELF = 1 / (1 + (D / D_h)^2) with D = tau/2 - |grad rho|^2 / (8 rho) and the spin-resolved uniform-gas
		//D_h = 2^(2/3) c_F rho^(5/3), c_F = (3/10) (3 pi^2)^(2/3); the prefactor is written out so a wrong constant fails
		const double c_f = 0.3 * std::pow(3 * constants::PI * constants::PI, 2.0 / 3.0);
		const double elf_ref = 1 / (1 + std::pow((0.5 * tau - 0.125 * grad2 / rho) / (std::cbrt(4.0) * c_f * std::pow(rho, 5.0 / 3.0)), 2));
		const double eli_ref = 0.5 * rho * std::pow(48 / (rho * tau - 0.25 * grad2), 3.0 / 8.0);
		EXPECT_NEAR(elf, elf_ref, 1e-4 * (1 + elf_ref));
		EXPECT_NEAR(eli, eli_ref, 1e-4 * (1 + eli_ref));
	}

	//the four Laplacian sources agree with each other and with the second central difference of the density
	TEST(WfnOpsTests, LaplacianVariantsAgree)
	{
		WFN w = make_wfn();
		const double h = 1e-3;
		double rho = 0, ng = 0, elf = 0, eli = 0, lap = 0, hess[9]{};
		w.computeValues(probe, rho, ng, hess, elf, eli, lap);
		double elf3 = 0, eli3 = 0, lap3 = 0, eli4 = 0, lap4 = 0;
		w.computeLapELIELF(probe, elf3, eli3, lap3);
		w.computeLapELI(probe, eli4, lap4);
		double ref = 0;
		for (int k = 0; k < 3; k++)
		{
			d3 Pp = probe, Pm = probe;
			Pp[k] += h;
			Pm[k] -= h;
			ref += (w.compute_dens(Pp) - 2 * rho + w.compute_dens(Pm)) / (h * h);
		}
		EXPECT_NEAR(w.computeLap(probe), ref, 1e-4 * (1 + std::abs(ref)));
		EXPECT_NEAR(lap, ref, 1e-4 * (1 + std::abs(ref)));
		EXPECT_NEAR(lap3, lap, 1e-10);
		EXPECT_NEAR(lap4, lap, 1e-10);
	}

	//analytic density gradient and analytic ELI-D gradient equal central differences of compute_dens and computeELI
	TEST(WfnOpsTests, GradientsMatchFiniteDifferences)
	{
		WFN w = make_wfn();
		const double h = 1e-3;
		d3 grad{}, eli_grad{};
		double eli = 0;
		w.computeGrad(probe, grad);
		w.computeELIGrad(probe, eli, eli_grad);
		for (int k = 0; k < 3; k++)
		{
			d3 Pp = probe, Pm = probe;
			Pp[k] += h;
			Pm[k] -= h;
			const double g = (w.compute_dens(Pp) - w.compute_dens(Pm)) / (2 * h);
			EXPECT_NEAR(grad[k], g, 1e-5 * (1 + std::abs(g))) << "density gradient " << k;
			const double ge = (w.computeELI(Pp) - w.computeELI(Pm)) / (2 * h);
			EXPECT_NEAR(eli_grad[k], ge, 1e-4 * (1 + std::abs(ge))) << "ELI gradient " << k;
		}
	}

	//a normalised s Gaussian of exponent a is a Gaussian charge cloud, so the ESP of a He nucleus with N such
	//electrons is (2 - N erf(sqrt(2a) r)) / r; r = 1 uses the tabulated Boys function, r = 6 (T = 72) the asymptotic branch
	TEST(WfnOpsTests, EspMatchesGaussianChargeCloud)
	{
		const double a = 1.0, occ = 1.0;
		WFN w = make_he_ion(a, occ);
		const WFN::ESP_pairs pairs = w.build_ESP_pairs();
		for (const double r : { 1.0, 6.0 })
		{
			const d3 P{ 0.6 * r, 0.0, 0.8 * r };
			const double ref = (2.0 - occ * std::erf(std::sqrt(2 * a) * r)) / r;
			EXPECT_NEAR(w.computeESP(P, pairs), ref, 1e-7) << "r = " << r;
		}
	}

	//MO::hdr writes the fixed-width wfn MO line: occupation branches for 2, 0 and fractional, energy padding per decade
	TEST(WfnOpsMoTests, HeaderFormatsOccupationAndEnergy)
	{
		EXPECT_EQ(MO(1, 2.0, -0.5).hdr(), "MO    1                  OCC NO =    2.00000000 ORB. ENERGY =   -0.500000\n");
		EXPECT_EQ(MO(12, 0.0, -15.25).hdr(), "MO   12                  OCC NO =    0.00000000 ORB. ENERGY =  -15.250000\n");
		EXPECT_EQ(MO(123, 0.5, -150.0).hdr(), "MO  123                  OCC NO =    0.500000  ORB. ENERGY = -150.000000\n");
		EXPECT_EQ(MO(2, 2.0, 3.0).hdr(), "MO    2                  OCC NO =    2.00000000 ORB. ENERGY =    3.000000\n");
	}

	//erase_coef takes a 1-based index, refuses 0, past-the-end and a wrong primitive count, and drops that entry
	TEST(WfnOpsMoTests, EraseCoefChecksIndexAndCountThenRemoves)
	{
		MO mo(1, 2.0, -1.0);
		for (double c : { 0.1, 0.2, 0.3 })
			mo.push_back_coef(c);
		EXPECT_FALSE(mo.erase_coef(0, 2));
		EXPECT_FALSE(mo.erase_coef(4, 2));
		EXPECT_FALSE(mo.erase_coef(2, 3));
		EXPECT_TRUE(mo.erase_coef(2, 2));
		ASSERT_EQ(mo.get_primitive_count(), 2);
		EXPECT_DOUBLE_EQ(mo.get_coefficient(0), 0.1);
		EXPECT_DOUBLE_EQ(mo.get_coefficient(1), 0.3);
	}

	//the unchecked coefficient accessors, set_coefficient, assign_coefs and get_spin all see the same storage
	TEST(WfnOpsMoTests, CoefficientAccessorsShareStorage)
	{
		MO mo(3, 1.0, -0.2, 1);
		EXPECT_EQ(mo.get_spin(), 1);
		mo.assign_coefs(vec{ 0.5, -0.25 });
		EXPECT_TRUE(mo.set_coefficient(1, 0.75));
		EXPECT_DOUBLE_EQ(mo.get_coefficient_f(1), 0.75);
		EXPECT_DOUBLE_EQ(mo.get_coefficient_ptr()[0], 0.5);
		EXPECT_EQ(mo.get_coefficient_ptr(), mo.get_coefficients().data());
		EXPECT_EQ(mo.get_ptr_coef_vector().size(), 2u);
		mo.assign_coefficients_size(4);
		EXPECT_EQ(mo.get_primitive_count(), 4);
		EXPECT_DOUBLE_EQ(mo.get_coefficient(3), 0.0);
	}

	//alpha and beta electron counts split the occupations by spin, unlike the occupied-MO count
	TEST(WfnOpsMoTests, AlphaBetaCountsSplitBySpin)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("O", 0.0, 0.0, 0.0, 8);
		w.push_back_MO(1, 1.0, -1.0, 0);
		w.push_back_MO(2, 1.0, -0.8, 0);
		w.push_back_MO(3, 0.5, -0.5, 1);
		w.push_back_MO(4, 0.0, 0.5, 1);
		EXPECT_DOUBLE_EQ(w.count_alpha_electrons(), 2.0);
		EXPECT_DOUBLE_EQ(w.count_beta_electrons(), 0.5);
		EXPECT_DOUBLE_EQ(w.count_nr_electrons(), 2.5);
		EXPECT_EQ(w.get_nmo(true), 3);
		EXPECT_EQ(w.get_nmo(false), 4);
		EXPECT_EQ(w.get_MO_op(3), 1);
	}

	//the MO deletion helpers keep nmo and the stored MOs consistent in every path
	TEST(WfnOpsMoTests, DeletionHelpersKeepCountsConsistent)
	{
		WFN w = make_wfn();
		ASSERT_EQ(w.get_nmo(), 3);
		w.push_back_MO(MO(4, 0.0, 1.0));
		EXPECT_EQ(w.get_nmo(), 4);
		EXPECT_DOUBLE_EQ(w.get_MO(3).get_energy(), 1.0);
		w.pop_back_MO();
		EXPECT_EQ(w.get_nmo(), 3);
		w.delete_unoccupied_MOs();
		EXPECT_EQ(w.get_nmo(), 2);
		EXPECT_EQ(w.get_nmo(true), 2);
		EXPECT_DOUBLE_EQ(w.get_MO_energy(1), -0.3);
		w.delete_MO(0);
		EXPECT_EQ(w.get_nmo(), 1);
		EXPECT_DOUBLE_EQ(w.get_MO_occ(0), 1.0);
		EXPECT_DOUBLE_EQ(w.get_MO_energy(0), -0.3);
		w.clear_MOs();
		EXPECT_EQ(w.get_nmo(), 0);
	}

	//the largest coefficient magnitude is taken over occupied MOs by default and over all MOs on request
	TEST(WfnOpsMoTests, MaximumCoefficientRespectsOccupation)
	{
		WFN w = make_wfn();
		EXPECT_DOUBLE_EQ(w.get_maximum_MO_coefficient(true), 0.7);
		EXPECT_DOUBLE_EQ(w.get_maximum_MO_coefficient(false), 0.7);
		double big[] = { 0.0, 0.0, -3.0 };
		w.add_primitive(1, 1, 2.0, big);
		EXPECT_DOUBLE_EQ(w.get_maximum_MO_coefficient(true), 0.7);
		EXPECT_DOUBLE_EQ(w.get_maximum_MO_coefficient(false), 3.0);
	}

	//per-atom basis set accessors read back what was pushed and answer -1 / false to clearly invalid indices
	TEST(WfnOpsAtomTests, BasisSetAccessorsAndBounds)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		EXPECT_TRUE(w.push_back_atom_basis_set(0, 10.0, 0.5, 1, 0));
		EXPECT_TRUE(w.push_back_atom_basis_set(0, 2.0, 0.6, 1, 0));
		EXPECT_TRUE(w.push_back_atom_basis_set(0, 0.8, 1.0, 2, 1));
		EXPECT_FALSE(w.push_back_atom_basis_set(7, 0.8, 1.0, 2, 1));
		EXPECT_EQ(w.get_atom_primitive_count(0), 3);
		EXPECT_EQ(w.get_atom_primitive_count(7), -1);
		EXPECT_EQ(w.get_atom_basis_set_size(0), 3);
		EXPECT_DOUBLE_EQ(w.get_atom_basis_set_exponent(0, 1), 2.0);
		EXPECT_DOUBLE_EQ(w.get_atom_basis_set_coefficient(0, 2), 1.0);
		EXPECT_DOUBLE_EQ(w.get_atom_basis_set_exponent(0, -1), -1.0);
		EXPECT_DOUBLE_EQ(w.get_atom_basis_set_coefficient(7, 0), -1.0);
		EXPECT_TRUE(w.change_atom_basis_set_exponent(0, 1, 2.5));
		EXPECT_FALSE(w.change_atom_basis_set_exponent(7, 1, 2.5));
		EXPECT_TRUE(w.change_atom_basis_set_coefficient(0, 2, 0.9));
		EXPECT_DOUBLE_EQ(w.get_atom_basis_set_entry(0, 1).get_exponent(), 2.5);
		EXPECT_DOUBLE_EQ(w.get_atom_basis_set_entry(0, 2).get_coefficient(), 0.9);
		EXPECT_TRUE(w.get_modified());
		EXPECT_EQ(w.get_basis_set_shell(0, 2), 1);
		EXPECT_EQ(w.get_basis_set_shell(7, 0), -1);
		EXPECT_EQ(w.get_atom_primitive_type(0, 2), 2);
		EXPECT_EQ(w.get_atom_primitive_type(0, 3), -1);
		EXPECT_EQ(w.get_nr_basis_set_loaded(), 1);
		EXPECT_TRUE(w.get_atom_basis_set_loaded(0));
		EXPECT_FALSE(w.get_atom_basis_set_loaded(-1));
		EXPECT_FALSE(w.erase_atom_primitive(0, 5));
		EXPECT_TRUE(w.erase_atom_primitive(0, 0));
		EXPECT_DOUBLE_EQ(w.get_atom_basis_set_exponent(0, 0), 2.5);
		EXPECT_TRUE(w.delete_basis_set());
		EXPECT_EQ(w.get_atom_primitive_count(0), 0);
		EXPECT_EQ(w.get_nr_basis_set_loaded(), 0);
	}

	//shell bookkeeping: starts and ends inside the atom's basis list, the start in the primitive list counts
	//3 entries per p shell, the shell centre is the centre of that first primitive
	TEST(WfnOpsAtomTests, ShellAccessorsCountPrimitivesPerShell)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		w.push_back_atom("H", 2.0, 0.0, 0.0, 1);
		w.push_back_MO(1, 2.0, -1.0);
		w.push_back_atom_basis_set(0, 10.0, 0.5, 1, 0);
		w.push_back_atom_basis_set(0, 2.0, 0.6, 1, 0);
		w.push_back_atom_basis_set(0, 0.8, 1.0, 2, 1);
		w.push_back_atom_basis_set(1, 1.2, 1.0, 1, 0);
		double c = 0.1;
		for (double e : { 10.0, 2.0 })
			w.add_primitive(1, 1, e, &c);
		for (int t = 2; t <= 4; t++)
			w.add_primitive(1, t, 0.8, &c);
		w.add_primitive(2, 1, 1.2, &c);
		EXPECT_EQ(w.get_atom_shell_count(0), 2);
		EXPECT_EQ(w.get_atom_shell_count(1), 1);
		EXPECT_EQ(w.get_atom_shell_count(9), -1);
		EXPECT_EQ(w.get_atom_shell_primitives(0, 0), 2);
		EXPECT_EQ(w.get_atom_shell_primitives(0, 1), 1);
		EXPECT_EQ(w.get_shell_type(0, 0), 1);
		EXPECT_EQ(w.get_shell_type(0, 1), 2);
		EXPECT_EQ(w.get_shell_type(0, 9), -1);
		EXPECT_EQ(w.get_shell_start(0, 0), 0);
		EXPECT_EQ(w.get_shell_start(0, 1), 2);
		EXPECT_EQ(w.get_shell_start(0, 9), -1);
		EXPECT_EQ(w.get_shell_end(0, 0), 1);
		EXPECT_EQ(w.get_shell_end(0, 1), 2);
		EXPECT_EQ(w.get_shell_end(0, 9), -1);
		EXPECT_EQ(w.get_shell_start_in_primitives(0, 1), 2);
		EXPECT_EQ(w.get_shell_start_in_primitives(1, 0), 5);
		EXPECT_EQ(w.get_shell_start_in_primitives(0, 9), -1);
		EXPECT_EQ(w.get_shell_center(0, 1), 1);
		EXPECT_EQ(w.get_shell_center(1, 0), 2);
		EXPECT_EQ(w.get_shell_center(0, 9), -1);
		EXPECT_EQ(w.get_atom_shell_count(9), -1);
	}

	//labels, charges and masses come from the atom list; out-of-range atoms give "?", -1 and 0
	TEST(WfnOpsAtomTests, LabelsChargesAndMasses)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C1", 0.0, 0.0, 0.0, 6);
		w.push_back_atom("Ra", 1.0, 0.0, 0.0, 88);
		EXPECT_EQ(w.get_atom_label(0u), "C1");
		EXPECT_EQ(w.get_atom_label(5u), "?");
		EXPECT_EQ(w.get_atom_label(1), "Ra");
		EXPECT_EQ(w.get_atom_charge(0), 6);
		EXPECT_EQ(w.get_atom_charge(-1), -1);
		EXPECT_EQ(w.get_atom_integer_mass(0), 12u);
		EXPECT_NEAR(w.get_atom_real_mass(0), 12.011, 1e-6);
		EXPECT_EQ(w.get_atom_integer_mass(1), 0u);
		EXPECT_DOUBLE_EQ(w.get_atom_real_mass(1), 0.0);
		w.set_atom_label(0, "C7");
		EXPECT_EQ(w.get_atom(0).get_label(), "C7");
		EXPECT_EQ(w.get_atom(-1).get_charge(), 0);
		EXPECT_EQ(w.get_atoms().size(), 2u);
		EXPECT_EQ(w.get_atoms_ptr()->size(), 2u);
		w.print_atom_long(0);
		atom copy = w.get_atom(0);
		EXPECT_TRUE(w.get_id_for_atom(0) == copy.get_ID());
		EXPECT_TRUE(w.erase_atom(1));
		EXPECT_EQ(w.get_ncen(), 1);
		EXPECT_EQ(w.get_atoms().size(), 1u);
	}

	//the three ECP tables fill core electron counts per element, set_ECPs overrides by atomic number
	TEST(WfnOpsAtomTests, EcpModesFillCoreElectrons)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("Rb", 0.0, 0.0, 0.0, 37);
		w.push_back_atom("C", 1.0, 0.0, 0.0, 6);
		w.set_has_ECPs(true, true, 1);
		EXPECT_TRUE(w.get_has_ECPs());
		EXPECT_EQ(w.get_ECP_mode(), 1);
		EXPECT_EQ(w.get_atom_ECP_electrons(0), 28);
		EXPECT_EQ(w.get_atom_ECP_electrons(1), 0);
		w.set_has_ECPs(true, true, 2);
		EXPECT_EQ(w.get_atom_ECP_electrons(0), 36);
		EXPECT_EQ(w.get_atom_ECP_electrons(1), 2);
		EXPECT_EQ(w.get_nr_ECP_electrons(), 38u);
		//pTB: Rb keeps the 28-electron def2 core, C gets its 1s pair
		w.set_has_ECPs(true, true, 3);
		EXPECT_EQ(w.get_atom_ECP_electrons(0), 28);
		EXPECT_EQ(w.get_atom_ECP_electrons(1), 2);
		ivec nr{ 6 }, els{ 4 };
		w.set_ECPs(nr, els);
		EXPECT_EQ(w.get_atom_ECP_electrons(1), 4);
		EXPECT_EQ(w.get_atom_ECP_electrons(0), 28);
		w.set_has_ECPs(false, false);
		EXPECT_FALSE(w.get_has_ECPs());
	}

	//the tables stop at Rn: an actinide used to index past the end of one and come back with
	//Z = -543649293 core electrons, and the run segfaulted rather than saying anything
	TEST(WfnOpsAtomTests, EcpTablesEndAtRadonAndAnActinideGetsNoCore)
	{
		EXPECT_EQ(constants::heaviest_ECP_element, 86);
		EXPECT_EQ(constants::ECP_core_electrons(constants::ECP_electrons, 86), 60);
		for (const int Z : { 87, 90, 92, 103, 118, 1000, -1 })
		{
			EXPECT_EQ(constants::ECP_core_electrons(constants::ECP_electrons, Z), 0) << "Z = " << Z;
			EXPECT_EQ(constants::ECP_core_electrons(constants::ECP_electrons_xTB, Z), 0) << "Z = " << Z;
			EXPECT_EQ(constants::ECP_core_electrons(constants::ECP_electrons_pTB, Z), 0) << "Z = " << Z;
		}
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("U", 0.0, 0.0, 0.0, 92);
		w.push_back_atom("H", 1.0, 0.0, 0.0, 1);
		std::ostringstream captured;
		std::streambuf *const previous = std::cout.rdbuf(captured.rdbuf());
		w.set_has_ECPs(true, true, 1);
		std::cout.rdbuf(previous);
		EXPECT_EQ(w.get_atom_ECP_electrons(0), 0);
		EXPECT_EQ(w.get_atom_ECP_electrons(1), 0);
		//and it says so, rather than quietly counting 92 nuclear charges against a valence basis
		EXPECT_NE(captured.str().find("Z = 92"), std::string::npos) << captured.str();
	}

	//charge = sum Z - sum occ, -1000 for an untyped atom; the multiplicity guess is 1 for even and 2 for odd electron counts
	TEST(WfnOpsAtomTests, ChargeAndMultiplicityFromOccupations)
	{
		WFN w = make_wfn();
		std::ostringstream log;
		EXPECT_EQ(w.calculate_charge(), 0);
		EXPECT_EQ(w.calculate_charge(log), 0);
		w.assign_charge(1);
		EXPECT_EQ(w.get_nr_electrons(), 2u);
		EXPECT_TRUE(w.guess_multiplicity(log));
		EXPECT_EQ(w.get_multi(), 1);
		w.assign_charge(0);
		EXPECT_TRUE(w.guess_multiplicity(log));
		EXPECT_EQ(w.get_multi(), 2);
		EXPECT_NE(log.str().find("multiplicity 1"), std::string::npos);
		EXPECT_NE(log.str().find("multiplicity 2"), std::string::npos);
		w.push_back_atom(atom());
		EXPECT_EQ(w.calculate_charge(), -1000);
		EXPECT_EQ(w.calculate_charge(log), -1000);
		EXPECT_NE(log.str().find("Atomtype misunderstanding"), std::string::npos);
	}

	//delete_Qs removes the Z = 119 dummy atoms and must shift only the centres of atoms behind them
	TEST(WfnOpsAtomTests, DeleteQsShiftsOnlyLaterCentres)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		w.push_back_atom("Q", 1.0, 0.0, 0.0, 119);
		w.push_back_atom("H", 2.0, 0.0, 0.0, 1);
		w.push_back_MO(1, 2.0, -1.0);
		double c = 0.5;
		w.add_primitive(1, 1, 1.0, &c);
		w.add_primitive(3, 1, 1.0, &c);
		w.delete_Qs();
		EXPECT_EQ(w.get_ncen(), 2);
		EXPECT_EQ(w.get_atom_charge(1), 1);
		EXPECT_EQ(w.get_center(0), 1);
		EXPECT_EQ(w.get_center(1), 2);
	}

	//a basis set per atom in gaussian, tonto or ORCA p order is recognised; missing basis sets and
	//d shells without a preceding p order are refused with -1
	TEST(WfnOpsOrderTests, CheckOrderRecognisesGaussianTontoAndOrca)
	{
		struct
		{
			ivec types;
			int expected;
		} cases[] = {
			{ { 1, 2, 2, 3, 3, 4, 4, 5, 6, 7, 8, 9, 10 }, 1 },
			{ { 1, 2, 3, 4, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, 2 },
			{ { 1, 4, 2, 3, 4, 2, 3, 5, 6, 7, 8, 9, 10 }, 3 },
			{ { 1, 3, 4, 2, 3, 4, 2, 5, 6, 7, 8, 9, 10 }, -1 },
		};
		for (const auto& cs : cases)
		{
			WFN w(e_origin::NOT_YET_DEFINED);
			w.push_back_atom("C", 0.0, 0.0, 0.0, 6);
			w.push_back_MO(1, 2.0, -1.0);
			w.push_back_atom_basis_set(0, 5.0, 1.0, 1, 0);
			w.push_back_atom_basis_set(0, 1.0, 1.0, 2, 1);
			w.push_back_atom_basis_set(0, 0.5, 1.0, 2, 1);
			w.push_back_atom_basis_set(0, 0.8, 1.0, 3, 2);
			double c = 0.1;
			for (int t : cs.types)
				w.add_primitive(1, t, 1.0, &c);
			EXPECT_EQ(w.check_order(true), cs.expected) << "first p type " << cs.types[1];
		}
		WFN nobasis(e_origin::NOT_YET_DEFINED);
		nobasis.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		EXPECT_EQ(nobasis.check_order(false), -1);
		WFN d_only(e_origin::NOT_YET_DEFINED);
		d_only.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		d_only.push_back_MO(1, 2.0, -1.0);
		d_only.push_back_atom_basis_set(0, 5.0, 1.0, 1, 0);
		d_only.push_back_atom_basis_set(0, 0.8, 1.0, 3, 1);
		double c = 0.1;
		for (int t : { 1, 5, 6, 7, 8, 9, 10 })
			d_only.add_primitive(1, t, 1.0, &c);
		EXPECT_EQ(d_only.check_order(false), -1);
	}

	//gaussian s/p/d order with two-primitive shells, coefficients distinct per primitive
	static WFN make_gaussian_ordered(const bool with_f)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		w.push_back_MO(1, 2.0, -1.0);
		w.push_back_MO(2, 1.0, -0.5);
		w.push_back_atom_basis_set(0, 5.0, 1.0, 1, 0);
		w.push_back_atom_basis_set(0, 1.0, 1.0, 2, 1);
		w.push_back_atom_basis_set(0, 0.5, 1.0, 2, 1);
		w.push_back_atom_basis_set(0, 0.9, 1.0, 3, 2);
		w.push_back_atom_basis_set(0, 0.4, 1.0, 3, 2);
		ivec types{ 1 };
		vec exps{ 5.0 };
		for (int t = 2; t <= 4; t++)
			for (double e : { 1.0, 0.5 })
				types.push_back(t), exps.push_back(e);
		for (int t = 5; t <= 10; t++)
			for (double e : { 0.9, 0.4 })
				types.push_back(t), exps.push_back(e);
		if (with_f)
		{
			w.push_back_atom_basis_set(0, 0.7, 1.0, 4, 3);
			for (int t : { 11, 12, 13, 17, 14, 15, 18, 19, 16, 20 })
				types.push_back(t), exps.push_back(0.7);
		}
		for (size_t j = 0; j < types.size(); j++)
		{
			double c[2]{ 0.1 * (j + 1), -0.03 * (j + 1) };
			w.add_primitive(1, types[j], exps[j], c);
		}
		w.set_exp_cutoff();
		return w;
	}

	//each new primitive must carry the coefficients its (type, exponent) had before sorting
	static void expect_coefficients_follow_primitives(const WFN& before, const WFN& after)
	{
		ASSERT_EQ(before.get_nex(), after.get_nex());
		for (int i = 0; i < after.get_nex(); i++)
		{
			int match = -1;
			for (int j = 0; j < before.get_nex(); j++)
				if (before.get_type(j) == after.get_type(i) && before.get_exponent(j) == after.get_exponent(i))
					match = j;
			ASSERT_GE(match, 0) << "primitive " << i;
			for (int mo = 0; mo < after.get_nmo(); mo++)
				EXPECT_DOUBLE_EQ(after.get_MO_coef(mo, i), before.get_MO_coef(mo, match)) << "primitive " << i << " MO " << mo;
		}
	}

	//sort_wfn(1) interleaves gaussian p and d shells into tonto order and moves the coefficients along,
	//so the density is unchanged and check_order sees tonto afterwards
	TEST(WfnOpsOrderTests, SortGaussianToTontoKeepsDensity)
	{
		WFN before = make_gaussian_ordered(false);
		WFN w = before;
		ASSERT_EQ(w.check_order(false), 1);
		EXPECT_TRUE(w.sort_wfn(1, true));
		const ivec expected{ 1, 2, 3, 4, 2, 3, 4, 5, 6, 7, 8, 9, 10, 5, 6, 7, 8, 9, 10 };
		for (size_t j = 0; j < expected.size(); j++)
			EXPECT_EQ(w.get_type((int)j), expected[j]) << "primitive " << j;
		EXPECT_DOUBLE_EQ(w.get_exponent(1), 1.0);
		EXPECT_DOUBLE_EQ(w.get_exponent(4), 0.5);
		EXPECT_DOUBLE_EQ(w.get_exponent(7), 0.9);
		EXPECT_DOUBLE_EQ(w.get_exponent(13), 0.4);
		expect_coefficients_follow_primitives(before, w);
		EXPECT_EQ(w.check_order(false), 2);
		for (const d3& P : { d3{ 0.3, -0.2, 0.4 }, d3{ -0.5, 0.1, 0.2 } })
			EXPECT_NEAR(w.compute_dens(P), before.compute_dens(P), 1e-12);
	}

	//sort_wfn(3) rotates every ORCA p triple z x y into x y z with its coefficients, again density preserving
	TEST(WfnOpsOrderTests, SortOrcaRotatesPShellKeepsDensity)
	{
		WFN before(e_origin::NOT_YET_DEFINED);
		before.push_back_atom("N", 0.0, 0.0, 0.0, 7);
		before.push_back_MO(1, 2.0, -1.0);
		before.push_back_atom_basis_set(0, 5.0, 1.0, 1, 0);
		before.push_back_atom_basis_set(0, 1.0, 1.0, 2, 1);
		before.push_back_atom_basis_set(0, 0.5, 1.0, 2, 1);
		const ivec types{ 1, 4, 2, 3, 4, 2, 3 };
		const vec exps{ 5.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5 };
		for (size_t j = 0; j < types.size(); j++)
		{
			double c = 0.1 * (j + 1);
			before.add_primitive(1, types[j], exps[j], &c);
		}
		before.set_exp_cutoff();
		WFN w = before;
		ASSERT_EQ(w.check_order(false), 3);
		EXPECT_TRUE(w.sort_wfn(3, false));
		const ivec expected{ 1, 2, 3, 4, 2, 3, 4 };
		for (size_t j = 0; j < expected.size(); j++)
			EXPECT_EQ(w.get_type((int)j), expected[j]) << "primitive " << j;
		expect_coefficients_follow_primitives(before, w);
		EXPECT_EQ(w.check_order(false), 2);
		const d3 P{ 0.3, -0.2, 0.4 };
		EXPECT_NEAR(w.compute_dens(P), before.compute_dens(P), 1e-12);
		EXPECT_TRUE(w.sort_wfn(2, true));
		EXPECT_TRUE(w.sort_wfn(42, true));
		EXPECT_FALSE(w.sort_wfn(5, false));
		EXPECT_FALSE(w.sort_wfn(52, false));
	}

	//sort_wfn(11) puts a gaussian f shell into natural 11..20 order with its coefficients
	TEST(WfnOpsOrderTests, SortGaussianFShellKeepsDensity)
	{
		WFN before = make_gaussian_ordered(true);
		WFN w = before;
		ASSERT_EQ(w.check_order(false), 11);
		EXPECT_TRUE(w.sort_wfn(11, false));
		for (int t = 11; t <= 20; t++)
			EXPECT_EQ(w.get_type(19 + t - 11), t) << "f primitive " << t;
		expect_coefficients_follow_primitives(before, w);
		const d3 P{ 0.3, -0.2, 0.4 };
		EXPECT_NEAR(w.compute_dens(P), before.compute_dens(P), 1e-12);
	}

	//the xyz writer emits the atom count, a comment and one label + three coordinates per atom in Angstrom
	TEST(WfnOpsIoTests, WriteXyzRoundTrip)
	{
		WFN w = make_wfn();
		const std::filesystem::path tmp = std::filesystem::temp_directory_path() / "nosphera2_wfnops_roundtrip.xyz";
		ASSERT_TRUE(w.write_xyz(tmp));
		std::ifstream in(tmp);
		int n = 0;
		std::string line;
		in >> n;
		std::getline(in, line);
		std::getline(in, line);
		EXPECT_EQ(n, 2);
		EXPECT_NE(line.find("NoSpherA2"), std::string::npos);
		for (int i = 0; i < n; i++)
		{
			std::string label;
			double x = 0, y = 0, z = 0;
			in >> label >> x >> y >> z;
			EXPECT_EQ(label, w.get_atom_label(i));
			const double f = w.get_isBohr() ? constants::bohr2ang(1.0) : 1.0;
			EXPECT_NEAR(x, f * w.get_atom_coordinate(i, 0), 1e-7);
			EXPECT_NEAR(y, f * w.get_atom_coordinate(i, 1), 1e-7);
			EXPECT_NEAR(z, f * w.get_atom_coordinate(i, 2), 1e-7);
		}
		in.close();
		std::filesystem::remove(tmp);
	}

	//a written wfn reads back with the same primitives and density; occupied = true drops the empty MO
	TEST(WfnOpsIoTests, WriteWfnRoundTripKeepsDensity)
	{
		WFN w = make_wfn();
		const std::filesystem::path tmp = std::filesystem::temp_directory_path() / "nosphera2_wfnops_roundtrip.wfn";
		ASSERT_TRUE(w.write_wfn(tmp, false, true));
		WFN back(tmp, false);
		EXPECT_EQ(back.get_ncen(), w.get_ncen());
		EXPECT_EQ(back.get_nex(), w.get_nex());
		EXPECT_EQ(back.get_nmo(), w.get_nmo(true));
		for (int j = 0; j < w.get_nex(); j++)
		{
			EXPECT_EQ(back.get_type(j), w.get_type(j));
			EXPECT_EQ(back.get_center(j), w.get_center(j));
		}
		EXPECT_NEAR(back.compute_dens(probe), w.compute_dens(probe), 1e-6 * (1 + w.compute_dens(probe)));
		std::filesystem::remove(tmp);
		EXPECT_FALSE(w.write_wfn(std::filesystem::temp_directory_path() / "nosphera2_wfnops_missing_dir" / "x.wfn", false, true));
	}

	//a written wfx reads back with the same density, and an unwritable path is reported as false
	TEST(WfnOpsIoTests, WriteWfxRoundTripAndFailure)
	{
		WFN w = make_wfn();
		MO extra(4, 0.0, 0.7);
		extra.assign_coefficients_size(w.get_nex());
		w.push_back_MO(extra);
		for (int j = 0; j < w.get_nex(); j++)
			EXPECT_TRUE(w.set_MO_coef(3, j, 0.01 * j));
		const std::filesystem::path tmp = std::filesystem::temp_directory_path() / "nosphera2_wfnops_roundtrip.wfx";
		ASSERT_TRUE(w.write_wfx(tmp, true));
		WFN back(tmp, false);
		EXPECT_EQ(back.get_ncen(), w.get_ncen());
		EXPECT_EQ(back.get_nmo(), w.get_nmo(true));
		EXPECT_NEAR(back.compute_dens(probe), w.compute_dens(probe), 1e-6 * (1 + w.compute_dens(probe)));
		std::filesystem::remove(tmp);
		EXPECT_FALSE(w.write_wfx(std::filesystem::temp_directory_path() / "nosphera2_wfnops_missing_dir" / "x.wfx", true));
	}

	//cubes attached to a wavefunction: write through the wrapper, read back by path, elementwise operations
	//through the wrappers, threshold and masks
	TEST(WfnOpsIoTests, CubeWrappersRoundTripAndOperate)
	{
		WFN w = make_wfn();
		const std::filesystem::path tmp = std::filesystem::temp_directory_path() / "nosphera2_wfnops_cube.cube";
		cube c({ 2, 3, 4 }, 0, true);
		for (int i = 0; i < 3; i++)
		{
			c.set_origin(i, -1.0 + 0.5 * i);
			c.set_vector(i, i, 0.5);
		}
		c.set_comment1("wfnops");
		c.set_comment2("cube");
		for (int x = 0; x < 2; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 4; z++)
					c.set_value(x, y, z, 1.0 + x + 10 * y + 100 * z);
		w.push_back_cube(c);
		EXPECT_EQ(w.get_cube_count(), 1);
		EXPECT_TRUE(w.get_cube_loaded(0));
		w.write_cube_file(0, tmp, false);
		EXPECT_EQ(w.get_cube_path(0), tmp);
		ASSERT_TRUE(w.push_back_cube(tmp.string(), true));
		EXPECT_EQ(w.get_cube_count(), 2);
		EXPECT_EQ(w.get_cube(1).get_size(2), 4);
		EXPECT_NEAR(w.get_cube(1).get_value(1, 2, 3), 322.0, 1e-6);
		EXPECT_NEAR(w.get_cube_ptr(1)->get_origin(1), -0.5, 1e-6);
		EXPECT_TRUE(w.cube_add(0, w.get_cube(1)));
		EXPECT_NEAR(w.get_cube(0).get_value(1, 2, 3), 644.0, 1e-6);
		EXPECT_TRUE(w.cube_subtract(0, w.get_cube(1)));
		EXPECT_TRUE(w.cube_multiply(0, w.get_cube(1)));
		EXPECT_NEAR(w.get_cube(0).get_value(0, 0, 1), 101.0 * 101.0, 1e-6);
		EXPECT_TRUE(w.cube_divide(0, w.get_cube(1)));
		EXPECT_NEAR(w.get_cube(0).get_value(0, 0, 1), 101.0, 1e-6);
		EXPECT_TRUE(w.apply_cube_thresh(0, 300.0));
		EXPECT_DOUBLE_EQ(w.get_cube(0).get_value(0, 0, 1), 0.0);
		EXPECT_NEAR(w.get_cube(0).get_value(0, 0, 3), 301.0, 1e-6);
		EXPECT_TRUE(w.apply_cube_negative_mask(1, w.get_cube(0)));
		EXPECT_DOUBLE_EQ(w.get_cube(1).get_value(0, 0, 3), 0.0);
		EXPECT_NEAR(w.get_cube(1).get_value(0, 0, 1), 101.0, 1e-6);
		EXPECT_TRUE(w.apply_cube_mask(1, w.get_cube(0)));
		EXPECT_DOUBLE_EQ(w.get_cube(1).get_value(0, 0, 1), 0.0);
		EXPECT_TRUE(w.apply_cube_thresh(0, w.get_cube(1), 1.0));
		EXPECT_DOUBLE_EQ(w.get_cube(0).get_value(0, 0, 3), 0.0);
		w.pop_back_cube();
		EXPECT_EQ(w.get_cube_count(), 1);
		ASSERT_TRUE(w.push_back_cube(tmp.string(), false));
		EXPECT_FALSE(w.get_cube_loaded(1));
		EXPECT_TRUE(w.read_cube(1, true, true));
		EXPECT_NEAR(w.get_cube(1).get_value(1, 2, 3), 322.0, 1e-6);
		std::filesystem::remove(tmp);
	}

	//the density cube of a normalised s orbital integrates to its occupation, both through calc_rho_cube on a
	//hand-made grid and through write_rho_cube, whose file lands next to the wavefunction path
	TEST(WfnOpsIoTests, RhoCubeIntegratesToElectronCount)
	{
		const double occ = 1.5;
		WFN w = make_he_ion(1.0, occ);
		const int n = 20;
		cube c({ n, n, n }, 0, true);
		for (int i = 0; i < 3; i++)
		{
			c.set_origin(i, -4.5);
			c.set_vector(i, i, 9.0 / n);
		}
		c.calc_dv();
		w.calc_rho_cube(c);
		EXPECT_NEAR(c.sum(), occ, 1e-4);
		EXPECT_NEAR(c.get_value(n / 2, n / 2, n / 2), w.compute_dens(c.get_pos(n / 2, n / 2, n / 2)), 1e-12);
		const std::filesystem::path base = std::filesystem::temp_directory_path() / "nosphera2_wfnops_rho.wfn";
		const std::filesystem::path out = std::filesystem::temp_directory_path() / "nosphera2_wfnops_rho_rho.cube";
		w.set_path(base);
		w.write_rho_cube(2.5, 0.25);
		ASSERT_TRUE(std::filesystem::exists(out));
		cube back(out, true, w, std::cout);
		EXPECT_NEAR(back.sum(), occ, 1e-3);
		std::filesystem::remove(out);
	}

	//print_primitive lists the coefficient of that primitive in every MO
	TEST(WfnOpsIoTests, PrintPrimitiveListsCoefficientsOfThatPrimitive)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
		w.push_back_MO(1, 2.0, -1.0);
		w.push_back_MO(2, 0.0, 0.5);
		double c0[] = { 0.111, 0.222 }, c1[] = { 0.333, 0.444 };
		w.add_primitive(1, 1, 1.0, c0);
		w.add_primitive(1, 1, 0.5, c1);
		std::ostringstream captured;
		std::streambuf* old = std::cout.rdbuf(captured.rdbuf());
		w.print_primitive(1);
		std::cout.rdbuf(old);
		EXPECT_NE(captured.str().find("0.333"), std::string::npos);
		EXPECT_NE(captured.str().find("0.444"), std::string::npos);
		EXPECT_EQ(captured.str().find("0.222"), std::string::npos);
	}

	//density and spin density matrix containers: push, resize, set, get and their out-of-range answers
	TEST(WfnOpsIoTests, DensityMatrixContainers)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_DM(1.5);
		w.resize_DM(3, 0.25);
		EXPECT_EQ(w.get_DM_size(), 3);
		EXPECT_TRUE(w.set_DM(2, -2.0));
		EXPECT_FALSE(w.set_DM(3, 1.0));
		EXPECT_DOUBLE_EQ(w.get_DM(0), 1.5);
		EXPECT_DOUBLE_EQ(w.get_DM(1), 0.25);
		EXPECT_DOUBLE_EQ(w.get_DM(2), -2.0);
		EXPECT_DOUBLE_EQ(w.get_DM(-1), -1.0);
		EXPECT_DOUBLE_EQ(w.get_DM(3), -1.0);
		w.push_back_SDM(0.5);
		w.resize_SDM(2, 0.75);
		EXPECT_EQ(w.get_SDM_size(), 2);
		EXPECT_TRUE(w.set_SDM(1, 4.0));
		EXPECT_FALSE(w.set_SDM(-1, 4.0));
		EXPECT_DOUBLE_EQ(w.get_SDM(0), 0.5);
		EXPECT_DOUBLE_EQ(w.get_SDM(1), 4.0);
		EXPECT_DOUBLE_EQ(w.get_SDM(2), -1.0);
	}

	//the small inline setters and flags on WFN store what they are given and keep nex in step with the arrays
	TEST(WfnOpsIoTests, InlineFlagsAndArraySetters)
	{
		WFN w = make_wfn();
		w.change_basis_set_name("def2-SVP");
		EXPECT_EQ(w.get_basis_set_name(), "def2-SVP");
		w.set_dist_switch(true);
		EXPECT_TRUE(w.get_dist_switch());
		w.set_dist_switch(false);
		EXPECT_FALSE(w.get_dist_switch());
		w.set_dist_switch();
		EXPECT_TRUE(w.get_dist_switch());
		w.set_d_f_switch(true);
		EXPECT_TRUE(w.get_d_f_switch());
		EXPECT_FALSE(w.get_modified());
		w.set_modified();
		EXPECT_TRUE(w.get_modified());
		w.set_origin(e_origin::wfn);
		EXPECT_EQ(w.get_origin(), e_origin::wfn);
		w.set_method("hf");
		EXPECT_EQ(w.get_method(), "hf");
		WFN u(e_origin::NOT_YET_DEFINED);
		u.push_back_atom("He", 0.0, 0.0, 0.0, 2);
		u.set_center(ivec{ 1, 1 });
		EXPECT_EQ(u.get_nex(), 2);
		u.set_types(ivec{ 1, 2 });
		u.set_exponents(vec{ 1.5, 0.5 });
		EXPECT_EQ(u.get_types()[1], 2);
		EXPECT_DOUBLE_EQ(u.get_exponents()[0], 1.5);
		EXPECT_EQ(u.get_centers()[1], 1);
		EXPECT_EQ(u.get_ptr_types()[1], 2);
		EXPECT_DOUBLE_EQ(u.get_ptr_exponents()[1], 0.5);
		EXPECT_EQ(u.get_ptr_centers()[0], 1);
		u.set_atoms(w.get_atoms());
		EXPECT_EQ(u.get_atoms().size(), 2u);
		u.set_ncen(2);
		u.clear_atom_basis_set(0);
		EXPECT_EQ(u.get_atom_basis_set_size(0), 0);
	}
}
