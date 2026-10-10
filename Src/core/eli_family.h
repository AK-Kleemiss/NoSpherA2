#pragma once
//Kohout's ELI family as point functions on density_source.h: ELI-D for same-spin and triplet-coupled
//pairs, ELI-q for same-spin pairs.  Conventions follow DGrid 5.2 [DG] where it differs from the papers.
//[I] Kohout, Pernal, Wagner, Grin, Theor. Chem. Acc. 112 (2004) 453 (parallel spin)
//[II] same authors, Theor. Chem. Acc. 113 (2005) 287 (antiparallel spin, ELIA)
//[III] Kohout, Wagner, Grin, Theor. Chem. Acc. 119 (2008) 413 (singlet and triplet pairs)
//[q] Kohout, Faraday Discuss. 135 (2007) 43 (q-restricted partitioning)
#include "wfn_class.h"
#include "density_source.h"

namespace eli_family
{
	//Per spin s (0 alpha, 1 beta): rho_s = sum_{i in s} n_i |phi_i|^2, grad rho_s, and
	//T_s = sum_{i in s} n_i |grad phi_i|^2 = 2 tau_s (DGrid's "tau" field is T_s/2)
	struct SpinFields
	{
		double rho[2]{ 0.0, 0.0 };
		d3 grad[2]{ {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
		double T[2]{ 0.0, 0.0 };
	};
	enum Spin { alpha = 0, beta = 1 };

	//Occupation split over the spins, read off the orbitals.  unrestricted: each MO in its own channel
	//(the spin flags count only while no MO holds more than one electron).
	//restricted_open: one set, occupations 0/1/2 with at least one 1 and multiplicity (if stated) = singles + 1;
	//a single is alpha, a double one of each.  halves: any other single set, occ/2 each.
	//Serial loops only: computeELISpinGrad calls it once per point of a climb
	enum class SpinSplit { unrestricted, restricted_open, halves };
	SpinSplit spin_split(const WFN& wave);
	inline void mo_spin_occupations(const double occ, const int op, const SpinSplit how, double n[2])
	{
		if (how == SpinSplit::unrestricted) { n[0] = op == 1 ? 0.0 : occ; n[1] = op == 1 ? occ : 0.0; }
		else if (how == SpinSplit::restricted_open) { n[0] = std::min(occ, 1.0); n[1] = occ - n[0]; }
		else n[0] = n[1] = 0.5 * occ;
	}
	//Broken-symmetry singlet test: occupied alpha and beta MOs compared in order (1e-4 of the largest
	//coefficient, sign free).  An occupied-space rotation leaves rho_beta = rho_alpha, so a difference is
	//confirmed on the spin density near each nucleus
	bool alpha_beta_orbitals_differ(const WFN& wave);

	void spin_fields(const WFN& wave, const d3& p, SpinFields& f);

	//Same-spin pair-density curvature, [I] Eq. 27 / [III] Eq. 29: g^s = rho_s T_s - 1/4 |grad rho_s|^2,
	//the O(s^2) coefficient of the spherically averaged rho_2^ss(r,s) ~ (s^2/6) g^s(r)
	inline double g_same_spin(const double rho, const d3& grad, const double T)
	{
		return rho * T - 0.25 * (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
	}
	inline double g_same_spin(const SpinFields& f, const int s) { return g_same_spin(f.rho[s], f.grad[s], f.T[s]); }

	//Triplet pair curvature, [III] Eq. 51: rho_2^(t) = rho_2^aa + rho_2^bb + rho_2^(t,0) with
	//rho_2^(t,0) = 1/2[rho_a(1)rho_b(2) + rho_b(1)rho_a(2)] - gamma_a(1,2) gamma_b(2,1).  Expanded to O(s^2)
	//and sphere-averaged (Hessian terms cancel):
	//g^(t) = g^a + g^b + 1/2 (rho_a T_b + rho_b T_a) - 1/4 grad rho_a . grad rho_b, = 3 g^a for a closed shell
	inline double g_triplet(const SpinFields& f)
	{
		const double dot = f.grad[0][0] * f.grad[1][0] + f.grad[0][1] * f.grad[1][1] + f.grad[0][2] * f.grad[1][2];
		return g_same_spin(f, 0) + g_same_spin(f, 1)
			+ 0.5 * (f.rho[0] * f.T[1] + f.rho[1] * f.T[0]) - 0.25 * dot;
	}

	//Micro-cell volume holding a fixed pair number, [I] Eq. 39 / [III] Eq. 52: V~_D = (12/g)^(3/8)
	inline double pair_volume(const double g) { return g > 0.0 ? std::pow(12.0 / g, constants::c_38) : 0.0; }

	//ELI-D, [I] Eq. 40: Y_D^s = rho_s (12/g^s)^(3/8); equals WFN::computeELI for a closed shell
	inline double eli_d(const double rho, const double g) { return rho > 0.0 && g > 0.0 ? rho * pair_volume(g) : 0.0; }

	//ELI-q [q], charge per micro-cell fixed instead of pair number: Y_q^s = g^s / (12 rho_s^(8/3)) = (Y_D^s)^(-8/3),
	//a strictly decreasing function of ELI-D, so the same topology inverted
	inline double eli_q(const double rho, const double g)
	{
		return rho > 0.0 && g > 0.0 ? g / (12.0 * std::pow(rho, 8.0 / 3.0)) : 0.0;
	}

	//rho^(t)/rho in Y_D^(t) = rho^(t) (12/g^(t))^(3/8), [III] Eq. 52.  DGrid's 1 - N_beta / (2(N-1)), not
	//[III] Eq. 5's (3/4)(N-2)/(N-1).  A constant, so closed-shell ELI-D(triplet) is a multiple of ELI-D(aa)
	double triplet_density_factor(const WFN& wave);

	//ELIA [II] / singlet ELI-q [III] Eq. 53.  For a single determinant rho_2^ab(1,2) = rho_a(1) rho_b(2)
	//exactly, so Y_q^(s) = 4 rho_a rho_b / rho^2 = 1 - zeta^2, zeta = (rho_a - rho_b)/rho: a spin-polarisation
	//map with no pair information, never a basin target.  A meaningful ELIA needs a correlated 2-RDM,
	//which no format read here carries
	inline double elia_singlet_eli_q(const SpinFields& f)
	{
		const double rho = f.rho[0] + f.rho[1];
		return rho > 0.0 ? 4.0 * f.rho[0] * f.rho[1] / (rho * rho) : 0.0;
	}
	const char* elia_status();

	//Members worth computing, decided from the orbitals (beta MO set, occupations), never from a stated
	//multiplicity, which only raises a warning.  Restricted: ELI-D(aa) only, since bb is identical and the
	//triplet a constant multiple.  Unrestricted: aa, bb and triplet, which have different basins.  ELI-q:
	//always computable, never a basin target (its maxima are ELI-D's minima).  ELIA: unrestricted only
	enum class Member { eli_d_aa, eli_d_bb, eli_d_triplet, eli_q_aa, eli_q_bb, elia_singlet };
	//Every member comes from one spin_fields pass, so there is no per-member evaluator
	const char* member_name(const Member m);
	const char* member_column(const Member m);
	bool basins_independent(const Member m);          //may a steepest ascent be run on this member?
	std::vector<Member> eli_variants_for(const WFN& wave, std::string* warning = nullptr);

	//-eli_family: N_alpha/N_beta, triplet factor and ELIA status; with a file of "x y z" lines (bohr)
	//one line of field values per point
	void report(const std::filesystem::path& wfn_path, const std::filesystem::path& points_file);

	//Uniform electron gas: g^s = (3/5)(6 pi^2)^(2/3) rho_s^(8/3), so Y_D = [12 / ((3/5)(6 pi^2)^(2/3))]^(3/8)
	//~ 1.10848 independent of density, and Y_q = Y_D^(-8/3)
	inline double ueg_g_coefficient() { return 0.6 * std::pow(6.0 * constants::PI2, 2.0 / 3.0); }
	inline double ueg_eli_d() { return std::pow(12.0 / ueg_g_coefficient(), constants::c_38); }
	inline double ueg_eli_q() { return ueg_g_coefficient() / 12.0; }
}

inline double calculate_eli_d(const WFN& w, const d3& p, const int spin)
{
	eli_family::SpinFields f;
	eli_family::spin_fields(w, p, f);
	return eli_family::eli_d(f.rho[spin], eli_family::g_same_spin(f, spin));
}
inline double calculate_eli_q(const WFN& w, const d3& p, const int spin)
{
	eli_family::SpinFields f;
	eli_family::spin_fields(w, p, f);
	return eli_family::eli_q(f.rho[spin], eli_family::g_same_spin(f, spin));
}
//density_factor depends on N and N_beta only: hoist triplet_density_factor(w) out of a grid loop
inline double calculate_eli_d_triplet(const WFN& w, const d3& p, const double density_factor)
{
	eli_family::SpinFields f;
	eli_family::spin_fields(w, p, f);
	return eli_family::eli_d(density_factor * (f.rho[0] + f.rho[1]), eli_family::g_triplet(f));
}
inline double calculate_eli_d_triplet(const WFN& w, const d3& p) { return calculate_eli_d_triplet(w, p, eli_family::triplet_density_factor(w)); }

//A density-only source has no spin resolution: PC07 orbital-free ELI-D of calculate_eli() and its
//ELI-q partner; no triplet, which needs N
template<PointDensity S> double calculate_eli_d(const S& s, const d3& p, const int) { return calculate_eli(s, p); }
template<PointDensity S> double calculate_eli_q(const S& s, const d3& p, const int)
{
	const double y = calculate_eli(s, p);
	return y > 0.0 ? std::pow(y, -8.0 / 3.0) : 0.0;
}
