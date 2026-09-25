#pragma once
//The ELI family of Kohout's electron localizability indicators, as pointwise field functions on
//density_source.h's concepts.  ELI-D for same-spin (alpha-alpha, beta-beta) and triplet-coupled pairs,
//and ELI-q for same-spin pairs.  Nothing here is wired into the cube writer or the basin code.
//
//References (the definitions below are quoted from these, not remembered):
//  [I]   M. Kohout, K. Pernal, F. R. Wagner, Yu. Grin, "Electron localizability indicator for correlated
//        wavefunctions. I: parallel-spin pairs", Theor. Chem. Acc. 112 (2004) 453.
//  [II]  M. Kohout, K. Pernal, F. R. Wagner, Yu. Grin, "... II: antiparallel-spin pairs", Theor. Chem.
//        Acc. 113 (2005) 287.                                    <- ELIA, see the note at the bottom
//  [III] M. Kohout, F. R. Wagner, Yu. Grin, "... III: singlet and triplet pairs", Theor. Chem. Acc. 119
//        (2008) 413-420, DOI 10.1007/s00214-007-0396-1.          <- triplet ELI-D, Eq. 52; ELI-q, Eq. 53
//  [q]   M. Kohout, Faraday Discuss. 135 (2007) 43.              <- the q-restricted partitioning
//  [DG]  M. Kohout, DGrid 5.2 (the reference implementation every convention here was checked against).
//
//Every convention below was pinned voxel-wise (ratio 1.0 to ~1e-8) against DGrid 5.2 on six systems:
//H2O, F2, epoxide (closed shell, N = 10/14/18) and three open-shell F atoms (N=9 S=1/2, N=8 S=1,
//N=7 S=3/2).  Where DGrid and the paper disagree, the disagreement is written out at the definition.
#include "wfn_class.h"
#include "density_source.h"

namespace eli_family
{
	//---------------------------------------------------------------------------------------------
	//The spin-resolved ingredients.  Everything in the family is built from these three per channel:
	//  rho_s   = sum_{i in s} n_i |phi_i|^2
	//  grad_s  = grad rho_s
	//  T_s     = sum_{i in s} n_i |grad phi_i|^2        (= 2 * the positive-definite kinetic energy
	//                                                     density; DGrid's own "tau" field is T_s/2,
	//                                                     its header reads 1/2*sum|grad-phi|2[r])
	//Index 0 = alpha, 1 = beta.
	//---------------------------------------------------------------------------------------------
	struct SpinFields
	{
		double rho[2]{ 0.0, 0.0 };
		d3 grad[2]{ {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
		double T[2]{ 0.0, 0.0 };
	};
	enum Spin { alpha = 0, beta = 1 };

	//One pass over the primitives, alpha and beta accumulated separately.  A restricted wavefunction
	//(no MO carries op 1) splits every occupation evenly between the two channels, which reproduces
	//rho_alpha = rho_beta = rho/2 exactly.
	void spin_fields(const WFN& wave, const d3& p, SpinFields& f);

	//Same-spin pair-density curvature, [I] Eq. 27 / [III] Eq. 29:
	//    g^s(r) = rho_s T_s - 1/4 |grad rho_s|^2
	//This is the O(s^2) coefficient of the spherically averaged same-spin pair density,
	//rho_2^{ss}(r,s) ~ (s^2/6) g^s(r).
	inline double g_same_spin(const double rho, const d3& grad, const double T)
	{
		return rho * T - 0.25 * (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
	}
	inline double g_same_spin(const SpinFields& f, const int s) { return g_same_spin(f.rho[s], f.grad[s], f.T[s]); }

	//Triplet-coupled pair curvature.  Spin-adapting the 2-matrix of a single determinant ([III] Sect. 2,
	//Eq. 51) splits the opposite-spin block into a singlet and an M_s = 0 triplet part, so
	//    rho_2^(t)(1,2) = rho_2^aa + rho_2^bb + rho_2^(t,0),
	//    rho_2^(t,0)    = 1/2[rho_a(1)rho_b(2) + rho_b(1)rho_a(2)] - gamma_a(1,2) gamma_b(2,1)
	//(the singlet partner carries +gamma_a gamma_b; the two add back to rho_a(1)rho_b(2)+rho_b(1)rho_a(2),
	//which is the whole opposite-spin pair density of a determinant).  Expanding each term about
	//r +- s/2 to O(s^2) and averaging over the sphere - the Hessian terms cancel - gives
	//    g^(t) = g^a + g^b + 1/2 (rho_a T_b + rho_b T_a) - 1/4 grad rho_a . grad rho_b
	//which collapses to 3 g^a for a closed shell - and DGrid's pair-volume ratio for the triplet is
	//exactly 3^(-3/8) on H2O, F2 and epoxide, to 8e-9.
	inline double g_triplet(const SpinFields& f)
	{
		const double dot = f.grad[0][0] * f.grad[1][0] + f.grad[0][1] * f.grad[1][1] + f.grad[0][2] * f.grad[1][2];
		return g_same_spin(f, 0) + g_same_spin(f, 1)
			+ 0.5 * (f.rho[0] * f.T[1] + f.rho[1] * f.T[0]) - 0.25 * dot;
	}

	//Pair-volume function, [I] Eq. 39 / [III] Eq. 52: the volume of the micro-cell that holds a fixed
	//pair number, V~_D = (12/g)^(3/8).  DGrid writes it as the field "pair-volume-function ...".
	inline double pair_volume(const double g) { return g > 0.0 ? std::pow(12.0 / g, constants::c_38) : 0.0; }

	//ELI-D, [I] Eq. 40:  Y_D^s(r) = rho_s(r) V~_D^s(r) = rho_s (12/g^s)^(3/8).
	//Identical to what WFN::computeELI already returns for a closed shell: there Rho = 2 rho_a and
	//tau = 2 T_a, so Rho*tau - 1/4|grad Rho|^2 = 4 g^a and 0.5*Rho*(48/4g)^(3/8) = rho_a (12/g^a)^(3/8).
	inline double eli_d(const double rho, const double g) { return rho > 0.0 && g > 0.0 ? rho * pair_volume(g) : 0.0; }

	//ELI-q for same-spin pairs, [q] / [DG] field "ELI-q alpha-alpha": the q-restricted partitioning
	//fixes the charge per micro-cell instead of the pair number, so the sampled quantity is the pair
	//number itself,
	//    Y_q^s(r) = g^s(r) / (12 rho_s(r)^(8/3)) = [Y_D^s(r)]^(-8/3).
	//Verified against DGrid to 1e-8 on all six systems.  For a single determinant ELI-q is therefore an
	//exact, strictly decreasing function of ELI-D - same topology, inverted, no new information.
	inline double eli_q(const double rho, const double g)
	{
		return rho > 0.0 && g > 0.0 ? g / (12.0 * std::pow(rho, 8.0 / 3.0)) : 0.0;
	}

	//The triplet density that ELI-D(triplet) samples, [III] Eq. 52: Y_D^(t) = rho^(t) (12/g^(t))^(3/8).
	//[III] Eq. 5 gives rho^(t) = (3/4)(N-2)/(N-1) rho for the unpolarized case; DGrid 5.2 does NOT use
	//that.  Its "rho triplet" field is, exactly on all six test systems,
	//    rho^(t) = rho * [1 - N_beta / (2(N-1))]
	//(13/18 for H2O, 19/26 for F2, 25/34 for epoxide, 3/4, 11/14, 5/6 for the three F atoms; Kohout's
	//Eq. 5 is off by 13/12 already for H2O).  We follow DGrid, the reference implementation.
	//Note this factor is a pure constant: rho^(t) is always proportional to rho, so for a closed shell
	//ELI-D(triplet) = 3^(-3/8) * (rho^(t)/rho_alpha) * ELI-D(alpha-alpha) - a constant times ELI-D.
	double triplet_density_factor(const WFN& wave);

	//---------------------------------------------------------------------------------------------
	//ELIA, the antiparallel-spin member ([II]), and the singlet ELI-q of [III] Eq. 53.
	//
	//Both sample the opposite-spin pair density, whose leading term at the coalescence point is the
	//singlet on-top density rho_2^(s)(r,r).  The derivation that matters here is one line: for a single
	//determinant the opposite-spin block of the 2-matrix has no exchange partner, because alpha and beta
	//spin-orbitals are orthogonal by spin, so
	//    rho_2^(ab)(1,2) = rho_a(1) rho_b(2)     exactly, at every separation,
	//and the singlet on-top density is rho_2^(s)(r,r) = rho_a(r) rho_b(r).  Substituted into Eq. 53,
	//    Y_q^(s)(r) = 4 rho_2^(s)(r,r) / rho(r)^2 = 4 rho_a rho_b / (rho_a + rho_b)^2 = 1 - zeta(r)^2,
	//    zeta = (rho_a - rho_b)/rho   (the spin polarisation).
	//So it is not merely uninformative, it is *exactly degenerate*: an algebraic function of the two spin
	//densities and nothing else, carrying no pair information at all.  For a closed-shell determinant
	//zeta == 0 and it is identically 1 everywhere - which is what the factor 4 in Eq. 53 is for.
	//
	//That expression is implemented below, because a caller asking for "the antiparallel member" deserves
	//a number and a reason rather than a missing symbol.  It is NOT a second opinion on the electron pair
	//structure; use it as a spin-polarisation map or not at all, and never run a basin ascent on it.
	//
	//DGrid 5.2 agrees that nothing better is available from this input: fed a single-determinant molden it
	//writes an identically zero field for both "ELIA both" and "ELIA alpha-beta" (0 nonzero of 49248
	//voxels on H2O, checked), with its own diagnostics reading "Singlet matrix not given" and "error in
	//read_dm_matrices ==> ALPHA-BETA 2-matrix element missing".
	//
	//What input would make it meaningful: the alpha-beta block of a *correlated* 2-matrix, i.e. the
	//genuine rho_2^(ab)(1,2) of a CI/CC/MCSCF wavefunction.  No format NoSpherA2 reads carries one -
	//wfn, wfx, molden, fchk and gbw all stop at orbitals plus occupations, and correlated *natural*
	//orbitals give a 1-RDM, from which the 2-RDM cannot be recovered.  DGrid takes it as separate .dm2
	//matrices.  Wiring a 2-matrix interface here would be dead code behind an input that never arrives.
	//---------------------------------------------------------------------------------------------
	inline double elia_singlet_eli_q(const SpinFields& f)
	{
		const double rho = f.rho[0] + f.rho[1];
		return rho > 0.0 ? 4.0 * f.rho[0] * f.rho[1] / (rho * rho) : 0.0;
	}
	const char* elia_status();

	//---------------------------------------------------------------------------------------------
	//Which members are worth computing for a given wavefunction, and which of them have a topology of
	//their own.  The decision is read off the wavefunction itself - whether a beta MO set exists
	//(get_MO_op_count(1)) and the actual occupations - never from a stated multiplicity.  A stated
	//multiplicity is used only as a cross-check that can raise a warning: a restricted file labelled as
	//a triplet is still restricted, and inventing a spin split from the label is how a spin heuristic
	//once doubled every delocalisation index while leaving the populations looking right.
	//
	//The rule, and the measurement behind each line:
	//  restricted (no beta MO set)   -> ELI-D(alpha-alpha) only.
	//      ELI-D(beta-beta) is the same field: rho_b = rho_a and T_b = T_a by construction.
	//      ELI-D(triplet) is a *constant* multiple of it, 3^(-3/8) * rho^(t)/rho_alpha - measured on
	//      H2O as 0.95670 against the predicted 2*(13/18)*3^(-3/8) = 0.956716, with DGrid's basin
	//      decomposition returning the same five basins, the same attractors and the same populations
	//      (1.7968 / 1.6792 / 1.6775 / 2.2420 / 2.2651).  Listing it separately would be five decimal
	//      places of nothing.  (Note DGrid itself has no beta orbital set for a restricted file: its
	//      tau_beta is 0, so its g_beta = -1/4|grad rho_b|^2 is negative and its elid_r_b_bb /
	//      eliq_r_b_bb fields are artefacts - eliq_r_b_bb comes out as one constant, 7.47e-13, over
	//      all 15625 voxels of the F(-) grid.  We return ELI-D_bb = ELI-D_aa, which is correct.)
	//  unrestricted                 -> ELI-D(alpha-alpha), ELI-D(beta-beta), ELI-D(triplet).
	//      All three are genuinely different fields with genuinely different basins.  Measured on the
	//      F atom doublet: DGrid finds four basins for both alpha-alpha and the triplet, but at
	//      different attractors and with different populations (1.2829/2.9960/2.1733/1.8770 against
	//      1.2829/1.2968/2.5635/3.1860), and the core maxima ratio is 1.03362 where the closed-shell
	//      relation would demand 0.99353.  The topology itself changes; the closed-shell shortcut does
	//      not survive spin polarisation.
	//  ELI-q                        -> computable always, basins never.
	//      Y_q = Y_D^(-8/3) exactly, a strictly decreasing function, so its maxima are ELI-D's minima.
	//      DGrid's own ELID_core ascent on ELI-q returns ~60 (H2O) / ~100 (F doublet) tail basins
	//      holding 0.0000-0.0005 e each.  basins_independent() refuses it.
	//  ELIA / singlet ELI-q         -> only for an unrestricted wavefunction, and only as 1 - zeta^2.
	//      Identically 1 for a restricted one, so there is nothing to report; see above.
	//---------------------------------------------------------------------------------------------
	enum class Member { eli_d_aa, eli_d_bb, eli_d_triplet, eli_q_aa, eli_q_bb, elia_singlet };
	//member_name() is the value to record in a method column, member_column() the header of the column
	//-eli_family's point dump writes it under.  There is deliberately no switch per member: one
	//-eli_family run evaluates the spin fields once and emits every member for that point, which is the
	//batched shape the family wants - g^aa, g^bb and g^(t) are three contractions of the same seven
	//invariants, so asking for them one at a time would re-do the primitive loop N times.
	const char* member_name(const Member m);
	const char* member_column(const Member m);
	bool basins_independent(const Member m);          //may a steepest ascent be run on this member?
	std::vector<Member> eli_variants_for(const WFN& wave, std::string* warning = nullptr);

	//CLI entry point for -eli_family: reports N_alpha/N_beta, the triplet density factor and the ELIA
	//status for one wavefunction and, given a file of "x y z" lines (bohr, '#' comments skipped),
	//writes one line of field values per point. That second mode is also the DGrid comparison tool.
	void report(const std::filesystem::path& wfn_path, const std::filesystem::path& points_file);

	//---------------------------------------------------------------------------------------------
	//Uniform electron gas limits, used as the analytic assertion in the tests.
	//  T_s^UEG  = (3/5)(6 pi^2)^(2/3) rho_s^(5/3)          (grad rho_s = 0)
	//  g^s,UEG  = (3/5)(6 pi^2)^(2/3) rho_s^(8/3)
	//  Y_D^UEG  = [12 / ((3/5)(6 pi^2)^(2/3))]^(3/8)       ~ 1.10848, density independent
	//  Y_q^UEG  = (Y_D^UEG)^(-8/3)
	//---------------------------------------------------------------------------------------------
	inline double ueg_g_coefficient() { return 0.6 * std::pow(6.0 * constants::PI2, 2.0 / 3.0); }
	inline double ueg_eli_d() { return std::pow(12.0 / ueg_g_coefficient(), constants::c_38); }
	inline double ueg_eli_q() { return ueg_g_coefficient() / 12.0; }
}

//---------------------------------------------------------------------------------------------
//The density_source.h-style point functions.  calculate_eli() there is ELI-D(alpha-alpha) already;
//these add the spin argument and the other family members.
//---------------------------------------------------------------------------------------------
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
//ELI-D for triplet-coupled pairs.  The density factor depends on N and N_beta only, so a caller
//sweeping a grid should hoist eli_family::triplet_density_factor(w) out of the loop.
inline double calculate_eli_d_triplet(const WFN& w, const d3& p, const double density_factor)
{
	eli_family::SpinFields f;
	eli_family::spin_fields(w, p, f);
	return eli_family::eli_d(density_factor * (f.rho[0] + f.rho[1]), eli_family::g_triplet(f));
}
inline double calculate_eli_d_triplet(const WFN& w, const d3& p) { return calculate_eli_d_triplet(w, p, eli_family::triplet_density_factor(w)); }

//A source that is only a density (a fitted Gaussian density, an atom model) has no orbitals and no
//spin resolution: it falls back to the PC07 orbital-free ELI-D of calculate_eli(), i.e. the
//spin-restricted value, and to its exact ELI-q partner.  No triplet member - that needs N.
template<PointDensity S> double calculate_eli_d(const S& s, const d3& p, const int) { return calculate_eli(s, p); }
template<PointDensity S> double calculate_eli_q(const S& s, const d3& p, const int)
{
	const double y = calculate_eli(s, p);
	return y > 0.0 ? std::pow(y, -8.0 / 3.0) : 0.0;
}
