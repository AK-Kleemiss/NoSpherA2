#pragma once
//The exact 2-RDM of a determinant CI that ORCA printed (%casscf PrintWF det, TPrintWF 1e-14): no gbw or
//.densities file carries a CI vector or a 2-RDM, so the "[2ud0...] coef" block of <stem>.out/.log next to the
//gbw is the only source.  Active orbitals only: the internal ones are doubly occupied in every determinant,
//so every RDM element that touches them is the mean-field one of the 1-RDM.
#include "wfn_class.h"

struct CiRdm2
{
	int first = 0;   //WFN MO index of active orbital 0 (ORCA's internal count)
	int n = 0;       //active orbitals
	vec g1[2];       //gamma^s_tu = <a+_ts a_us>, s = alpha, beta; [t*n+u]
	vec G;           //spin-summed Gamma_tuvw = sum_st <a+_ts a+_vt a_wt a_us>, [((t*n+u)*n+v)*n+w]; sum_tv G_ttvv = N(N-1)
	vec Gab;         //alpha-beta block <a+_ta a+_vb a_wb a_ua>; rho_2^ab(1,2) = sum Gab_tuvw phi_t phi_u(1) phi_v phi_w(2)
	std::filesystem::path source;
	size_t at(const int t, const int u, const int v, const int w) const { return ((size_t(t) * n + u) * n + v) * n + w; }
};

//Reads <stem>.out or <stem>.log next to wave.get_path().  false, with the reason on log, when there is no
//single-root det block, the WFN is unrestricted, or the CI 1-RDM does not reproduce the WFN occupations
//(1e-4: catches a truncated print, a phase error and a CI in other orbitals than the gbw holds)
bool load_ci_rdm2(const WFN& wave, CiRdm2& r, std::ostream& log);

//On-top alpha-beta pair density rho_2^ab(r,r) at p, the internal orbitals doubly occupied in every determinant;
//rho[0], rho[1] get rho_alpha, rho_beta of the same 2-RDM, so a single determinant gives exactly rho_a rho_b
double ci_on_top_ab(const WFN& wave, const CiRdm2& r, const d3& p, double rho[2]);

//Exported for the tests: dets as (ORCA string, coefficient), n orbitals, into r.g1/G/Gab in that orbital basis
void ci_rdm2_from_dets(const std::vector<std::pair<std::string, double>>& dets, CiRdm2& r);
