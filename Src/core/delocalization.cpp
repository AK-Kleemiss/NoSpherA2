#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "b2c.h"
#include "citations.h"
#include "eli_family.h"
#include "ci_rdm2.h"

delocalization_result delocalization_indices(const WFN &wavy, const basin_overlaps &ovl, const CiRdm2 *rdm)
{
	delocalization_result r;
	const int nb = static_cast<int>(ovl.S.size());
	const int n = ovl.nmo;
	r.lambda.assign(nb, 0.0);
	r.population.assign(nb, 0.0);
	r.outside_half.assign(nb, 0.0);
	if (nb == 0 || n == 0) return r;
	//Spin-orbital occupations as the ELI family splits them: an ROHF singly occupied orbital is all alpha
	const eli_family::SpinSplit how = eli_family::spin_split(wavy);
	vec na(n), nbeta(n);
	ivec spin(n);
	for (int i = 0; i < n; i++) {
		double ns[2];
		eli_family::mo_spin_occupations(wavy.get_MO_occ(ovl.mo_index[i]), wavy.get_MO_op(ovl.mo_index[i]), how, ns);
		na[i] = ns[0], nbeta[i] = ns[1];
		//Spatial orbitals stand for both spins and all of them exchange with one another; only
		//a genuinely spin-resolved set has an alpha and a beta block that must not mix
		spin[i] = how == eli_family::SpinSplit::unrestricted ? wavy.get_MO_op(ovl.mo_index[i]) : 0;
	}
	for (int b = 0; b < nb; b++)
		for (int i = 0; i < n; i++)
			r.population[b] += (na[i] + nbeta[i]) * ovl.at(b, i, i);
	//The basins and outside region together resolve the orbital metric.
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= i; j++) {
			if (spin[i] != spin[j]) continue;
			double s = 0.0;
			for (int b = 0; b < nb; b++) s += ovl.at(b, i, j);
			if (!ovl.outside.empty()) s += ovl.outside[basin_overlaps::packed(i, j)];
			r.identity_error = std::max(r.identity_error, std::abs(s - (i == j ? 1.0 : 0.0)));
		}
	if (rdm) {
		//Mean-field part over the natural orbitals, both spins in one occupation
		auto mean_field = [&](const vec &a, const vec &b) {
			double s = 0.0;
			for (int i = 0; i < n; i++)
				for (int j = 0; j <= i; j++) {
					const size_t k = basin_overlaps::packed(i, j);
					const double t = (na[i] + nbeta[i]) * (na[j] + nbeta[j]) * a[k] * b[k];
					s += i == j ? t : 2.0 * t;
				}
			return s;
		};
		//The cumulant L as an m^2 x m^2 matrix, so C(A,B) = s_A^T L s_B with s_A[t*m+u] = S^A_tu
		const int m = rdm->n;
		const size_t m2 = size_t(m) * m;
		vec L(rdm->G);
		for (int t = 0; t < m; t++)
			for (int v = 0; v < m; v++) {
				const double nn = wavy.get_MO_occ(rdm->first + t) * wavy.get_MO_occ(rdm->first + v);
				L[rdm->at(t, t, v, v)] -= nn;
				L[rdm->at(t, v, v, t)] += 0.5 * nn;
			}
		//An active orbital below the AOM's 1e-8 occupation cut has no overlaps, and its cumulant is as small
		ivec slot(m, -1);
		for (int i = 0; i < n; i++) {
			const int t = ovl.mo_index[i] - rdm->first;
			if (t >= 0 && t < m) slot[t] = i;
		}
		const bool out = !ovl.outside.empty();
		vec2 s(nb + 1, vec(m2, 0.0)), Ls(nb + 1, vec(m2, 0.0));
		for (int b = 0; b < nb + (out ? 1 : 0); b++) {
			const vec &S = b < nb ? ovl.S[b] : ovl.outside;
			for (int t = 0; t < m; t++)
				for (int u = 0; u < m; u++)
					if (slot[t] >= 0 && slot[u] >= 0) s[b][size_t(t) * m + u] = S[basin_overlaps::packed(slot[t], slot[u])];
			for (size_t i = 0; i < m2; i++)
				Ls[b][i] = std::inner_product(s[b].begin(), s[b].end(), L.begin() + i * m2, 0.0);
		}
		auto C = [&](const int a, const int b) { return std::inner_product(s[a].begin(), s[a].end(), Ls[b].begin(), 0.0); };
		r.exact = true;
		for (int b = 0; b < nb; b++) {
			r.lambda[b] = 0.5 * mean_field(ovl.S[b], ovl.S[b]) - C(b, b);
			if (out) r.outside_half[b] = 0.5 * mean_field(ovl.S[b], ovl.outside) - C(b, nb);
		}
		for (int a = 0; a < nb; a++)
			for (int b = a + 1; b < nb; b++) {
				r.pairs.push_back({ a, b });
				r.di.push_back(mean_field(ovl.S[a], ovl.S[b]) - 2.0 * C(a, b));
			}
		return r;
	}
	//The pair sum is symmetric in i and j, so the triangle is taken once and doubled off the
	//diagonal.  Mueller (Buijse-Baerends) exchange-correlation hole, sqrt(n_i n_j) per spin: n_i n_j for
	//0/1 occupations, and for natural orbitals it keeps lambda(A) + sum_B delta(A,B)/2 = N(A), which n_i n_j
	//breaks.  Relaxed densities' small negative n are clipped inside the root only.
	auto pair_sum = [&](const vec &a, const vec &b) {
		double s = 0.0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j <= i; j++) {
				if (spin[i] != spin[j]) continue;
				const size_t k = basin_overlaps::packed(i, j);
				const double t = (std::sqrt(std::max(na[i], 0.0) * std::max(na[j], 0.0)) +
					std::sqrt(std::max(nbeta[i], 0.0) * std::max(nbeta[j], 0.0))) * a[k] * b[k];
				s += i == j ? t : 2.0 * t;
			}
		return s;
	};
	for (int b = 0; b < nb; b++) {
		r.lambda[b] = pair_sum(ovl.S[b], ovl.S[b]);
		if (!ovl.outside.empty()) r.outside_half[b] = pair_sum(ovl.S[b], ovl.outside);
	}
	for (int a = 0; a < nb; a++)
		for (int b = a + 1; b < nb; b++) {
			r.pairs.push_back({ a, b });
			r.di.push_back(2.0 * pair_sum(ovl.S[a], ovl.S[b]));
		}
	return r;
}

void report_delocalization(const WFN &wavy, const basin_overlaps &ovl, const svec &labels, std::ostream &log, const double threshold)
{
	if (ovl.S.empty() || ovl.nmo == 0) return;
	//A correlated 1-RDM: the exact pair density of a CI printed next to the wavefunction, Mueller's otherwise
	const bool correlated = !eli_family::noon_summary(wavy).idempotent;
	CiRdm2 rdm;
	std::ostringstream why;
	const bool exact = correlated && load_ci_rdm2(wavy, rdm, why);
	const delocalization_result r = delocalization_indices(wavy, ovl, exact ? &rdm : nullptr);
	const int nb = static_cast<int>(r.lambda.size());
	auto name = [&](const int b) { return b < static_cast<int>(labels.size()) ? labels[b] : std::to_string(b + 1); };
	log << "\nDelocalization indices (" << ovl.nmo << " occupied orbitals):\n";
	citations::cite(citations::Method::LIDI, log);
	if (correlated)
		log << why.str() << (exact ? "  pair density: exact, from the CI 2-RDM\n"
			: "  pair density: Mueller, a 1-RDM approximation; exact indices need the CI 2-RDM\n");
	log << "  sum over basins and outside of S^A - identity: " << std::scientific << std::setprecision(2) << r.identity_error
		<< std::fixed << "   (the quadrature's own error; AIMAll's integrations reach ~1e-3)\n";
	//delta(A,B) summed over B is the count an atom shares with everything else; with lambda(A)
	//it has to give the population back, and the residual says which basin the grid missed
	log << "\n  basin  label                 N(A)     lambda(A)   sum_B delta(A,B)/2   residual\n";
	log << "  B includes the region beyond rho = " << std::scientific << std::setprecision(0) << basin_density_cutoff << std::fixed << std::setprecision(2) << " e/bohr^3.\n";
	for (int b = 0; b < nb; b++) {
		double half = 0.0;
		for (size_t p = 0; p < r.pairs.size(); p++)
			if (r.pairs[p][0] == b || r.pairs[p][1] == b) half += 0.5 * r.di[p];
		half += r.outside_half[b];
		log << std::setw(7) << b + 1 << "  " << std::left << std::setw(18) << name(b) << std::right << std::fixed << std::setprecision(4)
			<< std::setw(11) << r.population[b] << std::setw(12) << r.lambda[b] << std::setw(18) << half
			<< std::setw(12) << r.lambda[b] + half - r.population[b] << "\n";
	}
	ivec order(r.di.size());
	std::iota(order.begin(), order.end(), 0);
	std::sort(order.begin(), order.end(), [&](const int a, const int b) { return r.di[a] > r.di[b]; });
	log << "\n  basin pair                                delta(A,B)\n";
	int shown = 0;
	for (const int p : order) {
		if (r.di[p] < threshold) break;
		log << "  " << std::left << std::setw(18) << name(r.pairs[p][0]) << std::setw(18) << name(r.pairs[p][1])
			<< std::right << std::fixed << std::setprecision(4) << std::setw(12) << r.di[p] << "\n";
		shown++;
	}
	if (!shown) log << "  none above " << threshold << "\n";
	else log << "  " << static_cast<int>(r.di.size()) - shown << " further pairs below " << std::setprecision(2) << threshold << "\n";
	for (int a = 0; a < wavy.get_ncen(); a++)
		if (wavy.get_atom_ECP_electrons(a) > 0) {
			log << "  An ECP took core electrons out of the orbitals: N(A) here is the valence count\n"
				<< "  and is short of the basin population above by the core the Thakkar fill added.\n";
			break;
		}
}
