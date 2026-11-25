#include "pch.h"
#include "tuning.h"
#include "wfn_class.h"
#include "convenience.h"
#include "bondwise_analysis.h"
#include "properties.h"
#include "b2c.h"
#include "eli_family.h"
#include "topology.h"
#include "gaussian_atom.h"
#include "spherical_density.h"
#include "citations.h"
#include "section_log.h"
#include "ibo.h"
#include "nbo.h"
#include "basis_set.h"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>

//The fitted density in place of the orbitals for rho and the QTAIM basins when the options
//ask for it: -ri_fit <basis> before the analysis flag on a wavefunction, or -SALTED
//<model-dir> on an xyz. rho and its gradient then come from one loop over the auxiliary
//functions, which is what makes 0.05 A grids affordable. nullptr for the orbital route.
//ELI-D keeps the orbitals: the density-only estimate (PC07 alpha, the -eli cubes) is a
//constant wherever alpha is switched off, i.e. over the bonds and lone pairs, and has no
//basins to find there.
static std::unique_ptr<Gaussian_Molecule> fitted_source(const WFN &wavy, options &opt, density_field &field, std::ostream &log)
{
	if (!Gaussian_Molecule::requested(wavy, opt)) {
		err_checkf(wavy.get_nmo() > 0, "No orbitals in " + wavy.get_path().string() + "; give -SALTED <model-dir> before the analysis flag to analyse a predicted density", log);
		return nullptr;
	}
	auto fit = std::make_unique<Gaussian_Molecule>(wavy, opt);
	log << "Density source: " << (wavy.get_nmo() == 0 ? "SALTED prediction from " + opt.salted_model_dir.string() : "RI fit of " + std::to_string(wavy.get_nmo()) + " orbitals")
		<< ", " << fit->n_aux() << " auxiliary functions.\n"
		<< "  rho and its gradient come from one pass over the auxiliary functions, N_aux per point instead of\n"
		<< "  N_basis x N_MO, so grids of 0.05 A and finer are affordable. The QTAIM basins and populations below\n"
		<< "  are those of the fitted density: it matches the orbital density to ~1e-3 e/bohr^3 at bond critical\n"
		<< "  points, and its tail (rho < 0.01) is not reliable. ELI-D "
		<< (wavy.get_nmo() == 0 ? "needs orbitals and is skipped." : "is taken from the orbitals.") << std::endl;
	const Gaussian_Molecule *f = fit.get();
	field.rho = [f](const d3 &p) { return f->rho(p); };
	field.grad = [f](const d3 &p, d3 &g) { double lap; f->values(p, g, lap); };
	return fit;
}

//IBO and NBO labels number atoms from 1; the table uses the 0-based numbers of the QTAIM and ELI-D labels
static std::string zero_based(const std::string &l) {
	std::string s;
	for (size_t i = 0; i < l.size();) {
		if (!std::isdigit(static_cast<unsigned char>(l[i])) || i == 0 || !std::isalpha(static_cast<unsigned char>(l[i - 1]))) { s += l[i++]; continue; }
		size_t j = i;
		while (j < l.size() && std::isdigit(static_cast<unsigned char>(l[j]))) j++;
		s += std::to_string(std::stoi(l.substr(i, j - i)) - 1);
		i = j;
	}
	return s;
}

//repeated labels get (1), (2), ... in their order: after the tag of an orbital (LP(1) O0), after an ELI-D LP or bond (O0 LP(1));
//cores stay one label, the table sums them
static void number_repeats(svec &lab, const bool orbital) {
	std::map<std::string, int> n, seen;
	for (const std::string &l : lab) n[l]++;
	const auto ends = [](const std::string &l, const std::string &e) { return l.size() >= e.size() && l.compare(l.size() - e.size(), e.size(), e) == 0; };
	for (std::string &l : lab) {
		if (n[l] < 2 || (orbital ? l.rfind("CR ", 0) == 0 : !ends(l, " LP") && !ends(l, " bond"))) continue;
		const std::string k = "(" + std::to_string(++seen[l]) + ")";
		const size_t sp = l.find(' ');
		l = orbital && sp != std::string::npos ? l.substr(0, sp) + k + l.substr(sp) : l + k;
	}
}

//Kohout's spin-resolved ELI-D (eli_family.h). For an open shell the spin-summed ELI-D is neither spin's pair
//function; the alpha-alpha, beta-beta and triplet members are, each with its own maxima and basins.
//Members come from the orbitals (eli_variants_for), never from the stated multiplicity.
static void spin_eli_analysis(const WFN &l_w, const options &opt, const std::vector<atom> &atoms, const double shell_dist, const double shell_tol)
{
	using eli_family::Member;
	std::string warn;
	const std::vector<Member> members = eli_family::eli_variants_for(l_w, &warn);
	auto has = [&](const Member m) { return std::find(members.begin(), members.end(), m) != members.end(); };
	if (!opt.spin_eli) { std::cout << "\nSpin-resolved ELI-D: off (-no_spin_eli)." << std::endl; return; }
	const eli_family::SpinSplit how = eli_family::spin_split(l_w);
	const bool unrestricted = how == eli_family::SpinSplit::unrestricted;
	double N[2]{ 0.0, 0.0 };
	int norb[2]{ 0, 0 };
	bool open_restricted = false;
	for (int mo = 0; mo < l_w.get_nmo(); mo++) {
		const double occ = l_w.get_MO_occ(mo);
		if (occ == 0.0) continue;
		double n[2];
		eli_family::mo_spin_occupations(occ, l_w.get_MO_op(mo), how, n);
		for (int s = 0; s < 2; s++) if (n[s] != 0.0) { N[s] += n[s]; norb[s]++; }
		if (how == eli_family::SpinSplit::halves && std::abs(occ - 2.0) > 1e-6) open_restricted = true;
	}
	if (!has(Member::eli_d_bb)) {
		std::cout << "\nSpin-resolved ELI-D: skipped, "
			<< (open_restricted ? "fractional occupations in one MO set without spin labels (natural orbitals) - alpha and beta densities are not resolved"
				: !unrestricted ? "restricted wavefunction - ELI-D(alpha-alpha) = ELI-D(beta-beta) is the ELI-D above and the triplet member a constant multiple of it"
				: N[1] < 1e-8 ? "the beta orbital set holds no electrons"
				: "N_alpha = N_beta and the beta orbitals are the alpha ones - not spin-polarised")
			<< "." << std::endl;
		return;
	}
	if (opt.basin_cube) { std::cout << "\nSpin-resolved ELI-D: skipped, it runs on the analytic field only and -basin_cube was given." << std::endl; return; }
	citations::cite(citations::Method::ELIFamily, std::cout);
	const double tf = eli_family::triplet_density_factor(l_w);
	std::cout << "\nSpin-resolved ELI-D (Kohout): N_alpha = " << std::fixed << std::setprecision(4) << N[0] << ", N_beta = " << N[1]
		<< (how == eli_family::SpinSplit::restricted_open ? " (restricted open shell: singly occupied MOs alpha, doubly occupied one of each)"
			: std::abs(N[0] - N[1]) < 1e-8 ? " (broken symmetry: N_alpha = N_beta, alpha and beta orbitals differ)" : "")
		<< ". The ELI-D above is the spin-summed field, not a pair function of either spin." << std::endl;
	struct member_basins { std::vector<d4> maxima; svec labels; vec pop, vol; vec2 spin; double outside = 0.0; bool done = false; };
	const char *names[3] = { "alpha-alpha", "beta-beta", "triplet" };
	member_basins res[3];
	const int nf = has(Member::eli_d_triplet) ? 3 : 2;
	for (int f = 0; f < nf; f++) {
		//A channel carried by one orbital has g_s = 0 everywhere: no field, only round-off basins
		if (f < 2 && norb[f] < 2) {
			std::cout << "\nELI-D " << names[f] << ": skipped, the " << (f ? "beta" : "alpha") << " electrons sit in a single orbital, so g_s = 0 and the field is undefined." << std::endl;
			continue;
		}
		basin_stage_timer clock;
		const eli_spin_field eval = [&l_w, f, tf](const d3 &p, double &y, d3 &g, double *aux) { l_w.computeELISpinGrad(p, f, tf, y, g, aux); };
		member_basins &r = res[f];
		const std::vector<d4> all = analytic_eli_maxima(l_w, opt.debug, &eval);
		r.maxima = all;
		cubei none;
		ivec core_map, shell_map;
		std::vector<shell_split> split;
		unify_core_basins(none, r.maxima, atoms, &core_map);
		unify_shell_basins(none, r.maxima, &shell_map, shell_dist, shell_tol, &atoms, &l_w, &eval, &split);
		if (!shell_map.empty())
			for (size_t b = 1; b < core_map.size(); b++) core_map[b] = shell_map[core_map[b]];
		ivec edge_map;
		if (unify_boundary_basins(r.maxima, l_w, eval, &edge_map) > 0)
			for (size_t b = 1; b < core_map.size(); b++) core_map[b] = edge_map[core_map[b]];
		r.labels = assign_labels_to_basins(r.maxima, atoms, opt.debug, 1);
		number_repeats(r.labels, false);
		clock.lap(std::string("ELI-D ") + names[f] + " maxima");
		r.pop = integrate_basins_on_atomic_grids(nullptr, nullptr, all, l_w, opt.accuracy, true, r.vol, r.outside, nullptr, nullptr, opt.basin_grid, nullptr, nullptr, &core_map, &eval, &r.spin, &split);
		r.done = true;
		clock.lap(std::string("ELI-D ") + names[f] + " basins");
		double tot = 0.0, ts[2]{ 0.0, 0.0 };
		for (size_t b = 0; b < r.pop.size(); b++) {
			tot += r.pop[b];
			ts[0] += r.spin[b][0];
			ts[1] += r.spin[b][1];
		}
		//ponytail: rim seeds find empty surface attractors of the bounded field; they count in the totals but are not listed
		int empty = 0;
		for (size_t b = r.pop.size(); b-- > 0;)
			if (r.pop[b] < 5e-5) {
				r.pop.erase(r.pop.begin() + b); r.vol.erase(r.vol.begin() + b); r.spin.erase(r.spin.begin() + b);
				r.maxima.erase(r.maxima.begin() + b); r.labels.erase(r.labels.begin() + b);
				empty++;
			}
		//all.size() also counts the maxima unified into cores/shells and the empty rim attractors: say which
		const size_t merged = all.size() - r.maxima.size() - empty;
		std::cout << "\nELI-D " << names[f] << " Analysis (atomic quadrature grids), " << all.size() << " maxima";
		if (merged) std::cout << ", " << merged << " unified into core/shell basins";
		if (empty) std::cout << ", " << empty << " empty surface attractors (below 5e-5 e, not listed)";
		std::cout << ", " << r.maxima.size() << " basins";
		std::cout << ":\n"
			<< "  basin  label               electrons    N_alpha     N_beta       spin  <ELI-q>      volume         maximum        x          y          z\n";
		for (size_t b = 0; b < r.pop.size(); b++) {
			std::cout << std::setw(7) << b + 1 << "  " << std::left << std::setw(18) << r.labels[b] << std::right << std::fixed << std::setprecision(4)
				<< std::setw(11) << r.pop[b] << std::setw(11) << r.spin[b][0] << std::setw(11) << r.spin[b][1] << std::setw(11) << r.spin[b][0] - r.spin[b][1];
			//<ELI-q_s> = integral of rho_s Y_q over the basin / N_s; the triplet has no ELI-q partner
			//below 1e-3 e of the channel the average is a tail ratio of two vanishing integrals
			if (f < 2 && r.spin[b][f] > 1e-3) std::cout << std::setw(9) << r.spin[b][2] / r.spin[b][f];
			else std::cout << std::setw(9) << "-";
			std::cout << std::setw(12) << r.vol[b] << std::setw(16) << r.maxima[b][3]
				<< std::setprecision(3) << std::setw(11) << r.maxima[b][0] << std::setw(11) << r.maxima[b][1] << std::setw(11) << r.maxima[b][2] << "\n";
		}
		std::cout << std::setprecision(4) << "  total in basins: " << tot << "   N_alpha " << ts[0] << " of " << N[0] << ", N_beta " << ts[1] << " of " << N[1]
			<< "   outside every basin: " << r.outside << std::endl;
	}
	//aa vs bb: only basins with the same label pair. A label held once pairs outright (a unified core or shell
	//basin's maximum is an arbitrary point on its shell); a repeated label pairs by mutual nearest maxima within
	//a bohr. An alpha basin left without a partner is where the unpaired alpha electrons sit
	if (res[0].done && res[1].done) {
		const svec &LA = res[0].labels, &LB = res[1].labels;
		auto nearest = [](const d4 &m, const std::vector<d4> &in, const svec &lin, const std::string &l, double &d) {
			int best = -1;
			d = 1e30;
			for (size_t k = 0; k < in.size(); k++) {
				if (lin[k] != l) continue;
				const double q = std::sqrt(std::pow(m[0] - in[k][0], 2) + std::pow(m[1] - in[k][1], 2) + std::pow(m[2] - in[k][2], 2));
				if (q < d) { d = q; best = static_cast<int>(k); }
			}
			return best;
		};
		const std::vector<d4> &A = res[0].maxima, &B = res[1].maxima;
		std::vector<char> b_used(B.size(), 0);
		std::cout << "\nELI-D alpha-alpha <-> beta-beta basins (same label; a repeated label by mutual nearest maxima within 1 bohr):\n"
			<< "  aa basin                bb basin              distance  N_alpha(aa)  N_beta(bb)   N_a - N_b\n";
		for (size_t a = 0; a < A.size(); a++) {
			double d, back;
			const int b = nearest(A[a], B, LB, LA[a], d);
			const bool unique = std::count(LA.begin(), LA.end(), LA[a]) == 1 && std::count(LB.begin(), LB.end(), LA[a]) == 1;
			const bool pair = b >= 0 && (unique || (d < 1.0 && nearest(B[b], A, LA, LA[a], back) == static_cast<int>(a)));
			std::cout << std::setw(4) << a + 1 << " " << std::left << std::setw(18) << res[0].labels[a] << std::right;
			if (pair) {
				b_used[b] = 1;
				std::cout << std::setw(5) << b + 1 << " " << std::left << std::setw(18) << res[1].labels[b] << std::right << std::fixed << std::setprecision(3) << std::setw(8) << d
					<< std::setprecision(4) << std::setw(13) << res[0].spin[a][0] << std::setw(12) << res[1].spin[b][1] << std::setw(12) << res[0].spin[a][0] - res[1].spin[b][1] << "\n";
			}
			else
				std::cout << std::setw(5) << "-" << " " << std::left << std::setw(18) << "(no partner)" << std::right << std::setw(8) << "-"
					<< std::fixed << std::setprecision(4) << std::setw(13) << res[0].spin[a][0] << std::setw(12) << "-" << std::setw(12) << res[0].spin[a][0] << "\n";
		}
		for (size_t b = 0; b < B.size(); b++)
			if (!b_used[b])
				std::cout << std::setw(4) << "-" << " " << std::left << std::setw(18) << "(no partner)" << std::right << std::setw(5) << b + 1 << " " << std::left << std::setw(18) << res[1].labels[b] << std::right
					<< std::setw(8) << "-" << std::setw(13) << "-" << std::fixed << std::setprecision(4) << std::setw(12) << res[1].spin[b][1] << std::setw(12) << -res[1].spin[b][1] << "\n";
	}
	std::cout << "\nELIA (antiparallel pairs): not computed. For a single determinant the on-top pair density is exactly"
		" rho_alpha * rho_beta, so ELIA carries no pair information; it needs a correlated 2-matrix (-eli_family prints the details)." << std::endl;
	if (!warn.empty()) std::cout << "  WARNING: " << warn << std::endl;
}

//-eli_basin_cube: ELI-D inside the "A-B bond" basins, zero outside, on a box cropped to them. The voxels climb the
//gridded field to a terminal (a grid maximum, or a box-face maximum), only the terminals and the voxels on the basin
//surface walk the analytic field with classify (the integrator's own climb), so the surface is the integrator's.
//That runs on a coarse box only (a basin touching a face grows that face); the requested spacing is reached by halving.
static void write_eli_basin_cube(WFN &l_w, const std::vector<d4> &maxima, const svec &labels, const svec &want, const double spacing_ang,
	const std::function<int(const d3 &)> &classify, std::ostream &log) {
	auto matches = [](std::string lab, const std::string &w) {
		if (lab == w) return true;
		lab.erase(std::remove_if(lab.begin(), lab.end(), [](const char c) { return std::isdigit(static_cast<unsigned char>(c)); }), lab.end());
		return lab == w;
	};
	std::vector<char> target(labels.size() + 2, 0);
	std::vector<int> hits;
	for (size_t b = 0; b < labels.size(); b++) {
		const std::string &l = labels[b];
		const size_t dash = l.find('-'), sp = l.find(" bond");
		if (dash == std::string::npos || sp == std::string::npos || sp < dash) continue;
		const std::string x = l.substr(0, dash), y = l.substr(dash + 1, sp - dash - 1);
		if ((matches(x, want[0]) && matches(y, want[1])) || (matches(x, want[1]) && matches(y, want[0]))) { target[b + 1] = 1; hits.push_back(static_cast<int>(b)); }
	}
	if (hits.empty()) { log << "\n-eli_basin_cube: no " << want[0] << "-" << want[1] << " bond basin, no cube written." << std::endl; return; }
	log << "\n-eli_basin_cube: " << hits.size() << " basin(s):";
	for (const int b : hits) log << " " << b + 1 << " (" << labels[b] << ")";
	log << std::endl;

	d3 lo{ 1e300, 1e300, 1e300 }, hi{ -1e300, -1e300, -1e300 };
	for (const int b : hits) for (int i = 0; i < 3; i++) { lo[i] = std::min(lo[i], maxima[b][i] - 3.0); hi[i] = std::max(hi[i], maxima[b][i] + 3.0); }
	std::array<int, 3> n{};
	vec val, rho;
	std::vector<char> in;
	double h = 0.1;
	auto idx = [&](const int x, const int y, const int z) { return (static_cast<size_t>(x) * n[1] + y) * n[2] + z; };
	auto pos = [&](const size_t i) { return d3{ lo[0] + (i / (static_cast<size_t>(n[1]) * n[2])) * h, lo[1] + ((i / n[2]) % n[1]) * h, lo[2] + (i % n[2]) * h }; };
	auto inside = [&](const d3 &p) { const int b = classify(p); return static_cast<char>(b > 0 && b < static_cast<int>(target.size()) && target[b]); };
	//Labels the box at spacing h; returns which faces (-x,+x,-y,+y,-z,+z) the basin touches
	//ponytail: a grown box re-evaluates the whole box; keep the old voxels if growth ever dominates
	auto label_box = [&]() {
		for (int i = 0; i < 3; i++) n[i] = static_cast<int>(std::ceil((hi[i] - lo[i]) / h)) + 1;
		const size_t N = static_cast<size_t>(n[0]) * n[1] * n[2];
		val.assign(N, 0.0); rho.assign(N, 0.0);
#pragma omp parallel for schedule(dynamic, 1024)
		for (long long i = 0; i < static_cast<long long>(N); i++) l_w.computeRhoELI(pos(i), rho[i], val[i]);
		//Steepest ascent on the grid: up = the neighbour with the largest rise per bohr, -1 a terminal, -2 below the density cutoff
		std::vector<long long> up(N, -2);
#pragma omp parallel for schedule(static)
		for (long long i = 0; i < static_cast<long long>(N); i++) {
			if (rho[i] < basin_density_cutoff) continue;
			const int x = static_cast<int>(i / (static_cast<long long>(n[1]) * n[2])), y = static_cast<int>((i / n[2]) % n[1]), z = static_cast<int>(i % n[2]);
			double best = 0.0;
			long long to = -1;
			for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
				const int a = x + dx, b = y + dy, c = z + dz;
				if ((!dx && !dy && !dz) || a < 0 || b < 0 || c < 0 || a >= n[0] || b >= n[1] || c >= n[2]) continue;
				const size_t j = idx(a, b, c);
				if (rho[j] < basin_density_cutoff) continue;
				const double rise = (val[j] - val[i]) / std::sqrt(static_cast<double>(dx * dx + dy * dy + dz * dz));
				if (rise > best) { best = rise; to = static_cast<long long>(j); }
			}
			up[i] = to;
		}
		//Every voxel to its terminal, with path compression
		std::vector<long long> root(N, -2), path;
		for (size_t i = 0; i < N; i++) {
			if (up[i] == -2 || root[i] != -2) continue;
			long long j = static_cast<long long>(i);
			path.clear();
			while (root[j] == -2 && up[j] >= 0) { path.push_back(j); j = up[j]; }
			const long long r = root[j] != -2 ? root[j] : j;
			root[j] = r;
			for (const long long k : path) root[k] = r;
		}
		std::vector<long long> terms;
		for (size_t i = 0; i < N; i++) if (root[i] == static_cast<long long>(i)) terms.push_back(static_cast<long long>(i));
		std::vector<char> t_in(terms.size(), 0);
#pragma omp parallel for schedule(dynamic)
		for (long long t = 0; t < static_cast<long long>(terms.size()); t++) t_in[t] = inside(pos(terms[t]));
		std::vector<char> term_in(N, 0);
		for (size_t t = 0; t < terms.size(); t++) term_in[terms[t]] = t_in[t];
		in.assign(N, 0);
		for (size_t i = 0; i < N; i++) if (root[i] >= 0) in[i] = term_in[root[i]];
		//The basin surface walks the analytic field until no voxel flips
		std::vector<char> done(N, 0);
		for (int pass = 0; pass < 200; pass++) {
			std::vector<long long> edge;
			for (int x = 0; x < n[0]; x++) for (int y = 0; y < n[1]; y++) for (int z = 0; z < n[2]; z++) {
				const size_t i = idx(x, y, z);
				if (done[i] || rho[i] < basin_density_cutoff) continue;
				const bool e = (x > 0 && in[idx(x - 1, y, z)] != in[i]) || (x + 1 < n[0] && in[idx(x + 1, y, z)] != in[i])
					|| (y > 0 && in[idx(x, y - 1, z)] != in[i]) || (y + 1 < n[1] && in[idx(x, y + 1, z)] != in[i])
					|| (z > 0 && in[idx(x, y, z - 1)] != in[i]) || (z + 1 < n[2] && in[idx(x, y, z + 1)] != in[i]);
				if (e) edge.push_back(static_cast<long long>(i));
			}
			if (edge.empty()) break;
			std::vector<char> e_in(edge.size(), 0);
#pragma omp parallel for schedule(dynamic)
			for (long long k = 0; k < static_cast<long long>(edge.size()); k++) e_in[k] = inside(pos(edge[k]));
			size_t flips = 0;
			for (size_t k = 0; k < edge.size(); k++) { done[edge[k]] = 1; flips += in[edge[k]] != e_in[k]; in[edge[k]] = e_in[k]; }
			log << "  h " << std::setprecision(4) << h << " bohr, pass " << pass + 1 << ": " << terms.size() << " grid terminals, " << edge.size() << " surface voxels walked, " << flips << " flipped" << std::endl;
			if (!flips) break;
		}
		std::array<bool, 6> face{};
		for (size_t i = 0; i < N; i++) {
			if (!in[i]) continue;
			const int x = static_cast<int>(i / (static_cast<size_t>(n[1]) * n[2])), y = static_cast<int>((i / n[2]) % n[1]), z = static_cast<int>(i % n[2]);
			face[0] |= x == 0; face[1] |= x == n[0] - 1; face[2] |= y == 0; face[3] |= y == n[1] - 1; face[4] |= z == 0; face[5] |= z == n[2] - 1;
		}
		return face;
	};
	//The whole box is labelled only at 4 fine spacings. Each halving gives a new point the label of its enclosing
	//cell when all its corners agree and walks the analytic field otherwise, so only the cells the surface cuts are
	//climbed, and ELI-D is evaluated only inside.
	//ponytail: a protrusion thinner than a coarse cell (0.076 bohr at 0.01 A) between agreeing corners is lost
	const int levels = 2, stride = 1 << levels;
	const double hf = constants::ang2bohr(spacing_ang);
	h = hf * stride;
	for (int grow = 0; grow < 6; grow++) {
		const std::array<bool, 6> face = label_box();
		bool any = false;
		for (int i = 0; i < 3; i++) {
			if (face[2 * i]) { lo[i] -= 1.0; any = true; }
			if (face[2 * i + 1]) { hi[i] += 1.0; any = true; }
		}
		if (!any) break;
		log << "  the basin reaches the box face, growing it by 1 bohr" << std::endl;
	}
	//The fine grid spans the coarse basin's index box plus one coarse cell, aligned so every coarse voxel is a fine one
	std::array<int, 3> c0{ n[0], n[1], n[2] }, c1{ -1, -1, -1 };
	for (int x = 0; x < n[0]; x++) for (int y = 0; y < n[1]; y++) for (int z = 0; z < n[2]; z++)
		if (in[idx(x, y, z)]) { const int c[3] = { x, y, z }; for (int k = 0; k < 3; k++) { c0[k] = std::min(c0[k], c[k]); c1[k] = std::max(c1[k], c[k]); } }
	err_checkf(c1[0] >= 0, "-eli_basin_cube: the basin holds no voxel", log);
	for (int k = 0; k < 3; k++) { c0[k] = std::max(0, c0[k] - 1); c1[k] = std::min(n[k] - 1, c1[k] + 1); }
	const std::array<int, 3> nc = n;
	const std::vector<char> coarse = std::move(in);
	for (int k = 0; k < 3; k++) { lo[k] += c0[k] * h; n[k] = (c1[k] - c0[k]) * stride + 1; }
	h = hf;
	const size_t NF = static_cast<size_t>(n[0]) * n[1] * n[2];
	std::vector<signed char> lab(NF, -1);
	for (int x = c0[0]; x <= c1[0]; x++) for (int y = c0[1]; y <= c1[1]; y++) for (int z = c0[2]; z <= c1[2]; z++)
		lab[idx((x - c0[0]) * stride, (y - c0[1]) * stride, (z - c0[2]) * stride)] = coarse[(static_cast<size_t>(x) * nc[1] + y) * nc[2] + z];
	for (int s = stride / 2; s >= 1; s /= 2) {
		const int S = 2 * s;
		std::vector<long long> pend;
		for (int x = 0; x < n[0]; x += s) for (int y = 0; y < n[1]; y += s) for (int z = 0; z < n[2]; z += s) {
			const size_t i = idx(x, y, z);
			if (lab[i] >= 0) continue;
			const int x0 = x / S * S, y0 = y / S * S, z0 = z / S * S, ex = x % S ? S : 0, ey = y % S ? S : 0, ez = z % S ? S : 0;
			int ones = 0, tot = 0;
			for (int dx = 0; dx <= ex; dx += S) for (int dy = 0; dy <= ey; dy += S) for (int dz = 0; dz <= ez; dz += S) { tot++; ones += lab[idx(x0 + dx, y0 + dy, z0 + dz)]; }
			if (ones == tot) lab[i] = 1;
			else if (ones == 0) lab[i] = 0;
			else pend.push_back(static_cast<long long>(i));
		}
#pragma omp parallel for schedule(dynamic)
		for (long long k = 0; k < static_cast<long long>(pend.size()); k++) lab[pend[k]] = inside(pos(pend[k]));
		log << "  h " << std::setprecision(4) << h * s << " bohr: " << pend.size() << " points in cells cut by the surface walked the field" << std::endl;
	}
	val.assign(NF, 0.0); rho.assign(NF, 0.0);
#pragma omp parallel for schedule(dynamic, 1024)
	for (long long i = 0; i < static_cast<long long>(NF); i++) if (lab[i] == 1) l_w.computeRhoELI(pos(i), rho[i], val[i]);
	in.assign(lab.begin(), lab.end());
	double vol = 0.0, pop = 0.0;
	const double dv = h * h * h;
	cube out(n, l_w.get_ncen(), true);
	out.give_parent_wfn(l_w);
	for (int i = 0; i < 3; i++) { out.set_origin(i, lo[i]); for (int j = 0; j < 3; j++) out.set_vector(i, j, i == j ? h : 0.0); }
	out.calc_dv();
	for (int x = 0; x < n[0]; x++) for (int y = 0; y < n[1]; y++) for (int z = 0; z < n[2]; z++) {
		const size_t i = idx(x, y, z);
		out.set_value(x, y, z, in[i] ? val[i] : 0.0);
		if (in[i]) { vol += dv; pop += rho[i] * dv; }
	}
	std::string names;
	for (const int b : hits) names += (names.empty() ? "" : ", ") + labels[b];
	out.set_comment1("ELI-D inside the ELI-D basin(s) " + names + ", zero outside");
	out.set_comment2("NoSpherA2 -eli_basin_cube " + want[0] + " " + want[1] + ", spacing " + std::to_string(spacing_ang) + " Angstrom");
	std::filesystem::path p = l_w.get_path();
	p.replace_filename(p.stem().string() + "_eli_" + want[0] + "-" + want[1] + ".cube");
	out.set_path(p);
	log << "  " << n[0] << " x " << n[1] << " x " << n[2] << " voxels of " << std::setprecision(4) << spacing_ang << " A: volume " << vol << " bohr^3, "
		<< pop << " electrons (voxel sums; the quadrature table above is the reference)" << std::endl;
	log << "Writing " << p.string() << std::endl;
	err_checkf(out.write_file(true), "-eli_basin_cube: failed to write the cube", log);
}

//-eli_mask_cube: ELI-D on the nonzero voxels of a given grid, zero elsewhere - a region carried over from another
//complex where this one has no basin. The region cube's atoms are not read (a parent without atoms takes none).
static void write_eli_mask_cube(WFN &l_w, const std::filesystem::path &mask_path, std::ostream &log) {
	WFN none(e_origin::cub);
	const cube mask(mask_path, true, none, log);
	const i3 n = mask.get_sizes();
	cube out(n, l_w.get_ncen(), true);
	out.give_parent_wfn(l_w);
	d3 o;
	std::array<d3, 3> v;
	for (int i = 0; i < 3; i++) { o[i] = mask.get_origin(i); out.set_origin(i, o[i]); for (int j = 0; j < 3; j++) { v[i][j] = mask.get_vector(i, j); out.set_vector(i, j, v[i][j]); } }
	out.calc_dv();
	const double dv = out.get_dv();
	double vol = 0.0, pop = 0.0;
#pragma omp parallel for reduction(+:vol, pop) schedule(dynamic)
	for (int x = 0; x < n[0]; x++) for (int y = 0; y < n[1]; y++) for (int z = 0; z < n[2]; z++) {
		if (mask.get_value(x, y, z) == 0.0) continue;
		d3 p;
		for (int k = 0; k < 3; k++) p[k] = o[k] + x * v[0][k] + y * v[1][k] + z * v[2][k];
		double rho, eli;
		l_w.computeRhoELI(p, rho, eli);
		out.set_value(x, y, z, eli);
		vol += dv; pop += rho * dv;
	}
	out.set_comment1("ELI-D inside the region of " + mask_path.filename().string() + ", zero outside");
	out.set_comment2("NoSpherA2 -eli_mask_cube");
	std::filesystem::path p = l_w.get_path();
	p.replace_filename(p.stem().string() + "_eli_in_" + mask_path.stem().string() + ".cube");
	out.set_path(p);
	log << "-eli_mask_cube: " << n[0] << " x " << n[1] << " x " << n[2] << " voxels, region volume " << vol << " bohr^3, " << pop << " electrons" << std::endl;
	log << "Writing " << p.string() << std::endl;
	err_checkf(out.write_file(true), "-eli_mask_cube: failed to write the cube", log);
}

namespace {
//One column of the -eli_contributions table: c[b][k] electrons of ELI-D basin b from contributor k
struct contrib_column {
	std::string title, skipped;
	vec2 c;
	svec lab;
};
}

//p_k(b) = 2 sum_ij d_i d_j S_ij(b) for orbital k = sum_m coef[k][m] MO m of wavy; exact for a closed-shell
//determinant, so orbitals spanning the occupied space split each basin's population without remainder
static vec2 orbital_contributions(const WFN &wavy, const basin_overlaps &ovl, const vec2 &coef) {
	//wavy MO -> l_w MO (delete_unoccupied_MOs drops occ == 0 and keeps the order) -> slot in the overlaps
	ivec slot(wavy.get_nmo(), -1), occupied;
	for (int m = 0; m < wavy.get_nmo(); m++)
		if (wavy.get_MO_occ(m) != 0.0) occupied.push_back(m);
	for (int i = 0; i < ovl.nmo; i++) slot[occupied[ovl.mo_index[i]]] = i;
	const int nb = static_cast<int>(ovl.S.size()), nk = static_cast<int>(coef.size());
	vec2 c(nb, vec(nk, 0.0));
#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < nk; k++) {
		vec d(ovl.nmo, 0.0);
		for (int m = 0; m < wavy.get_nmo(); m++)
			if (slot[m] >= 0) d[slot[m]] = coef[k][m];
		for (int b = 0; b < nb; b++) {
			double s = 0.0;
			for (int i = 0; i < ovl.nmo; i++) {
				if (d[i] == 0.0) continue;
				const double *row = ovl.S[b].data() + (size_t)i * (i + 1) / 2;
				double t = 0.0;
				for (int j = 0; j < i; j++) t += d[j] * row[j];
				s += d[i] * (2.0 * t + d[i] * row[i]);
			}
			c[b][k] = 2.0 * s;
		}
	}
	return c;
}

//Per ELI-D basin its population, then each column's contributors of at least 0.05 e, largest first, and what
//the column leaves unlisted; the listed ones and the rest add up to the population
static void print_eli_contributions(const vec &pop, const svec &eli_labels, const std::vector<contrib_column> &cols, std::ostream &out) {
	constexpr double show = 0.05;
	constexpr int w = 22;
	out << "\nELI-D basin contributions (electrons; contributors >= " << show << " e, largest first)\n"
		<< "  atoms numbered from 0 in every column, as in the QTAIM and ELI-D tables; (rest) is what a column leaves unlisted\n";
	for (const contrib_column &col : cols)
		if (!col.skipped.empty()) out << "  " << col.title << " skipped: " << col.skipped << "\n";
	out << "  " << std::left << std::setw(30) << "ELI-D basin";
	for (const contrib_column &col : cols)
		if (col.skipped.empty()) out << " | " << std::setw(w) << col.title;
	out << std::right << "\n";
	for (size_t b = 0; b < pop.size(); b++) {
		std::vector<svec> cells;
		size_t rows = 0;
		for (const contrib_column &col : cols) {
			if (!col.skipped.empty()) continue;
			//an atom's core orbitals are one entry, "CR U1 x43"; a heavy core would fill the column otherwise
			std::vector<std::pair<std::string, double>> entry;
			std::map<std::string, std::pair<size_t, int>> core;
			double listed = 0.0, sum = 0.0;
			for (size_t k = 0; k < col.c[b].size(); k++) {
				sum += col.c[b][k];
				if (col.lab[k].rfind("CR ", 0) != 0) { entry.emplace_back(col.lab[k], col.c[b][k]); continue; }
				auto it = core.find(col.lab[k]);
				if (it == core.end()) { core[col.lab[k]] = { entry.size(), 1 }; entry.emplace_back(col.lab[k], col.c[b][k]); }
				else { entry[it->second.first].second += col.c[b][k]; it->second.second++; }
			}
			for (const auto &[name, at] : core)
				if (at.second > 1) entry[at.first].first += " x" + std::to_string(at.second);
			// rank on the printed 3 decimals, label on ties, so goldens do not swap equal rows between platforms
			std::sort(entry.begin(), entry.end(), [](const auto &x, const auto &y) {
				const long long a = std::llround(x.second * 1000.0), b = std::llround(y.second * 1000.0);
				return a != b ? a > b : x.first < y.first;
				});
			svec cell;
			for (const auto &[name, v] : entry) {
				if (v < show) break;
				std::ostringstream s;
				s << std::left << std::setw(w - 8) << name << std::right << std::fixed << std::setprecision(3) << std::setw(8) << v;
				cell.push_back(s.str());
				listed += v;
			}
			std::ostringstream rest;
			rest << std::left << std::setw(w - 8) << "(rest)" << std::right << std::fixed << std::setprecision(3) << std::setw(8) << sum - listed;
			cell.push_back(rest.str());
			rows = std::max(rows, cell.size());
			cells.push_back(cell);
		}
		std::ostringstream head;
		head << std::setw(4) << b + 1 << " " << std::left << std::setw(16) << eli_labels[b] << std::right << std::fixed << std::setprecision(3) << std::setw(9) << pop[b];
		for (size_t r = 0; r < rows; r++) {
			out << "  " << std::left << std::setw(30) << (r == 0 ? head.str() : "");
			for (const svec &cell : cells) {
				//the rest rows line up at the bottom of every column
				const size_t lead = rows - cell.size();
				out << " | " << std::setw(w) << (r < lead ? "" : cell[r - lead]);
			}
			out << std::right << "\n";
		}
		out << "  " << std::string(30 + cells.size() * (3 + w), '-') << "\n";
	}
}

void ELI_analysis(const WFN &wavy, options &opt) {
	err_checkf(wavy.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
	std::cout << (opt.qtaim_only ? "Analysing QTAIM basins" : "Analysing ELI basins") << " in the wavefunction..." << std::endl;
	if (!opt.qtaim_only) citations::cite(citations::Method::ELID, std::cout);
	density_field field;
	const std::unique_ptr<Gaussian_Molecule> fit = fitted_source(wavy, opt, field, std::cout);
	const density_field *fld = fit ? &field : nullptr;
	//A fitted source has no orbitals, so ELI comes from the density alone.
	if (fld)
		citations::cite(citations::Method::ELIOrbitalFree, std::cout);

	const double radius = opt.properties.radius;
	const double grid_spacing = opt.properties.resolution;
	properties_options prop_opt = opt.properties;
	WFN l_w = wavy;
	l_w.delete_unoccupied_MOs();
	if (!opt.eli_mask_cube.empty()) { write_eli_mask_cube(l_w, opt.eli_mask_cube, std::cout); return; }
	//g = rho tau - |grad rho|^2 / 4 vanishes identically when one orbital carries the density: ELI-D is 0/0 and
	//its basins are round-off (DGrid shatters the same way). That is the definition, not the integrator, so warn.
	if (l_w.get_nmo() < 2 && !opt.qtaim_only)
		std::cout << "WARNING: this wavefunction has a single occupied orbital, so ELI-D's pair"
			" density g = rho*tau - |grad rho|^2/4 is identically zero and the field is undefined."
			" The basins below are round-off structure, not chemistry - expect thousands of them and"
			" do not quote their populations. The QTAIM basins are unaffected." << std::endl;
	//By default both basin sets are found and integrated on the analytic field (critical-point search,
	//analytic_eli_maxima, quadrature points walked up the field). The cube is built only for -basin_cube,
	//a fitted density (its attractors are not the orbitals' critical points) and a wavefunction without orbitals.
	const bool stream_qtaim = !opt.basin_cube && !fld && l_w.get_nmo() > 0;
	const bool stream_eli = !opt.basin_cube;
	if (!opt.eli_basin_atoms.empty() && !stream_eli) std::cout << "WARNING: -eli_basin_cube needs the analytic ELI-D basins and is ignored with -basin_cube." << std::endl;
	const bool need_cube = (!stream_qtaim && !opt.eli_only) || (!stream_eli && !opt.qtaim_only);
	const std::vector<atom> atoms = wavy.get_atoms();
	cube rho, eli_cube;
	basin_stage_timer T;
	if (need_cube) {
		readxyzMinMax_fromWFN(wavy, prop_opt);
		rho = cube(prop_opt.NbSteps, l_w.get_ncen(), true);
		eli_cube = cube(prop_opt.NbSteps, l_w.get_ncen(), true);
		rho.give_parent_wfn(l_w);
		eli_cube.give_parent_wfn(l_w);
		std::cout << "Calcualting grid of size " << prop_opt.NbSteps[0] << " x " << prop_opt.NbSteps[1] << " x " << prop_opt.NbSteps[2] << "..." << std::endl;
		std::cout << "Number of points: " << prop_opt.NbSteps[0] * prop_opt.NbSteps[1] * prop_opt.NbSteps[2] << std::endl;
		std::cout << "Grid parameters:\n";
		std::cout << "  Min: (" << prop_opt.MinMax[0] << ", " << prop_opt.MinMax[1] << ", " << prop_opt.MinMax[2] << ")\n";
		std::cout << "  Max: (" << prop_opt.MinMax[3] << ", " << prop_opt.MinMax[4] << ", " << prop_opt.MinMax[5] << ")\n";
		vec stepsizes{ (prop_opt.MinMax[3] - prop_opt.MinMax[0]) / prop_opt.NbSteps[0],
					  (prop_opt.MinMax[4] - prop_opt.MinMax[1]) / prop_opt.NbSteps[1],
					  (prop_opt.MinMax[5] - prop_opt.MinMax[2]) / prop_opt.NbSteps[2] };
		for (int i = 0; i < 3; i++) {
			rho.set_origin(i, prop_opt.MinMax[i]);
			rho.set_vector(i, i, stepsizes[i]);
			eli_cube.set_origin(i, prop_opt.MinMax[i]);
			eli_cube.set_vector(i, i, stepsizes[i]);
		}
		rho.calc_dv();
		eli_cube.calc_dv();
		Calc_RhoEli(rho, eli_cube, l_w, radius, fld);
		T.lap("rho and ELI-D cube");
	}

	// print table of atom positions
	std::cout << "Atom positions:\n";
	for (int a = 0; a < atoms.size(); a++) {
		std::cout << "  Atom " << a << ": (" << atoms[a].get_coordinate(0) << ", " << atoms[a].get_coordinate(1) << ", " << atoms[a].get_coordinate(2) << ")\n";
	}

	//An ECP took the core electrons out of the density. The QTAIM basins get them back from
	//Thakkar's spherical core densities, the fill the Hirshfeld grids and the scattering
	//factors apply: the nucleus is a cusp again and its basin holds the atom's full count.
	//ELI-D stays on the valence density the wavefunction has.
	ivec ecp_atoms;
	std::vector<Thakkar> ecp_cores;
	double ecp_electrons = 0.0;
	for (int a = 0; a < l_w.get_ncen(); a++)
		if (l_w.get_atom_ECP_electrons(a) > 0) {
			const int mode = l_w.get_ECP_mode() > 0 ? l_w.get_ECP_mode() : 1;
			ecp_atoms.push_back(a);
			ecp_cores.emplace_back(l_w.get_atom_charge(a), mode);
			ecp_electrons += l_w.get_atom_ECP_electrons(a);
		}
	std::function<double(const d3&)> core_density = [&](const d3 &p) {
		double s = 0.0;
		for (size_t k = 0; k < ecp_atoms.size(); k++) {
			const d3 ap = l_w.get_atom_pos(ecp_atoms[k]);
			const double d = std::sqrt(std::pow(p[0] - ap[0], 2) + std::pow(p[1] - ap[1], 2) + std::pow(p[2] - ap[2], 2));
			s += ecp_cores[k].get_core_density(d, l_w.get_atom_ECP_electrons(ecp_atoms[k]));
		}
		return s;
	};
	std::function<void(const d3&, d3&)> core_gradient = [&](const d3 &p, d3 &g) {
		const double h = 1e-4;
		for (int k = 0; k < 3; k++) {
			d3 a = p, b = p;
			a[k] += h; b[k] -= h;
			g[k] = (core_density(a) - core_density(b)) / (2.0 * h);
		}
	};
	const bool fill_cores = !ecp_atoms.empty();
	if (fill_cores) {
		std::cout << "ECP cores of " << ecp_atoms.size() << " atoms filled with Thakkar densities: " << std::fixed << std::setprecision(1) << ecp_electrons << " electrons added for the QTAIM basins" << std::endl;
		for (int x = 0; x < rho.get_size(0); x++)
			for (int y = 0; y < rho.get_size(1); y++)
				for (int z = 0; z < rho.get_size(2); z++)
					if (rho.get_value(x, y, z) > 0.0) rho.set_value(x, y, z, rho.get_value(x, y, z) + core_density(rho.get_pos(x, y, z)));
	}
	if (need_cube && opt.debug) {
		rho.set_path("rho.cube");
		eli_cube.set_path("eli.cube");
		rho.write_file(true);
		eli_cube.write_file(true);
	}

	//The critical points are refined on the orbitals (V, G and K need them) even when the
	//basins follow the fitted density; a SALTED prediction has none, so they are skipped. -topology's Newton
	//search on the analytic Hessian, checked against Poincare-Hopf: no grid decides which points exist
	std::vector<critical_point> density_critical_points;
	if (opt.eli_only) std::cout << "-eli_only: density critical points and QTAIM basins skipped." << std::endl;
	else if (l_w.get_nmo() > 0) {
		//-topology's default 1.3 x (r_cov + r_cov) bond seeding misses hydrogen bonds and long contacts while
		//Poincare-Hopf still closes; 2.5 reaches O-H...O and C-H...O, and seeding and fragment count share it
		topology::options topo_opt;
		topo_opt.bond_scale = 2.5;
		const topology::result top = topology::analyze_topology(l_w, topology::nuclei_of(l_w), topo_opt);
		for (const topology::cp &p : top.points) {
			if (p.kind != topology::cp_kind::attractor && p.density < basin_density_cutoff) continue;
			//a saddle within a bohr of an ECP-cored nucleus sits in the pseudo-density's core hole and is no bond
			if (p.kind != topology::cp_kind::attractor && p.nearest_nucleus_distance < 1.0
				&& std::find(top.coreless_nuclei.begin(), top.coreless_nuclei.end(), p.nearest_nucleus) != top.coreless_nuclei.end()) continue;
			const critical_point_seed seed{ i3{ 0, 0, 0 }, p.position, p.density, p.gradient_norm, p.kind == topology::cp_kind::attractor && !p.is_nna, p.nearest_nucleus };
			density_critical_points.push_back(evaluate_critical_point(seed, p.position, l_w, p.iterations, true));
		}
		const section_log::tee cp_tee("qtaim", "Critical points of the electron density");
		std::cout << "Density critical points from the analytic field: NCP " << top.n_attractor << ", BCP " << top.n_bond << ", RCP "
			<< top.n_ring << ", CCP " << top.n_cage << "; Poincare-Hopf sum " << top.sum << " against " << top.target
			<< (top.complete ? " (COMPLETE)" : " (INCOMPLETE: " + top.diagnosis + ")") << std::endl;
	}
	else std::cout << "No orbitals: critical points (Hessian, V, G, K) need a wavefunction and are skipped." << std::endl;
	T.lap("density critical points");
	//Core shells make critical points of their own and an ECP atom a whole sphere of them,
	//none of which says anything about bonding and none of which any two machines find at
	//the same spots; only the nuclear attractor survives inside an atom's core radius. Sorted
	//by type, density and position so the listing reads the same everywhere.
	{
		std::vector<critical_point> kept;
		for (const critical_point &cp : density_critical_points) {
			bool core = false, nuclear = false;
			for (int a = 0; a < l_w.get_ncen(); a++) {
				const d3 apos = l_w.get_atom_pos(a);
				const double d2 = std::pow(cp.position[0] - apos[0], 2) + std::pow(cp.position[1] - apos[1], 2) + std::pow(cp.position[2] - apos[2], 2);
				if (d2 < 0.01) nuclear = true;
				else if (d2 < std::pow(core_shell_radius(l_w.get_atom_charge(a)), 2)) core = true;
			}
			if (!core || nuclear) kept.push_back(cp);
		}
		std::sort(kept.begin(), kept.end(), [](const critical_point &a, const critical_point &b) {
			if (a.type != b.type) return a.type < b.type;
			if (std::abs(a.density - b.density) > 1e-6 * std::max(1.0, std::abs(a.density))) return a.density > b.density;
			for (int k = 0; k < 3; k++)
				if (std::abs(a.position[k] - b.position[k]) > 1e-4) return a.position[k] < b.position[k];
			return false;
		});
		density_critical_points.swap(kept);
	}
	if (!opt.eli_only) {
	const section_log::tee cp_tee("qtaim");
	std::cout << "Density Critical Points";
	if (!density_critical_points.empty())
		std::cout << " (" << density_critical_points.size() << " found)";
	std::cout << ":\n";
	if (density_critical_points.empty()) {
		std::cout << "  No critical points found.\n";
	}
	else {
		constexpr int nw = 15; // numeric field width — wide enough for large core densities
		std::vector<topology::bond_path> bond_paths(density_critical_points.size());
		{
			const std::vector<topology::nucleus> nuclei = topology::nuclei_of(l_w);
#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < (int)density_critical_points.size(); i++)
				if (density_critical_points[i].type == "bond")
					bond_paths[i] = topology::trace_bond_path(l_w, density_critical_points[i].position, nuclei);
		}
		if (opt.fba) write_bond_path_profiles(l_w, bond_paths, opt.wfn.string(), nullptr); //beside the straight-line profiles
		for (size_t i = 0; i < density_critical_points.size(); i++) {
			const critical_point &cp = density_critical_points[i];

			std::string ownership_info;
			if (cp.type == "attractor" || cp.type == "bond") {
				double min_dist1 = std::numeric_limits<double>::max();
				double min_dist2 = std::numeric_limits<double>::max();
				int atom_index1 = -1;
				int atom_index2 = -1;
				for (int a = 0; a < l_w.get_ncen(); a++) {
					const d3 apos = l_w.get_atom_pos(a);
					const double dx = cp.position[0] - apos[0];
					const double dy = cp.position[1] - apos[1];
					const double dz = cp.position[2] - apos[2];
					const double dist2 = dx * dx + dy * dy + dz * dz;
					if (dist2 < min_dist1) {
						min_dist2 = min_dist1;
						atom_index2 = atom_index1;
						min_dist1 = dist2;
						atom_index1 = a;
					}
					else if (dist2 < min_dist2) {
						min_dist2 = dist2;
						atom_index2 = a;
					}
				}

				if (cp.type == "attractor" && atom_index1 >= 0) {
					ownership_info = "   " + l_w.get_atom_label(atom_index1) + std::to_string(atom_index1);
				}
				else if (cp.type == "bond" && atom_index1 >= 0 && atom_index2 >= 0) {
					ownership_info = "   "
						+ l_w.get_atom_label(atom_index1) + std::to_string(atom_index1)
						+ "-"
						+ l_w.get_atom_label(atom_index2) + std::to_string(atom_index2)
						+ " bond";
				}
			}

			std::cout << "\n  CP " << std::right << std::setw(3) << i + 1
				<< "  [" << std::left << std::setw(12) << cp.type << "]  "
				<< (cp.converged ? "converged" : "NOT converged")
				<< ownership_info << "\n";
			std::cout << std::fixed << std::setprecision(4) << std::right;
			std::cout << "    Position  :"
				<< std::setw(nw) << cp.position[0]
				<< std::setw(nw) << cp.position[1]
				<< std::setw(nw) << cp.position[2] << "\n";
			std::cout << "    Rho       :" << std::setw(nw) << cp.density << "\n";
			std::cout << "    GradRho   :"
				<< std::setw(nw) << cp.gradient[0]
				<< std::setw(nw) << cp.gradient[1]
				<< std::setw(nw) << cp.gradient[2]
				<< "   |GradRho| :" << std::setw(nw) << cp.gradient_norm << "\n";
			std::cout << "    HessRho_EigVals:"
				<< std::scientific << std::setprecision(4)
				<< std::setw(nw) << cp.hessian_eigenvalues[0]
				<< std::setw(nw) << cp.hessian_eigenvalues[1]
				<< std::setw(nw) << cp.hessian_eigenvalues[2]
				<< std::fixed << std::setprecision(4) << "\n";
			//The eigenvectors of a degenerate pair are any two in their plane and their signs
			//are free; both differ from machine to machine, so they are for -debug
			if (opt.debug) {
				std::cout << "    HessRho_EigVecs v1:"
					<< std::setw(nw) << cp.hessian_eigenvectors[0][0]
					<< std::setw(nw) << cp.hessian_eigenvectors[0][1]
					<< std::setw(nw) << cp.hessian_eigenvectors[0][2] << "\n";
				std::cout << "                    v2:"
					<< std::setw(nw) << cp.hessian_eigenvectors[1][0]
					<< std::setw(nw) << cp.hessian_eigenvectors[1][1]
					<< std::setw(nw) << cp.hessian_eigenvectors[1][2] << "\n";
				std::cout << "                    v3:"
					<< std::setw(nw) << cp.hessian_eigenvectors[2][0]
					<< std::setw(nw) << cp.hessian_eigenvectors[2][1]
					<< std::setw(nw) << cp.hessian_eigenvectors[2][2] << "\n";
			}
			std::cout << "    DelSqRho  :" << std::setw(nw) << cp.laplacian << "\n";
			if (std::isfinite(cp.ellipticity))
				std::cout << "    Bond Ellipticity:" << std::setw(nw) << cp.ellipticity << "\n";
			if (const topology::bond_path &bp = bond_paths[i]; bp.a >= 0) {
				//the traced path's ends, which for a long contact need not be the two nearest atoms
				std::cout << "    Bond Path :   " << l_w.get_atom_label(bp.a) << bp.a << "-" << l_w.get_atom_label(bp.b) << bp.b
					<< "   BPL :" << std::setw(nw) << bp.bpl << "   GPL :" << std::setw(nw) << bp.gpl
					<< "   BPL-GPL :" << std::setw(nw) << bp.bpl - bp.gpl << " bohr";
				if (bp.gap_a > 0.0 || bp.gap_b > 0.0)
					std::cout << "   end maxima off nuclei by " << bp.gap_a << " / " << bp.gap_b;
				std::cout << "\n";
			}
			if (std::isfinite(cp.virial_field) || std::isfinite(cp.kinetic_lagrangian) || std::isfinite(cp.kinetic_hamiltonian) || std::isfinite(cp.lagrangian_density)) {
				std::cout << "   V         :" << std::scientific << std::setprecision(4) << std::setw(nw) << cp.virial_field
					<< "   G         :" << std::fixed << std::setprecision(4) << std::setw(nw) << cp.kinetic_lagrangian
					<< "   K         :" << std::scientific << std::setprecision(4) << std::setw(nw) << cp.kinetic_hamiltonian
					<< "   L         :" << std::scientific << std::setprecision(4) << std::setw(nw) << cp.lagrangian_density
					<< std::fixed << std::setprecision(4) << "\n";
			}
		}
		std::cout << "\n";
	}
	}

	//Every nucleus is a maximum of the density, whatever the grid says
	std::vector<d3> nuclei;
	for (const atom &a : atoms) nuclei.push_back(a.get_pos());
	//A fitted density oscillates about zero in the tail, and every positive island out there is
	//a maximum with no neighbour to merge into (2000 of them on epoxide at 0.05 A, and the
	//persistence merge is quadratic in their number), so the search stops where the fit is no
	//longer trusted; its ~1e-3 error at a bond critical point also leaves a bump there that
	//5e-3 persistence keeps, hence the looser merge
	const double floor = basin_density_cutoff, persistence = fld ? 2e-2 : 5e-3;
	//Density attractors come from the critical-point search already run, and the quadrature walks the field.
	//A fitted density keeps the cube: those critical points are the orbitals', not the fit's attractors.
	std::pair<cubei, std::vector<d4>> qtaim_results;
	if (opt.eli_only) {}
	else if (stream_qtaim) {
		qtaim_results.second = streaming_density_attractors(l_w, density_critical_points, fill_cores ? &core_density : nullptr, fill_cores ? &core_gradient : nullptr, opt.debug);
		std::cout << "QTAIM attractors from the analytic field: " << qtaim_results.second.size() << " (" << l_w.get_ncen() << " nuclei, " << qtaim_results.second.size() - l_w.get_ncen() << " non-nuclear)" << std::endl;
	}
	else
		qtaim_results = topological_cube_analysis(&rho, atoms, opt.debug, true, floor, 1e-10, radius, persistence, &nuclei, &l_w, fill_cores ? &core_density : nullptr, fill_cores ? &core_gradient : nullptr, fld);
	T.lap("QTAIM attractors");
	svec labels = assign_labels_to_basins(qtaim_results.second, atoms, opt.debug);

	//Two integrations of the density over each basin set: the voxel sum, which is what the cube
	//resolution buys, and the atom-centred quadrature grids with the boundary decided by the
	//field itself, which is the number to compare with AIMAll and DGrid. The ELI-D basins follow
	//the orbitals' ELI-D and integrate the orbital density; only the QTAIM set uses the fit
	//A streaming ELI-D walks to every core shell's own maximum while the report keeps one core basin per atom,
	//so the integrator gets the unmerged maxima and the map to the basins
	std::vector<d4> eli_maxima_all;
	ivec eli_core_map;
	std::vector<shell_split> eli_split;
	//-eli_contributions: the ELI-D integration also keeps its overlap matrices and every point's share
	basin_overlaps eli_ovl;
	std::vector<basin_point> eli_points;
	vec eli_pop;
	bool contrib_orbitals = false;
	auto report = [&](const char *title, const std::pair<cubei, std::vector<d4>> &res, svec &lab, const bool eli, const bool stream) {
		//The voxel sum needs the basin cube; streaming has none
		if (!stream) {
			std::cout << "\n" << title << " (voxel sum):\n";
			integrate_values_in_basins(&rho, &(res.first), lab, opt.debug);
			T.lap(std::string(title) + " voxel sum");
		}
		vec vol;
		double outside = 0.0;
		//The overlap matrices cost one triangle per basin per thread; past the limit the delocalization indices are dropped, not the run
		basin_overlaps ovl;
		const bool orbitals = !eli && !fld && l_w.get_nmo() > 0;
		const size_t aom_bytes = orbitals ? (size_t)l_w.get_nmo() * (l_w.get_nmo() + 1) / 2 * res.second.size() * omp_get_max_threads() * sizeof(double) : 0;
		const bool want_aom = orbitals && aom_bytes < (size_t)512 * 1024 * 1024;
		if (orbitals && !want_aom)
			std::cout << "  Delocalization indices skipped: " << l_w.get_nmo() << " orbitals over " << res.second.size()
				<< " basins on " << omp_get_max_threads() << " threads would need " << aom_bytes / (1024 * 1024) << " MB of overlap matrices.\n";
		const bool mapped = eli && stream && !eli_core_map.empty();
		//-eli_basin_cube labels its voxels with the integrator's own climb, after the integration has filled the memo
		const std::function<void(const std::function<int(const d3 &)> &)> probe = [&](const std::function<int(const d3 &)> &classify) {
			write_eli_basin_cube(l_w, res.second, lab, opt.eli_basin_atoms, grid_spacing, classify, std::cout);
		};
		const bool want_probe = eli && stream && opt.eli_basin_atoms.size() == 2;
		//-eli_only with -eli_basin_cube: the cube alone, the basin populations are not integrated
		const bool cube_only = want_probe && opt.eli_only;
		const bool contrib = eli && opt.eli_contributions && !cube_only;
		const vec pop = integrate_basins_on_atomic_grids(stream ? nullptr : &rho, stream ? nullptr : &(res.first), mapped ? eli_maxima_all : res.second, l_w, opt.accuracy, eli, vol, outside, fill_cores && !eli ? &core_density : nullptr, fill_cores && !eli ? &core_gradient : nullptr, opt.basin_grid, eli ? nullptr : fld, want_aom ? &ovl : (contrib && contrib_orbitals ? &eli_ovl : nullptr), mapped ? &eli_core_map : nullptr, nullptr, nullptr, mapped ? &eli_split : nullptr, want_probe ? &probe : nullptr, cube_only, contrib ? &eli_points : nullptr);
		if (cube_only) { T.lap("ELI-D basin cube"); return; }
		if (eli) eli_pop = pop;
		const section_log::tee basin_tee(eli ? "eli" : "qtaim", eli ? "ELI-D basins: populations, volumes and maxima"
			: "QTAIM basins: populations, charges, volumes and maxima");
		std::cout << "\n" << title << " (atomic quadrature grids):\n";
		if (!eli) citations::cite(citations::Method::QTAIM, std::cout);
		//The maximum column is 16 wide: it carries rho at the attractor, ~3e8 at a uranium nucleus
		std::cout << "  basin  label               electrons" << (eli ? "" : "     charge") << "      volume         maximum        x          y          z\n";
		double total = 0.0;
		for (size_t b = 0; b < pop.size(); b++) {
			total += pop[b];
			std::cout << std::setw(7) << b + 1 << "  " << std::left << std::setw(18) << lab[b] << std::right << std::fixed
				<< std::setprecision(4) << std::setw(11) << pop[b];
			if (!eli) {
				//The label names the atom the maximum sits on; a non-nuclear attractor has no charge
				double Z = -1.0;
				for (int a = 0; a < l_w.get_ncen(); a++)
					if (lab[b] == l_w.get_atom_label(a) + std::to_string(a)) Z = l_w.get_atom_charge(a);
				if (Z < 0) std::cout << std::setw(11) << "-";
				else std::cout << std::setw(11) << Z - pop[b];
			}
			std::cout << std::setw(12) << vol[b] << std::setw(16) << res.second[b][3]
				<< std::setprecision(3) << std::setw(11) << res.second[b][0] << std::setw(11) << res.second[b][1] << std::setw(11) << res.second[b][2] << "\n";
		}
		std::cout << "  total in basins: " << std::setprecision(4) << total << "   outside every basin: " << outside << "\n";
		//A metal's core and its outer core shell (eli_core_radius), each summed; a core far from its
		//closed-shell count (eli_core_electrons, less what an ECP took) is flagged
		if (eli)
			for (int a = 0; a < l_w.get_ncen(); a++) {
				const std::string name = l_w.get_atom_label(a) + std::to_string(a);
				double core = 0.0, shell = 0.0;
				int n_core = 0, n_shell = 0;
				for (size_t b = 0; b < pop.size(); b++) {
					if (lab[b] == name + " core") { core += pop[b]; n_core++; }
					else if (lab[b] == name + " shell") { shell += pop[b]; n_shell++; }
				}
				const int expected = std::max(0, eli_core_electrons(static_cast<int>(l_w.get_atom_charge(a))) - static_cast<int>(l_w.get_atom_ECP_electrons(a)));
				if (n_shell) std::cout << "  " << name << ": core " << std::setprecision(4) << core << " e (closed shells " << expected << "), outer core shell " << shell << " e in "
					<< n_shell << " basins, together " << core + shell << " e\n";
				if (n_core && std::abs(core - expected) > std::max(1.0, 0.15 * expected))
					std::cout << "  WARNING: " << name << " core holds " << std::setprecision(2) << core << " e, its closed shells " << expected
					<< ": the core boundary sits in the wrong shell minimum\n";
			}
		if (want_aom) {
			std::cout.flush();
			section_log::heading("qtaim", "Localization and delocalization indices");
			report_delocalization(l_w, ovl, lab, std::cout);
			//-v: the overlap matrices to qtaim.aom, so an exact 2-RDM from outside (a CI vector) can be
			//contracted on the same basins: nmo nbasins, then MO index + occupation, then one packed triangle per basin
			if (opt.debug) {
				std::ofstream aom("qtaim.aom");
				aom << std::setprecision(17) << ovl.nmo << " " << ovl.S.size() << "\n";
				for (int i = 0; i < ovl.nmo; i++) aom << ovl.mo_index[i] << " " << l_w.get_MO_occ(ovl.mo_index[i]) << "\n";
				for (size_t b = 0; b < ovl.S.size(); b++) {
					aom << lab[b] << "\n";
					for (const double s : ovl.S[b]) aom << s << "\n";
				}
			}
		}
	};
	if (l_w.get_nmo() == 0 || opt.qtaim_only) {
		if (!opt.eli_only) report("QTAIM Analysis", qtaim_results, labels, false, stream_qtaim);
		if (!opt.qtaim_only) std::cout << "\nNo orbitals: the ELI-D basins are skipped." << std::endl;
		return;
	}

	std::pair<cubei, std::vector<d4>> eli_results;
	if (stream_eli) {
		//Gradient ascent on computeELIGrad from shells of seeds around every atom: no grid, and no Hessian (see below)
		eli_results.second = analytic_eli_maxima(l_w, opt.debug);
		eli_maxima_all = eli_results.second;
		T.lap("ELI-D maxima");
		std::cout << "ELI-D maxima from the analytic field: " << eli_results.second.size() << std::endl;
	}
	else {
	//ELI-D is undefined in the density tail, so the cube ends at the density isosurface.
	for (int x = 0; x < eli_cube.get_size(0); x++)
		for (int y = 0; y < eli_cube.get_size(1); y++)
			for (int z = 0; z < eli_cube.get_size(2); z++)
				if (rho.get_value(x, y, z) < basin_density_cutoff) eli_cube.set_value(x, y, z, 0.0);
	T.lap("ELI-D tail crop");
	//The only attractors found on the grid here (the QTAIM set is seeded from the nuclei), so for ELI-D the
	//resolution is an accuracy parameter, not a performance knob: coarser grids shift core populations and retype basins.
	if (grid_spacing > 0.05)
		std::cout << "WARNING: the ELI-D attractors are searched on the " << grid_spacing
			<< " A grid. Coarser than 0.05 A this basin set is less reliable: core populations have"
			" been seen to drift by 0.3 electrons and basins to appear or be retyped as lone pairs by"
			" 0.2 A. The QTAIM basins are unaffected." << std::endl;
	//The persistence merge absorbs a basin into its highest neighbour across their highest shared saddle. In a flat
	//valence shell every saddle is about as deep as the one to the core, so a loose threshold chains the shell into
	//the core basin. 3e-4 leaves only grid noise to it; unify_shell_basins merges shell maxima by distance (1.2 bohr)
	//and cannot chain into a core.
	//TODO: Na and Al cores come out too small; core_shell_radius does not reach the outer core shell.
	eli_results = topological_cube_analysis(&eli_cube, atoms, opt.debug, false, 0.0, 1e-10, radius, 3e-4);
	T.lap("ELI-D cube topology");
	}
	//There is no analytic ELI-D Hessian, only computeELIGrad, so maxima are not Newton-tested. Basin boundaries
	//come from walking every quadrature point up computeELIGrad, at any electron count.
	std::cout << "ELI-D basin boundaries from " << (stream_eli ? "the analytic field (the default)" : "the cube (-basin_cube)") << "." << std::endl;
	//The shells of a heavy atom's core structure ELI-D into several basins each; one core
	//basin per atom is what a bonding analysis wants, and what DGrid's ELIDcore gives
	const int core_merged = unify_core_basins(eli_results.first, eli_results.second, atoms, stream_eli ? &eli_core_map : nullptr);
	if (core_merged) std::cout << "Unified " << core_merged << " core-shell basins into their atoms' cores, " << eli_results.second.size() << " ELI-D basins remain." << std::endl;
	//Outside the cores the same sphere of maxima needs unify_shell_basins. NOS_ELI_SHELL_DIST / _TOL tune it;
	//distance 0 turns the merge off.
	double shell_dist = 1.2, shell_tol = 0.05;
	if (const char *e = tuning("NOS_ELI_SHELL_DIST")) { const double v = std::atof(e); if (v >= 0.0 && v < 10.0) shell_dist = v; }
	if (const char *e = tuning("NOS_ELI_SHELL_TOL")) { const double v = std::atof(e); if (v >= 0.0 && v < 1.0) shell_tol = v; }
	ivec eli_shell_map;
	const int shell_merged = unify_shell_basins(eli_results.first, eli_results.second, stream_eli ? &eli_shell_map : nullptr, shell_dist, shell_tol, &atoms, &l_w, nullptr, &eli_split);
	if (shell_merged) std::cout << "Unified " << shell_merged << " shattered shell basins, " << eli_results.second.size() << " ELI-D basins remain." << std::endl;
	//The integrator walks to one of eli_maxima_all and then reads eli_core_map, so the second merge
	//has to be composed into that map rather than replacing it
	if (!eli_shell_map.empty())
		for (size_t b = 1; b < eli_core_map.size(); b++) eli_core_map[b] = eli_shell_map[eli_core_map[b]];
	svec eli_labels = assign_labels_to_basins(eli_results.second, atoms, opt.debug, 1);
	number_repeats(eli_labels, false);
	//The orbital columns need one closed-shell determinant and the localisers' preconditions; checked here so a
	//column is skipped with its reason instead of the run exiting inside ibo.cpp
	std::string orbital_skip;
	if (opt.eli_contributions) {
		bool closed = true;
		for (int m = 0; m < wavy.get_nmo(); m++)
			if (wavy.get_MO_occ(m) != 0.0 && std::abs(wavy.get_MO_occ(m) - 2.0) > 1e-6) closed = false;
		const size_t bytes = (size_t)l_w.get_nmo() * (l_w.get_nmo() + 1) / 2 * eli_results.second.size() * omp_get_max_threads() * sizeof(double);
		if (wavy.get_is_unrestricted() || !closed) orbital_skip = "needs a closed-shell wavefunction (every occupation 0 or 2)";
		else if (bytes > (size_t)512 * 1024 * 1024) orbital_skip = "the overlap matrices would need " + std::to_string(bytes / (1024 * 1024)) + " MB";
		contrib_orbitals = orbital_skip.empty();
	}
	if (!opt.eli_only) report("QTAIM Analysis", qtaim_results, labels, false, stream_qtaim);
	report("ELI-D Analysis", eli_results, eli_labels, true, stream_eli);
	if (opt.eli_contributions && !eli_points.empty()) {
		std::vector<contrib_column> cols(4);
		cols[0].title = "QTAIM (Raub-Jansen)";
		cols[1].title = "PM-IBO";
		cols[2].title = "Boys-IBO";
		cols[3].title = "NLMO";
		const int nb = static_cast<int>(eli_pop.size());
		//Raub-Jansen: each ELI-D point's electrons go to the QTAIM atom its position climbs to
		if (opt.eli_only) cols[0].skipped = "-eli_only has no QTAIM basins";
		else {
			const int nq = static_cast<int>(qtaim_results.second.size());
			ivec at(eli_points.size(), 0);
			const std::function<void(const std::function<int(const d3 &)> &)> classify_points = [&](const std::function<int(const d3 &)> &classify) {
#pragma omp parallel for schedule(dynamic, 256)
				for (long long i = 0; i < (long long)eli_points.size(); i++) at[i] = classify(eli_points[i].p);
			};
			vec vol;
			double outside = 0.0;
			integrate_basins_on_atomic_grids(stream_qtaim ? nullptr : &rho, stream_qtaim ? nullptr : &(qtaim_results.first), qtaim_results.second, l_w, opt.accuracy, false, vol, outside, fill_cores ? &core_density : nullptr, fill_cores ? &core_gradient : nullptr, opt.basin_grid, fld, nullptr, nullptr, nullptr, nullptr, nullptr, &classify_points, true);
			cols[0].c.assign(nb, vec(nq + 1, 0.0));
			for (size_t i = 0; i < eli_points.size(); i++)
				cols[0].c[eli_points[i].basin - 1][at[i] > 0 ? at[i] - 1 : nq] += eli_points[i].q;
			cols[0].lab = labels;
			cols[0].lab.push_back("unassigned");
			T.lap("Raub-Jansen classification");
		}
		if (!contrib_orbitals)
			for (int c = 1; c < 4; c++) cols[c].skipped = orbital_skip;
		else {
			std::string ibo_skip;
			if (!origin_has_orca_pure_phases(wavy.get_origin())) ibo_skip = "IBOs need ORCA's pure-function phases (.gbw or ORCA .molden)";
			const std::shared_ptr<BasisSet> minao = BasisSetLibrary::get_basis_set("minao");
			for (int a = 0; a < wavy.get_ncen() && ibo_skip.empty(); a++)
				if (!minao->has_element(wavy.get_atom_charge(a)))
					ibo_skip = std::string("the MINAO reference basis has no ") + constants::atnr2letter(wavy.get_atom_charge(a));
			for (const bool boys : { false, true }) {
				contrib_column &col = cols[boys ? 2 : 1];
				if (!ibo_skip.empty()) { col.skipped = ibo_skip; continue; }
				const IBOResult ibo = intrinsic_bond_orbitals(wavy, boys);
				const int nk = static_cast<int>(ibo.mos.size());
				vec2 coef(nk, vec(wavy.get_nmo(), 0.0));
				for (int k = 0; k < nk; k++) {
					for (int i = 0; i < nk; i++) coef[k][ibo.mos[i]] = ibo.U(i, k);
					col.lab.push_back(zero_based(ibo_label(ibo, wavy, k)));
				}
				number_repeats(col.lab, true);
				col.c = orbital_contributions(wavy, eli_ovl, coef);
				T.lap(boys ? "Boys IBO contributions" : "PM IBO contributions");
			}
			if (wavy.get_nr_basis_set_loaded() != wavy.get_ncen() || wavy.get_d_f_switch())
				cols[3].skipped = "NBO needs the basis set on every atom and pure d/f functions";
			else {
				NboOptions nbo_opt;
				nbo_opt.threads = std::max(opt.threads, 0);
				nbo_opt.nrt_parent = opt.properties.nbo_lewis;
				//Under -fba the NBO thread holds the cores the search barriers would wait for
				if (opt.fba) nbo_opt.search_threads = 1;
				NboOrbitalSpace space;
				WFN copy(wavy);
				std::ostringstream nbo_log;
				const NboResults r = native_nbo(copy, nbo_opt, nbo_log, &space);
				std::vector<std::pair<int, std::string>> picked;
				vec2 coef;
				nbo_cube_wfn(wavy, r, space, "lewis", true, picked, nbo_log, &coef);
				for (const auto &p : picked) cols[3].lab.push_back(zero_based(p.second));
				number_repeats(cols[3].lab, true);
				cols[3].c = orbital_contributions(wavy, eli_ovl, coef);
				T.lap("NLMO contributions");
			}
		}
		const section_log::tee contrib_tee("eli", "ELI-D basin contributions");
		print_eli_contributions(eli_pop, eli_labels, cols, std::cout);
		if (fill_cores) std::cout << "  ECP core electrons count in the ELI-D core populations but in no column's (sum).\n";
	}
	if (opt.eli_only && stream_eli && opt.eli_basin_atoms.size() == 2) return;
	const section_log::tee spin_tee("eli", "Spin-resolved ELI-D");
	spin_eli_analysis(l_w, opt, atoms, shell_dist, shell_tol);
}

// ---------------------------------------------------------------------------
// QTAIM_ELI_mask
// ---------------------------------------------------------------------------

void QTAIM_ELI_mask(
	cube& rho,
	cube& eli,
	WFN& parent_wfn,
	const std::vector<atom>& atoms,
	const ivec& selected_indices,
	double background_value,
	const std::filesystem::path& output_path,
	bool debug,
	std::ostream& log
) {
	err_checkf(!atoms.empty(), "QTAIM_ELI_mask: no atoms available.", log);
	err_checkf(rho.get_size(0) == eli.get_size(0) &&
			   rho.get_size(1) == eli.get_size(1) &&
			   rho.get_size(2) == eli.get_size(2),
			   "QTAIM_ELI_mask: rho and eli grids have different dimensions.", log);

	// Validate atom indices
	for (int idx : selected_indices) {
		err_checkf(idx >= 0 && idx < (int)atoms.size(),
			"QTAIM_ELI_mask: atom index " + std::to_string(idx) +
			" is out of range (0.." + std::to_string((int)atoms.size() - 1) + ").", log);
	}

	log << "Running QTAIM topological analysis on density grid..." << std::endl;
	auto [basin_cube, maxima] = topological_cube_analysis(&rho, atoms, debug, false, basin_density_cutoff);
	svec labels = assign_labels_to_basins(maxima, atoms, debug, 0);

	// Map selected atom indices → set of 1-based basin IDs
	// Label format: element + atom_index, e.g. "C0", "H3"
	std::set<int> selected_basins;
	for (int b = 0; b < (int)labels.size(); b++) {
		const std::string& lbl = labels[b];
		// Extract trailing digit sequence (the atom index)
		auto it = std::find_if(lbl.rbegin(), lbl.rend(),
							   [](char c) { return !std::isdigit(static_cast<unsigned char>(c)); });
		if (it == lbl.rend()) continue; // entire string is digits — skip
		const std::string num_str(it.base(), lbl.end());
		if (num_str.empty()) continue;
		int atom_idx = std::stoi(num_str);
		if (std::find(selected_indices.begin(), selected_indices.end(), atom_idx)
				!= selected_indices.end()) {
			selected_basins.insert(b + 1); // basin IDs in cubei are 1-indexed
			if (debug)
				log << "  Basin " << (b + 1) << " (\"" << lbl << "\") selected.\n";
		}
	}

	if (selected_basins.empty()) {
		log << "Warning: no QTAIM basins matched the requested atom indices. "
			   "Output will contain only the background value.\n";
	}

	const int nx = eli.get_size(0);
	const int ny = eli.get_size(1);
	const int nz = eli.get_size(2);

	// Pass 1: find tight bounding box of selected voxels
	int xmin = nx, xmax = -1, ymin = ny, ymax = -1, zmin = nz, zmax = -1;
	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++) {
				if (selected_basins.count(basin_cube.get_value(x, y, z))) {
					xmin = std::min(xmin, x); xmax = std::max(xmax, x);
					ymin = std::min(ymin, y); ymax = std::max(ymax, y);
					zmin = std::min(zmin, z); zmax = std::max(zmax, z);
				}
			}

	if (xmax < 0) {
		// No selected voxels — write a 1×1×1 cube with background value
		xmin = xmax = 0; ymin = ymax = 0; zmin = zmax = 0;
	}

	log << "Bounding box of selected basins: ["
		<< xmin << ".." << xmax << "] x ["
		<< ymin << ".." << ymax << "] x ["
		<< zmin << ".." << zmax << "]\n";

	// Pass 2: build shrunk cube
	const std::array<int, 3> new_size = {xmax - xmin + 1, ymax - ymin + 1, zmax - zmin + 1};
	cube shrunk(new_size, (int)atoms.size(), true);
	shrunk.give_parent_wfn(parent_wfn);

	const auto new_origin = eli.get_pos(xmin, ymin, zmin);
	for (int i = 0; i < 3; i++) {
		shrunk.set_origin(i, new_origin[i]);
		for (int j = 0; j < 3; j++)
			shrunk.set_vector(i, j, eli.get_vector(i, j));
	}
	shrunk.calc_dv();

	for (int x = xmin; x <= xmax; x++)
		for (int y = ymin; y <= ymax; y++)
			for (int z = zmin; z <= zmax; z++) {
				double val = selected_basins.count(basin_cube.get_value(x, y, z))
							 ? eli.get_value(x, y, z)
							 : background_value;
				shrunk.set_value(x - xmin, y - ymin, z - zmin, val);
			}

	shrunk.set_comment1("QTAIM-masked ELI");
	shrunk.set_comment2("Selected atoms: " + [&]() {
		std::string s;
		for (int i = 0; i < (int)selected_indices.size(); i++) {
			if (i) s += ',';
			s += std::to_string(selected_indices[i]);
		}
		return s;
	}());
	shrunk.set_path(output_path);

	log << "Writing masked ELI cube to " << output_path.string() << " ..." << std::endl;
	err_checkf(shrunk.write_file(true), "QTAIM_ELI_mask: failed to write output cube.", log);
	log << "Done.\n";
}

// ---------------------------------------------------------------------------
// run_QTAIM_ELI_mask  (dispatch: cube-files mode vs WFN mode)
// ---------------------------------------------------------------------------

void run_QTAIM_ELI_mask(
	const std::filesystem::path& rho_or_wfn,
	const std::filesystem::path& eli_path,
	const ivec& selected_indices,
	double background_value,
	options& opt,
	std::ostream& log
) {
	if (!eli_path.empty()) {
		// ---- Cube-files mode ----
		log << "Reading density cube: " << rho_or_wfn.string() << std::endl;
		WFN dummy;
		cube rho(rho_or_wfn, true, dummy, log);
		log << "Reading ELI cube:     " << eli_path.string() << std::endl;
		cube eli(eli_path, true, dummy, log);
		rho.give_parent_wfn(dummy);
		eli.give_parent_wfn(dummy);

		const std::vector<atom> atoms = rho.get_parent_wfn_atoms();
		const std::filesystem::path out =
			eli_path.parent_path() / (eli_path.stem().string() + "_qtaim_masked.cube");

		QTAIM_ELI_mask(rho, eli, dummy, atoms, selected_indices, background_value, out, opt.debug, log);
	} else {
		// ---- WFN mode: compute rho and eli cubes first ----
		log << "Loading wavefunction: " << rho_or_wfn.string() << std::endl;
		WFN wavy(rho_or_wfn, opt.debug);
		wavy.delete_unoccupied_MOs();
		density_field field;
		const std::unique_ptr<Gaussian_Molecule> fit = fitted_source(wavy, opt, field, log);
		err_checkf(wavy.get_nmo() > 0, "ELI-D needs orbitals; -qtaim_eli cannot run on a SALTED density", log);

		properties_options prop_opt = opt.properties;
		readxyzMinMax_fromWFN(wavy, prop_opt);

		log << "Calculating density and ELI grid ("
			<< prop_opt.NbSteps[0] << " x " << prop_opt.NbSteps[1] << " x " << prop_opt.NbSteps[2]
			<< ") ..." << std::endl;

		cube rho(prop_opt.NbSteps, wavy.get_ncen(), true);
		cube eli(prop_opt.NbSteps, wavy.get_ncen(), true);
		rho.give_parent_wfn(wavy);
		eli.give_parent_wfn(wavy);

		const vec stepsizes{
			(prop_opt.MinMax[3] - prop_opt.MinMax[0]) / prop_opt.NbSteps[0],
			(prop_opt.MinMax[4] - prop_opt.MinMax[1]) / prop_opt.NbSteps[1],
			(prop_opt.MinMax[5] - prop_opt.MinMax[2]) / prop_opt.NbSteps[2]
		};
		for (int i = 0; i < 3; i++) {
			rho.set_origin(i, prop_opt.MinMax[i]);
			rho.set_vector(i, i, stepsizes[i]);
			eli.set_origin(i, prop_opt.MinMax[i]);
			eli.set_vector(i, i, stepsizes[i]);
		}
		rho.calc_dv();
		eli.calc_dv();

		Calc_RhoEli(rho, eli, wavy, prop_opt.radius, fit ? &field : nullptr);

		const std::vector<atom> atoms = wavy.get_atoms();
		const std::filesystem::path out =
			rho_or_wfn.parent_path() / "eli_qtaim_masked.cube";

		QTAIM_ELI_mask(rho, eli, wavy, atoms, selected_indices, background_value, out, opt.debug, log);
	}
}
