#include "pch.h"
#include "tuning.h"
#include "wfn_class.h"
#include "convenience.h"
#include "constants.h"
#include "b2c.h"

double core_shell_radius(const int Z)
{
	if (Z <= 2) return 0.0;
	if (Z <= 10) return 0.25;
	if (Z <= 36) return 1.0;
	if (Z <= 54) return 1.4;
	return 1.8;
}

//A d-block metal's outer core shell, (n-1)s,p with the d electrons, stands apart from its inner core in ELI-D and
//DGrid's ELIDcore keeps it so (Ca's 3s3p likewise); core_shell_radius reaches past it, these radii sit in the gap.
//Zn, Cd and Hg keep core_shell_radius: a closed d10 core shows no separate shell.
double eli_core_radius(const int Z)
{
	if (Z >= 19 && Z <= 29) return 0.5;  // K-Cu
	if (Z >= 39 && Z <= 47) return 0.6;  // Y-Ag
	if (Z == 57 || (Z >= 72 && Z <= 79)) return 0.8;  // La, Hf-Au
	return core_shell_radius(Z);
}

//The ELI-D N/O minimum of a free 5d atom holds less than 60 electrons, whatever the charge, because part of the 4f
//tail lies beyond it; hence a relative tolerance
int eli_core_electrons(const int Z)
{
	if (Z >= 19 && Z <= 29) return 10;
	if (Z >= 39 && Z <= 47) return 28;
	if (Z == 57) return 46;
	if (Z >= 72 && Z <= 79) return 60;
	//the largest closed core below Z: a noble gas, or one with a filled (n-1)d and 4f beneath it
	int core = 0;
	for (const int c : {2, 10, 18, 28, 36, 46, 54, 78, 86})
		if (c < Z) core = c;
	return core;
}

static bool in_outer_core_shell(const d4 &m, const std::vector<atom> &atoms)
{
	for (const atom &a : atoms) {
		const int Z = a.get_charge();
		const double inner = eli_core_radius(Z), outer = core_shell_radius(Z);
		if (inner >= outer) continue;
		const d3 p = a.get_pos();
		const double d2 = std::pow(m[0] - p[0], 2) + std::pow(m[1] - p[1], 2) + std::pow(m[2] - p[2], 2);
		if (d2 >= inner * inner && d2 < outer * outer) return true;
	}
	return false;
}

//Both merges below only decide which maxima belong together; this keeps the highest of each group, renumbers the cube
//and reports where each old basin went. keeper[b] is the group's lowest member, or -1 for a maximum standing alone.
static int collapse_maxima_groups(cubei &basin_cube, std::vector<d4> &maxima, const ivec &keeper, ivec *basin_map)
{
	const int nb = static_cast<int>(maxima.size());
	ivec target(nb + 1);
	for (int b = 0; b <= nb; b++) target[b] = b;
	for (int b = 0; b < nb; b++) {
		if (keeper[b] < 0 || keeper[b] == b) continue;
		target[b + 1] = target[keeper[b] + 1];
		const int keep = target[b + 1] - 1;
		//Symmetry-equivalent maxima tie up to rounding; the tie goes to the lexicographically
		//larger position so the surviving maximum is the same on every platform
		const bool tie = std::abs(maxima[b][3] - maxima[keep][3]) < 1e-8 * std::abs(maxima[keep][3]);
		if (tie ? maxima[b] > maxima[keep] : maxima[b][3] > maxima[keep][3]) std::swap(maxima[b], maxima[keep]);
	}
	ivec renumber(nb + 1, 0);
	std::vector<d4> kept;
	for (int b = 1; b <= nb; b++)
		if (target[b] == b) { kept.push_back(maxima[b - 1]); renumber[b] = static_cast<int>(kept.size()); }
	for (int x = 0; x < basin_cube.get_size(0); x++)
		for (int y = 0; y < basin_cube.get_size(1); y++)
			for (int z = 0; z < basin_cube.get_size(2); z++) {
				const int b = basin_cube.get_value(x, y, z);
				if (b > 0) basin_cube.set_value(x, y, z, renumber[target[b]]);
			}
	//The swap above only ever exchanges two maxima of one group, and every member of that group
	//shares one target, so renumber[target[b]] is the same number before and after it
	if (basin_map) {
		basin_map->assign(nb + 1, 0);
		for (int b = 1; b <= nb; b++) (*basin_map)[b] = renumber[target[b]];
	}
	const int merged = nb - static_cast<int>(kept.size());
	maxima.swap(kept);
	return merged;
}

int unify_core_basins(cubei &basin_cube, std::vector<d4> &maxima, const std::vector<atom> &atoms, ivec *basin_map)
{
	const int nb = static_cast<int>(maxima.size());
	ivec owner(nb, -1);
	for (int b = 0; b < nb; b++)
		for (size_t a = 0; a < atoms.size(); a++) {
			const d3 ap = atoms[a].get_pos();
			const double r = eli_core_radius(atoms[a].get_charge());
			if (std::pow(maxima[b][0] - ap[0], 2) + std::pow(maxima[b][1] - ap[1], 2) + std::pow(maxima[b][2] - ap[2], 2) < r * r) { owner[b] = static_cast<int>(a); break; }
		}
	//The atom's core keeps the highest of its maxima; the merged ones are dropped
	ivec keeper(nb, -1);
	for (int b = 0; b < nb; b++) {
		if (owner[b] < 0) continue;
		for (int c = 0; c < b; c++)
			if (owner[c] == owner[b]) { keeper[b] = (keeper[c] < 0 ? c : keeper[c]); break; }
	}
	return collapse_maxima_groups(basin_cube, maxima, keeper, basin_map);
}

int unify_shell_basins(cubei &basin_cube, std::vector<d4> &maxima, ivec *basin_map, double max_dist, double rel_tol, const std::vector<atom> *atoms, const WFN *wavy, const eli_spin_field *eli, std::vector<shell_split> *split)
{
	const int nb = static_cast<int>(maxima.size());
	if (split) split->clear();
	if (nb < 2 || max_dist <= 0.0) return 0;
	const double d2 = max_dist * max_dist;
	std::vector<char> shell(nb, 0);
	if (atoms)
		for (int b = 0; b < nb; b++) shell[b] = in_outer_core_shell(maxima[b], *atoms);
	ivec root(nb);
	for (int b = 0; b < nb; b++) root[b] = b;
	auto find = [&root](int b) { while (root[b] != b) b = root[b] = root[root[b]]; return b; };
	for (int b = 1; b < nb; b++)
		for (int c = 0; c < b; c++) {
			//Two outer-core-shell maxima never fold by distance: the fold chained genuine shell maxima through low-contrast
			//shoulders, transitively (5.4-5.7 e merged basins, Tauchert set). They join by persistence below; a shell-to-
			//valence pair folds at half the length
			if (shell[b] && shell[c]) continue;
			if (std::pow(maxima[b][0] - maxima[c][0], 2) + std::pow(maxima[b][1] - maxima[c][1], 2) + std::pow(maxima[b][2] - maxima[c][2], 2) > (shell[b] || shell[c] ? d2 / 4.0 : d2)) continue;
			const double hi = std::max(std::abs(maxima[b][3]), std::abs(maxima[c][3]));
			if (hi > 0.0 && std::abs(maxima[b][3] - maxima[c][3]) > rel_tol * hi) continue;
			const int rb = find(b), rc = find(c);
			//union by lower index, so a group's root is always its lowest member
			if (rb != rc) root[std::max(rb, rc)] = std::min(rb, rc);
		}
	//Saddle criterion (NOS_ELI_SHELL_SADDLE = phi, default 0.99, 0 turns it off): two maxima nearest the same atom, at radii
	//within 0.3 bohr and less than 120 deg apart about it, join when the radial maximum of the field along the arc between
	//them never drops below phi times the lower of the two
	double phi = 0.99;
	if (const char *e = tuning("NOS_ELI_SHELL_SADDLE")) { const double v = std::atof(e); if (v >= 0.0 && v < 1.0) phi = v; }
	//Shell persistence (NOS_ELI_SHELL_PERSIST = eps, default 5e-4, 0 turns it off): two shell maxima of the same atom join,
	//highest arc first, when the arc between them dips less than eps below the lower of their two groups' maxima (elder
	//rule). Contrast alone separates the two kinds: spurious shoulders dip <= 2.9e-4, genuine shell saddles >= 7.5e-4
	//(as little as 0.05 % relative), so neither a distance nor a relative cut works (eli_merged validation, 24 Oct, 8 systems).
	//Measured safe window 5e-4 to 6e-4 over 9 TM references + 20 Tauchert systems: [PtCl4]2- sigma_h dimple dies at 4e-4,
	//the NiZn PMe3 shell saddle at 7e-4 (eli_tm_reference + eli_merged, 24 Oct)
	double eps = 5e-4;
	if (const char *e = tuning("NOS_ELI_SHELL_PERSIST")) {
		const double v = std::atof(e);
		if (v >= 0.0) eps = v;
		else std::cout << "NOS_ELI_SHELL_PERSIST=" << e << " ignored (must be >= 0), eps stays " << eps << std::endl;
	}
	if (wavy && (phi > 0.0 || eps > 0.0)) {
		auto field = [&](const d3 &p) { double f; d3 g; if (eli) (*eli)(p, f, g, nullptr); else wavy->computeELIGrad(p, f, g); return f; };
		ivec na(nb, 0);
		for (int b = 0; b < nb; b++) {
			const d3 pb{ maxima[b][0], maxima[b][1], maxima[b][2] };
			for (int a = 1; a < wavy->get_ncen(); a++)
				if (array_length(pb, wavy->get_atom_pos(a)) < array_length(pb, wavy->get_atom_pos(na[b]))) na[b] = a;
			//a spin field's surface attractor sits on the rho cutoff, where the ray reaches past the domain and never dips
			if (wavy->compute_dens(pb) < 2.0 * basin_density_cutoff) na[b] = -1;
		}
		//The largest value along the ray from o through w within 0.25 bohr of r: 11 samples and a parabola
		auto radial_max = [&](const d3 &o, const d3 &w, const double r) {
			double s[11];
			int j = 0;
			for (int k = 0; k < 11; k++) {
				const double rr = r + (k - 5) * 0.05;
				s[k] = field(d3{ o[0] + rr * w[0], o[1] + rr * w[1], o[2] + rr * w[2] });
				if (s[k] > s[j]) j = k;
			}
			if (j == 0 || j == 10) return s[j];
			const double den = s[j + 1] - 2.0 * s[j] + s[j - 1];
			return den < 0.0 ? s[j] - std::pow(s[j + 1] - s[j - 1], 2) / (8.0 * den) : s[j];
		};
		std::vector<std::pair<int, int>> pairs;
		for (int b = 1; b < nb; b++)
			for (int c = 0; c < b; c++) {
				if (na[b] < 0 || na[b] != na[c] || find(b) == find(c)) continue;
				if (shell[b] && shell[c] ? !(eps > 0.0) : !(phi > 0.0)) continue;
				const d3 o = wavy->get_atom_pos(na[b]);
				const d3 u{ maxima[b][0] - o[0], maxima[b][1] - o[1], maxima[b][2] - o[2] }, v{ maxima[c][0] - o[0], maxima[c][1] - o[1], maxima[c][2] - o[2] };
				const double ru = array_length(u), rv = array_length(v);
				if (!(ru > 0.2) || !(rv > 0.2) || std::abs(ru - rv) > 0.3) continue;
				if ((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (ru * rv) < -0.5) continue;
				pairs.emplace_back(b, c);
			}
		//the lowest radial maximum along each arc, and the lower endpoint
		vec arc(pairs.size(), 0.0), end(pairs.size(), 0.0);
		//direction about the atom of each arc's lowest interior sample: the saddle a split shoulder exits through
		std::vector<d3> sad(pairs.size());
		//ponytail: every candidate pair gets its own arc, 16 x 11 field calls at worst; a merge tree over a sphere of
		//directions if a shell ever has hundreds of maxima
#pragma omp parallel for schedule(dynamic)
		for (long long q = 0; q < static_cast<long long>(pairs.size()); q++) {
			const d4 &mb = maxima[pairs[q].first], &mc = maxima[pairs[q].second];
			const d3 o = wavy->get_atom_pos(na[pairs[q].first]);
			const d3 u{ mb[0] - o[0], mb[1] - o[1], mb[2] - o[2] }, v{ mc[0] - o[0], mc[1] - o[1], mc[2] - o[2] };
			const double ru = array_length(u), rv = array_length(v);
			const d3 uh{ u[0] / ru, u[1] / ru, u[2] / ru }, vh{ v[0] / rv, v[1] / rv, v[2] / rv };
			end[q] = std::min(radial_max(o, uh, ru), radial_max(o, vh, rv));
			double lo = end[q], low = std::numeric_limits<double>::max();
			for (int k = 1; k <= 16; k++) {
				const double t = k / 17.0;
				d3 w{ (1.0 - t) * uh[0] + t * vh[0], (1.0 - t) * uh[1] + t * vh[1], (1.0 - t) * uh[2] + t * vh[2] };
				const double wl = array_length(w);
				if (!(wl > 1e-6)) { lo = -1.0; break; }
				w = d3{ w[0] / wl, w[1] / wl, w[2] / wl };
				const double s = radial_max(o, w, (1.0 - t) * ru + t * rv);
				if (s < low) { low = s; sad[q] = w; }
				lo = std::min(lo, s);
			}
			arc[q] = lo;
		}
		const bool log = tuning("NOS_ELI_SHELL_SADDLE_LOG") != nullptr;
		int joined = 0, persisted = 0;
		ivec order;
		for (size_t q = 0; q < pairs.size(); q++) {
			if (shell[pairs[q].first] && shell[pairs[q].second]) { order.push_back(static_cast<int>(q)); continue; }
			if (arc[q] < phi * end[q]) continue;
			joined++;
			const d4 &mb = maxima[pairs[q].first], &mc = maxima[pairs[q].second];
			if (log) std::cout << "  saddle join about atom " << na[pairs[q].first] << ": " << mb[0] << " " << mb[1] << " " << mb[2] << " (" << mb[3] << ") + " << mc[0] << " " << mc[1] << " " << mc[2] << " (" << mc[3] << ")" << std::endl;
			const int rb = find(pairs[q].first), rc = find(pairs[q].second);
			if (rb != rc) root[std::max(rb, rc)] = std::min(rb, rc);
		}
		//group maximum = the highest logged member, as collapse_maxima_groups keeps it
		vec ghi(nb);
		for (int b = 0; b < nb; b++) ghi[b] = maxima[b][3];
		for (int b = 0; b < nb; b++) { const int r = find(b); ghi[r] = std::max(ghi[r], ghi[b]); }
		std::stable_sort(order.begin(), order.end(), [&](const int x, const int y) { return arc[x] > arc[y]; });
		for (const int q : order) {
			const int rb = find(pairs[q].first), rc = find(pairs[q].second);
			if (rb == rc || arc[q] < std::min(ghi[rb], ghi[rc]) - eps) continue;
			persisted++;
			//the death to 3 digits: the log's 4 decimals alone carry +-1e-4 against a 3e-4 wide eps window
			if (log) {
				char death[16];
				std::snprintf(death, sizeof death, "%.2e", std::min(ghi[rb], ghi[rc]) - arc[q]);
				std::cout << "  shell persistence join about atom " << na[pairs[q].first] << ": arc " << arc[q] << " vs " << std::min(ghi[rb], ghi[rc])
						  << ", death " << death << std::endl;
			}
			root[std::max(rb, rc)] = std::min(rb, rc);
			ghi[std::min(rb, rc)] = std::max(ghi[rb], ghi[rc]);
		}
		//The elder rule hands a dying shoulder whole to one of several equally reachable groups (Ni(CO)4: all four ligand-
		//direction shoulders into -x, 3.7 e against 2.3). Its exits are the surviving group maxima it reaches within eps;
		//over two or more groups the integrator splits it by saddle. Exits only to group tops: an arc to another shoulder
		//runs over an axis maximum and would name the wrong group
		int splits = 0;
		if (split && eps > 0.0)
			for (int b = 0; b < nb; b++) {
				if (!shell[b] || na[b] < 0 || !(maxima[b][3] < ghi[find(b)])) continue;
				shell_split s{ d3{ maxima[b][0], maxima[b][1], maxima[b][2] }, wavy->get_atom_pos(na[b]), {} };
				ivec groups;
				for (size_t q = 0; q < pairs.size(); q++) {
					const int c = pairs[q].first == b ? pairs[q].second : pairs[q].second == b ? pairs[q].first : -1;
					if (c < 0 || !shell[c] || maxima[c][3] < ghi[find(c)] || arc[q] < maxima[b][3] - eps) continue;
					s.exits.emplace_back(sad[q], d3{ maxima[c][0], maxima[c][1], maxima[c][2] });
					if (log) std::cout << "  shoulder " << b << " exit to " << c << " saddle dir " << sad[q][0] << " " << sad[q][1] << " " << sad[q][2] << ", arc " << arc[q] << std::endl;
					if (std::find(groups.begin(), groups.end(), find(c)) == groups.end()) groups.push_back(find(c));
				}
				if (groups.size() < 2) continue;
				splits++;
				split->push_back(std::move(s));
			}
		if (splits) std::cout << "Split " << splits << " shell shoulders between their exit basins by saddle direction." << std::endl;
		std::cout << "Shell saddle merge (phi " << phi << ", shell persistence " << eps << "): " << pairs.size() << " same-atom pairs tested, "
				  << joined << " saddle-joined, " << persisted << " shell pairs joined." << std::endl;
	}
	ivec keeper(nb);
	for (int b = 0; b < nb; b++) { const int r = find(b); keeper[b] = (r == b ? -1 : r); }
	return collapse_maxima_groups(basin_cube, maxima, keeper, basin_map);
}

int unify_boundary_basins(std::vector<d4> &maxima, const WFN &wavy, const eli_spin_field &eli, ivec *basin_map)
{
	const int nb = static_cast<int>(maxima.size());
	std::vector<char> edge(nb, 0);
	for (int b = 0; b < nb; b++) edge[b] = wavy.compute_dens(d3{ maxima[b][0], maxima[b][1], maxima[b][2] }) < 1.05 * basin_density_cutoff;
	//No barrier on the arc from b to c about the atom nearest b, at their radii: the field stays within 1e-3 of the two
	//maxima at four points between. A degenerate sphere is flat there; equivalent lone pairs or bonds dip between them
	auto flat = [&](const d4 &b, const d4 &c) {
		const d3 pb{ b[0], b[1], b[2] }, pc{ c[0], c[1], c[2] };
		int na = 0;
		for (int a = 1; a < wavy.get_ncen(); a++)
			if (array_length(pb, wavy.get_atom_pos(a)) < array_length(pb, wavy.get_atom_pos(na))) na = a;
		const d3 o = wavy.get_atom_pos(na);
		const d3 u{ pb[0] - o[0], pb[1] - o[1], pb[2] - o[2] }, v{ pc[0] - o[0], pc[1] - o[1], pc[2] - o[2] };
		const double ru = array_length(u), rv = array_length(v);
		if (!(ru > 1e-6) || !(rv > 1e-6)) return false;
		const double lo = (1.0 - 1e-3) * std::min(b[3], c[3]);
		for (int k = 1; k <= 4; k++) {
			const double t = k / 5.0, r = (1.0 - t) * ru + t * rv;
			//normalised linear interpolation of the two directions: the arc is short, 3 bohr at most
			d3 w{ (1.0 - t) * u[0] / ru + t * v[0] / rv, (1.0 - t) * u[1] / ru + t * v[1] / rv, (1.0 - t) * u[2] / ru + t * v[2] / rv };
			const double wl = array_length(w);
			if (!(wl > 1e-6)) return false;
			const d3 p{ o[0] + r * w[0] / wl, o[1] + r * w[1] / wl, o[2] + r * w[2] / wl };
			double f;
			d3 g;
			eli(p, f, g, nullptr);
			if (!(f >= lo)) return false;
		}
		return true;
	};
	ivec root(nb);
	std::iota(root.begin(), root.end(), 0);
	auto find = [&root](int b) { while (root[b] != b) b = root[b] = root[root[b]]; return b; };
	//ponytail: equal value to 1e-3 and 3 bohr single linkage; surface maxima link without the arc test (equivalent
	//open-shell sites on the isosurface closer than 3 bohr would merge), interior ones only over a flat arc
	for (int b = 1; b < nb; b++)
		for (int c = 0; c < b; c++) {
			if (std::pow(maxima[b][0] - maxima[c][0], 2) + std::pow(maxima[b][1] - maxima[c][1], 2) + std::pow(maxima[b][2] - maxima[c][2], 2) > 9.0) continue;
			if (std::abs(maxima[b][3] - maxima[c][3]) > 1e-3 * std::max(maxima[b][3], maxima[c][3])) continue;
			const int rb = find(b), rc = find(c);
			if (rb == rc) continue;
			if ((edge[b] && edge[c]) || flat(maxima[b], maxima[c])) root[std::max(rb, rc)] = std::min(rb, rc);
		}
	ivec keeper(nb);
	for (int b = 0; b < nb; b++) { const int r = find(b); keeper[b] = (r == b ? -1 : r); }
	cubei none;
	return collapse_maxima_groups(none, maxima, keeper, basin_map);
}
