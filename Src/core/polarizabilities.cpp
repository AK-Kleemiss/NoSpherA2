#include "pch.h"
#include "wfn_class.h"
#include "properties.h"
#ifdef NOSPHERA2_USE_GPU
#include "aux_density_gpu.h"
#include "spherical_density_gpu.h"
#endif
#include "convenience.h"
#include "spherical_density.h"
#include "cube.h"
#include "constants.h"
#include "ibo.h"
#include "nbo.h"
#include "GridManager.h"
#include "isosurface.h"
#include "SALTED_predictor.h"
#include "crystal_energies.h"
#include "b2c.h"
#include "density_source.h"
#include "citations.h"
#include "eli_family.h"

namespace {

//Keith's bond graph: the pairs TFVC calls neighbours (no third atom inside the sphere on the pair, a Gabriel graph, so
//always connected), each with its bond point on the axis where make_chi found the density minimum (the midpoint where
//it found none) and the density there
struct BondGraph {
	std::vector<std::array<int, 2>> pair;
	std::vector<d3> point;
	vec rho;
};

BondGraph make_bond_graph(const WFN &w)
{
	const int n = w.get_ncen();
	BondGraph g;
	if (n < 2) return g;
	const vec chi = make_chi(w, 40, true, false);
	for (int a = 0; a < n; a++)
		for (int b = a + 1; b < n; b++) {
			const d3 A = w.get_atom_pos(a), B = w.get_atom_pos(b);
			const d3 mid{ (A[0] + B[0]) / 2, (A[1] + B[1]) / 2, (A[2] + B[2]) / 2 };
			const double half = array_length(A, B) / 2;
			bool gabriel = true;
			for (int c = 0; c < n && gabriel; c++)
				// on the sphere counts as inside: benzene para pairs have all four other C on it, and
				// an optimized geometry is off by ~1e-5 relative, so a strict < kept one of the three
				if (c != a && c != b && array_length(mid, w.get_atom_pos(c)) <= half * (1 + 1e-3)) gabriel = false;
			if (!gabriel) continue;
			const double x = chi.empty() ? 1.0 : chi[(size_t)a * n + b], t = x / (1 + x);
			g.pair.push_back({ a, b });
			g.point.push_back({ A[0] + t * (B[0] - A[0]), A[1] + t * (B[1] - A[1]), A[2] + t * (B[2] - A[2]) });
		}
	const int nb = (int)g.point.size();
	vec x(nb), y(nb), z(nb);
	for (int e = 0; e < nb; e++) { x[e] = g.point[e][0]; y[e] = g.point[e][1]; z[e] = g.point[e][2]; }
	g.rho.resize(nb);
	if (nb > 0) density_batch(w)(nb, x.data(), y.data(), z.data(), g.rho.data());
	return g;
}

//Keith's bond charges: q[e] is what pair[e][0] gets from bond e and -q[e] what pair[e][1] gets; an atom's bonds sum to
//its charge, and around every ring sum lambda q = 0 with PolaBer's lambda = 1/rho at the bond point (eqs 5-9). That is
//Kirchhoff on a network of conductance rho per bond: L phi = Q, q = rho (phi_a - phi_b), the smallest sum lambda q^2.
//Q is taken less its mean, which no bond can carry (an ion, or the integration error)
vec bond_charges(const BondGraph &g, const vec &Q)
{
	const int n = (int)Q.size(), m = n - 1; //the last atom is grounded
	vec q(g.pair.size(), 0.0);
	if (n < 2) return q;
	const double mean = std::accumulate(Q.begin(), Q.end(), 0.0) / n;
	vec L((size_t)m * m, 0.0), phi(n, 0.0), rhs(m);
	for (int a = 0; a < m; a++) rhs[a] = Q[a] - mean;
	for (size_t e = 0; e < g.pair.size(); e++) {
		const int a = g.pair[e][0], b = g.pair[e][1];
		const double c = std::max(g.rho[e], 1e-12);
		if (a < m) L[(size_t)a * m + a] += c;
		if (b < m) L[(size_t)b * m + b] += c;
		if (a < m && b < m) { L[(size_t)a * m + b] -= c; L[(size_t)b * m + a] -= c; }
	}
	err_checkf(solve_spd_system(L, (size_t)m, rhs) == 0, "bond charges: the bond graph is not connected", std::cout);
	std::copy(rhs.begin(), rhs.end(), phi.begin());
	for (size_t e = 0; e < g.pair.size(); e++)
		q[e] = std::max(g.rho[e], 1e-12) * (phi[g.pair[e][0]] - phi[g.pair[e][1]]);
	return q;
}

//An every-scheme grid keeps all atoms on helper grids, indexed by atom
vec3 &pol_grids(GridManager &gm, const ivec *&npts)
{
	GridData &gd = gm.getGridData();
	npts = gm.getNeedsHelper() ? &gd.helper_num_points_per_atom : &gd.num_points_per_atom;
	return gm.getNeedsHelper() ? gd.helper_grids : gd.atomic_grids;
}

//Puts the density of w on the grid built for another field: the frozen boundary keeps the zero-field weights
void fill_density(GridManager &gm, const WFN &w)
{
	const ivec *npts;
	vec3 &grids = pol_grids(gm, npts);
	const DensityBatch rho = density_batch(w);
	for (size_t a = 0; a < grids.size(); a++)
		rho((*npts)[a], grids[a][GridData::X].data(), grids[a][GridData::Y].data(), grids[a][GridData::Z].data(),
			grids[a][GridData::WFN_DENSITY].data());
}

//Per atom {mu_p x, y, z, mu_c x, y, z, Q} for one weight column. mu_p = -int w rho (r - R) and Q = Z - N as in PolaBer
//eq 3; mu_c = sum over bonds (R - R_b) q, plus the mean charge about the nuclear-charge centre, so the atoms add up to
//the molecular dipole about that centre
vec2 pol_dipoles(GridManager &gm, const WFN &w, const BondGraph &g, GridData::GridIndex col)
{
	const ivec *npts;
	const vec3 &grids = pol_grids(gm, npts);
	const int n = w.get_ncen();
	vec2 out(n, vec(7, 0.0));
	vec Q(n);
#pragma omp parallel for
	for (int a = 0; a < n; a++) {
		const vec2 &p = grids[a];
		const d3 R = w.get_atom_pos(a);
		double N = 0, mx = 0, my = 0, mz = 0;
		for (int i = 0; i < (*npts)[a]; i++) {
			const double f = p[GridData::WFN_DENSITY][i] * p[col][i];
			N += f;
			mx -= f * (p[GridData::X][i] - R[0]);
			my -= f * (p[GridData::Y][i] - R[1]);
			mz -= f * (p[GridData::Z][i] - R[2]);
		}
		Q[a] = w.get_atom_charge(a) - N - w.get_atom_ECP_electrons(a);
		out[a][0] = mx; out[a][1] = my; out[a][2] = mz; out[a][6] = Q[a];
	}
	const vec q = bond_charges(g, Q);
	const double mean = std::accumulate(Q.begin(), Q.end(), 0.0) / n;
	d3 Rc{ 0, 0, 0 };
	double Z = 0;
	for (int a = 0; a < n; a++) {
		Z += w.get_atom_charge(a);
		for (int k = 0; k < 3; k++) Rc[k] += w.get_atom_charge(a) * w.get_atom_pos(a)[k];
	}
	for (int k = 0; k < 3; k++) Rc[k] /= Z;
	for (int a = 0; a < n; a++)
		for (int k = 0; k < 3; k++) out[a][3 + k] = mean * (w.get_atom_pos(a)[k] - Rc[k]);
	for (size_t e = 0; e < g.pair.size(); e++)
		for (int s = 0; s < 2; s++) {
			const int a = g.pair[e][s];
			for (int k = 0; k < 3; k++)
				out[a][3 + k] += (w.get_atom_pos(a)[k] - g.point[e][k]) * (s == 0 ? q[e] : -q[e]);
		}
	return out;
}

}

//Distributed atomic polarizabilities from seven finite-field wavefunctions (0, +x, -x, +y, -y, +z, -z) for every
//real-space partition, each with the boundary frozen at zero field and moving with the field
void polarizabilities(options &opt, std::ostream &log2)
{
	using namespace std;
	std::vector<WFN> wavy;
	for (int i = 0; i < 7; i++)
	{
		wavy.emplace_back(e_origin::NOT_YET_DEFINED);
		wavy[i].read_known_wavefunction_format(opt.pol_wfns[i], log2, opt.debug);
		//compute_dens sizes its scratch by get_nmo(true); a gbw carries the virtuals
		wavy[i].delete_unoccupied_MOs();
		err_checkf(wavy[i].get_ncen() == wavy[0].get_ncen(), "The seven wavefunctions must hold the same atoms", log2);
	}
	const int n = wavy[0].get_ncen();
	const double eps = opt.efield;
	static const char *scheme[5] = { "Becke", "TFVC", "Hirshfeld", "MBIS", "EMBIS" };
	static const GridData::GridIndex col[5] = { GridData::BECKE_WEIGHT, GridData::TFVC_WEIGHT, GridData::HIRSH_WEIGHT,
		GridData::MBIS_WEIGHT, GridData::EMBIS_WEIGHT };
	static const char *mode[2] = { "frozen", "moving" };
	GridConfiguration cfg;
	cfg.accuracy = opt.accuracy;
	cfg.all_charges = true;
	ivec atoms(n);
	iota(atoms.begin(), atoms.end(), 0);
	std::ostream quiet(nullptr);
	std::ostream &glog = opt.debug ? log2 : quiet;

	//mu[mode][scheme][field] = [atom]{mu_p, mu_c, Q}
	vector<vector<vec3>> mu(2, vector<vec3>(5, vec3(7)));
	GridManager gm0(cfg);
	gm0.setup3DGridsForMolecule(wavy[0], atoms, {}, cell(), false, glog);
	const BondGraph g0 = make_bond_graph(wavy[0]);
	for (int s = 0; s < 5; s++)
		mu[0][s][0] = mu[1][s][0] = pol_dipoles(gm0, wavy[0], g0, col[s]);
	for (int f = 1; f < 7; f++) {
		log2 << "Partitioning field " << f << " of 6..." << endl;
		GridManager gm(cfg);
		gm.setup3DGridsForMolecule(wavy[f], atoms, {}, cell(), false, glog);
		const BondGraph g = make_bond_graph(wavy[f]);
		err_checkf(g.pair == g0.pair, "The bond graph changed with the field", log2);
		fill_density(gm0, wavy[f]);
		for (int s = 0; s < 5; s++) {
			mu[0][s][f] = pol_dipoles(gm0, wavy[f], g0, col[s]);
			mu[1][s][f] = pol_dipoles(gm, wavy[f], g, col[s]);
		}
	}

	//Programs differ in the sign of the field (Gaussian's lowers the dipole along +E); the molecular trace is positive
	double trace = 0;
	for (int a = 0; a < n; a++)
		for (int i = 0; i < 3; i++)
			trace += mu[0][0][1 + 2 * i][a][i] + mu[0][0][1 + 2 * i][a][3 + i] - mu[0][0][2 + 2 * i][a][i] - mu[0][0][2 + 2 * i][a][3 + i];
	const double sgn = trace < 0 ? -1.0 : 1.0;
	//part 0 = alpha_p, 1 = alpha_c, -1 = both; column j is the field direction
	auto alpha = [&](int m, int s, int a, int part, int i, int j) {
		double v = 0;
		for (int k = (part < 0 ? 0 : part); k <= (part < 0 ? 1 : part); k++)
			v += mu[m][s][1 + 2 * j][a][3 * k + i] - mu[m][s][2 + 2 * j][a][3 * k + i];
		return sgn * v / (2 * eps);
	};
	auto sym = [&](int m, int s, int a, int part, int i, int j) { return (alpha(m, s, a, part, i, j) + alpha(m, s, a, part, j, i)) / 2; };
	static const int ij[6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2} };
	auto label = [&](int a) { return constants::atnr2letter(wavy[0].get_atom_charge(a)) + to_string(a + 1); };

	log2 << "\nAtomic polarizabilities in a.u. (bohr^3) from fields of " << eps << " a.u." << (sgn < 0 ? " (field sign of the input reversed)" : "") << "\n"
		<< "alpha = alpha_p + alpha_c: alpha_p from the atom's own dipole, alpha_c from Keith's charge translation through the\n"
		<< "bonds (bond charges weighted 1/rho at the bond point, as in PolaBer); tensors symmetrized.\n"
		<< "frozen: partition weights and bond points of the zero-field density; moving: rebuilt for every field.\n"
		<< "\nBond graph (zero field): atoms, distance of the bond point from the first atom (bohr), rho there (a.u.)\n";
	for (size_t e = 0; e < g0.pair.size(); e++)
		log2 << setw(6) << label(g0.pair[e][0]) << setw(6) << label(g0.pair[e][1]) << fixed << setprecision(4)
			<< setw(10) << array_length(wavy[0].get_atom_pos(g0.pair[e][0]), g0.point[e]) << setw(10) << g0.rho[e] << "\n";

	const filesystem::path dat = wavy[0].get_path().parent_path() / "polarizabilities.dat";
	ofstream table(dat);
	table << "# NoSpherA2 atomic polarizabilities, a.u., symmetrized; a = total, p = polarization, c = charge translation\n"
		<< "# boundary scheme atom label Q a_xx a_yy a_zz a_xy a_xz a_yz p_xx p_yy p_zz p_xy p_xz p_yz c_xx c_yy c_zz c_xy c_xz c_yz\n";
	for (int m = 0; m < 2; m++)
		for (int s = 0; s < 5; s++) {
			log2 << "\n" << scheme[s] << ", " << mode[m] << " boundary\n"
				<< "  atom       Q        xx       yy       zz       xy       xz       yz      iso    iso_p    iso_c\n";
			for (int a = 0; a <= n; a++) {
				const bool sum = a == n;
				double t[6] = {}, iso[3] = {};
				for (int b = sum ? 0 : a; b < (sum ? n : a + 1); b++) {
					for (int c = 0; c < 6; c++) t[c] += sym(m, s, b, -1, ij[c][0], ij[c][1]);
					for (int k = 0; k < 3; k++) iso[0] += sym(m, s, b, -1, k, k) / 3, iso[1] += sym(m, s, b, 0, k, k) / 3, iso[2] += sym(m, s, b, 1, k, k) / 3;
				}
				log2 << setw(6) << (sum ? string("sum") : label(a)) << fixed << setprecision(4)
					<< setw(9) << (sum ? 0.0 : mu[m][s][0][a][6]);
				for (int c = 0; c < 6; c++) log2 << setw(9) << t[c];
				log2 << setw(9) << iso[0] << setw(9) << iso[1] << setw(9) << iso[2] << "\n";
				if (sum) continue;
				table << mode[m] << " " << scheme[s] << " " << a + 1 << " " << label(a) << scientific << setprecision(8) << " " << mu[m][s][0][a][6];
				for (int part = -1; part < 2; part++)
					for (int c = 0; c < 6; c++) table << " " << sym(m, s, a, part, ij[c][0], ij[c][1]);
				table << "\n";
			}
		}
	log2 << "\nTable for plotting: " << dat.string() << endl;
	std::cout << "\n\nProperties calculation done!" << std::endl;
}
