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

std::vector<Thakkar> make_thakkar_interpolators()
{
	// atom_models[atom.charge - 1] is the lookup convention (see
	// promolecular_fragment_densities_at), so index a must hold atomic
	// number a+1, not a=0 (which would build a bogus atomic-number-0 model
	// and read Thakkar_ns/np/nd/nf out of bounds for every real element).
	std::vector<Thakkar> atom_models;
	atom_models.reserve(92);
	for (int a = 0; a < 92; a++)
	{
		atom_models.emplace_back(a + 1);
		atom_models[a].make_interpolator(1.005 * 1.005 * 1.005, 1E-7);
	}
	return atom_models;
}

void print_time(_time_point &start, _time_point &end, std::ostream &file) {
	if (constants::hide_timings)
		return;
	//fixed/setprecision(0) must not leak into the caller's later prints
	const std::ios_base::fmtflags flags = file.flags();
	const std::streamsize prec = file.precision();
	if (get_sec(start, end) < 60)
		file << "Time to calculate Values: " << std::fixed << std::setprecision(0) << get_sec(start, end) << " s" << std::endl;
	else if (get_sec(start, end) < 3600)
		file << "Time to calculate Values: " << std::fixed << std::setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << std::endl;
	else
		file << "Time to calculate Values: " << std::fixed << std::setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << std::endl;
	file.flags(flags);
	file.precision(prec);
}

namespace {

bool is_within_radius(const d3 &pos, const std::vector<atom> &atoms, double radius_bohr)
{
	for (const atom &entry : atoms)
		if (array_length(pos, entry.get_pos()) < radius_bohr)
			return true;
	return false;
}

double sanitize_finite(double value)
{
	if (std::isnan(value) || std::isinf(value))
		return 0.0;
	return value;
}

// The in-radius points in one batch: the ESP walks the whole primitive-pair table per point, so the set is worth
// handing to a device at once. Points outside stay zero as in evaluate_cube_in_radius; wrap sums periodic images per cell.
template <typename BatchFn>
void evaluate_cube_in_radius_batched(
	cube &target,
	bool wrap,
	const std::vector<atom> &atoms,
	double radius_bohr,
	BatchFn &&batch)
{
	std::vector<d3> points;
	std::vector<i3> cells;
	target.evaluate_on_grid(
		[&](const d3 &pos, const i3 &, const i3 &mapped) {
			if (!is_within_radius(pos, atoms, radius_bohr))
				return 0.0;
#pragma omp critical(cube_gather)
			{
				points.push_back(pos);
				cells.push_back(mapped);
			}
			return 0.0;
		},
		wrap);
	vec values(points.size());
	batch(points, values.data());
	for (size_t i = 0; i < points.size(); i++)
	{
		const i3 &c = cells[i];
		target.set_value(c[0], c[1], c[2], target.get_value(c[0], c[1], c[2]) + values[i]);
	}
}

template <typename EvalFn>
void evaluate_cube_in_radius(
	cube &target,
	bool wrap,
	const std::vector<atom> &atoms,
	double radius_bohr,
	EvalFn &&evaluate_inside)
{
	target.evaluate_on_grid(
		[&](const d3 &pos) {
			if (!is_within_radius(pos, atoms, radius_bohr))
				return 0.0;
			return evaluate_inside(pos);
		},
		wrap);
}

template <typename EvalFn>
void evaluate_cube_in_radius_mapped(
	cube &target,
	bool wrap,
	const std::vector<atom> &atoms,
	double radius_bohr,
	EvalFn &&evaluate_inside_mapped)
{
	target.evaluate_on_grid(
		[&](const d3 &pos, const i3 &raw_idx, const i3 &mapped_idx) {
			if (!is_within_radius(pos, atoms, radius_bohr))
				return 0.0;
			return evaluate_inside_mapped(pos, raw_idx, mapped_idx);
		},
		wrap);
}

template <typename EvalFn>
void evaluate_cube_near_atom_mapped(
	cube &target,
	bool wrap,
	const d3 &atom_pos,
	double radius_bohr,
	EvalFn &&evaluate_inside_mapped)
{
	target.evaluate_on_grid(
		[&](const d3 &pos, const i3 &raw_idx, const i3 &mapped_idx) {
			if (array_length(pos, atom_pos) >= radius_bohr)
				return 0.0;
			return evaluate_inside_mapped(pos, raw_idx, mapped_idx);
		},
		wrap);
}

//Per B^3 block of cube points (and `extra` more layers on the high side), the atoms (in input order) whose centre lies
//within reach(b, a) of some point of block b
template <class Reach>
std::vector<std::vector<int>> block_atom_lists(const cube &c, const int B, const std::vector<d3> &pos, const Reach &reach, const int extra = 0)
{
	const i3 n = c.get_sizes();
	const i3 nb = { (n[0] + B - 1) / B, (n[1] + B - 1) / B, (n[2] + B - 1) / B };
	std::vector<std::vector<int>> lists((size_t)nb[0] * nb[1] * nb[2]);
#pragma omp parallel for schedule(dynamic)
	for (int b = 0; b < (int)lists.size(); b++)
	{
		const int bx = b / (nb[1] * nb[2]), by = (b / nb[2]) % nb[1], bz = b % nb[2];
		const int lo[3] = { bx * B, by * B, bz * B };
		const int hi[3] = { std::min(lo[0] + B + extra, n[0]) - 1, std::min(lo[1] + B + extra, n[1]) - 1, std::min(lo[2] + B + extra, n[2]) - 1 };
		//the farthest point of the block's parallelepiped from its centre is a corner
		d3 corner[8], centre = { 0, 0, 0 };
		for (int k = 0; k < 8; k++)
		{
			corner[k] = c.get_pos(k & 1 ? hi[0] : lo[0], k & 2 ? hi[1] : lo[1], k & 4 ? hi[2] : lo[2]);
			for (int d = 0; d < 3; d++) centre[d] += corner[k][d] / 8;
		}
		double rb = 0;
		for (int k = 0; k < 8; k++) rb = std::max(rb, array_length(corner[k], centre));
		for (int a = 0; a < (int)pos.size(); a++)
			if (array_length(centre, pos[a]) <= reach(b, a) + rb + 1E-6)
				lists[b].push_back(a);
	}
	return lists;
}

} // namespace

void Calc_Spherical_Dens(
	cube &CubeSpher,
	const WFN &wavy,
	double radius,
	std::ostream &file,
	bool wrap)
{
	using namespace std;
	_time_point start = get_time();

	const vector<Thakkar> atom_models = make_thakkar_interpolators();
	const double radius_bohr = constants::ang2bohr(radius);
	const vector<atom> wavy_atoms = wavy.get_atoms();

#ifdef NOSPHERA2_USE_GPU
	//Only the elements present travel to the device; wrap sums periodic images with atomics and stays on the host
	if (!wrap && aux_density_gpu_enabled())
	{
		const int ncen = wavy.get_ncen();
		vector<double> ax(ncen), ay(ncen), az(ncen);
		vector<int> at_tab(ncen), tab_of_z(92, -1), tab_off{ 0 };
		vector<double> r_tab, rho_tab;
		for (int a = 0; a < ncen; a++)
		{
			const int z = wavy.get_atom_charge(a) - 1;
			if (tab_of_z[z] < 0)
			{
				tab_of_z[z] = (int)tab_off.size() - 1;
				const vec &r = atom_models[z].get_radial_dist(), &rho = atom_models[z].get_radial_density_table();
				r_tab.insert(r_tab.end(), r.begin(), r.end());
				rho_tab.insert(rho_tab.end(), rho.begin(), rho.end());
				tab_off.push_back((int)r_tab.size());
			}
			at_tab[a] = tab_of_z[z];
			const d3 pos = wavy.get_atom_pos(a);
			ax[a] = pos[0], ay[a] = pos[1], az[a] = pos[2];
		}
		const i3 n = CubeSpher.get_sizes();
		double origin[3], vectors[9];
		for (int i = 0; i < 3; i++)
		{
			origin[i] = CubeSpher.get_origin(i);
			for (int j = 0; j < 3; j++)
				vectors[3 * i + j] = CubeSpher.get_vector(i, j);
		}
		vector<double> out((size_t)n[0] * n[1] * n[2]);
		if (spherical_density_gpu_eval(n[0], n[1], n[2], origin, vectors, ncen, ax.data(), ay.data(), az.data(), at_tab.data(),
									   (int)tab_off.size() - 1, tab_off.data(), r_tab.data(), rho_tab.data(),
									   atom_models[0].get_lincr(), atom_models[0].get_start(), radius_bohr, out.data()))
		{
			if (!constants::hide_gpu_notes)
				file << "GPU in use: spherical density grid" << endl;
			const double *v = out.data();
			for (int x = 0; x < n[0]; x++)
				for (int y = 0; y < n[1]; y++)
					for (int z = 0; z < n[2]; z++)
						CubeSpher.set_value(x, y, z, *v++);
			_time_point end = get_time();
			print_time(start, end, file);
			return;
		}
	}
#endif

	if (!wrap)
	{
		//Per 8^3 block of points only the atoms that can reach it: Thakkar is exactly 0 past its table end, so a
		//dropped atom only ever added +0.0 and the in-radius test sees every atom within radius_bohr; values bit-identical.
		constexpr int B = 8;
		const i3 n = CubeSpher.get_sizes();
		const i3 nb = { (n[0] + B - 1) / B, (n[1] + B - 1) / B, (n[2] + B - 1) / B };
		vector<d3> apos(wavy_atoms.size());
		vec reach(wavy_atoms.size());
		for (size_t a = 0; a < wavy_atoms.size(); a++)
		{
			apos[a] = wavy_atoms[a].get_pos();
			reach[a] = std::max(atom_models[wavy_atoms[a].get_charge() - 1].get_radial_dist().back(), radius_bohr);
		}
		const vector<vector<int>> block_atoms = block_atom_lists(CubeSpher, B, apos, [&](int, int a) { return reach[a]; });
		CubeSpher.evaluate_on_grid(
			[&](const d3 &pos, const i3 &idx, const i3 &) {
				const vector<int> &nearby = block_atoms[((size_t)(idx[0] / B) * nb[1] + idx[1] / B) * nb[2] + idx[2] / B];
				bool inside = false;
				for (const int a : nearby)
					if (array_length(pos, wavy_atoms[a].get_pos()) < radius_bohr) { inside = true; break; }
				if (!inside)
					return 0.0;
				double dens_all = 0.0;
				for (const int a : nearby)
					dens_all += atom_models[wavy_atoms[a].get_charge() - 1].get_interpolated_density(array_length(pos, wavy_atoms[a].get_pos()));
				return dens_all;
			});
	}
	else
		evaluate_cube_in_radius(
			CubeSpher,
			wrap,
			wavy_atoms,
			radius_bohr,
			[&](const d3 &pos) {
				vector<double> dists(wavy.get_ncen(), 0.0);
				for (int a = 0; a < wavy.get_ncen(); a++)
					dists[a] = array_length(pos, wavy.get_atom_pos(a));

				double dens_all = 0.0;
				for (int a = 0; a < wavy.get_ncen(); a++)
					dens_all += atom_models[wavy.get_atom_charge(a) - 1].get_interpolated_density(dists[a]);
				return dens_all;
			});

	_time_point end = get_time();
	print_time(start, end, file);
};

void Calc_Static_Def(
	cube &CubeDEF,
	cube &CubeRho,
	const WFN &wavy,
	double radius,
	std::ostream &file,
	bool wrap)
{
	using namespace std;
	_time_point start = get_time();
	vector<Thakkar> atoms;
	atoms.reserve(wavy.get_ncen());
	for (int a = 0; a < wavy.get_ncen(); a++)
		atoms.emplace_back(wavy.get_atom_charge(a));
	const double radius_bohr = constants::ang2bohr(radius);
	const vector<atom> wavy_atoms = wavy.get_atoms();

	evaluate_cube_in_radius_mapped(
		CubeDEF,
		wrap,
		wavy_atoms,
		radius_bohr,
		[&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
			vector<double> dists(wavy.get_ncen(), 0.0);
			for (int a = 0; a < wavy.get_ncen(); a++)
				dists[a] = array_length(pos, wavy.get_atom_pos(a));

			double dens_all = 0.0;
			for (int a = 0; a < wavy.get_ncen(); a++)
				dens_all += atoms[a].get_radial_density(dists[a]);

			dens_all -= CubeRho.get_value(mapped_idx[0], mapped_idx[1], mapped_idx[2]);
			return -dens_all;
		});

	_time_point end = get_time();
	print_time(start, end, file);
};

void Calc_Static_Def(
	std::vector<cube> &Cubes,
	const WFN &wavy,
	double radius,
	std::ostream &file,
	bool wrap)
{
	using namespace std;
	_time_point start = get_time();
	const double radius_bohr = constants::ang2bohr(radius);
	const vector<atom> wavy_atoms = wavy.get_atoms();

	evaluate_cube_in_radius_mapped(
		Cubes[cube_type::DEF],
		wrap,
		wavy_atoms,
		radius_bohr,
		[&](const d3 &, const i3 &, const i3 &mapped_idx) {
			const int x = mapped_idx[0];
			const int y = mapped_idx[1];
			const int z = mapped_idx[2];
			const double rho = Cubes[cube_type::Rho].get_value(x, y, z);
			const double spher = Cubes[cube_type::spherical_density].get_value(x, y, z);
			return rho - spher;
		});

	_time_point end = get_time();
	print_time(start, end, file);
};

void Calc_Hirshfeld(
	std::vector<cube> &Cubes,
	const WFN &wavy,
	double radius,
	int ignore_atom,
	std::ostream &file,
	bool wrap)
{
	using namespace std;
	_time_point start = get_time();

	vector<Thakkar> atoms;
	atoms.reserve(wavy.get_ncen());
	for (int a = 0; a < wavy.get_ncen(); a++)
		atoms.emplace_back(wavy.get_atom_charge(a));
	const double radius_bohr = constants::ang2bohr(radius);
	const vector<atom> focus_atom{ wavy.get_atoms()[ignore_atom] };

	evaluate_cube_in_radius_mapped(
		Cubes[cube_type::HDEF],
		wrap,
		focus_atom,
		radius_bohr,
		[&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
			double dens_choice = 0.0;
			double dens_all = 0.0;
			for (int a = 0; a < wavy.get_ncen(); a++)
			{
				const double dist = array_length(pos, wavy.get_atom_pos(a));
				const double temp = atoms[a].get_radial_density(dist);
				if (ignore_atom == a)
					dens_choice = temp;
				dens_all += temp;
			}

			if (std::abs(dens_all) < 1E-20)
				return 0.0;

			const double rho = Cubes[cube_type::Rho].get_value(mapped_idx[0], mapped_idx[1], mapped_idx[2]);
			return dens_choice / dens_all * rho - dens_choice;
		});

	_time_point end = get_time();
	print_time(start, end, file);
};

void Calc_Hirshfeld(
	cube &CubeHDEF,
	cube &CubeRho,
	cube &CubeSpherical,
	const WFN &wavy,
	double radius,
	int ignore_atom,
	std::ostream &file,
	bool wrap)
{
	using namespace std;
	_time_point start = get_time();
	Thakkar atom_model(wavy.get_atom_charge(ignore_atom));
	const double radius_bohr = constants::ang2bohr(radius);

	evaluate_cube_near_atom_mapped(
		CubeHDEF,
		wrap,
		wavy.get_atom_pos(ignore_atom),
		radius_bohr,
		[&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
			const double dist = array_length(pos, wavy.get_atom_pos(ignore_atom));
			const int x = mapped_idx[0];
			const int y = mapped_idx[1];
			const int z = mapped_idx[2];
			const double spherical = CubeSpherical.get_value(x, y, z);
			if (std::abs(spherical) < 1E-20)
				return 0.0;

			const double dens_choice = atom_model.get_radial_density(dist);
			return (dens_choice / spherical * CubeRho.get_value(x, y, z)) - dens_choice;
		});

	_time_point end = get_time();
	print_time(start, end, file);
};

void Calc_Hirshfeld_atom(
	std::vector<cube> &Cubes,
	const WFN &wavy,
	double radius,
	int ignore_atom,
	std::ostream &file,
	bool wrap)
{
	using namespace std;
	_time_point start = get_time();
	Thakkar atom_model(wavy.get_atom_charge(ignore_atom));
	const double radius_bohr = constants::ang2bohr(radius);

	evaluate_cube_near_atom_mapped(
		Cubes[cube_type::Hirsh],
		wrap,
		wavy.get_atom_pos(ignore_atom),
		radius_bohr,
		[&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
			const double dist = array_length(pos, wavy.get_atom_pos(ignore_atom));
			const double dens_choice = atom_model.get_radial_density(dist);
			const int x = mapped_idx[0];
			const int y = mapped_idx[1];
			const int z = mapped_idx[2];
			const double spherical = Cubes[cube_type::spherical_density].get_value(x, y, z);
			if (std::abs(spherical) < 1E-20)
				return 0.0;

			return dens_choice / spherical * Cubes[cube_type::Rho].get_value(x, y, z);
		});

	_time_point end = get_time();
	print_time(start, end, file);
};


void Calc_RhoEli(
	cube &CubeRho,
	cube &CubeEli,
	const WFN &wavy,
	double radius,
	const density_field *field)
{
	using namespace std;
	err_checkf(CubeRho.get_size(0) == CubeEli.get_size(0) && CubeRho.get_size(1) == CubeEli.get_size(1) && CubeRho.get_size(2) == CubeEli.get_size(2), "Cube sizes do not match", std::cout);

	const double radius_bohr = constants::ang2bohr(radius);
	const vector<atom> atoms = wavy.get_atoms();
	CubeEli.set_zero();

	evaluate_cube_in_radius_mapped(
		CubeRho,
		false,
		atoms,
		radius_bohr,
		[&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
			double rho = 0.0;
			double eli = 0.0;
			if (!field) wavy.computeRhoELI(pos, rho, eli);
			else {
				rho = field->rho(pos);
				if (wavy.get_nmo() > 0) eli = wavy.computeELI(pos);
			}
			CubeEli.set_value(mapped_idx[0], mapped_idx[1], mapped_idx[2], eli);
			return rho;
		});

};

void Calc_Rho_spherical_harmonics(
	cube &CubeRho,
	const WFN &wavy,
	std::ostream &file)
{
	using namespace std;
	_time_point start = get_time();

	ProgressBar *progress = new ProgressBar(CubeRho.get_size(0), 50, "=", " ", "Calculating Rho");

#pragma omp parallel shared(CubeRho)
	{
		vec2 d(wavy.get_ncen());
		for (int i = 0; i < wavy.get_ncen(); i++)
			d[i].resize(16, 0.0);
		const int n = wavy.get_nmo(true);
		vec phi(n, 0.0);
#pragma omp for schedule(dynamic)
		for (int i = 0; i < CubeRho.get_size(0); i++)
		{
			for (int j = 0; j < CubeRho.get_size(1); j++)
				for (int k = 0; k < CubeRho.get_size(2); k++)
					CubeRho.set_value(i, j, k, wavy.compute_dens(CubeRho.get_pos(i, j, k), d, phi));
			progress->update();
		}
	}
	delete (progress);

	_time_point end = get_time();
	print_time(start, end, file);
};

//The grid points within r of a centre, strips of 4x4 lines gathered in 4x4x4 tiles, so each 64-point block of a batch
//kernel is a compact cube of space; f(idx, x, y, z) per 65536 points or so. With wrap the points of the 27 periodic
//images, idx mapped into the cell
template <typename F>
static void tiled_points(const cube &C, const WFN &wavy, const double radius, const bool wrap, F &&f)
{
	const double r = constants::ang2bohr(radius);
	const std::vector<d3> centres = source_positions(wavy);
	const i3 n = C.get_sizes(), o = wrap ? n : i3{ 0, 0, 0 }, N = wrap ? i3{ 3 * n[0], 3 * n[1], 3 * n[2] } : n;
	const int nt = (N[1] + 3) / 4, U = (N[0] + 3) / 4 * nt, G = std::max(1, 65536 / std::max(1, 16 * N[2]));
	std::vector<std::vector<i3>> keep(G);
	std::vector<i3> idx;
	vec x, y, z;
	for (int s0 = 0; s0 < U; s0 += G) {
		const int ng = std::min(G, U - s0);
#pragma omp parallel for schedule(dynamic)
		for (int g = 0; g < ng; g++) {
			keep[g].clear();
			const int a = 4 * ((s0 + g) / nt), b = 4 * ((s0 + g) % nt);
			for (int c = 0; c < N[2]; c += 4)
				for (int i = a; i < std::min(a + 4, N[0]); i++)
					for (int j = b; j < std::min(b + 4, N[1]); j++)
						for (int k = c; k < std::min(c + 4, N[2]); k++)
							if (near_any(C.get_pos(i - o[0], j - o[1], k - o[2]), centres, r)) keep[g].push_back({ i, j, k });
		}
		idx.clear(), x.clear(), y.clear(), z.clear();
		for (int g = 0; g < ng; g++)
			for (int p = 0; p < (int)keep[g].size(); p++) {
				const i3 &t = keep[g][p];
				const d3 pos = C.get_pos(t[0] - o[0], t[1] - o[1], t[2] - o[2]);
				x.push_back(pos[0]), y.push_back(pos[1]), z.push_back(pos[2]);
				idx.push_back({ t[0] % n[0], t[1] % n[1], t[2] % n[2] });
			}
		f(idx, x, y, z);
	}
}

void Calc_Rho(cube &CubeRho, const WFN &wavy, double radius, std::ostream &file, bool wrap)
{
	_time_point start = get_time();
	CubeRho.set_zero();
	vec rho;
	tiled_points(CubeRho, wavy, radius, wrap, [&](const std::vector<i3> &idx, const vec &x, const vec &y, const vec &z) {
		rho.resize(idx.size());
		wavy.compute_dens_batch((int)idx.size(), x.data(), y.data(), z.data(), rho.data());
		for (int p = 0; p < (int)idx.size(); p++) {
			const i3 &i = idx[p];
			CubeRho.set_value(i[0], i[1], i[2], CubeRho.get_value(i[0], i[1], i[2]) + rho[p]);
		}
	});
	_time_point end = get_time();
	print_time(start, end, file);
}

void Calc_MO_spherical_harmonics(
	cube &CubeMO,
	const WFN &wavy,
	int MO,
	std::ostream &file,
	bool nodate)
{
	using namespace std;
	_time_point start = get_time();
	ProgressBar *progress = NULL;
	if (!nodate)
		progress = new ProgressBar(CubeMO.get_size(0), 50, "=", " ", "Calculating Values");

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < CubeMO.get_size(0); i++)
	{
		for (int j = 0; j < CubeMO.get_size(1); j++)
			for (int k = 0; k < CubeMO.get_size(2); k++)
				CubeMO.set_value(i, j, k, wavy.compute_MO_spherical(CubeMO.get_pos(i, j, k), MO));
		if (!nodate)
			progress->update();
	}
	if (!nodate)
	{
		delete (progress);

		_time_point end = get_time();
		print_time(start, end, file);
	}
};

void Calc_S_Rho(
	cube &Cube_S_Rho,
	const WFN &wavy,
	std::ostream &file,
	bool &nodate)
{
	using namespace std;
	_time_point start = get_time();
	ProgressBar *progress = NULL;
	if (!nodate)
		progress = new ProgressBar(Cube_S_Rho.get_size(0), 50, "=", " ", "Calculating Values");

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < Cube_S_Rho.get_size(0); i++)
	{
		vec2 d(wavy.get_ncen(), vec(16, 0.0));
		vec phi(wavy.get_nmo(), 0.0);
		for (int j = 0; j < Cube_S_Rho.get_size(1); j++)
			for (int k = 0; k < Cube_S_Rho.get_size(2); k++)
				Cube_S_Rho.set_value(i, j, k, wavy.compute_spin_dens(Cube_S_Rho.get_pos(i, j, k), d, phi));
		if (!nodate)
			progress->update();
	}
	if (!nodate)
	{
		delete (progress);

		_time_point end = get_time();
		print_time(start, end, file);
	}
};

//ELI-D(aa), ELI-D(bb) and the triplet ELI-D of [III] at every point from one eli_family::spin_fields call; a member not in
//want stays unloaded and unwritten. Where rho_s or g is not positive (a spin channel without electrons there, far tails)
//eli_d gives 0, and sanitize_finite keeps that cube's values finite as -eli does
void Calc_Eli_Spin(
	std::vector<cube> &cubes,
	const WFN &wavy,
	std::ostream &file,
	bool &nodate)
{
	using namespace std;
	_time_point start = get_time();
	const double factor = eli_family::triplet_density_factor(wavy);
	cube &aa = cubes[cube_type::Eli_aa], &bb = cubes[cube_type::Eli_bb], &tr = cubes[cube_type::Eli_triplet];
	ProgressBar *progress = NULL;
	if (!nodate)
		progress = new ProgressBar(aa.get_size(0), 50, "=", " ", "Calculating Values");
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < aa.get_size(0); i++)
	{
		eli_family::SpinFields f;
		for (int j = 0; j < aa.get_size(1); j++)
			for (int k = 0; k < aa.get_size(2); k++)
			{
				eli_family::spin_fields(wavy, aa.get_pos(i, j, k), f);
				aa.set_value(i, j, k, sanitize_finite(eli_family::eli_d(f.rho[0], eli_family::g_same_spin(f, 0))));
				if (bb.get_loaded()) bb.set_value(i, j, k, sanitize_finite(eli_family::eli_d(f.rho[1], eli_family::g_same_spin(f, 1))));
				if (tr.get_loaded()) tr.set_value(i, j, k, sanitize_finite(eli_family::eli_d(factor * (f.rho[0] + f.rho[1]), eli_family::g_triplet(f))));
			}
		if (!nodate)
			progress->update();
	}
	if (!nodate)
	{
		delete (progress);
		_time_point end = get_time();
		print_time(start, end, file);
	}
}

void Calc_Prop(
	std::vector<cube> &Cubes,
	const WFN &wavy,
	double radius,
	std::ostream &file,
	bool test,
	bool wrap)
{
	_time_point start = get_time();
	const bool rdg = Cubes[cube_type::RDG].get_loaded(), lap = Cubes[cube_type::Lap].get_loaded(), elf = Cubes[cube_type::Elf].get_loaded(), eli = Cubes[cube_type::Eli].get_loaded();
	const bool hess = rdg || lap;
	const int s = hess ? 11 : 5;
	cube rho_contrib(Cubes[cube_type::Rho]);
	rho_contrib.set_zero();
	vec F;
	tiled_points(rho_contrib, wavy, radius, wrap, [&](const std::vector<i3> &idx, const vec &x, const vec &y, const vec &z) {
		const int np = (int)idx.size();
		F.resize((size_t)s * np);
		wavy.compute_fields_batch(np, x.data(), y.data(), z.data(), hess, F.data());
		for (int p = 0; p < np; p++) {
			const double *q = F.data() + (size_t)s * p, rho = q[0], g2 = q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
			const i3 &i = idx[p];
			auto add = [&](cube &c, const double v) { c.set_value(i[0], i[1], i[2], c.get_value(i[0], i[1], i[2]) + sanitize_finite(v)); };
			if (elf)
				add(Cubes[cube_type::Elf], rho > 0 ? 1 / (1 + pow(constants::ctelf * pow(rho, constants::c_m53) * (q[4] * 0.5 - 0.125 * g2 / rho), 2)) : 0.0);
			if (eli) {
				const double g = rho * q[4] - 0.25 * g2;
				add(Cubes[cube_type::Eli], rho > 0 && g > 0 ? 0.5 * rho * pow(48 / g, constants::c_38) : 0.0);
			}
			if (lap)
				add(Cubes[cube_type::Lap], q[5] + q[6] + q[7]);
			if (rdg) {
				add(Cubes[cube_type::RDG], rho > 0 ? constants::alpha_coef * sqrt(g2) / pow(rho, constants::c_43) : 0.0);
				double H[9] = { q[5], q[8], q[9], q[8], q[6], q[10], q[9], q[10], q[7] };
				add(rho_contrib, get_lambda_1(H) < 0 ? -rho : rho);
			}
		}
	});
	// Calc_Rho already filled Rho; only the RDG run replaces it by sign(lambda2)*rho.
	if (rdg)
		finish_signed_rho(Cubes, rho_contrib);

	if (!test)
	{
		_time_point end = get_time();
		print_time(start, end, file);
	}
};

void Calc_ESP(
	cube &CubeESP,
	const WFN &wavy,
	double radius,
	bool no_date,
	std::ostream &file,
	bool wrap)
{
	using namespace std;
	_time_point start = get_time();

	const WFN::ESP_pairs pairs = wavy.build_ESP_pairs();
	const double radius_bohr = constants::ang2bohr(radius);
	const vector<atom> atoms = wavy.get_atoms();

	evaluate_cube_in_radius_batched(
		CubeESP,
		wrap,
		atoms,
		radius_bohr,
		[&](const std::vector<d3> &points, double *out) {
			wavy.computeESP_batch(points, pairs, out);
		});

	if (!no_date)
	{
		_time_point end = get_time();
		print_time(start, end, file);
	}
};

void Calc_MO(
	cube &CubeMO,
	int mo,
	const WFN &wavy,
	double radius,
	std::ostream &file,
	bool wrap)
{
	using namespace std;
	err_checkf(mo >= 0 && mo < wavy.get_nmo(), "MO index " + to_string(mo) + " out of range, the wavefunction has " + to_string(wavy.get_nmo()) + " MOs (0 to " + to_string(wavy.get_nmo() - 1) + ")!", file);
	_time_point start = get_time();
	const double radius_bohr = constants::ang2bohr(radius);
	const vector<atom> atoms = wavy.get_atoms();

	evaluate_cube_in_radius(
		CubeMO,
		wrap,
		atoms,
		radius_bohr,
		[&](const d3 &pos) {
			return wavy.computeMO(pos, mo);
		});

	_time_point end = get_time();
	print_time(start, end, file);
};

bool find_frontier_orbitals(
	const WFN &wavy,
	int &homo,
	int &lumo,
	bool &unrestricted)
{
	homo = -1;
	lumo = -1;
	unrestricted = wavy.get_is_unrestricted();

	const int n = wavy.get_nmo();
	if (n <= 0)
		return false;

	// Some formats (notably plain .wfn) carry no orbital energies at all. Only
	// trust energy ordering if at least one energy is non-zero.
	bool have_energies = false;
	for (int i = 0; i < n; i++) {
		if (wavy.get_MO_energy(i) != 0.0) {
			have_energies = true;
			break;
		}
	}

	// An occupation of exactly zero marks a virtual orbital. Use a small
	// threshold so that near-zero natural-orbital occupations from correlated
	// methods are still treated as occupied rather than virtual.
	constexpr double occ_threshold = 1.0e-6;

	double best_occ_energy = 0.0;
	double best_vir_energy = 0.0;
	for (int i = 0; i < n; i++) {
		const double occ = wavy.get_MO_occ(i);
		const double ener = wavy.get_MO_energy(i);
		if (occ > occ_threshold) {
			if (homo == -1 || (have_energies ? (ener > best_occ_energy) : (i > homo))) {
				homo = i;
				best_occ_energy = ener;
			}
		}
		else {
			if (lumo == -1 || (have_energies ? (ener < best_vir_energy) : (i < lumo))) {
				lumo = i;
				best_vir_energy = ener;
			}
		}
	}

	return (homo != -1 && lumo != -1);
}

void Calc_Fukui(
	std::vector<cube> &Cubes,
	const WFN &wavy,
	int homo,
	int lumo,
	double radius,
	std::ostream &file,
	bool wrap)
{
	using namespace std;
	_time_point start = get_time();
	const double radius_bohr = constants::ang2bohr(radius);
	const vector<atom> atoms = wavy.get_atoms();

	if (!Cubes[cube_type::Fukui_plus].get_loaded())
		return;

	// The caller is expected to have obtained these from find_frontier_orbitals,
	// which returns false rather than a negative index when no frontier pair
	// exists. Guard anyway: a negative index reaches computeMO as an
	// out-of-bounds orbital lookup, which is undefined rather than merely wrong.
	err_checkf(homo >= 0 && homo < wavy.get_nmo(), "Invalid HOMO index " + std::to_string(homo) + " for Fukui calculation!", file);
	err_checkf(lumo >= 0 && lumo < wavy.get_nmo(), "Invalid LUMO index " + std::to_string(lumo) + " for Fukui calculation!", file);

	// Deliberately TWO grid passes, each returning its value, rather than one
	// pass that fills all four cubes from a single lambda.
	//
	// The tempting single-pass version writes the other three cubes from inside
	// the functor with a get_value/set_value read-modify-write. That is a data
	// race in the wrapped (periodic) path: evaluate_on_grid parallelises over
	// the outer index and several raw indices wrap onto the SAME voxel, so two
	// threads can read-modify-write one value concurrently. Only the functor's
	// RETURN value is protected, by the "#pragma omp atomic" in evaluate_on_grid.
	//
	// Two passes cost exactly the same: one computeMO per point each, versus two
	// per point in a single pass. The orbital evaluation dominates, so this is
	// free, and f0 and the dual descriptor then follow elementwise from the two
	// finished cubes where nothing is shared between threads.
	evaluate_cube_in_radius(
		Cubes[cube_type::Fukui_plus],
		wrap,
		atoms,
		radius_bohr,
		[&](const d3 &pos) {
			// computeMO returns the orbital AMPLITUDE psi, so square it here.
			const double psi_lumo = wavy.computeMO(pos, lumo);
			return sanitize_finite(psi_lumo * psi_lumo);
		});

	evaluate_cube_in_radius(
		Cubes[cube_type::Fukui_minus],
		wrap,
		atoms,
		radius_bohr,
		[&](const d3 &pos) {
			const double psi_homo = wavy.computeMO(pos, homo);
			return sanitize_finite(psi_homo * psi_homo);
		});

	// f0 = (f+ + f-) / 2 and df = f+ - f-, elementwise. Each voxel is written by
	// exactly one thread here, so no synchronisation is needed.
	const i3 sizes = Cubes[cube_type::Fukui_plus].get_sizes();
#pragma omp parallel for schedule(dynamic)
	for (int x = 0; x < sizes[0]; x++)
		for (int y = 0; y < sizes[1]; y++)
			for (int z = 0; z < sizes[2]; z++)
			{
				const double f_plus = Cubes[cube_type::Fukui_plus].get_value(x, y, z);
				const double f_minus = Cubes[cube_type::Fukui_minus].get_value(x, y, z);
				Cubes[cube_type::Fukui_zero].set_value(x, y, z, 0.5 * (f_plus + f_minus));
				Cubes[cube_type::Dual_Descriptor].set_value(x, y, z, f_plus - f_minus);
			}
	// All four cubes were allocated together with loaded = true, so no flag needs
	// setting here; evaluate_on_grid has already re-set it for f+ and f-.

	_time_point end = get_time();
	print_time(start, end, file);
};

namespace {

// Replace the WFN_DENSITY column of a built grid by |psi_mo(r)|^2: calculatePartitionedCharges() integrates
// that column against all five partition weights without recomputing it, giving the condensed Fukui functions.
// wavy must be the original wavefunction with its virtual orbitals; the grid is built from a pruned copy.
void fill_density_column_with_orbital(GridManager &gm, const WFN &wavy, int mo)
{
	GridData &gd = gm.getGridData();
	const bool helper = gm.getNeedsHelper();
	vec3 &grids = helper ? gd.helper_grids : gd.atomic_grids;
	const ivec &npts = helper ? gd.helper_num_points_per_atom : gd.num_points_per_atom;

	const int n_grids = static_cast<int>(grids.size());
	for (int g = 0; g < n_grids; g++) {
		vec2 &atom_grid = grids[g];
		const int num_points = npts[g];
		const double *x = atom_grid[GridData::GridIndex::X].data();
		const double *y = atom_grid[GridData::GridIndex::Y].data();
		const double *z = atom_grid[GridData::GridIndex::Z].data();
		double *dens = atom_grid[GridData::GridIndex::WFN_DENSITY].data();
#pragma omp parallel for schedule(dynamic, 64)
		for (int p = 0; p < num_points; p++) {
			const double psi = wavy.computeMO({ x[p], y[p], z[p] }, mo);
			dens[p] = psi * psi;
		}
	}
}

}

CondensedFukuiResults Calc_Condensed_Fukui(
	const WFN &wavy,
	int homo,
	int lumo,
	const cell &unit_cell,
	int accuracy,
	std::ostream &file)
{
	CondensedFukuiResults out;
	const int ncen = wavy.get_ncen();
	if (ncen <= 0 || homo < 0 || lumo < 0)
		return out;

	GridConfiguration config;
	config.accuracy = accuracy;
	config.pbc = 0;
	config.debug = false;
	// Forces calculatePartitionedCharges down its five-wide branch, so one pass
	// yields Becke, Hirshfeld, TFVC, MBIS and EMBIS together.
	config.all_charges = true;

	ivec atom_list(ncen);
	for (int i = 0; i < ncen; i++)
		atom_list[i] = i;

	out.labels.resize(ncen);
	const std::vector<atom> atoms = wavy.get_atoms();
	for (int i = 0; i < ncen; i++)
		out.labels[i] = atoms[i].get_label();

	GridManager gm(config);

	// Build the grid and the five weight sets from a copy with the virtual
	// orbitals removed. Two reasons, and the first is a crash rather than a
	// preference: compute_dens sizes its scratch array by get_nmo(true), so
	// leaving several hundred virtuals in place overruns it. Second, the
	// partition weights are properties of the GROUND-STATE density - MBIS and
	// EMBIS are refined self-consistently against it - so they must be built
	// from rho, not from a frontier orbital. Only the integrand changes below.
	WFN temp = wavy;
	temp.delete_unoccupied_MOs();
	gm.setup3DGridsForMolecule(temp, atom_list, {}, unit_cell, false, file);

	// ECP cores are added to the populations by calculatePartitionedCharges.
	// That is right for a charge and wrong for a Fukui function, which knows
	// nothing about core electrons, so it is subtracted back out below.
	const bool has_ecps = temp.get_has_ECPs();

	auto accumulate = [&](int mo, vec2 &dest) {
		fill_density_column_with_orbital(gm, wavy, mo);
		PartitionResults res = gm.calculatePartitionedCharges(temp, unit_cell);
		dest.resize(5);
		for (int s = 0; s < 5; s++) {
			dest[s].resize(ncen, 0.0);
			for (int a = 0; a < ncen && a < static_cast<int>(res.atom_charges[s].size()); a++) {
				double v = res.atom_charges[s][a];
				if (has_ecps)
					v -= temp.get_atom_ECP_electrons(a);
				dest[s][a] = v;
			}
		}
	};

	accumulate(lumo, out.f_plus);
	accumulate(homo, out.f_minus);
	out.valid = true;
	return out;
}

void print_condensed_fukui(
	const CondensedFukuiResults &r,
	std::ostream &file)
{
	if (!r.valid)
		return;
	using namespace std;
	// Column order matches PartitionResults::CHARGE_ORDER.
	static const char *scheme[5] = { "Becke", "TFVC", "Hirshfeld", "MBIS", "EMBIS" };
	static const int order[5] = {
		PartitionResults::CHARGE_ORDER::S_HIRSH,
		PartitionResults::CHARGE_ORDER::S_BECKE,
		PartitionResults::CHARGE_ORDER::S_TFVC,
		PartitionResults::CHARGE_ORDER::S_MBIS,
		PartitionResults::CHARGE_ORDER::S_EMBIS
	};
	static const char *order_name[5] = { "Hirshfeld", "Becke", "TFVC", "MBIS", "EMBIS" };
	(void)scheme;

	const size_t n = r.labels.size();
	file << "\nCondensed (atom-summed) Fukui functions" << endl;
	file << "  f+_A = integral of w_A(r) |psi_LUMO(r)|^2   (A is attacked by a nucleophile -> electrophilic site)" << endl;
	file << "  f-_A = integral of w_A(r) |psi_HOMO(r)|^2   (A is attacked by an electrophile -> nucleophilic site)" << endl;
	file << "  df_A = f+_A - f-_A                          (> 0 electrophilic, < 0 nucleophilic)" << endl;
	file << "  Each column uses a different atomic partition w_A. Sum over atoms is 1 for each of f+ and f-." << endl;

	for (int quantity = 0; quantity < 3; quantity++) {
		const char *qname = (quantity == 0) ? "f+" : (quantity == 1) ? "f-" : "dual descriptor (f+ - f-)";
		file << "\n  " << qname << "\n";
		file << "  " << setw(8) << left << "Atom" << right;
		for (int s = 0; s < 5; s++)
			file << setw(13) << order_name[s];
		file << "\n";
		file << "  " << string(8 + 5 * 13, '-') << "\n";
		vec totals(5, 0.0);
		for (size_t a = 0; a < n; a++) {
			file << "  " << setw(8) << left << r.labels[a] << right;
			for (int s = 0; s < 5; s++) {
				const int idx = order[s];
				double v = 0.0;
				if (quantity == 0)
					v = r.f_plus[idx][a];
				else if (quantity == 1)
					v = r.f_minus[idx][a];
				else
					v = r.f_plus[idx][a] - r.f_minus[idx][a];
				totals[s] += v;
				file << setw(13) << fixed << setprecision(5) << v;
			}
			file << "\n";
		}
		file << "  " << string(8 + 5 * 13, '-') << "\n";
		file << "  " << setw(8) << left << "sum" << right;
		for (int s = 0; s < 5; s++)
			file << setw(13) << fixed << setprecision(5) << totals[s];
		file << "\n";
	}
	file << endl;
}

void fukui_analysis(options &opt, std::ostream &log2)
{
	using namespace std;
	log2 << NoSpherA2_message(opt.no_date);
	if (!opt.no_date)
		log2 << build_date;

	err_checkf(opt.wfn != "", "Error, no wfn file specified! Use -fukui_analysis <wfn> or -wfn <wfn>.", log2);
	WFN wavy(opt.wfn);
	log2 << "\nConceptual-DFT reactivity analysis of " << opt.wfn.string() << endl;
	citations::cite(citations::Method::Fukui, log2);
	log2 << "Read " << wavy.get_ncen() << " atoms and " << wavy.get_nmo()
		 << " molecular orbitals (" << wavy.get_nmo(true) << " occupied)." << endl;

	int homo = -1, lumo = -1;
	bool unrestricted = false;
	if (!find_frontier_orbitals(wavy, homo, lumo, unrestricted))
	{
		log2 << "\nERROR: no HOMO/LUMO pair found (HOMO index " << homo
			 << ", LUMO index " << lumo << ")." << endl;
		if (lumo == -1)
			log2 << "This wavefunction stores no virtual orbitals, so no Fukui function can be\n"
				 << "formed from it. Plain .wfn files and many .wfx files keep only the occupied\n"
				 << "orbitals; use a .gbw, .molden or .fchk instead." << endl;
		return;
	}

	log2 << "HOMO = MO " << homo << "  energy " << fixed << setprecision(6) << wavy.get_MO_energy(homo)
		 << "  occupation " << setprecision(3) << wavy.get_MO_occ(homo) << endl;
	log2 << "LUMO = MO " << lumo << "  energy " << setprecision(6) << wavy.get_MO_energy(lumo)
		 << "  occupation " << setprecision(3) << wavy.get_MO_occ(lumo) << endl;
	const double gap = wavy.get_MO_energy(lumo) - wavy.get_MO_energy(homo);
	log2 << "HOMO-LUMO gap: " << setprecision(6) << gap << " Hartree ("
		 << setprecision(3) << gap * constants::keV_per_hartree * 1000.0 << " eV)" << endl;
	if (unrestricted)
	{
		// Say which manifold actually supplied the pair, and what each manifold's
		// own frontier is. "Taken across both manifolds" on its own is not
		// actionable: for an open-shell system the answer usually comes entirely
		// from one spin, and which one changes the chemistry being described.
		static const char *spin_name[2] = { "alpha", "beta" };
		const int homo_op = wavy.get_MO_op(homo);
		const int lumo_op = wavy.get_MO_op(lumo);
		log2 << "\nUnrestricted wavefunction." << endl;
		log2 << "  The frontier pair is the globally highest occupied and globally lowest virtual\n"
			 << "  spin orbital, i.e. the electron that is easiest to remove and the level an added\n"
			 << "  electron would enter, regardless of spin. That is the spin-UNRESOLVED Fukui\n"
			 << "  function; it is not a spin-polarised (f_N / f_S) treatment." << endl;
		log2 << "  Chosen HOMO is " << spin_name[homo_op == 1 ? 1 : 0]
			 << ", chosen LUMO is " << spin_name[lumo_op == 1 ? 1 : 0] << "." << endl;

		// Per-manifold frontiers, so a suspicious global answer can be checked.
		for (int op = 0; op < 2; op++)
		{
			int h = -1, l = -1;
			double he = 0.0, le = 0.0;
			for (int i = 0; i < wavy.get_nmo(); i++)
			{
				if (wavy.get_MO_op(i) != op)
					continue;
				const double o = wavy.get_MO_occ(i);
				const double e = wavy.get_MO_energy(i);
				if (o > 1.0e-6) {
					if (h == -1 || e > he) { h = i; he = e; }
				}
				else {
					if (l == -1 || e < le) { l = i; le = e; }
				}
			}
			if (h == -1 && l == -1)
				continue;
			log2 << "  " << setw(5) << spin_name[op] << ": HOMO = MO ";
			if (h >= 0) log2 << h << " (" << fixed << setprecision(6) << he << ")"; else log2 << "none";
			log2 << ", LUMO = MO ";
			if (l >= 0) log2 << l << " (" << fixed << setprecision(6) << le << ")"; else log2 << "none";
			log2 << "\n";
		}
		if (homo_op == lumo_op)
			log2 << "  Both come from the " << spin_name[homo_op == 1 ? 1 : 0]
				 << " manifold, so this is effectively a " << spin_name[homo_op == 1 ? 1 : 0]
				 << "-only Fukui function." << endl;
		else
			log2 << "  They come from DIFFERENT manifolds, so f+ and f- describe different spin\n"
				 << "  channels. Interpret the dual descriptor with care." << endl;
		log2 << "  For a high-spin system the lowest virtual is often the spatial partner of a\n"
			 << "  singly-occupied orbital, which makes f+ and f- describe nearly the same region\n"
			 << "  and drives the dual descriptor toward zero. That is a property of the\n"
			 << "  approximation, not a statement about the molecule." << endl;
	}

	const CondensedFukuiResults condensed =
		Calc_Condensed_Fukui(wavy, homo, lumo, cell(), opt.accuracy, log2);
	if (!condensed.valid)
	{
		log2 << "Condensed Fukui functions could not be calculated." << endl;
		return;
	}
	print_condensed_fukui(condensed, log2);

	// The "so what" line. Ranked on the Hirshfeld column because that is the
	// partition the condensed-Fukui literature uses, so it is the one a reader
	// can compare against published numbers.
	const int H = PartitionResults::CHARGE_ORDER::S_HIRSH;
	const size_t n = condensed.labels.size();
	if (n > 0)
	{
		size_t most_electrophilic = 0, most_nucleophilic = 0;
		double best_pos = -1e300, best_neg = 1e300;
		for (size_t a = 0; a < n; a++)
		{
			const double df = condensed.f_plus[H][a] - condensed.f_minus[H][a];
			if (df > best_pos) { best_pos = df; most_electrophilic = a; }
			if (df < best_neg) { best_neg = df; most_nucleophilic = a; }
		}
		log2 << "Summary (Hirshfeld partition):" << endl;
		log2 << "  Most ELECTROPHILIC site (a nucleophile attacks here): "
			 << condensed.labels[most_electrophilic]
			 << "   df = " << fixed << setprecision(5) << best_pos << endl;
		log2 << "  Most NUCLEOPHILIC site (an electrophile attacks here): "
			 << condensed.labels[most_nucleophilic]
			 << "   df = " << fixed << setprecision(5) << best_neg << endl;
		log2 << "\nNote: these atom-summed values, not the cube maxima, are what carry the\n"
			 << "chemistry. A point-wise Fukui function peaks near the heaviest nucleus\n"
			 << "regardless of where the molecule actually reacts." << endl;
	}

	// Machine-readable copy beside the wavefunction, same format as the -fukui run.
	const std::filesystem::path summary_path =
		(wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_fukui.dat";
	ofstream summary(summary_path, ios::out);
	summary << "# NoSpherA2 condensed Fukui analysis\n";
	summary << "# Frontier-orbital approximation; no cubes were calculated.\n";
	summary << fixed;
	summary << "HOMO_index " << homo << "\n";
	summary << "HOMO_energy " << setprecision(6) << wavy.get_MO_energy(homo) << "\n";
	summary << "LUMO_index " << lumo << "\n";
	summary << "LUMO_energy " << setprecision(6) << wavy.get_MO_energy(lumo) << "\n";
	summary << "unrestricted " << (unrestricted ? 1 : 0) << "\n";
	static const int ord[5] = {
		PartitionResults::CHARGE_ORDER::S_HIRSH,
		PartitionResults::CHARGE_ORDER::S_BECKE,
		PartitionResults::CHARGE_ORDER::S_TFVC,
		PartitionResults::CHARGE_ORDER::S_MBIS,
		PartitionResults::CHARGE_ORDER::S_EMBIS
	};
	static const char *ord_name[5] = { "hirshfeld", "becke", "tfvc", "mbis", "embis" };
	summary << "# condensed Fukui: atom f+ f- df, one block per partition\n";
	for (int s = 0; s < 5; s++)
	{
		summary << "partition " << ord_name[s] << "\n";
		for (size_t a = 0; a < n; a++)
		{
			const double fp = condensed.f_plus[ord[s]][a];
			const double fm = condensed.f_minus[ord[s]][a];
			summary << "  " << condensed.labels[a]
					<< " " << setprecision(6) << fp
					<< " " << setprecision(6) << fm
					<< " " << setprecision(6) << (fp - fm) << "\n";
		}
	}
	summary.close();
	log2 << "\nWrote " << summary_path.filename().string() << endl;
}

namespace {

struct PromolecularFragmentDensities {
	double sum = 0.0;      // over all fragments
	double dominant = 0.0; // largest single-fragment contribution

	double total() const
	{
		return sum;
	}
};

struct PromolecularAtom {
	d3 pos;
	int charge = 1;
	int fragment = 0;
};

// nearby: indices into atoms, ascending; atoms left out must be past their spline_reach() (exactly 0).
// reach2[a]: spline_reach()^2 of atom a, past which it is skipped here too
PromolecularFragmentDensities promolecular_fragment_densities_at(
	const d3 &pos,
	const std::vector<PromolecularAtom> &atoms,
	const std::vector<int> &nearby,
	const std::vector<Thakkar> &atom_models,
	const vec &reach2)
{
	PromolecularFragmentDensities result;
	// atoms are grouped by fragment (add_promolecular_atoms appends whole fragments),
	// so a running per-fragment sum needs no per-fragment storage
	double fragment_sum = 0.0;
	int current_fragment = nearby.empty() ? 0 : atoms[nearby.front()].fragment;
	for (const int a : nearby)
	{
		const PromolecularAtom &atom = atoms[a];
		if (atom.fragment != current_fragment)
		{
			result.dominant = std::max(result.dominant, fragment_sum);
			fragment_sum = 0.0;
			current_fragment = atom.fragment;
		}
		// the block list holds atoms that reach some point of the block; past spline_reach() this one adds exactly 0
		const double dx = pos[0] - atom.pos[0], dy = pos[1] - atom.pos[1], dz = pos[2] - atom.pos[2];
		const double d2 = dx * dx + dy * dy + dz * dz;
		if (d2 > reach2[a])
			continue;
		// Table lookup for the mask pass only; lambda2 and the RDG come from the
		// analytic Thakkar derivatives in promolecular_derivatives_at(); the exact
		// Slater sums here cost several times the run for a handful of kept points
		const double contribution = atom_models[atom.charge - 1].get_interpolated_density_spline(std::sqrt(d2));
		fragment_sum += contribution;
		result.sum += contribution;
	}
	result.dominant = std::max(result.dominant, fragment_sum);
	return result;
}

bool is_promolecular_nci_point(
	const PromolecularFragmentDensities &densities,
	double total_density,
	double dominant_density_cutoff,
	double fragment_sum_cutoff)
{
	total_density = std::abs(total_density);
	if (total_density <= 1E-20)
		return false;

	const double fragment_density = densities.total();
	if (fragment_density <= 1E-20)
		return false;

	if (densities.dominant >= fragment_density * dominant_density_cutoff)
		return false;

	return fragment_density >= total_density * fragment_sum_cutoff;
}

// rho, rho', rho'' of one Thakkar atom on [r0, r0 + n h) as quintic Hermite pieces through the exact rho, rho', rho''
// at the nodes; rho' and rho'' are the piece's own derivatives
struct RadialQuintic
{
	double r0 = 0.0, inv_h = 0.0;
	std::vector<std::array<double, 6>> c; // per piece, in t = (r - r_i) / h

	bool eval(const double r, double &rho, double &d1, double &d2) const
	{
		const double s = (r - r0) * inv_h;
		if (!(s >= 0.0) || s >= (double)c.size())
			return false;
		const int i = (int)s;
		const double t = s - i;
		const std::array<double, 6> &k = c[i];
		if (std::isnan(k[0]))
			return false;
		rho =k[0] + t * (k[1] + t * (k[2] + t * (k[3] + t * (k[4] + t * k[5]))));
		d1 = (k[1] + t * (2 * k[2] + t * (3 * k[3] + t * (4 * k[4] + t * 5 * k[5])))) * inv_h;
		d2 = (2 * k[2] + t * (6 * k[3] + t * (12 * k[4] + t * 20 * k[5]))) * inv_h * inv_h;
		return true;
	}
};

// The table on [r0, r1). A piece that strays more than tol * env(r) from the exact rho, rho' or rho'' at t = 1/4, 1/2, 3/4
// (env: the sup envelope of the RDG pass, sampled every env_dr) is NaN, and atoms there take the Slater sums: Thakkar
// drops each exponent's term past z r = 46.5, a step no smooth piece follows. h is the step with the fewest NaN pieces;
// coarse steps fail on truncation, fine ones on rounding, since rho'' comes from differences of rho over h^2
RadialQuintic make_radial_quintic(const Thakkar &model, const double r0, const double r1, const vec &env, const double env_dr, const double tol)
{
	RadialQuintic best;
	double best_bad = 2.0;
	for (double h = 0.05; h >= 1E-4; h /= 2)
	{
		RadialQuintic q;
		q.r0 = r0;
		q.inv_h = 1.0 / h;
		const int n = (int)std::ceil((r1 - r0) / h);
		q.c.resize(n);
		double d1, d2;
		double f0 = model.get_radial_density(r0, d1, d2), m0 = h * d1, a0 = h * h * d2;
		int bad = 0;
		for (int i = 0; i < n; i++)
		{
			const double f1 = model.get_radial_density(r0 + (i + 1) * h, d1, d2), m1 = h * d1, a1 = h * h * d2;
			const double A = f1 - f0 - m0 - a0 / 2, B = m1 - m0 - a0, C = a1 - a0;
			q.c[i] = { f0, m0, a0 / 2, 10 * A - 4 * B + C / 2, -15 * A + 7 * B - C, 6 * A - 3 * B + C / 2 };
			f0 = f1, m0 = m1, a0 = a1;
			for (const double t : { 0.25, 0.5, 0.75 })
			{
				const double r = r0 + (i + t) * h, bound = tol * env[std::min((size_t)(r / env_dr), env.size() - 1)];
				double e1, e2, p, p1, p2;
				const double e = model.get_radial_density(r, e1, e2);
				q.eval(r, p, p1, p2);
				if (!(std::abs(p - e) <= bound && std::abs(p1 - e1) <= bound && std::abs(p2 - e2) <= bound))
				{
					q.c[i][0] = std::numeric_limits<double>::quiet_NaN();
					bad++;
					break;
				}
			}
		}
		const double frac = (double)bad / n;
		if (frac > best_bad)
			break; // past the bottom
		if (frac < best_bad)
			best = std::move(q), best_bad = frac;
		if (bad == 0)
			break;
	}
	return best;
}

// Atoms past from[slot] (per element, this block's far radius) are summed from their table instead of the Slater sums
struct FarAtoms
{
	const double *from;
	const std::vector<int> &slot;
	const std::vector<RadialQuintic> &tables;
};

// rho, grad rho and the Hessian of the promolecule at pos, each Thakkar atom's analytic rho', rho''
// summed through Centred<Thakkar>; replaces the finite-difference stencils on the rho cube, whose
// error scaled with the grid step and whose edge points were clamped. nearby: indices into atoms, ascending
double promolecular_derivatives_at(
	const d3 &pos,
	const std::vector<PromolecularAtom> &atoms,
	const std::vector<int> &nearby,
	const std::vector<Thakkar> &atom_models,
	d3 &grad,
	double *hessian,
	const FarAtoms *distant = nullptr)
{
	double rho = 0.0;
	grad = { 0.0, 0.0, 0.0 };
	std::fill(hessian, hessian + 9, 0.0);
	for (const int a : nearby)
	{
		const PromolecularAtom &atom = atoms[a];
		d3 g;
		double H[9];
		if (distant)
		{
			// calculate_hessian(Centred<A>) regrouped: grad = rho'/r d, H = (rho'' - rho'/r) d d / r^2 + rho'/r delta
			const int s = distant->slot[atom.charge - 1];
			const d3 d{ pos[0] - atom.pos[0], pos[1] - atom.pos[1], pos[2] - atom.pos[2] };
			const double r = array_length(d);
			double f, d1, d2;
			if (r >= distant->from[s] && r > 1E-4 && distant->tables[s].eval(r, f, d1, d2))
			{
				const double inv_r = 1.0 / r, a1 = d1 * inv_r, b = (d2 - a1) * inv_r * inv_r;
				rho += f;
				for (int i = 0; i < 3; i++)
				{
					grad[i] += a1 * d[i];
					for (int j = 0; j < 3; j++)
						hessian[3 * i + j] += b * d[i] * d[j] + (i == j ? a1 : 0.0);
				}
				continue;
			}
		}
		const Centred<Thakkar> source{ atom_models[atom.charge - 1], atom.pos };
		rho += calculate_hessian(source, pos, g, H);
		for (int k = 0; k < 3; k++)
			grad[k] += g[k];
		for (int k = 0; k < 9; k++)
			hessian[k] += H[k];
	}
	return rho;
}

double reduced_density_gradient(const double rho, const d3 &grad)
{
	if (std::abs(rho) <= 1E-20)
		return 0.0;
	const double rdg_factor = 2.0 * std::pow(3.0 * constants::PI * constants::PI, 1.0 / 3.0);
	return array_length(grad) / (rdg_factor * std::pow(std::abs(rho), 4.0 / 3.0));
}

void add_promolecular_atoms(
	const WFN &fragment,
	int fragment_id,
	std::vector<PromolecularAtom> &atoms)
{
	for (int i = 0; i < fragment.get_ncen(); i++)
	{
		PromolecularAtom atom_entry;
		atom_entry.pos = fragment.get_atom_pos(i);
		atom_entry.charge = fragment.get_atom_charge(i);
		atom_entry.fragment = fragment_id;
		atoms.push_back(atom_entry);
	}
}

void add_atoms_to_combined_wfn(const WFN &fragment, WFN &combined)
{
	for (int i = 0; i < fragment.get_ncen(); i++)
		combined.push_back_atom(
			fragment.get_atom_label(i),
			fragment.get_atom_coordinate(i, 0),
			fragment.get_atom_coordinate(i, 1),
			fragment.get_atom_coordinate(i, 2),
			fragment.get_atom_charge(i));
}

std::string tcl_quote_path(const std::filesystem::path &path)
{
	const std::string input = path.generic_string();
	std::string output = "\"";
	for (const char c : input)
	{
		if (c == '\\' || c == '"' || c == '$' || c == '[' || c == ']')
			output += '\\';
		output += c;
	}
	output += '"';
	return output;
}

void write_promolecular_nci_vmd(
	const pathvec &xyz_files,
	const std::filesystem::path &output_base,
	const properties_options &opts,
	std::ostream &log)
{
	const std::filesystem::path signed_rho_path = output_base.string() + "_signed_rho.cube";
	const std::filesystem::path rdg_path = output_base.string() + "_rdg.cube";
	const std::filesystem::path vmd_path = output_base.string() + "_nci.vmd";

	std::ofstream vmd_file(vmd_path, std::ios::out);
	err_checkf(vmd_file.good(), "Could not open " + vmd_path.string() + " for writing.", log);

	vmd_file
		<< "# VMD visualization for promolecular NCI/RDG analysis\n"
		<< "# Load with: vmd -e " << vmd_path.filename().string() << "\n"
		<< "display projection Orthographic\n"
		<< "display depthcue off\n"
		<< "axes location Off\n"
		<< "color Display Background white\n"
		<< "color scale method BGR\n"
		<< "\n";
	for (size_t f = 0; f < xyz_files.size(); f++)
		vmd_file
			<< "mol new " << tcl_quote_path(std::filesystem::absolute(xyz_files[f])) << " type xyz waitfor all\n"
			<< "set frag" << f + 1 << " [molinfo top]\n"
			<< "mol delrep 0 $frag" << f + 1 << "\n"
			<< "mol representation CPK\n"
			<< "mol color Name\n"
			<< "mol selection all\n"
			<< "mol material Opaque\n"
			<< "mol addrep $frag" << f + 1 << "\n"
			<< "\n";
	vmd_file
		<< "mol new " << tcl_quote_path(std::filesystem::absolute(signed_rho_path)) << " type cube first 0 last -1 step 1 waitfor 1 volsets {0 }\n"
		<< "set nci [molinfo top]\n"
		<< "mol delrep 0 $nci\n"
		<< "mol addfile " << tcl_quote_path(std::filesystem::absolute(rdg_path)) << " type cube first 0 last -1 step 1 waitfor 1 volsets {0 } $nci\n"
		<< "while {[molinfo $nci get numreps] > 0} {\n"
		<< "    mol delrep 0 $nci\n"
		<< "}\n"
		<< "mol addrep $nci\n"
		<< "mol modselect 0 $nci all\n"
		<< "mol modstyle 0 $nci Isosurface 0.5 1 0 0 1 1\n"
		<< "mol modcolor 0 $nci Volume 0\n"
		<< "mol modmaterial 0 $nci Opaque\n"
		<< "mol colupdate 0 $nci on\n"
		<< "mol scaleminmax $nci 0 " << -opts.promol_nci_colour_max << " " << opts.promol_nci_colour_max << "\n"
		<< "color scale method BGR\n"
		<< "color scale midpoint 0.5\n"
		<< "\n"
		<< "display resetview\n";

	log << "Wrote " << vmd_path << std::endl;
}

void write_promolecular_nci_plot_script(
	const std::filesystem::path &output_base,
	const properties_options &opts,
	std::ostream &log)
{
	const std::filesystem::path values_path = output_base.string() + "_values.dat";
	const std::filesystem::path plot_path = output_base.string() + "_plot.py";

	std::ofstream plot_file(plot_path, std::ios::out);
	err_checkf(plot_file.good(), "Could not open " + plot_path.string() + " for writing.", log);

	plot_file
		<< "from pathlib import Path\n"
		<< "import sys\n"
		<< "\n"
		<< "import matplotlib.pyplot as plt\n"
		<< "import numpy as np\n"
		<< "from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm\n"
		<< "\n"
		<< "def read_float_arg(name, default):\n"
		<< "    if name not in sys.argv[1:]:\n"
		<< "        return default\n"
		<< "    idx = sys.argv.index(name)\n"
		<< "    if idx + 1 >= len(sys.argv):\n"
		<< "        raise SystemExit(f\"{name} requires a numeric value\")\n"
		<< "    return float(sys.argv[idx + 1])\n"
		<< "\n"
		<< "dat_path = Path(__file__).with_name(" << std::quoted(values_path.filename().string()) << ")\n"
		<< "data = np.loadtxt(dat_path, comments=\"#\")\n"
		<< "if data.ndim == 1:\n"
		<< "    data = data.reshape(1, -1)\n"
		<< "\n"
		<< "signed_rho = data[:, 0]\n"
		<< "rdg = data[:, 1]\n"
		<< "xmax = " << std::setprecision(16) << opts.promol_nci_rho_abs_max << "\n"
		<< "ymax = " << std::setprecision(16) << opts.promol_nci_rdg_max << "\n"
		<< "xmax = read_float_arg(\"-xmax\", None if xmax < 0.0 else xmax)\n"
		<< "ymax = read_float_arg(\"-ymax\", None if ymax < 0.0 else ymax)\n"
		<< "mask = np.ones_like(signed_rho, dtype=bool)\n"
		<< "if xmax is not None:\n"
		<< "    mask &= np.abs(signed_rho) <= xmax\n"
		<< "if ymax is not None:\n"
		<< "    mask &= rdg <= ymax\n"
		<< "signed_rho = signed_rho[mask]\n"
		<< "rdg = rdg[mask]\n"
		<< "if signed_rho.size == 0:\n"
		<< "    raise SystemExit(\"No points remain after applying plot limits\")\n"
		<< "rho_min = float(np.min(signed_rho))\n"
		<< "rho_max = float(np.max(signed_rho))\n"
		<< "\n"
		<< "cmap = LinearSegmentedColormap.from_list(\"nci_bgr\", [\"blue\", \"green\", \"red\"])\n"
		<< "if rho_min < 0.0 < rho_max:\n"
		<< "    norm = TwoSlopeNorm(vmin=rho_min, vcenter=0.0, vmax=rho_max)\n"
		<< "else:\n"
		<< "    norm = None\n"
		<< "\n"
		<< "fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)\n"
		<< "scatter = ax.scatter(signed_rho, rdg, c=signed_rho, s=4, cmap=cmap, norm=norm, linewidths=0)\n"
		<< "ax.axvline(0.0, color=\"0.65\", linewidth=0.8)\n"
		<< "ax.set_xlabel(\"signed rho\")\n"
		<< "ax.set_ylabel(\"RDG\")\n"
		<< "if xmax is not None:\n"
		<< "    ax.set_xlim(-xmax, xmax)\n"
		<< "if ymax is not None:\n"
		<< "    ax.set_ylim(0.0, ymax)\n"
		<< "ax.set_title(dat_path.stem.replace(\"_values\", \"\"))\n"
		<< "cbar = fig.colorbar(scatter, ax=ax)\n"
		<< "cbar.set_label(\"signed rho\")\n"
		<< "png_path = dat_path.with_name(dat_path.stem.replace(\"_values\", \"_plot\") + \".png\")\n"
		<< "fig.savefig(png_path, dpi=300)\n"
		<< "print(f\"Wrote {png_path}\")\n"
		<< "if \"-show\" in sys.argv[1:]:\n"
		<< "    plt.show()\n"
		<< "else:\n"
		<< "    plt.close(fig)\n";

	log << "Wrote " << plot_path << std::endl;
}

} // namespace

void promolecular_nci_analysis(
	const pathvec &xyz_files,
	const properties_options &opts,
	std::ostream &log,
	const std::filesystem::path &cif)
{
	using namespace std;
	const _time_point t_start = get_time();

	err_checkf(xyz_files.size() >= 2, "Promolecular NCI needs at least two XYZ fragments.", log);
	citations::cite(citations::Method::NCI, log);

	// Output names join every fragment stem: a_b_c_values.dat etc.
	std::string joined_stems = xyz_files.front().stem().string();
	std::string joined_names = xyz_files.front().string();
	for (size_t f = 1; f < xyz_files.size(); f++)
	{
		joined_stems += "_" + xyz_files[f].stem().string();
		joined_names += (f + 1 == xyz_files.size() ? " and " : ", ") + xyz_files[f].string();
	}

	WFN combined(e_origin::xyz);
	combined.set_path(xyz_files.front().parent_path() / (joined_stems + ".xyz"));
	vector<PromolecularAtom> atoms;
	vector<std::array<double, 6>> fragment_boxes; // each fragment's atoms + radius, bohr
	for (size_t f = 0; f < xyz_files.size(); f++)
	{
		err_checkf(std::filesystem::exists(xyz_files[f]), "XYZ file does not exist: " + xyz_files[f].string(), log);
		WFN fragment(e_origin::xyz);
		fragment.read_xyz(xyz_files[f], log, false);
		properties_options fragment_opts = opts;
		readxyzMinMax_fromWFN(fragment, fragment_opts);
		fragment_boxes.push_back(fragment_opts.MinMax);
		add_atoms_to_combined_wfn(fragment, combined);
		add_promolecular_atoms(fragment, static_cast<int>(f) + 1, atoms);
	}

	properties_options local_opts = opts;
	// grid: a box around the fragments, or the unit cell when a cif is given (Olex2's xgrid spans the cell)
	vec2 cell_matrix(3, vec(3, 0.0));
	if (cif.empty())
	{
		readxyzMinMax_fromWFN(combined, local_opts);
		// an intermolecular point lies within radius of two fragments, so only the pairwise overlaps of the fragment boxes can hold one
		std::array<double, 6> overlaps = { 1E300, 1E300, 1E300, -1E300, -1E300, -1E300 };
		for (size_t f = 0; f < fragment_boxes.size(); f++)
			for (size_t g = f + 1; g < fragment_boxes.size(); g++)
			{
				std::array<double, 6> o;
				bool empty = false;
				for (int i = 0; i < 3; i++)
				{
					o[i] = std::max(fragment_boxes[f][i], fragment_boxes[g][i]);
					o[i + 3] = std::min(fragment_boxes[f][i + 3], fragment_boxes[g][i + 3]);
					empty |= o[i] >= o[i + 3];
				}
				if (!empty)
					for (int i = 0; i < 3; i++)
					{
						overlaps[i] = std::min(overlaps[i], o[i]);
						overlaps[i + 3] = std::max(overlaps[i + 3], o[i + 3]);
					}
			}
		err_checkf(overlaps[0] < overlaps[3], "The fragments are more than two -radius apart; no intermolecular region.", log);
		for (int i = 0; i < 3; i++)
		{
			local_opts.MinMax[i] = std::max(local_opts.MinMax[i], overlaps[i]);
			local_opts.MinMax[i + 3] = std::min(local_opts.MinMax[i + 3], overlaps[i + 3]);
			local_opts.NbSteps[i] = (int)ceil(constants::bohr2ang(local_opts.MinMax[i + 3] - local_opts.MinMax[i]) / local_opts.resolution);
			local_opts.NbSteps[i] += local_opts.NbSteps[i] % 2; // even, as readxyzMinMax_fromWFN
		}
		for (int i = 0; i < 3; i++)
			cell_matrix[i][i] = (local_opts.MinMax[i + 3] - local_opts.MinMax[i]) / local_opts.NbSteps[i];
	}
	else
		readxyzMinMax_fromCIF(cif, local_opts, cell_matrix);
	err_checkf(local_opts.NbSteps[0] > 1 && local_opts.NbSteps[1] > 1 && local_opts.NbSteps[2] > 1,
		"Promolecular NCI grid is too small; decrease -resolution or increase -radius.", log);

	const std::filesystem::path output_base = xyz_files.front().parent_path() / joined_stems;

	const vector<Thakkar> atom_models = make_thakkar_interpolators();

	cube rho_cube(local_opts.NbSteps, combined.get_ncen(), true);
	cube signed_rho_cube(local_opts.NbSteps, combined.get_ncen(), true);
	cube rdg_cube(local_opts.NbSteps, combined.get_ncen(), true);
	rho_cube.give_parent_wfn(combined);
	signed_rho_cube.give_parent_wfn(combined);
	rdg_cube.give_parent_wfn(combined);

	for (int i = 0; i < 3; i++)
	{
		rho_cube.set_origin(i, local_opts.MinMax[i]);
		signed_rho_cube.set_origin(i, local_opts.MinMax[i]);
		rdg_cube.set_origin(i, local_opts.MinMax[i]);
		for (int j = 0; j < 3; j++)
		{
			rho_cube.set_vector(i, j, cell_matrix[i][j]);
			signed_rho_cube.set_vector(i, j, cell_matrix[i][j]);
			rdg_cube.set_vector(i, j, cell_matrix[i][j]);
		}
	}
	rho_cube.calc_dv();
	signed_rho_cube.calc_dv();
	rdg_cube.calc_dv();

	rho_cube.set_comment1("Promolecular density using Thakkar spherical atoms");
	rho_cube.set_comment2("from " + joined_names);
	signed_rho_cube.set_comment1("Promolecular signed density using Thakkar spherical atoms");
	signed_rho_cube.set_comment2("from " + joined_names);
	rdg_cube.set_comment1("Promolecular reduced density gradient using Thakkar spherical atoms");
	rdg_cube.set_comment2("from " + joined_names);

	log << "Promolecular NCI analysis for " << xyz_files.size() << " fragments: " << joined_names << endl;
	log << "Grid points: " << local_opts.n_grid_points() << endl;
	log << "Dominant-fragment density discard cutoff: " << opts.promol_nci_rcut1 << endl;
	log << "Fragment-sum density keep cutoff: " << opts.promol_nci_rcut2 << endl;

	const _time_point t_mask = get_time();
	//Per 8^3 block only the atoms within spline_reach(): the rest add exactly 0, in order, so bit-identical
	constexpr int B = 8;
	vector<d3> apos(atoms.size());
	vec reach(atoms.size()), reach2(atoms.size());
	vector<int> all_atoms(atoms.size());
	for (size_t a = 0; a < atoms.size(); a++)
	{
		apos[a] = atoms[a].pos;
		reach[a] = atom_models[atoms[a].charge - 1].spline_reach();
		reach2[a] = reach[a] * reach[a];
		all_atoms[a] = (int)a;
	}
	const vector<vector<int>> block_atoms = block_atom_lists(rho_cube, B, apos, [&](int, int a) { return reach[a]; });
	const int nby = (rho_cube.get_size(1) + B - 1) / B, nbz = (rho_cube.get_size(2) + B - 1) / B;
	ProgressBar density_progress(rho_cube.get_size(0), 50, "=", " ", "Calculating promolecular rho");
#pragma omp parallel for schedule(dynamic)
	for (int x = 0; x < rho_cube.get_size(0); x++)
	{
		for (int y = 0; y < rho_cube.get_size(1); y++)
			for (int z = 0; z < rho_cube.get_size(2); z++)
			{
				const d3 pos = rho_cube.get_pos(x, y, z);
				const vector<int> &nearby = block_atoms[((size_t)(x / B) * nby + y / B) * nbz + z / B];
				const PromolecularFragmentDensities densities = promolecular_fragment_densities_at(pos, atoms, nearby, atom_models, reach2);
				const double rho = densities.total();
				rho_cube.set_value(x, y, z, rho);
				rdg_cube.set_value(
					x,
					y,
					z,
					is_promolecular_nci_point(densities, rho, opts.promol_nci_rcut1, opts.promol_nci_rcut2) ? 1.0 : 101.0);
			}
		density_progress.update();
	}
	std::cout << std::endl;

	i3 shrink_lower;
	i3 shrink_upper;
	if (rdg_cube.find_value_bounds(shrink_lower, shrink_upper, 101.0))
	{
		rho_cube.shrink_to_bounds(shrink_lower, shrink_upper);
		signed_rho_cube.shrink_to_bounds(shrink_lower, shrink_upper);
		rdg_cube.shrink_to_bounds(shrink_lower, shrink_upper);
		log << "Shrunk promolecular NCI work grid to "
			<< rdg_cube.get_size(0) << " x "
			<< rdg_cube.get_size(1) << " x "
			<< rdg_cube.get_size(2)
			<< " grid points, ignoring RDG mask value 101." << endl;
	}

	const _time_point t_rdg = get_time();
	ofstream values_file(output_base.string() + "_values.dat", ios::out);
	err_checkf(values_file.good(), "Could not open " + output_base.string() + "_values.dat for writing.", log);
	values_file << "# signed_rho rdg\n";
	values_file << "# rcut1 " << opts.promol_nci_rcut1 << " rcut2 " << opts.promol_nci_rcut2 << "\n";

	unsigned long long kept_points = 0;
#ifdef _OPENMP
	// schedule(dynamic) below makes the mapping of grid points to threads (and
	// therefore the row order in values_by_thread) non-deterministic between
	// runs; -promol_nci_single_thread pins this region to one thread so output
	// ordering is reproducible, e.g. for golden-file test generation.
	const int nci_write_threads = opts.promol_nci_single_threaded ? 1 : omp_get_max_threads();
#else
	const int nci_write_threads = 1;
#endif
	std::vector<std::ostringstream> values_by_thread(nci_write_threads);
	// masked points next to a kept one get their RDG too, so the mesh below ends a sheet in a face whose centre fails the
	// mask (dropped) instead of a wall along the mask's voxel staircase. Other masked points stay at 101, above any iso,
	// and the cubes are masked again after the mesh
	const int nx = rho_cube.get_size(0), ny = rho_cube.get_size(1), nz = rho_cube.get_size(2);
	const auto at = [&](int x, int y, int z) { return (static_cast<size_t>(x) * ny + y) * nz + z; };
	std::vector<char> masked(static_cast<size_t>(nx) * ny * nz), rim(masked.size(), 0);
	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++)
				masked[at(x, y, z)] = std::abs(rdg_cube.get_value(x, y, z) - 101.0) <= 1E-12;
	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++)
				if (!masked[at(x, y, z)])
					for (int i = std::max(x - 1, 0); i <= std::min(x + 1, nx - 1); i++)
						for (int j = std::max(y - 1, 0); j <= std::min(y + 1, ny - 1); j++)
							for (int k = std::max(z - 1, 0); k <= std::min(z + 1, nz - 1); k++)
								rim[at(i, j, k)] = masked[at(i, j, k)];
	// The analytic sums below, per 8^3 block plus the next grid layer (so a marching-cubes cell lies in the block of its
	// low corner), only over the atoms that can add 1E-14 rho_min / (4 N) there. env(r) = 2 sup over r' >= r of
	// |rho| + |rho'| + |rho''| + |rho'| / r', sampled every 0.01 bohr, bounds every term of rho, grad and Hessian, so the
	// atoms left out move each sum by < 1E-14 rho; rho_min is the block's smallest mask-pass spline rho (<= 0 keeps all)
	constexpr double env_dr = 0.01, inf = std::numeric_limits<double>::infinity();
	std::vector<int> slot(atom_models.size(), -1);
	std::vector<vec> envelope;
	for (const PromolecularAtom &atom : atoms)
	{
		if (slot[atom.charge - 1] >= 0)
			continue;
		slot[atom.charge - 1] = (int)envelope.size();
		vec &e = envelope.emplace_back();
		for (int k = 0; k < 40000 && (e.empty() || e.back() >= 1E-300); k++)
		{
			const double r = k * env_dr;
			double d1, d2;
			const double rho = atom_models[atom.charge - 1].get_radial_density(r, d1, d2);
			e.push_back(2 * (std::abs(rho) + std::abs(d1) + std::abs(d2) + (k > 0 ? std::abs(d1) / r : 0.0)));
		}
		for (int k = (int)e.size() - 2; k >= 0; k--)
			e[k] = std::max(e[k], e[k + 1]);
	}
	const int sby = (ny + B - 1) / B, sbz = (nz + B - 1) / B;
	const auto block_of = [&](int x, int y, int z) { return ((size_t)(x / B) * sby + y / B) * sbz + z / B; };
	vec rho_min((size_t)((nx + B - 1) / B) * sby * sbz, inf); // inf: nothing in the block is evaluated
#pragma omp parallel for schedule(dynamic)
	for (int b = 0; b < (int)rho_min.size(); b++)
	{
		const int bx = b / (sby * sbz) * B, by = b / sbz % sby * B, bz = b % sbz * B;
		for (int x = bx; x < std::min(bx + B + 1, nx); x++)
			for (int y = by; y < std::min(by + B + 1, ny); y++)
				for (int z = bz; z < std::min(bz + B + 1, nz); z++)
					if (!masked[at(x, y, z)] || rim[at(x, y, z)])
						rho_min[b] = std::min(rho_min[b], rho_cube.get_value(x, y, z));
	}
	const size_t nel = envelope.size();
	vec cut(rho_min.size() * nel);
	for (size_t b = 0; b < rho_min.size(); b++)
		for (size_t s = 0; s < nel; s++)
		{
			const double thr = 1E-14 * rho_min[b] / (4.0 * atoms.size());
			const auto it = std::partition_point(envelope[s].begin(), envelope[s].end(), [thr](double v) { return v >= thr; });
			cut[b * nel + s] = rho_min[b] == inf ? -inf : it == envelope[s].end() ? inf : (it - envelope[s].begin()) * env_dr;
		}
	// Past far_from (env < T rho_min, about 86 % of the terms on sucrose) an atom comes from its quintic table, which stays
	// within 1E-13 / (T N) env of the Slater sums: all N atoms together move each sum by < 1E-13 rho_min
	constexpr double T = 1E-7;
	vec far_from(cut.size(), inf);
	vector<RadialQuintic> tables(nel);
	for (size_t z = 0; z < slot.size(); z++)
	{
		const int s = slot[z];
		if (s < 0)
			continue;
		double lo = inf, hi = 0;
		for (size_t b = 0; b < rho_min.size(); b++)
		{
			if (rho_min[b] == inf)
				continue;
			const double thr = T * rho_min[b];
			const auto it = std::partition_point(envelope[s].begin(), envelope[s].end(), [thr](double v) { return v >= thr; });
			if (it != envelope[s].end())
				lo = std::min(lo, far_from[b * nel + s] = (it - envelope[s].begin()) * env_dr);
			if (cut[b * nel + s] != inf)
				hi = std::max(hi, cut[b * nel + s]);
		}
		if (lo < hi)
			tables[s] = make_radial_quintic(atom_models[z], lo, hi + 2.0, envelope[s], env_dr, 1E-13 / (T * atoms.size()));
	}
	const vector<vector<int>> rdg_atoms = block_atom_lists(rho_cube, B, apos, [&](int b, int a) { return cut[b * nel + slot[atoms[a].charge - 1]]; }, 1);
	const vector<vector<int>> face_atoms = block_atom_lists(rho_cube, B, apos, [&](int b, int a) { return rho_min[b] == inf ? -inf : reach[a]; }, 1);
	ProgressBar rdg_progress(rho_cube.get_size(0), 50, "=", " ", "Calculating promolecular RDG");
#pragma omp parallel reduction(+ : kept_points) num_threads(nci_write_threads)
	{
		int thread_id = 0;
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#endif
		std::ostringstream &local_values = values_by_thread[thread_id];
#pragma omp for schedule(dynamic)
	for (int x = 0; x < rho_cube.get_size(0); x++)
	{
		for (int y = 0; y < rho_cube.get_size(1); y++)
		{
			for (int z = 0; z < rho_cube.get_size(2); z++)
			{
				const bool is_masked = masked[at(x, y, z)];
				if (is_masked && !rim[at(x, y, z)])
				{
					signed_rho_cube.set_value(x, y, z, 0.0); // rdg is still the mask's 101
					continue;
				}

				d3 grad;
				double hessian[9];
				const size_t b = block_of(x, y, z);
				const FarAtoms distant{ &far_from[b * nel], slot, tables };
				const double rho = promolecular_derivatives_at(rho_cube.get_pos(x, y, z), atoms, rdg_atoms[b], atom_models, grad, hessian, &distant);
				const double lambda2 = get_lambda_1(hessian);
				const double signed_rho = lambda2 < 0.0 ? -rho : rho;
				const double rdg = sanitize_finite(reduced_density_gradient(rho, grad));

				signed_rho_cube.set_value(x, y, z, signed_rho);
				rdg_cube.set_value(x, y, z, rdg);
				if (is_masked)
					continue;

				if (opts.promol_nci_rho_abs_max >= 0.0 && std::abs(signed_rho) > opts.promol_nci_rho_abs_max)
					continue;
				if (opts.promol_nci_rdg_max >= 0.0 && rdg > opts.promol_nci_rdg_max)
					continue;

				local_values << scientific << setprecision(10) << signed_rho << " " << rdg << "\n";
				kept_points++;
			}
		}
		rdg_progress.update();
	}
	}
	for (const std::ostringstream &local_values : values_by_thread)
		values_file << local_values.str();
	const _time_point t_mesh = get_time();

	// marching cubes on the unmasked RDG; a face stays when its centre is intermolecular, coloured and written to the
	// .dat with sign(lambda2) rho there (VMD's BGR, +-colour_max)
	const std::vector<Triangle> triangles = marchingCubes(rdg_cube, opts.promol_nci_iso);
	std::vector<char> keep(triangles.size(), 0);
	vec centre_signed_rho(triangles.size());
	const vec2 to_index = crystal_energies::inverse3(cell_matrix);
#pragma omp parallel for schedule(dynamic, 256)
	for (long long i = 0; i < static_cast<long long>(triangles.size()); i++)
	{
		const d3 c = triangles[i].calc_center();
		// the cell holding the centre, by its low corner; a block with nothing evaluated (a centre rounded onto a face) keeps all atoms
		int idx[3];
		for (int d = 0; d < 3; d++)
		{
			double f = 0;
			for (int e = 0; e < 3; e++)
				f += to_index[d][e] * (c[e] - rho_cube.get_origin(e));
			idx[d] = std::clamp((int)std::floor(f), 0, rho_cube.get_size(d) - 1);
		}
		const size_t b = block_of(idx[0], idx[1], idx[2]);
		const PromolecularFragmentDensities densities = promolecular_fragment_densities_at(c, atoms, rho_min[b] == inf ? all_atoms : face_atoms[b], atom_models, reach2);
		if (!is_promolecular_nci_point(densities, densities.total(), opts.promol_nci_rcut1, opts.promol_nci_rcut2))
			continue;
		d3 grad;
		double hessian[9];
		const FarAtoms distant{ &far_from[b * nel], slot, tables };
		const double rho = promolecular_derivatives_at(c, atoms, rho_min[b] == inf ? all_atoms : rdg_atoms[b], atom_models, grad, hessian, rho_min[b] == inf ? nullptr : &distant);
		centre_signed_rho[i] = get_lambda_1(hessian) < 0.0 ? -rho : rho;
		keep[i] = 1;
	}
	std::vector<Triangle> faces;
	vec face_signed_rho;
	for (size_t i = 0; i < triangles.size(); i++)
	{
		if (!keep[i])
			continue;
		const Triangle &t = triangles[i];
		face_signed_rho.push_back(centre_signed_rho[i]);
		faces.emplace_back(t.get_v(1), t.get_v(3), t.get_v(2), // low RDG is inside: wind as the density surfaces
			mix_colour(centre_signed_rho[i], { { { 0, 0, 255 }, { 0, 255, 0 }, { 255, 0, 0 } } },
				-opts.promol_nci_colour_max, opts.promol_nci_colour_max));
	}
	writeColourObj(output_base.string() + "_nci.obj", faces);
	{
		ofstream face_file(output_base.string() + "_nci.dat");
		face_file << "# signed_rho\n" << scientific << setprecision(6);
		for (const double s : face_signed_rho)
			face_file << s << "\n";
	}
	log << "Intermolecular NCI surface: " << faces.size() << " of " << triangles.size() << " RDG = "
		<< opts.promol_nci_iso << " faces, " << output_base.string() << "_nci.obj" << endl;

	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++)
				if (rim[at(x, y, z)])
				{
					signed_rho_cube.set_value(x, y, z, 0.0);
					rdg_cube.set_value(x, y, z, 101.0);
				}
	const _time_point t_write = get_time();

	signed_rho_cube.set_path(output_base.string() + "_signed_rho.cube");
	rdg_cube.set_path(output_base.string() + "_rdg.cube");
	signed_rho_cube.write_file(true);
	rdg_cube.write_file(true);
	write_promolecular_nci_vmd(xyz_files, output_base, opts, log);
	write_promolecular_nci_plot_script(output_base, opts, log);
	if (!constants::hide_timings)
		log << "Promolecular NCI setup: " << get_msec(t_start, t_mask) << " ms, mask pass: " << get_msec(t_mask, t_rdg)
			<< " ms, RDG pass: " << get_msec(t_rdg, t_mesh) << " ms, surface: " << get_msec(t_mesh, t_write) << " ms, cube writing: " << get_msec(t_write, get_time()) << " ms" << endl;

	log << "Wrote " << signed_rho_cube.get_path() << endl;
	log << "Wrote " << rdg_cube.get_path() << endl;
	log << "Wrote " << output_base.string() + "_values.dat" << " with " << kept_points << " points" << endl;
}

void do_combine_mo(options &opt)
{
	using namespace std;
	WFN wavy1(e_origin::wfn);
	WFN wavy2(e_origin::wfn);
	WFN wavy3(e_origin::wfn);
	wavy1.read_wfn(opt.combine_mo[0], false, std::cout);
	wavy2.read_wfn(opt.combine_mo[1], false, std::cout);
	for (int i = 0; i < wavy1.get_ncen(); i++)
	{
		wavy3.push_back_atom(wavy1.get_atom(i));
	}
	for (int i = 0; i < wavy2.get_ncen(); i++)
	{
		wavy3.push_back_atom(wavy2.get_atom(i));
	}
	std::cout << "In total we have " << wavy3.get_ncen() << " atoms" << endl;

	readxyzMinMax_fromWFN(wavy1, opt.properties);
	properties_options prop2 = opt.properties;
	readxyzMinMax_fromWFN(wavy2, prop2);

	std::cout << "Read input\nCalculating for MOs ";
	for (int v1 = 0; v1 < opt.cmo1.size(); v1++)
	{
		std::cout << opt.cmo1[v1] << " ";
	}
	std::cout << "of fragment 1 and MOs ";
	for (int v1 = 0; v1 < opt.cmo2.size(); v1++)
	{
		std::cout << opt.cmo2[v1] << " ";
	}
	std::cout << "of fragment 2" << endl;
	std::array<double, 6> MinMax = { 100, 100, 100, -100, -100, -100 };
	std::array<int, 3> steps = { 0, 0, 0 };
	for (int i = 0; i < 3; i++)
	{
		if (opt.properties.MinMax[i] < MinMax[i])
			MinMax[i] = opt.properties.MinMax[i];
		if (opt.properties.MinMax[i + 3] > MinMax[i + 3])
			MinMax[i + 3] = opt.properties.MinMax[i + 3];
	}
	for (int i = 0; i < 3; i++)
	{
		if (prop2.MinMax[i] < MinMax[i])
			MinMax[i] = prop2.MinMax[i];
		if (prop2.MinMax[i + 3] > MinMax[i + 3])
			MinMax[i + 3] = prop2.MinMax[i + 3];
		steps[i] = (int)ceil(constants::bohr2ang(MinMax[i + 3] - MinMax[i]) / 0.1);
	}
	int counter = 0;
	cube total(steps, 0, true);
	cube MO1(steps, 0, true);
	MO1.give_parent_wfn(wavy3);
	MO1.set_na(wavy3.get_ncen());
	cube MO2(steps, 0, true);
	svec fns;
	for (int i = 0; i < 3; i++)
	{
		MO1.set_origin(i, MinMax[i]);
		MO1.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
		total.set_origin(i, MinMax[i]);
		total.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
		MO2.set_origin(i, MinMax[i]);
		MO2.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
	}
	for (int v1 = 0; v1 < opt.cmo1.size(); v1++)
	{
		MO1.set_zero();
		Calc_MO(MO1, opt.cmo1[v1] - 1, wavy1, 40, std::cout, false);
		for (int j = 0; j < opt.cmo2.size(); j++)
		{
			counter++;
			std::cout << "Running: " << counter << " of " << opt.cmo2.size() * opt.cmo1.size() << endl;
			string filename("");
			MO2.set_zero();
			Calc_MO(MO2, opt.cmo2[j] - 1, wavy2, 40, std::cout, false);
			std::cout << "writing files..." << flush;
			filename = wavy1.get_path().stem().string() + "_" + std::to_string(opt.cmo1[v1]) + "+" + wavy2.get_path().stem().string() + "_" + std::to_string(opt.cmo2[j]) + ".cube";
			fns.push_back(filename);
			total.set_zero();
			total = MO1;
			total += MO2;
			total.write_file(filename, false);
			filename = wavy1.get_path().stem().string() + "_" + std::to_string(opt.cmo1[v1]) + "-" + wavy2.get_path().stem().string() + "_" + std::to_string(opt.cmo2[j]) + ".cube";
			fns.push_back(filename);
			total.set_zero();
			total = MO1;
			total -= MO2;
			total.write_file(filename, false);
			std::cout << " ... done!" << endl;
		}
	}
	ofstream vmd("read_files.vmd");
	vmd << "mol addrep 0\nmol new {" + fns[0] + "} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 }\n";
	vmd << "animate style Loop\n";
	for (int i = 1; i < fns.size(); i++)
		vmd << "mol addfile {" + fns[i] + "} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n";
	vmd << "animate style Loop\ndisplay projection Orthographic\ndisplay depthcue off\n";
	vmd << "axes location Off\ndisplay rendermode GLSL\ncolor Display Background white\ncolor Element P purple\n";
	vmd << "color Element Ni green\ncolor Element C gray\nmol modstyle 0 0 CPK 1.000000 0.300000 12.000000 12.000000\n";
	vmd << "mol modcolor 0 0 Element\nmol color Element\nmol representation CPK 1.000000 0.300000 22.000000 22.000000\n";
	vmd << "mol selection all\nmol material Transparent\nmol addrep 0\nmol modstyle 1 0 Isosurface 0.020000 0 0 0 1 1\n";
	vmd << "mol modcolor 1 0 ColorID 0\nmol selection all\nmol material Transparent\nmol addrep 0\nmol modstyle 2 0 Isosurface -0.020000 0 0 0 1 1\nmol modcolor 2 0 ColorID 1\n";
	vmd << "mol selection all\nmol material Transparent\n";
	vmd.flush();
	vmd.close();
}

// end here
