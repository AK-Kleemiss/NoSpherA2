#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "constants.h"
#include "b2c.h"
#include "nos_math.h"
#include "citations.h"
#include "GridManager.h"
#include <limits>
#include <map>
#include <mutex>

namespace {

d3 get_axis_lengths(const cube *cub)
{
	return {
		std::sqrt(
			cub->get_vector(0, 0) * cub->get_vector(0, 0) +
			cub->get_vector(1, 0) * cub->get_vector(1, 0) +
			cub->get_vector(2, 0) * cub->get_vector(2, 0)),
		std::sqrt(
			cub->get_vector(0, 1) * cub->get_vector(0, 1) +
			cub->get_vector(1, 1) * cub->get_vector(1, 1) +
			cub->get_vector(2, 1) * cub->get_vector(2, 1)),
		std::sqrt(
			cub->get_vector(0, 2) * cub->get_vector(0, 2) +
			cub->get_vector(1, 2) * cub->get_vector(1, 2) +
			cub->get_vector(2, 2) * cub->get_vector(2, 2))
	};
}

// array_length(v) and array_length(a,b) are provided by convenience.h

double resolve_value_floor(const cube *cub, double value_floor)
{
	if (value_floor >= 0.0)
		return value_floor;
	return std::max(1e-8, cub->max_value() * 1e-6);
}

double resolve_gradient_cutoff(double gradient_epsilon)
{
	if (gradient_epsilon >= 0.0)
		return gradient_epsilon;
	return std::numeric_limits<double>::infinity();
}

d3 cube_gradient_at_point(const cube *cub, int x, int y, int z, const d3 &axis_lengths)
{
	return {
		(cub->get_value(x + 1, y, z) - cub->get_value(x - 1, y, z)) / (2.0 * axis_lengths[0]),
		(cub->get_value(x, y + 1, z) - cub->get_value(x, y - 1, z)) / (2.0 * axis_lengths[1]),
		(cub->get_value(x, y, z + 1) - cub->get_value(x, y, z - 1)) / (2.0 * axis_lengths[2])
	};
}

bool brackets_stationary_point(const cube *cub, int x, int y, int z)
{
	const double center = cub->get_value(x, y, z);
	const bool x_change = (center - cub->get_value(x - 1, y, z)) * (cub->get_value(x + 1, y, z) - center) <= 0.0;
	const bool y_change = (center - cub->get_value(x, y - 1, z)) * (cub->get_value(x, y + 1, z) - center) <= 0.0;
	const bool z_change = (center - cub->get_value(x, y, z - 1)) * (cub->get_value(x, y, z + 1) - center) <= 0.0;
	return x_change && y_change && z_change;
}

// Closed-form inverse of a row-major 3x3 matrix via the adjugate/determinant.
// Returns false (leaving inv untouched) if m is (numerically) singular.
// Used instead of a LAPACK call because this runs many times per Newton-Raphson
// iteration on a fixed-size 3x3 system, where solver call overhead would dominate.
bool invert_3x3(const double m[9], double inv[9])
{
	const double det =
		m[0] * (m[4] * m[8] - m[5] * m[7]) -
		m[1] * (m[3] * m[8] - m[5] * m[6]) +
		m[2] * (m[3] * m[7] - m[4] * m[6]);

	const double scale = std::max({ std::abs(m[0]), std::abs(m[1]), std::abs(m[2]),
									 std::abs(m[3]), std::abs(m[4]), std::abs(m[5]),
									 std::abs(m[6]), std::abs(m[7]), std::abs(m[8]) });
	const double threshold = std::max(1e-300, scale * scale * scale * 1e-12);
	if (!std::isfinite(det) || std::abs(det) < threshold)
		return false;

	const double inv_det = 1.0 / det;
	inv[0] = (m[4] * m[8] - m[5] * m[7]) * inv_det;
	inv[1] = (m[2] * m[7] - m[1] * m[8]) * inv_det;
	inv[2] = (m[1] * m[5] - m[2] * m[4]) * inv_det;
	inv[3] = (m[5] * m[6] - m[3] * m[8]) * inv_det;
	inv[4] = (m[0] * m[8] - m[2] * m[6]) * inv_det;
	inv[5] = (m[2] * m[3] - m[0] * m[5]) * inv_det;
	inv[6] = (m[3] * m[7] - m[4] * m[6]) * inv_det;
	inv[7] = (m[1] * m[6] - m[0] * m[7]) * inv_det;
	inv[8] = (m[0] * m[4] - m[1] * m[3]) * inv_det;
	return true;
}

d3 mat3_vec_mul(const double m[9], const d3 &v)
{
	return {
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
	};
}

bool is_seed_duplicate(const std::vector<critical_point_seed> &seeds, const critical_point_seed &candidate, double distance_tolerance)
{
	for (const critical_point_seed &seed : seeds) {
		if (array_length(seed.position, candidate.position) <= distance_tolerance)
			return true;
	}
	return false;
}

std::string classify_density_critical_point(int negative_count, int positive_count, int zero_count)
{
	if (zero_count > 0)
		return "degenerate";
	if (negative_count == 3)
		return "attractor";
	if (negative_count == 2 && positive_count == 1)
		return "bond";
	if (negative_count == 1 && positive_count == 2)
		return "ring";
	if (negative_count == 0 && positive_count == 3)
		return "cage";
	return "unknown";
}

} // namespace

critical_point evaluate_critical_point(
	const critical_point_seed &seed,
	const d3 &position,
	const WFN &wavy,
	int iterations,
	bool converged)
{
	critical_point result{};
	result.grid_index = seed.grid_index;
	result.seed_position = seed.position;
	result.position = position;
	result.seed_value = seed.value;
	result.iterations = iterations;
	result.converged = converged;
	result.ellipticity = std::numeric_limits<double>::quiet_NaN();
	result.virial_field = std::numeric_limits<double>::quiet_NaN();
	result.kinetic_lagrangian = std::numeric_limits<double>::quiet_NaN();
	result.kinetic_hamiltonian = std::numeric_limits<double>::quiet_NaN();
	result.lagrangian_density = std::numeric_limits<double>::quiet_NaN();

	d3 gradient{ 0.0, 0.0, 0.0 };
	wavy.computeGrad(position, gradient);
	result.gradient = gradient;
	result.gradient_norm = array_length(gradient);

	double rho = 0.0;
	double norm_grad = 0.0;
	double elf = 0.0;
	double eli = 0.0;
	double laplacian = 0.0;
	double hessian_data[9]{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	wavy.computeValues(position, rho, norm_grad, hessian_data, elf, eli, laplacian);
	result.density = rho;
	result.laplacian = laplacian;

	// Requested quantities from the rho field and its derivatives.
	// L is defined by user convention: L = K - G = (-1/4) * DelSqRho.
	const double grad2 = gradient[0] * gradient[0] + gradient[1] * gradient[1] + gradient[2] * gradient[2];
	if (rho > 1e-14 && std::isfinite(elf) && elf > 0.0 && elf < 1.0) {
		const double one_over_elf_minus_one = std::max(0.0, (1.0 / elf) - 1.0);
		const double d_term = std::pow(rho, constants::c_53) * std::sqrt(one_over_elf_minus_one) / constants::ctelf;
		const double tau = 2.0 * (d_term + 0.125 * grad2 / rho);
		const double g = 0.5 * tau;
		const double l = -0.25 * laplacian;
		const double k = g + l;
		//Local virial theorem, (1/4) DelSqRho = 2 G + V, with L = -(1/4) DelSqRho and K = G + L,
		//so V = -L - 2G = -(K + G). It used to read k - g, which is L again: every printed V
		//carried L's value, right magnitude at a bond critical point and the wrong sign
		const double v = -(k + g);
		result.kinetic_lagrangian = g;
		result.kinetic_hamiltonian = k;
		result.lagrangian_density = l;
		result.virial_field = v;
	}

	vec hessian_A(hessian_data, hessian_data + 9);
	vec hessian_W(3);
	if (try_make_Eigenvalues(hessian_A, hessian_W)) {
		const vec &eigenvalues = hessian_W;
		result.hessian_eigenvalues = { eigenvalues[0], eigenvalues[1], eigenvalues[2] };
		result.hessian_eigenvectors = {
			d3{ hessian_A[0 * 3 + 0], hessian_A[1 * 3 + 0], hessian_A[2 * 3 + 0] },
			d3{ hessian_A[0 * 3 + 1], hessian_A[1 * 3 + 1], hessian_A[2 * 3 + 1] },
			d3{ hessian_A[0 * 3 + 2], hessian_A[1 * 3 + 2], hessian_A[2 * 3 + 2] }
		};

		const double max_abs = std::max({ std::abs(eigenvalues[0]), std::abs(eigenvalues[1]), std::abs(eigenvalues[2]) });
		const double eigen_tolerance = std::max(1e-10, max_abs * 1e-8);
		result.negative_eigenvalues = 0;
		result.positive_eigenvalues = 0;
		result.zero_eigenvalues = 0;
		for (int i = 0; i < 3; i++) {
			if (eigenvalues[i] < -eigen_tolerance)
				result.negative_eigenvalues++;
			else if (eigenvalues[i] > eigen_tolerance)
				result.positive_eigenvalues++;
			else
				result.zero_eigenvalues++;
		}
		result.type = classify_density_critical_point(result.negative_eigenvalues, result.positive_eigenvalues, result.zero_eigenvalues);

		if (result.type == "bond" && std::abs(eigenvalues[1]) > eigen_tolerance)
			result.ellipticity = (eigenvalues[0] / eigenvalues[1]) - 1.0;
	}
	else {
		result.hessian_eigenvalues = { 0.0, 0.0, 0.0 };
		result.hessian_eigenvectors = { d3{ 0.0, 0.0, 0.0 }, d3{ 0.0, 0.0, 0.0 }, d3{ 0.0, 0.0, 0.0 } };
		result.negative_eigenvalues = 0;
		result.positive_eigenvalues = 0;
		result.zero_eigenvalues = 3;
		result.type = "unknown";
	}

	return result;
}

namespace {

bool try_merge_critical_point(std::vector<critical_point> &points, const critical_point &candidate, double distance_tolerance)
{
	for (critical_point &point : points) {
		if (array_length(point.position, candidate.position) > distance_tolerance)
			continue;

		const bool candidate_is_better = (candidate.converged && !point.converged) ||
			(candidate.converged == point.converged && candidate.gradient_norm < point.gradient_norm);
		if (candidate_is_better)
			point = candidate;
		return true;
	}
	return false;
}

} // namespace

bool b2c(const cube *cub, const std::vector<atom> &atoms, bool debug, bool bcp)
{
	using namespace std;
	iMatrix3 CP(cub->get_size(0), cub->get_size(1), cub->get_size(2));
	ivec2 Liste(3);
	double GradMax, xlength, ylength, zlength;
	xlength = std::sqrt(
		cub->get_vector(0, 0) * cub->get_vector(0, 0) +
		cub->get_vector(1, 0) * cub->get_vector(1, 0) +
		cub->get_vector(2, 0) * cub->get_vector(2, 0));

	ylength = std::sqrt(
		cub->get_vector(0, 1) * cub->get_vector(0, 1) +
		cub->get_vector(1, 1) * cub->get_vector(1, 1) +
		cub->get_vector(2, 1) * cub->get_vector(2, 1));

	zlength = std::sqrt(
		cub->get_vector(0, 2) * cub->get_vector(0, 2) +
		cub->get_vector(1, 2) * cub->get_vector(1, 2) +
		cub->get_vector(2, 2) * cub->get_vector(2, 2));
	if (debug)std::cout << "calculated lengths!" << endl;
	string s, s1, s2, s3;
	double xmin = cub->get_origin(0);//, xstep;
	double ymin = cub->get_origin(1);//, ystep;
	double zmin = cub->get_origin(2);//, zstep;

	if (debug)std::cout << "Resized CP and initialized it; GOING TO MAKE LISTE NOW" << endl;
	int iCP = 0, ListeMax;
	cubepoint Max = { 0,0,0,0.0 };
	vec Maxima;
	vec3 distances;
	svec Labels;
	ivec nrs;
	vector<vector<bvec>> border;
	ivec2 neighbours;
	vector<vector<cubepoint> > BCPs;
	distances.resize(3);
	for (int i = 0; i < 3; i++) {
		distances[i].resize(3);
		for (int j = 0; j < 3; j++) {
			distances[i][j].resize(3);
			for (int k = 0; k < 3; k++) distances[i][j][k] = sqrt(
				pow((1 - i) * xlength, 2)
				+ pow((1 - j) * ylength, 2)
				+ pow((1 - k) * zlength, 2));
		}
	}
	if (debug)std::cout << "calculated distances!" << endl;
	for (int x = 0; x < cub->get_size(0); x++)
		for (int y = 0; y < cub->get_size(1); y++)
			for (int z = 0; z < cub->get_size(2); z++) {
				if (CP(x, y, z) == 0) {
					ListeMax = 0;
					Max.x = x; Max.y = y; Max.z = z;
					int xs, ys, zs;
					do {
						ListeMax++;
						xs = Max.x, ys = Max.y, zs = Max.z;
						if (xs < cub->get_size(0) && ys < cub->get_size(1) && zs < cub->get_size(2))
							Max.value = (cub->get_value(xs, ys, zs));
						else
							Max.value = 0.0;
						if (Liste[0].size() <= ListeMax)
							Liste[0].push_back(Max.x);  //only append List when needed check if neccesarry
						else
							Liste[0][ListeMax - 1] = Max.x;
						if (Liste[1].size() <= ListeMax)
							Liste[1].push_back(Max.y);
						else
							Liste[1][ListeMax - 1] = Max.y;
						if (Liste[2].size() <= ListeMax)
							Liste[2].push_back(Max.z);
						else
							Liste[2][ListeMax - 1] = Max.z;
						GradMax = 0;
						for (int ix = xs - 1; ix < xs + 2; ix++)
							for (int iy = ys - 1; iy < ys + 2; iy++)
								for (int iz = zs - 1; iz < zs + 2; iz++) {
									if (ix == xs && iy == ys && iz == zs) continue;
									if (ix < 0 || iy < 0 || iz < 0) continue;
									if (ix >= cub->get_size(0) || iy >= cub->get_size(1) || iz >= cub->get_size(2)) continue;
									/*double dist=sqrt(
											 ((ix-xs)*xlength)*((ix-xs)*xlength)
											+((iy-ys)*ylength)*((iy-ys)*ylength)
											+((iz-zs)*zlength)*((iz-zs)*zlength));*/
											//    Tests if this voxel has maximum value compared to neighboring voxels
									if (((cub->get_value(ix, iy, iz) - (cub->get_value(xs, ys, zs))) / distances[1 + ix - xs][1 + iy - ys][1 + iz - zs]) > GradMax) {
										Max.x = ix; Max.y = iy; Max.z = iz;
										Max.value = cub->get_value(ix, iy, iz);
										GradMax = ((cub->get_value(ix, iy, iz)) - (cub->get_value(xs, ys, zs))) / distances[1 + ix - xs][1 + iy - ys][1 + iz - zs];
									}
								}
					} while (!(GradMax == 0 || CP(Max.x, Max.y, Max.z) > 0));
					if (CP(Max.x, Max.y, Max.z) > 0)
						for (int i = 0; i < ListeMax; i++)
							CP(Liste[0][i], Liste[1][i], Liste[2][i]) = CP(Max.x, Max.y, Max.z);
					else {
						const d3 pos = { Max.x * xlength + cub->get_origin(0),Max.y * ylength + cub->get_origin(1),Max.z * zlength + cub->get_origin(2) };
						if (debug) {
							std::cout << "DBUG: Position of CP: ";
							for (int i = 0; i < 3; i++)std::cout << pos[i] << " ";
							std::cout << "\n";
						}
						double min_dist, temp;
						unsigned int atom = 0;
						Maxima.push_back(Max.value);
						iCP++;
						min_dist = 10000;
						for (int i = 0; i < ListeMax; i++)
							CP(Liste[0][i], Liste[1][i], Liste[2][i]) = iCP;
						for (int i = 0; i < atoms.size(); i++) {
							temp = array_length(pos, atoms[i].get_pos());
							//if(debug)std::cout << "DBUG: Distance to atom " << i  << " with label " << atoms[i].get_label() << " is " << temp;
							if (min_dist > temp) {
								min_dist = temp;
								atom = i;
								//if(debug)std::cout << " *";
							}
							//if(debug)std::cout << endl;
						}
						Labels.push_back(atoms[atom].get_label());
						nrs.push_back(atom);
					}
					for (int i = 0; i < 3; i++)
						Liste[i].resize(0); //Revert Liste and free memory
				}
			}
	if (bcp) {
		neighbours.resize(iCP);
		BCPs.resize(iCP);
		border.resize(cub->get_size(0));
		for (int i = 0; i < cub->get_size(0); i++) {
			border[i].resize(cub->get_size(1));
			for (int j = 0; j < cub->get_size(1); j++) {
				border[i][j].resize(cub->get_size(2));
				for (int k = 0; k < cub->get_size(2); k++)
					border[i][j][k] = false;
			}
		}
		for (int x = 0; x < cub->get_size(0); x++)
			for (int y = 0; y < cub->get_size(1); y++)
				for (int z = 0; z < cub->get_size(2); z++)
					for (int ix = x - 1; ix < x + 2; ix++)
						for (int iy = y - 1; iy < y + 2; iy++)
							for (int iz = z - 1; iz < z + 2; iz++) {
								if (ix == x && iy == y && iz == z) continue;
								if (ix <= 0 || iy <= 0 || iz <= 0) continue;
								if (ix >= cub->get_size(0) || iy >= cub->get_size(1) || iz >= cub->get_size(2)) continue;
								//    Tests if this voxel has neighboring voxels, that are from a different basin and assigns list of neighbors for each basin
								if (CP(ix, iy, iz) != CP(x, y, z)) {
									border[x][y][z] = true;
									//CP holds 1-based basin ids, neighbours/BCPs are 0-based
									ivec& nb = neighbours[CP(x, y, z) - 1];
									const int other = CP(ix, iy, iz) - 1;
									if (find(nb.begin(), nb.end(), other) == nb.end())
										nb.push_back(other);
								}
							}
		//sanity check, all basins must be neighboring each other pairwise!
		for (int b = 0; b < iCP; b++)
			for (const int n : neighbours[b])
				if (find(neighbours[n].begin(), neighbours[n].end(), b) == neighbours[n].end()) {
					std::cout << "ERROR: Basins should be neighbours in a pairwise way! Basin " << b << " has neighbour " << n << ", but it does not appear to be the case the other way around!\n";
					return false;
				}
		//possibly better if we assigned max and min values for coordinates of basins to reduce amount of calculations, lets see...
		for (int b = 0; b < iCP; b++) {
			Max = { 0,0,0,0.0 };
			for (int n = 0; n < neighbours[b].size(); n++) {
				for (int x = 0; x < cub->get_size(0); x++)
					for (int y = 0; y < cub->get_size(1); y++)
						for (int z = 0; z < cub->get_size(2); z++) {
							//check that we are on A) A border B) in the basin we want C) the neighboring basin we look at is at this border D) the value is bigger than the previously found maximum
							if (!border[x][y][z])
								continue;
							if (CP(x, y, z) != b + 1)
								continue;
							bool found = false;
							for (int ix = x - 1; ix < x + 2; ix++)
								for (int iy = y - 1; iy < y + 2; iy++)
									for (int iz = z - 1; iz < z + 2; iz++) {
										if (ix == x && iy == y && iz == z) continue;
										if (ix <= 0 || iy <= 0 || iz <= 0) continue;
										if (ix >= cub->get_size(0) || iy >= cub->get_size(1) || iz >= cub->get_size(2)) continue;
										if (CP(ix, iy, iz) != neighbours[b][n] + 1) continue;
										else found = true;
									}
							if (found && cub->get_value(x, y, z) > Max.value)
								Max = { x,y,z,cub->get_value(x,y,z) };
						}
				BCPs[b].push_back(Max);
			}
		}
		for (int b = 0; b < iCP; b++)
			for (int n = 0; n < neighbours[b].size(); n++) {
				int back_reference = 0;
				for (int i = 0; i < neighbours[neighbours[b][n]].size(); i++)
					if (neighbours[neighbours[b][n]][i] == b)
						back_reference = i;
				if (BCPs[b][n].x - BCPs[neighbours[b][n]][back_reference].x > 2 || BCPs[b][n].y - BCPs[neighbours[b][n]][back_reference].y > 2 || BCPs[b][n].z - BCPs[neighbours[b][n]][back_reference].z > 2) {
					std::cout << "The BCPs on both sides of the ZFS are too far apart...\n";
				}
				else {
					std::cout << "BCP: " << Labels[b] << "-" << Labels[neighbours[b][n]] << " ED: " << (BCPs[b][n].value + BCPs[neighbours[b][n]][back_reference].value) / 2 << " Position: "
						<< (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2 * cub->get_vector(0, 0)
						+ (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2 * cub->get_vector(1, 0)
						+ (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2 * cub->get_vector(2, 0)
						+ xmin << " "
						<< (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2 * cub->get_vector(0, 1)
						+ (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2 * cub->get_vector(1, 1)
						+ (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2 * cub->get_vector(2, 1)
						+ ymin
						<< (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2 * cub->get_vector(0, 2)
						+ (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2 * cub->get_vector(1, 2)
						+ (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2 * cub->get_vector(2, 2)
						+ zmin
						<< "\n";
				}
			}
		if (debug)std::cout << "done with BCPs" << endl;
	}
	if (debug)
		std::cout << "done with liste, writing basins now!" << endl;
	ivec basins;
	vec EDS(iCP);
	vec VOL(iCP);
	double dv = cub->get_dv();
	if (debug)std::cout << "dv: " << dv << " iCP: " << iCP << endl;
	for (int x = 0; x < cub->get_size(0); x++)
		for (int y = 0; y < cub->get_size(1); y++)
			for (int z = 0; z < cub->get_size(2); z++)
				for (int a = 0; a < iCP; a++) if (a + 1 == CP(x, y, z)) {
					EDS[a] += cub->get_value(x, y, z) * dv;            //    Electrons in this voxel
					VOL[a] += dv;                                    //    Size of the voxel
				}
	if (debug) for (int a = 0; a < iCP; a++)std::cout << a << " eds: " << EDS[a] << " vol: " << VOL[a] << endl;
	unsigned int in = 9999999;
	std::cout << "I found " << iCP << " Basins." << endl;
	for (int a = 0; a < iCP; a++)std::cout << "Integrated value in Basin " << toString<int>(a + 1) << ": " << scientific << setw(14) << setprecision(7) << EDS[a]
		<< " Volume: " << VOL[a] << " Maximum: " << Maxima[a] << " possibly atom: " << Labels[a] << "_" << nrs[a] << "\n";
	if (debug) {
		std::cout << "DEBUG: Labels_size: " << Labels.size() << endl;
		for (int i = 0; i < Labels.size(); i++)std::cout << Labels[i] << " ";
		std::cout << endl;
	}

	std::cout << "Which of these do you want to include into one cube file?" << endl
		<< "Put 0 for end of input. If none are selected (first number is zero), all basins will be written into separate files" << endl;
	while (in != 0) {
		cin >> in;
		if (in != 0) basins.push_back(in);
	}
	string temp;
	string replace("E");
	const string base = (cub->get_path().parent_path() / cub->get_path().stem()).generic_string();
	ofstream logfile((base + ".b2c_log").c_str(), ios::out);
	logfile << "Number of Basins: " << iCP << endl;
	double Integral = 0.0;
	temp = base;
	if (basins.size() > 0) {
		temp += "_" + toString<int>((int)basins.size()) + "_basins";
		temp += ".cube";
		ofstream outfile(temp.c_str(), ios::out);
		outfile << "First comment line?" << endl << "second comment line" << endl << "   " << atoms.size();
		outfile << fixed;
		outfile << setprecision(6);
		outfile << setw(12);
		outfile << xmin << " " << ymin << " " << zmin << endl;
		for (int i = 0; i < 3; i++) {
			outfile << setprecision(0);
			outfile << setw(6);
			outfile << cub->get_size(i) << " ";
			outfile << setw(12);
			outfile << setprecision(6);
			for (int j = 0; j < 3; j++) outfile << cub->get_vector(i, j) << " ";
			outfile << "\n";
		}
		for (int i = 0; i < atoms.size(); i++) {
			outfile << "  " << atoms[i].get_charge();
			outfile << setw(12) << setprecision(6);
			outfile << atoms[i].get_charge() << " ";
			outfile << atoms[i].get_coordinate(0) << " " << atoms[i].get_coordinate(1) << " " << atoms[i].get_coordinate(2) << " ";
			outfile << "\n";
		}
		outfile.flush();
		stringstream stream;
		for (int x = 0; x < cub->get_size(0); x++) {
			for (int y = 0; y < cub->get_size(1); y++) {
				unsigned int r = 0;
				for (int z = 0; z < cub->get_size(2); z++) {
					bool include = false;
					for (int a = 0; a < basins.size(); a++)
						if (basins[a] == CP(x, y, z))
							include = true;
					stream << uppercase << scientific << setw(14) << setprecision(7) << cub->get_value(x, y, z) * include;
					outfile << stream.str();
					stream.str("");
					r++;
					if (r % 6 == 0) outfile << "\n";
				}
				if (r % 6 != 0) outfile << "\n";
			}
		}
		for (int a = 0; a < basins.size(); a++) {
			logfile << "Integrated value in Basin " << toString<int>(basins[a]) << ": " << scientific << setw(14) << setprecision(7) << EDS[basins[a] - 1]
				<< "Volume: " << VOL[basins[a] - 1] << " Atom: " << Labels[basins[a] - 1] << "_" << nrs[basins[a] - 1] << "\n";
			Integral += EDS[basins[a] - 1];
		}
		outfile.close();
	}
	else {
		for (int f = 0; f < iCP; f++) {
			temp = base + '_' + Labels[f] + '_' + toString<unsigned int>(nrs[f]) + '_' + toString<int>(f) + ".cube";
			ofstream outfile(temp.c_str(), ios::out);
			outfile << s1 << "\n" << s2 << "\n" << setw(5) << atoms.size();
			outfile << fixed << setprecision(6) << setw(12) << xmin << " " << ymin << " " << zmin << "\n";
			for (int i = 0; i < 3; i++)
				outfile << setprecision(0) << setw(5) << cub->get_size(i) << setw(12) << setprecision(6) << cub->get_vector(i, 0) << setw(12) << setprecision(6) << cub->get_vector(i, 1) << setw(12) << setprecision(6) << cub->get_vector(i, 2) << "\n";
			for (int i = 0; i < atoms.size(); i++) {
				outfile << setw(5) << atoms[i].get_charge() << setw(5) << atoms[i].get_charge() << ".000000";
				outfile << fixed << setw(12) << setprecision(6) << atoms[i].get_coordinate(0);
				outfile << fixed << setw(12) << setprecision(6) << atoms[i].get_coordinate(1);
				outfile << fixed << setw(12) << setprecision(6) << atoms[i].get_coordinate(2);
				outfile << "\n";
			}
			for (int x = 0; x < cub->get_size(0); x++) {
				for (int y = 0; y < cub->get_size(1); y++) {
					unsigned int r = 0;
					for (int z = 0; z < cub->get_size(2); z++) {
						outfile << uppercase << scientific << setw(13) << setprecision(5) << cub->get_value(x, y, z) * (CP(x, y, z) == (f + 1));
						r++;
						if (r % 6 == 0) outfile << "\n";
					}
					if (r % 6 != 0) outfile << "\n";
				}
			}
			logfile << "Integrated value in Basin " << toString<int>(f + 1) << ": " << scientific << setw(14) << setprecision(7) << EDS[f] << "Volume: " << VOL[f]
				<< " Atom: " << Labels[f] << "_" << nrs[f] << "\n";
			outfile.close();
			Integral += EDS[f];
		}
	}
	logfile << "Integral over all basins: " << Integral << endl;
	if (bcp)
		for (int b = 0; b < iCP; b++)
			for (int n = 0; n < neighbours[b].size(); n++) {
				int back_reference = 0;
				for (int i = 0; i < neighbours[neighbours[b][n]].size(); i++)
					if (neighbours[neighbours[b][n]][i] == b)
						back_reference = i;
				auto pos = cub->get_pos(
					(BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2,
					(BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2,
					(BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2
				);
				logfile << "BCP: " << Labels[b] << "-" << Labels[neighbours[b][n]] << " ED: " << (BCPs[b][n].value + BCPs[neighbours[b][n]][back_reference].value) / 2 << " Position: "
					<< pos[0] << " "
					<< pos[1] << " "
					<< pos[2] << endl;
			}
	logfile.flush();
	logfile.close();
	return true;
};

// Structure to store pre-computed gradient direction for each grid point
struct GradientDirection {
	int next_x, next_y, next_z;  // Next point to follow (-1 if local maximum or invalid)
	double gradient;              // Gradient magnitude
};

std::vector<critical_point_seed> find_cube_critical_point_seeds(const cube *cub, bool debug, double value_floor, double gradient_epsilon)
{
	const i3 sizes = cub->get_sizes();
	if (sizes[0] < 3 || sizes[1] < 3 || sizes[2] < 3)
		return {};

	const d3 axis_lengths = get_axis_lengths(cub);
	const double density_floor = resolve_value_floor(cub, value_floor);
	const double gradient_cutoff = resolve_gradient_cutoff(gradient_epsilon);
	const double distance_tolerance = 0.75 * std::min({ axis_lengths[0], axis_lengths[1], axis_lengths[2] });

	std::vector<critical_point_seed> seeds;
	for (int x = 1; x < sizes[0] - 1; x++) {
		for (int y = 1; y < sizes[1] - 1; y++) {
			for (int z = 1; z < sizes[2] - 1; z++) {
				const double value = cub->get_value(x, y, z);
				if (value <= density_floor)
					continue;
				if (!brackets_stationary_point(cub, x, y, z))
					continue;

				const d3 gradient = cube_gradient_at_point(cub, x, y, z, axis_lengths);
				const double gradient_norm = array_length(gradient);
				if (!std::isfinite(gradient_norm) || gradient_norm > gradient_cutoff)
					continue;

				// Accept this point as a seed if no axis-aligned neighbour has a strictly
				// smaller finite-difference gradient norm.  We deliberately check only the
				// 6 axis-aligned neighbours (not all 26) and use a small relative tolerance
				// so that equidistant nuclei (where adjacent grid points share the same
				// gradient norm) are not silently dropped.
				bool local_minimum = true;
				const double rel_tol = gradient_norm * 1e-4 + 1e-12;
				const int dx[6] = { -1, 1, 0, 0, 0, 0 };
				const int dy[6] = { 0, 0, -1, 1, 0, 0 };
				const int dz[6] = { 0, 0, 0, 0, -1, 1 };
				for (int n = 0; n < 6 && local_minimum; n++) {
					const int ix = x + dx[n], iy = y + dy[n], iz = z + dz[n];
					if (ix <= 0 || iy <= 0 || iz <= 0 || ix >= sizes[0] - 1 || iy >= sizes[1] - 1 || iz >= sizes[2] - 1)
						continue;
					if (cub->get_value(ix, iy, iz) <= density_floor)
						continue;
					const d3 neighbor_gradient = cube_gradient_at_point(cub, ix, iy, iz, axis_lengths);
					const double neighbor_norm = array_length(neighbor_gradient);
					if (neighbor_norm < gradient_norm - rel_tol)
						local_minimum = false;
				}
				if (!local_minimum)
					continue;

				critical_point_seed seed{
					{ x, y, z },
					cub->get_pos(x, y, z),
					value,
					gradient_norm
				};
				if (!is_seed_duplicate(seeds, seed, distance_tolerance))
					seeds.push_back(seed);
			}
		}
	}

	std::sort(seeds.begin(), seeds.end(), [](const critical_point_seed &left, const critical_point_seed &right) {
		if (left.gradient_norm != right.gradient_norm)
			return left.gradient_norm < right.gradient_norm;
		return left.value > right.value;
	});

	if (debug)
		std::cout << "Found " << seeds.size() << " cube-based critical-point seeds" << std::endl;
	return seeds;
}

std::vector<critical_point> refine_cube_critical_points(
	const cube *cub,
	const WFN &wavy,
	const std::vector<critical_point_seed> &seeds,
	bool debug,
	double value_floor,
	double gradient_tolerance,
	double step_tolerance,
	int max_iterations)
{
	const d3 axis_lengths = get_axis_lengths(cub);
	const double density_floor = resolve_value_floor(cub, value_floor);
	const double trust_radius = 1.5 * std::max({ axis_lengths[0], axis_lengths[1], axis_lengths[2] });
	const double nuclear_trust_radius = 0.25;
	const double nuclear_max_radius = 0.35;
	// Merge converged CPs within 0.1 bohr of each other (safely smaller than any bond length).
	const double merge_distance = 0.1;

	std::vector<critical_point> points;
	for (const critical_point_seed &seed : seeds) {
		// Check if this seed is a known nuclear seed (fallback to proximity check for legacy seeds).
		bool is_nuclear_seed = seed.is_nuclear_seed;
		int seed_nucleus = seed.nucleus_index;
		if (!is_nuclear_seed) {
			for (int a = 0; a < wavy.get_ncen(); a++) {
				if (array_length(seed.position, wavy.get_atom_pos(a)) < 0.2) {
					is_nuclear_seed = true;
					seed_nucleus = a;
					break;
				}
			}
		}
		d3 nuclear_center = seed_nucleus >= 0 ? wavy.get_atom_pos(seed_nucleus) : seed.position;
		// For nuclear seeds, use a much lower density floor (or none at all) since H nuclei have very low density.
		const double local_density_floor = is_nuclear_seed ? 1e-10 : density_floor;
		const double local_trust_radius = is_nuclear_seed ? nuclear_trust_radius : trust_radius;
		const int local_max_iterations = is_nuclear_seed ? std::max(max_iterations, 64) : max_iterations;

		d3 position = seed.position;
		d3 gradient{ 0.0, 0.0, 0.0 };
		wavy.computeGrad(position, gradient);
		double gradient_norm = array_length(gradient);

		bool converged = gradient_norm <= gradient_tolerance;
		int iterations = 0;
		for (; iterations < local_max_iterations && !converged; iterations++) {
			double rho = 0.0;
			double norm_grad = 0.0;
			double elf = 0.0;
			double eli = 0.0;
			double laplacian = 0.0;
			double hessian_data[9]{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			wavy.computeValues(position, rho, norm_grad, hessian_data, elf, eli, laplacian);
			if (rho <= local_density_floor)
				break;

			double hessian_inv[9];
			if (!invert_3x3(hessian_data, hessian_inv))
				break;

			d3 step = mat3_vec_mul(hessian_inv, gradient);
			step = { -step[0], -step[1], -step[2] };
			if (!std::isfinite(step[0]) || !std::isfinite(step[1]) || !std::isfinite(step[2]))
				break;
			const double step_norm = array_length(step);
			if (step_norm > local_trust_radius)
				step = { step[0] * local_trust_radius / step_norm,
						 step[1] * local_trust_radius / step_norm,
						 step[2] * local_trust_radius / step_norm };

			bool accepted = false;
			d3 accepted_position = position;
			d3 accepted_gradient = gradient;
			double accepted_norm = gradient_norm;
			double accepted_step_norm = 0.0;

			for (int attempt = 0; attempt < 8; attempt++) {
				const double damping = std::pow(0.5, attempt);
				const d3 trial_step{ step[0] * damping, step[1] * damping, step[2] * damping };
				d3 trial_position{
					position[0] + trial_step[0],
					position[1] + trial_step[1],
					position[2] + trial_step[2]
				};
				if (is_nuclear_seed && array_length(trial_position, nuclear_center) > nuclear_max_radius)
					continue;
				d3 trial_gradient{ 0.0, 0.0, 0.0 };
				wavy.computeGrad(trial_position, trial_gradient);
				const double trial_norm = array_length(trial_gradient);
				if (!std::isfinite(trial_norm))
					continue;
				const double trial_step_norm = array_length(trial_step);
				if (trial_norm <= gradient_norm || trial_step_norm <= step_tolerance) {
					accepted = true;
					accepted_position = trial_position;
					accepted_gradient = trial_gradient;
					accepted_norm = trial_norm;
					accepted_step_norm = trial_step_norm;
					break;
				}
			}

			if (!accepted)
				break;

			position = accepted_position;
			gradient = accepted_gradient;
			gradient_norm = accepted_norm;
			converged = gradient_norm <= gradient_tolerance;

			if (accepted_step_norm <= step_tolerance && gradient_norm <= 10.0 * gradient_tolerance)
				converged = true;
		}

		critical_point point = evaluate_critical_point(seed, position, wavy, iterations, converged);

		// Check if this CP is at a known atomic position (nuclear attractor).
		// Nuclear attractors should be kept even if density is very low (e.g., H atoms),
		// and even if not fully converged, as long as the Hessian indicates it's an attractor.
		bool is_nuclear_attractor = false;
		if (point.type == "attractor") {
			for (int a = 0; a < wavy.get_ncen(); a++) {
				if (array_length(point.position, wavy.get_atom_pos(a)) < 0.5) {
					is_nuclear_attractor = true;
					break;
				}
			}
		}

		// For non-nuclear CPs, require that they exceed the normal density floor
		if (!is_nuclear_attractor && point.density <= density_floor)
			continue;
		// Discard non-converged bond CPs; these are typically spurious basin-edge artifacts.
		if (point.type == "bond" && !point.converged)
			continue;
		if (!try_merge_critical_point(points, point, merge_distance))
			points.push_back(point);
	}

	std::sort(points.begin(), points.end(), [](const critical_point &left, const critical_point &right) {
		if (left.type != right.type)
			return left.type < right.type;
		return left.density > right.density;
	});

	if (debug)
		std::cout << "Refined " << points.size() << " unique critical points from cube seeds" << std::endl;
	return points;
}
std::vector<critical_point> analyze_cube_critical_points(
	const cube *cub,
	const WFN &wavy,
	bool debug,
	double value_floor,
	double gradient_epsilon,
	double gradient_tolerance,
	double step_tolerance,
	int max_iterations)
{
	// Grid seeds (non-nuclear CPs: bond, ring, cage)
	std::vector<critical_point_seed> seeds = find_cube_critical_point_seeds(cub, debug, value_floor, gradient_epsilon);

	// Invert the grid matrix once to convert WFN atom positions → grid indices.
	// cube::get_pos: pos = origin + V * idx  →  idx = V^{-1} * (pos − origin)
	double V[9];
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
			V[r * 3 + c] = cub->get_vector(r, c);
	double Vinv[9];
	err_checkf(invert_3x3(V, Vinv), "Cube grid matrix is singular!", std::cout);

	const i3 sizes = cub->get_sizes();
	const double density_floor = resolve_value_floor(cub, value_floor);
	const d3 axis_lengths = get_axis_lengths(cub);
	const double min_step = std::min({ axis_lengths[0], axis_lengths[1], axis_lengths[2] });
	const double seed_merge_dist = 0.75 * min_step;

	// Helper: project a real-space position onto the nearest interior grid index.
	auto pos_to_idx = [&](const d3 &pos) -> i3 {
		const d3 rhs{ pos[0] - cub->get_origin(0),
					  pos[1] - cub->get_origin(1),
					  pos[2] - cub->get_origin(2) };
		const d3 fidx = mat3_vec_mul(Vinv, rhs);
		return {
			std::max(1, std::min(sizes[0] - 2, (int)std::round(fidx[0]))),
			std::max(1, std::min(sizes[1] - 2, (int)std::round(fidx[1]))),
			std::max(1, std::min(sizes[2] - 2, (int)std::round(fidx[2])))
		};
	};

	auto add_seed_at = [&](const d3 &pos, bool skip_density_check = false, bool force_add = false, bool is_nuclear_seed = false, int nucleus_index = -1) {
		const i3 idx = pos_to_idx(pos);
		const double cv = cub->get_value(idx[0], idx[1], idx[2]);
		if (!skip_density_check && cv <= density_floor)
			return;
		// Compute approximate grid gradient to supply a gradient_norm hint.
		const d3 gv = cube_gradient_at_point(cub, idx[0], idx[1], idx[2], axis_lengths);
		critical_point_seed s{ idx, pos, cv, array_length(gv), is_nuclear_seed, nucleus_index };
		if (force_add || !is_seed_duplicate(seeds, s, seed_merge_dist))
			seeds.push_back(s);
	};

	// Add every WFN atom as an explicit nuclear-attractor seed.
	// Nuclear attractors should always exist, even if density is low (e.g., H atoms).
	for (int a = 0; a < wavy.get_ncen(); a++) {
		add_seed_at(wavy.get_atom_pos(a), true, true, true, a);  // force-add every nuclear seed
	}

	// Add bond-critical-point seeds along every atom pair that is plausibly bonded
	// (interatomic distance <= 1.3 * sum of CSD covalent radii).
	// We inject seeds at multiple positions along the internuclear vector to ensure
	// coverage even when the BCP is not exactly at the midpoint (e.g. polar bonds).
	for (int a = 0; a < wavy.get_ncen(); a++) {
		const d3 pa = wavy.get_atom_pos(a);
		const int za = wavy.get_atom_charge(a);
		const double ra = (za > 0 && za < 114) ? constants::covalent_radii[za] : 1.5;
		for (int b = a + 1; b < wavy.get_ncen(); b++) {
			const d3 pb = wavy.get_atom_pos(b);
			const int zb = wavy.get_atom_charge(b);
			const double rb = (zb > 0 && zb < 114) ? constants::covalent_radii[zb] : 1.5;
			const double dist = array_length(pa, pb);
			const double bond_threshold = constants::ang2bohr(1.3 * (ra + rb));
			if (dist > bond_threshold)
				continue;
			// Seeds at 10% to 90% along the bond (9 seeds per bond)
			for (int i = 1; i < 10; i++) {
				const double t = i * 0.1;
				const d3 bpos{
					pa[0] + t * (pb[0] - pa[0]),
					pa[1] + t * (pb[1] - pa[1]),
					pa[2] + t * (pb[2] - pa[2])
				};
				add_seed_at(bpos, false);  // skip_density_check=false for bond seeds
			}
		}
	}

	if (debug)
		std::cout << "Total seeds (grid + atom + bond): " << seeds.size() << std::endl;

	// Refine all the seed critical points
	auto points = refine_cube_critical_points(cub, wavy, seeds, debug, value_floor, gradient_tolerance, step_tolerance, max_iterations);

	return points;
}

//Assign gridded basins by near-grid steepest ascent.
std::pair<cubei, std::vector<d4>> topological_cube_analysis(const cube *cub, const std::vector<atom> &atoms, bool debug, bool bcp, double value_floor, double grad_epsilon, double assignment_radius, double merge_persistence, const std::vector<d3> *seeds, const WFN *field_wfn, const std::function<double(const d3&)> *core_density, const std::function<void(const d3&, d3&)> *core_gradient, const density_field *field)
{
	auto field_dens = [&](const d3 &p) { return (field ? field->rho(p) : field_wfn->compute_dens(p)) + (core_density ? (*core_density)(p) : 0.0); };
	auto field_grad = [&](const d3 &p, d3 &g) {
		if (field) field->grad(p, g); else field_wfn->computeGrad(p, g);
		if (core_gradient) { d3 c; (*core_gradient)(p, c); for (int k = 0; k < 3; k++) g[k] += c[k]; }
	};
	const int nx = cub->get_size(0), ny = cub->get_size(1), nz = cub->get_size(2);
	cubei basin_cube({ nx, ny, nz }, 0, true);
	d3 h;
	for (int i = 0; i < 3; i++)
		h[i] = std::sqrt(cub->get_vector(0, i) * cub->get_vector(0, i) + cub->get_vector(1, i) * cub->get_vector(1, i) + cub->get_vector(2, i) * cub->get_vector(2, i));
	const size_t n = static_cast<size_t>(nx) * ny * nz;
	auto lin = [&](int x, int y, int z) { return (static_cast<size_t>(x) * ny + y) * nz + z; };
	vec v(n);
	std::vector<unsigned char> valid(n, 0);
	const double r2 = assignment_radius > 0.0 ? std::pow(constants::ang2bohr(assignment_radius), 2) : -1.0;
#pragma omp parallel for collapse(3) schedule(static)
	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++) {
				const size_t i = lin(x, y, z);
				v[i] = cub->get_value(x, y, z);
				if (v[i] <= value_floor) continue;
				bool in = r2 <= 0.0;
				if (!in) {
					const d3 pos = cub->get_pos(x, y, z);
					for (const atom &a : atoms) {
						const d3 ap = a.get_pos();
						if (std::pow(pos[0] - ap[0], 2) + std::pow(pos[1] - ap[1], 2) + std::pow(pos[2] - ap[2], 2) <= r2) { in = true; break; }
					}
				}
				valid[i] = in ? 1 : 0;
			}
	auto ok = [&](int x, int y, int z) { return x >= 0 && y >= 0 && z >= 0 && x < nx && y < ny && z < nz && valid[lin(x, y, z)]; };
	//Measure sensitivity to the order of voxel updates.
	long long steep_calls = 0, steep_tied = 0, tied_paths = 0;
	bool path_tied = false;
	//Measure sensitivity to the ascent step.
	long long marginal_1e9 = 0, marginal_1e5 = 0, step_decisions = 0;
	double min_margin = 1.0;
	//And the one that decides where the attractors land: grad_epsilon is ABSOLUTE, so on a field whose
	//values are around 11 a top that varies in the eleventh digit counts as flat and the walk stops at
	//whichever voxel of it the path entered. flat_stops counts those stops and max_flat_best says how
	//much uphill was still there when the walk gave up; a count of zero would refute this too.
	long long flat_stops = 0;
	double max_flat_best = 0.0;
	//Highest 26-neighbour; false when none is higher
	auto steepest = [&](int x, int y, int z, int &bx, int &by, int &bz) {
		double best = 0.0;
		int n_best = 0;
		bx = x; by = y; bz = z;
		const double c = v[lin(x, y, z)];
		for (int ix = x - 1; ix <= x + 1; ix++)
			for (int iy = y - 1; iy <= y + 1; iy++)
				for (int iz = z - 1; iz <= z + 1; iz++) {
					if ((ix == x && iy == y && iz == z) || !ok(ix, iy, iz)) continue;
					const double d = std::sqrt(std::pow((ix - x) * h[0], 2) + std::pow((iy - y) * h[1], 2) + std::pow((iz - z) * h[2], 2));
					const double g = (v[lin(ix, iy, iz)] - c) / d;
					if (g > best) { best = g; bx = ix; by = iy; bz = iz; n_best = 1; }
					else if (n_best && g == best) n_best++;   //just as steep, reached later: the winner was picked by loop order alone
				}
		steep_calls++;
		if (n_best > 1) { steep_tied++; path_tied = true; }
		//Measure sensitivity to the grid spacing.
		if (best > 0.0 && best <= grad_epsilon) {
			flat_stops++;
			if (best > max_flat_best) max_flat_best = best;
		}
		return best > grad_epsilon;
	};
	ivec basin(n, 0);
	std::vector<d4> Maxima;
	std::vector<unsigned char> seeded(n, 0);
	ivec stamp(n, 0);
	int path_id = 0;
	if (seeds)
		for (const d3 &p : *seeds) {
			int c[3];
			bool inside = true;
			for (int d = 0; d < 3 && inside; d++) {
				c[d] = static_cast<int>(std::lround((p[d] - cub->get_origin(d)) / cub->get_vector(d, d)));
				inside = c[d] >= 0 && c[d] < (d == 0 ? nx : d == 1 ? ny : nz);
			}
			if (!inside || !valid[lin(c[0], c[1], c[2])] || basin[lin(c[0], c[1], c[2])]) continue;
			Maxima.push_back(d4{ p[0], p[1], p[2], v[lin(c[0], c[1], c[2])] });
			basin[lin(c[0], c[1], c[2])] = static_cast<int>(Maxima.size());
			seeded[lin(c[0], c[1], c[2])] = 1;
		}
	const int n_seeded = static_cast<int>(Maxima.size());
	ivec path;
	//The first pass stops at any assigned point; the refinement, given the interior mask, only
	//at a point no differently assigned neighbour touches
	auto ascend = [&](int x, int y, int z, const std::vector<unsigned char> *interior, const bool assign) {
		path.clear();
		path_id++;
		path_tied = false;
		d3 dr{ 0.0, 0.0, 0.0 };
		int cx = x, cy = y, cz = z;
		for (size_t guard = 0; guard < n; guard++) {
			const size_t ci = lin(cx, cy, cz);
			if (basin[ci] != 0 && (!interior || (*interior)[ci])) break;
			path.push_back(static_cast<int>(ci));
			stamp[ci] = path_id;
			//The analytic gradient when there is one, else central differences where both
			//sides exist and one-sided at the rim; in index units either way
			d3 s;
			if (field_wfn) {
				field_grad(cub->get_pos(cx, cy, cz), s);
				for (int d = 0; d < 3; d++) s[d] /= h[d];
			}
			else for (int d = 0; d < 3; d++) {
				const int px = cx + (d == 0), py = cy + (d == 1), pz = cz + (d == 2);
				const int mx = cx - (d == 0), my = cy - (d == 1), mz = cz - (d == 2);
				const bool hp = ok(px, py, pz), hm = ok(mx, my, mz);
				if (hp && hm) s[d] = (v[lin(px, py, pz)] - v[lin(mx, my, mz)]) / (2.0 * h[d]);
				else if (hp) s[d] = (v[lin(px, py, pz)] - v[ci]) / h[d];
				else if (hm) s[d] = (v[ci] - v[lin(mx, my, mz)]) / h[d];
				else s[d] = 0.0;
				s[d] /= h[d];
			}
			const double m = std::max({ std::abs(s[0]), std::abs(s[1]), std::abs(s[2]) });
			int nxp = cx, nyp = cy, nzp = cz;
			bool moved = false;
			if (m > 0.0) {
				int step[3];
				for (int d = 0; d < 3; d++) {
					const double f = s[d] / m;
					//how far is this component from the nearest lround() boundary at a half-integer?
					const double margin = std::abs(std::abs(f) - std::floor(std::abs(f)) - 0.5);
					step_decisions++;
					if (margin < min_margin) min_margin = margin;
					if (margin < 1e-9) marginal_1e9++;
					if (margin < 1e-5) marginal_1e5++;
					step[d] = static_cast<int>(std::lround(f));
					dr[d] += f - step[d];
					if (dr[d] > 0.5) { step[d]++; dr[d] -= 1.0; }
					else if (dr[d] < -0.5) { step[d]--; dr[d] += 1.0; }
				}
				nxp = cx + step[0]; nyp = cy + step[1]; nzp = cz + step[2];
				//The analytic gradient is trusted over the grid values, which a cusp sampled
				//at a tenth of an angstrom does not order; a step back onto this path means
				//the maximum lies between voxels and the current one stands for it
				moved = (step[0] || step[1] || step[2]) && ok(nxp, nyp, nzp)
					&& (field_wfn ? stamp[lin(nxp, nyp, nzp)] != path_id : v[lin(nxp, nyp, nzp)] > v[ci]);
			}
			if (!moved) {
				dr = { 0.0, 0.0, 0.0 };
				if (!steepest(cx, cy, cz, nxp, nyp, nzp)) {
					int id = basin[ci];
					if (id == 0) {
						const d3 pos = cub->get_pos(cx, cy, cz);
						//The voxel beside a seed can top the seed's own, and a nucleus whose core
						//an ECP removed wears a shell of maxima where its valence density peaks:
						//within a bohr of a seed a maximum is the seed's
						for (int m = 0; m < n_seeded && id == 0; m++)
							if (std::pow(pos[0] - Maxima[m][0], 2) + std::pow(pos[1] - Maxima[m][1], 2) + std::pow(pos[2] - Maxima[m][2], 2) < std::max(1.0, 2.25 * std::pow(std::max({ h[0], h[1], h[2] }), 2)))
								id = m + 1;
					}
					if (id == 0) {
						const d3 pos = cub->get_pos(cx, cy, cz);
						Maxima.push_back(d4{ pos[0], pos[1], pos[2], v[ci] });
						id = static_cast<int>(Maxima.size());
					}
					if (assign) for (const int q : path) basin[q] = id;
					if (path_tied) tied_paths++;
					return id;
				}
			}
			cx = nxp; cy = nyp; cz = nzp;
		}
		const int id = basin[lin(cx, cy, cz)];
		if (assign) for (const int q : path) basin[q] = id;
		if (path_tied) tied_paths++;
		return id;
	};
	if (field_wfn) {
		//A cusp basin two voxels across cannot be climbed voxel by voxel, so with the analytic
		//gradient every voxel sends a continuous trajectory instead: it ends within a voxel and
		//a half of a seed, in a voxel some earlier trajectory settled, or where the gradient
		//dies, and every voxel it crossed takes the answer. Threads share the answers as they
		//come; a stale read only makes a trajectory run a little further.
		std::cout << "Assigning basins along the density gradient..." << std::endl;
		const double hmax = std::max({ h[0], h[1], h[2] });
		const double catch2 = 2.25 * hmax * hmax;
		const double hmin = std::min({ h[0], h[1], h[2] });
		std::vector<d3> seed_pos;
		for (int m = 0; m < n_seeded; m++) seed_pos.push_back(d3{ Maxima[m][0], Maxima[m][1], Maxima[m][2] });
		//Densest voxels first, in chunks whose answers are applied together, so a trajectory
		//may stop in a voxel an earlier chunk settled and never in one its own chunk is still
		//deciding: the result does not depend on the thread count
		ivec order;
		for (size_t i = 0; i < n; i++) if (basin[i] == 0 && valid[i]) order.push_back(static_cast<int>(i));
		std::sort(order.begin(), order.end(), [&](int a, int b) { return v[a] > v[b] || (v[a] == v[b] && a < b); });
		const size_t chunk = std::max<size_t>(1, order.size() / 64 + 1);
		ivec result(n, 0);
		ivec unresolved;
		auto voxel_of = [&](const d3 &p, int *c) {
			for (int d = 0; d < 3; d++) {
				c[d] = static_cast<int>(std::lround((p[d] - cub->get_origin(d)) / cub->get_vector(d, d)));
				if (c[d] < 0 || c[d] >= (d == 0 ? nx : d == 1 ? ny : nz)) return false;
			}
			return true;
		};
		for (size_t c0 = 0; c0 < order.size(); c0 += chunk) {
			const size_t c1 = std::min(order.size(), c0 + chunk);
			unresolved.clear();
#pragma omp parallel
		{
			ivec crossed;
#pragma omp for schedule(dynamic, 64)
			for (long long oi = static_cast<long long>(c0); oi < static_cast<long long>(c1); oi++) {
				const size_t i = order[oi];
				crossed.clear();
				crossed.push_back(static_cast<int>(i));
				const int x = static_cast<int>(i / (static_cast<size_t>(ny) * nz)), y = static_cast<int>((i / nz) % ny), z = static_cast<int>(i % nz);
				d3 r = cub->get_pos(x, y, z), g;
				int id = 0;
				double last_rho = -1.0;
				for (int step = 0; step < 4000 && id == 0; step++) {
					for (int m = 0; m < n_seeded; m++)
						if (std::pow(r[0] - seed_pos[m][0], 2) + std::pow(r[1] - seed_pos[m][1], 2) + std::pow(r[2] - seed_pos[m][2], 2) < catch2) { id = m + 1; break; }
					if (id) break;
					double near2 = 1e300;
					for (const d3 &sp : seed_pos) near2 = std::min(near2, std::pow(r[0] - sp[0], 2) + std::pow(r[1] - sp[1], 2) + std::pow(r[2] - sp[2], 2));
					const double sl = (near2 < 1.0 ? 0.3 : 0.6) * hmin;
					//A nucleus an ECP hollowed out has a sphere of maxima around it instead
					//of a cusp; a trajectory that stops climbing has reached such a top
					const double rho_here = field_dens(r);
					if (rho_here <= last_rho) break;
					last_rho = rho_here;
					field_grad(r, g);
					double gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
					if (gn < 1e-14) break;
					d3 mid;
					for (int k = 0; k < 3; k++) mid[k] = r[k] + 0.5 * sl * g[k] / gn;
					field_grad(mid, g);
					gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
					if (gn < 1e-14) break;
					for (int k = 0; k < 3; k++) r[k] += sl * g[k] / gn;
					int c[3];
					if (!voxel_of(r, c)) break;
					const size_t ci = lin(c[0], c[1], c[2]);
					if (!valid[ci]) break;
					if (ci != static_cast<size_t>(crossed.back())) {
						const int b = basin[ci];
						if (b != 0 && std::find(crossed.begin(), crossed.end(), static_cast<int>(ci)) == crossed.end()) { id = b; break; }
						if (std::find(crossed.begin(), crossed.end(), static_cast<int>(ci)) == crossed.end()) crossed.push_back(static_cast<int>(ci));
					}
				}
				if (id == 0)
					//A top within a bohr of a seed is the seed's - the ECP shell again
					for (int m = 0; m < n_seeded && id == 0; m++)
						if (std::pow(r[0] - seed_pos[m][0], 2) + std::pow(r[1] - seed_pos[m][1], 2) + std::pow(r[2] - seed_pos[m][2], 2) < 1.0) id = m + 1;
				if (id == 0) {
					//Nothing caught it: a maximum between voxels, a non-nuclear one, or a
					//trajectory that left the region; the highest voxel crossed stands for it,
					//resolved in order once the chunk is done
					size_t top = crossed[0];
					for (const int q : crossed) if (v[q] > v[top]) top = q;
#pragma omp critical
					unresolved.push_back(static_cast<int>(i));
					result[i] = -static_cast<int>(top) - 1;
				}
				else result[i] = id;
			}
		}
			std::sort(unresolved.begin(), unresolved.end());
			for (const int i : unresolved) {
				const size_t top = static_cast<size_t>(-result[i] - 1);
				int id = basin[top];
				const int tx = static_cast<int>(top / (static_cast<size_t>(ny) * nz)), ty = static_cast<int>((top / nz) % ny), tz = static_cast<int>(top % nz);
				const d3 pos = cub->get_pos(tx, ty, tz);
				for (size_t m = 0; m < Maxima.size() && id == 0; m++)
					if (std::pow(pos[0] - Maxima[m][0], 2) + std::pow(pos[1] - Maxima[m][1], 2) + std::pow(pos[2] - Maxima[m][2], 2) < catch2) id = static_cast<int>(m) + 1;
				if (id == 0) {
					Maxima.push_back(d4{ pos[0], pos[1], pos[2], v[top] });
					id = static_cast<int>(Maxima.size());
				}
				result[i] = id;
			}
			for (size_t oi = c0; oi < c1; oi++) basin[order[oi]] = result[order[oi]];
		}
	}
	else {
		std::cout << "Assigning basins by near-grid ascent..." << std::endl;
		//An assigned voxel terminates ascent only after its basin is settled.
		const bool reverse_scan = std::getenv("NOS_BASIN_SCAN_REVERSE") != nullptr;
		if (reverse_scan) std::cout << "NOS_BASIN_SCAN_REVERSE is set: scanning seed voxels in the opposite order (diagnostic)" << std::endl;
		for (int xr = 0; xr < nx; xr++)
			for (int yr = 0; yr < ny; yr++)
				for (int zr = 0; zr < nz; zr++) {
					const int x = reverse_scan ? nx - 1 - xr : xr;
					const int y = reverse_scan ? ny - 1 - yr : yr;
					const int z = reverse_scan ? nz - 1 - zr : zr;
					const size_t i = lin(x, y, z);
					if (basin[i] == 0 && valid[i]) ascend(x, y, z, nullptr, true);
				}
	}
	//Refinement: a point with a differently assigned 6-neighbour is sent up again and may only
	//settle on an interior point or a maximum; the new labels apply once all are known
	std::vector<unsigned char> interior(n, 1);
	const int dx6[6] = { 1, -1, 0, 0, 0, 0 }, dy6[6] = { 0, 0, 1, -1, 0, 0 }, dz6[6] = { 0, 0, 0, 0, 1, -1 };
	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++) {
				const size_t i = lin(x, y, z);
				if (!valid[i]) { interior[i] = 0; continue; }
				for (int k = 0; k < 6; k++) {
					const int ix = x + dx6[k], iy = y + dy6[k], iz = z + dz6[k];
					if (ok(ix, iy, iz) && basin[lin(ix, iy, iz)] != basin[i]) { interior[i] = 0; break; }
				}
			}
	if (!field_wfn) {
		ivec refined(basin);
		for (size_t i = 0; i < n; i++) {
			if (!valid[i] || interior[i]) continue;
			const int x = static_cast<int>(i / (static_cast<size_t>(ny) * nz)), y = static_cast<int>((i / nz) % ny), z = static_cast<int>(i % nz);
			refined[i] = ascend(x, y, z, &interior, false);
		}
		basin.swap(refined);
		if (basin_timing_enabled()) {
			std::cout << "Near-grid ascent: " << steep_calls << " steepest-neighbour steps, " << steep_tied
				<< " of them tied (" << std::fixed << std::setprecision(3)
				<< (steep_calls ? 100.0 * static_cast<double>(steep_tied) / static_cast<double>(steep_calls) : 0.0)
				<< " %), " << tied_paths << " voxel paths went through at least one tie" << std::endl;
			std::cout << "Near-grid ascent: " << step_decisions << " rounded step decisions, " << marginal_1e5
				<< " within 1e-5 of flipping and " << marginal_1e9 << " within 1e-9, closest margin "
				<< std::scientific << std::setprecision(3) << min_margin << std::defaultfloat << std::endl;
			std::cout << "Near-grid ascent: " << flat_stops << " walks stopped on a top that grad_epsilon ("
				<< std::scientific << std::setprecision(1) << grad_epsilon << ") calls flat, largest uphill "
				<< std::setprecision(3) << max_flat_best << " still available there" << std::defaultfloat << std::endl;
		}
	}
	int nb = static_cast<int>(Maxima.size());
	std::cout << "I found " << nb << " Basins." << std::endl;
	//Exclude maxima on the density-crop rim.
	std::vector<char> rim(nb + 1, 0);
	int n_rim = 0;
	for (int b = 1; b <= nb; b++) {
		if (b <= n_seeded) continue;
		int c[3];
		bool inside = true;
		for (int d = 0; d < 3 && inside; d++) {
			c[d] = static_cast<int>(std::lround((Maxima[b - 1][d] - cub->get_origin(d)) / cub->get_vector(d, d)));
			inside = c[d] >= 0 && c[d] < (d == 0 ? nx : d == 1 ? ny : nz);
		}
		if (!inside) continue;
		for (int k = 0; k < 6 && !rim[b]; k++)
			if (!ok(c[0] + dx6[k], c[1] + dy6[k], c[2] + dz6[k])) rim[b] = 1;
		n_rim += rim[b];
	}
	//Persistence merge: the saddle between two basins is the highest of the lower values over
	//their shared faces; a maximum less than merge_persistence of its height above its highest
	//saddle is grid noise and joins the basin behind that saddle
	if ((merge_persistence > 0.0 || n_rim > 0) && nb > 1) {
		std::map<std::pair<int, int>, double> pass;
		for (int x = 0; x < nx; x++)
			for (int y = 0; y < ny; y++)
				for (int z = 0; z < nz; z++) {
					const size_t i = lin(x, y, z);
					const int a = basin[i];
					if (a == 0) continue;
					for (int k = 0; k < 6; k += 2) {
						const int ix = x + dx6[k], iy = y + dy6[k], iz = z + dz6[k];
						if (!ok(ix, iy, iz)) continue;
						const int b = basin[lin(ix, iy, iz)];
						if (b == 0 || b == a) continue;
						const double s = std::min(v[i], v[lin(ix, iy, iz)]);
						double &e = pass[{ std::min(a, b), std::max(a, b) }];
						if (s > e) e = s;
					}
				}
		ivec target(nb + 1);
		for (int b = 0; b <= nb; b++) target[b] = b;
		auto root = [&](int b) { while (target[b] != b) b = target[b]; return b; };
		for (;;) {
			//The least persistent basin first, so a chain of noise collapses in order
			int worst = -1, into = -1;
			double worst_rel = std::max(merge_persistence, 0.0);
			for (int b = 1; b <= nb; b++) {
				if (root(b) != b || b <= n_seeded) continue;
				double saddle = -1.0;
				int nb_into = -1;
				for (const auto &e : pass) {
					const int p = root(e.first.first), q = root(e.first.second);
					if (p == q || (p != b && q != b)) continue;
					const int other = p == b ? q : p;
					if (e.second > saddle && Maxima[other - 1][3] >= Maxima[b - 1][3]) { saddle = e.second; nb_into = other; }
				}
				if (nb_into < 0) continue;
				const double rel = (Maxima[b - 1][3] - saddle) / Maxima[b - 1][3];
				if (rel < worst_rel) { worst_rel = rel; worst = b; into = nb_into; }
			}
			if (worst < 0) break;
			if (debug) std::cout << "Merging basin " << worst << " into " << into << " (persistence " << worst_rel << ")\n";
			target[worst] = into;
		}
		//Report unresolved rim voxels outside the basins.
		int rim_dropped = 0;
		for (int b = 1; b <= nb; b++)
			if (rim[b] && root(b) == b) { target[b] = 0; rim_dropped++; }
		if (rim_dropped)
			std::cout << "Left " << rim_dropped << " maxima on the rim of the analysed region unresolved; their density is reported outside every basin." << std::endl;
		ivec renumber(nb + 1, 0);
		std::vector<d4> kept;
		for (int b = 1; b <= nb; b++)
			if (root(b) == b) { kept.push_back(Maxima[b - 1]); renumber[b] = static_cast<int>(kept.size()); }
		for (size_t i = 0; i < n; i++) if (basin[i]) basin[i] = renumber[root(basin[i])];
		if (kept.size() != Maxima.size()) std::cout << "Merged " << Maxima.size() - kept.size() << " noise maxima, " << kept.size() << " basins remain." << std::endl;
		Maxima.swap(kept);
	}
	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++)
				basin_cube.set_value(x, y, z, basin[lin(x, y, z)]);
	return { basin_cube, Maxima };
}

double core_shell_radius(const int Z)
{
	if (Z <= 2) return 0.0;
	if (Z <= 10) return 0.25;
	//Use element-dependent core radii for basin grouping.
	if (Z <= 36) return 1.0;
	if (Z <= 54) return 1.4;
	return 1.8;
}

//Both merges below decide only WHICH maxima belong together, and then do the same three things with
//the answer: keep the highest of each group, renumber the cube, and report where each old basin went.
//keeper[b] is the group's representative - its lowest member - or -1 for a maximum that stands alone.
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
			const double r = core_shell_radius(atoms[a].get_charge());
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

//Merge grid-fragmented ELI-D shells by persistence.
int unify_shell_basins(cubei &basin_cube, std::vector<d4> &maxima, ivec *basin_map, double max_dist, double rel_tol)
{
	const int nb = static_cast<int>(maxima.size());
	if (nb < 2 || max_dist <= 0.0) return 0;
	const double d2 = max_dist * max_dist;
	ivec root(nb);
	for (int b = 0; b < nb; b++) root[b] = b;
	auto find = [&root](int b) { while (root[b] != b) b = root[b] = root[root[b]]; return b; };
	for (int b = 1; b < nb; b++)
		for (int c = 0; c < b; c++) {
			if (std::pow(maxima[b][0] - maxima[c][0], 2) + std::pow(maxima[b][1] - maxima[c][1], 2) + std::pow(maxima[b][2] - maxima[c][2], 2) > d2) continue;
			const double hi = std::max(std::abs(maxima[b][3]), std::abs(maxima[c][3]));
			if (hi > 0.0 && std::abs(maxima[b][3] - maxima[c][3]) > rel_tol * hi) continue;
			const int rb = find(b), rc = find(c);
			//union by lower index, so a group's root is always its lowest member
			if (rb != rc) root[std::max(rb, rc)] = std::min(rb, rc);
		}
	ivec keeper(nb);
	for (int b = 0; b < nb; b++) { const int r = find(b); keeper[b] = (r == b ? -1 : r); }
	return collapse_maxima_groups(basin_cube, maxima, keeper, basin_map);
}

//Newton-Raphson onto the nearest critical point of the field, then the negative-definite test.
//Newton converges to whatever critical point is nearest, of any type, which is the point: a
//candidate sitting next to a saddle comes back rejected rather than dragged uphill to some
//maximum elsewhere. The step is damped until it lowers the gradient norm, so a bad quadratic
//model costs iterations and not a runaway.
bool converge_to_maximum(const scalar_field &field, d3 &p, const double step_limit, const int max_iterations, const double gradient_tolerance)
{
	//Central differences of the analytic gradient: the cancellation at 1e-3 bohr leaves about
	//1e-13 of noise on a curvature of order one, far below what the definiteness test asks
	constexpr double fd = 1e-3;
	d3 g;
	if (!std::isfinite(field(p, g))) return false;
	for (int it = 0; it < max_iterations; it++) {
		double H[9];
		for (int k = 0; k < 3; k++) {
			d3 a = p, b = p, ga, gb;
			a[k] += fd;
			b[k] -= fd;
			field(a, ga);
			field(b, gb);
			for (int j = 0; j < 3; j++) H[3 * j + k] = (ga[j] - gb[j]) / (2.0 * fd);
		}
		for (int i = 0; i < 3; i++)
			for (int j = i + 1; j < 3; j++) {
				const double m = 0.5 * (H[3 * i + j] + H[3 * j + i]);
				H[3 * i + j] = H[3 * j + i] = m;
			}
		const double g0 = array_length(g);
		if (!std::isfinite(g0)) return false;
		if (g0 <= gradient_tolerance) {
			vec A(H, H + 9), W(3);
			if (!try_make_Eigenvalues(A, W)) return false;
			const double max_abs = std::max({ std::abs(W[0]), std::abs(W[1]), std::abs(W[2]) });
			const double tol = std::max(1e-10, max_abs * 1e-8);
			return W[0] < -tol && W[1] < -tol && W[2] < -tol;
		}
		double inv[9];
		if (!invert_3x3(H, inv)) return false;
		d3 s = mat3_vec_mul(inv, g);
		for (int k = 0; k < 3; k++) s[k] = -s[k];
		const double n = array_length(s);
		if (!std::isfinite(n) || n == 0.0) return false;
		if (n > step_limit)
			for (int k = 0; k < 3; k++) s[k] *= step_limit / n;
		bool accepted = false;
		for (int att = 0; att < 10 && !accepted; att++) {
			const double d = std::pow(0.5, att);
			const d3 q{ p[0] + d * s[0], p[1] + d * s[1], p[2] + d * s[2] };
			d3 gq;
			if (!std::isfinite(field(q, gq))) continue;
			if (array_length(gq) < g0) {
				p = q;
				g = gq;
				accepted = true;
			}
		}
		if (!accepted) return false;
	}
	return false;
}

//Accept non-nuclear attractors only when the analytic field confirms a maximum.
std::vector<d4> streaming_density_attractors(const WFN &wavy, const std::vector<critical_point> &critical_points, const std::function<double(const d3&)> *core_density, const std::function<void(const d3&, d3&)> *core_gradient, const bool debug)
{
	auto rho = [&](const d3 &p) { return wavy.compute_dens(p) + (core_density ? (*core_density)(p) : 0.0); };
	scalar_field field = [&](const d3 &p, d3 &g) {
		wavy.computeGrad(p, g);
		if (core_gradient) {
			d3 c;
			(*core_gradient)(p, c);
			for (int k = 0; k < 3; k++) g[k] += c[k];
		}
		return rho(p);
	};
	std::vector<d4> maxima;
	for (int a = 0; a < wavy.get_ncen(); a++) {
		const d3 p = wavy.get_atom_pos(a);
		maxima.push_back(d4{ p[0], p[1], p[2], rho(p) });
	}
	const size_t nuclei = maxima.size();
	constexpr double nuclear_radius = 0.5;   //a critical point this close to a nucleus is that nucleus
	constexpr double duplicate_radius2 = 0.01;

	//Keep perturbations inside the candidate attractor basin.
	auto perturbation_for = [&](const d3 &at, const double cap) {
		double nearest2 = std::numeric_limits<double>::max();
		auto note = [&](const d3 &q) {
			const double d2 = std::pow(at[0] - q[0], 2) + std::pow(at[1] - q[1], 2) + std::pow(at[2] - q[2], 2);
			//1e-4 bohr^2 is 0.01 bohr, the floor below - anything nearer is the candidate itself or a
			//second seed that converged onto it, and neither bounds anything.
			if (d2 > 1e-4) nearest2 = std::min(nearest2, d2);
		};
		for (const critical_point &o : critical_points)
			if (o.converged) note(o.position);
		for (const d4 &m : maxima) note(d3{ m[0], m[1], m[2] });
		//Limit the displacement used to test a candidate maximum.
		return std::min(cap, std::max(0.01, std::sqrt(nearest2) / 3.0));
	};

	//A candidate is an attractor if the analytic field has a maximum there: the point itself must
	//converge, and three starts displaced along the axes must all come back to the same place.
	auto accept = [&](const d3 &candidate, d3 &converged, const double cap = 0.1) {
		for (size_t a = 0; a < nuclei; a++)
			if (array_length(candidate, d3{ maxima[a][0], maxima[a][1], maxima[a][2] }) < nuclear_radius) return false;
		const double delta = perturbation_for(candidate, cap);
		const double back = 0.5 * delta;
		converged = candidate;
		int fail_t = -1; double fail_dist = 0.0; const char *why = "";
		for (int t = 0; t < 4; t++) {
			//t == 0 is the point itself; the three after it start delta off along each axis and
			//have to come back to the same place
			d3 q = candidate;
			if (t > 0) q[t - 1] += delta;
			if (!converge_to_maximum(field, q)) { fail_t = t; why = "the Newton iteration did not reach a maximum"; }
			else if (t == 0) { converged = q; continue; }
			else {
				fail_dist = array_length(q, converged);
				if (fail_dist < back) continue;
				fail_t = t; why = "it came back somewhere else";
			}
			//Record whether a candidate failed to converge or reached another maximum.
			if (debug) std::cout << "Dropped a non-nuclear attractor candidate that is not a maximum of the analytic field at "
				<< candidate[0] << " " << candidate[1] << " " << candidate[2]
				<< " (start " << fail_t << ": " << why << ", " << fail_dist
				<< " bohr away, displaced by " << delta << ")" << std::endl;
			return false;
		}
		for (const d4 &m : maxima)
			if (std::pow(converged[0] - m[0], 2) + std::pow(converged[1] - m[1], 2) + std::pow(converged[2] - m[2], 2) < duplicate_radius2) return false;
		return true;
	};

	auto keep = [&](const d3 &converged, const char *from) {
		maxima.push_back(d4{ converged[0], converged[1], converged[2], rho(converged) });
		if (debug) std::cout << "Kept a non-nuclear attractor at " << converged[0] << " " << converged[1] << " "
			<< converged[2] << " with rho " << rho(converged) << " (from " << from << ")" << std::endl;
	};

	d3 converged;
	for (const critical_point &cp : critical_points) {
		if (!cp.converged || cp.type != "attractor") continue;
		if (accept(cp.position, converged)) keep(converged, "the seed cube");
	}

	//Sample bond lines for maxima absent from the seed cube.
	const double line_step = 0.02;   //bohr. CCH's bump is 3.1e-5 e/bohr^3 over 0.19 bohr with an axial
	//curvature of -0.0249, so a step this size falls 5e-6 across it - six times finer than the bump it
	//has to resolve. It is an absolute length on purpose: a fraction of the bond would make the
	//nomination depend on the bond, which is the mistake being fixed here. Cost is one density per
	//point, about 5400 for all of sucrose's bonds against millions in the integration itself.
	for (int a = 0; a < wavy.get_ncen(); a++) {
		const d3 pa = wavy.get_atom_pos(a);
		const int za = wavy.get_atom_charge(a);
		const double ra = (za > 0 && za < 114) ? constants::covalent_radii[za] : 1.5;
		for (int b = a + 1; b < wavy.get_ncen(); b++) {
			const d3 pb = wavy.get_atom_pos(b);
			const int zb = wavy.get_atom_charge(b);
			const double rb = (zb > 0 && zb < 114) ? constants::covalent_radii[zb] : 1.5;
			const double dist = array_length(pa, pb);
			//The same 1.3 x sum of CSD covalent radii the cube's own bond seeds use
			if (dist > constants::ang2bohr(1.3 * (ra + rb)) || dist < 2.0 * line_step) continue;
			const int n = (int)(dist / line_step);
			const d3 u{ (pb[0] - pa[0]) / dist, (pb[1] - pa[1]) / dist, (pb[2] - pa[2]) / dist };
			auto at = [&](const int i) { return d3{ pa[0] + i * line_step * u[0], pa[1] + i * line_step * u[1], pa[2] + i * line_step * u[2] }; };
			std::vector<double> prof((size_t)n + 1);
			for (int i = 0; i <= n; i++) prof[(size_t)i] = rho(at(i));
			for (int i = 1; i < n; i++) {
				if (prof[i] <= prof[i - 1] || prof[i] <= prof[i + 1]) continue;
				//How far this maximum's basin reaches along the line, from the profile and nothing else:
				//the nearest turning point on either side. The cube's critical points would otherwise be
				//what bounds the displacement, and a candidate the bond line found precisely because the
				//cube was too coarse must not then be judged against what the cube found.
				int lo = i, hi = i;
				while (lo > 0 && prof[lo - 1] < prof[lo]) lo--;
				while (hi < n && prof[hi + 1] < prof[hi]) hi++;
				const double reach = line_step * std::min(i - lo, hi - i);
				if (accept(at(i), converged, std::min(0.1, std::max(0.01, reach / 3.0)))) keep(converged, "a bond line");
			}
		}
	}
	return maxima;
}

std::vector<d4> analytic_eli_maxima(const WFN &wavy, const bool debug)
{
	auto sphere = [](const int n) {
		std::vector<d3> d(n);
		for (int i = 0; i < n; i++) {
			const double z = 1.0 - 2.0 * (i + 0.5) / n;
			const double s = std::sqrt(std::max(0.0, 1.0 - z * z));
			const double phi = 2.39996322972865332 * i;
			d[i] = d3{ s * std::cos(phi), s * std::sin(phi), z };
		}
		return d;
	};
	const std::vector<d3> inner = sphere(26), outer = sphere(110);
	//bohr; the core shells of a fourth-row atom sit inside 0.1, a lone pair or bond maximum 0.6-2.5 out
	const double radii[] = { 0.03, 0.06, 0.1, 0.15, 0.22, 0.3, 0.4, 0.55, 0.75, 1.0, 1.3, 1.7, 2.2, 2.8, 3.5 };
	const int ncen = wavy.get_ncen();
	std::vector<d3> seeds;
	for (int a = 0; a < ncen; a++) {
		const d3 pa = wavy.get_atom_pos(a);
		seeds.push_back(pa);
		for (const double r : radii)
			for (const d3 &u : (r < 0.5 ? inner : outer)) {
				const d3 p{ pa[0] + r * u[0], pa[1] + r * u[1], pa[2] + r * u[2] };
				//Every seed belongs to the atom nearest to it, so no region is seeded twice
				bool own = true;
				for (int b = 0; b < ncen && own; b++)
					if (b != a && array_length(p, wavy.get_atom_pos(b)) < r) own = false;
				if (own) seeds.push_back(p);
			}
	}
	//Monotone ascent along the normalised gradient: a step is taken only if ELI-D rises, grows
	//after a success and halves after a failure. There is deliberately no Hessian test at the end
	//- a hydrogen's valence maximum sits on the cusp of the nucleus, and the Newton test ate six
	//real H basins in UH6 (see converge_to_maximum)
	const int ns = static_cast<int>(seeds.size());
	std::vector<d4> ends(ns, d4{ 0.0, 0.0, 0.0, -1.0 });
#pragma omp parallel for schedule(dynamic, 16)
	for (int s = 0; s < ns; s++) {
		d3 p = seeds[s], g;
		if (wavy.compute_dens(p) < basin_density_cutoff) continue;
		double f;
		wavy.computeELIGrad(p, f, g);
		if (!std::isfinite(f)) continue;
		bool inside = true;
		double step = 0.05;
		for (int it = 0; it < 500 && step > 1e-5; it++) {
			const double gn = array_length(g);
			if (!std::isfinite(gn) || gn < 1e-12) break;
			const d3 t{ p[0] + step * g[0] / gn, p[1] + step * g[1] / gn, p[2] + step * g[2] / gn };
			double ft;
			d3 gt;
			wavy.computeELIGrad(t, ft, gt);
			if (std::isfinite(ft) && ft > f) {
				p = t; f = ft; g = gt;
				step = std::min(1.5 * step, 0.3);
				//ELI-D is undefined in the tail and the basins end at the same isosurface
				if (wavy.compute_dens(p) < basin_density_cutoff) { inside = false; break; }
			}
			else step *= 0.5;
		}
		if (inside) ends[s] = d4{ p[0], p[1], p[2], f };
	}
	//Highest first; an end within 0.1 bohr of a kept one is that maximum reached again
	std::vector<int> order;
	for (int s = 0; s < ns; s++) if (ends[s][3] >= 0.0) order.push_back(s);
	std::stable_sort(order.begin(), order.end(), [&](const int x, const int y) { return ends[x][3] > ends[y][3]; });
	std::vector<d4> maxima;
	for (const int s : order) {
		bool seen = false;
		for (const d4 &m : maxima)
			if (std::pow(m[0] - ends[s][0], 2) + std::pow(m[1] - ends[s][1], 2) + std::pow(m[2] - ends[s][2], 2) < 0.01) { seen = true; break; }
		if (!seen) maxima.push_back(ends[s]);
	}
	if (debug)
		std::cout << "ELI-D maxima from " << ns << " seeds by analytic ascent: " << order.size() << " climbs finished, " << maxima.size() << " distinct maxima" << std::endl;
	return maxima;
}

static bool g_beta_spheres = true;
void beta_spheres_set_enabled(const bool on) { g_beta_spheres = on; }
bool beta_spheres_enabled() { return g_beta_spheres; }
//The adaptive step uses field turning to grow on straight trajectories.
static constexpr double adp_cap_default = 8.0;    //at most this many times the validated floor step
static constexpr double adp_grow_default = 0.99999;   //midpoint cosine that earns a doubling
static constexpr double adp_keep_default = 0.999;   //below this the step is thrown away and retaken at the floor
static constexpr double adp_reach_default = 0.25;  //fraction of the distance to the nearest maximum
static double g_adp_cap = adp_cap_default;
static double g_adp_grow = adp_grow_default;
static double g_adp_keep = adp_keep_default;
static double g_adp_reach = adp_reach_default;
//Measure adaptive-step settings against basin populations.
static double env_double(const char *name, const double fallback, const double upper = 0.0)
{
	const char *v = std::getenv(name); // Flawfinder: ignore - parsed as one positive double
	if (v == nullptr || *v == '\0') return fallback;
	try {
		const double d = std::stod(v);
		if (d > 0.0 && std::isfinite(d) && (upper <= 0.0 || d <= upper)) return d;
		std::cout << "Ignoring " << name << "=" << v << ": not a positive finite number";
		if (upper > 0.0) std::cout << " of at most " << upper;
		std::cout << std::endl;
	}
	catch (const std::exception &) {
		std::cout << "Ignoring " << name << "=" << v << ": not a number" << std::endl;
	}
	return fallback;
}
static void adp_knobs_from_env()
{
	//From the defaults every time, so clearing the variables puts the validated numbers back
	g_adp_cap = env_double("NOS_ADP_CAP", adp_cap_default);
	g_adp_grow = env_double("NOS_ADP_GROW", adp_grow_default);
	g_adp_keep = env_double("NOS_ADP_KEEP", adp_keep_default);
	g_adp_reach = env_double("NOS_ADP_REACH", adp_reach_default);
}
//Keep beta spheres within the smallest sampled safe radius.
static constexpr double beta_margin_default = 0.9;
double basin_beta_margin() { return env_double("NOS_BETA_MARGIN", beta_margin_default, 1.0); }
void basin_adaptive_step_knobs(double &cap, double &grow, double &keep, double &reach)
{
	cap = g_adp_cap; grow = g_adp_grow; keep = g_adp_keep; reach = g_adp_reach;
}
//Keep adaptive step growth optional.
static bool g_adaptive_step = false;
//Count rejected adaptive-step proposals separately.
static std::atomic<long long> g_adp_steps{ 0 }, g_adp_tries{ 0 }, g_adp_turn{ 0 }, g_adp_fall{ 0 }, g_adp_shrink{ 0 };
//Trajectories that ran their step budget out. Counted separately from the stalls because they are
//the opposite failure - a walk that never gave up rather than one that gave up too early - and
//because a streaming walk that reaches the cap has no cube to fall back on.
static std::atomic<long long> g_adp_exhaust{ 0 };
//Steps a single trajectory may take before it is given up on. 2000 was chosen for a walk whose step
//had a hard floor; the shrink took that floor away, so whether it is still enough is a measurement,
//and NOS_BASIN_STEP_CAP is how that measurement gets taken without a rebuild per value.
static int g_step_cap = 2000;
//Recover from a shortened step by doubling toward the base step.
static double g_step_relax = 1.5;
//Whether a floor step that fails to rise may halve and try again. On by default - it is the fix, and
//3716 of CCH's 3716 stalls were this - and NOS_BASIN_SHRINK=0 turns it off, which is how the test that
//asserts the stall count is zero gets to show the count it is asserting against.
static bool g_step_shrink = true;
//Account separately for density that never reaches a basin.
static std::atomic<long long> g_stall_vacuum{ 0 }, g_stall_field{ 0 }, g_stall_far{ 0 };
static std::atomic<double> g_stall_far_rho{ 0.0 };
//Bucket stalled gradients to distinguish tail and interior failures.
static constexpr int g_stall_basins = 8;
static std::atomic<long long> g_stall_gn[4] = {};
//The >=1e-2 bucket again, but cleared when an integration starts instead of when it prints, so it
//still holds the last field's count after the call returns. That is what basin_stalls_on_a_slope
//hands to the test; the rest of these are read once by the log line and reset in the same breath.
static std::atomic<long long> g_stall_slope{ 0 };
static std::atomic<long long> g_stall_to[g_stall_basins + 1] = {};
//A negative distance marks a stall below the density cutoff.
static void basin_stall_seen(const double dist, const double rho, const double gn = -1.0)
{
	if (dist < 0.0) { g_stall_vacuum.fetch_add(1, std::memory_order_relaxed); return; }
	g_stall_field.fetch_add(1, std::memory_order_relaxed);
	if (gn >= 0.0) {
		const int b = gn < 1e-6 ? 0 : gn < 1e-4 ? 1 : gn < 1e-2 ? 2 : 3;
		g_stall_gn[b].fetch_add(1, std::memory_order_relaxed);
		if (b == 3) g_stall_slope.fetch_add(1, std::memory_order_relaxed);
	}
	if (dist <= 1.0) return;
	g_stall_far.fetch_add(1, std::memory_order_relaxed);
	double cur = g_stall_far_rho.load(std::memory_order_relaxed);
	while (rho > cur && !g_stall_far_rho.compare_exchange_weak(cur, rho, std::memory_order_relaxed)) {}
}
long long basin_stalls_on_a_slope()
{
	return g_stall_slope.load();
}
static void basin_stall_gave_to(const int basin)
{
	g_stall_to[basin >= 0 && basin < g_stall_basins ? basin : g_stall_basins].fetch_add(1, std::memory_order_relaxed);
}
static void basin_stall_counters(long long &vacuum, long long &in_field, long long &beyond_a_bohr, double &worst_rho)
{
	vacuum = g_stall_vacuum.load();
	in_field = g_stall_field.load();
	beyond_a_bohr = g_stall_far.load();
	worst_rho = g_stall_far_rho.load();
}
static void basin_stall_counters_reset()
{
	g_stall_vacuum = 0; g_stall_field = 0; g_stall_far = 0; g_stall_far_rho = 0.0;
	for (auto &c : g_stall_gn) c = 0;
	for (auto &c : g_stall_to) c = 0;
}
static inline void adp_count(std::atomic<long long> &c) { c.fetch_add(1, std::memory_order_relaxed); }
void basin_adaptive_step_counters(long long &steps, long long &proposed, long long &turned_back, long long &fell_back)
{
	steps = g_adp_steps.load();
	proposed = g_adp_tries.load();
	turned_back = g_adp_turn.load();
	fell_back = g_adp_fall.load();
}
void basin_adaptive_step_counters_reset()
{
	g_adp_steps = 0; g_adp_tries = 0; g_adp_turn = 0; g_adp_fall = 0; g_adp_shrink = 0; g_adp_exhaust = 0;
}
void basin_adaptive_step_set_enabled(const bool on) { g_adaptive_step = on; if (on) adp_knobs_from_env(); }
bool basin_adaptive_step_enabled() { return g_adaptive_step; }
static double g_basin_step_scale = 1.0;
void basin_step_scale_set(const double f) { g_basin_step_scale = f > 0.0 ? f : 1.0; }
double basin_step_scale() { return g_basin_step_scale; }
static bool g_basin_timing = false;
void basin_timing_set_enabled(const bool on) { g_basin_timing = on; }
bool basin_timing_enabled() { return g_basin_timing; }
void basin_stage_timer::lap(const std::string &what) {
	const auto now = std::chrono::steady_clock::now();
	const double s = std::chrono::duration<double>(now - t).count();
	t = now;
	if (g_basin_timing) std::cout << "  [timing] " << what << ": " << std::fixed << std::setprecision(2) << s << " s" << std::endl;
}

//Integrate basin populations on atom-centred quadrature grids.
vec integrate_basins_on_atomic_grids(const cube *cub, const cubei *basin_cube, const std::vector<d4> &maxima, const WFN &wavy, const int accuracy, const bool eli_field, vec &volumes, double &outside, const std::function<double(const d3&)> *core_density, const std::function<void(const d3&, d3&)> *core_gradient, const int grid_boost, const density_field *field, basin_overlaps *ovl, const ivec *maximum_basin)
{
	//The filled core steers the trajectories only. An ECP atom's grid is built for its
	//valence basis and cannot integrate a 1s at Z = 80, so the core electrons are added to
	//the nucleus's basin by count once the valence density is integrated; a Thakkar core
	//lies whole inside its atom's basin
	auto valence = [&](const d3 &p) { return field ? field->rho(p) : wavy.compute_dens(p); };
	//Streaming: no cube and no basin cube, the maxima are the whole topology and every point
	//finds its basin by walking the field
	const bool streaming = cub == nullptr || basin_cube == nullptr;
	//Cleared here and not at the print, so the count belongs to this call and survives it
	g_stall_slope = 0;
	if (const char *e = std::getenv("NOS_BASIN_STEP_CAP")) {
		const int v = std::atoi(e);
		//A cap below a few hundred steps would make every trajectory unfinished and the measurement
		//meaningless, so a typo reads as "leave it alone" rather than as a setting.
		if (v >= 100) g_step_cap = v;
	}
	if (const char *e = std::getenv("NOS_BASIN_STEP_RELAX")) {
		const double v = std::atof(e);
		//Below 1.0 the step would shrink on a *successful* step too, which is a different algorithm
		//and not one anybody asked for; above 4.0 it is a cycle with a longer period.
		if (v >= 1.0 && v <= 4.0) g_step_relax = v;
	}
	if (const char *e = std::getenv("NOS_BASIN_SHRINK")) g_step_shrink = std::string(e) != "0";
	//The basin a maximum belongs to, 1-based; without a map that is the maximum's own index
	auto basin_of = [&](const size_t m) { return maximum_basin ? (*maximum_basin)[m + 1] : static_cast<int>(m) + 1; };
	int nb = 0;
	if (!streaming) nb = basin_cube->max_value();
	else if (!maximum_basin) nb = static_cast<int>(maxima.size());
	else for (size_t m = 0; m < maxima.size(); m++) nb = std::max(nb, basin_of(m));
	//The overlap matrices ride along on the same points and the same weights as the populations:
	//the density a point contributes is sum_i occ_i phi_i^2, so the diagonal of what is
	//accumulated here sums to exactly the population below and the two can never disagree
	if (ovl && field) ovl = nullptr;
	if (ovl) {
		ovl->mo_index.clear();
		for (int m = 0; m < wavy.get_nmo(); m++)
			if (std::abs(wavy.get_MO_occ(m)) > 1e-8) ovl->mo_index.push_back(m);
		ovl->nmo = static_cast<int>(ovl->mo_index.size());
		ovl->S.assign(nb, vec(ovl->triangle(), 0.0));
		ovl->outside.assign(ovl->triangle(), 0.0);
	}
	vec pop(nb, 0.0);
	volumes.assign(nb, 0.0);
	outside = 0.0;
	double cutoff_outside = 0.0, unresolved_outside = 0.0;
	const int nx = streaming ? 0 : cub->get_size(0), ny = streaming ? 0 : cub->get_size(1), nz = streaming ? 0 : cub->get_size(2);
	//Without a cube there is nothing to read a spacing off, so the trajectory keeps the step of
	//a 0.1 A grid: the integrator is then the one the gridded path has been validated against
	//and only the basin bookkeeping changes
	d3 h{ constants::ang2bohr(0.1), constants::ang2bohr(0.1), constants::ang2bohr(0.1) };
	if (cub)
		for (int i = 0; i < 3; i++)
			h[i] = std::sqrt(cub->get_vector(0, i) * cub->get_vector(0, i) + cub->get_vector(1, i) * cub->get_vector(1, i) + cub->get_vector(2, i) * cub->get_vector(2, i));
	//A third of a voxel with a midpoint step near a nucleus: the Euler step at half a voxel
	//put the N-H boundary of NH3BH3 0.03 e off AIMAll, this is within 0.006. Beyond 1.5 bohr
	//of every nucleus the field is smooth enough for a whole voxel.
	const double voxel = std::min({ h[0], h[1], h[2] });
	const std::vector<atom> atoms = wavy.get_atoms();
	const double sscale = basin_step_scale();
	auto step_at = [&](const d3 &p) {
		double d2 = std::numeric_limits<double>::max();
		for (const atom &at : atoms) {
			const d3 ap = at.get_pos();
			d2 = std::min(d2, std::pow(p[0] - ap[0], 2) + std::pow(p[1] - ap[1], 2) + std::pow(p[2] - ap[2], 2));
		}
		if (d2 < 2.25) return sscale * 0.3 * voxel;
		//The density cutoff bounds the tail outside the basins.
		if (streaming && d2 > 36.0) return sscale * std::min(0.25 * std::sqrt(d2), 4.0);
		return sscale * voxel;
	};
	const double step = 0.3 * voxel;
	//Cube cell of a position and the position within it; false outside the cube
	auto cell = [&](const d3 &p, int *c, d3 &f) {
		if (streaming) return false;
		const int sz[3] = { nx, ny, nz };
		for (int d = 0; d < 3; d++) {
			const double t = (p[d] - cub->get_origin(d)) / cub->get_vector(d, d);
			c[d] = static_cast<int>(std::floor(t));
			if (c[d] < 0 || c[d] >= sz[d] - 1) return false;
			f[d] = t - c[d];
		}
		return true;
	};
	//Basin of the nearest grid point, and whether every voxel within band of the cell agrees.
	//Three voxels: the grid's own boundary can be off by one or two, and a point that close
	//to it is cheap to send up the field
	const int band = 1;
	auto lookup = [&](const d3 &p, bool &settled) {
		int c[3]; d3 f;
		settled = false;
		if (!cell(p, c, f)) return 0;
		const int b = basin_cube->get_value(c[0] + (f[0] > 0.5), c[1] + (f[1] > 0.5), c[2] + (f[2] > 0.5));
		settled = b != 0;
		for (int dx = -band; dx < band + 2 && settled; dx++)
			for (int dy = -band; dy < band + 2 && settled; dy++)
				for (int dz = -band; dz < band + 2; dz++) {
					const int x = c[0] + dx, y = c[1] + dy, z = c[2] + dz;
					if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) continue;
					if (basin_cube->get_value(x, y, z) != b) { settled = false; break; }
				}
		return b;
	};
	//The basin whose maximum lies within reach of p, 0 when none does. Two voxels serve an
	//ELI-D maximum, which is broad; a nucleus gets a tenth of a bohr, since a hydroxyl
	//hydrogen's basin is 0.4 bohr thick and a wider net catches the oxygen's electrons
	const double catch2 = eli_field ? std::pow(2.0 * std::max({ h[0], h[1], h[2] }), 2) : 0.01;
	//Store each attractor centre and squared beta-sphere radius.
	vec beta2(maxima.size(), 0.0);
	std::vector<d3> bcen(maxima.size());
	for (size_t m = 0; m < maxima.size(); m++) bcen[m] = d3{ maxima[m][0], maxima[m][1], maxima[m][2] };
	auto at_maximum = [&](const d3 &p) {
		for (size_t m = 0; m < maxima.size(); m++) {
			if (std::pow(p[0] - maxima[m][0], 2) + std::pow(p[1] - maxima[m][1], 2) + std::pow(p[2] - maxima[m][2], 2) < catch2)
				return basin_of(m);
			if (beta2[m] > 0.0 && std::pow(p[0] - bcen[m][0], 2) + std::pow(p[1] - bcen[m][1], 2) + std::pow(p[2] - bcen[m][2], 2) < beta2[m])
				return basin_of(m);
		}
		return 0;
	};
	//The basin of the maximum nearest p, 0 when the nearest is further than reach. A stalled
	//trajectory has nowhere else to go: the field it was climbing has run out of slope, and the
	//point still has to belong to somebody
	auto nearest_maximum = [&](const d3 &p, const double reach) {
		int best = 0;
		double d2 = reach * reach;
		for (size_t m = 0; m < maxima.size(); m++) {
			const double q = std::pow(p[0] - maxima[m][0], 2) + std::pow(p[1] - maxima[m][1], 2) + std::pow(p[2] - maxima[m][2], 2);
			if (q < d2) { d2 = q; best = basin_of(m); }
		}
		return best;
	};
	//Half the distance to the nearest maximum. A step longer than the validated one must still
	//land inside the catch radius (or the beta sphere) of the attractor it is walking into, since
	//at_maximum only ever looks at where the walk landed. Only a grown step pays for this loop.
	auto reach_limit = [&](const d3 &p) {
		double d2 = std::numeric_limits<double>::max();
		for (size_t m = 0; m < maxima.size(); m++)
			d2 = std::min(d2, std::pow(p[0] - maxima[m][0], 2) + std::pow(p[1] - maxima[m][1], 2) + std::pow(p[2] - maxima[m][2], 2));
		return g_adp_reach * std::sqrt(d2);
	};
	//Evaluate density and gradient in the same orbital pass.
	auto gradient = [&](const d3 &p, d3 &g, double *val = nullptr) {
		if (!eli_field) {
			if (field) { field->grad(p, g); if (val) *val = field->rho(p); }
			else wavy.computeGrad(p, g, val);
			if (core_gradient) { d3 c; (*core_gradient)(p, c); for (int k = 0; k < 3; k++) g[k] += c[k]; }
			if (val && core_density) *val += (*core_density)(p);
			return;
		}
		double e;
		wavy.computeELIGrad(p, e, g);
		if (val) *val = e;
	};
	//The field's value at p and its gradient in one call, which is what the climb needs to see
	//that it has stopped rising: ELI-D's value costs nothing beside its gradient, computeELIGrad
	//building both from the same orbital pass
	auto value_and_gradient = [&](const d3 &p, d3 &g) {
		double v = 0.0;
		gradient(p, g, &v);
		return v;
	};
	//A beta sphere terminates trajectories that cannot leave its attractor basin.
	basin_stage_timer T;
	const std::string fieldname = eli_field ? "ELI-D " : "QTAIM ";
	double margin = 0.0;
	if (streaming && beta_spheres_enabled() && !maxima.empty()) {
		constexpr int ndir = 302;
		static const std::vector<d3> dirs = [] {
			std::vector<d3> d(ndir);
			for (int i = 0; i < ndir; i++) {
				const double z = 1.0 - 2.0 * (i + 0.5) / ndir;
				const double s = std::sqrt(std::max(0.0, 1.0 - z * z));
				const double phi = 2.39996322972865332 * i;   //golden angle: no two samples line up
				d[i] = d3{ s * std::cos(phi), s * std::sin(phi), z };
			}
			return d;
		}();
		const double march = 0.05;   //bohr; the radius is only ever needed to within a step
		const int nm = static_cast<int>(maxima.size());
		//Once, outside the parallel loop: every sphere has to be drawn at the same margin
		margin = basin_beta_margin();
#pragma omp parallel for schedule(dynamic)
		for (int m = 0; m < nm; m++) {
			//Ascend onto the attractor first. A cube maximum is a voxel centre, so half a voxel
			//out on the far side of the true top the radial derivative already points back in:
			//the march stops at its very first sample and the sphere collapses to nothing. That
			//is the whole reason the density gained (its maxima are nuclei, exact to the last
			//digit) and ELI-D gained not one trajectory.
			d3 c{ maxima[m][0], maxima[m][1], maxima[m][2] };
			{
				d3 g;
				double f = value_and_gradient(c, g), sl = 0.5 * voxel;
				for (int it = 0; it < 60 && sl > 1e-4; it++) {
					const double gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
					if (gn < 1e-12) break;
					const d3 t{ c[0] + sl * g[0] / gn, c[1] + sl * g[1] / gn, c[2] + sl * g[2] / gn };
					d3 gt;
					const double ft = value_and_gradient(t, gt);
					if (ft > f) { c = t; f = ft; g = gt; }
					else sl *= 0.5;
				}
			}
			//Walked further than a voxel and a half: that is not this maximum refined any more,
			//it is a different attractor, and a sphere around it would answer for the wrong
			//basin. Such a maximum keeps the plain catch radius and no sphere.
			if (std::pow(c[0] - maxima[m][0], 2) + std::pow(c[1] - maxima[m][1], 2) + std::pow(c[2] - maxima[m][2], 2) > std::pow(1.5 * voxel, 2))
				continue;
			double cap = 3.0;
			for (int n = 0; n < nm; n++) {
				if (n == m || basin_of(n) == basin_of(m)) continue;
				const double d = std::sqrt(std::pow(c[0] - maxima[n][0], 2) + std::pow(c[1] - maxima[n][1], 2) + std::pow(c[2] - maxima[n][2], 2));
				cap = std::min(cap, 0.45 * d);
			}
			double r = cap;
			for (const d3 &u : dirs) {
				double rr = march;
				for (; rr <= cap + 1e-12; rr += march) {
					const d3 q{ c[0] + rr * u[0], c[1] + rr * u[1], c[2] + rr * u[2] };
					d3 g;
					gradient(q, g);
					if (g[0] * u[0] + g[1] * u[1] + g[2] * u[2] >= 0.0) break;
				}
				r = std::min(r, rr - march);
				if (r <= march) break;
			}
			if (r <= 2.0 * march) continue;
			bcen[m] = c;
			beta2[m] = std::pow(margin * r, 2);
		}
	}
	T.lap(fieldname + "beta spheres");
	if (basin_timing_enabled() && margin > 0.0) {
		size_t with = 0;
		for (const double b : beta2) if (b > 0.0) with++;
		std::cout << "  [timing] " << fieldname << "beta spheres: margin " << margin << ", "
			<< with << " of " << beta2.size() << " maxima carry one" << std::endl;
	}
	//Use at least level 3; strongly ionic centres can require level 4 for 0.01 e agreement.
	GridConfiguration config;
	config.accuracy = std::max(accuracy, 3);
	config.alpha_max_scale = static_cast<double>(grid_boost) * grid_boost;
	config.radial_step_scale = grid_boost;
	config.angular_boost = grid_boost - 1;
	config.partition_type = PartitionType::Becke;
	config.no_density_eval = true;
	GridManager grids(config);
	ivec every_atom(wavy.get_ncen());
	std::iota(every_atom.begin(), every_atom.end(), 0);
	grids.setup3DGridsForMolecule(wavy, every_atom);
	const GridData &gd = grids.getGridData();
	T.lap(fieldname + "atomic quadrature grids");
	//The rho = 1e-4 isosurface bounds gridded and streaming basins.
	const double stall_reach = 1e30;
	//The trajectory and quadrature share the density isosurface.
	const double stall_floor = basin_density_cutoff;
	//Resolve a stalled density trajectory from its starting side of the separatrix.
	auto stalled = [&](const d3 &start, const d3 &r, const double gn = -1.0) {
		const double rho = valence(r);
		if (rho < stall_floor) { basin_stall_seen(-1.0, 0.0); return 0; }
		double d2 = std::numeric_limits<double>::max();
		for (size_t m = 0; m < maxima.size(); m++)
			d2 = std::min(d2, std::pow(r[0] - maxima[m][0], 2) + std::pow(r[1] - maxima[m][1], 2) + std::pow(r[2] - maxima[m][2], 2));
		basin_stall_seen(std::sqrt(d2), rho, gn);
		const int b = nearest_maximum(eli_field ? r : start, stall_reach);
		basin_stall_gave_to(b);
		return b;
	};
	auto climb = [&](const d3 &p, long long &lb, long long &ll) {
		if (valence(p) < basin_density_cutoff) return 0;
		bool settled;
		int b = lookup(p, settled);
		//A gridded ELI-D takes a settled cell straight from the cube and leaves the crop outside;
		//streaming has no cube to ask and walks from every point
		if (eli_field && !streaming && (settled || b == 0)) return b;
		int c[3]; d3 f;
		//Off the cube there is nothing to integrate; streaming has no cube to be off
		if (!streaming && !cell(p, c, f)) return 0;
		//Already inside a beta sphere (or on a maximum): the answer needs no trajectory, so it is
		//not counted as one either
		if (const int m0 = at_maximum(p)) return m0;
		lb++;
		d3 r = p, g;
		double last_value = -1.0;
		//Grow only steps whose field direction stays nearly straight.
		double mult = 1.0;
		//Shorten a base step if it would decrease the field.
		double shrink = 1.0;
		//Was the step that reached r actually longer than the floor? mult is raised at the end of a
		//step, so mult > 1 at the top of the next iteration says "the next step may be grown", not
		//"the last one was" - and only the latter is grounds for throwing a point away.
		bool grown_last = false;
		d3 r_prev = p;
		const bool grow = basin_adaptive_step_enabled();
		double value_prev = -1.0;
		for (int s = 0; s < g_step_cap; s++) {
			if (eli_field && valence(r) < basin_density_cutoff) return 0;
			if (grow) adp_count(g_adp_steps);
			const int m = at_maximum(r);
			if (m) return m;
			const double floor_step = step_at(r);
			//The gridded ELI-D climb is steered by the cube below and never asked whether it is
			//still rising; every streaming walk is, since nothing else can stop it
			double here = last_value;
			if (!eli_field || streaming) {
				here = value_and_gradient(r, g);
				if (here <= last_value) {
					//Retry a grown step at the base length after excessive turning.
					if (grown_last) { adp_count(g_adp_fall); r = r_prev; last_value = value_prev; mult = 1.0; grown_last = false; continue; }
					mult = 1.0;
					if (!eli_field && g_step_shrink && shrink > 0.0625) {
						adp_count(g_adp_shrink);
						shrink *= 0.5;
						r = r_prev;
						last_value = value_prev;
						grown_last = false;
						continue;
					}
					const int n = stalled(p, r, std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]));
					if (n) return n;
					break;
				}
			}
			else gradient(r, g);
			double gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
			if (gn < 1e-12) {
				//A critical point, and now the only way to reach the stall path: 3716 of CCH's 3716
				//old stalls carried a gradient of 1e-2 or more, so none of them was ever this. A
				//streaming density walk has no cube behind it, so breaking here reports the point
				//outside every basin - 0.0102 e of CCH - where the tie-break has a definite answer.
				if (!eli_field) { const int n = stalled(p, r, gn); if (n) return n; }
				break;
			}
			//Preserve the midpoint expression order used for validated basin populations.
			const d3 g0{ g[0], g[1], g[2] };
			const double gn0 = gn;
			const d3 dir{ g0[0] / gn0, g0[1] / gn0, g0[2] / gn0 };
			//Retry a step at the base length when its field direction turns too far.
			const double base = floor_step * shrink;
			double sl = base;
			double cosine = 0.0;
			bool stepped = false;
			for (int attempt = 0; attempt < 2 && !stepped; attempt++) {
				sl = mult > 1.0 ? std::min(base * mult, std::max(base, reach_limit(r))) : base;
				if (mult > 1.0) adp_count(g_adp_tries);
				d3 mid;
				for (int k = 0; k < 3; k++) mid[k] = r[k] + 0.5 * sl * g0[k] / gn0;
				gradient(mid, g);
				gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
				if (gn < 1e-12) break;
				cosine = (g[0] * dir[0] + g[1] * dir[1] + g[2] * dir[2]) / gn;
				if (mult > 1.0 && cosine < g_adp_keep) { adp_count(g_adp_turn); mult = 1.0; continue; }
				stepped = true;
			}
			//Same case one level in: the attempt loop gives up only when the midpoint gradient
			//vanished too.
			if (!stepped) {
				if (!eli_field) { const int n = stalled(p, r, gn); if (n) return n; }
				break;
			}
			r_prev = r;
			value_prev = last_value;
			last_value = here;
			grown_last = sl > base;
			for (int k = 0; k < 3; k++) r[k] += sl * g[k] / gn;
			if (shrink < 1.0) shrink = std::min(1.0, shrink * g_step_relax);
			else if (cosine > g_adp_grow && mult < g_adp_cap && grow) mult *= 2.0;
			if (!streaming) {
				const int b2 = lookup(r, settled);
				if (b2 == 0 && eli_field) { ll++; break; }
				if (b2) b = b2;
				if (settled && eli_field) break;
			}
		}
		//Assign a stalled streaming trajectory by its starting point.
		if (streaming) { adp_count(g_adp_exhaust); return stalled(p, r); }
		return b;
	};
	//Bisect a radial cell only when its edge trajectories reach different basins.
	constexpr int bisections = 6;
	long long boundary_points = 0, lost = 0;
	for (size_t a = 0; a < gd.atomic_grids.size(); a++) {
		const vec &X = gd.atomic_grids[a][GridData::X], &Y = gd.atomic_grids[a][GridData::Y], &Z = gd.atomic_grids[a][GridData::Z], &W = gd.atomic_grids[a][GridData::BECKE_WEIGHT];
		const int np = static_cast<int>(X.size());
		const d3 centre = atoms[a].get_pos();
		vec radius(np), shells;
		for (int i = 0; i < np; i++) radius[i] = std::sqrt(std::pow(X[i] - centre[0], 2) + std::pow(Y[i] - centre[1], 2) + std::pow(Z[i] - centre[2], 2));
		shells = radius;
		std::sort(shells.begin(), shells.end());
		shells.erase(std::unique(shells.begin(), shells.end(), [](double x, double y) { return std::abs(x - y) < 1e-8; }), shells.end());
		//The cell of point i as a shell segment: the two radii halfway to the neighbouring
		//shells. False for a cell with no radial extent, which is given to one basin whole.
		auto cell_edges = [&](const int i, double &in, double &out) {
			const size_t k = std::lower_bound(shells.begin(), shells.end(), radius[i] - 1e-8) - shells.begin();
			in = k > 0 ? 0.5 * (shells[k - 1] + shells[k]) : 0.0;
			out = k + 1 < shells.size() ? 0.5 * (shells[k] + shells[k + 1]) : shells[k];
			return radius[i] > 1e-8 && out > in + 1e-8;
		};
		//A point at radius r along point i's direction
		auto along_i = [&](const int i, const double r) {
			const double f = r / radius[i];
			return d3{ centre[0] + (X[i] - centre[0]) * f, centre[1] + (Y[i] - centre[1]) * f, centre[2] + (Z[i] - centre[2]) * f };
		};
		//Reuse probes only for adjacent shells with matching directions.
		ivec partner(np, -1);
		{
			ivec start;
			for (int i = 0; i < np; i++)
				if (i == 0 || std::abs(radius[i] - radius[i - 1]) > 1e-8) start.push_back(i);
			start.push_back(np);
			for (size_t s = 1; s + 1 < start.size(); s++) {
				const int a0 = start[s - 1], a1 = start[s], a2 = start[s + 1];
				if (a1 - a0 != a2 - a1) continue;
				//Only a cell whose inner edge is the previous cell's outer edge, which needs the
				//two shells to be neighbours in shells[] - a pair closer than the 1e-8 the unique
				//pass merges on shares one entry and has no edge between them
				double i0, o0, i1, o1;
				if (!cell_edges(a0, i0, o0) || !cell_edges(a1, i1, o1) || o0 != i1) continue;
				for (int j = 0; j < a1 - a0; j++) {
					const int lo = a0 + j, hi = a1 + j;
					const double dot = ((X[lo] - centre[0]) * (X[hi] - centre[0]) + (Y[lo] - centre[1]) * (Y[hi] - centre[1]) + (Z[lo] - centre[2]) * (Z[hi] - centre[2])) / (radius[lo] * radius[hi]);
					if (dot > 1.0 - 1e-12) partner[hi] = lo;
				}
			}
		}
		//Every cell's outer probe, climbed once. A gridded ELI-D takes most cells from the cube
		//without probing at all, so it keeps the old on-demand path and this pass is skipped.
		ivec outer_probe(np, -1);
		if (!(eli_field && !streaming)) {
#pragma omp parallel
			{
				long long lb = 0, ll = 0;
#pragma omp for schedule(dynamic, 16)
				for (int i = 0; i < np; i++) {
					double in, out;
					if (W[i] == 0.0 || !cell_edges(i, in, out)) continue;
					outer_probe[i] = climb(along_i(i, out), lb, ll);
				}
#pragma omp critical
				{ boundary_points += lb; lost += ll; }
			}
		}
#pragma omp parallel
		{
			vec lp(nb, 0.0), lv(nb, 0.0);
			double lo = 0.0, lc = 0.0, lu = 0.0;
			long long lb = 0, ll = 0;
			//ponytail: one triangle per basin per thread, nb * nmo^2 / 2 doubles each; the caller
			//sizes the job, a molecule big enough to hurt here has other limits first
			vec2 ls;
			vec ls_out;
			vec phi;
			vec2 dbuf;
			if (ovl) {
				ls.assign(nb, vec(ovl->triangle(), 0.0));
				ls_out.assign(ovl->triangle(), 0.0);
				phi.resize(wavy.get_nmo(), 0.0);
				dbuf.assign(wavy.get_ncen(), vec(16, 0.0));
			}
			//Rank-1 update of one basin's triangle: the point's share of the quadrature weight
			//times the outer product of the orbitals it sees
			auto accumulate = [&](const int b, const double wq) {
				if (!ovl || wq == 0.0) return;
				double *Sb = b > 0 ? ls[b - 1].data() : ls_out.data();
				const int *idx = ovl->mo_index.data();
				for (int a2 = 0; a2 < ovl->nmo; a2++) {
					const double pa = wq * phi[idx[a2]];
					if (pa == 0.0) continue;
					double *row = Sb + (size_t)a2 * (a2 + 1) / 2;
					for (int b2 = 0; b2 <= a2; b2++) row[b2] += pa * phi[idx[b2]];
				}
			};
#pragma omp for schedule(dynamic, 16)
			for (int i = 0; i < np; i++) {
				const double w = W[i];
				if (w == 0.0) continue;
				const d3 p{ X[i], Y[i], Z[i] };
				bool settled;
				int b = lookup(p, settled);
				//The same density, taken from the orbital pass that also hands out phi
				const double rho = ovl ? wavy.compute_dens(p, dbuf, phi) : valence(p);
				//A basin's share of the cell's quadrature weight; the weight itself stays with
				//the rule, only who gets it is decided here
				auto give = [&](const int bb, const double fr) {
					if (fr <= 0.0) return;
					if (bb == 0) {
						lo += w * rho * fr;
						if (rho < basin_density_cutoff) lc += w * rho * fr; else lu += w * rho * fr;
						accumulate(0, w * fr);
						return;
					}
					lp[bb - 1] += w * rho * fr;
					lv[bb - 1] += w * fr;
					accumulate(bb, w * fr);
				};
				if (rho < basin_density_cutoff) { give(0, 1.0); continue; }
				//For a gridded ELI-D a cell whose neighbourhood agrees is taken from the grid, as
				//before; streaming has no grid to take it from and every cell is refined
				if (eli_field && !streaming && (b == 0 || settled)) { give(b, 1.0); continue; }
				//Climb the cell centre only when the edge probes do not decide its basin.
				int bc = -1;
				auto centre_basin = [&]() { if (bc < 0) bc = climb(p, lb, ll); return bc; };
				double inner, outer;
				if (!cell_edges(i, inner, outer)) { give(centre_basin(), 1.0); continue; }
				auto along = [&](const double r) { return along_i(i, r); };
				//A failed edge probe defers to the cell centre.
				const int pin = partner[i];
				int bi = pin >= 0 && outer_probe[pin] >= 0 ? outer_probe[pin] : climb(along(inner), lb, ll);
				int bo = outer_probe[i] >= 0 ? outer_probe[i] : climb(along(outer), lb, ll);
				if (bi == 0) bi = centre_basin();
				if (bo == 0) bo = centre_basin();
				if (bi == bo) { give(bi, 1.0); continue; }
				//ponytail: one crossing per cell. Three basins meeting inside a single quadrature
				//cell is a smaller thing than the rule's own error; bisect for more if it is not
				double lo_r = inner, hi_r = outer;
				for (int it = 0; it < bisections; it++) {
					const double mid = 0.5 * (lo_r + hi_r);
					int bm = climb(along(mid), lb, ll);
					if (bm == 0) bm = centre_basin();
					if (bm == bi) lo_r = mid; else hi_r = mid;
				}
				const double rc = 0.5 * (lo_r + hi_r);
				//Each side gets the density it carries over its own part of the shell segment,
				//whose volume goes as r^3
				const double wi = valence(along(0.5 * (inner + rc))) * (rc * rc * rc - inner * inner * inner);
				const double wo = valence(along(0.5 * (rc + outer))) * (outer * outer * outer - rc * rc * rc);
				const double sum = wi + wo;
				if (sum > 0.0) { give(bi, wi / sum); give(bo, wo / sum); }
				else give(centre_basin(), 1.0);
			}
#pragma omp critical
			{
				for (int b = 0; b < nb; b++) { pop[b] += lp[b]; volumes[b] += lv[b]; }
				if (ovl) {
					for (int b = 0; b < nb; b++)
						for (size_t t = 0; t < ovl->S[b].size(); t++) ovl->S[b][t] += ls[b][t];
					for (size_t t = 0; t < ovl->outside.size(); t++) ovl->outside[t] += ls_out[t];
				}
				outside += lo;
				cutoff_outside += lc;
				unresolved_outside += lu;
				boundary_points += lb;
				lost += ll;
			}
		}
	}
	if (core_density)
		for (int a = 0; a < wavy.get_ncen(); a++) {
			const int ncore = wavy.get_atom_ECP_electrons(a);
			if (ncore <= 0) continue;
			bool settled;
			//Streaming has no basin cube to read the nucleus out of; it is one of the maxima by
			//construction, so the maximum it sits on is its basin
			const int b = streaming ? at_maximum(wavy.get_atom_pos(a)) : lookup(wavy.get_atom_pos(a), settled);
			if (b > 0) pop[b - 1] += ncore;
			else outside += ncore;
		}
	T.lap(fieldname + "point loop");
	{
		//Read and reset unconditionally so the counts never carry from one field or one molecule
		//into the next, but print only under -basin_timing: the golden files are line-by-line
		//captures of this console log, so an extra line here shifts every basin row below it.
		long long sv, sf, sfar; double srho;
		basin_stall_counters(sv, sf, sfar, srho);
		long long sgn[4], sto[g_stall_basins + 1];
		for (int i = 0; i < 4; i++) sgn[i] = g_stall_gn[i].load();
		for (int i = 0; i <= g_stall_basins; i++) sto[i] = g_stall_to[i].load();
		const long long sx = g_adp_exhaust.exchange(0);
		basin_stall_counters_reset();
		//Built in its own stream on purpose. Scientific notation is the only readable form for a
		//density that spans 1e-6 to 1e-1, and setting it on std::cout leaves the precision behind
		//for whatever prints next - the basin table is three lines below and the golden files are
		//captures of it.
		if (g_basin_timing && sf + sv + sx > 0) {
			std::ostringstream line;
			line << std::scientific << std::setprecision(3);
			line << "  " << fieldname << "trajectories that stalled: " << sf << " in the field, "
				<< sv << " below " << stall_floor << " e/bohr^3; " << sfar
				<< " further than a bohr from every attractor, the densest at " << srho << " e/bohr^3";
			line << std::defaultfloat << "\n  " << fieldname << "stall gradients: ";
			const char *edge[4] = { "<1e-6", "<1e-4", "<1e-2", ">=1e-2" };
			for (int i = 0; i < 4; i++) line << (i ? ", " : "") << sgn[i] << " " << edge[i];
			line << "; " << sx << " ran the step budget out; handed to basin";
			for (int i = 0; i <= g_stall_basins; i++)
				if (sto[i]) line << " " << (i == g_stall_basins ? std::string("other") : std::to_string(i)) << ":" << sto[i];
			std::cout << line.str() << std::endl;
		}
	}
	if (g_basin_timing) std::cout << "  [timing] " << fieldname << "outside the density isosurface: " << cutoff_outside << " e, unresolved inside: " << unresolved_outside << " e" << std::endl;
	if (g_basin_timing && g_adaptive_step) {
		const long long st = g_adp_steps.exchange(0), tr = g_adp_tries.exchange(0);
		const long long tu = g_adp_turn.exchange(0), fa = g_adp_fall.exchange(0);
		const long long sh = g_adp_shrink.exchange(0);
		std::cout << "  [timing] " << fieldname << "grown step: " << st << " steps, " << tr
			<< " proposals, " << tu << " turned back, " << fa << " fell back, " << sh << " shrunk below the floor, "
			<< std::fixed << std::setprecision(1)
			<< (tr ? 100.0 * static_cast<double>(tu + fa) / static_cast<double>(tr) : 0.0)
			<< " % of proposals wasted" << std::endl;
		//A step count without the knobs it was taken at is not a measurement of anything
		std::cout << "  [timing] " << fieldname << "grown step knobs: cap " << std::setprecision(4) << g_adp_cap
			<< ", grow " << std::setprecision(8) << g_adp_grow << ", keep " << g_adp_keep
			<< ", reach " << std::setprecision(4) << g_adp_reach
			<< ", exp cutoff " << constants::exp_cutoff << std::endl;
	}
	std::cout << "Quadrature points sent along the field: " << boundary_points << ", left the grid: " << lost << std::endl;
	return pop;
}

//Partition the determinant pair density over basin overlap matrices.
delocalization_result delocalization_indices(const WFN &wavy, const basin_overlaps &ovl)
{
	delocalization_result r;
	const int nb = static_cast<int>(ovl.S.size());
	const int n = ovl.nmo;
	r.lambda.assign(nb, 0.0);
	r.population.assign(nb, 0.0);
	r.outside_half.assign(nb, 0.0);
	if (nb == 0 || n == 0) return r;
	//Infer spin resolution from orbital occupations.
	double max_occ = 0.0;
	for (int i = 0; i < n; i++) max_occ = std::max(max_occ, std::abs(wavy.get_MO_occ(ovl.mo_index[i])));
	const bool restricted = max_occ > 1.0 + 1e-6;
	const double m = restricted ? 2.0 : 1.0;
	vec occ(n);
	ivec spin(n);
	for (int i = 0; i < n; i++) {
		occ[i] = wavy.get_MO_occ(ovl.mo_index[i]) / m;
		//Spatial orbitals stand for both spins and all of them exchange with one another; only
		//a genuinely spin-resolved set has an alpha and a beta block that must not mix
		spin[i] = restricted ? 0 : wavy.get_MO_op(ovl.mo_index[i]);
	}
	for (int b = 0; b < nb; b++)
		for (int i = 0; i < n; i++)
			r.population[b] += m * occ[i] * ovl.at(b, i, i);
	//The basins and outside region together resolve the orbital metric.
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= i; j++) {
			if (spin[i] != spin[j]) continue;
			double s = 0.0;
			for (int b = 0; b < nb; b++) s += ovl.at(b, i, j);
			if (!ovl.outside.empty()) s += ovl.outside[basin_overlaps::packed(i, j)];
			r.identity_error = std::max(r.identity_error, std::abs(s - (i == j ? 1.0 : 0.0)));
		}
	//The pair sum is symmetric in i and j, so the triangle is taken once and doubled off the
	//diagonal
	auto pair_sum = [&](const vec &a, const vec &b) {
		double s = 0.0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j <= i; j++) {
				if (spin[i] != spin[j]) continue;
				const size_t k = basin_overlaps::packed(i, j);
				const double t = occ[i] * occ[j] * a[k] * b[k];
				s += i == j ? t : 2.0 * t;
			}
		return s;
	};
	for (int b = 0; b < nb; b++) {
		r.lambda[b] = m * pair_sum(ovl.S[b], ovl.S[b]);
		if (!ovl.outside.empty()) r.outside_half[b] = m * pair_sum(ovl.S[b], ovl.outside);
	}
	for (int a = 0; a < nb; a++)
		for (int b = a + 1; b < nb; b++) {
			r.pairs.push_back({ a, b });
			r.di.push_back(2.0 * m * pair_sum(ovl.S[a], ovl.S[b]));
		}
	return r;
}

void report_delocalization(const WFN &wavy, const basin_overlaps &ovl, const svec &labels, std::ostream &log, const double threshold)
{
	const delocalization_result r = delocalization_indices(wavy, ovl);
	const int nb = static_cast<int>(r.lambda.size());
	if (nb == 0 || ovl.nmo == 0) return;
	auto name = [&](const int b) { return b < static_cast<int>(labels.size()) ? labels[b] : std::to_string(b + 1); };
	log << "\nDelocalization indices (" << ovl.nmo << " occupied orbitals):\n";
	citations::cite(citations::Method::LIDI, log);
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

vec integrate_values_in_basins(const cube *cub, const cubei *basin_cube, svec& basin_label, bool debug)
{
	const int basin_count = basin_cube->max_value();
	vec EDS(basin_count);
	vec VOL(basin_count);
	double dv = cub->get_dv();
	err_checkf(basin_cube->get_size(0) == cub->get_size(0) && basin_cube->get_size(1) == cub->get_size(1) && basin_cube->get_size(2) == cub->get_size(2), "Basin cube and original cube must have the same dimensions!", std::cout);
	err_checkf(basin_cube->get_dv() - cub->get_dv() < 1E-10, "Basin must have reasonable size!", std::cout);
	if (debug)
		std::cout << "dv: " << dv << " iCP: " << basin_count << std::endl;

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
	for (int x = 0; x < cub->get_size(0); x++)
		for (int y = 0; y < cub->get_size(1); y++)
			for (int z = 0; z < cub->get_size(2); z++) {
				const int basin_index = basin_cube->get_value(x, y, z) - 1;
				if (basin_index < 0 || basin_index >= basin_count)
					continue;
				const double contribution = cub->get_value(x, y, z) * dv;
#ifdef _OPENMP
#pragma omp atomic
#endif
				EDS[basin_index] += contribution; // Value in this voxel
#ifdef _OPENMP
#pragma omp atomic
#endif
				VOL[basin_index] += dv; // Volume of this voxel
			}
	for (int a = 0; a < basin_count; a++) {
		if (EDS[a] > 0.001)
			std::cout << "basin: " << std::setw(4) << a << " label: " << std::setw(16) << basin_label[a] << " integrated value: " << std::setw(8) << std::setprecision(4) << std::fixed << EDS[a] << " volume:" << std::setw(10) << VOL[a] << "\n";
	}
	for (int a = 0; a < basin_count; a++) {
		if (EDS[a] < 0.001)
			std::cout << "WARNING: Integrated value in basin " << a + 1 << " is very low (" << EDS[a] << ") and might be inaccurate due to numerical errors!\n";
	}
	return EDS;
};

svec assign_labels_to_basins(const std::vector<d4> &Maxima, const std::vector<atom> &atoms, bool debug, int type_switch)
{
	using namespace std;
	svec result(Maxima.size());
	//type switch defines the field that was analyzed and will define how we assign labels to the basins
	switch (type_switch) {
		case 0: //QTAIM case
		for (size_t i = 0; i < Maxima.size(); i++) {
			const d3 pos = {Maxima[i][0], Maxima[i][1], Maxima[i][2]};
			double min_dist = std::numeric_limits<double>::max();
			int atom_index = -1;
			for (size_t j = 0; j < atoms.size(); j++) {
				const d3 apos = atoms[j].get_pos();
				const double dx = pos[0] - apos[0];
				const double dy = pos[1] - apos[1];
				const double dz = pos[2] - apos[2];
				const double dist2 = dx * dx + dy * dy + dz * dz;
				if (dist2 < min_dist) {
					min_dist = dist2;
					atom_index = j;
				}
			}
			err_checkf(atom_index >= 0, "No atom found for basin " + toString<size_t>(i) + " at position (" + toString<double>(pos[0]) + ", " + toString<double>(pos[1]) + ", " + toString<double>(pos[2]) + ")!", std::cout);
			result[i] = atoms[atom_index].get_label() + to_string(atom_index);
			//A maximum off every nucleus is a non-nuclear attractor, real or a grid bump on a
			//flat bond; it must not carry the atom's label into the charge column
			if (min_dist > 0.25) result[i] = "NNA near " + result[i];
		}
		break;
		case 1: //ELI case, assign core basins based on atom within, valence basins as the connected basin and bonds as the two closest atoms.
			for (size_t i = 0; i < Maxima.size(); i++) {
				const d3 pos = {Maxima[i][0], Maxima[i][1], Maxima[i][2]};
				double min_dist1 = std::numeric_limits<double>::max();
				double min_dist2 = std::numeric_limits<double>::max();
				int atom_index1 = -1;
				int atom_index2 = -1;
				for (size_t j = 0; j < atoms.size(); j++) {
					const d3 apos = atoms[j].get_pos();
					const double dx = pos[0] - apos[0];
					const double dy = pos[1] - apos[1];
					const double dz = pos[2] - apos[2];
					const double dist2 = dx * dx + dy * dy + dz * dz;
					if (dist2 < min_dist1) {
						min_dist2 = min_dist1;
						atom_index2 = atom_index1;
						min_dist1 = dist2;
						atom_index1 = j;
					}
					else if (dist2 < min_dist2) {
						min_dist2 = dist2;
						atom_index2 = j;
					}
				}
				err_checkf(atom_index1 >= 0, "No atom found for basin " + toString<size_t>(i) + " at position (" + toString<double>(pos[0]) + ", " + toString<double>(pos[1]) + ", " + toString<double>(pos[2]) + ")!", std::cout);
				//the charge is read after the index is known to be one: atoms[-1] was being indexed to
				//compute core_dist one line above the check that atom_index1 exists at all
				const double core_dist = std::pow(core_shell_radius(atoms[atom_index1].get_charge()), 2);
				//Assign every basin of a monoatomic system to its nucleus.
				const bool lone_atom = atom_index2 < 0;
				const double ratio = lone_atom ? 0.0 : std::max(1e-5, min_dist1) / std::max(1e-5, min_dist2);
				//Label a bond basin only for a bonded pair that brackets its maximum.
				bool between = false;
				err_checkf(atom_index2 >= 0 || atoms.size() == 1, "Only one atom found for basin " + toString<size_t>(i) + " at position (" + toString<double>(pos[0]) + ", " + toString<double>(pos[1]) + ", " + toString<double>(pos[2]) + ")!", std::cout);
				if (!lone_atom) {
					const d3 p1 = atoms[atom_index1].get_pos();
					const d3 p2 = atoms[atom_index2].get_pos();
					const int z1 = atoms[atom_index1].get_charge();
					const int z2 = atoms[atom_index2].get_charge();
					const double r1 = (z1 > 0 && z1 < 114) ? constants::covalent_radii[z1] : 1.5;
					const double r2 = (z2 > 0 && z2 < 114) ? constants::covalent_radii[z2] : 1.5;
					const double dAB = std::sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));
					between = dAB <= constants::ang2bohr(1.3 * (r1 + r2)) &&
						(std::sqrt(min_dist1) + std::sqrt(min_dist2)) <= 1.25 * std::max(1e-5, dAB);
				}
				if (atoms[atom_index1].get_charge() == 1 && min_dist1 < 0.36) // The basin holding a proton: its maximum sits within 0.6 bohr of the nucleus
					result[i] = atoms[atom_index1].get_label() + to_string(atom_index1);
				else if (lone_atom) // one atom in the molecule: inside its core radius a core shell, outside it the valence shell
					result[i] = atoms[atom_index1].get_label() + to_string(atom_index1) + (min_dist1 < core_dist ? " core" : " LP");
				else if (min_dist1 < core_dist && atoms[atom_index1].get_charge() > 2) // If the maximum is very close to an atom, we assume it's a core basin and label it with that atom
					result[i] = atoms[atom_index1].get_label() + to_string(atom_index1) + " core";
				else if ((ratio < 0.333 || ratio > 3) && atoms[atom_index1].get_charge() > 2) // If the maximum is significantly closer to one atom than to the other, we assume it's a valence basin and label it with the closest atom
					result[i] = atoms[atom_index1].get_label() + to_string(atom_index1) + " LP";
				else if (between) // between its two nearest nuclei: a bond basin, labelled with both atoms
					result[i] = atoms[atom_index1].get_label() + to_string(atom_index1) + "-" + atoms[atom_index2].get_label() + to_string(atom_index2) + " bond";
				else // not between them: it belongs to the nearest atom alone, core inside the core radius
					result[i] = atoms[atom_index1].get_label() + to_string(atom_index1) + (min_dist1 < core_dist && atoms[atom_index1].get_charge() > 2 ? " core" : " LP");
			}
			break;
		default:
			err_not_impl_f("Label assignment type " + toString<int>(type_switch) + " is not implemented!", std::cout);
	}
	return result;
}
