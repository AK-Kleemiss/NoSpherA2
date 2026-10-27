#include "pch.h"
#include "structure_factors.h"
#include "AtomGrid.h"
#include "GridManager.h"
#include "spherical_density.h"
#include "constants.h"
#include "tuning.h"

//TFVC, MBIS and EMBIS shares on XCW's grid from the density the SCF carries, see structure_factors::xcw_grid. The fits are
//the stand-alone ones, make_chi, make_MBIS_vectors and make_EMBIS_tensors, on build_I's points, each MBIS/EMBIS fit
//starting from the last. Their response to rho is left out of the perturbation: within an evaluation the shares are fixed
struct xcw_weights {
	PartitionType part = PartitionType::Hirshfeld;
	WFN atoms;
	vec3 grid; //[atom][GridData::GridIndex] what the fits read: points, w0 and rho
	ivec np;
	occ::qm::AOBasis ao;
	vec x, y, z, chi;
	ivec Z;
	std::vector<std::pair<vec, vec>> mbis;
	std::vector<std::pair<vec2, vec>> embis;
	int fits = 0, sweeps = 0;
	double ms = 0;
	std::string said;
};

void structure_factors::setup_grid_weights(const occ::qm::AOBasis& ao)
{
	const int ng = static_cast<int>(XG.w0.size());
	XG.ws = std::make_shared<xcw_weights>();
	xcw_weights& W = *XG.ws;
	W.part = opt->partition_type;
	W.ao = ao;
	W.x.resize(ng); W.y.resize(ng); W.z.resize(ng); W.Z.resize(ng);
	W.grid.assign(ng, vec2(GridData::GridIndex::WFN_DENSITY + 1));
	W.np.resize(ng);
	for (int g = 0; g < ng; g++) {
		atom a;
		for (int i = 0; i < 3; i++) a.set_coordinate(i, asym_atoms[g].pos[i]);
		a.set_charge(asym_atoms[g].type);
		W.atoms.push_back_atom(a);
		W.x[g] = asym_atoms[g].pos[0]; W.y[g] = asym_atoms[g].pos[1]; W.z[g] = asym_atoms[g].pos[2]; W.Z[g] = asym_atoms[g].type;
		const int n = W.np[g] = static_cast<int>(XG.w0[g].size());
		vec2& G = W.grid[g];
		G[GridData::GridIndex::X].resize(n); G[GridData::GridIndex::Y].resize(n); G[GridData::GridIndex::Z].resize(n);
		for (int p = 0; p < n; p++) {
			G[GridData::GridIndex::X][p] = W.x[g] + XG.d1[g][p];
			G[GridData::GridIndex::Y][p] = W.y[g] + XG.d2[g][p];
			G[GridData::GridIndex::Z][p] = W.z[g] + XG.d3[g][p];
		}
		G[GridData::GridIndex::WEIGHT] = XG.w0[g];
		G[GridData::GridIndex::WFN_DENSITY].resize(n);
	}
	XG.rho_s.assign(ng, vec());
	std::cout << "XCW " << (W.part == PartitionType::TFVC ? "TFVC" : W.part == PartitionType::MBIS ? "MBIS" : "EMBIS")
		<< " shares on the grid follow the density; Hirshfeld's until the first one" << std::endl;
}

void structure_factors::update_grid_weights(const vec2& rho, const dMatrix2& D)
{
	xcw_weights& W = *XG.ws;
	const auto t0 = std::chrono::steady_clock::now();
	const int ng = static_cast<int>(W.np.size());
	const char* tl = tuning("NOSPHERA2_XCW_WEIGHT_TOL");
	const double tol = tl ? std::atof(tl) : 1e-6;
	if (W.part == PartitionType::TFVC) {
		//chi from the line densities of D itself, rho = 2 chi^T D chi as grid_rho has it
		const int n = static_cast<int>(D.extent(0));
		occ::Mat Dm(n, n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) Dm(i, j) = D(i, j);
		const DensityBatch line = [&](const int m, const double* x, const double* y, const double* z, double* out) {
			occ::Mat3N p(3, m);
			for (int i = 0; i < m; i++) { p(0, i) = x[i]; p(1, i) = y[i]; p(2, i) = z[i]; }
			const occ::Mat phi = occ::gto::evaluate_basis(W.ao, p, 0).phi;
			const occ::Mat t = phi * Dm;
			for (int i = 0; i < m; i++) out[i] = 2.0 * t.row(i).dot(phi.row(i));
		};
		//make_chi's warnings once each, not at every fit
		std::stringstream said;
		std::streambuf* const old = std::cout.rdbuf(said.rdbuf());
		W.chi = make_chi(W.atoms, 40, true, false, line);
		std::cout.rdbuf(old);
		if (said.str() != W.said) std::cout << said.str();
		W.said = said.str();
		vec dab(static_cast<size_t>(ng) * ng);
		for (int a = 0; a < ng; a++)
			for (int b = 0; b < ng; b++) dab[static_cast<size_t>(a) * ng + b] = std::hypot(W.x[b] - W.x[a], W.y[b] - W.y[a], W.z[b] - W.z[a]);
#pragma omp parallel
		{
			vec pa_b(ng), pa_tv(ng);
			for (int g = 0; g < ng; g++) {
				const vec2& G = W.grid[g];
#pragma omp for schedule(dynamic, 256) nowait
				for (int p = 0; p < W.np[g]; p++)
					XG.w[g][p] = XG.w0[g][p] * get_integration_weights(ng, W.Z.data(), W.x.data(), W.y.data(), W.z.data(), g,
						G[GridData::GridIndex::X][p], G[GridData::GridIndex::Y][p], G[GridData::GridIndex::Z][p], pa_b, pa_tv, W.chi, dab.data())[1];
			}
		}
	}
	else {
		for (int g = 0; g < ng; g++)
			for (int p = 0; p < W.np[g]; p++) W.grid[g][GridData::GridIndex::WFN_DENSITY][p] = rho[g][p] / XG.w0[g][p];
		std::ostream none(nullptr);
		int k = 0;
		if (W.part == PartitionType::MBIS || W.embis.empty()) {
			W.mbis = make_MBIS_vectors(W.atoms, W.grid, W.np, false, none, W.mbis, tol, &k);
			W.sweeps += k;
		}
		if (W.part == PartitionType::EMBIS) {
			W.embis = make_EMBIS_tensors(W.atoms, W.grid, W.np, false, W.mbis, none, W.embis, tol, &k);
			W.sweeps += k;
		}
		std::vector<MBIS_Atom> m;
		std::vector<EMBIS_Atom> e;
		for (int a = 0; a < ng; a++)
			if (W.part == PartitionType::MBIS) m.emplace_back(W.Z[a], W.mbis[a].first, W.mbis[a].second);
			else e.emplace_back(W.Z[a], W.embis[a].first, W.embis[a].second);
		//the share of g at its own points, rho_g / Sum_a rho_a of the fitted pro-atoms out to the fits' reach
#pragma omp parallel
		for (int g = 0; g < ng; g++) {
			const vec2& G = W.grid[g];
#pragma omp for schedule(static) nowait
			for (int p = 0; p < W.np[g]; p++) {
				double own = 0, all = 0;
				for (int a = 0; a < ng; a++) {
					const d3 d = { G[GridData::GridIndex::X][p] - W.x[a], G[GridData::GridIndex::Y][p] - W.y[a], G[GridData::GridIndex::Z][p] - W.z[a] };
					const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
					if (r2 > constants::far_away_sq) continue;
					const double v = W.part == PartitionType::MBIS ? m[a].get_radial_density(std::sqrt(r2)) : e[a].get_density(d);
					all += v;
					if (a == g) own = v;
				}
				XG.w[g][p] = all > 0.0 ? XG.w0[g][p] * own / all : 0.0;
			}
		}
	}
	W.fits++;
	W.ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

std::string structure_factors::grid_weight_report()
{
	if (!XG.ws) return "";
	xcw_weights& W = *XG.ws;
	std::ostringstream s;
	s << (W.part == PartitionType::TFVC ? "TFVC" : W.part == PartitionType::MBIS ? "MBIS" : "EMBIS") << " shares: " << W.fits << " fits";
	if (W.part != PartitionType::TFVC) s << ", " << W.sweeps << " sweeps";
	s << ", " << std::fixed << std::setprecision(2) << W.ms / 1000.0 << " s";
	W.fits = W.sweeps = 0;
	W.ms = 0;
	return s.str();
}
