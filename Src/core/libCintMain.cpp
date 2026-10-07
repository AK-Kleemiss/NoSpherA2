#include "pch.h"
#include "constants.h"
#include "libCintMain.h"
#include "libCintKernels.h"
#include "tuning.h"
//
#if defined(__APPLE__)
// On macOS we are using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#elif !defined(NSA2_OPENBLAS)
// Linux/Windows with oneMKL
#include <mkl.h>
#endif

// Function to compute three-center two-electron integrals (eri3c)
template <typename Kernel>
void computeEri3c(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c)
{
	int nQM = param1.get_nbas();
	int nAux = param2.get_nbas();

	Int_Params combined(param1, param2);
	// combined.print_data("combined");

	ivec bas = combined.get_bas();
	ivec atm = combined.get_atm();
	vec env = combined.get_env();

	ivec shl_slice = {
		0,
		nQM,
		0,
		nQM,
		nQM,
		nQM + nAux,
	};

	int nat = combined.get_natoms();
	int nbas = combined.get_nbas();

	assert(shl_slice[1] <= nbas);
	assert(shl_slice[3] <= nbas);
	assert(shl_slice[5] <= nbas);

	ivec aoloc = Kernel::gen_loc(bas, nbas);
	int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
	int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];
	int naok = aoloc[shl_slice[5]] - aoloc[shl_slice[4]];

	libcint::CINTOpt* opty = nullptr;
	Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());

	// Compute integrals
	vec res((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);
	eri3c.resize((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);

	Kernel::drv(res.data(),
		1,
		shl_slice.data(),
		aoloc.data(),
		opty,
		atm.data(), nat,
		bas.data(), nbas,
		env.data());

	// FOR TESTING PURPOSES!!!!
	// GTOnr3c_drv(int3c2e_sph, res.data(), 1, shl_slice.data(), aoloc.data(), NULL, atm.data(), nat, bas.data(), nbas, env.data());

	// res is in fortran order, write the result in regular ordering
	for (int k = 0; k < naok; k++)
	{
		for (int j = 0; j < naoj; j++)
		{
			for (int i = 0; i < naoi; i++)
			{
				std::size_t idx_F = i + j * (size_t)naoi + k * ((size_t)naoi * (size_t)naoj);
				std::size_t idx_C = i * ((size_t)naoj * (size_t)naok) + j * (size_t)naok + k;
				eri3c[idx_C] = res[idx_F];
			}
		}
	}
}
template void computeEri3c<Coulomb3C_SPH>(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c);

template <typename Kernel>
void compute2C(Int_Params &params, vec &ret) {
	ivec bas = params.get_bas();
	ivec atm = params.get_atm();
	vec env = params.get_env();

	int nbas = params.get_nbas();
	int nat = params.get_natoms();

	ivec shl_slice = { 0, nbas, 0, nbas };
	ivec aoloc = Kernel::gen_loc(bas, nbas);

	int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
	int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];

	libcint::CINTOpt* opty = nullptr;
	Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());

	// Compute integrals
	vec res((size_t)naoi * (size_t)naoj, 0.0);
	ret.resize((size_t)naoi * (size_t)naoj, 0.0);
	Kernel::drv(res.data(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());

	// res is in fortran order, write the result in regular ordering
	for (int i = 0; i < naoi; i++)
	{
		for (int j = 0; j < naoj; j++)
		{
			ret[(size_t)j * (size_t)naoi + i] = res[(size_t)i * (size_t)naoj + j];
		}
	}
}
template void compute2C<Coulomb2C_SPH>(Int_Params &params, vec &ret);
template void compute2C<Coulomb2C_CRT>(Int_Params &params, vec &ret);
template void compute2C<Overlap2C_SPH>(Int_Params &params, vec &ret);
template void compute2C<Overlap2C_CRT>(Int_Params &params, vec &ret);


extern "C" libcint::CINTIntegralFunction int3c1e_sph;

//Per-shell integrals behind each 3C kernel. pair and aux give the Schwarz factors of
//|(ab|P)| <= sqrt((ab|ab)) sqrt((P|P)); nullptr where the metric has none.
//ponytail: the overlap metric is not screened, its bound needs int4c1e which this libcint build may lack
template <typename K> struct Shell3C;
template <> struct Shell3C<Coulomb3C_SPH> {
	static constexpr libcint::CINTIntegralFunction *three = libcint::int3c2e_sph, *pair = libcint::int2e_sph, *aux = libcint::int2c2e_sph;
};
template <> struct Shell3C<Coulomb3C_CRT> {
	static constexpr libcint::CINTIntegralFunction *three = libcint::int3c2e_cart, *pair = libcint::int2e_cart, *aux = libcint::int2c2e_cart;
};
template <> struct Shell3C<Overlap3C_SPH> {
	static constexpr libcint::CINTIntegralFunction *three = int3c1e_sph, *pair = nullptr, *aux = nullptr;
};

//rho_P = sum_ab w D_ab (ab|P) over orbital shell pairs a >= b (w = 2 off the diagonal), one shell
//triplet at a time. A triplet is skipped when w max|D_ab| Q_ab Q_P < NOS_RI_SCREEN, the
//density-weighted Schwarz screen XCW uses for its stored ERIs. This replaces an atom-pair overlap
//screen that rode on constants::exp_cutoff, which the ELI-D tail correction (16 Sep 2026) made
//so tight that it stopped screening.
//Tried 7 Oct 2026 on sucrose and dropped: QVl distance screening (Hollman, Schaefer, Valeev, JCP 142,
//154106 (2015)) lies on the same time/error curve as a tighter Schwarz threshold (1e-13: 3.14 s, 1.1e-7
//vs 3.24 s, 9e-8), and libcint's PTR_EXPCUTOFF at its floor of 40 instead of 60 changes nothing.
template <typename Kernel>
void computeRho(
	const Int_Params &normal_basis,
	const Int_Params &aux_basis,
	const dMatrix2 &dm,
	vec &rho,
	const std::optional<ivec> asym_atm_list)
{
	using S = Shell3C<Kernel>;
	Int_Params combined(normal_basis, aux_basis);

	ivec bas = combined.get_bas();
	ivec atm = combined.get_atm();
	vec  env = combined.get_env();

	const int nQM = normal_basis.get_nbas();
	const int nAux = aux_basis.get_nbas();
	const int nat = combined.get_natoms();
	const int nbas = combined.get_nbas();

	ivec aoloc = Kernel::gen_loc(bas, nbas);
	const int aux0 = aoloc[nQM];
	const int naux = aoloc[nQM + nAux] - aux0;
	rho.assign(naux, 0.0);

	double thr = constants::ri_screen_threshold;
	if (const char *e = tuning("NOS_RI_SCREEN")) thr = std::atof(e);
	const bool bound = S::pair != nullptr && thr > 0.0;

	int dmax_orb = 0, dmax_aux = 0;
	for (int s = 0; s < nbas; s++) {
		int &d = s < nQM ? dmax_orb : dmax_aux;
		d = std::max(d, aoloc[s + 1] - aoloc[s]);
	}
	//libcint mallocs scratch per call unless handed one; the size query is the call without output.
	//{s,s,s,s} over every shell is the bound pyscf's GTOmax_cache_size uses
	CACHE_SIZE_T ncache = 0;
	for (int s = 0; s < nbas; s++) {
		int shls[4] = { s, s, s, s };
		ncache = std::max(ncache, S::three(nullptr, nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), nullptr, nullptr));
		if (bound) ncache = std::max({ ncache,
			S::pair(nullptr, nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), nullptr, nullptr),
			S::aux(nullptr, nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), nullptr, nullptr) });
	}

	//Q_P = sqrt(max_p (p|p)) per aux shell
	vec qaux(nAux, 1.0);
	double qaux_max = 1.0;
	if (bound) {
		vec buf(static_cast<size_t>(dmax_aux) * dmax_aux), cache(ncache);
		qaux_max = 0.0;
		for (int P = 0; P < nAux; P++) {
			int shls[2] = { nQM + P, nQM + P };
			const int dp = aoloc[nQM + P + 1] - aoloc[nQM + P];
			double m = 0.0;
			if (S::aux(buf.data(), nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), nullptr, cache.data()))
				for (int f = 0; f < dp; f++) m = std::max(m, std::abs(buf[f + dp * f]));
			qaux[P] = std::sqrt(m);
			qaux_max = std::max(qaux_max, qaux[P]);
		}
	}

	libcint::CINTOpt *opty = nullptr, *opt2 = nullptr;
	Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());
	if (bound) libcint::int2e_optimizer(&opt2, atm.data(), nat, bas.data(), nbas, env.data());

	//Every shell a sums into its own vector, and those are added to rho in a fixed order as soon as all before
	//them are done: rho is bitwise the same whichever thread ran what, and only the out-of-order window is held.
	//Per-thread vectors added in arrival order moved the sucrose tsc by 1e-6 between identical runs (the metric
	//is ill-conditioned). Same speed as before (sucrose/TZVP 4.19 s vs 4.23 s); heaviest-a-first was slower.
	std::vector<vec> part(nQM);
	int merged = 0;

	long long done = 0;
	{ //the bar ends its line when it goes out of scope, before the summary below
		ProgressBar pb(nQM, 60, "#", " ", "Calculating Eri3c Matrix");
#pragma omp parallel reduction(+:done)
		{
			vec cache(ncache), dblk(static_cast<size_t>(dmax_orb) * dmax_orb);
			vec buf3(static_cast<size_t>(dmax_orb) * dmax_orb * dmax_aux), buf2(static_cast<size_t>(dmax_orb) * dmax_orb * dmax_orb * dmax_orb);
#pragma omp for schedule(dynamic)
			for (int a = 0; a < nQM; a++) {
				vec loc(naux, 0.0);
				const int a0 = aoloc[a], da = aoloc[a + 1] - a0;
				for (int b = 0; b <= a; b++) {
					const int b0 = aoloc[b], db = aoloc[b + 1] - b0, dab = da * db;
					const double w = a == b ? 1.0 : 2.0;
					//same element order as the integral block: i (shell a) fastest, then j (shell b)
					double dm_max = 0.0;
					for (int j = 0; j < db; j++)
						for (int i = 0; i < da; i++) {
							dblk[i + da * j] = w * dm(b0 + j, a0 + i);
							dm_max = std::max(dm_max, std::abs(dblk[i + da * j]));
						}
					if (dm_max == 0.0) continue;
					double qab = 1.0;
					if (bound) {
						int shls[4] = { a, b, a, b };
						double m = 0.0;
						if (S::pair(buf2.data(), nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), opt2, cache.data()))
							for (int j = 0; j < db; j++)
								for (int i = 0; i < da; i++) m = std::max(m, std::abs(buf2[i + da * j + dab * (i + da * j)]));
						qab = std::sqrt(m);
						if (dm_max * qab * qaux_max < thr) continue;
					}
					for (int P = 0; P < nAux; P++) {
						if (bound && dm_max * qab * qaux[P] < thr) continue;
						int shls[3] = { a, b, nQM + P };
						done++;
						if (!S::three(buf3.data(), nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), opty, cache.data())) continue;
						const int p0 = aoloc[nQM + P] - aux0, dp = aoloc[nQM + P + 1] - aoloc[nQM + P];
						for (int k = 0; k < dp; k++) {
							const double *col = buf3.data() + static_cast<size_t>(dab) * k;
							double sum = 0.0;
							for (int ij = 0; ij < dab; ij++) sum += col[ij] * dblk[ij];
							loc[p0 + k] += sum;
						}
					}
				}
				pb.update();
#pragma omp critical
				{
					part[a] = std::move(loc);
					for (; merged < nQM && !part[merged].empty(); merged++) {
						for (int p = 0; p < naux; p++) rho[p] += part[merged][p];
						vec().swap(part[merged]);
					}
				}
			}
		}
	}
	libcint::CINTdel_optimizer(&opty);
	if (opt2) libcint::CINTdel_optimizer(&opt2);
	const long long total = static_cast<long long>(nQM) * (nQM + 1) / 2 * nAux;
	std::cout << "Schwarz screening (" << thr << ") kept " << done << " of " << total << " shell triplets" << std::endl;
}
template void computeRho<Coulomb3C_SPH>(
	const Int_Params &normal_basis,
	const Int_Params &aux_basis,
	const dMatrix2 &dm,
	vec &rho,
	const std::optional<ivec> asym_atm_list);
template void computeRho<Coulomb3C_CRT>(
	const Int_Params &normal_basis,
	const Int_Params &aux_basis,
	const dMatrix2 &dm,
	vec &rho,
	const std::optional<ivec> asym_atm_list);
template void computeRho<Overlap3C_SPH>(
	const Int_Params &normal_basis,
	const Int_Params &aux_basis,
	const dMatrix2 &dm,
	vec &rho,
	const std::optional<ivec> asym_atm_list);


template <typename Kernel>
void compute3C(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c) {
	int nQM = param1.get_nbas();
	int nAux = param2.get_nbas();
	Int_Params combined(param1, param2);
	ivec bas = combined.get_bas();
	ivec atm = combined.get_atm();
	vec env = combined.get_env();
	int nat = combined.get_natoms();
	int nbas = combined.get_nbas();
	//ivec aoloc = make_loc(bas, nbas);
	ivec aoloc = Kernel::gen_loc(bas, nbas);
	unsigned long long int naoi = aoloc[nQM] - aoloc[0];
	unsigned long long int naoj = aoloc[nQM] - aoloc[0];
	unsigned long long int naok = aoloc[nQM + nAux] - aoloc[nQM];
	eri3c.resize(naoi * naoj * naok, 0.0);
	libcint::CINTOpt* opty = nullptr;
	Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());
	ivec shl_slice = { 0, nQM, 0, nQM, nQM, nQM + nAux };
	Kernel::drv(eri3c.data(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());

}
template void compute3C<Coulomb3C_SPH>(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c);
template void compute3C<Overlap3C_SPH>(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c);



//Cartesian to real spherical transformation matrix
//
//Kwargs :
//normalized:
//How the Cartesian GTOs are normalized.  'sp' means the s and p
//functions are normalized(this is the convention used by libcint
//    library).
dMatrix2 cart2sph(const int l, const bool normalized) {
	int n_cart = (l + 1) * (l + 2) / 2;

	dMatrix2 c_tensor(n_cart, n_cart);
	for (int i = 0; i < n_cart; i++) {
		c_tensor(i, i) = 1.0;
	}

	if (l == 0 || l == 1) { //For s and p functions, the transformation is trivial
		if (l == 1) { //libcint is built with PYPZPX: spherical p comes as py, pz, px
			std::fill(c_tensor.container().begin(), c_tensor.container().end(), 0.0);
			c_tensor(1, 0) = c_tensor(2, 1) = c_tensor(0, 2) = 1.0;
		}
		if (normalized) {
			return c_tensor;
		}
		else {
			double norm_factor = (l == 0) ? 0.282094791773878143 : 0.488602511902919921;
			for (auto &val : c_tensor.container()) {
				val *= norm_factor;
			}
			return c_tensor;
		}
	}

	err_checkf(l <= 15, "cart2sph_matrix: l must be <= 15", std::cout);

	int n_sph = 2 * l + 1;
	vec c_sph(static_cast<size_t>(n_sph) * n_cart, 0.0);

	libcint::CINTc2s_ket_sph(c_sph.data(), n_cart, c_tensor.data(), l);
	//Transform back to row-major order
	dMatrix2 c_sph_RM(n_cart, n_sph);
	for (int i = 0; i < n_cart; i++) {
		for (int j = 0; j < n_sph; j++) {
			c_sph_RM(i, j) = c_sph[j * n_cart + i];
		}
	}
	return c_sph_RM;
}


//Returns a n_cart*n_sph matrix used in transforming a cartesian DM to a spherical DM
dMatrix2 get_cart2sph_matrix(const WFN &cart_wfn, const bool normalized) {
	//First collect the complete number of spherical and cartesian functions used in the wavefunction
	int max_l = 0;
	int n_cart = 0, n_sph = 0;
	for (const atom &a : cart_wfn.get_atoms()) {
		int prim = 0;
		for (int shell = 0; shell < a.get_shellcount_size(); shell++) {
			const int type = a.get_basis_set_type(prim) - 1;

			n_cart += ((type + 1) * (type + 2)) / 2;
			n_sph += 2 * type + 1;
			if (type > max_l) {
				max_l = type;
			}

			prim += a.get_shellcount(shell);
		}
	}

	std::vector<dMatrix2> conversion_matrices(max_l + 1);
	for (int l = 0; l <= max_l; l++) {
		conversion_matrices[l] = cart2sph(l, normalized);
	}

	std::cout << "Number of cartesian functions: " << n_cart << ", number of spherical functions: " << n_sph << std::endl;
	dMatrix2 c_dm(n_cart, n_sph);
	int cart_idx = 0, sph_idx = 0;
	for (const atom &a : cart_wfn.get_atoms()) {
		int prim = 0;
		for (int shell = 0; shell < a.get_shellcount_size(); shell++) {
			const int type = a.get_basis_set_type(prim) - 1;
			const dMatrix2 &c_sph = conversion_matrices[type];
			const int n_cart_shell = ((type + 1) * (type + 2)) / 2;
			const int n_sph_shell = 2 * type + 1;
			//Fill in the appropriate block in the c_dm matrix
			for (int i = 0; i < n_cart_shell; i++) {
				for (int j = 0; j < n_sph_shell; j++) {
					c_dm(cart_idx + i, sph_idx + j) = c_sph(i, j);
				}
			}
			cart_idx += n_cart_shell;
			sph_idx += n_sph_shell;
			prim += a.get_shellcount(shell);
		}
	}
	////print all of c_dm
	//for (int i = 0; i < c_dm.extent(0); i++) {
	//    for (int j = 0; j < c_dm.extent(1); j++) {
	//        std::cout << c_dm(i, j) << " ";
	//    }
	//    std::cout << std::endl;
	//}

	return c_dm;
}


vec eval_GTO_sph(Int_Params& params, vec2& grid, ivec& shl_slice) {
	ivec bas = params.get_bas();
	ivec atm = params.get_atm();
	vec env = params.get_env();


	//grid = numpy.asarray(grid, dtype = numpy.double, order = 'F'); grid holds the three coordinate rows
	const int ngrid = (int)grid[0].size();
	vec fortran_grid((size_t)ngrid * 3);
	for (int i = 0; i < ngrid; i++) {
		for (int j = 0; j < 3; j++) {
			fortran_grid[(size_t)i + (size_t)j * ngrid] = grid[j][i];
		}
	}

	int nbas = params.get_nbas();
	int nat = params.get_natoms();

	if (shl_slice.size() == 0) {
		shl_slice = { 0, nbas};
	}
	//ivec aoloc = Kernel::gen_loc(bas, nbas);
	ivec aoloc = make_loc<COORDINATE_TYPE::SPH>(bas, nbas);
	int nao = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];

	//non0tab = numpy.ones(((ngrids+BLKSIZE-1)//BLKSIZE,nbas),dtype = numpy.uint8)
	std::vector<uint8_t> non0table(static_cast<size_t>((ngrid + 56 - 1) / 56) * nbas, 1);

	// Compute integrals
	vec res((size_t)nao * (size_t)ngrid, 0.0);
	GTOval_sph(ngrid, shl_slice.data(), aoloc.data(), res.data(), fortran_grid.data(), non0table.data(), atm.data(), nat, bas.data(), nbas, env.data());

	vec ret((size_t)nao * (size_t)ngrid, 0.0);
	// res is in fortran order, write the result in regular ordering
	for (int i = 0; i < nao; i++)
	{
		for (int j = 0; j < ngrid; j++)
		{
			ret[(size_t)j * (size_t)nao + i] = res[(size_t)i * (size_t)ngrid + j];
		}
	}
	return ret;
}