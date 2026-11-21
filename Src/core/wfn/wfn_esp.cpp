#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "mo_class.h"
#include "cube.h"
#include "constants.h"
#include "fchk.h"
#include "basis_set.h"
#include "nos_math.h"
#include "libCintMain.h"
#include "eli_family.h"
#include "integrator.h"
#include "cell.h"
#ifdef NOSPHERA2_USE_GPU
#include "aux_density_gpu.h"
#include "esp_gpu.h"
#include "basin_field_gpu.h"
#endif

// Boys function F_m(T) by a 6-term Taylor expansion around a tabulated grid (step 0.1 up to
// T = 30, error < 1E-10), asymptotic beyond (1E-14); expn = exp(-T), which the caller has anyway
// Pairs stop at g x g, so MaxFn <= 8 and the Taylor reads m .. m + 5: rows of 15, not 31, keep a row in two cache lines
static constexpr int boys_mmax = 8 + 6, boys_nT = 301;
static constexpr double boys_step = 0.1;
// 1/k! factors as multiplies: a double division is a software sequence on a device
static constexpr double boys_inv_k[6] = { 0.0, 1.0, 0.5, 1.0 / 3.0, 0.25, 0.2 };
// 1 / (2n - 1) for the downward F_n recursion, n = 1 .. MaxFn
static constexpr double boys_inv_odd[9] = { 0.0, 1.0, 1.0 / 3.0, 0.2, 1.0 / 7.0, 1.0 / 9.0, 1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0 };
// File scope: no static-init guard in the hot loop, and esp_boys_table() hands it to the GPU
static const vec boys_tab = []()
{
	vec F((size_t)boys_nT * (boys_mmax + 1));
	for (int i = 0; i < boys_nT; i++)
	{
		const double T0 = i * boys_step, e = exp(-T0);
		// F_mmax by its series e^-T sum_k (2T)^k / ((2m+1)...(2m+2k+1)), then downward
		double term = 1.0 / (2 * boys_mmax + 1), sum = term;
		for (int k = 1; term > 1E-17 * sum; k++)
			term *= 2 * T0 / (2 * boys_mmax + 2 * k + 1), sum += term;
		double *row = F.data() + (size_t)i * (boys_mmax + 1);
		row[boys_mmax] = e * sum;
		for (int n = boys_mmax - 1; n >= 0; n--)
			row[n] = (2 * T0 * row[n + 1] + e) / (2 * n + 1);
	}
	return F;
}();

const double *esp_boys_table(int &nT, int &stride, double &step)
{
	nT = boys_nT, stride = boys_mmax + 1, step = boys_step;
	return boys_tab.data();
}

static double boys(const int m, const double T, const double expn)
{
	if (T >= (boys_nT - 1) * boys_step)
	{
		double f = 0.5 * sqrt(constants::PI / T); // F_0, then upward, stable at this T
		for (int n = 1; n <= m; n++)
			f = ((2 * n - 1) * f - expn) / (2 * T);
		return f;
	}
	const int i = (int)(T * (1.0 / boys_step) + 0.5);
	const double dT = i * boys_step - T, *row = boys_tab.data() + (size_t)i * (boys_mmax + 1) + m;
	double f = row[0], pw = 1;
	for (int k = 1; k <= 5; k++)
		pw *= dT * boys_inv_k[k], f += row[k] * pw;
	return f;
}
const double WFN::fj(int &j, int &l, int &m, double &aa, double &bb) const
{
	double temp = 0.0;
	double temp2 = 0.0;
	int a = 0, b = 0;
	for (int i = std::max(0, j - m); i <= std::min(j, l); i++)
	{
		// pre = factorial[l] / factorial[l - i] / factorial[i] * factorial[m] / factorial[m - j + i] / factorial[j - i];
		temp2 = static_cast<double>(pre[j][l][m][i]);
		a = l - i;
		b = m + i - j;
		if (a != 0)
			temp2 *= pow(aa, a);
		if (b != 0)
			temp2 *= pow(bb, b);
		temp += temp2;
	}
	return temp;
};

// (l, r, s) terms of one axis with l_i + l_j = L: sum over l <= L, r <= l / 2 of (l - 2 r) / 2 + 1; L <= 8 (g x g)
static int esp_axis_terms(const int L)
{
	static constexpr int n[9] = { 1, 2, 5, 8, 14, 20, 30, 40, 55 };
	return n[L];
}

WFN::ESP_pairs WFN::build_ESP_pairs() const
{
	ESP_pairs t;
	const int MO = get_nmo(true), stride = get_nmo(false);
	const int nprim = get_nex();
	const double *coef = get_coef_primitive_major(); // [prim][mo]
	int l_i[3], l_j[3];
	double Pi[3], Pj[3];
	// Tables are sized for g x g (Afac_pre [9][5][9], Fn [9]): an h primitive would index past them
	// silently, and the OCC reader accepts l up to 10
	for (int iprim = 0; iprim < nprim; iprim++)
	{
		constants::type2vector(get_type(iprim), l_i);
		if (l_i[0] < 0 || l_i[0] + l_i[1] + l_i[2] > 4)
			err_not_impl_f("ESP of a primitive beyond g functions (type " + std::to_string(get_type(iprim)) + ")", std::cout);
	}
	t.off.push_back(0);
	for (int iprim = 0; iprim < nprim; iprim++)
	{
		const int iat = get_center(iprim) - 1;
		constants::type2vector(get_type(iprim), l_i);
		const double iex = get_exponent(iprim);
		for (int jprim = iprim; jprim < nprim; jprim++)
		{
			const int jat = get_center(jprim) - 1;
			const double jex = get_exponent(jprim);
			const double ex_sum = iex + jex;
			double sqd = 0;
			for (int k = 0; k < 3; k++)
				sqd += pow(atoms[iat].get_coordinate(k) - atoms[jat].get_coordinate(k), 2);
			const double prefac = constants::TWO_PI / ex_sum * exp(-iex * jex * sqd / ex_sum);
			if (prefac < 1E-10)
				continue;
			double dij = 0;
			for (int mo = 0; mo < MO; mo++)
				dij += get_MO_occ(mo) * coef[iprim * stride + mo] * coef[jprim * stride + mo];
			// screening: the dropped pairs move the far-field ESP of sucrose by 3E-6 at most (18 Sep 2026)
			if (abs(prefac * dij) < 1E-8)
				continue;
			constants::type2vector(get_type(jprim), l_j);
			t.ex_sum.push_back(ex_sum);
			t.weight.push_back(prefac * dij * (iprim != jprim ? 2.0 : 1.0));
			d3 P;
			for (int k = 0; k < 3; k++)
			{
				P[k] = (atoms[iat].get_coordinate(k) * iex + atoms[jat].get_coordinate(k) * jex) / ex_sum;
				Pi[k] = P[k] - atoms[iat].get_coordinate(k);
				Pj[k] = P[k] - atoms[jat].get_coordinate(k);
			}
			t.P.push_back(P);
			t.L.push_back({ l_i[0] + l_j[0], l_i[1] + l_j[1], l_i[2] + l_j[2] });
			// Afac without its PC power: sign * fj * (1/4g)^(r+s) / (r! s! (l-2r-2s)!)
			for (int k = 0; k < 3; k++)
				for (int l = 0; l <= l_i[k] + l_j[k]; l++)
				{
					const double fjtmp = (l % 2 ? -1.0 : 1.0) * fj(l, l_i[k], l_j[k], Pi[k], Pj[k]);
					for (int r = 0; r <= l / 2; r++)
						for (int s = 0; s <= (l - 2 * r) / 2; s++)
						{
							t.coef.push_back((s % 2 ? -fjtmp : fjtmp) * pow(0.25 / ex_sum, r + s) / Afac_pre[l][r][s]);
							t.pc_pow.push_back((unsigned char)(l - 2 * r - 2 * s));
							t.fn_idx.push_back((unsigned char)(l - 2 * r - s));
						}
				}
			t.off.push_back((int)t.coef.size());
		}
	}
	return t;
}

// Pair loop over VL points: every index belongs to the pair, so lanes carry only arithmetic and vectorise.
// Each lane replays the scalar operations in scalar order, so any VL is bitwise identical (EspTests pins it)
template <int VL>
static void esp_lanes(const WFN &wfn, const WFN::ESP_pairs &t,
					  const double *px, const double *py, const double *pz, double *out)
{
	double ESP[VL];
	for (int v = 0; v < VL; v++)
		ESP[v] = 0.0;
	for (int iat = 0; iat < wfn.get_ncen(); iat++)
	{
		const double cx = wfn.get_atom_coordinate(iat, 0), cy = wfn.get_atom_coordinate(iat, 1), cz = wfn.get_atom_coordinate(iat, 2);
		// ECP/xTB/pTB: MOs hold only valence electrons, so the core is removed from the nuclear charge
		const double Z = wfn.get_atom_charge(iat) - wfn.get_atom_ECP_electrons(iat);
		for (int v = 0; v < VL; v++)
		{
			const double dx = px[v] - cx, dy = py[v] - cy, dz = pz[v] - cz;
			ESP[v] += Z / sqrt((dx * dx + dy * dy) + dz * dz);
		}
	}

	// B[axis][k]: one axis summed by the F index its (l, r, s) terms feed, so the product loop indexes Fn directly
	double Fn[9][VL], B[3][9][VL], pw[9][VL], PC[3][VL], T[VL], expc[VL], term[VL]; // MaxFn <= 8, L[k] <= 8 for g x g
	const int npairs = (int)t.weight.size();
	for (int p = 0; p < npairs; p++)
	{
		const double ex_sum = t.ex_sum[p];
		const std::array<int, 3> &L = t.L[p];
		const int MaxFn = L[0] + L[1] + L[2];
		const double Px = t.P[p][0], Py = t.P[p][1], Pz = t.P[p][2], w = t.weight[p];
		for (int v = 0; v < VL; v++)
		{
			PC[0][v] = Px - px[v], PC[1][v] = Py - py[v], PC[2][v] = Pz - pz[v];
			T[v] = ex_sum * ((PC[0][v] * PC[0][v] + PC[1][v] * PC[1][v]) + PC[2][v] * PC[2][v]);
		}
		int c = t.off[p];
		// s-s: one term per axis, PC power 0, and F_0 never reads exp(-T): bitwise the general path without the exp
		if (MaxFn == 0)
		{
			const double cc = (t.coef[c] * t.coef[c + 1]) * t.coef[c + 2];
			for (int v = 0; v < VL; v++)
				ESP[v] -= w * (cc * boys(0, T[v], 0.0));
			continue;
		}
		// Beyond T = 60, exp(-T) is below the last bit of F_0, so the libm call is skipped
		for (int v = 0; v < VL; v++)
		{
			expc[v] = T[v] < 60.0 ? exp(-T[v]) : 0.0;
			Fn[MaxFn][v] = boys(MaxFn, T[v], expc[v]);
		}
		for (int nu = MaxFn - 1; nu >= 0; nu--)
			for (int v = 0; v < VL; v++)
				Fn[nu][v] = (expc[v] + 2 * T[v] * Fn[nu + 1][v]) * boys_inv_odd[nu + 1];

		// Each (l, r, s) term feeds one F index, so summing each axis into B[k] first factorises the
		// nl * nm * nn triple product exactly into (L0+1)(L1+1)(L2+1) terms with index kx + ky + kz
		for (int k = 0; k < 3; k++)
		{
			for (int v = 0; v < VL; v++)
				pw[0][v] = 1.0;
			for (int n = 1; n <= L[k]; n++)
				for (int v = 0; v < VL; v++)
					pw[n][v] = pw[n - 1][v] * PC[k][v];
			for (int i = 0; i <= L[k]; i++)
				for (int v = 0; v < VL; v++)
					B[k][i][v] = 0.0;
			for (int i = esp_axis_terms(L[k]); i--; c++)
			{
				const double cf = t.coef[c];
				const double *const p_pw = pw[t.pc_pow[c]];
				double *const p_B = B[k][t.fn_idx[c]];
				for (int v = 0; v < VL; v++)
					p_B[v] += cf * p_pw[v];
			}
		}

		for (int v = 0; v < VL; v++)
			term[v] = 0.0;
		for (int kx = 0; kx <= L[0]; kx++)
			for (int ky = 0; ky <= L[1]; ky++)
				for (int kz = 0; kz <= L[2]; kz++)
				{
					const double *const f = Fn[kx + ky + kz];
					for (int v = 0; v < VL; v++)
						term[v] += (B[0][kx][v] * B[1][ky][v]) * B[2][kz][v] * f[v];
				}
		for (int v = 0; v < VL; v++)
			ESP[v] -= w * term[v];
	}
	for (int v = 0; v < VL; v++)
		out[v] = ESP[v];
}

const double WFN::computeESP(const d3 &PosGrid, const ESP_pairs &t) const
{
	double ESP = 0;
	esp_lanes<1>(*this, t, &PosGrid[0], &PosGrid[1], &PosGrid[2], &ESP);
	return ESP;
};

// Host path: VL points share one pass over the pair table, leftovers go scalar.  The loop is bound by
// fetching the table, so wider passes win until the lane state no longer fits in L1
template <int VL>
static void esp_host_pass(const WFN &wfn, const WFN::ESP_pairs &t, const std::vector<d3> &points, double *out)
{
	const int np = (int)points.size(), nblk = np / VL;
#pragma omp parallel for schedule(static)
	for (int b = 0; b < nblk; b++)
	{
		double px[VL], py[VL], pz[VL];
		for (int v = 0; v < VL; v++)
			px[v] = points[b * VL + v][0], py[v] = points[b * VL + v][1], pz[v] = points[b * VL + v][2];
		esp_lanes<VL>(wfn, t, px, py, pz, out + (size_t)b * VL);
	}
	for (int i = nblk * VL; i < np; i++)
		out[i] = wfn.computeESP(points[i], t);
}

// Every point walks the whole O(nprim^2) pair table, so a whole set goes to the device at once; host
// fallback without a device or when the set is too small to pay the copies
void WFN::computeESP_batch(const std::vector<d3> &points, const ESP_pairs &t, double *out) const
{
	const int np = (int)points.size();
	if (np == 0)
		return;
#ifdef NOSPHERA2_USE_GPU
	const int ncen = get_ncen(), npairs = (int)t.weight.size();
	// The device buffers stay allocated between calls; each call uploads the atoms and the pair table, O(npairs) next to O(np * npairs)
	if (npairs > 0 && (long long)np * npairs >= (1LL << 22) && aux_density_gpu_enabled())
	{
		vec ax(ncen), ay(ncen), az(ncen), q(ncen);
		for (int a = 0; a < ncen; a++)
		{
			ax[a] = atoms[a].get_coordinate(0), ay[a] = atoms[a].get_coordinate(1), az[a] = atoms[a].get_coordinate(2);
			q[a] = get_atom_charge(a) - atoms[a].get_ECP_electrons();
		}
		int nT = 0, stride = 0;
		double step = 0;
		const double *tab = esp_boys_table(nT, stride, step);
		// One launch for the whole set: the kernel is latency bound and needs every SM occupied, so
		// split launches lose more in tail occupancy than overlapping the host cores gains
		if (esp_gpu_eval(ncen, ax.data(), ay.data(), az.data(), q.data(),
						 npairs, t.ex_sum.data(), t.weight.data(), t.P[0].data(), t.L[0].data(),
						 t.off.data(), t.coef.data(), t.pc_pow.data(), t.fn_idx.data(),
						 nT, stride, step, tab, np, points[0].data(), out))
			return;
	}
#endif
	// Never fewer VL blocks than threads: idle cores cost more than the wide pass gains.  All widths are bitwise identical
	const int nthr = omp_get_max_threads();
	if (np >= 64 * nthr)
		esp_host_pass<64>(*this, t, points, out);
	else if (np >= 8 * nthr)
		esp_host_pass<8>(*this, t, points, out);
	else
		esp_host_pass<1>(*this, t, points, out);
}
