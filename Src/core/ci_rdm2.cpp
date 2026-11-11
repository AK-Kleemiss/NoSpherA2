#include "pch.h"
#include "ci_rdm2.h"
#include <bit>
#include <unordered_map>
#include <numeric>

namespace
{
	using Det = std::array<uint64_t, 2>; //alpha, beta occupation bits
	struct DetHash
	{
		size_t operator()(const Det& d) const { return std::hash<uint64_t>()(d[0] * 0x9E3779B97F4A7C15ull ^ d[1]); }
	};

	//a+ or a of orbital t, spin s, in the alpha-first order |a1 a2 ... b1 b2 ...>; 0 when it kills the state
	int apply(Det& d, const int s, const int t, const bool create)
	{
		const uint64_t bit = 1ull << t;
		if (bool(d[s] & bit) == create) return 0;
		const int before = std::popcount(d[s] & (bit - 1)) + (s == 1 ? std::popcount(d[0]) : 0);
		d[s] ^= bit;
		return before & 1 ? -1 : 1;
	}

	//T_{..t..} -> sum_t U(t,j) T_{..t..} on every index in turn; U(t,j) at [t*n+j]
	void rotate(vec& T, const int rank, const int n, const vec& U)
	{
		size_t pre = 1, post = T.size();
		for (int k = 0; k < rank; k++)
		{
			post /= n;
			vec out(T.size(), 0.0);
			for (size_t a = 0; a < pre; a++)
				for (int t = 0; t < n; t++)
					for (int j = 0; j < n; j++)
					{
						const double u = U[size_t(t) * n + j];
						if (u == 0.0) continue;
						const double* src = &T[(a * n + t) * post];
						double* dst = &out[(a * n + j) * post];
						for (size_t b = 0; b < post; b++) dst[b] += u * src[b];
					}
			T.swap(out);
			pre *= n;
		}
	}
}

void ci_rdm2_from_dets(const std::vector<std::pair<std::string, double>>& dets, CiRdm2& r)
{
	const int n = r.n = static_cast<int>(dets.front().first.size());
	std::unordered_map<Det, double, DetHash> c;
	for (const auto& [str, coef] : dets)
	{
		//ORCA orders spin orbitals 1a 1b 2a 2b ...: to alpha-first each alpha passes the betas below it
		Det d{ 0, 0 };
		int swaps = 0;
		for (int t = 0; t < n; t++)
		{
			if (str[t] == '2' || str[t] == 'u') { d[0] |= 1ull << t; swaps += std::popcount(d[1]); }
			if (str[t] == '2' || str[t] == 'd') d[1] |= 1ull << t;
		}
		c[d] = swaps & 1 ? -coef : coef;
	}
	const size_t n2 = size_t(n) * n;
	r.g1[0].assign(n2, 0.0);
	r.g1[1].assign(n2, 0.0);
	r.G.assign(n2 * n2, 0.0);
	r.Gab.assign(n2 * n2, 0.0);
	//<I| a+_P a+_R a_S a_Q |J> over the printed J; an I that was not printed has c_I = 0.
	//ponytail: serial, N_det n_el^2 n^2 hash lookups; a CAS(12,12) of 10^6 determinants wants OpenMP over J
	for (const auto& [J, cJ] : c)
		for (int sq = 0; sq < 2; sq++)
			for (int q = 0; q < n; q++)
			{
				Det J1 = J;
				const int f1 = apply(J1, sq, q, false);
				if (!f1) continue;
				for (int p = 0; p < n; p++)
				{
					Det I = J1;
					const int f2 = apply(I, sq, p, true);
					const auto it = f2 ? c.find(I) : c.end();
					if (it != c.end()) r.g1[sq][size_t(p) * n + q] += it->second * cJ * f1 * f2;
				}
				for (int ss = 0; ss < 2; ss++)
					for (int s = 0; s < n; s++)
					{
						Det J2 = J1;
						const int f2 = apply(J2, ss, s, false);
						if (!f2) continue;
						for (int rr = 0; rr < n; rr++)
						{
							Det J3 = J2;
							const int f3 = apply(J3, ss, rr, true);
							if (!f3) continue;
							for (int p = 0; p < n; p++)
							{
								Det I = J3;
								const int f4 = apply(I, sq, p, true);
								const auto it = f4 ? c.find(I) : c.end();
								if (it == c.end()) continue;
								const double v = it->second * cJ * f1 * f2 * f3 * f4;
								r.G[r.at(p, q, rr, s)] += v;
								if (sq == 0 && ss == 1) r.Gab[r.at(p, q, rr, s)] += v;
							}
						}
					}
			}
}

double ci_on_top_ab(const WFN& wave, const CiRdm2& r, const d3& p, double rho[2])
{
	//One spin's share of the internal density: alpha-beta pairs are core-core, core-active and active-active
	double core = 0.0;
	for (int i = 0; i < r.first; i++)
	{
		const double phi = wave.computeMO(p, i);
		core += phi * phi;
	}
	const int n = r.n;
	const size_t n2 = size_t(n) * n;
	vec P(n2);
	{
		vec phi(n);
		for (int t = 0; t < n; t++) phi[t] = wave.computeMO(p, r.first + t);
		for (int t = 0; t < n; t++)
			for (int u = 0; u < n; u++) P[size_t(t) * n + u] = phi[t] * phi[u];
	}
	const double act[2] = { std::inner_product(P.begin(), P.end(), r.g1[0].begin(), 0.0),
		std::inner_product(P.begin(), P.end(), r.g1[1].begin(), 0.0) };
	double ot = 0.0;
	for (size_t i = 0; i < n2; i++)
		if (P[i] != 0.0) ot += P[i] * std::inner_product(P.begin(), P.end(), r.Gab.begin() + i * n2, 0.0);
	rho[0] = core + act[0];
	rho[1] = core + act[1];
	return core * core + core * act[1] + act[0] * core + ot;
}

bool load_ci_rdm2(const WFN& wave, CiRdm2& r, std::ostream& log)
{
	std::filesystem::path out;
	for (const char* ext : { ".out", ".log" })
	{
		std::filesystem::path p = wave.get_path();
		p.replace_extension(ext);
		if (std::filesystem::exists(p)) { out = p; break; }
	}
	const std::string how = "(ORCA: %casscf PrintWF det, TPrintWF 1e-14 end, output kept as <stem>.out next to the gbw)";
	if (out.empty())
	{
		log << "  no 2-RDM: no " << wave.get_path().stem().string() << ".out/.log next to the wavefunction " << how << "\n";
		return false;
	}
	std::ifstream in(out);
	std::vector<std::pair<std::string, double>> dets;
	int internal = 0, active = -1, roots = 0;
	bool block = false, seen = false;
	std::string line;
	while (std::getline(in, line))
	{
		if (line.find("Spin-Determinant CI Printing") != std::string::npos)
		{
			block = seen = true;
			dets.clear();
			roots = 0;
			continue;
		}
		std::istringstream s(line);
		std::string w;
		s >> w;
		if (!block)
		{
			//"Internal 0 - -1 ( 0 orbitals)", "Active 0 - 18 ( 19 orbitals)"; the last one counts
			int a, b;
			std::string dash;
			if ((w == "Internal" || w == "Active") && (s >> a >> dash >> b) && dash == "-") (w == "Internal" ? internal : active) = b - a + 1;
			continue;
		}
		if (w == "ROOT") { roots++; continue; }
		if (roots == 0) continue;
		if (w.size() > 2 && w.front() == '[' && w.back() == ']')
		{
			double coef;
			if (roots == 1 && (s >> coef)) dets.emplace_back(w.substr(1, w.size() - 2), coef);
			continue;
		}
		if (!w.empty()) block = false;
	}
	auto refuse = [&](const std::string& why) {
		log << "  no 2-RDM from " << out.filename().string() << ": " << why << "\n";
		return false;
	};
	if (!seen || dets.empty()) return refuse("no determinant CI print " + how);
	if (roots != 1) return refuse(std::to_string(roots) + " roots printed; a state average is no one state's 2-RDM");
	const int n = static_cast<int>(dets.front().first.size());
	if (n > 64) return refuse(std::to_string(n) + " active orbitals, at most 64");
	if (n != active) return refuse("determinants over " + std::to_string(n) + " orbitals, active space " + std::to_string(active));
	double norm = 0.0;
	for (const auto& d : dets)
	{
		if (d.first.size() != size_t(n)) return refuse("determinant strings of different lengths");
		norm += d.second * d.second;
	}
	if (std::abs(norm - 1.0) > 1e-6) return refuse("the printed determinants carry a norm of " + std::to_string(norm) + ", TPrintWF was too loose");
	for (int mo = 0; mo < wave.get_nmo(); mo++)
		if (wave.get_MO_op(mo) == 1) return refuse("an unrestricted wavefunction; the CI belongs to one restricted orbital set");
	if (internal + n > wave.get_nmo()) return refuse("internal + active orbitals exceed the " + std::to_string(wave.get_nmo()) + " MOs read");
	r.first = internal;
	ci_rdm2_from_dets(dets, r);
	r.source = out;
	//The gbw of an ORCA full CI keeps its guess orbitals; read_gbw turned them into natural orbitals C U, and
	//the CI, written in the gbw orbitals, has to follow.  U spans all MOs, so only an all-active CI maps
	const vec& U = wave.get_natural_rotation();
	if (!U.empty())
	{
		if (internal != 0 || n != wave.get_nmo())
			return refuse("the gbw orbitals were rotated to natural orbitals, which maps only a CI over every orbital");
		rotate(r.g1[0], 2, n, U);
		rotate(r.g1[1], 2, n, U);
		rotate(r.G, 4, n, U);
		rotate(r.Gab, 4, n, U);
	}
	double dev = 0.0;
	for (int mo = 0; mo < wave.get_nmo(); mo++)
	{
		const int t = mo - internal;
		if (t < 0) dev = std::max(dev, std::abs(wave.get_MO_occ(mo) - 2.0));
		else if (t >= n) dev = std::max(dev, std::abs(wave.get_MO_occ(mo)));
		else
			for (int u = 0; u < n; u++)
				dev = std::max(dev, std::abs(r.g1[0][size_t(t) * n + u] + r.g1[1][size_t(t) * n + u] - (t == u ? wave.get_MO_occ(mo) : 0.0)));
	}
	if (dev > 1e-4)
		return refuse("its 1-RDM differs from the orbital occupations by " + std::to_string(dev) + ": the CI is not in these orbitals");
	std::ostringstream s;
	s << "  2-RDM: " << dets.size() << " determinants from " << out.filename().string() << ", " << internal << " internal + " << n
		<< " active orbitals, |1 - norm| " << std::scientific << std::setprecision(1) << std::abs(1.0 - norm)
		<< ", 1-RDM vs occupations " << dev << "\n";
	log << s.str();
	return true;
}
