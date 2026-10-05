#include "pch.h"
#include "nbo.h"
#include "wfn_class.h"
#include "constants.h"
#include "convenience.h"

#include <Eigen/Dense>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace
{
	//copies of the helpers in nao.cpp's anonymous namespace
	MatrixXd to_eigen(const dMatrix2& m)
	{
		const int r = static_cast<int>(m.extent(0)), c = static_cast<int>(m.extent(1));
		MatrixXd out(r, c);
		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++)
				out(i, j) = m(i, j);
		return out;
	}

	dMatrix2 to_dmatrix(const MatrixXd& m)
	{
		dMatrix2 out(m.rows(), m.cols());
		for (int i = 0; i < m.rows(); i++)
			for (int j = 0; j < m.cols(); j++)
				out(i, j) = m(i, j);
		return out;
	}

	MatrixXd sym_power(const MatrixXd& M, const double p, const double rel_floor = 1e-10)
	{
		//an empty spin channel
		if (M.rows() == 0 || M.cols() == 0)
			return M;
		Eigen::SelfAdjointEigenSolver<MatrixXd> es(M);
		VectorXd w = es.eigenvalues();
		const double cut = rel_floor * std::max(w.maxCoeff(), 1e-300);
		VectorXd f(w.size());
		for (int i = 0; i < w.size(); i++)
			f(i) = (w(i) > cut) ? std::pow(w(i), p) : (p < 0.0 ? 0.0 : std::pow(std::max(w(i), 0.0), p));
		return es.eigenvectors() * f.asDiagonal() * es.eigenvectors().transpose();
	}

	//body of a $SECTION between the keyword and $END
	std::string section_body(const std::string& text, const std::string& name)
	{
		const size_t start = text.find(name);
		if (start == std::string::npos) return std::string();
		const size_t after = start + name.size();
		const size_t end = text.find("$END", after);
		return text.substr(after, (end == std::string::npos ? text.size() : end) - after);
	}

	//every number of a section body, skipping the "CENTER =" tags
	vec numbers_of(const std::string& body)
	{
		vec out;
		size_t i = 0;
		while (i < body.size()) {
			const char c = body[i];
			const bool numeric = std::isdigit(static_cast<unsigned char>(c)) ||
				((c == '-' || c == '+' || c == '.') && i + 1 < body.size() &&
					(std::isdigit(static_cast<unsigned char>(body[i + 1])) || body[i + 1] == '.'));
			if (!numeric) { i++; continue; }
			size_t used = 0;
			try {
				out.push_back(std::stod(body.substr(i, 32), &used));
			}
			catch (const std::exception&) {
				i++;
				continue;
			}
			i += std::max<size_t>(used, 1);
		}
		return out;
	}

	//FILE47 packs the upper triangle column-major, i.e. the lower triangle row-major
	dMatrix2 unpack(const vec& values, const size_t offset, const int n)
	{
		dMatrix2 M(n, n);
		size_t k = offset;
		for (int j = 0; j < n; j++)
			for (int i = 0; i <= j; i++, k++) {
				M(i, j) = values[k];
				M(j, i) = values[k];
			}
		return M;
	}

	//NBO LABEL codes in write_nbo() order: a shell's components run m = 0, +1, -1, +2, -2, ...
	bool decode_label(const int label, int& l, int& component)
	{
		static const ivec2 codes = {
			{ 1 },
			{ 103, 101, 102 },
			{ 255, 252, 253, 254, 251 },
			{ 351, 352, 353, 354, 355, 356, 357 },
			{ 451, 452, 453, 454, 455, 456, 457, 458, 459 } };
		for (int ll = 0; ll < static_cast<int>(codes.size()); ll++)
			for (int c = 0; c < static_cast<int>(codes[ll].size()); c++)
				if (codes[ll][c] == label) { l = ll; component = c; return true; }
		return false;
	}

	//NBO's angular label for component c of shell l, in FILE47's m = 0, +1, -1, ... order
	std::string lang_label(const int l, const int c)
	{
		static const char* p[3] = { "pz", "px", "py" };
		static const char* d[5] = { "dz2", "dxz", "dyz", "dx2y2", "dxy" };
		if (l == 0) return "s";
		if (l == 1) return p[c];
		if (l == 2) return d[c];
		const int m = (c == 0) ? 0 : ((c + 1) / 2) * ((c % 2) ? 1 : -1);
		const char letter = "spdfghik"[std::min(l, 7)];
		if (m == 0) return std::string(1, letter) + "(0)";
		return std::string(1, letter) + "(" + (m > 0 ? "c" : "s") + std::to_string(std::abs(m)) + ")";
	}

	//NBO's NAO table order: s; px, py, pz; dxy, dxz, dyz, dx2y2, dz2; f(0), f(c1), f(s1), ...
	//so for p and d not the FILE47 order
	int lang_rank(const int l, const int c)
	{
		static const int p[3] = { 2, 0, 1 };        //pz, px, py -> 3rd, 1st, 2nd
		static const int d[5] = { 4, 1, 2, 3, 0 };  //dz2, dxz, dyz, dx2y2, dxy
		if (l == 1) return p[c];
		if (l == 2) return d[c];
		return c;
	}

	std::string pad(const std::string& s, const size_t w)
	{
		return s.size() >= w ? s : std::string(w - s.size(), ' ') + s;
	}

	std::string pad(const int v, const size_t w) { return pad(std::to_string(v), w); }

	//compare_nbo_results() matches orbitals on nbo_run.cpp's whitespace-collapsed description,
	//so the native side must produce exactly that string
	std::string collapse(const std::string& s)
	{
		std::string out;
		bool space = false;
		for (const char c : s) {
			if (std::isspace(static_cast<unsigned char>(c))) { space = !out.empty(); continue; }
			if (space) out.push_back(' ');
			space = false;
			out.push_back(c);
		}
		return out;
	}

	//"BD (1) O 1- H 2" / "LP (2) O 1 s" - the layout of the NBO hybrid table.
	std::string orbital_description(const NboFunction& f, const std::vector<NAOAtom>& atoms)
	{
		std::string s = f.type + " (" + std::to_string(f.multiplicity) + ") ";
		for (size_t k = 0; k < f.centers.size(); k++) {
			if (k) s += "-";
			s += pad(constants::atnr2letter(atoms[f.centers[k]].Z), 2) + pad(f.centers[k] + 1, 3);
		}
		if (f.centers.size() == 1) s += " s";
		return collapse(s);
	}

	//"LP ( 1)Cl 2" / "BD*( 1) C 1- C 6" - the E2 table uses a different layout for the same thing.
	std::string e2_label(const NboFunction& f, const std::vector<NAOAtom>& atoms)
	{
		std::string s = f.type;
		while (s.size() < 3) s += ' ';
		s += "(" + pad(f.multiplicity, 2) + ")";
		for (size_t k = 0; k < f.centers.size(); k++) {
			if (k) s += "-";
			s += pad(constants::atnr2letter(atoms[f.centers[k]].Z), 2) + pad(f.centers[k] + 1, 3);
		}
		return collapse(s);
	}

	struct AtomIndices {
		ivec core, valence, rydberg, all;
	};

	std::vector<AtomIndices> atom_indices(const NAOResult& nao)
	{
		std::vector<AtomIndices> out(nao.atoms.size());
		for (size_t i = 0; i < nao.orbitals.size(); i++) {
			const NAO& o = nao.orbitals[i];
			AtomIndices& a = out[o.atom];
			a.all.push_back(static_cast<int>(i));
			if (o.type == NAOClass::Core) a.core.push_back(static_cast<int>(i));
			else if (o.type == NAOClass::Valence) a.valence.push_back(static_cast<int>(i));
			else a.rydberg.push_back(static_cast<int>(i));
		}
		return out;
	}

	//Leading eigenpair of the block of G over the given indices, embedded back into the full space.
	double leading_block(const MatrixXd& G, const ivec& idx, VectorXd& v)
	{
		const int k = static_cast<int>(idx.size());
		MatrixXd B(k, k);
		for (int i = 0; i < k; i++)
			for (int j = 0; j < k; j++)
				B(i, j) = G(idx[i], idx[j]);
		Eigen::SelfAdjointEigenSolver<MatrixXd> es(B);
		v = VectorXd::Zero(G.rows());
		for (int i = 0; i < k; i++) v(idx[i]) = es.eigenvectors()(i, k - 1);
		return es.eigenvalues()(k - 1);
	}

	//minority-centre share below which a two-centre candidate is no bond
	double bond_minority_floor()
	{
		static const double v = [] {
			const char* e = std::getenv("NBO_BOND_FLOOR");
			return e ? std::atof(e) : 0.15;
		}();
		return v;
	}

	double weight_on(const VectorXd& v, const ivec& idx)
	{
		double w = 0.0;
		for (const int i : idx) w += v(i) * v(i);
		return w;
	}

	MatrixXd corner(const MatrixXd& M, const ivec& idx)
	{
		const int m = static_cast<int>(idx.size());
		MatrixXd out(m, m);
		for (int p = 0; p < m; p++)
			for (int q = 0; q < m; q++)
				out(p, q) = M(idx[p], idx[q]);
		return out;
	}

	void scatter_corner(MatrixXd& M, const MatrixXd& block, const ivec& idx)
	{
		const int m = static_cast<int>(idx.size());
		for (int p = 0; p < m; p++)
			for (int q = 0; q < m; q++)
				M(idx[p], idx[q]) = block(p, q);
	}

	VectorXd gather(const VectorXd& v, const ivec& idx)
	{
		const int m = static_cast<int>(idx.size());
		VectorXd out(m);
		for (int p = 0; p < m; p++) out(p) = v(idx[p]);
		return out;
	}
}

dMatrix2 nao_density(const dMatrix2& P, const dMatrix2& S, const dMatrix2& C)
{
	const MatrixXd Ce = to_eigen(C), Se = to_eigen(S);
	return to_dmatrix(MatrixXd(Ce.transpose() * Se * to_eigen(P) * Se * Ce));
}

dMatrix2 nao_operator(const dMatrix2& F, const dMatrix2& C)
{
	const MatrixXd Ce = to_eigen(C);
	return to_dmatrix(MatrixXd(Ce.transpose() * to_eigen(F) * Ce));
}

NboInput read_file47(const std::filesystem::path& file)
{
	std::ifstream in(file);
	err_checkf(in.good(), "Cannot read FILE47: " + file.string(), std::cout);
	const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

	NboInput res;
	const size_t head = text.find("$GENNBO");
	err_checkf(head != std::string::npos, "No $GENNBO line in " + file.string(), std::cout);
	const std::string header = text.substr(head, text.find('\n', head) - head);
	const size_t nbas = header.find("NBAS=");
	err_checkf(nbas != std::string::npos, "No NBAS in the $GENNBO line of " + file.string(), std::cout);
	res.n = std::stoi(header.substr(nbas + 5));
	res.open_shell = header.find("OPEN") != std::string::npos;

	//$BASIS: one CENTER and one LABEL per basis function, shells contiguous
	const std::string basis = section_body(text, "$BASIS");
	const size_t label_at = basis.find("LABEL");
	err_checkf(label_at != std::string::npos, "No LABEL array in the $BASIS of " + file.string(), std::cout);
	const vec centers = numbers_of(basis.substr(0, label_at));
	const vec labels = numbers_of(basis.substr(label_at));
	err_checkf(static_cast<int>(centers.size()) == res.n && static_cast<int>(labels.size()) == res.n,
			   "The $BASIS of " + file.string() + " does not hold NBAS centres and labels", std::cout);
	res.ao.resize(res.n);
	std::map<std::pair<int, int>, int> shells_seen;
	for (int k = 0; k < res.n;) {
		int l = 0, component = 0;
		err_checkf(decode_label(static_cast<int>(std::llround(labels[k])), l, component),
				   "Unknown FILE47 LABEL code " + std::to_string(std::llround(labels[k])), std::cout);
		err_checkf(component == 0, "FILE47 shell does not start at its first component", std::cout);
		const int atom = static_cast<int>(std::llround(centers[k])) - 1;
		const int shell = shells_seen[{ atom, l }]++;
		for (int c = 0; c <= 2 * l; c++, k++) {
			err_checkf(k < res.n, "FILE47 ends inside a shell", std::cout);
			res.ao[k].atom = atom;
			res.ao[k].l = l;
			res.ao[k].shell = shell;
			res.ao[k].m = c;
		}
	}

	const size_t packed = static_cast<size_t>(res.n) * (res.n + 1) / 2;
	const vec ovlp = numbers_of(section_body(text, "$OVERLAP"));
	err_checkf(ovlp.size() >= packed, "$OVERLAP of " + file.string() + " is short", std::cout);
	res.overlap = unpack(ovlp, 0, res.n);
	const int blocks = res.open_shell ? 2 : 1;
	const vec dens = numbers_of(section_body(text, "$DENSITY"));
	err_checkf(dens.size() >= packed * blocks, "$DENSITY of " + file.string() + " is short", std::cout);
	for (int b = 0; b < blocks; b++) res.density.push_back(unpack(dens, packed * b, res.n));
	const vec fock = numbers_of(section_body(text, "$FOCK"));
	if (fock.size() >= packed * blocks)
		for (int b = 0; b < blocks; b++) res.fock.push_back(unpack(fock, packed * b, res.n));
	return res;
}

bvec2 bondable_pairs(const std::vector<atom>& atoms, const double scale)
{
	const size_t na = atoms.size();
	bvec2 out(na, bvec(na, false));
	for (size_t a = 0; a < na; a++)
		for (size_t b = a + 1; b < na; b++) {
			double d2 = 0.0;
			for (unsigned int k = 0; k < 3; k++) {
				const double d = atoms[a].get_coordinate(k) - atoms[b].get_coordinate(k);
				d2 += d * d;
			}
			const int za = atoms[a].get_charge(), zb = atoms[b].get_charge();
			const double ra = za > 0 && za < 114 ? constants::covalent_radii[za] : 1.5;
			const double rb = zb > 0 && zb < 114 ? constants::covalent_radii[zb] : 1.5;
			const double cut = constants::ang2bohr(scale * (ra + rb));
			out[a][b] = out[b][a] = d2 <= cut * cut;
		}
	return out;
}

namespace {
	struct ProfClock {
		std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
		double lap() {
			const auto now = std::chrono::steady_clock::now();
			const double d = std::chrono::duration<double>(now - t0).count();
			t0 = now;
			return d;
		}
	};
}

NboLewis nbo_search(const NAOResult& nao, const dMatrix2& gamma, const bvec2& bondable,
					const int n_pairs, const double scale, const NboOptions& options, std::ostream& log)
{
	ProfClock prof;
	double p_ladder = 0, p_scf = 0, p_owso = 0, p_anti = 0, p_comp = 0, p_final = 0;
	double p_pivot = 0, p_ryd = 0;
	int n_sweeps = 0, n_eig1 = 0, n_eig2 = 0, n_levels_used = 0;
	vec ch_trace;  //convergence history, printed under -debug
	const int n = static_cast<int>(gamma.extent(0));
	const int natoms = static_cast<int>(nao.atoms.size());
	const MatrixXd G0 = to_eigen(gamma);
	MatrixXd G = G0;
	const std::vector<AtomIndices> idx = atom_indices(nao);
#ifdef _OPENMP
	const int nthreads = options.search_threads > 0 ? options.search_threads
					   : options.threads > 0        ? options.threads
													: omp_get_max_threads();
#else
	const int nthreads = 1;
#endif

	NboLewis res;
	res.gamma = gamma;
	res.topo.assign(natoms, ivec(natoms, 0));
	std::vector<VectorXd> vectors;
	auto accept = [&](const VectorXd& v, const std::string& type, const ivec& centers) {
		NboFunction f;
		f.type = type;
		f.centers = centers;
		f.multiplicity = 1;
		for (const NboFunction& o : res.orbitals)
			if (o.type == type && o.centers == centers) f.multiplicity++;
		f.occupancy = v.dot(G0 * v);
		res.orbitals.push_back(f);
		vectors.push_back(v);
		//deplete, so the next block search sees only what is left
		const double occ = v.dot(G * v);
		G -= occ * v * v.transpose();
	};

	//A core NAO is already a one-centre orbital of occupancy ~2, kept as it is.
	for (int a = 0; a < natoms; a++)
		for (const int i : idx[a].core) {
			VectorXd v = VectorXd::Zero(n);
			v(i) = 1.0;
			accept(v, "CR", { a });
			res.topo[a][a]++;
		}

	//threshold ladder: one-centre, then two-centre blocks at each threshold
	static const vec ladder = { 1.90, 1.80, 1.70, 1.60, 1.50, 1.40, 1.30, 1.20, 1.10,
								1.00, 0.90, 0.80, 0.70, 0.60, 0.50 };
	ivec used(natoms, 0);
	ivec cap(natoms, 0);
	for (int a = 0; a < natoms; a++) cap[a] = static_cast<int>(idx[a].valence.size());
	for (int relax = 0; relax < 2; relax++) {
	  if (relax && static_cast<int>(res.orbitals.size()) >= n_pairs) break;
	  for (const double t0 : ladder) {
		const double t = std::min(t0, options.occupancy_threshold) * scale / 2.0;
		bool progress = true;
		while (progress && static_cast<int>(res.orbitals.size()) < n_pairs) {
			progress = false;
			//Candidate eigensolves at one threshold are independent.
			vec lam1(natoms, 0.0);
			ivec have1(natoms, 0);
			std::vector<VectorXd> cand1(natoms);
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
			for (int a = 0; a < natoms; a++) {
				if (idx[a].valence.empty()) continue;
				if (!relax && used[a] >= cap[a]) continue;
				lam1[a] = leading_block(G, idx[a].all, cand1[a]);
				have1[a] = 1;
			}
			double best = t;
			VectorXd bv;
			int ba = -1;
			for (int a = 0; a < natoms; a++) {
				if (!have1[a]) continue;
				n_eig1++;
				if (lam1[a] > best) { best = lam1[a]; bv = cand1[a]; ba = a; }
			}
			if (ba >= 0) {
				accept(bv, "LP", { ba });
				res.topo[ba][ba]++;
				used[ba]++;
				res.threshold = t;
				progress = true;
				continue;
			}
			double bestp = t;
			VectorXd bvp;
			int pa = -1, pb = -1;
			//O(N^2) small diagonalisations per accepted bond
			ivec2 plist;
			for (int a = 0; a < natoms; a++) {
				if (idx[a].valence.empty()) continue;
				if (!relax && used[a] >= cap[a]) continue;
				for (int b = a + 1; b < natoms; b++) {
					if (idx[b].valence.empty()) continue;
					if (!relax && used[b] >= cap[b]) continue;
					if (!bondable.empty() && !bondable[a][b]) continue;
					plist.push_back({ a, b });
				}
			}
			const int np = static_cast<int>(plist.size());
			vec lam2(np, 0.0);
			std::vector<VectorXd> cand2(np);
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
			for (int c = 0; c < np; c++) {
				ivec pair = idx[plist[c][0]].all;
				pair.insert(pair.end(), idx[plist[c][1]].all.begin(), idx[plist[c][1]].all.end());
				lam2[c] = leading_block(G, pair, cand2[c]);
			}
			for (int c = 0; c < np; c++) {
				n_eig2++;
				if (lam2[c] <= bestp) continue;
				const double wa = weight_on(cand2[c], idx[plist[c][0]].all);
				const double wb = weight_on(cand2[c], idx[plist[c][1]].all);
				if (std::min(wa, wb) < bond_minority_floor() * (wa + wb)) continue;
				bestp = lam2[c];
				bvp = cand2[c];
				pa = plist[c][0];
				pb = plist[c][1];
			}
			if (pa >= 0) {
				accept(bvp, "BD", { pa, pb });
				res.topo[pa][pb]++;
				res.topo[pb][pa]++;
				used[pa]++;
				used[pb]++;
				res.threshold = t;
				progress = true;
			}
		}
		if (static_cast<int>(res.orbitals.size()) >= n_pairs) break;
	  }
	}
	res.n_lewis = static_cast<int>(res.orbitals.size());
	{
		double tr = 0.0;
		for (int i = 0; i < n; i++) tr += G0(i, i);
		err_checkf(res.n_lewis > 0,
				   "NBO: the search found no orbital at all - this spin's density holds only " +
				   std::to_string(tr) + " electrons over " + std::to_string(n) +
				   " NAOs, so there is no Lewis structure in it.  Check the wavefunction the .47 "
				   "was built from: an unnormalised or truncated density reads like this one.",
				   std::cout);
	}
	p_ladder = prof.lap();

	//Refine the Lewis orbitals self-consistently after the threshold search.
	{
		//The orbital centres and their density corners are fixed during the sweep.
		std::vector<ivec> sub(res.n_lewis);
		std::vector<MatrixXd> G0_corner(res.n_lewis);
		for (int j = 0; j < res.n_lewis; j++) {
			for (const int a : res.orbitals[j].centers)
				sub[j].insert(sub[j].end(), idx[a].all.begin(), idx[a].all.end());
			G0_corner[j] = corner(G0, sub[j]);
			//Corner updates require vectors supported only on their assigned centres.
			bvec on(n, false);
			for (const int i : sub[j]) on[i] = true;
			for (int i = 0; i < n; i++)
				err_checkf(on[i] || vectors[j](i) == 0.0,
						   "NBO: a search vector reaches outside its own centres", std::cout);
		}
		vec occ(res.n_lewis, 0.0);
		MatrixXd sum = MatrixXd::Zero(n, n);
		for (int j = 0; j < res.n_lewis; j++) {
			occ[j] = vectors[j].dot(G0 * vectors[j]);
			const VectorXd vs = gather(vectors[j], sub[j]);
			MatrixXd block = corner(sum, sub[j]);
			block += occ[j] * vs * vs.transpose();
			scatter_corner(sum, block, sub[j]);
		}
		//orbitals on disjoint centres update disjoint corners, so one level runs in parallel
		ivec level(res.n_lewis, 0);
		int n_levels = 0;
		for (int j = 0; j < res.n_lewis; j++) {
			//a core orbital is a single NAO and stays one, so it never enters a sweep
			if (res.orbitals[j].type == "CR") continue;
			for (int i = 0; i < j; i++) {
				if (res.orbitals[i].type == "CR") continue;
				bool share = false;
				for (const int a : res.orbitals[j].centers)
					for (const int b : res.orbitals[i].centers)
						if (a == b) share = true;
				if (share) level[j] = std::max(level[j], level[i] + 1);
			}
			n_levels = std::max(n_levels, level[j] + 1);
		}
		std::vector<ivec> by_level(n_levels);
		n_levels_used = n_levels;
		for (int j = 0; j < res.n_lewis; j++)
			if (res.orbitals[j].type != "CR") by_level[level[j]].push_back(j);
		//per-orbital slot, max taken afterwards: MSVC's OpenMP 2.0 has no max reduction
		vec moved_at(res.n_lewis, 0.0);
		for (int sweep = 0; sweep < 200; sweep++) {
			for (int L = 0; L < n_levels; L++) {
				const ivec& lv = by_level[L];
				const int n_lv = static_cast<int>(lv.size());
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
				for (int t = 0; t < n_lv; t++) {
					const int j = lv[t];
					const ivec& s = sub[j];
					const int m = static_cast<int>(s.size());
					//one gather and scatter serve both rank-one updates and the eigensolve's difference
					MatrixXd block = corner(sum, s);
					VectorXd vs = gather(vectors[j], s);
					block -= occ[j] * vs * vs.transpose();
					Eigen::SelfAdjointEigenSolver<MatrixXd> es(MatrixXd(G0_corner[j] - block));
					VectorXd v = VectorXd::Zero(n);
					for (int p = 0; p < m; p++) v(s[p]) = es.eigenvectors()(p, m - 1);
					if (v.dot(vectors[j]) < 0.0) v = -v;
					moved_at[j] = (v - vectors[j]).norm();
					vectors[j] = v;
					vs = gather(v, s);
					//v vanishes outside s, so its occupancy against G0 is the corner's
					occ[j] = vs.dot(G0_corner[j] * vs);
					block += occ[j] * vs * vs.transpose();
					scatter_corner(sum, block, s);
				}
			}
			double change = 0.0;
			for (int j = 0; j < res.n_lewis; j++) change = std::max(change, moved_at[j]);
			n_sweeps++;
			ch_trace.push_back(change);
			if (change < 1e-10) break;
		}
		p_scf = prof.lap();
		//Occupancy-weighted orthogonalization removes overlap among Lewis orbitals.
		MatrixXd VL(n, res.n_lewis);
		for (int j = 0; j < res.n_lewis; j++) VL.col(j) = vectors[j];
		VectorXd w(res.n_lewis);
		for (int j = 0; j < res.n_lewis; j++) w(j) = std::max(occ[j], 1e-6);
		const MatrixXd S = VL.transpose() * VL;
		VL = VL * w.asDiagonal() *
			 sym_power(MatrixXd(w.asDiagonal() * S * w.asDiagonal()), -0.5);
		for (int j = 0; j < res.n_lewis; j++) vectors[j] = VL.col(j);
		p_owso = prof.lap();
	}

	//A bond c_A h_A + c_B h_B has the antibond c_B h_A - c_A h_B; constructed, not searched, as it
	//has almost no occupancy to find it by.
	std::vector<VectorXd> non_lewis;
	std::vector<NboFunction> non_lewis_fn;
	for (int j = 0; j < res.n_lewis; j++) {
		if (res.orbitals[j].type != "BD") continue;
		const int a = res.orbitals[j].centers[0], b = res.orbitals[j].centers[1];
		VectorXd va = VectorXd::Zero(n), vb = VectorXd::Zero(n);
		for (const int i : idx[a].all) va(i) = vectors[j](i);
		for (const int i : idx[b].all) vb(i) = vectors[j](i);
		const double ca = va.norm(), cb = vb.norm();
		if (ca < 1e-8 || cb < 1e-8) continue;
		VectorXd w = (cb / ca) * va - (ca / cb) * vb;
		const double nrm = w.norm();
		if (nrm < 1e-8) continue;
		NboFunction f;
		f.type = "BD*";
		f.centers = { a, b };
		f.multiplicity = res.orbitals[j].multiplicity;
		non_lewis.push_back(w / nrm);
		non_lewis_fn.push_back(f);
	}
	//Each antibond is orthogonal to its own bond only; nbo_e2 and the occupancies need an orthonormal set,
	//so take the (orthonormal) Lewis span out of them and orthonormalize them symmetrically
	if (!non_lewis.empty()) {
		MatrixXd VL(n, res.n_lewis), A(n, static_cast<int>(non_lewis.size()));
		for (int j = 0; j < res.n_lewis; j++) VL.col(j) = vectors[j];
		for (size_t j = 0; j < non_lewis.size(); j++) A.col(static_cast<int>(j)) = non_lewis[j];
		A -= VL * (VL.transpose() * A);
		A = A * sym_power(MatrixXd(A.transpose() * A), -0.5);
		for (size_t j = 0; j < non_lewis.size(); j++) non_lewis[j] = A.col(static_cast<int>(j));
	}
	p_anti = prof.lap();

	//Project the Rydberg space out of the Lewis and antibond subspaces.
	{
		const int k = res.n_lewis + static_cast<int>(non_lewis.size());
		MatrixXd M(n, k);
		for (int j = 0; j < res.n_lewis; j++) M.col(j) = vectors[j];
		for (size_t j = 0; j < non_lewis.size(); j++)
			M.col(res.n_lewis + static_cast<int>(j)) = non_lewis[j];
		const MatrixXd Q = MatrixXd::Identity(n, n) -
			M * sym_power(MatrixXd(M.transpose() * M), -1.0) * M.transpose();
		//Use local basis vectors for deterministic Rydberg labels.
		p_comp = prof.lap();
		std::vector<VectorXd> extra;
		ivec owner;
		std::vector<VectorXd> residual(n);
		for (int i = 0; i < n; i++) residual[i] = Q.col(i);
		bvec used(n, false);
		//Orthogonalize the residuals without repeated full projections.
		vec nrm(n, 0.0);
		for (int step = 0; step < n - k; step++) {
			int pivot = -1;
			double best = 0.0;
#pragma omp parallel for schedule(static) num_threads(nthreads)
			for (int i = 0; i < n; i++)
				if (!used[i]) nrm[i] = residual[i].norm();
			for (int i = 0; i < n; i++) {
				if (used[i]) continue;
				if (nrm[i] > best) { best = nrm[i]; pivot = i; }
			}
			err_checkf(pivot >= 0 && best > 1e-8,
					   "NBO: the residual space is smaller than the orbital count demands", std::cout);
			used[pivot] = true;
			const VectorXd u = residual[pivot] / best;
			extra.push_back(u);
			owner.push_back(nao.orbitals[pivot].atom);
#pragma omp parallel for schedule(static) num_threads(nthreads)
			for (int i = 0; i < n; i++)
				if (!used[i]) residual[i] -= u.dot(residual[i]) * u;
		}
		p_pivot = prof.lap();
		for (int a = 0; a < natoms; a++) {
			ivec mine;
			for (size_t j = 0; j < extra.size(); j++)
				if (owner[j] == a) mine.push_back(static_cast<int>(j));
			if (mine.empty()) continue;
			const int kk = static_cast<int>(mine.size());
			MatrixXd B(n, kk);
			for (int j = 0; j < kk; j++) B.col(j) = extra[mine[j]];
			Eigen::SelfAdjointEigenSolver<MatrixXd> e2(MatrixXd(B.transpose() * G0 * B));
			for (int j = 0; j < kk; j++) {
				const VectorXd v = B * e2.eigenvectors().col(kk - 1 - j);
				NboFunction f;
				//still in the atom's valence shell: lone vacancy; pushed out of it: Rydberg
				f.type = (weight_on(v, idx[a].valence) > 0.5) ? "LV" : "RY";
				f.centers = { a };
				non_lewis.push_back(v);
				non_lewis_fn.push_back(f);
			}
		}
		p_ryd = prof.lap();
		p_comp += p_pivot + p_ryd;
	}

	for (size_t j = 0; j < non_lewis.size(); j++) {
		NboFunction f = non_lewis_fn[j];
		//occupancy is computed below from the same vector
		if (f.type != "BD*") {
			f.multiplicity = 1;
			for (const NboFunction& o : res.orbitals)
				if (o.type == f.type && o.centers == f.centers) f.multiplicity++;
		}
		res.orbitals.push_back(f);
		vectors.push_back(non_lewis[j]);
	}

	//Compute occupations and polarization from the final NBO vectors.
	const int n_orb = static_cast<int>(res.orbitals.size());
#pragma omp parallel for schedule(static) num_threads(nthreads)
	for (int j = 0; j < n_orb; j++) {
		NboFunction& f = res.orbitals[j];
		f.occupancy = vectors[j].dot(G0 * vectors[j]);
		f.coefficients.assign(vectors[j].data(), vectors[j].data() + n);
		for (const int a : f.centers) {
			double wa = 0.0;
			vec lchar(4, 0.0);
			for (const int i : idx[a].all) {
				const double c2 = vectors[j](i) * vectors[j](i);
				wa += c2;
				lchar[std::min(nao.orbitals[i].l, 3)] += c2;
			}
			f.center_weight.push_back(wa);
			for (double& x : lchar) x = (wa > 1e-300) ? x / wa : 0.0;
			f.center_lchar.push_back(lchar);
		}
	}

	double lewis_density = 0.0;
	for (int j = 0; j < res.n_lewis; j++) lewis_density += res.orbitals[j].occupancy;
	double total = 0.0;
	for (int i = 0; i < n; i++) total += G0(i, i);
	res.rho_nl = total - lewis_density;
	p_final = prof.lap();
	if (options.debug) {
		log << "NBO search: n=" << n << " atoms=" << natoms << " lewis=" << res.n_lewis
				  << " orbitals=" << res.orbitals.size() << " sweeps=" << n_sweeps
				  << " levels=" << n_levels_used << " eig=" << n_eig1 + n_eig2
				  << " | ladder " << p_ladder << " sweep " << p_scf << " owso " << p_owso
				  << " antibonds " << p_anti << " complement " << p_comp << " (pivot "
				  << p_pivot << " rydberg " << p_ryd << ") coefficients " << p_final
				  << std::endl;
		//Report a sweep that reaches the iteration cap.
		log << "NBO search: sweep change";
		for (size_t i = 0; i < ch_trace.size(); i++)
			if (i < 3 || i + 3 >= ch_trace.size())
				log << " " << i << ":" << ch_trace[i];
		log << std::endl;
	}
	return res;
}

std::vector<NboE2Entry> nbo_e2(NboLewis& lewis, const dMatrix2& fock_nao,
							   const double threshold_kcal)
{
	std::vector<NboE2Entry> out;
	if (fock_nao.extent(0) == 0) return out;
	const int n = static_cast<int>(lewis.orbitals.size());
	const MatrixXd F = to_eigen(fock_nao);
	MatrixXd V(F.rows(), n);
	for (int j = 0; j < n; j++)
		for (int i = 0; i < F.rows(); i++)
			V(i, j) = lewis.orbitals[j].coefficients[i];
	const MatrixXd Fnbo = V.transpose() * F * V;
	for (int i = 0; i < n; i++) lewis.orbitals[i].energy = Fnbo(i, i);
	for (int i = 0; i < lewis.n_lewis; i++) {
		for (int j = lewis.n_lewis; j < n; j++) {
			const double de = Fnbo(j, j) - Fnbo(i, i);
			if (std::abs(de) < 1e-8) continue;
			const double f = Fnbo(i, j);
			//E(2) = -q_i F(i,j)^2 / (eps_j - eps_i), printed as a stabilisation in kcal/mol
			const double e2 = lewis.orbitals[i].occupancy * f * f / de * constants::kcal_mol_per_hartree;
			if (e2 < threshold_kcal) continue;
			NboE2Entry entry;
			entry.donor_index = i + 1;
			entry.acceptor_index = j + 1;
			entry.energy_kcal = e2;
			entry.e_diff = de;
			entry.fij = f;
			out.push_back(entry);
		}
	}
	return out;
}

namespace
{
	//One spin's worth of the analysis, appended to the results.
	void analyse_spin(NboResults& res, const NboInput& in, const NAOResult& nao,
					  const bvec2& bondable, const bvec2& nrt_bondable, const dMatrix2& density,
					  const dMatrix2& fock,
					  const int n_pairs, const double scale, const std::string& spin,
					  const NboOptions& options, std::vector<NboLewis>& out_lewis, std::ostream& log)
	{
		const dMatrix2 gamma = nao_density(density, in.overlap, nao.C);
		const dMatrix2 fock_nao = fock.extent(0) ? nao_operator(fock, nao.C) : dMatrix2();
		const auto clock = [] { return std::chrono::steady_clock::now(); };
		auto t = clock();
		NboLewis lewis = nbo_search(nao, gamma, bondable, n_pairs, scale, options, log);
		res.search_seconds += std::chrono::duration<double>(clock() - t).count();
		t = clock();
		std::vector<NboE2Entry> e2 = nbo_e2(lewis, fock_nao, options.e2_threshold_kcal);
		res.e2_seconds += std::chrono::duration<double>(clock() - t).count();
		//before the renumbering below, while e2's indices still point into lewis.orbitals
		if (options.nrt)
			//log, not std::cout: under -fba this runs on the NBO thread and cout is RGBI's
			native_nrt(res.nrt, nao, lewis, e2, nrt_bondable, options, spin, scale, log);

		//NBO numbers by type group: the Lewis set in found order, then LV, BD*, RY
		ivec order;
		static const char* groups[] = { "CR", "LP", "BD", "LV", "BD*", "RY" };
		for (const char* g : groups)
			for (size_t j = 0; j < lewis.orbitals.size(); j++)
				if (lewis.orbitals[j].type == g) order.push_back(static_cast<int>(j));
		err_checkf(order.size() == lewis.orbitals.size(), "NBO: an orbital carries an unknown type",
				   std::cout);
		ivec position(lewis.orbitals.size(), 0);
		const int base = static_cast<int>(res.orbitals.size());
		for (size_t k = 0; k < order.size(); k++) position[order[k]] = static_cast<int>(k);

		for (const int j : order) {
			const NboFunction& f = lewis.orbitals[j];
			NboOrbital o;
			o.index = base + position[j] + 1;
			o.type = f.type;
			o.spin = spin;
			o.description = orbital_description(f, nao.atoms);
			o.occupancy = f.occupancy;
			o.energy = f.energy;
			for (size_t k = 0; k < f.centers.size(); k++) {
				NboHybrid h;
				h.center = f.centers[k] + 1;
				h.element = constants::atnr2letter(nao.atoms[f.centers[k]].Z);
				h.weight_percent = 100.0 * f.center_weight[k];
				//NBO's sign: first centre positive, so an antibond's second centre carries the minus
				const double c = std::sqrt(std::max(f.center_weight[k], 0.0));
				h.coefficient = (k == 1 && f.type == "BD*") ? -c : c;
				h.s = 100.0 * f.center_lchar[k][0];
				h.p = 100.0 * f.center_lchar[k][1];
				h.d = 100.0 * f.center_lchar[k][2];
				h.f = 100.0 * f.center_lchar[k][3];
				o.centers.push_back(f.centers[k] + 1);
				o.hybrids.push_back(h);
			}
			res.orbitals.push_back(o);
		}
		for (NboE2Entry& e : e2) {
			e.spin = spin;
			e.donor = e2_label(lewis.orbitals[e.donor_index - 1], nao.atoms);
			e.acceptor = e2_label(lewis.orbitals[e.acceptor_index - 1], nao.atoms);
			e.donor_index = base + position[e.donor_index - 1] + 1;
			e.acceptor_index = base + position[e.acceptor_index - 1] + 1;
			res.e2.push_back(e);
		}
		out_lewis.push_back(lewis);
	}

	//NAOs ordered by atom, l and NBO's component order
	void fill_nao_table(std::vector<NboNao>& out, const NAOResult& nao, const vec& occupancy,
						const dMatrix2& fock_nao)
	{
		ivec order(nao.orbitals.size());
		for (size_t i = 0; i < order.size(); i++) order[i] = static_cast<int>(i);
		std::stable_sort(order.begin(), order.end(), [&](const int a, const int b) {
			const NAO& x = nao.orbitals[a];
			const NAO& y = nao.orbitals[b];
			if (x.atom != y.atom) return x.atom < y.atom;
			if (x.l != y.l) return x.l < y.l;
			const int rx = lang_rank(x.l, x.m), ry = lang_rank(y.l, y.m);
			if (rx != ry) return rx < ry;
			return x.shell < y.shell;
		});
		int index = 0;
		for (const int i : order) {
			const NAO& o = nao.orbitals[i];
			NboNao n;
			n.index = ++index;
			n.center = o.atom + 1;
			n.element = constants::atnr2letter(nao.atoms[o.atom].Z);
			n.lang = lang_label(o.l, o.m);
			n.type = (o.type == NAOClass::Core) ? "Cor"
				   : (o.type == NAOClass::Valence) ? "Val" : "Ryd";
			n.shell = std::to_string(o.n) + std::string(1, "spdfghik"[std::min(o.l, 7)]);
			n.occupancy = occupancy[i];
			if (fock_nao.extent(0)) n.energy = fock_nao(i, i);
			out.push_back(n);
		}
	}
}

NboResults native_nbo(WFN& wavy, const NboOptions& options, std::ostream& log)
{
	const auto clock = [] { return std::chrono::steady_clock::now(); };
	const auto t_f47 = clock();
	std::filesystem::path f47 = options.file47;
	bool temporary = false;
	if (f47.empty()) {
		f47 = std::filesystem::temp_directory_path() /
			  ("nbo_native_" + std::to_string(static_cast<unsigned long long>(
				   std::chrono::steady_clock::now().time_since_epoch().count())) + ".47");
		err_checkf(wavy.write_nbo(f47, options.debug, options.debug ? &log : nullptr, ""),
				   "Could not write the FILE47 the native NBO analysis reads", log);
		temporary = !options.keep_file47;
	}
	const NboInput in = read_file47(f47);
	if (temporary) std::filesystem::remove(f47);

	NboResults res;
	res.file47_seconds = std::chrono::duration<double>(clock() - t_f47).count();
	res.name = wavy.get_path().stem().string();
	res.source = wavy.get_path().string();
	res.version = "NoSpherA2 native NBO";
	res.open_shell = in.open_shell;
	res.e2_threshold_kcal = options.e2_threshold_kcal;

	std::vector<atom> atoms = wavy.get_atoms();
	const bvec2 bondable = bondable_pairs(atoms);
	const bvec2 nrt_bondable = bondable_pairs(atoms, options.nrt_bond_scale);
	ivec ecp(atoms.size(), 0);
	if (wavy.get_has_ECPs())
		for (size_t a = 0; a < atoms.size(); a++)
			ecp[a] = wavy.get_atom_ECP_electrons(static_cast<int>(a));

	std::vector<NboLewis> lewis;
	if (!in.open_shell) {
		const auto t_nao = clock();
		const NAOResult nao = build_naos(in.density[0], in.overlap, in.ao, atoms, ecp);
		res.nao_seconds = std::chrono::duration<double>(clock() - t_nao).count();
		double electrons = 0.0;
		for (const NAOAtom& a : nao.atoms) electrons += a.population;
		const int n_pairs = static_cast<int>(std::llround(electrons / 2.0));
		analyse_spin(res, in, nao, bondable, nrt_bondable, in.density[0],
					 in.fock.empty() ? dMatrix2() : in.fock[0], n_pairs, 2.0, "", options, lewis, log);
		vec occ(nao.orbitals.size(), 0.0);
		for (size_t i = 0; i < occ.size(); i++) occ[i] = nao.orbitals[i].occupation;
		fill_nao_table(res.nao, nao, occ, in.fock.empty() ? dMatrix2() : nao_operator(in.fock[0], nao.C));
		for (const NAOAtom& a : nao.atoms) {
			NboAtomPopulation p;
			p.element = constants::atnr2letter(a.Z);
			p.index = a.index + 1;
			p.charge = a.charge;
			p.core = a.core;
			p.valence = a.valence;
			p.rydberg = a.rydberg;
			p.total = a.population;
			res.npa.push_back(p);
		}
	}
	else {
		//NPA tables share total-density NAOs; NBO searches use separate spin densities.
		const auto t_nao = clock();
		const dMatrix2 total_density = to_dmatrix(MatrixXd(to_eigen(in.density[0]) + to_eigen(in.density[1])));
		const dMatrix2 spin_density = to_dmatrix(MatrixXd(to_eigen(in.density[0]) - to_eigen(in.density[1])));
		const NAOResult total_nao = build_naos(total_density, in.overlap, in.ao, atoms, ecp);
		const NAOResult a_nao = build_naos(in.density[0], in.overlap, in.ao, atoms, ecp);
		const NAOResult b_nao = build_naos(in.density[1], in.overlap, in.ao, atoms, ecp);
		res.nao_seconds = std::chrono::duration<double>(clock() - t_nao).count();
		for (int s = 0; s < 2; s++) {
			const NAOResult& nao = s ? b_nao : a_nao;
			double electrons = 0.0;
			for (const NAOAtom& a : nao.atoms) electrons += a.population;
			analyse_spin(res, in, nao, bondable, nrt_bondable, in.density[s],
						 static_cast<int>(in.fock.size()) > s ? in.fock[s] : dMatrix2(),
						 static_cast<int>(std::llround(electrons)), 1.0, s ? "beta" : "alpha",
						 options, lewis, log);
		}
		vec occ(total_nao.orbitals.size(), 0.0);
		for (size_t i = 0; i < occ.size(); i++)
			occ[i] = total_nao.orbitals[i].occupation;
		fill_nao_table(res.nao, total_nao, occ,
					   in.fock.empty() ? dMatrix2() : nao_operator(in.fock[0], total_nao.C));
		for (int s = 0; s < 2; s++) {
			const dMatrix2 spin_nao = nao_density(in.density[s], in.overlap, total_nao.C);
			vec spin_occ(total_nao.orbitals.size(), 0.0);
			for (size_t i = 0; i < spin_occ.size(); i++) spin_occ[i] = spin_nao(i, i);
			fill_nao_table(s ? res.nao_beta : res.nao_alpha, total_nao, spin_occ,
						   static_cast<int>(in.fock.size()) > s ? nao_operator(in.fock[s], total_nao.C)
																: dMatrix2());
		}
		const dMatrix2 spin_nao = nao_density(spin_density, in.overlap, total_nao.C);
		for (const NAOAtom& x : total_nao.atoms) {
			NboAtomPopulation p;
			p.element = constants::atnr2letter(x.Z);
			p.index = x.index + 1;
			p.core = x.core;
			p.valence = x.valence;
			p.rydberg = x.rydberg;
			p.total = x.population;
			p.charge = x.charge;
			for (size_t i = 0; i < total_nao.orbitals.size(); i++)
				if (total_nao.orbitals[i].atom == x.index) p.spin_density += spin_nao(i, i);
			p.has_spin_density = true;
			res.npa.push_back(p);
		}
	}
	if (options.debug) {
		for (const NboLewis& l : lewis)
			log << "native NBO: " << l.n_lewis << " Lewis orbitals, rho(NL) = " << l.rho_nl
				<< ", threshold " << l.threshold << std::endl;
		log << "native NBO phases: file47 " << res.file47_seconds << " s, NAO " << res.nao_seconds
			<< " s, search " << res.search_seconds << " s, E2 " << res.e2_seconds
			<< " s (NRT reports its own)" << std::endl;
	}
	return res;
}

namespace {
	//"BD O1-H2" from the hybrids, which a parsed gennbo result carries too; the description is the fallback
	std::string nbo_label(const NboOrbital& o)
	{
		if (o.hybrids.empty()) return o.description;
		std::string s = o.type + " ";
		for (size_t k = 0; k < o.hybrids.size(); k++)
			s += (k ? "-" : "") + o.hybrids[k].element + std::to_string(o.hybrids[k].center);
		return s;
	}

	//The spins a table carries, in the order they first appear; one empty entry for a closed shell.
	template <class Row>
	svec spins_of(const std::vector<Row>& rows)
	{
		svec spins;
		for (const Row& x : rows)
			if (std::find(spins.begin(), spins.end(), x.spin) == spins.end()) spins.push_back(x.spin);
		return spins;
	}

	//"Resonance weights (alpha spin):" - the spin goes into the title, not into a column of every row
	std::string title(const std::string& what, const std::string& spin, std::string note = "")
	{
		if (!spin.empty())
			note += (note.empty() ? "" : ", ") + spin + (spin == "alpha" || spin == "beta" ? " spin" : "");
		return "\n" + what + (note.empty() ? "" : " (" + note + ")") + ":\n";
	}

	//a nearly empty orbital comes out at either sign, and "-0.00000" reads as a negative occupancy
	double printable(const double x)
	{
		return std::abs(x) < 5e-6 ? 0.0 : x;
	}
}

void print_nrt(const NboResults& r, std::ostream& out)
{
	using namespace std;
	const NboNrt& n = r.nrt;
	if (!n.present) return;
	const ostream_format_guard restore_format(out);
	//bond orders carry atom indices only; the populations table is where the elements are
	std::map<int, std::string> element;
	for (const NboAtomPopulation& p : r.npa) element[p.index] = p.element;
	for (const NboValency& v : n.valencies) element[v.atom] = v.element;
	const auto label = [&element](const int a) {
		const auto it = element.find(a);
		return (it == element.end() ? std::string("?") : it->second) + std::to_string(a);
	};

	out << "\nNatural resonance theory (in house):\n"
		<< "  " << n.structures_used << " of " << n.structures_found
		<< " resonance structures carry the fit, D(0) = " << fixed << setprecision(5) << n.d_0
		<< ", D(w) = " << n.d_w << "\n";
	//an open shell runs the search once per spin, and both runs push the same budget line
	svec seen;
	for (const std::string& s : n.notes)
		if (std::find(seen.begin(), seen.end(), s) == seen.end()) {
			seen.push_back(s);
			out << "  " << s << "\n";
		}

	for (const std::string& spin : spins_of(n.weights)) {
		out << title("Resonance weights", spin) << "  structure   weight %   change from the parent\n";
		for (const NboResonanceWeight& w : n.weights)
			if (w.spin == spin)
				out << setw(11) << w.structure << setprecision(2) << setw(11) << w.weight_percent << "   "
					<< w.changes << "\n";
	}
	//the bond-order diagonal is the atom's lone-pair count, so it goes into the valency table
	std::map<std::pair<std::string, int>, double> lone_pairs;
	for (const NboBondOrder& b : n.bond_orders)
		if (b.diagonal) lone_pairs[{ b.spin, b.atom1 }] = b.total;
	for (const std::string& spin : spins_of(n.bond_orders)) {
		out << title("Natural bond orders", spin) << "  atom   atom        total  covalent     ionic\n"
			<< setprecision(4);
		for (const NboBondOrder& b : n.bond_orders)
			if (b.spin == spin && !b.diagonal)
				out << "  " << left << setw(7) << label(b.atom1) << setw(7) << label(b.atom2) << right
					<< setw(10) << b.total << setw(10) << b.covalent << setw(10) << b.ionic << "\n";
	}
	for (const std::string& spin : spins_of(n.valencies)) {
		out << title("Natural atomic valencies", spin)
			<< "  atom      valency  covalency  electrovalency  electrons  lone pairs\n" << setprecision(4);
		for (const NboValency& v : n.valencies) {
			if (v.spin != spin) continue;
			out << "  " << left << setw(7) << (v.element + std::to_string(v.atom)) << right << setw(10)
				<< v.valency << setw(11) << v.covalency << setw(16) << v.electrovalency << setw(11)
				<< v.electron_count;
			const auto lp = lone_pairs.find({ v.spin, v.atom });
			if (lp != lone_pairs.end()) out << setw(12) << lp->second;
			out << "\n";
		}
	}
	for (const std::string& s : n.symmetry_forms) out << "  " << s << "\n";
}

void print_nbo(const NboResults& r, std::ostream& out)
{
	using namespace std;
	//fixed/setprecision would otherwise outlive this table
	const ostream_format_guard restore_format(out);
	if (!r.npa.empty()) {
		const bool spin = r.npa.front().has_spin_density;
		out << "\nNatural population analysis (in house):\n"
			<< "  atom       charge      core   valence   Rydberg     total" << (spin ? "  spin density\n" : "\n")
			<< fixed << setprecision(5);
		double charge_sum = 0.0, total_sum = 0.0;
		for (const NboAtomPopulation& p : r.npa) {
			out << "  " << left << setw(7) << (p.element + std::to_string(p.index)) << right << setw(10)
				<< p.charge << setw(10) << p.core << setw(10) << p.valence << setw(10) << p.rydberg
				<< setw(10) << p.total;
			if (p.has_spin_density) out << setw(14) << p.spin_density;
			out << "\n";
			charge_sum += p.charge;
			total_sum += p.total;
		}
		//charges sum to the molecular charge, populations to the electron count
		out << "  " << left << setw(7) << "total" << right << setw(10) << printable(charge_sum) << setw(40) << total_sum
			<< "\n";
	}


	//diffuse or polarisation functions leave dozens of near-empty Rydberg NBOs per atom: counted, not
	//listed; the JSON keeps every one
	constexpr double rydberg_floor = 1e-4;
	const auto hybrid = [&out](const NboHybrid& h) {
		out << "   " << left << setw(6) << (h.element + std::to_string(h.center)) << right << setprecision(2)
			<< setw(9) << h.weight_percent << setw(8) << h.s << setw(8) << h.p << setw(8) << h.d << setw(8)
			<< h.f;
		if (h.sp_exponent() > 0.0) out << setw(8) << h.sp_exponent();
		out << "\n";
	};
	for (const std::string& spin : spins_of(r.orbitals)) {
		out << title("Natural bond orbitals", spin, "in house")
			<< "    nbo  type  centres       occupancy     energy   hybrid weight %     s %     p %     d %     f %    sp^x\n";
		int hidden = 0;
		double hidden_occupancy = 0.0;
		for (const NboOrbital& o : r.orbitals) {
			if (o.spin != spin) continue;
			if (o.type.rfind("RY", 0) == 0 && o.occupancy < rydberg_floor) {
				hidden++;
				hidden_occupancy += o.occupancy;
				continue;
			}
			const std::string centres = nbo_label(o).substr(o.hybrids.empty() ? 0 : o.type.size() + 1);
			out << setw(7) << o.index << "  " << left << setw(6) << o.type << setw(12) << centres << right
				<< fixed << setprecision(5) << setw(11) << printable(o.occupancy) << setw(11) << o.energy;
			if (o.hybrids.empty()) out << "\n";
			for (size_t k = 0; k < o.hybrids.size(); k++) {
				if (k) out << string(49, ' ');
				hybrid(o.hybrids[k]);
			}
		}
		if (hidden)
			out << "  " << hidden << " further Rydberg NBOs below " << scientific << setprecision(0)
				<< rydberg_floor << " e, together " << fixed << setprecision(5) << printable(hidden_occupancy)
				<< " e\n";
	}

	if (!r.e2.empty()) {
		//"4 BD O1-H2": the number finds the orbital in the table above, the label says what it is
		std::map<std::pair<std::string, int>, const NboOrbital*> by_index;
		for (const NboOrbital& o : r.orbitals) by_index[{ o.spin, o.index }] = &o;
		const auto name = [&by_index](const std::string& spin, const int index, const std::string& fallback) {
			const auto it = by_index.find({ spin, index });
			return it == by_index.end() ? fallback : std::to_string(index) + " " + nbo_label(*it->second);
		};
		for (const std::string& spin : spins_of(r.e2)) {
			out << title("Second-order donor-acceptor energies in the NBO basis", spin, "in house") << "  "
				<< left << setw(21) << "donor" << setw(21) << "acceptor" << right << setw(14) << "E(2) kcal/mol"
				<< setw(15) << "E(j)-E(i) Eh" << setw(12) << "F(i,j) Eh" << "\n";
			for (const NboE2Entry& e : r.e2)
				if (e.spin == spin)
					out << "  " << left << setw(21) << name(spin, e.donor_index, e.donor) << setw(21)
						<< name(spin, e.acceptor_index, e.acceptor) << right << fixed << setprecision(2)
						<< setw(14) << e.energy_kcal << setprecision(3) << setw(15) << e.e_diff << setw(12)
						<< e.fij << "\n";
		}
	}
	print_nrt(r, out);
}
