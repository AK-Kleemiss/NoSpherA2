#include "pch.h"
#include "tuning.h"
#include "nbo.h"
#include "wfn_class.h"
#include "constants.h"
#include "convenience.h"

#include <Eigen/Dense>
#include <fstream>
#include <nlohmann/json.hpp>

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
			const char* e = tuning("NBO_BOND_FLOOR");
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
	//full n x n per block, row = MO; only -nbo_cube reads it, so a file without one is no error here
	const vec mo = numbers_of(section_body(text, "$LCAOMO"));
	const size_t square = static_cast<size_t>(res.n) * res.n;
	if (mo.size() >= square * blocks)
		for (int b = 0; b < blocks; b++) {
			dMatrix2 C(res.n, res.n);
			for (int k = 0; k < res.n; k++)
				for (int i = 0; i < res.n; i++) C(k, i) = mo[square * b + static_cast<size_t>(k) * res.n + i];
			res.lcaomo.push_back(C);
		}
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
			bool ok = d2 <= cut * cut;
			//No bond through an atom lying on the line between them (A-X-B above 150 degrees): the loose NRT
			//scale otherwise admits Ni(CO)4's Ni-O across the linear Ni-C-O
			for (size_t x = 0; x < na && ok; x++) {
				if (x == a || x == b) continue;
				double ax[3], bx[3];
				for (unsigned int k = 0; k < 3; k++) {
					ax[k] = atoms[a].get_coordinate(k) - atoms[x].get_coordinate(k);
					bx[k] = atoms[b].get_coordinate(k) - atoms[x].get_coordinate(k);
				}
				const double la = ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2];
				const double lb = bx[0] * bx[0] + bx[1] * bx[1] + bx[2] * bx[2];
				if (la >= d2 || lb >= d2) continue;
				const double c = (ax[0] * bx[0] + ax[1] * bx[1] + ax[2] * bx[2]) / std::sqrt(la * lb);
				ok = c > -0.8660254037844386;  //cos 150
			}
			out[a][b] = out[b][a] = ok;
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

//Every structure that moves a lone pair on a along an alternating path a-x=y-...-b onto a bond to b, which has
//a vacancy: LP(a) -> a=x, x=y -> x-y, ..., the last atom gains a bond to b. t is valence-only (no cores on the
//diagonal), sigma is the bond graph, at most max_pi pi bonds are shifted.
static void alternating_shifts(ivec2& t, const ivec2& sigma, const bvec& vacant, const int cur, const int max_pi,
							   bvec& seen, std::vector<ivec2>& out)
{
	const int n = static_cast<int>(t.size());
	for (int x = 0; x < n; x++) {
		if (seen[x] || sigma[cur][x] < 1) continue;
		if (vacant[x]) {
			ivec2 u = t;
			u[cur][x]++;
			u[x][cur]++;
			out.push_back(u);
		}
		if (max_pi == 0) continue;
		for (int y = 0; y < n; y++) {
			if (seen[y] || y == x || t[x][y] < 2) continue;
			t[cur][x]++; t[x][cur]++; t[x][y]--; t[y][x]--;
			seen[x] = seen[y] = true;
			alternating_shifts(t, sigma, vacant, y, max_pi - 1, seen, out);
			seen[x] = seen[y] = false;
			t[cur][x]--; t[x][cur]--; t[x][y]++; t[y][x]++;
		}
	}
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
	//d-block metals: their cap is -nbo_tm_pairs and the first relax pass keeps it. Only with -nbo_tm_pairs: holding
	//the default (valence NAO) metal cap while the ligands' caps lift booked Ni(CO)4's d pairs as three Ni-C "bonds"
	//97 % on Ni (NBO 7: one dative Ni-C and five d lone pairs)
	std::vector<char> tm(natoms, 0);
	for (int a = 0; a < natoms; a++) {
		const int z = nao.atoms[a].Z;
		tm[a] = (z >= 21 && z <= 30) || (z >= 39 && z <= 48) || z == 57 || (z >= 72 && z <= 80) || z == 89 || (z >= 104 && z <= 112);
		if (tm[a] && options.tm_pairs > 0) cap[a] = options.tm_pairs;
	}
	const bool any_tm = options.tm_pairs > 0 && std::any_of(tm.begin(), tm.end(), [](const char c) { return c != 0; });
	//$CHOOSE: the ladder only looks at blocks the chosen structure still lacks, so fed its own structure
	//it accepts the same orbitals in the same order
	ivec2 want = options.parent_lewis;
	std::vector<ivec> want_multi = options.parent_multi;
	const bool chosen = !want.empty();
	//-nbo_centres: connected sets of 3..max_centres atoms over the bondable graph (a chain or a ring through each
	//member), searched only where one- and two-centre blocks fail; under $CHOOSE only the chosen sets
	std::vector<ivec> msets;
	{
		std::set<ivec> sets;
		if (chosen) sets.insert(want_multi.begin(), want_multi.end());
		else if (options.max_centres > 2) {
			const auto bonded = [&](const int a, const int b) { return a != b && (bondable.empty() || bondable[a][b]); };
			for (int b = 0; b < natoms; b++)
				for (int a = 0; a < natoms; a++)
					for (int c = a + 1; c < natoms; c++)
						if (bonded(a, b) && bonded(b, c) && !idx[a].valence.empty() && !idx[b].valence.empty() &&
							!idx[c].valence.empty()) {
							ivec s{ a, b, c };
							std::sort(s.begin(), s.end());
							sets.insert(s);
						}
			for (int k = 4; k <= std::min(options.max_centres, 4); k++) {
				std::set<ivec> grown;
				for (const ivec& s : sets) {
					if (static_cast<int>(s.size()) != k - 1) continue;
					for (int x = 0; x < natoms; x++) {
						if (idx[x].valence.empty() || std::count(s.begin(), s.end(), x)) continue;
						if (std::none_of(s.begin(), s.end(), [&](const int a) { return bonded(a, x); })) continue;
						ivec t = s;
						t.push_back(x);
						std::sort(t.begin(), t.end());
						grown.insert(t);
					}
				}
				sets.insert(grown.begin(), grown.end());
			}
			for (const ivec& s : options.skip_multi) sets.erase(s);
			//-nrt_atoms: every centre inside the subspace, else ligand ring pi reads as allyl-type 3c units that do
			//pay their 0.05 e (18 of them, all C-C-C / N-C-C, on NiLi_Empty)
			const ivec& sub = options.nrt_subspace;
			if (!sub.empty())
				std::erase_if(sets, [&](const ivec& s) {
					return std::any_of(s.begin(), s.end(), [&](const int a) { return std::find(sub.begin(), sub.end(), a + 1) == sub.end(); });
				});
		}
		msets.assign(sets.begin(), sets.end());
	}
	if (chosen) {
		int total = static_cast<int>(res.orbitals.size() + want_multi.size());
		for (int a = 0; a < natoms && a < static_cast<int>(want.size()); a++)
			for (int b = a; b < natoms; b++) total += want[a][b];
		err_checkf(static_cast<int>(want.size()) == natoms && total == n_pairs,
				   "NBO: the chosen Lewis structure fills " + std::to_string(total) + " orbitals over " +
				   std::to_string(want.size()) + " atoms, this density " + std::to_string(n_pairs) +
				   " over " + std::to_string(natoms), log);
	}
	//relax 0: every cap; 1: main-group caps lifted; 2: metal caps too, logged
	int relax = 0;
	const auto full = [&](const int a) { return used[a] >= cap[a] && (relax == 0 || (relax == 1 && any_tm && tm[a])); };
	for (; relax < 3; relax++) {
	  if (relax && static_cast<int>(res.orbitals.size()) >= n_pairs) break;
	  if (relax == 2 && !any_tm) break;
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
				if (full(a)) continue;
				if (chosen && want[a][a] == 0) continue;
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
				if (chosen) want[ba][ba]--;
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
				if (full(a)) continue;
				for (int b = a + 1; b < natoms; b++) {
					if (idx[b].valence.empty()) continue;
					if (full(b)) continue;
					if (chosen ? want[a][b] == 0 : !bondable.empty() && !bondable[a][b]) continue;
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
				if (chosen) {
					want[pa][pb]--;
					want[pb][pa]--;
				}
				used[pa]++;
				used[pb]++;
				res.threshold = t;
				progress = true;
				continue;
			}
			//k-centre blocks: a superset of a failed pair always scores higher, so every centre must hold
			//10 % of the orbital or a bond with a tail would pass as a three-centre bond
			std::vector<const ivec*> mlist;
			for (const ivec& s : msets) {
				if (std::any_of(s.begin(), s.end(), full)) continue;
				if (chosen && std::find(want_multi.begin(), want_multi.end(), s) == want_multi.end()) continue;
				mlist.push_back(&s);
			}
			const int nm = static_cast<int>(mlist.size());
			vec lamm(nm, 0.0);
			std::vector<VectorXd> candm(nm);
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
			for (int c = 0; c < nm; c++) {
				ivec block;
				for (const int a : *mlist[c]) block.insert(block.end(), idx[a].all.begin(), idx[a].all.end());
				lamm[c] = leading_block(G, block, candm[c]);
			}
			double bestm = t;
			int pm = -1;
			for (int c = 0; c < nm; c++) {
				if (lamm[c] <= bestm) continue;
				if (std::any_of(mlist[c]->begin(), mlist[c]->end(), [&](const int a) { return weight_on(candm[c], idx[a].all) < 0.10; })) continue;
				bestm = lamm[c];
				pm = c;
			}
			if (pm >= 0) {
				const ivec& s = *mlist[pm];
				accept(candm[pm], std::to_string(s.size()) + "C", s);
				res.multi.push_back(s);
				if (chosen) want_multi.erase(std::find(want_multi.begin(), want_multi.end(), s));
				for (const int a : s) used[a]++;
				res.threshold = t;
				progress = true;
			}
		}
		if (static_cast<int>(res.orbitals.size()) >= n_pairs) break;
	  }
	}
	for (int a = 0; a < natoms; a++)
		if (tm[a] && used[a] > cap[a])
			log << "NBO: " << constants::atnr2letter(nao.atoms[a].Z) << a + 1 << " holds " << used[a] << " Lewis pairs, above its "
				<< cap[a] << ": no other placement filled the " << n_pairs << " pairs. Its bond orders are on a wider valence "
				<< "space than the other metals'; compare -nbo_tm_pairs runs" << std::endl;
	//What thresholds, the minority floor and the caps left of a chosen structure (a dative bond, a lone
	//pair in a Rydberg-heavy block), best block first; the count check above ends the loop.
	while (chosen && static_cast<int>(res.orbitals.size()) < n_pairs) {
		double best = -1.0;
		VectorXd bv;
		int ba = -1, bb = -1;
		for (int a = 0; a < natoms; a++)
			for (int b = a; b < natoms; b++) {
				if (want[a][b] == 0) continue;
				ivec block = idx[a].all;
				if (b != a) block.insert(block.end(), idx[b].all.begin(), idx[b].all.end());
				VectorXd v;
				const double lam = leading_block(G, block, v);
				if (lam > best) { best = lam; bv = v; ba = a; bb = b; }
			}
		int bm = -1;
		for (int c = 0; c < static_cast<int>(want_multi.size()); c++) {
			ivec block;
			for (const int a : want_multi[c]) block.insert(block.end(), idx[a].all.begin(), idx[a].all.end());
			VectorXd v;
			const double lam = leading_block(G, block, v);
			if (lam > best) { best = lam; bv = v; bm = c; }
		}
		if (bm >= 0) {
			accept(bv, std::to_string(want_multi[bm].size()) + "C", want_multi[bm]);
			res.multi.push_back(want_multi[bm]);
			want_multi.erase(want_multi.begin() + bm);
			continue;
		}
		accept(bv, ba == bb ? "LP" : "BD", ba == bb ? ivec{ ba } : ivec{ ba, bb });
		res.topo[ba][bb]++;
		want[ba][bb]--;
		if (ba != bb) {
			res.topo[bb][ba]++;
			want[bb][ba]--;
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
	bvec2 allowed;
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
		//An antibond holding over half an orbital means bond and antibond are two lone pairs rotated
		//into each other (a four-electron pair), which the ladder accepted from the leading eigenvalue.
		//A chosen parent ($CHOOSE) places its bonds regardless of bondable, so banning one there would recurse forever.
		if (options.parent_lewis.empty() && w.dot(G0 * w) > 0.5 * scale * nrm * nrm) {
			if (allowed.empty()) allowed = bondable.empty() ? bvec2(natoms, bvec(natoms, true)) : bondable;
			allowed[a][b] = allowed[b][a] = false;
			if (options.debug)
				log << "NBO search: " << nao.atoms[a].label << "-" << nao.atoms[b].label
					<< " is no bond after refinement, searching again without it" << std::endl;
		}
		NboFunction f;
		f.type = "BD*";
		f.centers = { a, b };
		f.multiplicity = res.orbitals[j].multiplicity;
		non_lewis.push_back(w / nrm);
		non_lewis_fn.push_back(f);
	}
	//every pass removes at least one pair, so this ends
	if (!allowed.empty()) return nbo_search(nao, gamma, allowed, n_pairs, scale, options, log);
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
	//A k-centre bond sum_i c_i h_i leaves k - 1 combinations of its own hybrids. Without what is already placed
	//(Lewis, antibonds, earlier units: two units on one set share hybrids), G0 diagonalised in them: the least
	//occupied is kC*, the rest kCn.  ponytail: sequential, not symmetric like the antibonds; order = Lewis order
	for (int j = 0; j < res.n_lewis; j++) {
		const ivec& cs = res.orbitals[j].centers;
		const int k = static_cast<int>(cs.size());
		if (k < 3) continue;
		MatrixXd H = MatrixXd::Zero(n, k);
		VectorXd c(k);
		for (int i = 0; i < k; i++) {
			for (const int x : idx[cs[i]].all) H(x, i) = vectors[j](x);
			c(i) = H.col(i).norm();
			if (c(i) > 1e-12) H.col(i) /= c(i);
		}
		Eigen::SelfAdjointEigenSolver<MatrixXd> pc(MatrixXd(MatrixXd::Identity(k, k) - c * c.transpose() / c.squaredNorm()));
		MatrixXd A = H * pc.eigenvectors().rightCols(k - 1);
		MatrixXd P(n, res.n_lewis + static_cast<int>(non_lewis.size()));
		for (int q = 0; q < res.n_lewis; q++) P.col(q) = vectors[q];
		for (size_t q = 0; q < non_lewis.size(); q++) P.col(res.n_lewis + static_cast<int>(q)) = non_lewis[q];
		A -= P * (P.transpose() * A);
		Eigen::SelfAdjointEigenSolver<MatrixXd> sa(MatrixXd(A.transpose() * A));
		ivec keep;
		for (int q = 0; q < k - 1; q++)
			if (sa.eigenvalues()(q) > 1e-6) keep.push_back(q);
		if (keep.empty()) continue;
		MatrixXd B(n, static_cast<int>(keep.size()));
		for (size_t q = 0; q < keep.size(); q++)
			B.col(static_cast<int>(q)) = A * sa.eigenvectors().col(keep[q]) / std::sqrt(sa.eigenvalues()(keep[q]));
		Eigen::SelfAdjointEigenSolver<MatrixXd> eo(MatrixXd(B.transpose() * G0 * B));
		for (int q = static_cast<int>(keep.size()) - 1; q >= 0; q--) {
			NboFunction f;
			f.type = std::to_string(k) + (q ? "Cn" : "C*");
			f.centers = cs;
			non_lewis.push_back(B * eo.eigenvectors().col(q));
			non_lewis_fn.push_back(f);
		}
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
	//The ladder is greedy: in a six-ring it can take two pi bonds that leave the third pi pair on two atoms that
	//are not neighbours, a half-filled LP on one and a half-filled LV on the other (a mis-paired Kekule structure;
	//E2 between them then moves more than one electron). Shift the pair along an alternating path onto a bond to
	//the vacancy under $CHOOSE and keep the lowest rho(NL), the search's own criterion. Each kept step lowers
	//rho(NL), so the loop ends.
	while (!chosen) {
		bvec half_lp(natoms, false), half_lv(natoms, false);
		for (const NboFunction& f : res.orbitals) {
			if (f.type == "LP" && f.occupancy < 0.6 * scale) half_lp[f.centers[0]] = true;
			if (f.type == "LV" && f.occupancy > 0.3 * scale) half_lv[f.centers[0]] = true;
		}
		ivec2 base = res.topo;
		for (int a = 0; a < natoms; a++) base[a][a] -= static_cast<int>(idx[a].core.size());
		std::vector<ivec2> tries;
		for (int a = 0; a < natoms; a++) {
			if (!half_lp[a]) continue;
			base[a][a]--;
			bvec seen(natoms, false);
			seen[a] = true;
			alternating_shifts(base, res.topo, half_lv, a, 3, seen, tries);
			base[a][a]++;
		}
		std::sort(tries.begin(), tries.end());
		tries.erase(std::unique(tries.begin(), tries.end()), tries.end());
		NboLewis best;
		best.rho_nl = res.rho_nl - 1e-4;
		bool found = false;
		for (const ivec2& t : tries) {
			NboOptions o = options;
			o.parent_lewis = t;
			o.parent_multi = res.multi;
			std::ostringstream quiet;
			NboLewis r = nbo_search(nao, gamma, bondable, n_pairs, scale, o, quiet);
			if (r.rho_nl < best.rho_nl) { best = std::move(r); found = true; }
		}
		if (!found) break;
		log << "NBO search: a half-filled lone pair and vacancy re-paired along an alternating path, rho(NL) "
			<< std::fixed << std::setprecision(5) << res.rho_nl << " -> " << best.rho_nl << " e ("
			<< tries.size() << " shifts tried)" << std::defaultfloat << std::endl;
		res = std::move(best);
	}
	//-nbo_centres: every k-centre unit pays for itself, the search without it (leave-one-out) loses 0.05 e of
	//rho(NL), else a delocalised 2c bond (BD + E2) or a slice of a ring pi system reads as a fake 3c unit.  A unit
	//the re-search swaps in for a dropped one is tested in turn.  The inner searches carry skip_multi, so they do
	//not gate again.  The 10 % per centre also holds for the final orbital: the ladder checks it on the block
	//eigenvector, and orthogonalisation turned Ni-P-P units into a P->Ni bond with 0.4 % on the second P (NiZn,
	//PdZn, PdLi MeCN).  ponytail: one search per unit tested, greedy in list order
	if (!chosen && options.max_centres > 2 && options.skip_multi.empty() && !res.multi.empty()) {
		const double need = 0.05 * (scale / 2.0);
		const double rho_all = res.rho_nl;
		NboOptions o = options;
		std::set<ivec> tested;
		//per centre the best of the unit's orbitals: a 3c-4e nonbonding orbital has a node on the middle atom
		const auto thin = [&](const ivec& u) {
			vec h(u.size(), 0.0);
			for (int j = 0; j < res.n_lewis; j++) {
				const NboFunction& f = res.orbitals[j];
				if (f.centers != u) continue;
				for (size_t t = 0; t < h.size() && t < f.center_weight.size(); t++) h[t] = std::max(h[t], f.center_weight[t]);
			}
			return *std::min_element(h.begin(), h.end()) < 0.10;
		};
		for (bool again = true; again;) {
			again = false;
			for (const ivec& u : std::vector<ivec>(res.multi)) {
				if (!tested.insert(u).second) continue;
				o.skip_multi.push_back(u);
				std::ostringstream quiet;
				NboLewis r = nbo_search(nao, gamma, bondable, n_pairs, scale, o, quiet);
				if (thin(u) || r.rho_nl - res.rho_nl < need) {
					res = std::move(r);
					again = true;
					break;
				}
				o.skip_multi.pop_back();
			}
		}
		log << "NBO search: " << tested.size() << " k-centre unit(s) tested, " << res.multi.size()
			<< " kept (each costs >= " << need << " e of rho(NL) when left out, 10 % on every centre), rho(NL) " << std::fixed
			<< std::setprecision(5) << rho_all << " -> " << res.rho_nl << " e" << std::defaultfloat << std::endl;
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

//which NRT structures an E2 average used, for the log and the printout
static std::string e2_avg_window(const NboE2Average& a)
{
	std::ostringstream s;
	if (a.min_percent > 0.0) s << "of weight >= " << a.min_percent << " %";
	else s << "within " << a.tol_percent << " % of the leading weight";
	return s.str();
}

namespace
{
	//-nrt_e2_avg: one Lewis structure's E2 says as much about the structure as about the molecule. Every NRT
	//structure above the cut is placed under $CHOOSE, its E2 summed per channel (type + centres), and the
	//channels averaged with the NRT weights renormalised over the set; a channel a structure lacks counts 0.
	void e2_nrt_average(NboResults& res, const NAOResult& nao, const dMatrix2& gamma, const dMatrix2& fock_nao,
						const bvec2& bondable, const int n_pairs, const double scale, const NboOptions& options,
						std::ostream& log)
	{
		if (scale != 2.0 || !fock_nao.extent(0)) {
			log << "-nrt_e2_avg: closed shell with a Fock matrix only, skipped" << std::endl;
			return;
		}
		const auto t0 = std::chrono::steady_clock::now();
		NboE2Average& avg = res.e2_avg;
		avg.tol_percent = options.nrt_e2_avg_percent;
		avg.min_percent = options.nrt_e2_avg_min_percent;
		std::vector<const NboNrtCandidate*> picked;
		double wsum = 0.0, wmax = 0.0;
		for (const auto& w : res.nrt.weights) wmax = std::max(wmax, w.weight_percent);
		const double floor = avg.min_percent > 0.0 ? avg.min_percent : wmax * (1.0 - options.nrt_e2_avg_percent / 100.0);
		for (size_t i = 0; i < res.nrt.weights.size() && i < res.nrt.candidates.size(); i++)
			if (res.nrt.weights[i].weight_percent >= floor) {
				picked.push_back(&res.nrt.candidates[i]);
				avg.structures.push_back(res.nrt.weights[i].structure);
				avg.weights.push_back(res.nrt.weights[i].weight_fraction);
				wsum += res.nrt.weights[i].weight_fraction;
			}
		if (picked.empty()) {
			log << "-nrt_e2_avg: no NRT structures" << std::endl;
			return;
		}
		for (double& w : avg.weights) w /= wsum;
		const auto atom = [&](const int c) { return constants::atnr2letter(nao.atoms[c].Z) + std::to_string(c + 1); };
		const auto key = [&](const NboFunction& f) {
			std::string s = f.type + " ";
			for (size_t k = 0; k < f.centers.size(); k++) s += (k ? "-" : "") + atom(f.centers[k]);
			return s;
		};
		//One structure per thread, each search single-threaded: a search alone is mostly serial (ladder picks,
		//Rydberg pivots, Eigen eigensolves), so the structures, which share nothing, are the parallel axis.
		const int ns = static_cast<int>(picked.size());
		avg.rho_nl.assign(ns, 0.0);
		//per entry {kcal, electrons}: second order moves q_i (F_ij/dE)^2 = E2/dE into the acceptor, which diverges as
		//dE -> 0 (a 2 kcal/mol pair read as 2000 e); the 2x2 mixing q_i sin^2(theta) = q_i (1 - |dE|/sqrt(dE^2 + 4F^2))/2
		//has the same small-F limit and stops at half the donor
		const auto moved = [](const NboE2Entry& e) {
			if (e.fij == 0.0) return 0.0;
			const double occ = e.energy_kcal / constants::kcal_mol_per_hartree * std::abs(e.e_diff) / (e.fij * e.fij);
			return 0.5 * occ * (1.0 - std::abs(e.e_diff) / std::sqrt(e.e_diff * e.e_diff + 4.0 * e.fij * e.fij));
		};
		std::vector<std::vector<std::pair<std::pair<std::string, std::string>, std::pair<double, double>>>> found(ns);
#ifdef _OPENMP
		const int nthreads = std::min(ns, options.threads > 0 ? options.threads : omp_get_max_threads());
#endif
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
		for (int s = 0; s < ns; s++) {
			NboOptions chosen = options;
			chosen.parent_lewis = picked[s]->topo;
			chosen.parent_multi = picked[s]->multi;
			chosen.search_threads = 1;
			std::ostringstream quiet;
			NboLewis lewis = nbo_search(nao, gamma, bondable, n_pairs, scale, chosen, quiet);
			avg.rho_nl[s] = lewis.rho_nl;
			for (const NboE2Entry& e : nbo_e2(lewis, fock_nao, options.e2_threshold_kcal))
				found[s].push_back({ { key(lewis.orbitals[e.donor_index - 1]), key(lewis.orbitals[e.acceptor_index - 1]) },
									 { e.energy_kcal, moved(e) } });
		}
		//per channel, the summed E2 and electrons of each structure
		std::map<std::pair<std::string, std::string>, std::pair<vec, vec>> sums;
		for (int s = 0; s < ns; s++)
			for (const auto& [k, e] : found[s]) {
				auto& [v, q] = sums[k];
				v.resize(ns, 0.0);
				q.resize(ns, 0.0);
				v[s] += e.first;
				q[s] += e.second;
			}
		for (const auto& [k, vq] : sums) {
			const vec& v = vq.first;
			const vec& q = vq.second;
			NboE2Channel c;
			c.donor = k.first;
			c.acceptor = k.second;
			c.min_kcal = *std::min_element(v.begin(), v.end());
			c.max_kcal = *std::max_element(v.begin(), v.end());
			c.min_e = *std::min_element(q.begin(), q.end());
			c.max_e = *std::max_element(q.begin(), q.end());
			for (size_t s = 0; s < v.size(); s++) {
				c.mean_kcal += avg.weights[s] * v[s];
				c.mean_e += avg.weights[s] * q[s];
				if (v[s] > 0.0) c.present += avg.weights[s];
			}
			avg.channels.push_back(c);
		}
		std::stable_sort(avg.channels.begin(), avg.channels.end(),
						 [](const NboE2Channel& a, const NboE2Channel& b) { return a.mean_kcal > b.mean_kcal; });
		avg.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		log << "-nrt_e2_avg: " << picked.size() << " NRT structure(s) " << e2_avg_window(avg)
			<< ", covering " << std::fixed << std::setprecision(1) << 100.0 * wsum << " % of the weight, "
			<< avg.channels.size() << " channels, " << std::setprecision(2) << avg.seconds << " s"
			<< std::defaultfloat << std::endl;
	}

	//One spin's worth of the analysis, appended to the results.
	void analyse_spin(NboResults& res, const NboInput& in, const NAOResult& nao,
					  const bvec2& bondable, const bvec2& nrt_bondable, const vec2& xyz, const dMatrix2& density,
					  const dMatrix2& fock,
					  const int n_pairs, const double scale, const std::string& spin,
					  const NboOptions& options, std::vector<NboLewis>& out_lewis, std::ostream& log)
	{
		const dMatrix2 gamma = nao_density(density, in.overlap, nao.C);
		const dMatrix2 fock_nao = fock.extent(0) ? nao_operator(fock, nao.C) : dMatrix2();
		const auto clock = [] { return std::chrono::steady_clock::now(); };
		auto t = clock();
		NboLewis lewis = nbo_search(nao, gamma, bondable, n_pairs, scale, options, log);
		if (!options.parent_lewis.empty()) {
			//the free search once more, so the log can say what the chosen structure costs
			NboOptions free = options;
			free.parent_lewis.clear();
			free.parent_multi.clear();
			std::ostringstream quiet;
			const NboLewis own = nbo_search(nao, gamma, bondable, n_pairs, scale, free, quiet);
			std::string diff;
			int n_diff = 0;
			for (size_t a = 0; a < own.topo.size(); a++)
				for (size_t b = a; b < own.topo.size(); b++) {
					if (own.topo[a][b] == lewis.topo[a][b]) continue;
					const std::string la = constants::atnr2letter(nao.atoms[a].Z) + std::to_string(a + 1);
					//valence lone pairs, as the file gives them: the diagonal also holds the cores
					const int core = a != b ? 0 : static_cast<int>(std::count_if(nao.orbitals.begin(), nao.orbitals.end(),
						[&](const NAO& o) { return o.atom == static_cast<int>(a) && o.type == NAOClass::Core; }));
					if (n_diff++ < 8)
						diff += " " + (a == b ? "LP(" + la + ")" : la + "-" + constants::atnr2letter(nao.atoms[b].Z) +
								std::to_string(b + 1)) + " " + std::to_string(own.topo[a][b] - core) + "->" +
								std::to_string(lewis.topo[a][b] - core);
				}
			log << "NRT parent from " << (options.nrt_parent.empty() ? std::string("the caller") : options.nrt_parent.string())
				<< ": rho(NL) = " << std::fixed << std::setprecision(5) << lewis.rho_nl << " e against "
				<< own.rho_nl << " e for the search's own best structure, "
				<< (n_diff ? std::to_string(n_diff) + " entries differ (own->chosen):" + diff + (n_diff > 8 ? " ..." : "")
						   : std::string("the same structure"))
				<< std::defaultfloat << std::endl;
		}
		res.search_seconds += std::chrono::duration<double>(clock() - t).count();
		t = clock();
		std::vector<NboE2Entry> e2 = nbo_e2(lewis, fock_nao, options.e2_threshold_kcal);
		res.e2_seconds += std::chrono::duration<double>(clock() - t).count();
		//before the renumbering below, while e2's indices still point into lewis.orbitals
		if (options.nrt)
			//log, not std::cout: under -fba this runs on the NBO thread and cout is RGBI's
			native_nrt(res.nrt, nao, lewis, e2, nrt_bondable, options, spin, scale, log, &xyz);
		if (options.nrt && (options.nrt_e2_avg_percent > 0.0 || options.nrt_e2_avg_min_percent > 0.0))
			e2_nrt_average(res, nao, gamma, fock_nao, bondable, n_pairs, scale, options, log);

		//NBO numbers by type group: the Lewis set in found order, then LV, BD*, RY
		ivec order;
		static const char* groups[] = { "CR", "LP", "BD", "3C", "4C", "LV", "BD*", "3Cn", "4Cn", "3C*", "4C*", "RY" };
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
		lewis.printed = order;
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

std::string read_lewis_structure(const std::filesystem::path& file, const NAOResult& nao, const int n_pairs,
								 const std::string& wavefunction, ivec2& topo, std::ostream& log, std::vector<ivec>* multi)
{
	using nlohmann::json;
	std::ifstream in(file);
	if (!in.good()) return "cannot open " + file.string();
	json root;
	try {
		root = json::parse(in);
	}
	catch (const json::exception& e) {
		return std::string("not valid JSON: ") + e.what();
	}
	//a whole .native.nbo.json carries the structure under "parent_lewis"
	const json& j = root.is_object() && root.contains("parent_lewis") ? root["parent_lewis"] : root;
	if (!j.is_object()) return "expected a JSON object with \"lone_pairs\" and \"bonds\"";
	for (const auto& [key, value] : j.items()) {
		static const std::set<std::string> known = { "format", "version", "wavefunction", "atoms", "lone_pairs",
													 "bonds", "three_centre", "multi_centre", "alpha", "beta" };
		if (!known.count(key)) return "unknown key \"" + key + "\"";
	}
	if (j.contains("format") && j["format"] != "NoSpherA2 Lewis structure")
		return "\"format\" is " + j["format"].dump() + ", not \"NoSpherA2 Lewis structure\"";
	if (j.contains("version") && (!j["version"].is_number_integer() || j["version"].get<int>() != 1))
		return "\"version\" " + j["version"].dump() + " is not 1, the only version this build reads";
	if (j.contains("alpha") || j.contains("beta"))
		return "per-spin \"alpha\"/\"beta\" structures are not supported: -nrt_parent takes a closed-shell wavefunction";
	for (const char* k : { "three_centre", "multi_centre" })
		if (!multi && j.contains(k) && !(j[k].is_array() && j[k].empty()))
			return std::string("\"") + k + "\" bonds are not supported here: the caller places one- and two-centre orbitals only";
	const int na = static_cast<int>(nao.atoms.size());
	if (j.contains("atoms") && !(j["atoms"].is_number_integer() && j["atoms"].get<int>() == na))
		return "\"atoms\" is " + j["atoms"].dump() + ", the wavefunction has " + std::to_string(na);
	if (j.contains("wavefunction") && j["wavefunction"].is_string() &&
		j["wavefunction"].get<std::string>() != wavefunction)
		log << "Warning: -nrt_parent " << file.string() << " was written for " << j["wavefunction"].get<std::string>()
			<< ", this is " << wavefunction << "; the atom count and elements match, so it is used" << std::endl;
	//"Pd1", or a bare 1-based index as a number or a string
	const auto atom_of = [&](const json& x, int& a) -> std::string {
		std::string s = x.is_string() ? x.get<std::string>() : x.is_number_integer() ? std::to_string(x.get<int>()) : "";
		if (s.empty()) return "atom " + x.dump() + " is neither a label like \"C1\" nor a 1-based index";
		const size_t d = s.find_first_of("0123456789");
		if (d == std::string::npos || s.find_first_not_of("0123456789", d) != std::string::npos || d > 2)
			return "atom " + x.dump() + " is neither a label like \"C1\" nor a 1-based index";
		a = std::stoi(s.substr(d)) - 1;
		if (a < 0 || a >= na) return "atom " + x.dump() + " is outside 1.." + std::to_string(na);
		const std::string el = constants::atnr2letter(nao.atoms[a].Z), given = s.substr(0, d);
		const auto same = [](char p, char q) { return std::tolower(p) == std::tolower(q); };
		if (d > 0 && (el.size() != d || !std::equal(el.begin(), el.end(), given.begin(), same)))
			return "atom " + x.dump() + ": atom " + std::to_string(a + 1) + " is " + el + ", not " + given;
		return "";
	};
	topo.assign(na, ivec(na, 0));
	if (!j.contains("lone_pairs") && !j.contains("bonds")) return "holds neither \"lone_pairs\" nor \"bonds\"";
	if (j.contains("lone_pairs")) {
		if (!j["lone_pairs"].is_object()) return "\"lone_pairs\" must be an object like {\"O1\": 2}";
		bvec seen(na, false);
		for (const auto& [key, value] : j["lone_pairs"].items()) {
			int a = -1;
			const std::string why = atom_of(json(key), a);
			if (!why.empty()) return "lone_pairs: " + why;
			if (seen[a]) return "lone_pairs: atom " + std::to_string(a + 1) + " is listed twice";
			seen[a] = true;
			if (!value.is_number_integer() || value.get<int>() < 0)
				return "lone_pairs: \"" + key + "\": " + value.dump() + " is not a lone-pair count >= 0";
			topo[a][a] = value.get<int>();
		}
	}
	if (j.contains("bonds")) {
		if (!j["bonds"].is_array()) return "\"bonds\" must be an array of [A, B, order]";
		for (const json& b : j["bonds"]) {
			if (!b.is_array() || b.size() != 3) return "bond " + b.dump() + " is not [A, B, order]";
			int x = -1, y = -1;
			std::string why = atom_of(b[0], x);
			if (why.empty()) why = atom_of(b[1], y);
			if (!why.empty()) return "bond " + b.dump() + ": " + why;
			if (x == y) return "bond " + b.dump() + " joins an atom to itself";
			if (!b[2].is_number_integer() || b[2].get<int>() < 1 || b[2].get<int>() > 3)
				return "bond " + b.dump() + ": the order must be 1, 2 or 3";
			if (topo[x][y]) return "bond " + b.dump() + " is listed twice";
			topo[x][y] = topo[y][x] = b[2].get<int>();
		}
	}
	//one electron pair per unit over 3 or 4 atoms; "three_centre" is the older name of the same list
	if (multi) multi->clear();
	for (const char* k : { "three_centre", "multi_centre" }) {
		if (!multi || !j.contains(k)) continue;
		if (!j[k].is_array()) return std::string("\"") + k + "\" must be an array of [A, B, C(, D)]";
		for (const json& u : j[k]) {
			if (!u.is_array() || u.size() < 3 || u.size() > 4) return "unit " + u.dump() + " is not [A, B, C(, D)]";
			ivec s;
			for (const json& x : u) {
				int a = -1;
				const std::string why = atom_of(x, a);
				if (!why.empty()) return "unit " + u.dump() + ": " + why;
				s.push_back(a);
			}
			std::sort(s.begin(), s.end());
			if (std::adjacent_find(s.begin(), s.end()) != s.end()) return "unit " + u.dump() + " names an atom twice";
			multi->push_back(s);
		}
	}
	int lps = 0, bonds = multi ? static_cast<int>(multi->size()) : 0, cores = 0;
	for (int a = 0; a < na; a++) {
		lps += topo[a][a];
		for (int b = a + 1; b < na; b++) bonds += topo[a][b];
	}
	for (const NAO& o : nao.orbitals) cores += o.type == NAOClass::Core;
	const int have = 2 * (bonds + lps + cores), need = 2 * n_pairs;
	if (have != need)
		return std::to_string(bonds) + " bond pairs + " + std::to_string(lps) + " lone pairs + " + std::to_string(cores) +
			   " core pairs hold " + std::to_string(have) + " electrons, the wavefunction has " + std::to_string(need) +
			   ": " + std::to_string(std::abs(have - need) / 2) + " pair(s) too " + (have > need ? "many" : "few");
	return "";
}

void refuse_unsupported_nbo_source(const WFN& wavy, std::ostream& log)
{
	err_checkf(wavy.get_nr_basis_set_loaded() == wavy.get_ncen(),
		"NBO needs the contracted basis set; a primitive-only source (.wfn/.wfx) does not carry one."
		" Use the .gbw, .fchk or .molden of the same calculation.", log);
	err_checkf(!wavy.get_d_f_switch(),
		"NBO needs a spherical basis; this wavefunction carries Cartesian d/f shells. "
		"Convert through a spherical .gbw/.molden of the same calculation.", log);
}

NboResults native_nbo(WFN& wavy, const NboOptions& options, std::ostream& log, NboOrbitalSpace* space)
{
	const auto clock = [] { return std::chrono::steady_clock::now(); };
	const auto t_f47 = clock();
	std::filesystem::path f47 = options.file47;
	bool temporary = false;
	if (f47.empty()) {
		f47 = std::filesystem::temp_directory_path() /
			  ("nbo_native_" + std::to_string(static_cast<unsigned long long>(
				   std::chrono::steady_clock::now().time_since_epoch().count())) + ".47");
		refuse_unsupported_nbo_source(wavy, log);
		if (!wavy.write_nbo(f47, options.debug, options.debug ? &log : nullptr, "")) {
			std::error_code ec;
			std::filesystem::remove(f47, ec);
			err_checkf(false, "Could not write the FILE47 the native NBO analysis reads", log);
		}
		temporary = !options.keep_file47;
	}
	const NboInput in = read_file47(f47);
	if (temporary) {
		//Windows refuses to delete a file another process holds open (the CI runner's virus scanner on the
		//freshly written .47); a stale temp file must not throw away an analysis that has its input already
		std::error_code ec;
		std::filesystem::remove(f47, ec);
	}

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
	vec2 xyz(atoms.size(), vec(3, 0.0));
	for (size_t a = 0; a < atoms.size(); a++)
		for (int x = 0; x < 3; x++) xyz[a][x] = atoms[a].get_coordinate(x);
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
		NboOptions chosen = options;
		if (!options.nrt_parent.empty()) {
			const std::string why = read_lewis_structure(options.nrt_parent, nao, n_pairs,
														 wavy.get_path().filename().string(), chosen.parent_lewis, log,
														 &chosen.parent_multi);
			err_checkf(why.empty(), "-nrt_parent " + options.nrt_parent.string() + ": " + why, log);
			//NRT rebuilds the parent's orbitals over its own bondable pairs and cannot place one outside them
			for (size_t a = 0; a < atoms.size() && options.nrt; a++)
				for (size_t b = a + 1; b < atoms.size(); b++)
					err_checkf(!chosen.parent_lewis[a][b] || nrt_bondable[a][b],
							   "-nrt_parent " + options.nrt_parent.string() + ": bond " +
							   constants::atnr2letter(nao.atoms[a].Z) + std::to_string(a + 1) + "-" +
							   constants::atnr2letter(nao.atoms[b].Z) + std::to_string(b + 1) +
							   " is longer than -nrt_bond_scale " + std::to_string(options.nrt_bond_scale) +
							   " x the covalent radius sum; raise it or drop the bond", log);
		}
		analyse_spin(res, in, nao, bondable, nrt_bondable, xyz, in.density[0],
					 in.fock.empty() ? dMatrix2() : in.fock[0], n_pairs, 2.0, "", chosen, lewis, log);
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
		if (space) {
			space->nao = nao;
			space->overlap = in.overlap;
			space->lcaomo = in.lcaomo.empty() ? dMatrix2() : in.lcaomo[0];
			space->lewis = lewis[0];
		}
	}
	else {
		err_checkf(!space, "-nbo_cube/-nlmo_cube: closed-shell wavefunctions only; this one is open shell", log);
		err_checkf(options.nrt_parent.empty() && options.parent_lewis.empty(),
				   "-nrt_parent takes a closed-shell wavefunction; this one is open shell, whose alpha and beta "
				   "Lewis structures the search places apart", log);
		//one NAO basis from the total density for the NPA tables and both spins' searches, as in NBO 7's per-spin tables
		const auto t_nao = clock();
		const dMatrix2 total_density = to_dmatrix(MatrixXd(to_eigen(in.density[0]) + to_eigen(in.density[1])));
		const dMatrix2 spin_density = to_dmatrix(MatrixXd(to_eigen(in.density[0]) - to_eigen(in.density[1])));
		const NAOResult total_nao = build_naos(total_density, in.overlap, in.ao, atoms, ecp);
		res.nao_seconds = std::chrono::duration<double>(clock() - t_nao).count();
		const dMatrix2 spin_nao[2] = { nao_density(in.density[0], in.overlap, total_nao.C), nao_density(in.density[1], in.overlap, total_nao.C) };
		for (int s = 0; s < 2; s++) {
			double electrons = 0.0;
			for (size_t i = 0; i < spin_nao[s].extent(0); i++) electrons += spin_nao[s](i, i);
			analyse_spin(res, in, total_nao, bondable, nrt_bondable, xyz, in.density[s],
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
			vec spin_occ(total_nao.orbitals.size(), 0.0);
			for (size_t i = 0; i < spin_occ.size(); i++) spin_occ[i] = spin_nao[s](i, i);
			fill_nao_table(s ? res.nao_beta : res.nao_alpha, total_nao, spin_occ,
						   static_cast<int>(in.fock.size()) > s ? nao_operator(in.fock[s], total_nao.C)
																: dMatrix2());
		}
		const dMatrix2 spin_pop = nao_density(spin_density, in.overlap, total_nao.C);
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
				if (total_nao.orbitals[i].atom == x.index) p.spin_density += spin_pop(i, i);
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
	//the parent the search starts from: dense in column blocks up to 40 atoms, one row per atom beyond
	for (const NboTopo& t : n.parent_topo) {
		const int na = static_cast<int>(t.matrix.size());
		out << title("Parent Lewis structure", t.spin);
		if (na <= 40) {
			out << "  valence lone pairs on the diagonal, bond orders off it\n";
			for (int c0 = 0; c0 < na; c0 += 12) {
				out << "\n        ";
				for (int b = c0; b < std::min(na, c0 + 12); b++) out << setw(6) << label(b + 1);
				out << "\n";
				for (int a = 0; a < na; a++) {
					out << "  " << left << setw(6) << label(a + 1) << right;
					for (int b = c0; b < std::min(na, c0 + 12); b++)
						out << setw(6) << (t.matrix[a][b] ? std::to_string(t.matrix[a][b]) : ".");
					out << "\n";
				}
			}
		}
		else {
			out << "  atom   lone pairs   partner:bond order\n";
			for (int a = 0; a < na; a++) {
				out << "  " << left << setw(7) << label(a + 1) << right << setw(10) << t.matrix[a][a] << "  ";
				for (int b = 0; b < na; b++)
					if (b != a && t.matrix[a][b]) out << " " << label(b + 1) << ":" << t.matrix[a][b];
				out << "\n";
			}
		}
	}
	if (!n.parent_topo.empty())
		out << "\nParent Lewis structure as NRT input (JSON):\n" << parent_lewis_json(r) << "\n";

	out << "\nNatural resonance theory (in house):\n"
		<< "  " << n.structures_used << " of " << n.structures_found
		<< " resonance structures carry the fit, D(0) = " << fixed << setprecision(5) << n.d_0
		<< ", D(w) = " << n.d_w << "\n";
	for (const NboNrtD& d : n.spin_d)
		out << "  in NBO's units" << (d.spin.empty() ? "" : " (" + d.spin + ")") << ": D(0) = " << d.d_0 << ", D(w) = " << d.d_w << "\n";
	//an open shell runs the search once per spin, and both runs push the same budget line
	svec seen;
	for (const std::string& s : n.notes)
		if (std::find(seen.begin(), seen.end(), s) == seen.end()) {
			seen.push_back(s);
			out << "  " << s << "\n";
		}

	for (const std::string& spin : spins_of(n.weights)) {
		out << title("Resonance weights", spin) << "  structure   weight %   change from the parent\n";
		//structure 0 is the parent, listed only when the fit drops it (a retained one is its own row)
		const bool retained = std::any_of(n.weights.begin(), n.weights.end(), [&](const NboResonanceWeight& w) {
			return w.spin == spin && w.changes == "parent";
		});
		for (const NboNrtD& d : n.spin_d)
			if (!retained && d.spin == spin && d.parent_percent >= 0.0)
				out << setw(11) << 0 << setprecision(2) << setw(11) << d.parent_percent << "   parent\n";
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
		out << title("Natural bond orders", spin) << "  atom   atom        total  covalent     ionic   range at D(w) + 0.1 %\n"
			<< setprecision(4);
		for (const NboBondOrder& b : n.bond_orders)
			if (b.spin == spin && !b.diagonal) {
				out << "  " << left << setw(7) << label(b.atom1) << setw(7) << label(b.atom2) << right
					<< setw(10) << b.total << setw(10) << b.covalent << setw(10) << b.ionic;
				if (b.range_lo >= 0.0) out << "   " << b.range_lo << " - " << b.range_hi;
				out << "\n";
			}
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
	if (!r.e2_avg.structures.empty()) {
		const NboE2Average& a = r.e2_avg;
		out << "\nSecond-order E(2) averaged over NRT structures " << e2_avg_window(a)
			<< " (in house; channels sum orbitals of one type on the same centres, weights renormalised):\n";
		for (size_t i = 0; i < a.structures.size(); i++)
			out << "  structure " << setw(3) << a.structures[i] << "  weight " << setprecision(4) << setw(7) << a.weights[i]
				<< "  rho(NL) " << setprecision(5) << a.rho_nl[i] << "\n";
		out << "  " << left << setw(21) << "donor" << setw(21) << "acceptor" << right << setw(10) << "<E(2)>"
			<< setw(9) << "min" << setw(9) << "max" << setw(9) << "present" << setw(9) << "<q>/e" << "\n";
		for (const NboE2Channel& c : a.channels)
			out << "  " << left << setw(21) << c.donor << setw(21) << c.acceptor << right << setprecision(2)
				<< setw(10) << c.mean_kcal << setw(9) << c.min_kcal << setw(9) << c.max_kcal << setw(9) << c.present
				<< setprecision(4) << setw(9) << c.mean_e << "\n";
	}
	print_nrt(r, out);
}

WFN nbo_cube_wfn(const WFN& wavy, const NboResults& r, const NboOrbitalSpace& space, const std::string& spec,
				 const bool nlmo, std::vector<std::pair<int, std::string>>& picked, std::ostream& log)
{
	const char* opt = nlmo ? "-nlmo_cube" : "-nbo_cube";
	const NboLewis& lewis = space.lewis;
	const int n = static_cast<int>(space.nao.C.extent(0));
	err_checkf(!wavy.get_is_unrestricted(), std::string(opt) + ": closed-shell wavefunctions only", log);
	err_checkf(static_cast<int>(space.lcaomo.extent(0)) == n,
			   std::string(opt) + ": the FILE47 carries no $LCAOMO to reach the primitives through", log);
	//selection as printed numbers, 1-based
	ivec sel;
	std::stringstream ss(spec);
	std::string tok;
	while (std::getline(ss, tok, ',')) {
		if (tok.empty()) continue;
		if (tok == "lewis") {
			for (int k = 1; k <= lewis.n_lewis; k++) sel.push_back(k);
			continue;
		}
		if (std::all_of(tok.begin(), tok.end(), [](unsigned char c) { return std::isdigit(c); })) {
			const int k = std::stoi(tok);
			err_checkf(k >= 1 && k <= static_cast<int>(r.orbitals.size()),
					   std::string(opt) + ": NBO " + tok + " does not exist, there are " + std::to_string(r.orbitals.size()), log);
			sel.push_back(k);
			continue;
		}
		std::string label = tok;
		std::replace(label.begin(), label.end(), ':', ' ');
		const size_t before = sel.size();
		for (const NboOrbital& o : r.orbitals)
			if (nbo_label(o) == label) sel.push_back(o.index);
		//a label from another resonance structure (BD Pd1-Zn2 where this Lewis structure has LP Zn2) is skipped,
		//not fatal: one batch can then ask for every orbital any structure might hold
		if (sel.size() == before) log << "Warning: " << opt << ": no NBO labelled " << label << " in this Lewis structure, skipped" << std::endl;
	}
	err_checkf(!sel.empty(), std::string(opt) + ": none of the requested NBOs exist", log);
	err_checkf(static_cast<int>(sel.size()) <= wavy.get_nmo(), std::string(opt) + ": more orbitals than MO slots", log);
	ivec internal(sel.size());
	for (size_t k = 0; k < sel.size(); k++) {
		internal[k] = lewis.printed[sel[k] - 1];
		//the Lewis set comes first in lewis.orbitals and only it has NLMOs
		err_checkf(!nlmo || internal[k] < lewis.n_lewis, std::string(opt) + ": NBO " + std::to_string(sel[k]) + " (" +
				   nbo_label(r.orbitals[sel[k] - 1]) + ") is non-Lewis and has no NLMO", log);
	}

	//orbitals in the NAO basis, one column per selection
	MatrixXd V(n, sel.size());
	if (!nlmo)
		for (size_t k = 0; k < sel.size(); k++)
			for (int i = 0; i < n; i++) V(i, k) = lewis.orbitals[internal[k]].coefficients[i];
	else {
		//Gamma / 2 projects on the occupied space of a closed-shell determinant; tr P - tr P^2 says how far it is from one
		const MatrixXd P = to_eigen(lewis.gamma) / 2.0;
		log << opt << ": occupied-space projector defect tr(P) - tr(P^2) = " << std::scientific << std::setprecision(2)
			<< P.trace() - P.squaredNorm() << std::defaultfloat << std::endl;
		const int nl = lewis.n_lewis;
		MatrixXd B(n, nl);
		for (int j = 0; j < nl; j++)
			for (int i = 0; i < n; i++) B(i, j) = lewis.orbitals[j].coefficients[i];
		const MatrixXd M = P * B;
		const MatrixXd L = M * sym_power(MatrixXd(M.transpose() * M), -0.5);
		for (size_t k = 0; k < sel.size(); k++) V.col(k) = L.col(internal[k]);
	}

	//NAO -> AO -> MO expansion t = C_MO S x -> primitives
	const MatrixXd T = to_eigen(space.lcaomo) * to_eigen(space.overlap) * to_eigen(space.nao.C) * V;
	WFN w(wavy);
	const int nmo = std::min(wavy.get_nmo(), n);
	picked.clear();
	for (size_t k = 0; k < sel.size(); k++) {
		const NboOrbital& o = r.orbitals[sel[k] - 1];
		picked.emplace_back(sel[k], nbo_label(o));
		//NAO weights per atom, the largest first
		std::vector<std::pair<double, int>> share(space.nao.atoms.size());
		for (size_t a = 0; a < share.size(); a++) share[a] = { 0.0, static_cast<int>(a) };
		for (int i = 0; i < n; i++) share[space.nao.orbitals[i].atom].first += V(i, k) * V(i, k);
		std::sort(share.rbegin(), share.rend());
		double norm = 0.0;
		for (int m = 0; m < nmo; m++) norm += T(m, k) * T(m, k);
		log << (nlmo ? "NLMO " : "NBO ") << sel[k] << " " << picked.back().second << std::fixed << std::setprecision(2);
		if (nlmo) {
			double parent = 0.0;
			for (int i = 0; i < n; i++) parent += V(i, k) * lewis.orbitals[internal[k]].coefficients[i];
			log << ", " << 100.0 * parent * parent << " % parent NBO";
		}
		log << ", atoms (%):";
		for (size_t a = 0; a < share.size() && share[a].first >= 0.005; a++)
			log << " " << constants::atnr2letter(space.nao.atoms[share[a].second].Z) << share[a].second + 1 << " "
				<< 100.0 * share[a].first;
		log << std::scientific << std::setprecision(1) << ", norm lost in the MO expansion " << 1.0 - norm
			<< std::defaultfloat << std::endl;
		for (int p = 0; p < wavy.get_nex(); p++) {
			double v = 0.0;
			for (int m = 0; m < nmo; m++) v += T(m, k) * wavy.get_MO_coef(m, p);
			w.set_MO_coef(static_cast<int>(k), p, v);
		}
	}
	return w;
}
