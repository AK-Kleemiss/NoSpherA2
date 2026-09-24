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
    //--------------------------------------------------------------------------------------
    // small linear-algebra helpers
    //--------------------------------------------------------------------------------------
    //These are copies of the ones in nao.cpp's anonymous namespace.  Twenty lines duplicated is
    //cheaper than a third header that two files under active parallel development both include.

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
        //A spin with nothing in it is a legitimate input: the beta spin of a hydrogen atom has no
        //occupied NAO, so the Lewis set is empty and this is called on a 0 x 0 matrix.  Eigen's
        //maxCoeff() on an empty expression is undefined - in a release build it reads past the end
        //and segfaults, which is how -nbo_native on a one-electron wavefunction died.  The power of
        //an empty matrix is that matrix, and every product below is already well defined for it.
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

    //--------------------------------------------------------------------------------------
    // FILE47 reading
    //--------------------------------------------------------------------------------------

    //The whole body of a $SECTION, without the keyword and without the $END.
    std::string section_body(const std::string& text, const std::string& name)
    {
        const size_t start = text.find(name);
        if (start == std::string::npos) return std::string();
        const size_t after = start + name.size();
        const size_t end = text.find("$END", after);
        return text.substr(after, (end == std::string::npos ? text.size() : end) - after);
    }

    //Every number of a section body, ignoring the "CENTER =" style tags around them.
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

    //FILE47 packs a symmetric matrix as the upper triangle in column-major order, which reads as
    //the lower triangle in row-major order - the same n(n+1)/2 values either way.
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

    //NBO's LABEL codes, in the order write_nbo() emits them: the components of a shell run
    //m = 0, +1, -1, +2, -2, ... and the code order below is that order, not the numeric one.
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

    //--------------------------------------------------------------------------------------
    // labels
    //--------------------------------------------------------------------------------------

    //The angular label NBO prints for component c of a shell of angular momentum l, in the
    //m = 0, +1, -1, ... order the FILE47 uses.
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

    //Where NBO puts that component in its NAO table: s; px, py, pz; dxy, dxz, dyz, dx2y2, dz2;
    //f(0), f(c1), f(s1), ... - so for p and d it is not the FILE47 order.
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

    //The parser in nbo_run.cpp stores whitespace-collapsed descriptions, and compare_nbo_results()
    //matches orbitals on exactly that string, so the native side has to produce the collapsed form
    //of NBO's fixed-width layout rather than something merely equivalent.
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

    //--------------------------------------------------------------------------------------
    // the search
    //--------------------------------------------------------------------------------------

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

    //Smallest share of a two-centre candidate that may sit on its minor centre.  Tuned on the 22
    //reference molecules: 0.05 lets formate and ozone buy an electron pair with a "bond" that is
    //really a lone pair, 0.15 reproduces NBO's Lewis topology on both without moving any molecule
    //that already agreed.  Override with NBO_BOND_FLOOR if a system ever needs it retuned.
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

NboLewis nbo_search(const NAOResult& nao, const dMatrix2& gamma, const bvec2& bondable,
                    const int n_pairs, const double scale, const NboOptions& options)
{
    const int n = static_cast<int>(gamma.extent(0));
    const int natoms = static_cast<int>(nao.atoms.size());
    const MatrixXd G0 = to_eigen(gamma);
    MatrixXd G = G0;
    const std::vector<AtomIndices> idx = atom_indices(nao);

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
        //depletion: take the orbital's current occupancy out of the working density, so the next
        //block search sees only what is left
        const double occ = v.dot(G * v);
        G -= occ * v * v.transpose();
    };

    //1. cores.  A core NAO is already a one-centre orbital of occupancy ~2 and NBO keeps it as it
    //is; taking the leading eigenvector of the core block instead moves nothing measurable.
    for (int a = 0; a < natoms; a++)
        for (const int i : idx[a].core) {
            VectorXd v = VectorXd::Zero(n);
            v(i) = 1.0;
            accept(v, "CR", { a });
            res.topo[a][a]++;
        }

    //2. the threshold ladder.  At each threshold every one-centre block is searched before any
    //two-centre one - that preference is what makes a lone pair a lone pair rather than half of a
    //bond - and the accepted orbital is always the most occupied candidate found.
    //
    //An atom may not carry more valence orbitals - lone pairs plus bonds, each bond of a multiple
    //bond counted once - than it has valence NAOs: four for a main-group atom, nine for a transition
    //metal.  That cap is what makes NBO's Lewis structure for SF6 four bonds and two fluoride lone
    //pairs rather than six bonds: the sixth S-F candidate has occupancy 1.99 and wins on occupancy,
    //it just does not fit in sulphur's octet.  Without it nine of the 22 reference molecules came out
    //with a different topology - always one bond too many on the central atom - and with it two.
    //If the cap cannot be satisfied at all the ladder is run again without it, so a genuinely
    //hypervalent density still gets a complete Lewis set rather than a truncated one.
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
            double best = t;
            VectorXd bv;
            int ba = -1;
            for (int a = 0; a < natoms; a++) {
                if (idx[a].valence.empty()) continue;
                if (!relax && used[a] >= cap[a]) continue;
                VectorXd v;
                const double lam = leading_block(G, idx[a].all, v);
                if (lam > best) { best = lam; bv = v; ba = a; }
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
            //every bondable atom pair, O(N^2) small diagonalisations per accepted bond; at
            //reference scale the whole search is milliseconds
            for (int a = 0; a < natoms; a++) {
                if (idx[a].valence.empty()) continue;
                if (!relax && used[a] >= cap[a]) continue;
                for (int b = a + 1; b < natoms; b++) {
                    if (idx[b].valence.empty()) continue;
                    if (!relax && used[b] >= cap[b]) continue;
                    if (!bondable.empty() && !bondable[a][b]) continue;
                    ivec pair = idx[a].all;
                    pair.insert(pair.end(), idx[b].all.begin(), idx[b].all.end());
                    VectorXd v;
                    const double lam = leading_block(G, pair, v);
                    if (lam <= bestp) continue;
                    const double wa = weight_on(v, idx[a].all);
                    const double wb = weight_on(v, idx[b].all);
                    //A candidate that sits almost entirely on one centre is a lone pair, not a
                    //bond, and accepting it as a bond costs a whole electron pair of the Lewis
                    //structure: formate took a 93 %-on-oxygen "third C-O bond" at occupancy 1.977
                    //in place of the reference's third lone pair on O, because the ladder always
                    //accepts the most occupied candidate and that one wins at the top rung.  The
                    //floor is the polarity at which NBO stops calling something a bond.
                    if (std::min(wa, wb) < bond_minority_floor() * (wa + wb)) continue;
                    bestp = lam;
                    bvp = v;
                    pa = a;
                    pb = b;
                }
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

    //3. self consistency.  The ladder decides one orbital at a time out of a density that still
    //holds every orbital found after it, so each accepted orbital carries the bias of the ones it
    //did not know about.  Sweeping every orbital against the density with all the *other* accepted
    //ones removed, to convergence, takes that bias out; on water's O-H bond it moves the occupancy
    //from 1.9857 to within 1e-4 of NBO's 1.99963, which is the difference between failing and
    //passing the 2e-3 tolerance.
    {
        vec occ(res.n_lewis, 0.0);
        MatrixXd sum = MatrixXd::Zero(n, n);
        for (int j = 0; j < res.n_lewis; j++) {
            occ[j] = vectors[j].dot(G0 * vectors[j]);
            sum += occ[j] * vectors[j] * vectors[j].transpose();
        }
        for (int sweep = 0; sweep < 200; sweep++) {
            double change = 0.0;
            for (int j = 0; j < res.n_lewis; j++) {
                //a core orbital is a single NAO and stays one - NBO prints 2.00000 for it either way
                if (res.orbitals[j].type == "CR") continue;
                sum -= occ[j] * vectors[j] * vectors[j].transpose();
                ivec sub;
                for (const int a : res.orbitals[j].centers)
                    sub.insert(sub.end(), idx[a].all.begin(), idx[a].all.end());
                VectorXd v;
                leading_block(MatrixXd(G0 - sum), sub, v);
                if (v.dot(vectors[j]) < 0.0) v = -v;
                change = std::max(change, (v - vectors[j]).norm());
                vectors[j] = v;
                occ[j] = v.dot(G0 * v);
                sum += occ[j] * v * v.transpose();
            }
            if (change < 1e-10) break;
        }
        //4. the refined set is still not orthogonal - two bonds at the same atom overlap by a fifth.
        //An occupancy-weighted symmetric orthogonalisation (OWSO) spreads that symmetrically, which
        //is what keeps two equivalent bonds equivalent while letting the occupied ones keep their
        //shape at the expense of the empty ones.
        MatrixXd VL(n, res.n_lewis);
        for (int j = 0; j < res.n_lewis; j++) VL.col(j) = vectors[j];
        VectorXd w(res.n_lewis);
        for (int j = 0; j < res.n_lewis; j++) w(j) = std::max(occ[j], 1e-6);
        const MatrixXd S = VL.transpose() * VL;
        VL = VL * w.asDiagonal() *
             sym_power(MatrixXd(w.asDiagonal() * S * w.asDiagonal()), -0.5);
        for (int j = 0; j < res.n_lewis; j++) vectors[j] = VL.col(j);
    }

    //5. the valence antibonds.  A bond is c_A h_A + c_B h_B over two normalised hybrids, so its
    //antibond - the other vector of that same two-dimensional space - is c_B h_A - c_A h_B.  It is
    //constructed, not searched: it has almost no occupancy to find it by.
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

    //6. everything the Lewis set and the antibonds do not span.  Their hybrids reach a little way
    //into the extra-valence NAOs - the reference's own LP(1) on water's oxygen carries 0.17% d
    //character, and confining the search to the natural minimal basis costs 0.014 e of bond
    //occupancy - so the complement is taken over the whole NAO space.  It has to be taken with a
    //projector rather than by Schmidt: the antibonds are not orthogonal to the Lewis set, and
    //projecting against a non-orthonormal set silently produces vectors of the wrong length (it put
    //a 202%-weight "LV" of occupancy 1.03 into water).
    {
        const int k = res.n_lewis + static_cast<int>(non_lewis.size());
        MatrixXd M(n, k);
        for (int j = 0; j < res.n_lewis; j++) M.col(j) = vectors[j];
        for (size_t j = 0; j < non_lewis.size(); j++)
            M.col(res.n_lewis + static_cast<int>(j)) = non_lewis[j];
        const MatrixXd Q = MatrixXd::Identity(n, n) -
            M * sym_power(MatrixXd(M.transpose() * M), -1.0) * M.transpose();
        //Q's own eigenvectors would do as a basis of the complement, but they come out of a
        //degenerate eigenspace and are therefore arbitrary - nothing ties one of them to an atom.
        //Projecting the NAO unit vectors instead and taking them in order of the largest surviving
        //norm (a pivoted Gram-Schmidt) picks a basis each of whose members belongs to one NAO, and
        //so to one atom; that is what makes the printed RY orbitals one-centre.
        std::vector<VectorXd> extra;
        ivec owner;
        std::vector<VectorXd> residual(n);
        for (int i = 0; i < n; i++) residual[i] = Q.col(i);
        bvec used(n, false);
        for (int step = 0; step < n - k; step++) {
            int pivot = -1;
            double best = 0.0;
            for (int i = 0; i < n; i++) {
                if (used[i]) continue;
                const double nrm = residual[i].norm();
                if (nrm > best) { best = nrm; pivot = i; }
            }
            err_checkf(pivot >= 0 && best > 1e-8,
                       "NBO: the residual space is smaller than the orbital count demands", std::cout);
            used[pivot] = true;
            const VectorXd u = residual[pivot] / best;
            extra.push_back(u);
            owner.push_back(nao.orbitals[pivot].atom);
            for (int i = 0; i < n; i++)
                if (!used[i]) residual[i] -= u.dot(residual[i]) * u;
        }
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
                //a residual orbital that still lives in the atom's valence shell is a lone
                //vacancy; one pushed out of it is Rydberg
                f.type = (weight_on(v, idx[a].valence) > 0.5) ? "LV" : "RY";
                f.centers = { a };
                non_lewis.push_back(v);
                non_lewis_fn.push_back(f);
            }
        }
    }

    for (size_t j = 0; j < non_lewis.size(); j++) {
        NboFunction f = non_lewis_fn[j];
        f.occupancy = non_lewis[j].dot(G0 * non_lewis[j]);
        if (f.type != "BD*") {
            f.multiplicity = 1;
            for (const NboFunction& o : res.orbitals)
                if (o.type == f.type && o.centers == f.centers) f.multiplicity++;
        }
        res.orbitals.push_back(f);
        vectors.push_back(non_lewis[j]);
    }

    //8. coefficients, polarisation and l character.  The occupancies are taken here, from the final
    //vectors: the Lewis ones were refined and then orthogonalised after they were accepted.
    for (size_t j = 0; j < res.orbitals.size(); j++) {
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
                      const NboOptions& options, std::vector<NboLewis>& out_lewis)
    {
        const dMatrix2 gamma = nao_density(density, in.overlap, nao.C);
        const dMatrix2 fock_nao = fock.extent(0) ? nao_operator(fock, nao.C) : dMatrix2();
        //A stage nobody times is a stage nobody can make faster: sucrose spent 50 of its 54 s
        //somewhere in here while the only phase line in the log was NRT's own 1.4 s.
        const auto clock = [] { return std::chrono::steady_clock::now(); };
        auto t = clock();
        NboLewis lewis = nbo_search(nao, gamma, bondable, n_pairs, scale, options);
        res.search_seconds += std::chrono::duration<double>(clock() - t).count();
        t = clock();
        std::vector<NboE2Entry> e2 = nbo_e2(lewis, fock_nao, options.e2_threshold_kcal);
        res.e2_seconds += std::chrono::duration<double>(clock() - t).count();
        //before the renumbering below, while e2's indices still point into lewis.orbitals
        if (options.nrt)
            native_nrt(res.nrt, nao, lewis, e2, nrt_bondable, options, spin, scale, std::cout);

        //NBO numbers its NBOs by type group: the Lewis set in the order it was found, then the
        //non-Lewis set as LV, BD*, RY
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
                //NBO's sign convention: the first centre positive, the second one following the
                //bond/antibond pair, so an antibond carries the minus
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

    //The NAO table in NBO's own order: per atom, per l, components in NBO's printing order, and
    //inside one component the shells by descending occupancy.
    void fill_nao_table(NboResults& res, const NAOResult& nao, const vec& occupancy,
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
            res.nao.push_back(n);
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
                     in.fock.empty() ? dMatrix2() : in.fock[0], n_pairs, 2.0, "", options, lewis);
        vec occ(nao.orbitals.size(), 0.0);
        for (size_t i = 0; i < occ.size(); i++) occ[i] = nao.orbitals[i].occupation;
        fill_nao_table(res, nao, occ, in.fock.empty() ? dMatrix2() : nao_operator(in.fock[0], nao.C));
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
        //NBO analyses the two spin densities independently and prints a spin-summed NAO table
        const auto t_nao = clock();
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
                         options, lewis);
        }
        //the spin-summed NAO table needs one set of NAOs; the alpha set carries the labels and the
        //occupancies are the two spin occupancies of the same NAO index, which is an approximation
        //to what NBO prints (it re-derives a spin-averaged set)
        vec occ(a_nao.orbitals.size(), 0.0);
        for (size_t i = 0; i < occ.size(); i++)
            occ[i] = a_nao.orbitals[i].occupation + b_nao.orbitals[i].occupation;
        fill_nao_table(res, a_nao, occ,
                       in.fock.empty() ? dMatrix2() : nao_operator(in.fock[0], a_nao.C));
        for (size_t a = 0; a < a_nao.atoms.size(); a++) {
            const NAOAtom& x = a_nao.atoms[a];
            const NAOAtom& y = b_nao.atoms[a];
            NboAtomPopulation p;
            p.element = constants::atnr2letter(x.Z);
            p.index = x.index + 1;
            p.core = x.core + y.core;
            p.valence = x.valence + y.valence;
            p.rydberg = x.rydberg + y.rydberg;
            p.total = x.population + y.population;
            p.charge = x.Z_eff - p.total;
            p.spin_density = x.population - y.population;
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

//-nrt spends most of a run's time and used to report nowhere a reader looks: the resonance weights,
//the bond orders and the valencies went into <stem>.native.nbo.json and NoSpherA2.log carried only
//the two NRT citations, so a 6.6 s search on a 14-atom complex was indistinguishable from an
//ignored flag.  The layout follows gennbo's own headings on purpose, so the two routes can be read
//side by side; every number printed here is the one written to the JSON, unrounded there.
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

    out << "\n NATURAL RESONANCE THEORY ANALYSIS (in house):\n\n"
        << " " << n.structures_used << " of " << n.structures_found
        << " resonance structures carry the fit, D(0) = " << fixed << setprecision(5) << n.d_0
        << ", D(w) = " << n.d_w << "\n";
    for (const std::string& s : n.notes) out << "   " << s << "\n";

    if (!n.weights.empty()) {
        out << "\n  RS   Weight(%)   Added(Removed)\n"
            << " ---------------------------------------------------------------------------------\n";
        for (const NboResonanceWeight& w : n.weights)
            out << setw(4) << w.structure << setprecision(2) << setw(11) << w.weight_percent << "   "
                << (w.spin.empty() ? "" : w.spin + ": ") << w.changes << "\n";
    }
    if (!n.bond_orders.empty()) {
        out << "\n Natural Bond Order\n"
            << "   Atom  Atom      Total   Covalent      Ionic\n"
            << " ---------------------------------------------------------------------------------\n";
        for (const NboBondOrder& b : n.bond_orders) {
            out << "  " << left << setw(6) << label(b.atom1) << setw(6)
                << (b.diagonal ? std::string() : label(b.atom2)) << right << fixed << setprecision(4)
                << setw(11) << b.total;
            if (b.diagonal)
                out << "        ---        ---";
            else
                out << setw(11) << b.covalent << setw(11) << b.ionic;
            if (!b.spin.empty()) out << "  " << b.spin;
            out << "\n";
        }
    }
    if (!n.valencies.empty()) {
        out << "\n Natural Atomic Valencies\n"
            << "   Atom    Valency  Covalency  Electroval.   Electrons\n"
            << " ---------------------------------------------------------------------------------\n";
        for (const NboValency& v : n.valencies) {
            out << "  " << left << setw(6) << (v.element + std::to_string(v.atom)) << right << fixed
                << setprecision(4) << setw(11) << v.valency << setw(11) << v.covalency << setw(11)
                << v.electrovalency << setw(12) << v.electron_count;
            if (!v.spin.empty()) out << "  " << v.spin;
            out << "\n";
        }
    }
    for (const std::string& s : n.symmetry_forms) out << " " << s << "\n";
}

void print_nbo(const NboResults& r, std::ostream& out)
{
    using namespace std;
    //fixed/setprecision below stay on the stream after this table, so everything printed through it
    //afterwards would carry two decimals
    const ostream_format_guard restore_format(out);
    //The populations were JSON-only for the same reason the resonance tables were: nobody wrote the
    //branch.  NPA charges are the most-read line of an NBO run, and on -nbo_native they appeared
    //nowhere in NoSpherA2.log.  The per-NAO table stays in <stem>.native.nbo.json - it is 189 rows
    //on a 14-atom complex and the populations are its summary.
    if (!r.npa.empty()) {
        const bool spin = r.npa.front().has_spin_density;
        out << "\n NATURAL POPULATION ANALYSIS (in house):\n\n"
            << "   Atom      Charge       Core    Valence    Rydberg      Total"
            << (spin ? "   Spin dens.\n" : "\n")
            << " ---------------------------------------------------------------------------------\n";
        double charge_sum = 0.0, total_sum = 0.0;
        for (const NboAtomPopulation& p : r.npa) {
            out << "  " << left << setw(4) << (p.element + std::to_string(p.index)) << right << fixed
                << setprecision(5) << setw(12) << p.charge << setw(11) << p.core << setw(11)
                << p.valence << setw(11) << p.rydberg << setw(11) << p.total;
            if (p.has_spin_density) out << setw(13) << p.spin_density;
            out << "\n";
            charge_sum += p.charge;
            total_sum += p.total;
        }
        //the two sums are the check a reader can make on the spot: the charges add to the molecular
        //charge and the populations to the number of electrons the wavefunction carries
        out << "  " << left << setw(4) << "sum" << right << setw(12) << charge_sum << setw(44)
            << total_sum << "\n";
    }
    out << "\n NATURAL BOND ORBITAL ANALYSIS (in house):\n\n"
        << "                                                      Principal Delocalizations\n"
        << "  NBO                         Occupancy    Energy\n"
        << " ---------------------------------------------------------------------------------\n";
    for (const NboOrbital& o : r.orbitals) {
        out << setw(4) << o.index << ". " << left << setw(22) << o.description << right
            << fixed << setprecision(5) << setw(12) << o.occupancy << setw(12) << o.energy;
        if (!o.spin.empty()) out << "  " << o.spin;
        out << "\n";
        for (const NboHybrid& h : o.hybrids)
            out << "              " << fixed << setprecision(2) << setw(7) << h.weight_percent
                << "% " << setw(2) << h.element << setw(3) << h.center << "  s(" << setw(6) << h.s
                << "%)p" << setw(7) << h.p << "%  d" << setw(7) << h.d << "%  f" << setw(6) << h.f
                << "%\n";
    }
    if (!r.e2.empty()) {
        out << "\n SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS\n\n"
            << "     Donor NBO              Acceptor NBO            E(2)   E(NL)-E(L)  F(L,NL)\n"
            << " ---------------------------------------------------------------------------------\n";
        for (const NboE2Entry& e : r.e2)
            out << " " << left << setw(22) << e.donor << setw(22) << e.acceptor << right << fixed
                << setprecision(2) << setw(8) << e.energy_kcal << setprecision(3) << setw(11)
                << e.e_diff << setw(10) << e.fij << (e.spin.empty() ? "" : "  " + e.spin) << "\n";
    }
    print_nrt(r, out);
}
