#include "pch.h"
#include "nao.h"
#include "wfn_class.h"
#include "integration_params.h"
#include "libCintMain.h"
#include "constants.h"
#include "nos_math.h"
#include "citations.h"

#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace
{
    //--------------------------------------------------------------------------------------
    // small linear-algebra helpers
    //--------------------------------------------------------------------------------------

    MatrixXd to_eigen(const dMatrix2 &m)
    {
        const int n = static_cast<int>(m.extent(0)), c = static_cast<int>(m.extent(1));
        MatrixXd out(n, c);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < c; j++)
                out(i, j) = m(i, j);
        return out;
    }

    dMatrix2 to_dmatrix(const MatrixXd &m)
    {
        dMatrix2 out(m.rows(), m.cols());
        for (int i = 0; i < m.rows(); i++)
            for (int j = 0; j < m.cols(); j++)
                out(i, j) = m(i, j);
        return out;
    }

    //M^p for a symmetric positive semi-definite M.  Eigenvalues below rel_floor * max are
    //dropped for a negative power, which turns an inverse into a pseudo-inverse rather than
    //letting a near-linear dependence blow the transform up.
    MatrixXd sym_power(const MatrixXd &M, const double p, const double rel_floor = 1e-10)
    {
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(M);
        VectorXd w = es.eigenvalues();
        const double cut = rel_floor * std::max(w.maxCoeff(), 1e-300);
        VectorXd f(w.size());
        for (int i = 0; i < w.size(); i++)
            f(i) = (w(i) > cut) ? std::pow(w(i), p) : (p < 0.0 ? 0.0 : std::pow(std::max(w(i), 0.0), p));
        return es.eigenvectors() * f.asDiagonal() * es.eigenvectors().transpose();
    }

    //T = W (W S W)^-1/2, the occupancy-weighted symmetric orthogonalisation of Reed/Weinhold.
    //It reduces to S^-1/2 for equal weights and to the identity for S = 1, and it is what makes
    //the strongly occupied orbitals keep their shape while the diffuse ones absorb the
    //orthogonalisation tails.
    MatrixXd owso(const MatrixXd &S, const VectorXd &weights)
    {
        VectorXd w = weights;
        //ponytail: relative floor on the weights. Rydberg occupancies span many decades and a
        //weight of 1e-12 next to one of 1e-3 makes W S W numerically singular. Drop the floor
        //only if a case shows the clamping matters.
        const double wmax = std::max(w.maxCoeff(), 1e-300);
        for (int i = 0; i < w.size(); i++)
            w(i) = std::max(w(i), 1e-6 * wmax);
        const MatrixXd WSW = w.asDiagonal() * S * w.asDiagonal();
        //W S W carries the weights squared, so its eigenvalues span (wmin/wmax)^2 even for a
        //perfectly conditioned S, and a fixed relative floor in sym_power would delete exactly
        //the columns of smallest weight - which is how a full-rank basis lost a Rydberg
        //direction on nh3li and 3.5e-4 electrons with it.  Scale the floor with the weights.
        const double s = w.minCoeff() / wmax;
        return w.asDiagonal() * sym_power(WSW, -0.5, 1e-10 * s * s);
    }

    //--------------------------------------------------------------------------------------
    // periodic table bookkeeping
    //--------------------------------------------------------------------------------------

    int period_of(const int Z)
    {
        if (Z <= 2) return 1;
        if (Z <= 10) return 2;
        if (Z <= 18) return 3;
        if (Z <= 36) return 4;
        if (Z <= 54) return 5;
        if (Z <= 86) return 6;
        return 7;
    }

    //0 = s, 1 = p, 2 = d, 3 = f
    int block_of(const int Z)
    {
        if (Z == 1 || Z == 2 || Z == 3 || Z == 4 || Z == 11 || Z == 12 || Z == 19 || Z == 20 ||
            Z == 37 || Z == 38 || Z == 55 || Z == 56 || Z == 87 || Z == 88)
            return 0;
        if ((Z >= 57 && Z <= 70) || (Z >= 89 && Z <= 102)) return 3;
        if ((Z >= 21 && Z <= 30) || (Z >= 39 && Z <= 48) || (Z >= 71 && Z <= 80) || (Z >= 103 && Z <= 112))
            return 2;
        return 1;
    }

    //The internal spherical basis runs m = 0, +1, -1, +2, -2 (the ORCA/Gaussian order the
    //sph2cart tables and the .47 writer both assume), not libcint's -l .. +l.
    std::string shell_label(const int l, const int m)
    {
        static const char *lc = "spdfghik";
        const char letter = (l < 8) ? lc[l] : '?';
        if (l == 0) return std::string(1, letter);
        if (l == 1) {
            static const char *p[3] = { "z", "x", "y" };
            return std::string(1, letter) + p[m];
        }
        if (l == 2) {
            static const char *d[5] = { "z2", "xz", "yz", "x2y2", "xy" };
            return std::string(1, letter) + d[m];
        }
        const int mv = (m == 0) ? 0 : ((m + 1) / 2) * ((m % 2) ? 1 : -1);
        return std::string(1, letter) + "(" + (mv >= 0 ? "+" : "") + std::to_string(mv) + ")";
    }

    const char *class_label(const NAOClass c)
    {
        switch (c) {
        case NAOClass::Core: return "Cor";
        case NAOClass::Valence: return "Val";
        default: return "Ryd";
        }
    }
}

void natural_minimal_shells(const int Z, int (&n_shell)[4], int (&n_core)[4])
{
    const int p = period_of(Z), b = block_of(Z);
    n_shell[0] = p;
    n_shell[1] = std::max(0, (b == 1) ? p - 1 : p - 2);
    n_shell[2] = (p < 4) ? 0 : std::max(0, (b == 0) ? p - 4 : p - 3);
    n_shell[3] = (p < 6) ? 0 : std::max(0, (b == 0) ? p - 6 : p - 5);
    for (int l = 0; l < 4; l++) n_core[l] = n_shell[l];
    //the valence shell of a block is one s plus the shell the block is named after
    n_core[0] = std::max(0, n_core[0] - 1);
    if (b == 1) n_core[1] = std::max(0, n_core[1] - 1);
    if (b == 2 || b == 3) n_core[2] = std::max(0, n_core[2] - 1);
    if (b == 3) n_core[3] = std::max(0, n_core[3] - 1);
}

std::vector<NAOBasisFunction> spherical_ao_map(const WFN &wavy)
{
    err_checkf(!wavy.get_d_f_switch(),
               "NAO/NPA needs a spherical basis; this wavefunction carries Cartesian d/f shells. "
               "Convert through a spherical .gbw/.molden of the same calculation.",
               std::cout);
    Int_Params params(wavy);
    const ivec bas = params.get_bas();
    const size_t nbas = params.get_nbas();
    std::vector<NAOBasisFunction> map;
    map.reserve(params.get_nao());
    std::map<std::pair<int, int>, int> shell_counter;
    for (size_t s = 0; s < nbas; s++) {
        const int a = bas[8 * s + 0], l = bas[8 * s + 1];
        const int idx = shell_counter[{ a, l }]++;
        for (int m = 0; m <= 2 * l; m++)
            map.push_back(NAOBasisFunction{ a, l, idx, m });
    }
    return map;
}

NAOResult build_naos(const dMatrix2 &P_in, const dMatrix2 &S_in, const std::vector<NAOBasisFunction> &ao,
                     const std::vector<atom> &atoms, const ivec &ecp_electrons)
{
    const int nao = static_cast<int>(ao.size());
    err_checkf(static_cast<int>(P_in.extent(0)) == nao && static_cast<int>(S_in.extent(0)) == nao,
               "NAO: density and overlap do not match the basis map (" +
                   std::to_string(P_in.extent(0)) + "/" + std::to_string(S_in.extent(0)) + " vs " +
                   std::to_string(nao) + ")",
               std::cout);
    const MatrixXd P = to_eigen(P_in), S = to_eigen(S_in);
    const MatrixXd SPS = S * P * S;

    //index lists per (atom, l) group, split by shell and by m
    struct Group {
        int atom = 0, l = 0, nshell = 0;
        std::vector<ivec> idx;  //[shell][m]
    };
    std::map<std::pair<int, int>, Group> groups;
    for (int i = 0; i < nao; i++) {
        Group &g = groups[{ ao[i].atom, ao[i].l }];
        g.atom = ao[i].atom;
        g.l = ao[i].l;
        if (static_cast<int>(g.idx.size()) <= ao[i].shell) g.idx.resize(ao[i].shell + 1);
        g.idx[ao[i].shell].push_back(i);
    }

    //---------------------------------------------------------------- 1. pre-NAOs
    //Inside one atom, shells of different l are orthogonal by symmetry and so are different m,
    //so the only non-trivial block is (atom, l) averaged over m.  Solving
    //(S P S) c = w S c there gives orbitals with c^T S c = 1 whose occupancy is exactly w.
    //Using the atom's own block of P with its own block of S as the metric instead - the net
    //rather than the gross atomic population - moves half an electron per carbon in epoxide, and
    //keeping these orbitals but weighting the orthogonalisation below by the net population
    //c^T P^A c is worse again (0.4 e on epoxide's oxygen).  Gross it is, for both.
    MatrixXd C = MatrixXd::Zero(nao, nao);
    VectorXd pre_occ = VectorXd::Zero(nao);
    std::vector<NAO> orbitals(nao);
    int col = 0;
    std::map<std::pair<int, int>, ivec> group_columns;  //(atom, l) -> column indices, shell-major
    for (auto &kv : groups) {
        Group &g = kv.second;
        g.nshell = static_cast<int>(g.idx.size());
        const int nm = 2 * g.l + 1, ns = g.nshell;
        MatrixXd Sb = MatrixXd::Zero(ns, ns), Pb = MatrixXd::Zero(ns, ns);
        for (int s1 = 0; s1 < ns; s1++)
            for (int s2 = 0; s2 < ns; s2++)
                for (int m = 0; m < nm; m++) {
                    Sb(s1, s2) += S(g.idx[s1][m], g.idx[s2][m]) / nm;
                    Pb(s1, s2) += SPS(g.idx[s1][m], g.idx[s2][m]) / nm;
                }
        //(S P S)^A c = w S^A c becomes the ordinary eigenproblem X (S P S)^A X y = w y with
        //X = (S^A)^-1/2 and c = X y, which keeps c^T S c = 1
        const MatrixXd X = sym_power(Sb, -0.5);
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(X * Pb * X);
        //descending occupancy: the natural order of the shells of a given l
        for (int k = ns - 1; k >= 0; k--) {
            const VectorXd c = X * es.eigenvectors().col(k);
            const double w = es.eigenvalues()(k);
            const int shell = ns - 1 - k;
            for (int m = 0; m < nm; m++) {
                for (int s = 0; s < ns; s++)
                    C(g.idx[s][m], col) = c(s);
                pre_occ(col) = w;
                orbitals[col].atom = g.atom;
                orbitals[col].l = g.l;
                orbitals[col].m = m;
                orbitals[col].shell = shell;
                orbitals[col].n = g.l + 1 + shell;
                group_columns[{ g.atom, g.l }].push_back(col);
                col++;
            }
        }
    }
    err_checkf(col == nao, "NAO: lost basis functions while building pre-NAOs", std::cout);

    //---------------------------------------------------------------- 2. NMB / NRB partition
    for (auto &kv : groups) {
        const Group &g = kv.second;
        int shells[4] = { 0, 0, 0, 0 }, cores[4] = { 0, 0, 0, 0 };
        const int Z = atoms[g.atom].get_charge();
        natural_minimal_shells(Z, shells, cores);
        int n_nmb = (g.l < 4) ? shells[g.l] : 0;
        int n_core = (g.l < 4) ? cores[g.l] : 0;
        //an ECP removes the lowest shells of an l from the basis, so cap on what is there and
        //take the missing ones off the core count
        if (n_nmb > g.nshell) {
            n_core = std::max(0, n_core - (n_nmb - g.nshell));
            n_nmb = g.nshell;
        }
        for (int c : group_columns[{ g.atom, g.l }]) {
            const int sh = orbitals[c].shell;
            orbitals[c].type = (sh < n_core) ? NAOClass::Core
                             : (sh < n_nmb)  ? NAOClass::Valence
                                             : NAOClass::Rydberg;
            //an ECP core shell is gone from the basis, so the remaining ones are labelled from
            //the first shell the basis actually carries
            if (ecp_electrons[g.atom] > 0 && (g.l < 4)) {
                int full_shells[4] = { 0, 0, 0, 0 }, full_cores[4] = { 0, 0, 0, 0 };
                natural_minimal_shells(Z, full_shells, full_cores);
                const int missing = std::max(0, full_shells[g.l] - g.nshell);
                orbitals[c].n += missing;
            }
        }
    }

    //---------------------------------------------------------------- 3. orthogonalisation
    //Three sets in decreasing priority: core, valence, Rydberg.  Each is Schmidt-projected out
    //of everything above it - so a core keeps its shape exactly, which is what makes its
    //occupancy come out at 1.99995 rather than 1.9991 - and then occupancy-weighted
    //symmetrically orthogonalised among its own members.
    //Core and valence stay separate classes: orthogonalising the whole natural minimal basis in
    //one OWSO, the way the 1985 paper reads, moves epoxide and nh3bh3 further from NBO 7 (worst
    //|dq| 0.020 -> 0.022 and 0.019 -> 0.022) and benzene from 0.004 to 0.006.
    ivec cols_by_class[3];
    for (int i = 0; i < nao; i++)
        cols_by_class[static_cast<int>(orbitals[i].type)].push_back(i);
    //(atom, l) -> columns, shell-major, for the m-averaging in step 4 and in the weight update
    std::map<std::pair<int, int>, ivec> l_blocks;
    for (int i = 0; i < nao; i++)
        l_blocks[{ orbitals[i].atom, orbitals[i].l }].push_back(i);

    //The weights stay the pre-NAO occupancies.  Re-running steps 3 and 4 with the weights taken
    //from the orbitals they produce does converge, but to the wrong answer - epoxide's hydrogens
    //end at +0.66 - so the pre-NAO occupancy is the weight, not a first guess at one.
    MatrixXd done(nao, 0);  //everything orthonormalised so far, in S
    for (int cls = 0; cls < 3; cls++) {
        const ivec &cols = cols_by_class[cls];
        if (cols.empty()) continue;
        MatrixXd B(nao, static_cast<int>(cols.size()));
        VectorXd w(static_cast<int>(cols.size()));
        for (size_t j = 0; j < cols.size(); j++) {
            B.col(static_cast<int>(j)) = C.col(cols[j]);
            w(static_cast<int>(j)) = std::max(pre_occ(cols[j]), 0.0);
        }
        if (done.cols() > 0) {
            B -= done * (done.transpose() * S * B);
            for (int j = 0; j < B.cols(); j++)
                B.col(j) /= std::sqrt(std::max(B.col(j).dot(S * B.col(j)), 1e-300));
        }
        B = B * owso(MatrixXd(B.transpose() * S * B), w);
        //the weighted inverse square root leaves the near-zero-weight directions orthonormal only
        //to about 1e-5; one unweighted Loewdin on a matrix that is already I + O(1e-5) cleans that
        //up without moving the occupied orbitals
        B = B * sym_power(MatrixXd(B.transpose() * S * B), -0.5);
        for (size_t j = 0; j < cols.size(); j++) C.col(cols[j]) = B.col(static_cast<int>(j));
        const int old = static_cast<int>(done.cols());
        done.conservativeResize(nao, old + B.cols());
        done.rightCols(B.cols()) = B;
    }

    //---------------------------------------------------------------- 4. natural character
    //The set is orthonormal now but the orthogonalisation mixed the shells of one l, so
    //re-diagonalise inside every (atom, l) block - m-averaged again, which is what keeps the
    //components of a shell degenerate.  Valence and Rydberg go in one block: separating them
    //leaves the valence-Rydberg coupling in place and inflates the Rydberg occupancies tenfold.
    //The mixing is intra-atomic and unitary, so it moves no charge between atoms.
    const MatrixXd Porb = C.transpose() * SPS * C;
    for (auto &kv : l_blocks) {
        const ivec &cols = kv.second;
        const int l = kv.first.second, nm = 2 * l + 1;
        const int ns = static_cast<int>(cols.size()) / nm;
        if (ns <= 1) continue;
        //cols is shell-major with the m components contiguous
        MatrixXd B = MatrixXd::Zero(ns, ns);
        for (int s1 = 0; s1 < ns; s1++)
            for (int s2 = 0; s2 < ns; s2++)
                for (int m = 0; m < nm; m++)
                    B(s1, s2) += Porb(cols[s1 * nm + m], cols[s2 * nm + m]) / nm;
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(B);
        MatrixXd rot(nao, ns * nm);
        for (int k = ns - 1; k >= 0; k--) {
            const int shell = ns - 1 - k;
            for (int m = 0; m < nm; m++) {
                VectorXd v = VectorXd::Zero(nao);
                for (int s = 0; s < ns; s++) v += es.eigenvectors()(s, k) * C.col(cols[s * nm + m]);
                rot.col(shell * nm + m) = v;
            }
        }
        for (int j = 0; j < ns * nm; j++) C.col(cols[j]) = rot.col(j);
    }

    //---------------------------------------------------------------- results
    NAOResult res;
    const MatrixXd Pfin = C.transpose() * SPS * C;
    for (int i = 0; i < nao; i++) orbitals[i].occupation = Pfin(i, i);

    res.C = to_dmatrix(C);
    res.orbitals = orbitals;
    res.atoms.resize(atoms.size());
    for (size_t a = 0; a < atoms.size(); a++) {
        res.atoms[a].index = static_cast<int>(a);
        res.atoms[a].Z = atoms[a].get_charge();
        res.atoms[a].label = atoms[a].get_label();
        res.atoms[a].Z_eff = atoms[a].get_charge() - ecp_electrons[a];
    }
    for (const NAO &o : orbitals) {
        NAOAtom &at = res.atoms[o.atom];
        at.population += o.occupation;
        if (o.type == NAOClass::Core) at.core += o.occupation;
        else if (o.type == NAOClass::Valence) at.valence += o.occupation;
        else at.rydberg += o.occupation;
    }
    for (NAOAtom &at : res.atoms) {
        at.charge = at.Z_eff - at.population;
        res.population += at.population;
        res.core += at.core;
        res.valence += at.valence;
        res.rydberg += at.rydberg;
    }
    return res;
}

namespace
{
    NAOResult analyse(const dMatrix2 &P, const dMatrix2 &S, const std::vector<NAOBasisFunction> &ao,
                      const WFN &wavy)
    {
        const std::vector<atom> ats = wavy.get_atoms();
        ivec ecp(ats.size(), 0);
        if (wavy.get_has_ECPs())
            for (size_t a = 0; a < ats.size(); a++)
                ecp[a] = wavy.get_atom_ECP_electrons(static_cast<int>(a));
        return build_naos(P, S, ao, ats, ecp);
    }
}

dMatrix2 ao_overlap(const WFN &wavy)
{
    Int_Params params(wavy);
    vec S_flat;
    compute2C<Overlap2C_SPH>(params, S_flat);
    const size_t n = static_cast<size_t>(std::llround(std::sqrt(static_cast<double>(S_flat.size()))));
    dMatrix2 S = reshape<dMatrix2>(S_flat, Shape2D(n, n));
    //An ORCA-convention density (gbw, and a molden written from one) carries the opposite sign on the
    //|m| >= 3 components, so the overlap next to it has to take that sign as well - see
    //origin_has_orca_pure_phases.  This is the same correction the FILE47 writer applies; without it
    //Tr(P S) misses up to 0.3 e (SF6) and every NAO population inherits it.  Scanned against the
    //electron count on CuF2_i_func/71/calc.gbw (shells up to i): flipping every |m| >= 3 gives
    //46.99929 of 47, stopping at |m| <= 3 gives 46.99845, flipping nothing 46.95769 - so "all |m| >= 3"
    //it is, and the 7e-4 that remains is a separate high-l matter, identical for the gbw and the molden.
    if (origin_has_orca_pure_phases(wavy.get_origin())) {
        const ivec bas = params.get_bas();
        bvec flip(n, false);
        size_t k = 0;
        for (size_t s = 0; s < params.get_nbas() && k < n; s++) {
            const int l = bas[8 * s + 1];
            for (int m = -l; m <= l && k < n; m++, k++)
                flip[k] = std::abs(m) >= 3;
        }
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                if (flip[i] != flip[j]) S(i, j) = -S(i, j);
    }
    return S;
}

NPAResult natural_population_analysis(const WFN &wavy)
{
    const std::vector<NAOBasisFunction> ao = spherical_ao_map(wavy);
    const dMatrix2 S = ao_overlap(wavy);
    const dMatrix2 P = wavy.get_dm();
    err_checkf(P.extent(0) == ao.size(),
               "NAO/NPA: this wavefunction carries no density matrix over the contracted basis "
               "(a .wfn/.wfx holds primitives only) - use the .gbw, .fchk or .molden of the same run.",
               std::cout);

    NPAResult res;
    const dMatrix2 Pb = wavy.get_dm_beta();
    if (wavy.get_is_unrestricted() && Pb.extent(0) == ao.size()) {
        //different NAOs for different spins, as NBO does by default: the two spin densities are
        //analysed independently and the charges come from the sum of the two populations
        dMatrix2 Pa(ao.size(), ao.size());
        for (size_t i = 0; i < ao.size(); i++)
            for (size_t j = 0; j < ao.size(); j++)
                Pa(i, j) = P(i, j) - Pb(i, j);
        res.alpha = analyse(Pa, S, ao, wavy);
        res.beta = analyse(Pb, S, ao, wavy);
        res.spin_resolved = true;
        res.total = res.alpha;
        for (size_t a = 0; a < res.total.atoms.size(); a++) {
            NAOAtom &t = res.total.atoms[a];
            const NAOAtom &b = res.beta.atoms[a];
            t.population += b.population;
            t.core += b.core;
            t.valence += b.valence;
            t.rydberg += b.rydberg;
            t.charge = t.Z_eff - t.population;
            res.spin_population.push_back(res.alpha.atoms[a].population - b.population);
        }
        res.total.population = res.alpha.population + res.beta.population;
        res.total.core = res.alpha.core + res.beta.core;
        res.total.valence = res.alpha.valence + res.beta.valence;
        res.total.rydberg = res.alpha.rydberg + res.beta.rydberg;
        //the orbital table of a spin-resolved run belongs to the two spin sets
        res.total.orbitals.clear();
        res.total.C = dMatrix2();
    }
    else {
        err_checkf(!wavy.get_is_unrestricted(),
                   "NAO/NPA: this reader did not keep a beta density, so the open-shell case "
                   "cannot be resolved by spin. Read the calculation from a .gbw or .molden.",
                   std::cout);
        res.total = analyse(P, S, ao, wavy);
    }
    return res;
}

namespace
{
    //with_charge = false for one spin on its own, where Z_eff minus that spin's population is not
    //a charge and printing it invites the reader to add the two tables up
    void print_one(const NAOResult &r, const std::string &title, std::ostream &out,
                   const bool with_charge = true)
    {
        using namespace std;
        if (!r.orbitals.empty()) {
            out << "\n NATURAL POPULATIONS:  " << title << "\n\n"
                << "  NAO Atom No lang   Type(AO)    Occupancy\n"
                << " -------------------------------------------------\n";
            int last_atom = -1;
            for (size_t i = 0; i < r.orbitals.size(); i++) {
                const NAO &o = r.orbitals[i];
                if (o.atom != last_atom && last_atom >= 0) out << "\n";
                last_atom = o.atom;
                out << setw(4) << i + 1 << setw(5) << r.atoms[o.atom].label << setw(3) << o.atom + 1
                    << "  " << left << setw(7) << shell_label(o.l, o.m) << right
                    << class_label(o.type) << "(" << setw(2) << o.n
                    << string(1, "spdfghik"[std::min(o.l, 7)]) << ")" << setw(12) << fixed
                    << setprecision(5) << o.occupation << "\n";
            }
        }
        out << "\n Summary of Natural Population Analysis:\n\n"
            << "                                     Natural Population\n"
            << "             Natural    ---------------------------------------------\n"
            << "  Atom No    Charge        Core      Valence    Rydberg      Total\n"
            << " --------------------------------------------------------------------\n";
        for (const NAOAtom &a : r.atoms) {
            out << setw(5) << r.atoms[a.index].label << setw(3) << a.index + 1 << fixed
                << setprecision(5);
            if (with_charge) out << setw(11) << a.charge;
            else out << setw(11) << "-";
            out << setw(13) << a.core << setw(12) << a.valence << setw(11) << a.rydberg << setw(12)
                << a.population << "\n";
        }
        out << " ====================================================================\n"
            << " * Total * " << fixed << setprecision(5);
        if (with_charge)
            out << setw(9) << std::accumulate(r.atoms.begin(), r.atoms.end(), 0.0,
                                              [](double s, const NAOAtom &a) { return s + a.charge; });
        else
            out << setw(9) << "-";
        out << setw(13) << r.core << setw(12) << r.valence << setw(11) << r.rydberg << setw(12)
            << r.population << "\n";

        if (r.orbitals.empty()) return;
        out << "\n    Atom No         Natural Electron Configuration\n"
            << " ----------------------------------------------------------------------------\n";
        for (const NAOAtom &a : r.atoms) {
            std::map<std::pair<int, int>, double> shells;  //(n, l) -> occupancy
            bool has_core = false;
            for (const NAO &o : r.orbitals) {
                if (o.atom != a.index) continue;
                if (o.type == NAOClass::Core) { has_core = true; continue; }
                shells[{ o.n, o.l }] += o.occupation;
            }
            out << setw(7) << r.atoms[a.index].label << setw(3) << a.index + 1 << "      "
                << (has_core ? "[core]" : "      ");
            for (const auto &kv : shells) {
                if (kv.second < 0.005) continue;
                out << kv.first.first << "spdfghik"[std::min(kv.first.second, 7)] << "("
                    << fixed << setprecision(2) << setw(5) << kv.second << ")";
            }
            out << "\n";
        }
    }
}

void print_npa(const NPAResult &result, std::ostream &out)
{
    citations::cite(citations::Method::NAONPA, out);
    if (result.spin_resolved) {
        print_one(result.alpha, "alpha spin natural atomic orbital occupancies", out, false);
        print_one(result.beta, "beta spin natural atomic orbital occupancies", out, false);
        print_one(result.total, "both spins", out);
        out << "\n\n  Atom No     Natural Charge      Spin Population (alpha - beta)\n"
            << " ---------------------------------------------------------------\n";
        for (size_t a = 0; a < result.total.atoms.size(); a++)
            out << std::setw(7) << result.total.atoms[a].label << std::setw(3) << a + 1 << std::fixed
                << std::setprecision(5) << std::setw(16) << result.total.atoms[a].charge
                << std::setw(22) << result.spin_population[a] << "\n";
    }
    else {
        print_one(result.total, "Natural atomic orbital occupancies", out);
    }
    out << std::endl;
}
