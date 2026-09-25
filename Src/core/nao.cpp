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
    //Diagnostic knobs for the NAO construction, off unless the environment sets them.  They
    //exist to run the two arms of one experiment - see the comments at steps 3 and 4 - and are
    //deliberately not command-line options: build_naos takes no options struct and this is a
    //measurement, not a feature.
    bool nao_env(const char *name)
    {
        const char *v = std::getenv(name);
        return v && *v && *v != '0';
    }

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
    //Both of those numbers are NPA charges, and the class-split arm proved an NPA comparison blind
    //to the class partition: 0.695 e of benzene's Rydberg set moved while all 111 charges agreed to
    //1e-10.  NAO_PRENAO_NET=1 re-runs the net variant so the Rydberg metric can judge it instead.
    //It is worth re-running because of what fixes the class totals after step 3: the Rydberg set is
    //then the S-orthogonal complement of the span of the natural minimal pre-NAOs, so its total
    //population depends only on that span - on these eigenvectors and on which shells step 2 calls
    //minimal - and on nothing inside step 3.  NAO_OWSO_OFF=1 below is the check of that claim.
    //Measured (job 586936, all 25 wavefunctions): the net variant is not a near miss, it is a
    //different answer.  The Rydberg total after step 3 goes from 10.40557 to 194.78163 e over the
    //22 closed shells and the final one from 4.00571 to 56.45473 against gennbo's 2.74986, every
    //single molecule between 8x and 72x, and the pre-NAO occupancies then sum to 0.8 N instead of
    //1.6 N.  Gross is confirmed by the Rydberg metric and no longer rests on the charge argument.
    const bool prenao_net = nao_env("NAO_PRENAO_NET");
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
                    Pb(s1, s2) += (prenao_net ? P(g.idx[s1][m], g.idx[s2][m])
                                              : SPS(g.idx[s1][m], g.idx[s2][m])) / nm;
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
        //NAO_OWSO_OFF=1 drops the occupancy weighting and leaves the plain Loewdin below as the
        //whole of the within-class orthogonalisation.  That is not a candidate - it throws away the
        //point of the OWSO - it is the test of the invariance claimed above: the weighting picks
        //different vectors inside the class but cannot change the class's span, so every class
        //TOTAL after step 3 must come out unchanged while the final ones move.  If a step-3 total
        //does move, the span argument is wrong and the search cannot be narrowed by it.
        //Measured (job 586936): confirmed.  Every molecule's Rydberg total after step 3 is
        //unchanged to all five printed decimals (ethane 0.43210, benzene 1.12682, sf6 0.71526)
        //while the final ones move everywhere.  So no choice inside this loop can change a class
        //TOTAL at step 3 - but the distribution over (atom, l) blocks is NOT invariant, and that is
        //what the final answer sees: the no-valence-block excess went 1.39828 -> 1.14700 e and the
        //final excess 1.25585 -> 1.12478 e.  The weighting is therefore a lever on the final
        //numbers, just not on the step-3 totals, and "only step 3 crosses l" was the wrong
        //localisation.  It is still not a candidate, because the pre-registered gate in
        //tests/nbo_reference_v2/step3_stages.py fired on exactly this arm: the mean improves while
        //pf5 goes 0.93 -> 1.12, so2 0.94 -> 1.05 and sf6 0.98 -> 1.19 x gennbo, sf6 worse by
        //0.077 e on its own.  Flattening the three molecules that were right to improve the average
        //is the fudge-factor signature the gate was written to catch.
        if (!nao_env("NAO_OWSO_OFF"))
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
    //components of a shell degenerate.  Valence and Rydberg go in one block, and NAO_CLASS_SPLIT=1
    //runs the alternative so the claim can be checked instead of believed.  Measured over the same
    //22 wavefunctions against gennbo 7's own NAO tables: splitting the classes raises the mean
    //intra-atomic valence->Rydberg leak from 0.0123 to 0.0690 e over 111 atoms (benzene's Rydberg
    //set 0.432 -> 1.127 e), and it does so in all 22 - including pf5, so2 and sf6, the only three
    //whose leak has the opposite sign in the one-block form.  One block is the better of the two,
    //so whatever is left of the leak is upstream of here: step 3, or which n,l count as valence.
    //That last sentence was too strong and job 586936 corrected it.  The class totals step 4
    //inherits are fixed by the span of the natural minimal pre-NAOs alone (NAO_OWSO_OFF leaves them
    //identical to five decimals), and the only other input to that span, the net-density pre-NAO
    //variant, is catastrophic.  So there is nothing upstream left to change: the final numbers move
    //only through which vectors step 3 hands to each (atom, l) block and what this step then does
    //with them.  That last clause read "and NBO 7 prints no intermediate table, so neither side of
    //that can be arbitrated", which over-generalised from occupancies to everything: NBO 7 prints no
    //intermediate OCCUPANCY table, but `$NBO AONAO=W $END` writes its AO -> NAO matrix to lfn 33
    //(verified on this install, job 588007).  The vectors themselves are arbitrable; see NAO_DUMP_C.
    //The mixing is intra-atomic and unitary, so it moves no charge between atoms - and that is not
    //a reassurance, it is a warning.  The two arms differ by 0.695 e in benzene's Rydberg
    //population and by 1e-10 in every one of 111 NPA charges, so an NPA comparison cannot see this
    //error at all and agreement there says nothing about the class partition.
    const MatrixXd Porb = C.transpose() * SPS * C;
    //NAO_DUMP_STEP3: the occupancies step 4 inherits.  If the valence deficit against NBO 7 is
    //already visible here, step 4 is not the place to look for it.
    if (nao_env("NAO_DUMP_STEP3")) {
        std::cout << "STEP3 atom l shell class occ_per_component pre_occ_per_component" << std::endl;
        for (const auto &kv : l_blocks) {
            const int nm = 2 * kv.first.second + 1;
            const ivec &cols = kv.second;
            for (size_t sh = 0; sh * nm < cols.size(); sh++) {
                double occ = 0.0, pre = 0.0;
                for (int m = 0; m < nm; m++) {
                    occ += Porb(cols[sh * nm + m], cols[sh * nm + m]) / nm;
                    pre += pre_occ(cols[sh * nm + m]) / nm;
                }
                std::cout << "STEP3 " << kv.first.first << " " << kv.first.second << " " << sh
                          << " " << static_cast<int>(orbitals[cols[sh * nm]].type)
                          << " " << std::setprecision(8) << std::fixed << occ
                          << " " << pre << std::endl;
            }
        }
    }
    //NAO_CLASS_SPLIT runs the rejected variant: the shells of one (atom, l) are re-diagonalised
    //within their own class instead of all together.  Note what the one-block form is: the top
    //eigenvalue of the m-averaged block is an upper bound on any single shell's own diagonal, so
    //the valence shell can only come out with MORE population from one block than from a split
    //one.  Splitting therefore cannot reduce a Rydberg excess - it has to increase it - which is
    //what makes this arm a real test rather than a search for a better number.
    const bool class_split = nao_env("NAO_CLASS_SPLIT");
    std::vector<ivec> blocks;
    for (auto &kv : l_blocks) {
        const int l_of_block = kv.first.second, nm_of_block = 2 * l_of_block + 1;
        if (!class_split) { blocks.push_back(kv.second); continue; }
        ivec per_class[3];
        const ivec &all = kv.second;
        for (size_t j = 0; j * nm_of_block < all.size(); j++) {
            const int cls = static_cast<int>(orbitals[all[j * nm_of_block]].type);
            for (int m = 0; m < nm_of_block; m++)
                per_class[cls].push_back(all[j * nm_of_block + m]);
        }
        for (int cls = 0; cls < 3; cls++)
            if (!per_class[cls].empty()) blocks.push_back(per_class[cls]);
    }
    for (const ivec &block_cols : blocks) {
        const ivec &cols = block_cols;
        const int l = orbitals[cols[0]].l, nm = 2 * l + 1;
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

    //NAO_DUMP_C: the AO -> NAO transformation itself, one line per NAO.  This is the only
    //arbitrated INTERMEDIATE that exists.  NBO 7 prints no pre-NAO and no step-3 occupancy table,
    //but `$NBO AONAO=W $END` writes its own AO -> NAO matrix to lfn 33 at nine decimals - verified
    //on this install, job 588007, where AONAO=W48 is refused because 48 is reserved - and both
    //sides read the SAME .47, so the AO order, S and P are literally the same arrays.  A final
    //occupancy table can only say that a block came out with the wrong population; the overlap
    //c_native^T S c_nbo says which orbital has the wrong shape, which is what a fix needs.
    //
    //What it said (jobs 589060 + 589634, 8 molecules, 0 VOID; every gate passed, completeness to
    //2.74e-09).  Mean m-averaged mixing defect per shell, rank-paired inside each (atom, l) block:
    //Core 0.00000 (34 shells), Valence 0.00783 (72), Rydberg 0.37253 (273).  The cores being EXACTLY
    //right exonerates the AO read, S, P and the core partition in one number.  Weighted by occupancy
    //the ranking INVERTS, and that is the finding: of 1.33759 e of mis-shaped density, the valence
    //shells carry 0.94266 e - 0.51486 e of it outside their own (atom, l) block - against the whole
    //Rydberg set's 0.39467 e (0.19472 e outside), whose per-shell error is 48x larger.  None of the
    //valence figure is an ordering artefact: the ordering-insensitive arm agrees to five decimals.  So the per-shell error lives in the Rydberg
    //construction while the charge that moves is valence, which is what d(Val) = -d(Ryd) looks like
    //one level down - and only the electron-weighted number is commensurable with the NPA failure.
    //Do not quote the per-shell mean alone: it ranks a badly-shaped empty Rydberg shell above a
    //nearly-right doubly-occupied valence one.  Tables in tests/nbo_reference_v2/README.md.
    if (nao_env("NAO_DUMP_C")) {
        std::cout << "NAOC index atom l m shell class occ coefficients[" << nao << "]" << std::endl;
        for (int i = 0; i < nao; i++) {
            const NAO &o = orbitals[i];
            std::cout << "NAOC " << i << " " << o.atom << " " << o.l << " " << o.m << " " << o.shell
                      << " " << static_cast<int>(o.type) << " " << std::setprecision(10)
                      << std::fixed << o.occupation;
            for (int k = 0; k < nao; k++) std::cout << " " << C(k, i);
            std::cout << std::endl;
        }
    }

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
    //ORCA stores the |m| >= 3 components - f(+-3), g(+-3), g(+-4) - with the sign opposite to
    //libcint's, and the gbw reader keeps its convention in the density, so the overlap next to that
    //density has to take ORCA's sign as well.  This is the same correction the FILE47 writer
    //applies; without it Tr(P S) misses up to 0.3 e (SF6) and every NAO population inherits it.
    if (wavy.get_origin() == e_origin::gbw) {
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
