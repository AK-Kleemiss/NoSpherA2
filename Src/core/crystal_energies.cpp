#include "pch.h"
#include "crystal_energies.h"
#include "SALTED_predictor.h"
#include "basis_set.h"
#include "nos_math.h"
#include "npy.h"

namespace crystal_energies {
    vec2 real_sh_rotation(const int l, const vec2& R)
    {
        const int n = 2 * l + 1, K = 4 * n + 8;
        vec2 D(n, vec(n, 0.0));
        if (l == 0) { D[0][0] = 1.0; return D; }
        //Y_lm(R^-1 u) = sum_k D[k][m] Y_lk(u), least squares over a Fibonacci sphere, R^-1 = R^T
        vec2 Y(K, vec(n)), Yr(K, vec(n)), G(n, vec(n, 0.0));
        for (int k = 0; k < K; k++) {
            const double z = 1.0 - 2.0 * (k + 0.5) / K, r = std::sqrt(1.0 - z * z), phi = k * 2.399963229728653;
            const double u[3] = { r * std::cos(phi), r * std::sin(phi), z };
            double v[3];
            for (int x = 0; x < 3; x++) v[x] = R[0][x] * u[0] + R[1][x] * u[1] + R[2][x] * u[2];
            for (int m = 0; m < n; m++) Y[k][m] = constants::spherical_harmonic(l, m - l, u), Yr[k][m] = constants::spherical_harmonic(l, m - l, v);
        }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < K; k++) G[i][j] += Y[k][i] * Y[k][j];
        for (int m = 0; m < n; m++) {
            vec b(n, 0.0);
            for (int i = 0; i < n; i++)
                for (int k = 0; k < K; k++) b[i] += Y[k][i] * Yr[k][m];
            solve_linear_system(G, b);
            for (int i = 0; i < n; i++) D[i][m] = b[i];
        }
        return D;
    }

    void transform(WFN& aux, vec& coef, const vec2& R, const vec& t)
    {
        std::vector<atom> atoms = aux.get_atoms();
        std::map<int, vec2> D;
        int offset = 0;
        for (int a = 0; a < (int)atoms.size(); a++) {
            double x[3];
            for (int i = 0; i < 3; i++) x[i] = t[i] + R[i][0] * atoms[a].get_coordinate(0) + R[i][1] * atoms[a].get_coordinate(1) + R[i][2] * atoms[a].get_coordinate(2);
            for (int i = 0; i < 3; i++) atoms[a].set_coordinate(i, x[i]);
            int prim = 0;
            for (int shell = 0; shell < (int)atoms[a].get_shellcount_size(); shell++) {
                const int l = atoms[a].get_basis_set_entry(prim).get_type(), n = 2 * l + 1;
                if (D.find(l) == D.end()) D[l] = real_sh_rotation(l, R);
                err_checkf(offset + n <= (int)coef.size(), "More auxiliary functions than coefficients", std::cout);
                const vec c(coef.begin() + offset, coef.begin() + offset + n);
                for (int i = 0; i < n; i++) {
                    coef[offset + i] = 0.0;
                    for (int j = 0; j < n; j++) coef[offset + i] += D[l][i][j] * c[j];
                }
                offset += n, prim += atoms[a].get_shellcount(shell);
            }
        }
        err_checkf(offset == (int)coef.size(), "Coefficient count does not match the auxiliary basis", std::cout);
        aux.set_atoms(atoms);
    }

    std::vector<symop> symops(const cell& c)
    {
        const std::vector<ivec2> sym = c.get_sym();
        const std::vector<vec> trans = c.get_trans();
        std::vector<symop> ops(trans[0].size());
        for (int s = 0; s < (int)ops.size(); s++) {
            ops[s].rot.assign(3, ivec(3));
            ops[s].trans.resize(3);
            for (int i = 0; i < 3; i++) {
                ops[s].trans[i] = trans[i][s];
                for (int j = 0; j < 3; j++) ops[s].rot[i][j] = sym[j][i][s];
            }
        }
        return ops;
    }

    std::string symop_string(const symop& op, const ivec& n)
    {
        const char axis[3] = { 'x', 'y', 'z' };
        std::string s;
        for (int i = 0; i < 3; i++) {
            std::string comp;
            for (int j = 0; j < 3; j++) {
                if (op.rot[i][j] == 0) continue;
                comp += op.rot[i][j] < 0 ? "-" : comp.empty() ? "" : "+";
                if (std::abs(op.rot[i][j]) != 1) comp += std::to_string(std::abs(op.rot[i][j])) + "*";
                comp += axis[j];
            }
            const double t = op.trans[i] + n[i];
            if (std::abs(t) > 1e-6) {
                const int twelfths = (int)std::lround(t * 12.0);
                if (std::abs(twelfths - t * 12.0) < 1e-6) {
                    const int g = std::gcd(std::abs(twelfths), 12);
                    comp += (twelfths < 0 ? "-" : "+") + std::to_string(std::abs(twelfths) / g) + (g == 12 ? "" : "/" + std::to_string(12 / g));
                }
                else comp += (t < 0 ? "-" : "+") + std::to_string(std::abs(t));
            }
            s += (i ? "," : "") + comp;
        }
        return s;
    }

    vec2 cell_matrix(const cell& c)
    {
        vec2 M(3, vec(3));
        for (int j = 0; j < 3; j++) {
            const vec v = c.get_coords_cartesian(j == 0, j == 1, j == 2, true);
            for (int i = 0; i < 3; i++) M[i][j] = v[i];
        }
        return M;
    }

    vec2 inverse3(const vec2& M)
    {
        vec2 I(3, vec(3));
        const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
        err_checkf(std::abs(det) > 1e-12, "Singular cell matrix", std::cout);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                const int r1 = (j + 1) % 3, r2 = (j + 2) % 3, c1 = (i + 1) % 3, c2 = (i + 2) % 3;
                I[i][j] = (M[r1][c1] * M[r2][c2] - M[r1][c2] * M[r2][c1]) / det;
            }
        return I;
    }

    vec fitted_coefficients(const WFN& wavy, const std::filesystem::path& coef_file, WFN& aux, options& opt)
    {
        vec coef;
        if (!coef_file.empty()) {
            err_checkf(!opt.aux_basis.empty(), "No auxiliary basis set specified! Use -ri_fit BEFORE the interaction energy flag", std::cout);
            aux = generate_aux_wfn(wavy, opt.aux_basis);
            std::vector<unsigned long> shape; bool fortran_order;
            npy::LoadArrayFromNumpy(coef_file.string(), shape, fortran_order, coef);
        }
        else if (opt.SALTED) {
            err_checkf(!opt.salted_model_dir.empty(), "No SALTED model directory specified! Use -SALTED <model-dir> BEFORE the interaction energy flag", std::cout);
            SALTEDPredictor SP(wavy, opt);
            if (!SP.basis_set_loaded()) load_basis_into_WFN(SP.wavy, BasisSetLibrary::get_basis_set(SP.get_dfbasis_name()));
            coef = SP.gen_SALTED_densities();
            err_checkf(SP.wavy.get_ncen() == wavy.get_ncen(), "The SALTED model does not cover every atom of " + wavy.get_path().string(), std::cout);
            aux = SP.wavy;
            aux.set_origin(e_origin::NOT_YET_DEFINED);
        }
        else {
            err_checkf(!opt.aux_basis.empty(), "No auxiliary basis set specified! Use -ri_fit <basis> or -SALTED <model-dir> BEFORE the interaction energy flag", std::cout);
            err_checkf(wavy.get_nmo() > 0, "Only a wavefunction can be fitted; " + wavy.get_path().string() + " needs coefficients or -SALTED <model-dir>", std::cout);
            aux = generate_aux_wfn(wavy, opt.aux_basis);
            coef = DensityFitting::density_fit(wavy, aux, DensityFitting::config_from_options(opt));
        }
        return coef;
    }

    namespace {
        std::string key(const vec& x)
        {
            char s[64];
            snprintf(s, sizeof(s), "%.3f %.3f %.3f", x[0], x[1], x[2]);
            return s;
        }
        vec mul(const vec2& M, const vec& v)
        {
            vec r(3, 0.0);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) r[i] += M[i][j] * v[j];
            return r;
        }
        vec2 mul(const vec2& A, const vec2& B)
        {
            vec2 C(3, vec(3, 0.0));
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++) C[i][j] += A[i][k] * B[k][j];
            return C;
        }
    }

    std::vector<pair> contacts(const cell& c, const std::vector<WFN>& mol, const double cutoff)
    {
        const std::vector<symop> ops = symops(c);
        const vec2 M = cell_matrix(c), Mi = inverse3(M);
        const double cut = constants::ang2bohr(cutoff), cut2 = cut * cut;
        const int nm = (int)mol.size();
        std::vector<vec> cen(nm, vec(3, 0.0));
        vec rad(nm, 0.0), spacing(3);
        for (int i = 0; i < nm; i++) {
            const int na = mol[i].get_ncen();
            err_checkf(na > 0, "A molecule without atoms", std::cout);
            for (int a = 0; a < na; a++)
                for (int x = 0; x < 3; x++) cen[i][x] += mol[i].get_atom_coordinate(a, x) / na;
            for (int a = 0; a < na; a++) {
                double d2 = 0.0;
                for (int x = 0; x < 3; x++) d2 += std::pow(mol[i].get_atom_coordinate(a, x) - cen[i][x], 2);
                rad[i] = std::max(rad[i], std::sqrt(d2));
            }
        }
        //interplanar spacing along axis x is 1/|row x of M^-1|
        for (int x = 0; x < 3; x++) spacing[x] = 1.0 / std::sqrt(Mi[x][0] * Mi[x][0] + Mi[x][1] * Mi[x][1] + Mi[x][2] * Mi[x][2]);
        std::vector<pair> out;
        for (int i = 0; i < nm; i++) {
            const int na = mol[i].get_ncen();
            vec2 PA(na, vec(3));
            for (int a = 0; a < na; a++)
                for (int x = 0; x < 3; x++) PA[a][x] = mol[i].get_atom_coordinate(a, x);
            for (int j = i; j < nm; j++) {
                std::set<std::string> seen;
                const vec cA = mul(Mi, cen[i]);
                const int nb = mol[j].get_ncen();
                for (int s = 0; s < (int)ops.size(); s++) {
                    vec2 Rf(3, vec(3));
                    for (int x = 0; x < 3; x++)
                        for (int y = 0; y < 3; y++) Rf[x][y] = ops[s].rot[x][y];
                    const vec2 R = mul(M, mul(Rf, Mi));
                    //rotated atoms and centroid of B, the lattice translation is added per image
                    vec2 XB(nb, vec(3));
                    for (int b = 0; b < nb; b++) {
                        vec p(3);
                        for (int x = 0; x < 3; x++) p[x] = mol[j].get_atom_coordinate(b, x);
                        XB[b] = mul(R, p);
                    }
                    const vec cR = mul(R, cen[j]);
                    vec cB0 = mul(Rf, mul(Mi, cen[j]));
                    ivec base(3), range(3), n(3);
                    for (int x = 0; x < 3; x++) cB0[x] += ops[s].trans[x], base[x] = (int)std::lround(cA[x] - cB0[x]), range[x] = (int)std::ceil((rad[i] + rad[j] + cut) / spacing[x]) + 1;
                    for (n[0] = base[0] - range[0]; n[0] <= base[0] + range[0]; n[0]++)
                        for (n[1] = base[1] - range[1]; n[1] <= base[1] + range[1]; n[1]++)
                            for (n[2] = base[2] - range[2]; n[2] <= base[2] + range[2]; n[2]++) {
                                vec t(3), cX(3);
                                for (int x = 0; x < 3; x++) t[x] = ops[s].trans[x] + n[x];
                                const vec tc = mul(M, t);
                                for (int x = 0; x < 3; x++) cX[x] = cR[x] + tc[x];
                                double dc2 = 0.0;
                                for (int x = 0; x < 3; x++) dc2 += std::pow(cX[x] - cen[i][x], 2);
                                //no atom pair can be closer than the centroids minus both extents
                                if (std::sqrt(dc2) - rad[i] - rad[j] >= cut) continue;
                                double dmin2 = DBL_MAX;
                                for (int b = 0; b < nb; b++)
                                    for (int a = 0; a < na; a++) {
                                        double d2 = 0.0;
                                        for (int x = 0; x < 3; x++) d2 += std::pow(XB[b][x] + tc[x] - PA[a][x], 2);
                                        dmin2 = std::min(dmin2, d2);
                                    }
                                if (dmin2 >= cut2) continue;
                                if (i == j && dc2 < 1e-6) continue;
                                if (!seen.insert(key(cX)).second) continue;
                                if (i == j) {
                                    //the same pair seen from B: A is the image of B under the inverse operation
                                    vec d(3), cY(3, 0.0);
                                    for (int x = 0; x < 3; x++) d[x] = cen[i][x] - tc[x];
                                    for (int x = 0; x < 3; x++)
                                        for (int y = 0; y < 3; y++) cY[x] += R[y][x] * d[y];
                                    seen.insert(key(cY));
                                }
                                pair p;
                                p.A = i, p.B = j, p.op = s, p.n = n, p.distance = constants::bohr2ang(std::sqrt(dc2)), p.R = R, p.t = tc;
                                p.centroid_A = cen[i], p.centroid_B = cX;
                                for (int x = 0; x < 3; x++) p.centroid_A[x] = constants::bohr2ang(p.centroid_A[x]), p.centroid_B[x] = constants::bohr2ang(p.centroid_B[x]);
                                p.symop = symop_string(ops[s], n);
                                out.push_back(p);
                            }
                }
            }
        }
        std::sort(out.begin(), out.end(), [](const pair& a, const pair& b) { return a.A != b.A ? a.A < b.A : a.B != b.B ? a.B < b.B : a.distance < b.distance; });
        return out;
    }

    job read_job(const std::filesystem::path& file)
    {
        err_checkf(std::filesystem::exists(file), "The job file does not exist: " + file.string(), std::cout);
        const std::filesystem::path dir = file.parent_path();
        auto resolve = [&dir](const std::string& p) { const std::filesystem::path q(p); return q.is_absolute() ? q : dir / q; };
        job J;
        J.output = resolve("interaction_energies.txt");
        std::ifstream in(file);
        std::string line;
        while (getline_universal(in, line)) {
            const std::string entry = trim(line);
            if (entry.empty() || entry[0] == '#') continue;
            std::istringstream ss(entry);
            std::string key, a, b;
            ss >> key >> a;
            std::getline(ss, b);
            b = trim(b);
            if (key == "cif") J.cif = resolve(a);
            else if (key == "cutoff") J.cutoff = std::stod(a);
            else if (key == "output") J.output = resolve(a);
            else if (key == "molecule") J.structures.push_back(resolve(a)), J.coefficients.push_back(b.empty() ? std::filesystem::path() : resolve(b));
            else err_checkf(false, "Unknown keyword \"" + key + "\" in " + file.string(), std::cout);
        }
        err_checkf(!J.cif.empty(), "The job file " + file.string() + " names no cif", std::cout);
        err_checkf(!J.structures.empty(), "The job file " + file.string() + " names no molecule", std::cout);
        return J;
    }

    void write_table(const std::vector<pair>& pairs, const pathvec& structures, std::ostream& out)
    {
        const double kJ = constants::kcal_mol_per_hartree * 4.184;
        out << "# NoSpherA2 interaction energies in kJ/mol; B = symop(B) + n, centroids and R_AB in Angstrom\n";
        for (int i = 0; i < (int)structures.size(); i++) out << "# molecule " << i << " " << structures[i].string() << "\n";
        out << "# A B symop n1 n2 n3 R_AB xA yA zA xB yB zB E_elst E_pol E_disp E_rep E_total\n";
        for (const pair& p : pairs) {
            out << p.A << " " << p.B << " " << p.symop << " " << p.n[0] << " " << p.n[1] << " " << p.n[2] << std::fixed << std::setprecision(4) << " " << p.distance;
            for (int x = 0; x < 3; x++) out << " " << p.centroid_A[x];
            for (int x = 0; x < 3; x++) out << " " << p.centroid_B[x];
            out << std::setprecision(3) << " " << p.E.electrostatic() * kJ << " " << (p.E.pol_A + p.E.pol_B) * kJ << " " << p.E.disp * kJ << " " << p.E.rep * kJ << " " << p.E.total() * kJ << "\n";
        }
        out << std::flush;
    }

    int run(options& opt)
    {
        const job J = read_job(opt.interaction_energies_job);
        cell c(J.cif, std::cout, opt.debug, true);
        std::cout << std::endl;
        std::vector<WFN> aux;
        std::vector<vec> coef;
        for (int i = 0; i < (int)J.structures.size(); i++) {
            const WFN wavy(J.structures[i]);
            WFN a;
            coef.push_back(fitted_coefficients(wavy, J.coefficients[i], a, opt));
            aux.push_back(a);
        }
        std::vector<pair> pairs = contacts(c, aux, J.cutoff);
        std::cout << pairs.size() << " molecule pairs within " << J.cutoff << " A" << std::endl;
        const double kJ = constants::kcal_mol_per_hartree * 4.184;
        for (int p = 0; p < (int)pairs.size(); p++) {
            WFN B = aux[pairs[p].B];
            vec cB = coef[pairs[p].B];
            transform(B, cB, pairs[p].R, pairs[p].t);
            pairs[p].E = DensityFitting::interaction_energy(coef[pairs[p].A], aux[pairs[p].A], cB, B, opt.repulsion_overlap, opt.repulsion_exchange);
            std::cout << "PAIR " << p + 1 << "/" << pairs.size() << " " << pairs[p].A << " " << pairs[p].B << " " << pairs[p].symop << " n=" << pairs[p].n[0] << "," << pairs[p].n[1] << "," << pairs[p].n[2]
                      << std::fixed << std::setprecision(3) << " R=" << pairs[p].distance << " total=" << pairs[p].E.total() * kJ << " kJ/mol" << std::endl;
        }
        std::ofstream out(J.output);
        err_checkf(out.good(), "Could not write " + J.output.string(), std::cout);
        write_table(pairs, J.structures, out);
        write_table(pairs, J.structures, std::cout);
        std::cout << "Wrote " << J.output.string() << std::endl;
        return 0;
    }
}
