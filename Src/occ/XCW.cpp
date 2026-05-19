#include "XCW.h"
#include "../convenience.h"
#include "../scattering_factors.h"
#include "../npy.h"
#include "../integration_params.h"
#include "../libCintMain.h"
#include "../nos_math.h"


//#include <occ/io/xyz.h>
//#include <occ/core/molecule.h>
//#include <occ/qm/hf.h>
//#include <occ/qm/scf.h>
//#include <occ/qm/scf_impl.h>

void XCW::U_frac2U_rec(){
    // For now only second order ADPs are supported
    vec norm(3);
    vec2 rec_matrix(3, vec(3));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rec_matrix[i][j] = unit_cell.get_rcm(i, j);
        }
    }
    const double scale = constants::ang2bohr(1) / constants::TWO_PI;
    std::transform(rec_matrix.begin(), rec_matrix.end(), rec_matrix.begin(), [scale](std::vector<double>& vec) {
        std::transform(vec.begin(), vec.end(), vec.begin(), [scale](double x) { return x * scale; });
        return vec;});
    norm[0] = std::sqrt(rec_matrix[0][0] * rec_matrix[0][0] + rec_matrix[1][0] * rec_matrix[1][0] + rec_matrix[2][0] * rec_matrix[2][0]);
    norm[1] = std::sqrt(rec_matrix[0][1] * rec_matrix[0][1] + rec_matrix[1][1] * rec_matrix[1][1] + rec_matrix[2][1] * rec_matrix[2][1]);
    norm[2] = std::sqrt(rec_matrix[0][2] * rec_matrix[0][2] + rec_matrix[1][2] * rec_matrix[1][2] + rec_matrix[2][2] * rec_matrix[2][2]);
    vec2 transform(3, vec(6));
    transform[0][0] = norm[0] * norm[0];
    transform[0][1] = norm[1] * norm[1];
    transform[0][2] = norm[2] * norm[2];
    transform[0][3] = norm[0] * norm[1];
    transform[0][4] = norm[0] * norm[2];
    transform[0][5] = norm[1] * norm[2];
    for (int a = 0; a < ncen; a++) {
        if (wave.get_atom(a).get_ADPs()[0].size() > 0) {
            vec2 ADPs = wave.get_atom(a).get_ADPs();
            ADPs[0][0] *= transform[0][0];
            ADPs[0][1] *= transform[0][1];
            ADPs[0][2] *= transform[0][2];
            ADPs[0][3] *= transform[0][3];
            ADPs[0][4] *= transform[0][4];
            ADPs[0][5] *= transform[0][5];
            wave.set_atom_ADPs(a, ADPs);
        }
    }
}

void XCW::U_rec2U_cart(){
	// For now only second order ADPs are supported
    for (int a = 0; a < ncen; a++) {
        vec2 ADPs = wave.get_atom(a).get_ADPs();
        if (wave.get_atom(a).get_ADPs()[0].size() > 0) {
            vec2 ADP_matrix(3, vec(3));
            ADP_matrix[0][0] = ADPs[0][0];
            ADP_matrix[0][1] = ADPs[0][3];
            ADP_matrix[0][2] = ADPs[0][4];
            ADP_matrix[1][0] = ADPs[0][3];
            ADP_matrix[1][1] = ADPs[0][1];
            ADP_matrix[1][2] = ADPs[0][5];
            ADP_matrix[2][0] = ADPs[0][4];
            ADP_matrix[2][1] = ADPs[0][5];
            ADP_matrix[2][2] = ADPs[0][2];
            vec2 cart_matrix(3, vec(3));
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    cart_matrix[i][j] = unit_cell.get_cm(i, j);
                }
            }
            const double scale = constants::bohr2ang(1);
            std::transform(cart_matrix.begin(), cart_matrix.end(), cart_matrix.begin(), [scale](std::vector<double>& vec) {
                std::transform(vec.begin(), vec.end(), vec.begin(), [scale](double x) { return x * scale; });
                return vec;});
            ADP_matrix = self_dot(self_dot(cart_matrix, ADP_matrix, true, false), cart_matrix, false, false);
            ADPs[0][0] = ADP_matrix[0][0];
            ADPs[0][1] = ADP_matrix[1][1];
            ADPs[0][2] = ADP_matrix[2][2];
            ADPs[0][3] = ADP_matrix[1][0];
            ADPs[0][4] = ADP_matrix[0][2];
            ADPs[0][5] = ADP_matrix[1][2];
            wave.set_atom_ADPs(a, ADPs);
        }
    }
}

void XCW::eval_DW() {
    //Setup
    //Converts angstrom to bohr OR MORE IMPORTANTLY reciprocal bohr to reciprocal angstrom
    const double angstrom2bohr = constants::ang2bohr(1);
    const double scale = angstrom2bohr / constants::TWO_PI;    
    std::vector<int> level;
    level.reserve(ncen);
    //Figure out which level of anisotropic displacements parameters are avaialable
    for (int a = 0; a < ncen; a++) {
		vec2 ADPs = wave.get_atom(a).get_ADPs();
        if (ADPs[2].size() != 0) {
            level.emplace_back(3);
        }
        else if (ADPs[1].size() != 0) {
            level.emplace_back(2);
        }
        else if (ADPs[0].size() != 0) {
            level.emplace_back(1);
        }
        else {
            level.emplace_back(0);
        }
    }
    // Convert ADPs from fractional to Cartesian coordinates
    U_frac2U_rec();
    U_rec2U_cart();
    // Save multiplicative factors for each level
    // const vec mult3 = {1,1,1,3,3,3,3,3,3,6};
    // const vec mult4 = {1,1,1,4,4,4,4,4,4,6,6,6,12,12,12};
    //Eigen::VectorXd mult3(10);
    //mult3 << 1,1,1,3,3,3,3,3,3,6;
    //Eigen::VectorXd mult4(15);
    //mult4 << 1,1,1,4,4,4,4,4,4,6,6,6,12,12,12;
    // Precompute monomials
    // Initialize q for future use
    vec2 q(nr, vec(3));
    //vec qx(nr);
    //vec qy(nr);
    //vec qz(nr);
    for (int h = 0; h < nr; h++) {
        q[h][0] = k_pt[0][h];
        q[h][1] = k_pt[1][h];
        q[h][2] = k_pt[2][h];
        //qx[h] = k_pt[0][h];
        //qy[h] = k_pt[0][h];
        //qz[h] = k_pt[0][h];
    }
    std::transform(q.begin(), q.end(), q.begin(), [scale](std::vector<double>& vec) {
        std::transform(vec.begin(), vec.end(), vec.begin(), [scale](double x) { return x * scale; });
        return vec;});
    // I need to figure out how to get q or if I have to implement it
    //Eigen::MatrixXd m3(nr, 10);
    // vec2 m3(nr, vec(10));
    //Eigen::MatrixXd m4(nr, 15);
    // vec2 m4(nr, vec(15));
    // SHOULD BE MERGABLE WITH ABOVE FOR
     //for (int h = 0; h < nr; h++) {
     //    m3[h][0] = std::pow(qx[h], 3);
     //    m3[h][1] = std::pow(qy[h], 3);
     //    m3[h][2] = std::pow(qz[h], 3);
     //    m3[h][3] = std::pow(qx[h], 2) * qy[h];
     //    m3[h][4] = std::pow(qx[h], 2) * qz[h];
     //    m3[h][5] = std::pow(qy[h], 2) * qx[h];
     //    m3[h][6] = std::pow(qy[h], 2) * qz[h];
     //    m3[h][7] = std::pow(qz[h], 2) * qx[h];
     //    m3[h][8] = std::pow(qz[h], 2) * qy[h];
     //    m3[h][9] = qx[h] * qy[h] * qz[h];
     //    m4[h][0] = std::pow(qx[h], 4);
     //    m4[h][1] = std::pow(qy[h], 4);
     //    m4[h][2] = std::pow(qz[h], 4);
     //    m4[h][3] = std::pow(qx[h], 3) * qy[h];
     //    m4[h][4] = std::pow(qx[h], 3) * qz[h];
     //    m4[h][5] = std::pow(qy[h], 3) * qx[h];
     //    m4[h][6] = std::pow(qy[h], 3) * qz[h];
     //    m4[h][7] = std::pow(qz[h], 3) * qx[h];
     //    m4[h][8] = std::pow(qz[h], 3) * qy[h];
     //    m4[h][9] = std::pow(qx[h], 2) * std::pow(qy[h], 2);
     //    m4[h][10] = std::pow(qx[h], 2) * std::pow(qz[h], 2);
     //    m4[h][11] = std::pow(qy[h], 2) * std::pow(qz[h], 2);
     //    m4[h][12] = std::pow(qx[h], 2) * qy[h] * qz[h];
     //    m4[h][13] = std::pow(qy[h], 2) * qx[h] * qz[h];
     //    m4[h][14] = std::pow(qz[h], 2) * qx[h] * qy[h];
     //}
    //m3.col(0) = (qx.pow(3)).matrix();
    //m3.col(1) = (qy.pow(3)).matrix();
    //m3.col(2) = (qz.pow(3)).matrix();
    //m3.col(3) = (qx.square() * qy).matrix();
    //m3.col(4) = (qx.square() * qz).matrix();
    //m3.col(5) = (qy.square() * qx).matrix();
    //m3.col(6) = (qy.square() * qz).matrix();
    //m3.col(7) = (qz.square() * qx).matrix();
    //m3.col(8) = (qz.square() * qy).matrix();
    //m3.col(9) = (qx * qy * qz).matrix();
    //m4.col(0) = qx.pow(4);
    //m4.col(1) = qy.pow(4);
    //m4.col(2) = qz.pow(4);
    //m4.col(3) = qx.pow(3) * qy;
    //m4.col(4) = qx.pow(3) * qz;
    //m4.col(5) = qy.pow(3) * qx;
    //m4.col(6) = qy.pow(3) * qz;
    //m4.col(7) = qz.pow(3) * qx;
    //m4.col(8) = qz.pow(3) * qy;
    //m4.col(9) = qx.square() * qy.square();
    //m4.col(10) = qx.square() * qz.square();
    //m4.col(11) = qy.square() * qz.square();
    //m4.col(12) = qx.square() * qy * qz;
    //m4.col(13) = qy.square() * qx * qz;
    //m4.col(14) = qz.square() * qx * qy;
	// Compute Debye-Waller factors
	//Eigen::MatrixXcd DW(nr, ncen);
    vec2 DW(ncen, vec(nr));
    //Eigen::Matrix3d Uij;
    for (int a = 0; a < ncen; a++) {
        vec2 ADPs = wave.get_atom(a).get_ADPs();
        vec2 Uij;
        if (level[a] > 0) {
            Uij = { { ADPs[0][0], ADPs[0][3], ADPs[0][4] },
                         { ADPs[0][3], ADPs[0][1], ADPs[0][5] },
						 { ADPs[0][4], ADPs[0][5], ADPs[0][2] } };
            //Uij << ADPs[0][0], ADPs[0][3], ADPs[0][4],
            //    ADPs[0][3], ADPs[0][1], ADPs[0][5],
            //    ADPs[0][4], ADPs[0][5], ADPs[0][2];
        }
        switch(level[a]) {
        case 0: {
                // Isotropic
                double U = U_iso[a];
                double temp;
                for (int r = 0; r < nr; r++) {
                    temp = -0.5 * (constants::TWO_PI * constants::TWO_PI) * U * (q[r][0] * q[r][0] + q[r][1] * q[r][1] + q[r][2] * q[r][2]);
                    DW[a][r] = std::exp(temp);
                }
                break;
            }
        case 1: {
                // Anisotropic U_ij
                double temp;
                for (int h = 0; h < nr; h++) {
                    vec q_ = { q[h][0], q[h][1], q[h][2] };
                    temp = -0.5 * (constants::TWO_PI * constants::TWO_PI) * dot_BLAS(dot(Uij, q_, true), q_, false);
                    DW[a][h] = std::exp(temp);
                }
                break;
            }
    //    case 2: {
    //            // Anisotropic C_ijk
    //            Eigen::VectorXd GCC = toEigenVector(ADPs[1]);
    //            Eigen::VectorXd term2 = -0.5 * (q * Uij).cwiseProduct(q).rowwise().sum();
    //            Eigen::VectorXcd term3 = (-std::complex<double>(0, 1) / 6.0) * (m3 * (mult3.array() * GCC.array()).matrix()).cast<std::complex<double>>();
				//DW.row(a) = (term2 + term3).array().exp();
    //            break;
    //        }
    //    case 3: {
    //            // Anisotropic D_ijkl
    //            Eigen::VectorXd GCC = toEigenVector(ADPs[1]);
    //            Eigen::VectorXd GCD = toEigenVector(ADPs[2]);
    //            Eigen::VectorXd term2 = -0.5 * (q * Uij).cwiseProduct(q).rowwise().sum();
    //            Eigen::VectorXcd term3 = (-std::complex<double>(0, 1) / 6.0) * (m3 * (mult3.array() * GCC.array()).matrix()).cast<std::complex<double>>();
    //            Eigen::VectorXd term4 = (1.0 / 24.0) * (m4 * (mult4.array() * GCD.array()).matrix());
				//DW.row(a) = (term2 + term3 + term4).array().exp();
    //            break;
    //        }
        }
    }
    set_DW(DW);
    // closing void
}

void XCW::eval_phase() {
    const double bohr2angstrom = constants::bohr2ang(1);
    const double angstrom2bohr = constants::ang2bohr(1);
    const double scale = angstrom2bohr / constants::TWO_PI;
    cdouble exponent;
    cvec2 phase(ncen, cvec(nr));

    vec2 cm = { { unit_cell.get_cm(0,0), unit_cell.get_cm(0,1), unit_cell.get_cm(0,2)},
                            { unit_cell.get_cm(1,0), unit_cell.get_cm(1,1), unit_cell.get_cm(1,2)},
                            { unit_cell.get_cm(2,0), unit_cell.get_cm(2,1), unit_cell.get_cm(2,2)} };
    std::transform(cm.begin(), cm.end(), cm.begin(), [bohr2angstrom](std::vector<double>& vec) {
        std::transform(vec.begin(), vec.end(), vec.begin(), [bohr2angstrom](double x) { return x * bohr2angstrom; });
        return vec;});
    for (int at = 0; at < ncen; at++) {
        vec pos_frac = { asym_atoms[at].frac_pos[0], asym_atoms[at].frac_pos[1], asym_atoms[at].frac_pos[2] };
        vec new_pos_cart = dot(cm, pos_frac, true);
        for (int r = 0; r < nr; r++) {
            vec q = { k_pt[0][r], k_pt[1][r], k_pt[2][r] };
            std::transform(q.begin(), q.end(), q.begin(), [scale](double x) { return x * scale; });
            exponent = cdouble(0, constants::TWO_PI * dot_BLAS(q, new_pos_cart, false));
            phase[at][r] = std::exp(exponent);
        }
    }
    set_phase(phase);
}

void XCW::eval_translation_phase() {
    const double angstrom2bohr = constants::ang2bohr(1);
    const double bohr2angstrom = constants::bohr2ang(1);
    const double scale = angstrom2bohr / constants::TWO_PI;
    vec2 trans = unit_cell.get_trans();
    vec2 cm = { { unit_cell.get_cm(0,0), unit_cell.get_cm(0,1), unit_cell.get_cm(0,2)},
                                  { unit_cell.get_cm(1,0), unit_cell.get_cm(1,1), unit_cell.get_cm(1,2)},
                                  { unit_cell.get_cm(2,0), unit_cell.get_cm(2,1), unit_cell.get_cm(2,2)} };
    std::transform(cm.begin(), cm.end(), cm.begin(), [bohr2angstrom](std::vector<double>& vec) {
              std::transform(vec.begin(), vec.end(), vec.begin(), [bohr2angstrom](double x) { return x * bohr2angstrom; });
              return vec;});
    cvec2 phase(nr_small, cvec(trans[0].size(), 0));
    for (int r = 0; r < nr_small; r++) {
        ivec asym_list = generate_asym_lookup(r);
        vec q_temp = { k_pt[0][asym_list[0]], k_pt[1][asym_list[0]], k_pt[2][asym_list[0]] };
        std::transform(q_temp.begin(), q_temp.end(), q_temp.begin(), [scale](double x) { return x * scale; });
        for (int t = 0; t < trans[0].size(); t++) {
            vec trans_temp = { trans[0][t], trans[1][t], trans[2][t] };
            trans_temp = dot(cm, trans_temp, true);
            cdouble exponent(0, constants::TWO_PI * dot_BLAS(q_temp, trans_temp, false));
            phase[r][t] = std::exp(exponent);
        }
    }
    set_translation_phase(phase);
    // closing function
}

cvec3 XCW::eval_I() {
    // Debugging, calculating atomic scattering factors for single reflex
    //cvec3 atomic_scattering_factor(nmo, cvec2(nmo, cvec(ncen)));

    cvec3 I(nmo, cvec2(nmo, cvec(nr_small, 0)));
    int at = 0, mu = 0, nu = 0, r = 0, s = 0, r_asym = 0;

    // Do the setup of primitives here so it's only done once
    std::vector<ao_data> ao_data_shells;
    const std::map<int, LibCintBasis>& basis_sets_map = params.get_basis_sets();

    // Create primitives for each shell of each atom and store in ao_data_shells
    for (int atm = 0; atm < ncen; atm++) {
		const atom atom = wave.get_atom(asym_atom_list[atm]);
        d3 pos = atom.get_pos();
        int atomic_number = atom.get_charge();
        if (basis_sets_map.find(atomic_number) == basis_sets_map.end()) {
            continue;
        }
        const LibCintBasis& basis = basis_sets_map.at(atomic_number);
        int prim_idx = 0;
        for (int shell = 0; shell < basis.shellcount.size(); shell++) {
            const int shell_type = basis.shelltypes[shell];
            int shell_start = prim_idx;
            std::vector<primitive> tmp_prims;
            for (int prim = 0; prim < basis.shellcount[shell]; prim++, prim_idx++) {
                if (prim_idx >= basis.exponents.size()) break;
                const double alpha = basis.exponents[prim_idx];
                const double coeff = basis.coefficients[prim_idx];
                tmp_prims.emplace_back(0, shell_type, alpha, coeff);
            }
			const int l = shell_type;
            for (int m = -l; m <= l; m++) {
				ao_data_shells.push_back({ tmp_prims, pos, m });
            }
        }
    }
	cvec2 XCW_integral;

 //   vec ovlp;
 //   Int_Params wavy_params(wave);
 //   compute2C<Overlap2C_SPH>(wavy_params, ovlp);
	//dMatrixRef2 ovlp_mat(ovlp.data(), nmo, nmo);

    GridConfiguration config;
    config.accuracy = opt->accuracy;
    config.partition_type = opt->partition_type;
    config.no_density_eval = true;
    config.pbc = opt->pbc;
    config.debug = opt->debug;
    config.all_charges = opt->all_charges;
    GridManager grid_manager(config);
    WFN temp = wave;
    temp.delete_unoccupied_MOs();
    grid_manager.setup3DGridsForMolecule(temp, asym_atom_list, needs_grid, unit_cell);
    int cofs;
    // Computes I as a three dimensional tensor mu x nu x hkl
    bool equal = false;
    double cutoff = 0;
    std::vector<ivec> asym_lookup(nr_small);
    for (r = 0; r < nr_small; r++) {
        asym_lookup[r] = generate_asym_lookup(r);
    }
	for (mu = 0; mu < nmo; mu++) {
        vec2 mu_vals;
        const ao_data& mu_prims = ao_data_shells[mu];
		const std::vector<primitive>& mu_primitives = mu_prims.prims;
        const double& mp0 = mu_prims.pos[0];
        const double& mp1 = mu_prims.pos[1];
        const double& mp2 = mu_prims.pos[2];
        for (nu = mu; nu < nmo; nu++) {
            if (nu == mu) {
                equal = true;
            }
            else {
                equal = false;
            }
            const ao_data& nu_prims = ao_data_shells[nu];
            const std::vector<primitive>& nu_primitives = nu_prims.prims;
            const double& np0 = nu_prims.pos[0];
            const double& np1 = nu_prims.pos[1];
            const double& np2 = nu_prims.pos[2];
            // Screening
            const double dist0 = mp0 - np0;
            const double dist1 = mp1 - np1;
            const double dist2 = mp2 - np2;
            const double dist = std::sqrt(dist0 * dist0 + dist1 * dist1 + dist2 * dist2);
            if (!equal) {
                double c_tot = 0;
                for (cofs = 0; cofs < mu_primitives.size(); cofs++) {
					double alpha = mu_primitives[cofs].get_exp();
					double beta = nu_primitives[cofs].get_exp();
					double alpha_k_l = alpha * beta / (alpha + beta);
                    c_tot += std::abs(mu_primitives[cofs].get_coef()) * std::abs(nu_primitives[cofs].get_coef()) * std::pow(constants::PI / alpha_k_l, 1.5);
                }
                c_tot *= std::sqrt((2*mu_primitives[0].get_type() + 1) * (2 * nu_primitives[0].get_type() + 1) / (constants::TWO_PI * constants::TWO_PI));
                double alpha_min = mu_primitives[mu_primitives.size() - 1].get_exp();
                double beta_min = nu_primitives[nu_primitives.size() - 1].get_exp();
                double a_min = (alpha_min + beta_min) / (alpha_min * beta_min);
                const double screen = 0.001;
                cutoff = std::sqrt(a_min * std::log(c_tot / screen));
            }
            if (dist > cutoff && !equal) {
                continue;
            }

            cdouble temp1;
            double asym_fact;
            int idx;
            XCW_integral = calculateXCWintegral(grid_manager, mu_prims, nu_prims, asym_atom_list, k_pt, equal, mu_vals);
            //for (at = 0; at < ncen; at++) {
            //    atomic_scattering_factor[mu][nu][at] = XCW_integral[at][0];
            //}
            for (r = 0; r < nr_small; r++) {
				const cvec& translation_phase_r = translation_phase[r];
                const ivec& asym_list = asym_lookup[r];
                for (at = 0; at < ncen; at++) {
					asym_fact = asym_atoms[at].asym_fact;
					const cvec& XCW_at = XCW_integral[at];
					const vec& DW_at = DW_fact[at];
					const cvec& phase_at = phase_fact[at];
                    temp1 = 0;
                    for (r_asym = 0; r_asym < asym_list.size(); r_asym++) {
						idx = asym_list[r_asym];
                        temp1 += XCW_at[idx] * DW_at[idx] * phase_at[idx] * translation_phase_r[r_asym];
                    }
                    I[mu][nu][r] += temp1 * asym_fact;
                }
            }
        }
        // Print out that contributions from mu are done
		std::cout << "Finished contributions from mu = " << mu << std::endl;
    }
    dMatrix2 D = wave.get_dm();

     //Debugging, calculating atomic scattering factors
    //cvec scattering_factors(ncen);
    //for (int at = 0; at < ncen; at++) {
    //    for (mu = 0; mu < nmo; mu++) {
    //        for (nu = mu; nu < nmo; nu++) {
    //            double D_ = D(mu, nu);
    //            if (mu == nu) {
    //                scattering_factors[at] += atomic_scattering_factor[mu][nu][at] * D_;
    //            }
    //            else {
    //                scattering_factors[at] += 2.00 * atomic_scattering_factor[mu][nu][at] * D_;
    //            }
    //        }
    //    }
    //}
    //std::cout << "Atomic Scattering Factors" << std::endl;
    //for (at = 0; at < ncen; at++) {
    //    std::cout << scattering_factors[at] << std::endl;
    //}
    return I;
}

ivec XCW::generate_asym_lookup(const int r) {
        ivec asym_list;
        auto it = hkl.begin();
        std::advance(it, r);
        ivec3 rots = unit_cell.get_sym();
        i3 tempv;
        const i3& hkl_temp = *it;
        for (int s = 0; s < rots[0][0].size(); s++) {
            tempv = { 0, 0, 0 };
            for (int h = 0; h < 3; h++) {
                for (int j = 0; j < 3; j++) {
                    tempv[j] += hkl_temp[h] * rots[j][h][s];
                }
            }
            int idx_ = 0;
            auto idx = hkl_enlarged.find(tempv);
            if (idx != hkl_enlarged.end()) {
                idx_ = std::distance(hkl_enlarged.begin(), idx);
            }
            asym_list.push_back(idx_);
        }
    return asym_list;
    // closing function
}

cvec2 XCW::calculateXCWintegral(GridManager& grid_manager, const ao_data& mu_data, const ao_data& nu_data, const ivec& asym_atom_list, vec2& k_pt, bool& equal, vec2& mu_vals) {
    GridData& GD = grid_manager.getGridData();
    vec2* grids = grid_manager.getNeedsHelper() ? GD.helper_grids.data() : GD.atomic_grids.data();
    const int n_grids = grid_manager.getNeedsHelper() ? GD.helper_grids.size() : GD.atomic_grids.size();
    const int* points = grid_manager.getNeedsHelper() ? GD.helper_num_points_per_atom.data() : GD.num_points_per_atom.data();
    const double mp0 = mu_data.pos[0];
    const double mp1 = mu_data.pos[1];
    const double mp2 = mu_data.pos[2];
    const double np0 = nu_data.pos[0];
    const double np1 = nu_data.pos[1];
    const double np2 = nu_data.pos[2];
    const std::vector<primitive> mu_prims = mu_data.prims;
    const std::vector<primitive> nu_prims = nu_data.prims;
    if (equal) {
		mu_vals.resize(n_grids);
		for (int g = 0; g < n_grids; g++) {
			const int num_points = points[g];
			mu_vals[g].resize(num_points);
			vec2& atom_grid = grids[g];
			const double* x_ptr = atom_grid[GridData::GridIndex::X].data();
			const double* y_ptr = atom_grid[GridData::GridIndex::Y].data();
			const double* z_ptr = atom_grid[GridData::GridIndex::Z].data();
			double* overlap_ptr = atom_grid[GridData::GridIndex::WFN_DENSITY].data();
#pragma omp parallel for schedule(dynamic, 4)
			for (int p = 0; p < num_points; p++) {
				std::array<double, 4> d_mu;
				double* d_mu_ptr = d_mu.data();
				*d_mu_ptr = x_ptr[p] - mp0;
				*(d_mu_ptr + 1) = y_ptr[p] - mp1;
				*(d_mu_ptr + 2) = z_ptr[p] - mp2;
				*(d_mu_ptr + 3) = std::hypot(*(d_mu_ptr), *(d_mu_ptr + 1), *(d_mu_ptr + 2));
				mu_vals[g][p] = wave.eval_ao(d_mu, mu_prims, mu_data.m);
				overlap_ptr[p] = std::pow(mu_vals[g][p], 2);
			}
		}
    }
    else {
		for (int g = 0; g < n_grids; g++) {
			const int num_points = points[g];
			vec2& atom_grid = grids[g];
			const double* x_ptr = atom_grid[GridData::GridIndex::X].data();
			const double* y_ptr = atom_grid[GridData::GridIndex::Y].data();
			const double* z_ptr = atom_grid[GridData::GridIndex::Z].data();
			double* overlap_ptr = atom_grid[GridData::GridIndex::WFN_DENSITY].data();
			#pragma omp parallel for schedule(dynamic, 4)
			for (int p = 0; p < num_points; p++) {
				std::array<double, 4> d_mu, d_nu;
				double* d_mu_ptr = d_mu.data();
				double* d_nu_ptr = d_nu.data();
				*d_mu_ptr = x_ptr[p] - mp0;
				*(d_mu_ptr + 1) = y_ptr[p] - mp1;
				*(d_mu_ptr + 2) = z_ptr[p] - mp2;
				*(d_mu_ptr + 3) = std::hypot(*(d_mu_ptr), *(d_mu_ptr + 1), *(d_mu_ptr + 2));
				*d_nu_ptr = x_ptr[p] - np0;
				*(d_nu_ptr + 1) = y_ptr[p] - np1;
				*(d_nu_ptr + 2) = z_ptr[p] - np2;
				*(d_nu_ptr + 3) = std::hypot(*(d_nu_ptr), *(d_nu_ptr + 1), *(d_nu_ptr + 2));
                overlap_ptr[p] = (*d_mu_ptr + 3) > 12 || (*d_nu_ptr + 3) > 12 ? 0 : mu_vals[g][p] * wave.eval_ao(d_nu, nu_prims, nu_data.m);
			}
		}
    }
    vec2 d1, d2, d3, dens;
    std::vector<_time_point> time_points({ get_time() });
    _time_point end;
    grid_manager.getDensityVectors(wave, asym_atom_list, d1, d2, d3, dens);
    const int points_ = grid_manager.getTotalGridPoints();
    cvec2 sf;
    calc_SF(points_, k_pt, d1, d2, d3, dens, sf, std::cout, time_points.front(), end, opt->debug, true, true);
    return sf;
}

void XCW::parse_anom_atoms(std::vector<anom_atom>& anom_atoms) {
    std::ifstream file(opt->anom_disp_path);
    if (!file) {
		std::cout << "Could not open anomalous dispersion file. Continuing without anomalous dispersions." << std::endl;
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty())
            continue;
        std::istringstream iss(line);
        std::string symbol;
        double real_part, imag_part;
        if (iss >> symbol >> real_part >> imag_part) {
            if (!symbol.empty() && symbol[0] != '_' && symbol != "loop_") {
                anom_atoms.push_back({symbol, cdouble(real_part, imag_part)});
            }
        }
    }
}

cvec XCW::eval_anom_disp() {
    std::vector<anom_atom> anom_atoms;
    parse_anom_atoms(anom_atoms);
    cvec corr(nr_small, 0);
    int r, at, r_asym;
	for (int at = 0; at < ncen; at++) {
        const char* symbol = constants::atnr2letter(asym_atoms[at].type[0]);
		for (const anom_atom& anom_atom : anom_atoms) {
			if (symbol == anom_atom.identifier) {
				asym_atoms[at].anom = anom_atom.dispersion;
				break;
			}
		}
	}
    std::vector<ivec> asym_lookup(nr_small);
    for (r = 0; r < nr_small; r++) {
        asym_lookup[r] = generate_asym_lookup(r);
    }
    for (r = 0; r < nr_small; r++) {
        const ivec& lookup = asym_lookup[r];
        for (at = 0; at < ncen; at++) {
            cdouble temp1 = 0;
            for (r_asym = 0; r_asym < lookup.size(); r_asym++) {
                temp1 += phase_fact[at][lookup[r_asym]] * DW_fact[at][lookup[r_asym]] * translation_phase[r][r_asym];
            }
			corr[r] += temp1 * asym_atoms[at].asym_fact * asym_atoms[at].anom;
        }
    }
    return corr;
}

void XCW::calc_F_calc(const cvec3& I, const cvec& corr) {
    cvec F_calc(nr_small, 0);
    dMatrix2 D = wave.get_dm();
    int mu, nu, r;
    for (r = 0; r < nr_small; r++) {
        for (mu = 0; mu < nmo; mu++) {
            for (nu = mu; nu < nmo; nu++) {
                double D_ = D(mu, nu);
                if (mu == nu) {
                    F_calc[r] += I[mu][nu][r] * D_;
                }
                else {
                    F_calc[r] += 2.00 * I[mu][nu][r] * D_;
                }
            }
        }
		F_calc[r] += corr[r];
    }
    std::cout << "F_calc values" << std::endl;
    for (cdouble val : F_calc) {
        std::cout << std::fixed << std::setprecision(5) << std::pow(std::abs(val), 2) << std::endl;
    }
    std::cout << "Finished" << std::endl;
    // closing function
}

//bool set_occ_data_path(const std::filesystem::path& path)
//{
//#ifdef _WIN32
//    return _putenv_s("OCC_DATA_PATH", path.string().c_str()) == 0;
//#else
//    return setenv("OCC_DATA_PATH", path.string().c_str(), 1) == 0;
//#endif
//}

//void XCW::do_SCF() {
//    // Use occ to load a molecule from the water.xyz file and create a 6-31G basis set for it, then run Hartree-Fock on it and print the SCF energy
//
//    bool x = set_occ_data_path(std::filesystem::current_path());
//	occ::core::Molecule mol = occ::io::molecule_from_xyz_file("water.xyz");
//	occ::qm::AOBasis bs = occ::qm::AOBasis::load(mol.atoms(), "def2-TZVP");
//    occ::qm::HartreeFock hf(bs);
//    occ::qm::SCF scf(hf, occ::qm::SpinorbitalKind::Restricted);
//	scf.set_charge_multiplicity(0,1);
//
//    const double& energy_conv = scf.convergence_settings.energy_threshold;
//    const double& commutator_conv = scf.convergence_settings.commutator_threshold;
//    bool converged = false;
//	int iteration = 0;
//
//	// Initial SCF setup
//    bool incremental = false;
//    scf.update_occupied_orbital_count();
//    scf.compute_initial_guess();
//    scf.ctx.K = scf.m_procedure.compute_schwarz_ints();
//    occ::Mat D_diff = scf.ctx.mo.D;
//    occ::Mat D_last;
//    occ::Mat FD_comm = occ::Mat::Zero(scf.ctx.F.rows(), scf.ctx.F.cols());
//    occ::Mat F_diis;
//    double ehf_last;
//    scf.update_scf_energy(incremental);
//    //occ::log::info("starting {} scf iterations", scf.scf_kind());
//    //occ::log::debug("{} electrons total", scf.ctx.n_electrons);
//    //occ::log::debug("{} alpha electrons", scf.n_alpha());
//    //occ::log::debug("{} beta electrons", scf.n_beta());
//    //occ::log::debug("net charge {}", scf.charge());
//    scf.total_time = 0.0;
//
//    // SCF Loop
//    do {
//        const auto tstart = std::chrono::high_resolution_clock::now();
//		std::cout << "SCF Iteration " << scf.iter + 1 << std::endl;
//        scf.iter++;
//
//        SCF_iteration(scf, ehf_last, D_last, D_diff, incremental, F_diis);
//
//        const auto tstop = std::chrono::high_resolution_clock::now();
//        const std::chrono::duration<double> time_elapsed = tstop - tstart;
//
//        if (scf.iter == 1) {
//            occ::log::info("{:>4s} {: >20s} {: >12s} {: >12s}  {: >8s}", "#",
//                "E (Hartrees)", "|dE|/E", "max|FDS-SDF|", "T (s)");
//        }
//        occ::log::info("{:>4d} {:>20.12f} {:>12.5e} {:>12.5e}  {:>8.2e}", scf.iter,
//            scf.ctx.energy["total"], scf.ediff_rel, scf.diis_error, time_elapsed.count());
//        occ::log::flush();
//
//        SCF_convergence_check(scf, converged, energy_conv, commutator_conv, ehf_last);
//
//        scf.total_time += time_elapsed.count();
//    } while (!converged && (scf.iter < scf.maxiter));
//    occ::log::info("{} spinorbital SCF energy converged after {:.5f} seconds",
//        scf.scf_kind(), scf.total_time);
//    occ::log::info("{}", scf.ctx.energy.to_string());
//    scf.ctx.converged = true;
//    double e = scf.ctx.energy["total"];
//}

//void XCW::SCF_iteration(auto& scf, double& ehf_last, occ::Mat& D_last, occ::Mat& D_diff, bool& incremental, occ::Mat& F_diis) {
//    ehf_last = scf.ctx.energy["electronic"];
//	std::cout << std::setprecision(12) << "Ehf last: " << ehf_last << std::endl;
//    D_last = scf.ctx.mo.D;
//    scf.ctx.H = scf.ctx.T + scf.ctx.V; // +scf.ctx.Vecp + scf.ctx.V_ext; 
//    scf.m_procedure.update_core_hamiltonian(scf.ctx.mo, scf.ctx.H); 
//	incremental = true; // Memory for whether or not the current Fock matrux was built incrementally
//
//    // Managing incremental Fock building
//    //If the incremental Fock build has not started yet but conditions allow it
//    if (!scf.incremental_Fbuild_started &&
//        scf.convergence_settings.start_incremental_fock(scf.diis_error)) {
//        scf.incremental_Fbuild_started = true;
//        scf.reset_incremental_fock_formation = false;
//        scf.last_reset_iteration = scf.iter - 1;
//        scf.next_reset_threshold = scf.diis_error / 10;
//        occ::log::debug("starting incremental fock build");
//    }
//
//    //If the incremental Fock build needs to be reset or has not started yet
//    if (scf.reset_incremental_fock_formation ||
//        !scf.incremental_Fbuild_started) {
//        scf.ctx.F = scf.ctx.H;
//        D_diff = scf.ctx.mo.D;
//        incremental = false;
//    }
//
//    //If the incremental Fock build already started but a reset is due
//    if (scf.reset_incremental_fock_formation && scf.incremental_Fbuild_started) {
//        scf.reset_incremental_fock_formation = false;
//        scf.last_reset_iteration = scf.iter;
//        scf.next_reset_threshold = scf.diis_error / 10;
//        occ::log::debug("resetting incremental fock build");
//    }
//
//    // build a new Fock matrix
//    std::swap(scf.ctx.mo.D, D_diff);
//    scf.ctx.F += scf.m_procedure.compute_fock(scf.ctx.mo, scf.ctx.K);
//    std::swap(scf.ctx.mo.D, D_diff);
//
//    // compute HF energy with the non-extrapolated Fock matrix
//    scf.update_scf_energy(incremental);
//	std::cout << std::setprecision(12) << "Electronic energy after SCF iteration: " << scf.ctx.energy["electronic"] << std::endl;
//    scf.ediff_rel = std::abs((scf.ctx.energy["electronic"] - ehf_last) / scf.ctx.energy["electronic"]);
//
//    F_diis = scf.diis.update(scf.ctx.S, scf.ctx.mo.D, scf.ctx.F);
//    // double prev_error = diis_error;
//    scf.diis_error = scf.diis.max_error();
//    /*
//    bool use_ediis = (diis_error > 1e-1) || (prev_error /
//    diis.min_error() > 1.1);
//
//    Mat F_ediis = scf.ediis.update(scf.ctx.D, scf.ctx.F, scf.ctx.energy["electronic"]);
//    if(use_ediis) F_diis = F_ediis;
//    else if(diis_error > 1e-4) {
//        F_diis = (10 * diis_error) * F_ediis + (1 - 10 * diis_error)
//    * F_diis;
//    }
//    */
//
//    // Check for incremental Fock building
//    if (scf.diis_error < scf.next_reset_threshold || scf.iter - scf.last_reset_iteration >= 15)
//        scf.reset_incremental_fock_formation = true;
//
//    scf.ctx.orthogonalizer.orthogonalize_molecular_orbitals(scf.ctx.mo, F_diis);
//    D_diff = scf.ctx.mo.D - D_last;
//    //closing function
//}

//void XCW::SCF_convergence_check(auto& scf, bool& converged, const double& energy_conv, const double& commutator_conv, const double& ehf_last) {
//    if (scf.convergence_settings.energy_and_commutator_converged(scf.ediff_rel, scf.diis_error)) {
//        converged = true;
//    }
//    // closing function
//}


void XCW::calc_F_calc_fast(const cvec& corr) {
	cvec2 atomic_scattering_factors(ncen, cvec(nr, 0));
	cvec F_calc(nr_small, 0);

    // Calculate atomic scattering factors for each atom and symmetry generated reflexes
    GridConfiguration config;
    config.accuracy = opt->accuracy;
    config.partition_type = opt->partition_type;
    config.pbc = opt->pbc;
    config.no_density_eval = false;
    config.debug = opt->debug;
    config.all_charges = opt->all_charges;
    GridManager grid_manager(config);
    WFN temp = wave;
    vec2 d1, d2, d3, dens;
    std::vector<_time_point> time_points({ get_time() });
    _time_point end;
    temp.delete_unoccupied_MOs();
    grid_manager.setup3DGridsForMolecule(temp, asym_atom_list, needs_grid, unit_cell);
    grid_manager.calculateNonSphericalDensities(temp, unit_cell);
    svec time_descriptions;
    grid_manager.addTimingInfoToVecs(time_points, time_descriptions);
    PartitionResults results = grid_manager.calculatePartitionedCharges(temp, unit_cell);
    grid_manager.printChargeTable(labels, temp, asym_atom_list, std::cout, results);
    time_points.push_back(get_time());
    time_descriptions.push_back("calculate charges");
    grid_manager.getDensityVectors(temp, asym_atom_list, d1, d2, d3, dens);
    const int points = grid_manager.getTotalGridPoints();
    calc_SF(points, k_pt, d1, d2, d3, dens, atomic_scattering_factors, std::cout, time_points.front(), end, opt->debug, true, true);
    // Calculate F_calc
	auto it = hkl.begin();
	for (int r = 0; r < nr_small; r++) {
		const ivec& lookup = generate_asym_lookup(r);
		for (int at = 0; at < ncen; at++) {
			for (int r_asym = 0; r_asym < lookup.size(); r_asym++) {
				F_calc[r] += atomic_scattering_factors[at][lookup[r_asym]] * DW_fact[at][lookup[r_asym]] * phase_fact[at][lookup[r_asym]] * translation_phase[r][r_asym] * asym_atoms[at].asym_fact;
			}
		}
		// Add anomalous dispersion correction
        const i3& hkl_temp = *it;
        it++;
        F_calc[r] += corr[r];
		std::cout << "F_calc for hkl = " << hkl_temp[0] << ", " << hkl_temp[1] << ", " << hkl_temp[2] << ": " << std::fixed << std::setprecision(5) << std::pow(std::abs(F_calc[r]), 2) << std::endl;
	}
	//dump F_calc values as binary file called F_calc
	std::ofstream fout("F_calc.bin", std::ios::out | std::ios::binary);
	fout.write(reinterpret_cast<const char*>(F_calc.data()), F_calc.size() * sizeof(cdouble));
    fout.close();

}