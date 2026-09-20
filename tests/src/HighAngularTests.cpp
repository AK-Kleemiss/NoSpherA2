
#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/fchk.h"
#include "core/AtomGrid.h"
#include "core/SALTED_utilities.h"
#include "core/scattering_factors.h"
#include "core/nos_math.h"
#include "core/GridManager.h"
#include "core/atoms.h"
#include "core/tsc_block.h"
#include "core/cell.h"
#include "core/wfn_class.h"
#include "core/properties.h"
#include "core/integrator.h"
#include "core/integration_params.h"
#include "core/basis_set.h"
#include "core/geometry_aid.h"
#include "core/crystal_energies.h"
#include "core/NoSpherA2.h"
#include "core/isosurface.h"
#include "core/npy.h"
#include "core/libCintMain.h"
#include "core/throughput.h"
#include "core/i_tensor_stream.h"
#include "core/tsc_label_converter.h"
#include "core/sphere_lebedev_rule.h"
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#undef I
#ifdef NOSPHERA2_USE_GPU
#include "core/blas_gpu.h"
#include "core/aux_density_gpu.h"
#endif

static constexpr double PI_VAL = 3.14159265358979323846;


namespace NoSpherA2UnitTests
{
    TEST(GbwHighAngularTests, ReadsIFunctionFixtures)
    {
        const auto root = nos_test_repo_root();
        const std::array<std::string, 3> angles = { "71", "113", "149" };
        for (const auto& angle : angles) {
            const auto input = root / "tests" / "CuF2_i_func" / angle / "calc.gbw";
            if (!std::filesystem::exists(input)) GTEST_SKIP() << "Missing " << input;
            WFN wave(input, false);
            EXPECT_EQ(wave.get_nmo(), 670);
            EXPECT_GT(wave.get_nex(), 0);
            EXPECT_NE(std::find(wave.get_types().begin(), wave.get_types().end(), 56), wave.get_types().end());
            EXPECT_EQ(*std::max_element(wave.get_types().begin(), wave.get_types().end()), 84);
            const double density = wave.compute_dens({ 0.23, -0.41, 0.67 });
            EXPECT_TRUE(std::isfinite(density));
            EXPECT_GT(density, 0.0);
        }
    }

    //<x^a e^-ar^2 | x^b e^-br^2> along one axis, PA/PB the centre P = (aA+bB)/p relative to A and B
    static double prim_overlap_1d(const int a, const int b, const double PA, const double PB, const double p)
    {
        double s = 0;
        for (int i = 0; i <= a; i++)
            for (int j = i % 2; j <= b; j += 2)
            {
                double t = std::pow(PA, a - i) * std::pow(PB, b - j);
                for (int k = 1; k <= i; k++) t *= double(a - i + k) / k;
                for (int k = 1; k <= j; k++) t *= double(b - j + k) / k;
                for (int k = i + j - 1; k > 0; k -= 2) t *= k;
                s += t / std::pow(2 * p, (i + j) / 2);
            }
        return s * std::sqrt(constants::PI / p);
    }

    static vec2 primitive_overlap(const WFN& wave)
    {
        const int nex = wave.get_nex();
        vec2 S(nex, vec(nex));
#pragma omp parallel for schedule(dynamic)
        for (int a = 0; a < nex; a++)
        {
            int la[3], lb[3];
            constants::type2vector(wave.get_type(a), la);
            const double al = wave.get_exponent(a);
            for (int b = 0; b <= a; b++)
            {
                constants::type2vector(wave.get_type(b), lb);
                const double be = wave.get_exponent(b), p = al + be;
                double s = 1, AB2 = 0;
                for (int k = 0; k < 3; k++)
                {
                    const double A = wave.get_atom_coordinate(wave.get_center(a) - 1, k), B = wave.get_atom_coordinate(wave.get_center(b) - 1, k), P = (al * A + be * B) / p;
                    AB2 += (A - B) * (A - B);
                    s *= prim_overlap_1d(la[k], lb[k], P - A, P - B, p);
                }
                S[a][b] = S[b][a] = s * std::exp(-al * be / p * AB2);
            }
        }
        return S;
    }

    //ORCA stores every contraction normalised, so each shell's m = 0 function and every occupied MO of the i-function
    //GBWs has unit norm under the analytic primitive overlap; this pins the h and i tables where the grid integral cannot
    TEST(GbwHighAngularTests, OccupiedMOsAreNormalised_full)
    {
        if (const char* env = std::getenv("RUN_FULL_TEST"); !env || std::string(env) == "0" || std::string(env) == "false")
            GTEST_SKIP() << "Set RUN_FULL_TEST=1 to build the primitive overlap of the CuF2 i-function GBWs";
        const auto root = nos_test_repo_root();
        const std::array<std::string, 3> angles = { "71", "113", "149" };
        for (const auto& angle : angles) {
            const auto input = root / "tests" / "CuF2_i_func" / angle / "calc.gbw";
            if (!std::filesystem::exists(input)) GTEST_SKIP() << "Missing " << input;
            WFN wave(input, false);
            wave.delete_unoccupied_MOs();
            const int nex = wave.get_nex(), nmo = wave.get_nmo();
            const vec2 S = primitive_overlap(wave);
            for (int at = 0; at < wave.get_ncen(); at++)
            {
                const atom& A = wave.get_atom(at);
                int prim = 0;
                for (int sh = 0; sh < (int)A.get_shellcount_size(); sh++)
                {
                    const int n = A.get_shellcount(sh), l = A.get_basis_set_type(prim) - 1, nc = constants::n_cart(l), nsph = constants::n_spher(l);
                    double norm = 0;
                    for (int u = 0; u < n; u++)
                        for (int v = 0; v < n; v++)
                        {
                            const double p = A.get_basis_set_exponent(prim + u) + A.get_basis_set_exponent(prim + v);
                            double t = 0;
                            for (int c1 = 0; c1 < nc; c1++)
                                for (int c2 = 0; c2 < nc; c2++)
                                {
                                    int e1[3], e2[3]; constants::type2vector(constants::first_type[l] + c1, e1); constants::type2vector(constants::first_type[l] + c2, e2);
                                    double o = 1;
                                    for (int k = 0; k < 3; k++) o *= prim_overlap_1d(e1[k], e2[k], 0, 0, p);
                                    t += constants::sph2cart(l)[c1 * nsph] * constants::sph2cart(l)[c2 * nsph] * o;
                                }
                            norm += A.get_basis_set_coefficient(prim + u) * A.get_basis_set_coefficient(prim + v) * t;
                        }
                    EXPECT_NEAR(norm, 1.0, 1e-12) << angle << " atom " << at << " shell " << sh << " l " << l;
                    prim += n;
                }
            }
            double electrons = 0;
            for (int i = 0; i < nmo; i++)
            {
                double norm = 0;
                for (int a = 0; a < nex; a++)
                {
                    double t = 0;
                    for (int b = 0; b < nex; b++) t += S[a][b] * wave.get_MO_coef(i, b);
                    norm += wave.get_MO_coef(i, a) * t;
                }
                EXPECT_NEAR(norm, 1.0, 1e-8) << angle << " MO " << i;
                electrons += wave.get_MO_occ(i) * norm;
            }
            EXPECT_NEAR(electrons, wave.get_nr_electrons(), 1e-6) << angle;
        }
    }

    TEST(GbwHighAngularTests, IntegratesElectronCount_full)
    {
        if (const char* env = std::getenv("RUN_FULL_TEST"); !env || std::string(env) == "0" || std::string(env) == "false")
            GTEST_SKIP() << "Set RUN_FULL_TEST=1 to integrate the CuF2 i-function densities on a Becke grid";
        const auto root = nos_test_repo_root();
        const std::array<std::string, 3> angles = { "71", "113", "149" };
        for (const auto& angle : angles) {
            const auto input = root / "tests" / "CuF2_i_func" / angle / "calc.gbw";
            if (!std::filesystem::exists(input)) GTEST_SKIP() << "Missing " << input;
            WFN wave(input, false);
            wave.delete_unoccupied_MOs();
            GridConfiguration config;
            config.partition_type = PartitionType::Becke;
            config.accuracy = 4;
            GridManager gm(config);
            ivec atom_list(wave.get_ncen());
            for (int i = 0; i < wave.get_ncen(); i++) atom_list[i] = i;
            gm.setup3DGridsForMolecule(wave, atom_list);
            const PartitionResults res = gm.calculatePartitionedCharges(wave);
            EXPECT_NEAR(res.overall_charges[PartitionResults::S_BECKE], wave.get_nr_electrons(), 1e-3) << angle;
        }
    }
    //H2 in a made-up spherical basis with one primitive per l, converted through the OCC constructor: every MO keeps
    //unit norm and stays orthogonal to the others under the analytic primitive overlap, which pins the sph2cart tables,
    //the |m| phase flips and the type order for every l at once. OCC's Hartree-Fock provides the MOs up to h; beyond
    //that libcint's Rys quadrature breaks down on MSVC (12 roots and up), so the l = 10 MOs are OCC's overlap matrix
    //Loewdin-orthonormalised with the lowest one doubly occupied
    static WFN occ_h2(const int lmax)
    {
        spdlog::set_level(spdlog::level::err);
        const std::vector<occ::core::Atom> atoms{ { 1, 0.0, 0.0, -0.7 }, { 1, 0.0, 0.0, 0.7 } };
        std::vector<occ::gto::Shell> shells;
        for (const auto& at : atoms)
            for (int l = 0; l <= lmax; l++)
            {
                shells.emplace_back(l, std::vector<double>{ l == 0 ? 1.2 : 0.8 + 0.1 * l }, std::vector<vec>{ { 1.0 } }, std::array<double, 3>{ at.x, at.y, at.z });
                shells.back().kind = occ::gto::Shell::Kind::Spherical;
                shells.back().incorporate_shell_norm();
            }
        occ::gto::AOBasis basis(atoms, shells, "l" + std::to_string(lmax));
        basis.set_pure(true);
        occ::qm::HartreeFock hf(basis);
        if (lmax <= 5)
        {
            occ::qm::SCF<occ::qm::HartreeFock> scf(hf, occ::qm::SpinorbitalKind::Restricted);
            scf.set_charge_multiplicity(0, 1);
            scf.compute_initial_guess();
            scf.compute_scf_energy();
            return WFN(scf.wavefunction(), false);
        }
        occ::qm::Wavefunction wf;
        wf.basis = basis;
        wf.atoms = atoms;
        wf.nbf = (int)basis.nbf();
        wf.num_electrons = 2;
        wf.mo.n_alpha = wf.mo.n_beta = 1;
        wf.mo.n_ao = wf.nbf;
        wf.mo.C = Eigen::SelfAdjointEigenSolver<occ::Mat>(hf.compute_overlap_matrix()).operatorInverseSqrt();
        wf.mo.energies = occ::Vec::Zero(wf.nbf);
        wf.mo.update_occupied_orbitals();
        wf.mo.update_density_matrix();
        return WFN(wf, false);
    }

    static void expect_orthonormal(const WFN& wave)
    {
        const int nex = wave.get_nex(), nmo = wave.get_nmo();
        const vec2 S = primitive_overlap(wave);
        vec2 SC(nmo, vec(nex));
        for (int i = 0; i < nmo; i++)
            for (int a = 0; a < nex; a++)
                for (int b = 0; b < nex; b++) SC[i][a] += S[a][b] * wave.get_MO_coef(i, b);
        double electrons = 0;
        for (int i = 0; i < nmo; i++)
            for (int j = 0; j <= i; j++)
            {
                double o = 0;
                for (int a = 0; a < nex; a++) o += wave.get_MO_coef(i, a) * SC[j][a];
                EXPECT_NEAR(o, i == j ? 1.0 : 0.0, 1e-8) << "MOs " << i << " " << j;
                if (i == j) electrons += wave.get_MO_occ(i) * o;
            }
        EXPECT_NEAR(electrons, 2.0, 1e-8);
    }

    TEST(OccHighAngularTests, HartreeFockMOsAreOrthonormalToH)
    {
        const WFN wave = occ_h2(5);
        EXPECT_EQ(wave.get_nex(), 2 * 56);
        EXPECT_EQ(wave.get_nmo(), 2 * 36);
        expect_orthonormal(wave);
    }

    TEST(OccHighAngularTests, LoewdinMOsAreOrthonormalToL10)
    {
        const WFN wave = occ_h2(10);
        EXPECT_EQ(wave.get_nex(), 2 * 286);
        EXPECT_EQ(wave.get_nmo(), 2 * 121);
        expect_orthonormal(wave);
    }

    static void expect_same_density(const WFN& a, const WFN& b, const double rtol, const std::string& what)
    {
        ASSERT_EQ(a.get_ncen(), b.get_ncen()) << what;
        for (int i = 0; i < a.get_ncen(); i++)
            for (const double r : { 0.3, 1.1, 2.5 })
            {
                const d3 pos{ a.get_atom_coordinate(i, 0) + r, a.get_atom_coordinate(i, 1) + 0.5 * r, a.get_atom_coordinate(i, 2) - 0.7 * r };
                const double da = a.compute_dens(pos);
                EXPECT_NEAR(da, b.compute_dens(pos), rtol * std::max(1.0, da)) << what << " atom " << i << " r " << r;
            }
    }

    //every format is written through write_wfn and write_wfx and read back
    TEST(FormatConsistencyTests, SameDensityFromEveryFormat)
    {
        const auto root = nos_test_repo_root() / "tests";
        {
            WFN a(root / "cytidine_tonto" / "stdout_cyt", false), b(root / "cytidine_tonto" / "cyt.wfn", false);
            //the H positions in the tonto output carry three decimals, cyt.wfn eight
            expect_same_density(a, b, 1e-3, "cytidine tonto vs its wfn");
        }
        const std::filesystem::path inputs[] = { root / "epoxide_gbw" / "epoxide.gbw", root / "CuF2_i_func" / "71" / "calc.gbw", root / "molden_file" / "Sc_full.molden", root / "molden_file" / "Ce_full.molden", root / "NiP3_fchk" / "good.fchk", root / "alanine_occ" / "alanine.owf.fchk" };
        const auto tmp = std::filesystem::temp_directory_path() / "nosphera2_format_roundtrip.wfn", tmpx = std::filesystem::temp_directory_path() / "nosphera2_format_roundtrip.wfx";
        for (const auto& input : inputs)
        {
            WFN a(input, false);
            ASSERT_TRUE(a.write_wfn(tmp, false, true)) << input;
            ASSERT_TRUE(a.write_wfx(tmpx, true)) << input;
            WFN b(tmp, false), c(tmpx, false);
            expect_same_density(a, b, 1e-6, input.filename().string() + " vs its wfn");
            expect_same_density(a, c, 1e-6, input.filename().string() + " vs its wfx");
            EXPECT_EQ(a.get_nmo(true), c.get_nmo()) << input;
            EXPECT_EQ(a.get_nr_electrons(), c.get_nr_electrons()) << input;
        }
        std::filesystem::remove(tmp);
        std::filesystem::remove(tmpx);
    }

    //the primitive-expanded WFN back into a contracted cartesian fchk through free_fchk and read again: tonto with the
    //basis set from the library file, gbw and molden with the basis the file itself carries
    TEST(FormatConsistencyTests, FchkFullCircle)
    {
        const auto root = nos_test_repo_root() / "tests";
        const auto tmp = std::filesystem::temp_directory_path() / "nosphera2_fchk_full_circle.fchk";
        std::ostringstream log;
        for (const auto& input : { root / "NiP3_fchk" / "in.ffn", root / "epoxide_gbw" / "epoxide.gbw", root / "molden_file" / "Sc_full.molden", root / "NiP3_fchk" / "good.fchk" })
        {
            WFN a(input, false);
            a.set_basis_set_name((root / "NiP3_fchk" / "def2-TZVP").string());
            a.assign_charge(a.calculate_charge());
            ASSERT_TRUE(a.guess_multiplicity(log)) << input;
            ASSERT_TRUE(free_fchk(log, tmp, "", a, false, true)) << input;
            WFN b(tmp, false);
            expect_same_density(a, b, 1e-6, input.filename().string() + " vs its fchk");
            EXPECT_EQ(a.get_nr_electrons(), b.get_nr_electrons()) << input;
        }
        std::filesystem::remove(tmp);
    }

    //accuracy 5 is the largest Becke grid; the fixture with the highest angular momentum of every format has to integrate to its electron count
    TEST(FormatConsistencyTests, ElectronCountOnLargestGrid_full)
    {
        if (const char* env = std::getenv("RUN_FULL_TEST"); !env || std::string(env) == "0" || std::string(env) == "false")
            GTEST_SKIP() << "Set RUN_FULL_TEST=1 to integrate every format on the accuracy 5 grid";
        const auto root = nos_test_repo_root() / "tests";
        const std::filesystem::path inputs[] = { root / "CuF2_i_func" / "71" / "calc.gbw", root / "molden_file" / "Ce_full.molden", root / "molden_file" / "Co2.molden", root / "NiP3_fchk" / "good.fchk", root / "alanine_occ" / "alanine.owf.fchk", root / "sucrose_fchk_SF" / "sucrose.fchk", root / "molden_file" / "Ce_full.wfn", root / "grown" / "water.wfx", root / "cytidine_tonto" / "stdout_cyt", root / "ptb_H_file" / "wfn.xtb" };
        for (const auto& input : inputs)
        {
            WFN wave(input, false);
            wave.delete_unoccupied_MOs();
            double electrons = 0;
            for (int i = 0; i < wave.get_nmo(); i++)
                electrons += wave.get_MO_occ(i);
            GridConfiguration config;
            config.partition_type = PartitionType::Becke;
            config.accuracy = 5;
            GridManager gm(config);
            ivec atom_list(wave.get_ncen());
            std::iota(atom_list.begin(), atom_list.end(), 0);
            gm.setup3DGridsForMolecule(wave, atom_list);
            const PartitionResults res = gm.calculatePartitionedCharges(wave);
            EXPECT_NEAR(res.overall_charges[PartitionResults::S_BECKE], electrons, 1e-4 * electrons) << input;
        }
    }
    //OCC writes fchk in the Gaussian convention, so the reader's pure conventions are checked against the bridge for every l
    TEST(FormatConsistencyTests, OccFchkMatchesBridge)
    {
        spdlog::set_level(spdlog::level::err);
        for (int lmax = 1; lmax <= 4; lmax++)
        {
            const std::vector<occ::core::Atom> atoms{ { 1, 0.0, 0.0, -0.7 }, { 1, 0.0, 0.3, 0.7 } };
            std::vector<occ::gto::Shell> shells;
            for (const auto& at : atoms)
                for (int l = 0; l <= lmax; l++)
                {
                    shells.emplace_back(l, std::vector<double>{ 1.2, 0.5 + 0.1 * l }, std::vector<vec>{ { 0.4, 0.7 } }, std::array<double, 3>{ at.x, at.y, at.z });
                    shells.back().kind = occ::gto::Shell::Kind::Spherical;
                    shells.back().incorporate_shell_norm();
                }
            occ::gto::AOBasis basis(atoms, shells, "test");
            basis.set_pure(true);
            occ::qm::HartreeFock hf(basis);
            occ::qm::SCF<occ::qm::HartreeFock> scf(hf, occ::qm::SpinorbitalKind::Restricted);
            scf.set_charge_multiplicity(0, 1);
            scf.compute_initial_guess();
            scf.compute_scf_energy();
            occ::qm::Wavefunction wf = scf.wavefunction();
            const auto tmp = std::filesystem::temp_directory_path() / ("nosphera2_occ_roundtrip_l" + std::to_string(lmax) + ".fchk");
            {
                occ::io::FchkWriter writer(tmp.string());
                wf.save(writer);
                writer.write();
            }
            expect_same_density(WFN(wf, false), WFN(tmp, false), 1e-6, "occ pure fchk lmax " + std::to_string(lmax));
        }
    }
    TEST(OccHighAngularTests, IntegratesElectronCountToL10)
    {
        for (const int lmax : { 5, 10 })
        {
            WFN wave = occ_h2(lmax);
            wave.delete_unoccupied_MOs();
            GridConfiguration config;
            config.partition_type = PartitionType::Becke;
            config.accuracy = 4;
            GridManager gm(config);
            ivec atom_list{ 0, 1 };
            gm.setup3DGridsForMolecule(wave, atom_list);
            const PartitionResults res = gm.calculatePartitionedCharges(wave);
            EXPECT_NEAR(res.overall_charges[PartitionResults::S_BECKE], 2.0, 1e-4) << "lmax " << lmax;
        }
    }

} // namespace NoSpherA2UnitTests
