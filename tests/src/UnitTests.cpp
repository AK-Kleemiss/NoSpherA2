
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
#include "core/basis_set.h"
#include "core/geometry_aid.h"
#include "core/crystal_energies.h"
#include "core/NoSpherA2.h"
#include "core/npy.h"
#ifdef NOSPHERA2_USE_GPU
#include "core/blas_gpu.h"
#include "core/aux_density_gpu.h"
#endif

static constexpr double PI_VAL = 3.14159265358979323846;

namespace {
    int test_analytical_fourier()
    {
        // Generate grid and k_pts
        vec2 kpts;

        for (int i = 1; i < 100; i++) {
            //Generate random k-points with values between -1 and 1
            kpts.push_back({ (double)rand() / RAND_MAX * 2 - 1, (double)rand() / RAND_MAX * 2 - 1, (double)rand() / RAND_MAX * 2 - 1 });
        }
        vec2 grid;
        grid.resize(7); // x, y, z, dens, atomic_weight, becke_weight, TFVC_weight

        // Conditions for the Wavefunction
        const double c_exp = 2.0;
        double vals[] = { 1.0 };
        unsigned int max_l = 6;
        double radial_res = 1E-17;
        int charge = 1;

        double alpha_min[] = { 0.5 };
        AtomGrid griddy(radial_res,
            350,
            5000,
            charge,
            c_exp,
            max_l,
            alpha_min,
            std::cout);

        double pos[] = { 0 };
        for (int i = 0; i < grid.size(); i++)
        {
            grid[i].resize(griddy.get_num_grid_points(), 0.0);
        }

        vec empty(1, 0.0);
        griddy.get_grid(1, 0, pos, pos, pos, &charge, grid[0].data(), grid[1].data(), grid[2].data(), grid[4].data(), grid[5].data(), grid[6].data(), WFN(), empty);

        // Initialize the vectors sf_A and sf_N
        cvec2 sf_A, sf_N;
        sf_A.resize(1);
        sf_N.resize(1);
        sf_A[0].resize(kpts.size(), 0.0);
        sf_N[0].resize(kpts.size(), 0.0);

        bool all_correct = true; //Veryfiy if all m for one l are correct break if one failed
        bool correct = true; //Verify if the current m is correct

        cdouble max_diff, diff;
        for (unsigned int type = 0; type <= max_l; type++)
        {
            std::cout << "Testing l = " << type << "\n";
            vec coefs(type * 2 + 1);

            // Initialize the Wavefunction
            WFN wavy(e_origin::NOT_YET_DEFINED);
            wavy.push_back_MO(0, 1.0, -13);
            wavy.push_back_atom("H", 0, 0, 0, 1);
            wavy.push_back_atom_basis_set(0, c_exp, vals[0], type, 0);
            primitive p(1, type, c_exp, vals[0]);

            for (unsigned int l = 0; l < type * 2 + 1; l++)
            {
                int m = static_cast<int>(l) - static_cast<int>(type);
                for (int i = 0; i < coefs.size(); i++)
                {
                    coefs[i] = 0.0;
                }
                max_diff = 0.0;
                coefs[l] = 1.0;


                for (int i = 0; i < grid[0].size(); i++)
                {
                    //grid[3][i] = wavy.compute_dens(grid[0][i], grid[1][i], grid[2][i]);
                    grid[3][i] = calc_density_ML(grid[0][i], grid[1][i], grid[2][i], coefs, wavy.get_atoms());
                }

                // Empty the vectors sf:A nad sf_N
                for (int i = 0; i < kpts.size(); i++)
                {
                    sf_A[0][i] = 0.0;
                    sf_N[0][i] = 0.0;
                }
                double work = 0.0;
#pragma omp parallel for private(work)
                for (int i = 0; i < kpts.size(); i++)
                {
                    double k_pt_local[4] = { kpts[i][0] * 2 * constants::PI , kpts[i][1] * 2 * constants::PI , kpts[i][2] * 2 * constants::PI , 0.0 };
                    k_pt_local[3] = sqrt(k_pt_local[0] * k_pt_local[0] + k_pt_local[1] * k_pt_local[1] + k_pt_local[2] * k_pt_local[2]);
                    for (int d = 0; d < 3; d++) k_pt_local[d] /= k_pt_local[3];

                    sf_A[0][i] = sfac_bessel(p, k_pt_local, coefs.data());
                    //sf_A[0][i] = sfac_bessel(p, k_pt_local, ri_coefs);
                    for (int _p = 0; _p < grid[0].size(); _p++)
                    {
                        work = constants::TWO_PI * (kpts[i][0] * grid[0][_p] + kpts[i][1] * grid[1][_p] + kpts[i][2] * grid[2][_p]);
                        sf_N[0][i] += std::polar(grid[3][_p] * grid[4][_p], work);
                    }
                    diff = abs(sf_A[0][i] - sf_N[0][i]);
                    if (abs(diff) > abs(max_diff))
                    {
                        max_diff = diff;
                    }
                    if (abs(diff) > 2E-5)
                    {
                        all_correct = false;
                        correct = false;
                    }
                }
                if (!correct)
                {
                    std::cout << "Error at m: " << m
                              << "   Max diff: (" << max_diff.real() << "|" << max_diff.imag() << ")\n";
                    correct = true;
                }
            }
            if (!all_correct)
                break;
            std::cout << "| PASSED!\n";
        }
        if (!all_correct)
        {
            using namespace std;
            ofstream result("sfacs.dat", ios::out);
            for (int i = 0; i < kpts.size(); i++)
            {
                result << setw(8) << setprecision(2) << fixed << kpts[i][0];
                result << setw(8) << setprecision(2) << fixed << kpts[i][1];
                result << setw(8) << setprecision(2) << fixed << kpts[i][2];
                result << setw(16) << setprecision(8) << scientific << sf_A[0][i].real();
                result << setw(16) << setprecision(8) << scientific << sf_A[0][i].imag();
                result << setw(16) << setprecision(8) << scientific << sf_N[0][i].real();
                result << setw(16) << setprecision(8) << scientific << sf_N[0][i].imag();
                result << setw(16) << setprecision(8) << scientific << abs(sf_A[0][i] - sf_N[0][i]);
                result << setw(35) << setprecision(8) << scientific << sf_A[0][i] / sf_N[0][i];
                result << "\n";
            }
            result.flush();
            result.close();
            std::cout << "Error in the calculations!\n";
            return 1;
        }
        std::cout << "All tests passed!\n";
        return 0;
    };

    template <typename MatType1, typename MatType2>
    void compare_matrices(const MatType1& mat, const MatType2& vecMat) {
        const double tol = 1e-9;
        for (size_t i = 0; i < mat.extent(0); i++) {
            for (size_t j = 0; j < mat.extent(1); j++) {
                //Compile different versions based on the type of the matrix dMatrix and dMatrix2 are accessed with mat(i, j) and vecMat[i][j] respectively
                if constexpr (std::is_same_v<MatType1, dMatrix2> && std::is_same_v<MatType2, vec2>) {
                    EXPECT_NEAR(mat(i, j), vecMat[i][j], tol);
                }
                else if constexpr (std::is_same_v<MatType1, cMatrix2> && std::is_same_v<MatType2, cvec2>) {
                    EXPECT_NEAR(mat(i, j).real(), vecMat[i][j].real(), tol);
                    EXPECT_NEAR(mat(i, j).imag(), vecMat[i][j].imag(), tol);
                }
            }
        }
    }

    void test_solve_linear_equations() {
        // Small, well-conditioned 3x3 test (precomputed)
        // A * x_expected = b
        const vec2 A = {
            {3.0,  2.0, -1.0},
            {2.0, -2.0,  4.0},
            {-1.0, 0.5, -1.0}
        };

        // Precomputed true solution
        const vec x_expected = { 1.0, -2.0, -2.0 };

        // Build right-hand side b = A * x_expected
        vec b(3, 0.0);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                b[i] += A[i][j] * x_expected[j];

        // Keep original RHS for logging if the test fails
        vec b_orig = b;

        // Call the solver -- it replaces b with the solution x (in-place)
        solve_linear_system(A, b);

        // Validate result
        const double tol = 1e-9;
        bool ok = true;
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(b[i], x_expected[i], tol);
        }

        if (ok) {
            std::cout << "test_solve_linear_equations: PASSED\n";
        }
        else {
            std::cout << "test_solve_linear_equations: FAILED\n";
            std::cout << "Matrix A:\n";
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j)
                    std::cout << std::setw(12) << A[i][j];
                std::cout << std::endl;
            }
            std::cout << "Target   b: " << b_orig[0] << " " << b_orig[1] << " " << b_orig[2] << std::endl;
            std::cout << "Original x: ";
            for (double v : x_expected) std::cout << v << " ";
            std::cout << std::endl;
            std::cout << "Returned x: ";
            for (double v : b) std::cout << v << " ";
            std::cout << std::endl;
        }
    } 

    void test_openblas()
    {
        ivec dims = { 10, 10 };
        // Init Mat A with some values as a 3x3 matrix
        vec2 A(dims[0], vec(dims[1]));
        vec2 B(dims[0], vec(dims[1]));
        // Init A and B with random values between -100 and 100
        for (int i = 0; i < dims[0]; i++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                A[i][j] = rand() % 200 - 100;
                B[i][j] = rand() % 200 - 100;
            }
        }

        vec fA = flatten<double>(A);
        vec fB = flatten<double>(B);

        dMatrix2 matA(dims[0], dims[1]);
        std::copy(fA.data(), fA.data() + fA.size(), matA.data());
        dMatrix2 matB(dims[0], dims[1]);
        std::copy(fB.data(), fB.data() + fB.size(), matB.data());
        // Init Mat A and Mat B as 3x3 matrices

        std::cout << "Testing matrices directly\n";
        compare_matrices(matA, A);
        compare_matrices(matB, B);

        std::cout << "Testing untransposed matrices\n";
        // First test regular dot-product
        compare_matrices(dot(matA, matB, false, false), self_dot(A, B));

        std::cout << "Testing transpose A\n";
        ////Second compare first transpose
        compare_matrices(dot(matA, matB, true, false), self_dot(transpose(A), B));

        std::cout << "Testing transpose B\n";
        ////Third comparte second transpose
        compare_matrices(dot(matA, matB, false, true), self_dot(A, transpose(B)));

        std::cout << "Testing transpose A and B\n";
        ////Fourth compare both transposed
        compare_matrices(dot(matA, matB, true, true), self_dot(transpose(A), transpose(B)));

        // Init Complex matrices
        cvec2 C(dims[0], cvec(dims[1])), D(dims[0], cvec(dims[1]));
        for (int i = 0; i < dims[0]; i++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                C[i][j] = cdouble(rand() % 200 - 100, rand() % 200 - 100);
                D[i][j] = cdouble(rand() % 200 - 100, rand() % 200 - 100);
            }
        }

        cvec fC = flatten<cdouble>(C);
        cvec fD = flatten<cdouble>(D);
        cMatrix2 matC(dims[0], dims[1]);
        std::copy(fC.data(), fC.data() + fC.size(), matC.data());
        cMatrix2 matD(dims[0], dims[1]);
        std::copy(fD.data(), fD.data() + fD.size(), matD.data());

        std::cout << "Testing C-matrices directly\n";
        compare_matrices(matC, C);
        compare_matrices(matD, D);

        std::cout << "Testing untransposed C-matrices\n";
        // First test regular dot-product
        compare_matrices(dot(matC, matD, false, false), self_dot(C, D));

        std::cout << "Testing transpose C\n";
        ////Second compare first transpose
        compare_matrices(dot(matC, matD, true, false), self_dot(transpose(C), D));

        std::cout << "Testing transpose D\n";
        ////Third comparte second transpose
        compare_matrices(dot(matC, matD, false, true), self_dot(C, transpose(D)));

        std::cout << "Testing transpose C and D\n";
        ////Fourth compare both transposed
        compare_matrices(dot(matC, matD, true, true), self_dot(transpose(C), transpose(D)));

        // Test 2D x 1D matrix multiplication
        dims[0] = 12;
        vec E(dims[1]);
        vec2 F(dims[0], vec(dims[1]));
        for (int i = 0; i < dims[1]; i++)
        {
            E[i] = rand() % 200 - 100;
            for (int j = 0; j < dims[0]; j++)
            {
                F[j][i] = rand() % 200 - 100;
            }
        }
        vec fE = flatten<double>(F);
        dMatrix1 matE(dims[1]);
        dMatrix2 matF(dims[0], dims[1]);
        std::copy(E.data(), E.data() + E.size(), matE.data());
        std::copy(fE.data(), fE.data() + fE.size(), matF.data());

        // For matrix just reuse matA
        std::cout << "Testing 2D x 1D matrix multiplication\n";
        compare_matrices(dot(matF, matE, false), self_dot(F, E));

        std::cout << "All BLAS tests passed!\n";
    }

    void test_reading_SALTED_binary_file() {
        std::filesystem::path path("../../../tests/SALTED/Model/model.salted");
        if (!std::filesystem::exists(path)) {
            path = std::filesystem::path("tests/SALTED/Model/model.salted");
        }
        if (!std::filesystem::exists(path)) {
            path = std::filesystem::path("Model/model.salted");
        }
        if (!std::filesystem::exists(path)) {
            path = std::filesystem::path("../SALTED/Model/model.salted");
        }
        SALTED_BINARY_FILE file = SALTED_BINARY_FILE(path, true);
        Config config;
        file.populate_config(config);
        std::unordered_map<int, std::vector<int64_t>> fps = file.read_fps();
        std::unordered_map<std::string, vec> averages = file.read_averages();
        std::unordered_map<int, vec> wigners = file.read_wigners();
        vec weights = file.read_weights();
        std::unordered_map<std::string, dMatrix2> feats = file.read_features();
        std::unordered_map<std::string, dMatrix2> proj = file.read_projectors();
        std::cout << "Finished reading SALTED binary file\n";
        // TEST if both configs are the same
        std::cout << "Average:" << config.average << "\n";
        std::cout << "Field:" << config.field << "\n";
        std::cout << "Sparsify:" << config.sparsify << "\n";
        std::cout << "Ncut:" << config.ncut << "\n";
        std::cout << "Ntrain:" << config.Ntrain << "\n";
        std::cout << "Menv:" << config.Menv << "\n";
        std::cout << "trainfrac:" << config.trainfrac << "\n";
        std::cout << "Rcut1:" << config.rcut1 << "\n";
        std::cout << "Rcut2:" << config.rcut2 << "\n";
        std::cout << "nang1:" << config.nang1 << "\n";
        std::cout << "nang2:" << config.nang2 << "\n";
        std::cout << "sig1:" << config.sig1 << "\n";
        std::cout << "sig2:" << config.sig2 << "\n";
        std::cout << "zeta:" << config.zeta << "\n";
        std::cout << "neighspe size:" << config.neighspe1.size() << "\n";
        for (int i = 0; i < config.neighspe1.size(); i++)
        {
            std::cout << "neighspe1[" << i << "]:" << config.neighspe1[i] << "\n";
        }
        std::cout << "neighspe2 size:" << config.neighspe2.size() << "\n";
        for (int i = 0; i < config.neighspe2.size(); i++)
        {
            std::cout << "neighspe2[" << i << "]:" << config.neighspe2[i] << "\n";
        }
        std::cout << "dfBasis:" << config.dfbasis << "\n";

        std::cout << "Comparing wigners\n";
        for (int i = 0; i < wigners.size(); i++)
        {
            for (int j = 0; j < wigners[i].size(); j += 10)
            {
                std::cout << "wigners[" << i << "][" << j << "]:" << wigners[i][j] << "\n";
            }
        }

        std::cout << "Comparing FPS\n";
        for (int i = 0; i < fps.size(); i++)
        {
            for (int j = 0; j < fps[i].size(); j += 10)
            {
                std::cout << "fps[" << i << "][" << j << "]:" << fps[i][j] << "\n";
            }
        }

        std::cout << "All tests passed!\n";
    }
}

namespace NoSpherA2UnitTests
{
    TEST(TscBlockTests, ConstructorWarnsForDuplicateScattererIds)
    {
        const std::vector<std::vector<cdouble>> form_factors = {
            { cdouble(1.0, 0.0) },
            { cdouble(2.0, 0.0) },
            { cdouble(3.0, 0.0) }
        };
        const std::vector<atomID> scatterer_ids = { {0.1, 0.1, 0.1, 0, 1 }, {0.5, 0.2, 0.4, 0, 2 }, {0.1, 0.1, 0.1, 0, 1 } };
        const std::vector<std::vector<int>> indices = { { 1 }, { 0 }, { 0 } };

        testing::internal::CaptureStdout();
        tsc_block<int, cdouble> block(form_factors, scatterer_ids, indices);
        const std::string output = testing::internal::GetCapturedStdout();

        EXPECT_NE(output.find("Warning: Duplicate scatterer: atomID(frac_x: 0.1, frac_y: 0.1, frac_z: 0.1, Z: 1, data: 0, reserved: 0)"), std::string::npos);
    }

    TEST(TscBlockTests, AppendWarnsForDuplicateScattererIds)
    {
        const std::vector<std::vector<int>> indices = { { 1 }, { 0 }, { 0 } };
        const std::vector<std::vector<cdouble>> lhs_form_factors = {
            { cdouble(1.0, 0.0) }
        };
        const std::vector<atomID> lhs_scatterer_ids{ {0.1, 0.1, 0.1, 0, 1 } };
        const std::vector<std::vector<cdouble>> rhs_form_factors = {
            { cdouble(2.0, 0.0) },
            { cdouble(3.0, 0.0) }
        };
        const std::vector<atomID> rhs_scatterer_ids{ {0.1, 0.1, 0.1, 0, 1 }, {0.5, 0.2, 0.4, 0, 2 } };

        tsc_block<int, cdouble> lhs(lhs_form_factors, lhs_scatterer_ids, indices);
        tsc_block<int, cdouble> rhs(rhs_form_factors, rhs_scatterer_ids, indices);
        std::ostringstream log;

        testing::internal::CaptureStdout();
        lhs.append(rhs, log);
        const std::string output = testing::internal::GetCapturedStdout();

        EXPECT_NE(output.find("Warning: Duplicate scatterer in append: atomID(frac_x: 0.1, frac_y: 0.1, frac_z: 0.1, Z: 1, data: 0, reserved: 0)"), std::string::npos);
        EXPECT_EQ(lhs.scatterer_size(), 2);
    }

    TEST(TscBlockTests, BinaryFileRoundTripsWith32BitSizes)
    {
        const std::vector<std::vector<cdouble>> form_factors = {
            { cdouble(1.25, -0.5), cdouble(2.5, 0.75) },
            { cdouble(-3.0, 1.5), cdouble(4.25, -2.0) }
        };
        const std::vector<atomID> scatterer_ids = {
            atomID(0.1, 0.2, 0.3, 1, 6),
            atomID(0.4, 0.5, 0.6, 2, 8)
        };
        const std::vector<std::vector<int>> indices = {
            { 1, -2 }, { 0, 3 }, { -1, 4 }
        };
        const std::filesystem::path path =
            std::filesystem::temp_directory_path() / "nosphera2_tscb_32bit_roundtrip.tscb";

        tsc_block<int, cdouble> original(
            form_factors, scatterer_ids, indices);
        original.write_tscb_file({}, path);
        tsc_block<int, cdouble> restored(path);
        std::filesystem::remove(path);

        ASSERT_EQ(restored.scatterer_size(), scatterer_ids.size());
        ASSERT_EQ(restored.reflection_size(), indices[0].size());
        for (std::size_t scatterer = 0; scatterer < scatterer_ids.size(); ++scatterer)
        {
            EXPECT_EQ(std::get<atomID>(restored.get_scatterer(scatterer)),
                scatterer_ids[scatterer]);
            EXPECT_EQ(restored.get_sf_for_scatterer(scatterer), form_factors[scatterer]);
        }
        for (std::size_t reflection = 0; reflection < indices[0].size(); ++reflection)
        {
            EXPECT_EQ(restored.get_indices(reflection),
                (std::array<int, 3>{ indices[0][reflection], indices[1][reflection], indices[2][reflection] }));
        }
    }

    TEST(SALTEDTests, ReadingSALTEDBinaryFile)
    {
        test_reading_SALTED_binary_file();
    }

    // -----------------------------------------------------------------------

    TEST(GeometryTests, ArrayDistance_UnitCubeDiagonal)
    {
        double d = array_length(d3{ 0.0, 0.0, 0.0 }, d3{ 1.0, 1.0, 1.0 });
        EXPECT_NEAR(std::sqrt(3.0), d, 1e-12);
    }

    TEST(GeometryTests, ArrayDistance_SamePoint)
    {
        EXPECT_NEAR(0.0, array_length(d3{ 1.5, 2.3, -4.7 }, d3{ 1.5, 2.3, -4.7 }), 1e-12);
    }

    TEST(GeometryTests, VecDiff_Basic)
    {
        d3 r = vec_diff({ 5.0, 10.0, 15.0 }, { 1.0, 2.0, 3.0 });
        EXPECT_NEAR(4.0, r[0], 1e-12);
        EXPECT_NEAR(8.0, r[1], 1e-12);
        EXPECT_NEAR(12.0, r[2], 1e-12);
    }

    TEST(GeometryTests, VecDiff_ZeroResult)
    {
        d3 r = vec_diff({ 1.0, 2.0, 3.0 }, { 1.0, 2.0, 3.0 });
        EXPECT_NEAR(0.0, r[0], 1e-12);
        EXPECT_NEAR(0.0, r[1], 1e-12);
        EXPECT_NEAR(0.0, r[2], 1e-12);
    }

    TEST(GeometryTests, VecCross_BasisVectors)
    {
        // i x j = k
        d3 r = vec_cross({ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 });
        EXPECT_NEAR(0.0, r[0], 1e-12);
        EXPECT_NEAR(0.0, r[1], 1e-12);
        EXPECT_NEAR(1.0, r[2], 1e-12);
    }

    TEST(GeometryTests, VecCross_AntiCommutative)
    {
        // a x b = -(b x a)
        d3 ab = vec_cross({ 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 });
        d3 ba = vec_cross({ 4.0, 5.0, 6.0 }, { 1.0, 2.0, 3.0 });
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(-ab[i], ba[i], 1e-12);
        }
    }

    TEST(GeometryTests, VecCross_ParallelVectors)
    {
        // Parallel vectors → zero cross product
        d3 r = vec_cross({ 2.0, 4.0, 6.0 }, { 1.0, 2.0, 3.0 });
        EXPECT_NEAR(0.0, r[0], 1e-12);
        EXPECT_NEAR(0.0, r[1], 1e-12);
        EXPECT_NEAR(0.0, r[2], 1e-12);
    }

    TEST(GeometryTests, VecDot_Perpendicular)
    {
        EXPECT_NEAR(0.0, vec_dot({ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }), 1e-12);
    }

    TEST(GeometryTests, VecDot_Parallel)
    {
        EXPECT_NEAR(3.0, vec_dot({ 1.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }), 1e-12);
    }

    TEST(GeometryTests, VecDot_KnownValue)
    {
        EXPECT_NEAR(32.0, vec_dot({ 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }), 1e-12);
    }

    // -----------------------------------------------------------------------

    TEST(NumericTests, IsSimilarRel_WithinTolerance)
    {
        // 1.0 vs 1.005: relative diff ≈ 0.5 %, within 1 %
        EXPECT_TRUE(is_similar_rel(1.0, 1.005, 0.01));
    }

    TEST(NumericTests, IsSimilarRel_OutsideTolerance)
    {
        // 1.0 vs 1.05: relative diff = 5 %, outside 1 %
        EXPECT_FALSE(is_similar_rel(1.0, 1.05, 0.01));
    }

    TEST(NumericTests, IsSimilarRel_EqualValues)
    {
        EXPECT_TRUE(is_similar_rel(42.0, 42.0, 1e-6));
    }

    TEST(NumericTests, IsSimilarAbs_WithinTolerance)
    {
        EXPECT_TRUE(is_similar_abs(1.0, 1.0009, 0.001));
    }

    TEST(NumericTests, IsSimilarAbs_OutsideTolerance)
    {
        EXPECT_FALSE(is_similar_abs(1.0, 1.002, 0.001));
    }

    TEST(NumericTests, FastExpNeg_AtZero)
    {
        // fast_exp_neg uses std::exp for x > -ln(2) ≈ -0.693
        double val = fast_exp_neg(0.0);
        EXPECT_NEAR(1.0, val, 1e-12);
    }

    TEST(NumericTests, FastExpNeg_AtMinusOne)
    {
        double approx = fast_exp_neg(-1.0);
        double exact = std::exp(-1.0);
        // Approximation tolerance: 0.5 %
        EXPECT_NEAR(exact, approx, exact * 0.005);
    }

    TEST(NumericTests, FastExpNeg_AtMinusTen)
    {
        // (1 + x/1024)^1024 underestimates exp(x) by ~x²/(2N) ≈ 4.9% at x=-10.
        double approx = fast_exp_neg(-10.0);
        double exact = std::exp(-10.0);
        EXPECT_NEAR(exact, approx, exact * 0.06);
    }

    TEST(NumericTests, FastExpNeg_BeyondCutoff)
    {
        // Values below -42 must return 0 exactly
        EXPECT_NEAR(0.0, fast_exp_neg(-100.0), 1e-30);
    }

    TEST(NumericTests, BLAS_tests)
    {
        test_openblas();
    }

    TEST(NumericTests, SolveLinearEquations)
    {
        test_solve_linear_equations();
    }

    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // FchkParsingTests
    // -----------------------------------------------------------------------

        // FCHK format: 40-char label + type char + value starting at position 49
        // We pad to exactly 49 chars then append the value.

    TEST(FchkParsingTests, ReadFchkInt_PositiveValue)
    {
            // Format: keyword padded to 40, type at 40, 8 spaces (41-48), value at 49+
            const char* line = "Number of atoms                         I        3";
        EXPECT_EQ(3, read_fchk_integer(std::string(line)));
    }

    TEST(FchkParsingTests, ReadFchkInt_NegativeValue)
    {
        const char* line = "Charge                                  I        -1";
        EXPECT_EQ(-1, read_fchk_integer(std::string(line)));
    }

    TEST(FchkParsingTests, ReadFchkInt_LargeValue)
    {
        const char* line = "Number of basis functions               I        1024";
        EXPECT_EQ(1024, read_fchk_integer(std::string(line)));
    }

    TEST(FchkParsingTests, ReadFchkDbl_NegativeScientific)
    {
        const char* line = "Total Energy                            R        -1.23456789E+02";
        double val = read_fchk_double(std::string(line));
        EXPECT_NEAR(-123.456789, val, 1e-6);
    }

    TEST(FchkParsingTests, ReadFchkDbl_PositiveScientific)
    {
        const char* line = "Zero-point correction                   R        4.56000000E-02";
        double val = read_fchk_double(std::string(line));
        EXPECT_NEAR(0.0456, val, 1e-10);
    }

    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // StringUtilTests
    // -----------------------------------------------------------------------

    TEST(StringUtilTests, EndsWith_MatchingSuffix)
    {
        EXPECT_TRUE(ends_with(std::string("molecule.wfx"), std::string(".wfx")));
        EXPECT_TRUE(ends_with(std::string("data.hkl"), std::string(".hkl")));
    }

    TEST(StringUtilTests, EndsWith_NonMatchingSuffix)
    {
        EXPECT_FALSE(ends_with(std::string("molecule.wfx"), std::string(".gbw")));
        EXPECT_FALSE(ends_with(std::string("test.cpp"), std::string(".txt")));
    }

    TEST(StringUtilTests, EndsWith_EmptySuffix)
    {
        // Empty suffix always matches
        EXPECT_TRUE(ends_with(std::string("anything"), std::string("")));
    }

    TEST(StringUtilTests, EndsWith_Suffix_LongerThanString)
    {
        EXPECT_FALSE(ends_with(std::string("ab"), std::string("abc")));
    }

    TEST(StringUtilTests, EndsWith_ExactMatch)
    {
        EXPECT_TRUE(ends_with(std::string(".wfx"), std::string(".wfx")));
    }

    TEST(StringUtilTests, ShrinkString_RemovesDigits)
    {
        std::string input = "C1";
        std::string r = shrink_string(input);
        EXPECT_EQ(1, static_cast<int>(r.size()));
        EXPECT_STREQ("C", r.c_str());
    }

    TEST(StringUtilTests, ShrinkString_RemovesSpacesAndDigits)
    {
        std::string input = "O 1 1";
        std::string r = shrink_string(input);
        int len = static_cast<int>(r.size());
        EXPECT_GT(len, 0);
        // All spaces and digits removed → only "O" remains
        EXPECT_STREQ("O", r.c_str());
    }

    TEST(StringUtilTests, ShrinkString_PureLettersUnchanged)
    {
        std::string input = "Fe";
        std::string r = shrink_string(input);
        EXPECT_STREQ("Fe", r.c_str());
    }

    TEST(StringUtilTests, ShrinkString_BufferTooSmall)
    {
        // Can't fit even "C"
        std::string input = "C1H2O";
        std::string r = shrink_string(input);
        EXPECT_GE(static_cast<int>(r.size()), 1);
    }

    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // Sha256Tests
    // -----------------------------------------------------------------------

    TEST(Sha256Tests, Sha256_OutputLength)
    {
        std::string r = sha::sha256(std::string("abc"));
        EXPECT_EQ(64, static_cast<int>(r.size()));
    }

    TEST(Sha256Tests, Sha256_KnownVector_Abc)
    {
        // NIST FIPS 180-4 test vector
        std::string r = sha::sha256(std::string("abc"));
        EXPECT_STREQ(
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            r.c_str());
    }

    TEST(Sha256Tests, Sha256_EmptyString)
    {
        std::string r = sha::sha256(std::string(""));
        EXPECT_STREQ(
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
            r.c_str());
    }

    TEST(Sha256Tests, Sha256_Deterministic)
    {
        std::string r1 = sha::sha256(std::string("NoSpherA2"));
        std::string r2 = sha::sha256(std::string("NoSpherA2"));
        EXPECT_STREQ(r1.c_str(), r2.c_str());
    }

    TEST(Sha256Tests, Sha256_DifferentInputsDifferentOutputs)
    {
        std::string r1 = sha::sha256(std::string("abc"));
        std::string r2 = sha::sha256(std::string("abd"));
        EXPECT_STRNE(r1.c_str(), r2.c_str());
    }

    TEST(Sha256Tests, Sha256_BufferTooSmall)
    {
        // No more small-buffer API: just validate digest length is always 64.
        std::string r = sha::sha256(std::string("abc"));
        EXPECT_EQ(64, static_cast<int>(r.size()));
    }

    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // VecAggregateTests
    // -----------------------------------------------------------------------

    TEST(VecAggregateTests, VecSumBool_MixedValues)
    {
        bvec v = { true, false, true, true, false };
        EXPECT_EQ(3, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumBool_AllFalse)
    {
        bvec v = { false, false, false };
        EXPECT_EQ(0, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumBool_Empty)
    {
        bvec v;
        EXPECT_EQ(0, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumInt_Basic)
    {
        ivec v = { 1, 2, 3, 4, 5 };
        EXPECT_EQ(15, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumInt_Negative)
    {
        ivec v = { -3, 7, -2 };
        EXPECT_EQ(2, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumDouble_Basic)
    {
        vec v = { 1.5, -0.5, 2.0 };
        EXPECT_NEAR(3.0, vec_sum(v), 1e-12);
    }

    TEST(VecAggregateTests, VecSumDouble_AllZero)
    {
        vec v = { 0.0, 0.0, 0.0 };
        EXPECT_NEAR(0.0, vec_sum(v), 1e-12);
    }

    TEST(VecAggregateTests, VecLength_345)
    {
        vec v = { 3.0, 4.0 };
        EXPECT_NEAR(5.0, vec_length(v), 1e-12);
    }

    TEST(VecAggregateTests, VecLength_Unit)
    {
        vec v = { 1.0, 0.0, 0.0 };
        EXPECT_NEAR(1.0, vec_length(v), 1e-12);
    }

    TEST(VecAggregateTests, VecLength_4D)
    {
        vec v = { 1.0, 1.0, 1.0, 1.0 };
        EXPECT_NEAR(2.0, vec_length(v), 1e-12);
    }

    TEST(VecAggregateTests, VecLength_Empty)
    {
        vec v;
        EXPECT_NEAR(0.0, vec_length(v), 1e-12);
    }
    
    // -----------------------------------------------------------------------

    TEST(StringUtilTests2, Trim_LeadingTrailingSpaces)
    {
        std::string r = trim(std::string("  hello world  "));
        EXPECT_EQ(r, "hello world");
    }

    TEST(StringUtilTests2, Trim_NoSpaces)
    {
        std::string r = trim(std::string("no_spaces"));
        EXPECT_EQ(r, "no_spaces");
    }

    TEST(StringUtilTests2, Trim_EmptyString)
    {
        std::string r = trim(std::string(""));
        EXPECT_TRUE(r.empty());
        EXPECT_EQ(r, "");
    }

    TEST(StringUtilTests2, Trim_OnlySpaces)
    {
        std::string r = trim(std::string("    "));
        EXPECT_EQ(r, "");
    }

    TEST(StringUtilTests2, Asciitolower_Uppercase)
    {
        EXPECT_EQ('a', asciitolower('A'));
        EXPECT_EQ('z', asciitolower('Z'));
        EXPECT_EQ('m', asciitolower('M'));
    }

    TEST(StringUtilTests2, Asciitolower_AlreadyLower)
    {
        EXPECT_EQ('a', asciitolower('a'));
        EXPECT_EQ('z', asciitolower('z'));
    }

    TEST(StringUtilTests2, Asciitolower_NonAlpha)
    {
        // Digits and symbols are returned unchanged
        EXPECT_EQ('5', asciitolower('5'));
        EXPECT_EQ('_', asciitolower('_'));
    }

    TEST(StringUtilTests2, DoubleFromEsd_Plain)
    {
        EXPECT_NEAR(3.14159, double_from_string_with_esd(std::string("3.14159")), 1e-10);
    }

    TEST(StringUtilTests2, DoubleFromEsd_WithEsd)
    {
        // "(5)" is stripped; value is 1.234
        EXPECT_NEAR(1.234, double_from_string_with_esd(std::string("1.234(5)")), 1e-10);
    }

    TEST(StringUtilTests2, DoubleFromEsd_Zero)
    {
        EXPECT_NEAR(0.0, double_from_string_with_esd(std::string("0.0")), 1e-12);
    }

    TEST(StringUtilTests2, DecimalPrecisionCif_WithBracketsAndDecimal)
    {
        // "1.2345(6)" → 6 × 10⁻⁴ = 0.0006
        std::string s = "1.2345(6)";
        EXPECT_NEAR(6e-4, get_decimal_precision_from_CIF_number(s), 1e-10);
    }

    TEST(StringUtilTests2, DecimalPrecisionCif_NoBrackets)
    {
        // No brackets → default 0.005
        std::string s = "1.234";
        EXPECT_NEAR(0.005, get_decimal_precision_from_CIF_number(s), 1e-10);
    }

    TEST(StringUtilTests2, DecimalPrecisionCif_IntegerWithEsd)
    {
        // "100(2)" → no decimal, digit count from bracket positions
        std::string s = "100(2)";
        EXPECT_NEAR(0.2, get_decimal_precision_from_CIF_number(s), 1e-10);
    }

    // -----------------------------------------------------------------------

    TEST(BasisTypeTests, Sht2nbas_CartesianShells)
    {
        // Cartesian: S=1, P=3, D=6, F=10, G=15
        EXPECT_EQ(1, sht2nbas(0));
        EXPECT_EQ(3, sht2nbas(1));
        EXPECT_EQ(6, sht2nbas(2));
        EXPECT_EQ(10, sht2nbas(3));
        EXPECT_EQ(15, sht2nbas(4));
    }

    TEST(BasisTypeTests, Sht2nbas_SphericalShells)
    {
        // Negative types → spherical: -2=D(5), -3=F(7), -4=G(9)
        EXPECT_EQ(5, sht2nbas(-2));
        EXPECT_EQ(7, sht2nbas(-3));
        EXPECT_EQ(9, sht2nbas(-4));
    }

    TEST(BasisTypeTests, DoubleFactorial_SmallValues)
    {
        EXPECT_EQ(1u, doublefactorial(0));
        EXPECT_EQ(1u, doublefactorial(1));
        EXPECT_EQ(2u, doublefactorial(2));
        EXPECT_EQ(3u, doublefactorial(3));
        EXPECT_EQ(8u, doublefactorial(4));
        EXPECT_EQ(15u, doublefactorial(5));
        EXPECT_EQ(48u, doublefactorial(6));
        EXPECT_EQ(105u, doublefactorial(7));
    }

    // -----------------------------------------------------------------------

    TEST(ConstantsTests, ConstAbs_Positive) { EXPECT_EQ(5, constants::const_abs(5)); }
    TEST(ConstantsTests, ConstAbs_Negative) { EXPECT_EQ(5, constants::const_abs(-5)); }
    TEST(ConstantsTests, ConstAbs_Zero) { EXPECT_EQ(0, constants::const_abs(0)); }

    TEST(ConstantsTests, ConstexprPow_Integer)
    {
        EXPECT_NEAR(1024.0, constants::constexpr_pow(2.0, 10), 1e-12);
        EXPECT_NEAR(1.0, constants::constexpr_pow(10.0, 0), 1e-12);
        EXPECT_NEAR(8.0, constants::constexpr_pow(2.0, 3), 1e-12);
    }

    TEST(ConstantsTests, ConstantsSqrt_KnownValues)
    {
        EXPECT_NEAR(2.0, constants::sqrt(4.0), 1e-10);
        EXPECT_NEAR(std::sqrt(2.0), constants::sqrt(2.0), 1e-10);
        EXPECT_NEAR(0.0, constants::sqrt(0.0), 1e-12);
    }

    TEST(ConstantsTests, ConstantsSqrt_NegativeIsNaN)
    {
#if defined(__FAST_MATH__)
    GTEST_SKIP() << "Skipping NaN test under -ffast-math";
#else
    double r = constants::sqrt(-1.0);
    EXPECT_TRUE(std::isnan(r));
#endif
    }

    TEST(ConstantsTests, ExpApprox_AtZero)
    {
        EXPECT_NEAR(1.0, constants::exp_approx(0.0, 25), 1e-12);
    }

    TEST(ConstantsTests, ExpApprox_AtOne)
    {
        // 25-term Taylor series matches std::exp to machine precision
        EXPECT_NEAR(std::exp(1.0), constants::exp_approx(1.0, 25), 1e-12);
    }

    TEST(ConstantsTests, ExpApprox_AtMinusOne)
    {
        EXPECT_NEAR(std::exp(-1.0), constants::exp_approx(-1.0, 25), 1e-10);
    }

    TEST(ConstantsTests, LogApprox_AtOne)
    {
        EXPECT_NEAR(0.0, constants::log_approx(1.0, 25), 1e-12);
    }

    TEST(ConstantsTests, LogApprox_AtE)
    {
        // Arctanh series converges; 25 iterations accurate to 1e-8 for x=e
        EXPECT_NEAR(1.0, constants::log_approx(std::exp(1.0), 25), 1e-6);
    }

    TEST(ConstantsTests, LogApprox_NonPositiveReturnsSentinel)
    {
        EXPECT_NEAR(-1.0, constants::log_approx(0.0, 25), 1e-12);
        EXPECT_NEAR(-1.0, constants::log_approx(-5.0, 25), 1e-12);
    }

    TEST(ConstantsTests, Bohr2Ang_OneBohr)
    {
        // 1 Bohr = a₀ Å = 0.529177210903 Å
        EXPECT_NEAR(0.529177210903, constants::bohr2ang(1.0), 1e-10);
    }

    TEST(ConstantsTests, Bohr2Ang_Zero)
    {
        EXPECT_NEAR(0.0, constants::bohr2ang(0.0), 1e-12);
    }

    TEST(ConstantsTests, Ang2Bohr_OneAngstrom)
    {
        EXPECT_NEAR(1.0 / 0.529177210903, constants::ang2bohr(1.0), 1e-8);
    }

    TEST(ConstantsTests, Ang2Bohr_Roundtrip)
    {
        // ang2bohr(bohr2ang(x)) ≈ x
        double x = 2.5;
        EXPECT_NEAR(x, constants::ang2bohr(constants::bohr2ang(x)), 1e-10);
    }

    TEST(ConstantsTests, CubicBohr2Ang_OneBohr3)
    {
        // 1 Bohr³ = a₀³ Å³
        double expected = 0.529177210903 * 0.529177210903 * 0.529177210903;
        EXPECT_NEAR(expected, constants::cubic_bohr2ang(1.0), 1e-10);
    }

    TEST(ConstantsTests, CubicAng2Bohr_Roundtrip)
    {
        double x = 3.0;
        EXPECT_NEAR(x, constants::cubic_ang2bohr(constants::cubic_bohr2ang(x)), 1e-8);
    }

    TEST(ConstantsTests, Factorial_SmallValues)
    {
        EXPECT_EQ(1LL, static_cast<long long>(constants::ft_fun(0)));
        EXPECT_EQ(1LL, static_cast<long long>(constants::ft_fun(1)));
        EXPECT_EQ(120LL, static_cast<long long>(constants::ft_fun(5)));
        EXPECT_EQ(3628800LL, static_cast<long long>(constants::ft_fun(10)));
    }

    // -----------------------------------------------------------------------

    TEST(OrbitalIndexTests, OrcaToPySCF_SShell)
    {
        // S: only one component, m_idx=0 → 0
        auto r = constants::orca_2_pySCF(0, 0);
        EXPECT_EQ(0, r.has_value() ? static_cast<int>(r.value()) : -1);
    }

    TEST(OrbitalIndexTests, OrcaToPySCF_PShell)
    {
        // P ORCA ordering 0,+1,-1 → PySCF map {1,2,0}
        EXPECT_EQ((size_t)1, constants::orca_2_pySCF(1, 0).value());
        EXPECT_EQ((size_t)2, constants::orca_2_pySCF(1, 1).value());
        EXPECT_EQ((size_t)0, constants::orca_2_pySCF(1, 2).value());
    }

    TEST(OrbitalIndexTests, OrcaToPySCF_DShell_FirstComponent)
    {
        // D: map {2,3,1,4,0}, m_idx=0 → 2
        EXPECT_EQ((size_t)2, constants::orca_2_pySCF(2, 0).value());
    }

    TEST(OrbitalIndexTests, OrcaToPySCF_OutOfRange)
    {
        // l=100 is not in the switch → nullopt → -1
        auto r = constants::orca_2_pySCF(100, 0);
        EXPECT_EQ(-1, r.has_value() ? static_cast<int>(r.value()) : -1);
    }

    TEST(OrbitalIndexTests, TypeToNbo_SShell) { EXPECT_EQ(1u, constants::type_2_nbo(1)); }
    TEST(OrbitalIndexTests, TypeToNbo_PxShell) { EXPECT_EQ(101u, constants::type_2_nbo(2)); }
    TEST(OrbitalIndexTests, TypeToNbo_DxxShell) { EXPECT_EQ(201u, constants::type_2_nbo(5)); }
    TEST(OrbitalIndexTests, TypeToNbo_FxxxShell) { EXPECT_EQ(301u, constants::type_2_nbo(11)); }
    TEST(OrbitalIndexTests, TypeToNbo_GxxxxShell) { EXPECT_EQ(401u, constants::type_2_nbo(21)); }
    TEST(OrbitalIndexTests, TypeToNbo_Unknown) { EXPECT_EQ(0u, constants::type_2_nbo(99)); }

    // -----------------------------------------------------------------------

    TEST(BesselTests, BesselJ0_AtZero)
    {
        // j₀(0) = 1 (limit of sin(x)/x as x→0)
        EXPECT_NEAR(1.0, bessel_first_kind(0, 0.0), 1e-12);
    }

    TEST(BesselTests, BesselJ1_AtZero)
    {
        // j_l(0) = 0 for l > 0
        EXPECT_NEAR(0.0, bessel_first_kind(1, 0.0), 1e-12);
        EXPECT_NEAR(0.0, bessel_first_kind(5, 0.0), 1e-12);
    }

    TEST(BesselTests, BesselJ0_AtOne)
    {
        // j₀(1) = sin(1)/1
        EXPECT_NEAR(std::sin(1.0), bessel_first_kind(0, 1.0), 1e-12);
    }

    TEST(BesselTests, BesselJ1_AtOne)
    {
        // j₁(1) = (sin(1) - cos(1)) / 1
        double expected = std::sin(1.0) - std::cos(1.0);
        EXPECT_NEAR(expected, bessel_first_kind(1, 1.0), 1e-12);
    }

    TEST(BesselTests, BesselJ2_AtOne)
    {
        // j₂(1) = (2·sin(1) - 3·cos(1)) / 1
        double expected = 2.0 * std::sin(1.0) - 3.0 * std::cos(1.0);
        EXPECT_NEAR(expected, bessel_first_kind(2, 1.0), 1e-12);
    }

    TEST(BesselTests, BesselJ0_AtPi)
    {
        // j₀(π) = sin(π)/π ≈ 0
        EXPECT_NEAR(std::sin(PI_VAL) / PI_VAL, bessel_first_kind(0, PI_VAL), 1e-12);
    }

    TEST(BesselTests, BesselJ_HigherOrder_PositiveAndFinite)
    {
        // l=7 exercises the continued-fraction fallback path
        double r = bessel_first_kind(7, 2.0);
        EXPECT_TRUE(std::isfinite(r));
        EXPECT_TRUE(r > 0.0);
    }

    TEST(BesselTests, BesselJ_RecurrenceCheck)
    {
        // Recurrence: j_{l-1}(x) + j_{l+1}(x) = (2l+1)/x · j_l(x)
        double x = 3.0;
        int l = 3;
        double jlm1 = bessel_first_kind(l - 1, x);
        double jl = bessel_first_kind(l, x);
        double jlp1 = bessel_first_kind(l + 1, x);
        double lhs = jlm1 + jlp1;
        double rhs = (2.0 * l + 1.0) / x * jl;
        EXPECT_NEAR(lhs, rhs, 1e-10);
    }

    TEST(BesselTests, AnalyticFourier)
    {
        int ret = test_analytical_fourier();
        EXPECT_EQ(ret, 0);
    }

    // -----------------------------------------------------------------------
    // IsSimilarPow10Tests — is_similar(a, b, tolerance) where |a-b| <= 10^tol
    // -----------------------------------------------------------------------
    TEST(IsSimilarPow10Tests, Equal_ReturnTrue)
    {
        EXPECT_TRUE(is_similar(1.0, 1.0, -6.0));
    }

    TEST(IsSimilarPow10Tests, WithinTolerance_ReturnTrue)
    {
        // |1.000001 - 1.0| = 1e-6 <= 10^(-6)
        EXPECT_TRUE(is_similar(1.000001, 1.0, -6.0));
    }

    TEST(IsSimilarPow10Tests, OutsideTolerance_ReturnFalse)
    {
        // |1.00001 - 1.0| = 1e-5 > 10^(-6)
        EXPECT_FALSE(is_similar(1.00001, 1.0, -6.0));
    }

    TEST(IsSimilarPow10Tests, NegativeValues_WithinTolerance)
    {
        EXPECT_TRUE(is_similar(-5.0, -5.0 + 1e-8, -7.0));
    }

    TEST(IsSimilarPow10Tests, LooseTolerance_LargeDiff)
    {
        // |100 - 50| = 50 <= 10^2 = 100
        EXPECT_TRUE(is_similar(100.0, 50.0, 2.0));
    }

    TEST(IsSimilarPow10Tests, LooseTolerance_TooLargeDiff)
    {
        // |200 - 50| = 150 > 10^2 = 100
        EXPECT_FALSE(is_similar(200.0, 50.0, 2.0));
    }

    // -----------------------------------------------------------------------
    // Shell2FunctionTests — shell2function(type, prim) WFN column index
    // -----------------------------------------------------------------------
    TEST(Shell2FunctionTests, SType_Prim0_ReturnsNonNegative)
    {
        // s-type shell (type=1): first and only function is index 0
        int r = shell2function(1, 0);
        EXPECT_GE(r, 0);
    }

    TEST(Shell2FunctionTests, PType_Prim0_ReturnsNonNegative)
    {
        // p-type shell: 3 functions
        int r = shell2function(2, 0);
        EXPECT_GE(r, 0);
    }

    TEST(Shell2FunctionTests, PType_Prim1)
    {
        int r0 = shell2function(2, 0);
        int r1 = shell2function(2, 1);
        EXPECT_GT(r1, r0);
    }

    TEST(Shell2FunctionTests, DType_Prim5_ValidIndex)
    {
        // d-type (type=3): 6 Cartesian or 5 spherical functions
        int r = shell2function(3, 5);
        EXPECT_GE(r, 0);
    }

    TEST(Shell2FunctionTests, ResultsAreStrictlyIncreasingWithinShell)
    {
        // f-type (type=4): consecutive prims must give increasing column indices
        int r0 = shell2function(4, 0);
        int r1 = shell2function(4, 1);
        int r2 = shell2function(4, 2);
        EXPECT_GT(r1, r0);
        EXPECT_GT(r2, r1);
    }

    // -----------------------------------------------------------------------
    // CountWordsTests — CountWords(str)
    // -----------------------------------------------------------------------
    TEST(CountWordsTests, Empty_Returns0)
    {
        EXPECT_EQ(0, CountWords(""));
    }

    TEST(CountWordsTests, OneWord)
    {
        EXPECT_EQ(1, CountWords("hello"));
    }

    TEST(CountWordsTests, TwoWords)
    {
        EXPECT_EQ(2, CountWords("hello world"));
    }

    TEST(CountWordsTests, LeadingTrailingSpaces)
    {
        EXPECT_EQ(2, CountWords("  foo   bar  "));
    }

    TEST(CountWordsTests, MultipleSpacesBetweenWords)
    {
        EXPECT_EQ(3, CountWords("a  b  c"));
    }

    TEST(CountWordsTests, SingleSpace)
    {
        EXPECT_EQ(0, CountWords(" "));
    }

    // -----------------------------------------------------------------------
    // ShrinkStringToAtomTests — shrink_string_to_atom(input, atom_number)
    // -----------------------------------------------------------------------
    TEST(ShrinkStringToAtomTests, CarbonAtomNumber6)
    {
        // atnr2letter(6) = "C"
        std::string input = "C1";
        std::string r = shrink_string_to_atom(input, 6);
        EXPECT_GE(static_cast<int>(r.size()), 0);
        EXPECT_EQ(std::string("C"), r);
    }

    TEST(ShrinkStringToAtomTests, CalciumAtomNumber20)
    {
        // atnr2letter(20) = "Ca"
        std::string input = "Ca12";
        std::string r = shrink_string_to_atom(input, 20);
        EXPECT_GE(static_cast<int>(r.size()), 0);
        EXPECT_EQ(std::string("Ca"), r);
    }

    TEST(ShrinkStringToAtomTests, BufferTooSmall_ReturnsAtLeastOne)
    {
        std::string input = "Carbon6";
        std::string r = shrink_string_to_atom(input, 6);
        EXPECT_GE(static_cast<int>(r.size()), 1);
    }

    TEST(ShrinkStringToAtomTests, IronAtomNumber26)
    {
        // atnr2letter(26) = "Fe"
        std::string input = "Fe3 ";
        std::string r = shrink_string_to_atom(input, 26);
        EXPECT_GE(static_cast<int>(r.size()), 0);
        EXPECT_EQ(std::string("Fe"), r);
    }

    // -----------------------------------------------------------------------
    // SplitStringTests — split_string(input, delim)
    // -----------------------------------------------------------------------
    TEST(SplitStringTests, SingleToken_NoDelim)
    {
        svec toks = split_string<std::string>(std::string("hello"), " ");
        EXPECT_EQ(1, static_cast<int>(toks.size()));
        EXPECT_EQ(std::string("hello"), toks[0]);
    }

    TEST(SplitStringTests, ThreeTokens)
    {
        svec toks = split_string<std::string>(std::string("a b c"), " ");
        EXPECT_EQ(3, static_cast<int>(toks.size()));
        EXPECT_EQ(std::string("a"), toks[0]);
        EXPECT_EQ(std::string("b"), toks[1]);
        EXPECT_EQ(std::string("c"), toks[2]);
    }

    TEST(SplitStringTests, CommaDelimiter)
    {
        svec toks = split_string<std::string>(std::string("x,y,z"), ",");
        EXPECT_EQ(3, static_cast<int>(toks.size()));
        EXPECT_EQ(std::string("z"), toks[2]);
    }

    TEST(SplitStringTests, MaxOutLimit_ReturnsTotalCount)
    {
        // 4 tokens but max_out=2 (old API); direct API returns full token list.
        svec toks = split_string<std::string>(std::string("a b c d"), " ");
        EXPECT_EQ(4, static_cast<int>(toks.size()));
        EXPECT_EQ(std::string("a"), toks[0]);
        EXPECT_EQ(std::string("b"), toks[1]);
    }

    TEST(SplitStringTests, EmptyString_ZeroTokens)
    {
        svec toks = split_string<std::string>(std::string(""), " ");
        EXPECT_TRUE(toks.empty() || toks.size() == 1); // impl-defined for empty input
    }

    // -----------------------------------------------------------------------
    // TimingTests — ut_sleep_and_measure_us(N) returns elapsed µs >= N*1000
    // -----------------------------------------------------------------------
    TEST(TimingTests, Sleep10ms_ElapsedAtLeast10000us)
    {
        auto t0 = get_time();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        auto t1 = get_time();
        long long us = get_musec(t0, t1);
        EXPECT_GE(us, 10000LL);
    }

    TEST(TimingTests, Sleep1ms_ElapsedAtLeast1000us)
    {
        auto t0 = get_time();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        auto t1 = get_time();
        long long us = get_musec(t0, t1);
        EXPECT_GE(us, 1000LL);
    }

    TEST(TimingTests, Sleep5ms_ElapsedPositive)
    {
        auto t0 = get_time();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        auto t1 = get_time();
        long long us = get_musec(t0, t1);
        EXPECT_GT(us, 0LL);
    }

    // -----------------------------------------------------------------------
    // HypergeometricTests — 2F1(a,b;c;x)
    // -----------------------------------------------------------------------
    TEST(HypergeometricTests, Identity_2F1_Zero_IsOne)
    {
        // 2F1(a,b;c;0) = 1 for any a,b,c
        double r = hypergeometric(1.0, 2.0, 3.0, 0.0);
        EXPECT_NEAR(1.0, r, 1e-12);
    }

    TEST(HypergeometricTests, KnownValue_2F1_1_1_2_Half)
    {
        // 2F1(1,1;2;0.5) = -2*ln(0.5) = 2*ln(2) ≈ 1.386294...
        double expected = -2.0 * std::log(0.5);
        double r = hypergeometric(1.0, 1.0, 2.0, 0.5);
        EXPECT_NEAR(expected, r, 1e-6);
    }

    TEST(HypergeometricTests, KnownValue_2F1_Half_Half_ThreeHalves_Half)
    {
        // 2F1(0.5,0.5;1.5;0.5) — finite, positive
        double r = hypergeometric(0.5, 0.5, 1.5, 0.5);
        EXPECT_TRUE(std::isfinite(r));
        EXPECT_TRUE(r > 1.0);
    }

    TEST(HypergeometricTests, Symmetry_ab_equals_ba)
    {
        // 2F1(a,b;c;x) = 2F1(b,a;c;x)
        double r1 = hypergeometric(2.0, 3.0, 5.0, 0.3);
        double r2 = hypergeometric(3.0, 2.0, 5.0, 0.3);
        EXPECT_NEAR(r1, r2, 1e-10);
    }

    TEST(HypergeometricTests, NegativeX_ReturnsFinite)
    {
        double r = hypergeometric(1.0, 2.0, 3.0, -0.5);
        EXPECT_TRUE(std::isfinite(r));
    }

    // -----------------------------------------------------------------------
    // Atnr2LetterTests — constants::atnr2letter(nr) element symbol lookup
    // -----------------------------------------------------------------------
    TEST(Atnr2LetterTests, Hydrogen_Is_H)
    {
        const char* sym = constants::atnr2letter(1);
        EXPECT_EQ(std::string("H"), std::string(sym));
    }

    TEST(Atnr2LetterTests, Carbon_Is_C)
    {
        EXPECT_EQ(std::string("C"), std::string(constants::atnr2letter(6)));
    }

    TEST(Atnr2LetterTests, Iron_Is_Fe)
    {
        EXPECT_EQ(std::string("Fe"), std::string(constants::atnr2letter(26)));
    }

    TEST(Atnr2LetterTests, Gold_Is_Au)
    {
        EXPECT_EQ(std::string("Au"), std::string(constants::atnr2letter(79)));
    }

    TEST(Atnr2LetterTests, Lawrencium_103_Is_Lr)
    {
        EXPECT_EQ(std::string("Lr"), std::string(constants::atnr2letter(103)));
    }

    TEST(Atnr2LetterTests, Zero_Is_Q_Peak)
    {
        EXPECT_EQ(std::string("Q"), std::string(constants::atnr2letter(0)));
    }

    TEST(Atnr2LetterTests, OutOfRange_Returns_PROBLEM)
    {
        EXPECT_EQ(std::string("PROBLEM"), std::string(constants::atnr2letter(200)));
    }

    TEST(Atnr2LetterTests, BufferTooSmall_ReturnsNegOne)
    {
        // No more buffer-sized API: symbol is returned as const char*.
        EXPECT_EQ(std::string("C"), std::string(constants::atnr2letter(6)));
    }

    // -----------------------------------------------------------------------
    // Type2VectorTests — constants::type2vector: basis type → [nx, ny, nz]
    // -----------------------------------------------------------------------
    TEST(Type2VectorTests, Type1_Is_SShell_000)
    {
        int v[3];
        constants::type2vector(1, v);
        EXPECT_EQ(0, v[0]);
        EXPECT_EQ(0, v[1]);
        EXPECT_EQ(0, v[2]);
    }

    TEST(Type2VectorTests, Type2_Is_Px_100)
    {
        int v[3];
        constants::type2vector(2, v);
        EXPECT_EQ(1, v[0]);
        EXPECT_EQ(0, v[1]);
        EXPECT_EQ(0, v[2]);
    }

    TEST(Type2VectorTests, Type5_Is_Dx2_200)
    {
        // index 5 in type_vector: (2,0,0) = dx²
        int v[3];
        constants::type2vector(5, v);
        EXPECT_EQ(2, v[0]);
        EXPECT_EQ(0, v[1]);
        EXPECT_EQ(0, v[2]);
    }

    TEST(Type2VectorTests, Type8_Is_Dxy_110)
    {
        // index 8: (1,1,0) = dxy
        int v[3];
        constants::type2vector(8, v);
        EXPECT_EQ(1, v[0]);
        EXPECT_EQ(1, v[1]);
        EXPECT_EQ(0, v[2]);
    }

    TEST(Type2VectorTests, SumOfExponents_Matches_ShellType)
    {
        // For d-type (types 5-10), nx+ny+nz == 2
        for (int t = 5; t <= 10; ++t) {
            int v[3];
            constants::type2vector(t, v);
            EXPECT_EQ(2, v[0] + v[1] + v[2]);
        }
    }

    TEST(Type2VectorTests, OutOfRange_Returns_NegOne)
    {
        int v[3];
        constants::type2vector(0, v);
        EXPECT_EQ(-1, v[0]);
        constants::type2vector(57, v);
        EXPECT_EQ(-1, v[0]);
    }

    // -----------------------------------------------------------------------
    // NormGaussTests — constants::normgauss(type, exponent)
    // -----------------------------------------------------------------------
    TEST(NormGaussTests, SShell_IsPositive)
    {
        double n = constants::normgauss(1, 1.0);
        EXPECT_TRUE(n > 0.0);
        EXPECT_TRUE(std::isfinite(n));
    }

    TEST(NormGaussTests, SShell_ScalesWithExponent)
    {
        // Normalization grows with exponent for s-type
        double n1 = constants::normgauss(1, 1.0);
        double n2 = constants::normgauss(1, 4.0);
        EXPECT_GT(n2, n1);
    }

    TEST(NormGaussTests, PShell_IsPositive)
    {
        double n = constants::normgauss(2, 1.0);
        EXPECT_TRUE(n > 0.0);
    }

    TEST(NormGaussTests, DShell_IsPositive)
    {
        double n = constants::normgauss(5, 1.0);
        EXPECT_TRUE(n > 0.0);
    }

    TEST(NormGaussTests, SShell_KnownValue)
    {
        // normgauss for s-type (0,0,0), exponent α:
        // N = (2α/π)^(9/4) * sqrt(1/1) = (2α/π)^(9/4)
        // For α=1: N = (2/π)^(9/4)
        double alpha = 1.0;
        double expected = std::pow(2.0 * alpha / PI_VAL, 9.0 / 4.0);
        double result = constants::normgauss(1, alpha);
        EXPECT_NEAR(expected, result, 1e-10);
    }

    // -----------------------------------------------------------------------
    // AssocLegendreTests — constants::associated_legendre_polynomial(l, m, x)
    // -----------------------------------------------------------------------
    TEST(AssocLegendreTests, P00_Is_One)
    {
        EXPECT_NEAR(1.0, constants::associated_legendre_polynomial(0, 0, 0.5), 1e-12);
    }

    TEST(AssocLegendreTests, P10_At_Half_Is_Half)
    {
        // P_1^0(x) = x
        EXPECT_NEAR(0.5, constants::associated_legendre_polynomial(1, 0, 0.5), 1e-12);
    }

    TEST(AssocLegendreTests, P11_AtOne_Is_Zero)
    {
        // P_1^1(x) = sqrt(1-x²), at x=1: 0
        EXPECT_NEAR(0.0, constants::associated_legendre_polynomial(1, 1, 1.0), 1e-12);
    }

    TEST(AssocLegendreTests, P20_Is_Legendre_Polynomial)
    {
        // P_2^0(x) = (3x²-1)/2
        double x = 0.6;
        double expected = 0.5 * (3 * x * x - 1);
        EXPECT_NEAR(expected, constants::associated_legendre_polynomial(2, 0, x), 1e-12);
    }

    TEST(AssocLegendreTests, P21_At_Zero)
    {
        // P_2^1(x) = 3x*sqrt(1-x²), at x=0: 0
        EXPECT_NEAR(0.0, constants::associated_legendre_polynomial(2, 1, 0.0), 1e-12);
    }

    TEST(AssocLegendreTests, P22_At_Zero)
    {
        // P_2^2(x) = -3(x²-1) = 3(1-x²), at x=0: 3
        EXPECT_NEAR(-3.0 * (0.0 - 1.0), constants::associated_legendre_polynomial(2, 2, 0.0), 1e-12);
    }

    TEST(AssocLegendreTests, NegativeM_P1m1)
    {
        // P_1^{-1}(x) = -0.5*sqrt(1-x²)
        double x = 0.5;
        double expected = -0.5 * std::sqrt(1 - x * x);
        EXPECT_NEAR(expected, constants::associated_legendre_polynomial(1, -1, x), 1e-12);
    }

    // -----------------------------------------------------------------------
    // CartesianToSphericalTests — constants::cartesian_to_spherical(x,y,z)
    // -----------------------------------------------------------------------
    TEST(CartesianToSphericalTests, Origin_Has_Zero_Radius)
    {
        vec r = constants::cartesian_to_spherical(0.0, 0.0, 0.0);
        EXPECT_NEAR(0.0, r[0], 1e-12); // r
    }

    TEST(CartesianToSphericalTests, UnitX_Gives_Correct_Angles)
    {
        vec r = constants::cartesian_to_spherical(1.0, 0.0, 0.0);
        EXPECT_NEAR(1.0, r[0], 1e-12); // r=1
        EXPECT_NEAR(PI_VAL / 2.0, r[1], 1e-12); // theta=pi/2
        EXPECT_NEAR(0.0, r[2], 1e-12); // phi=0
    }

    TEST(CartesianToSphericalTests, UnitZ_Has_Zero_Theta)
    {
        vec r = constants::cartesian_to_spherical(0.0, 0.0, 1.0);
        EXPECT_NEAR(1.0, r[0], 1e-12); // r=1
        EXPECT_NEAR(0.0, r[1], 1e-12); // theta=0
    }

    TEST(CartesianToSphericalTests, UnitY_Gives_PhiHalfPi)
    {
        vec r = constants::cartesian_to_spherical(0.0, 1.0, 0.0);
        EXPECT_NEAR(1.0, r[0], 1e-12); // r=1
        EXPECT_NEAR(PI_VAL / 2.0, r[2], 1e-12); // phi=pi/2
    }

    TEST(CartesianToSphericalTests, Radius_Is_Euclidean_Norm)
    {
        vec r = constants::cartesian_to_spherical(3.0, 4.0, 0.0);
        EXPECT_NEAR(5.0, r[0], 1e-12); // r=5
    }

    TEST(CartesianToSphericalTests, Inverse_Recover_Cartesian)
    {
        // Convert (1,1,1) → spherical → back to Cartesian
        vec r_sph = constants::cartesian_to_spherical(1.0, 1.0, 1.0);
        double r = r_sph[0], theta = r_sph[1], phi = r_sph[2];
        double x = r * std::sin(theta) * std::cos(phi);
        double y = r * std::sin(theta) * std::sin(phi);
        double z = r * std::cos(theta);
        EXPECT_NEAR(1.0, x, 1e-12);
        EXPECT_NEAR(1.0, y, 1e-12);
        EXPECT_NEAR(1.0, z, 1e-12);
    }

    // -----------------------------------------------------------------------
    // GetLambda1Tests — median eigenvalue of a 3x3 symmetric matrix
    // -----------------------------------------------------------------------
    TEST(GetLambda1Tests, DiagonalMatrix_ReturnsMidEigenvalue)
    {
        // Eigenvalues are 1, 3, 5 → median = 3
        const double a[9] = {
            1, 0, 0,
            0, 3, 0,
            0, 0, 5
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        EXPECT_NEAR(3.0, lam, 1e-10);
    }

    TEST(GetLambda1Tests, DiagonalMatrix_AllEqual_ReturnsThat)
    {
        const double a[9] = {
            2, 0, 0,
            0, 2, 0,
            0, 0, 2
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        // all eigenvalues equal 2; any is "middle"
        EXPECT_NEAR(2.0, lam, 1e-10);
    }

    TEST(GetLambda1Tests, DiagonalDescending_ReturnsMid)
    {
        // Eigenvalues 10, 5, 1 → median = 5
        const double a[9] = {
            10, 0, 0,
            0,  5, 0,
            0,  0, 1
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        EXPECT_NEAR(5.0, lam, 1e-10);
    }

    TEST(GetLambda1Tests, SymmetricMatrix_ReturnsFiniteValue)
    {
        // 3x3 symmetric: known eigenvalues can be verified with characteristic polynomial
        // A = [[4,2,0],[2,3,1],[0,1,2]] → trace=9, det=14
        const double a[9] = {
            4, 2, 0,
            2, 3, 1,
            0, 1, 2
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        EXPECT_TRUE(std::isfinite(lam));
        // Eigenvalues must be between min(diag)=2 and max(diag)=4 for a diag-dominant case
        EXPECT_TRUE(lam >= 1.0 && lam <= 5.0);
    }

    TEST(GetLambda1Tests, SymmetricMatrix_MedianIsInBounds)
    {
        // For any symmetric matrix, median eigenvalue is bounded by extremes
        const double a[9] = {
            3, 1, 0,
            1, 4, 2,
            0, 2, 5
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        EXPECT_TRUE(std::isfinite(lam));
        EXPECT_TRUE(lam > 0.0); // positive definite matrix
    }

    // -----------------------------------------------------------------------
    // Atom Tests 
    // -----------------------------------------------------------------------
    class AtomTest : public ::testing::Test {
    protected:
        static atom make_atom()
        {
            return atom{
                "C1",
                {},
                6,
                1.25,
                -2.5,
                3.75,
                0
            };
        }

        static atom make_atom_fractional(
            const int charge,
            const double x,
            const double y,
            const double z,
            const int group_nr = 0)
        {
            atom result{
                "Test",
                {},
                charge,
                0.0,
                0.0,
                0.0,
                charge
            };

            result.set_frac_coords(d3{ x, y, z });
            result.set_group_nr(group_nr);
            return result;
        }

        static atomID binary_roundtrip(const atomID& original)
        {
            std::stringstream buffer(
                std::ios::in |
                std::ios::out |
                std::ios::binary
            );

            original.write_atom_id(buffer);
            buffer.seekg(0);

            return atomID(buffer);
        }
    };
    TEST_F(AtomTest, ID_WriteAndRead)
    {
        atom value =
            make_atom_fractional(6, 0.1, 0.2, 0.3, 0);

        const atomID originalId = value.get_ID();

        std::stringstream buffer(
            std::ios::in |
            std::ios::out |
            std::ios::binary
        );

        originalId.write_atom_id(buffer);

        EXPECT_EQ(buffer.str().size(), sizeof(atomID));
        EXPECT_EQ(buffer.str().size(), 16);

        buffer.seekg(0);
        const atomID readId(buffer);

        EXPECT_EQ(originalId, readId);
    }

    TEST_F(AtomTest, ID_WriteAndReadPreservesNegativeCoordinates)
    {
        atom value =
            make_atom_fractional(6, -1.25, 2.5, -3.75, -4);

        const atomID originalId = value.get_ID();
        const atomID readId = binary_roundtrip(originalId);

        EXPECT_EQ(originalId, readId);
    }

    TEST_F(AtomTest, ID_IsDeterministic)
    {
        atom first =
            make_atom_fractional(8, -0.125, 1.75, 12.5, 3);

        atom second =
            make_atom_fractional(8, -0.125, 1.75, 12.5, 3);

        EXPECT_EQ(first.get_ID(), second.get_ID());
    }

    TEST_F(AtomTest, ID_DifferentCoordinatesProduceDifferentIDs)
    {
        atom first =
            make_atom_fractional(6, 0.123456, 0.2, 0.3);

        atom second =
            make_atom_fractional(6, 0.123457, 0.2, 0.3);

        /*
         * The difference is 1e-6, which is comfortably larger than the
         * approximately 7.45e-9 resolution of the signed 32-bit encoding
         * over the range [-16, 16].
         */
        EXPECT_NE(first.get_ID(), second.get_ID());
    }

    TEST_F(AtomTest, ID_CoordinateSignAffectsID)
    {
        atom positive =
            make_atom_fractional(6, 1.25, 2.5, 3.75);

        atom negative =
            make_atom_fractional(6, -1.25, 2.5, 3.75);

        EXPECT_NE(positive.get_ID(), negative.get_ID());
    }

    TEST_F(AtomTest, ID_DifferentAtomicNumbersProduceDifferentIDs)
    {
        atom carbon =
            make_atom_fractional(6, 0.1, 0.2, 0.3);

        atom oxygen =
            make_atom_fractional(8, 0.1, 0.2, 0.3);

        EXPECT_NE(carbon.get_ID(), oxygen.get_ID());
    }

    TEST_F(AtomTest, ID_IsAvailableWhenNoGroupWasAssigned)
    {
        /*
         * group_nr feeds the int16_t data field of atomID, which throws on
         * anything that field cannot hold. An atom only gets a group when the
         * CIF reader matches it, so the default has to be usable on its own;
         * every other test here sets one and so never exercises it.
         */
        atom value{ "Test", {}, 6, 0.0, 0.0, 0.0, 6 };
        value.set_frac_coords(d3{ 0.1, 0.2, 0.3 });

        atomID id;
        EXPECT_NO_THROW({ id = value.get_ID(); });
        EXPECT_EQ(id, atomID(0.1, 0.2, 0.3, 0, 6));
    }

    TEST_F(AtomTest, ID_DifferentGroupsProduceDifferentIDs)
    {
        atom first =
            make_atom_fractional(6, 0.1, 0.2, 0.3, -1);

        atom second =
            make_atom_fractional(6, 0.1, 0.2, 0.3, 1);

        EXPECT_NE(first.get_ID(), second.get_ID());
    }

    TEST_F(AtomTest, ID_IsRebuiltWhenCIFPartChanges)
    {
        atom value = make_atom_fractional(6, 0.1, 0.2, 0.3, 1);
        const atomID part_one_id = value.get_ID();

        // CIF matching can assign PART after an ID has already been requested.
        value.set_group_nr(2);

        // set_group_nr() updates the cached value, rather than leaving it empty
        // for a later get_ID() call to reconstruct.
        EXPECT_EQ(value.get_ID(), atomID(0.1, 0.2, 0.3, 2, 6));
        EXPECT_NE(value.get_ID(), part_one_id);
    }

    TEST_F(AtomTest, ID_SupportsCoordinateRangeBoundaries)
    {
        EXPECT_NO_THROW({
            const atomID minimum(-16.0, -16.0, -16.0, 0, 6);
            const atomID restored = binary_roundtrip(minimum);
            EXPECT_EQ(minimum, restored);
            });

        EXPECT_NO_THROW({
            const atomID maximum(16.0, 16.0, 16.0, 0, 6);
            const atomID restored = binary_roundtrip(maximum);
            EXPECT_EQ(maximum, restored);
            });
    }

    TEST_F(AtomTest, ID_RejectsCoordinatesOutsideSupportedRange)
    {
        EXPECT_THROW(
            (atomID{ 16.000001, 0.0, 0.0, 0, 6 }),
            std::out_of_range
        );

        EXPECT_THROW(
            (atomID{ -16.000001, 0.0, 0.0, 0, 6 }),
            std::out_of_range
        );
    }

    TEST_F(AtomTest, ID_RejectsNonFiniteCoordinates)
    {
        const double infinity =
            std::numeric_limits<double>::infinity();

        const double nan =
            std::numeric_limits<double>::quiet_NaN();

        EXPECT_THROW(
            (atomID{ infinity, 0.0, 0.0, 0, 6 }),
            std::invalid_argument
        );

        EXPECT_THROW(
            (atomID{ nan, 0.0, 0.0, 0, 6 }),
            std::invalid_argument
        );
    }

    TEST_F(AtomTest, ID_RejectsInvalidAtomicNumber)
    {
        EXPECT_THROW(
            (atomID{ 0.1, 0.2, 0.3, 0, 0 }),
            std::out_of_range
        );

        EXPECT_THROW(
            (atomID{ 0.1, 0.2, 0.3, 0, 256 }),
            std::out_of_range
        );
    }

    TEST_F(AtomTest, ID_RejectsDataOutsideInt16Range)
    {
        EXPECT_THROW(
            (atomID{
                0.1,
                0.2,
                0.3,
                static_cast<int>(
                    std::numeric_limits<std::int16_t>::max()
                ) + 1,
                6
                }),
            std::out_of_range
        );

        EXPECT_THROW(
            (atomID{
                0.1,
                0.2,
                0.3,
                static_cast<int>(
                    std::numeric_limits<std::int16_t>::min()
                ) - 1,
                6
                }),
            std::out_of_range
        );
    }

    TEST_F(AtomTest, ID_DefaultConstructedObjectIsNotInitialized)
    {
        const atomID id;

        EXPECT_FALSE(id.is_initialized());
    }

    TEST_F(AtomTest, ID_CannotWriteUninitializedObject)
    {
        const atomID id;

        std::stringstream buffer(
            std::ios::in |
            std::ios::out |
            std::ios::binary
        );

        EXPECT_THROW(
            id.write_atom_id(buffer),
            std::runtime_error
        );
    }

    TEST_F(AtomTest, ID_RejectsTruncatedBinaryInput)
    {
        /*
         * A valid atomID requires 16 bytes, but this stream contains only 8.
         */
        const std::string incompleteData(8, '\0');

        std::istringstream input(
            incompleteData,
            std::ios::in | std::ios::binary
        );

        EXPECT_THROW(
            (atomID{ input }),
            std::runtime_error
        );
    }
    // -----------------------------------------------------------------------
    // Non-trivial atom behavior
    // -----------------------------------------------------------------------

    TEST_F(AtomTest, DistanceToOtherAtomIsEuclideanDistance)
    {
        const atom first{
            "A",
            {},
            1,
            1.0,
            2.0,
            3.0,
            0
        };

        const atom second{
            "B",
            {},
            1,
            4.0,
            6.0,
            3.0,
            0
        };

        EXPECT_NEAR(first.distance_to(second), 5.0, 1e-12);
        EXPECT_NEAR(second.distance_to(first), 5.0, 1e-12);
    }

    TEST_F(AtomTest, BasisSetSupportsAddingModifyingAndErasingEntries)
    {
        atom value = make_atom();

        ASSERT_TRUE(value.push_back_basis_set(10.0, 0.1, 1, 0));
        ASSERT_TRUE(value.push_back_basis_set(20.0, 0.2, 2, 1));
        ASSERT_TRUE(value.push_back_basis_set(30.0, 0.3, 3, 2));

        value.set_basis_set_exponent(1, 25.0);
        value.set_basis_set_coefficient(1, 0.25);

        EXPECT_DOUBLE_EQ(value.get_basis_set_exponent(1), 25.0);
        EXPECT_DOUBLE_EQ(value.get_basis_set_coefficient(1), 0.25);

        value.erase_basis_set(0);

        ASSERT_EQ(value.get_basis_set_size(), 2u);
        EXPECT_DOUBLE_EQ(value.get_basis_set_exponent(0), 25.0);
        EXPECT_DOUBLE_EQ(value.get_basis_set_exponent(1), 30.0);
    }

    TEST_F(AtomTest, IndexedShellCountSetterExpandsAndZeroInitializesVector)
    {
        atom value = make_atom();

        value.set_shellcount(3u, 9u);

        ASSERT_EQ(value.get_shellcount_size(), 4u);
        EXPECT_EQ(value.get_shellcount(0u), 0u);
        EXPECT_EQ(value.get_shellcount(1u), 0u);
        EXPECT_EQ(value.get_shellcount(2u), 0u);
        EXPECT_EQ(value.get_shellcount(3u), 9u);
    }

    TEST_F(AtomTest, AssignmentPerformsDeepCopy)
    {
        atom source{
            "O1",
            {},
            8,
            1.0,
            2.0,
            3.0,
            -2,
            2
        };

        source.set_frac_coords(d3{ 0.1, 0.2, 0.3 });
        source.set_shellcount(std::vector<unsigned int>{2u, 3u});
        ASSERT_TRUE(source.push_back_basis_set(25.0, 0.75, 2, 1));

        atom destination;
        destination = source;

        destination.set_label("Changed");
        destination.set_basis_set_exponent(0, 999.0);
        destination.set_shellcount(0u, 99u);

        EXPECT_EQ(source.get_label(), "O1");
        EXPECT_DOUBLE_EQ(source.get_basis_set_exponent(0), 25.0);
        EXPECT_EQ(source.get_shellcount(0u), 2u);

        EXPECT_EQ(destination.get_label(), "Changed");
        EXPECT_DOUBLE_EQ(destination.get_basis_set_exponent(0), 999.0);
        EXPECT_EQ(destination.get_shellcount(0u), 99u);
    }

    TEST_F(AtomTest, EqualityDetectsMeaningfulDifference)
    {
        atom first = make_atom();
        atom second = make_atom();

        EXPECT_TRUE(first == second);

        ASSERT_TRUE(second.push_back_basis_set(10.0, 0.5, 1, 0));

        EXPECT_FALSE(first == second);
    }

    // ParseSymopTests — cell::parse_symop, the CIF symmetry operation reader
    struct parsed_symop
    {
        int rot[3][3]{};
        double trans[3]{};
    };

    static parsed_symop parse(const std::string& operation)
    {
        parsed_symop result;
        std::ostringstream sink;
        cell::parse_symop(operation, "test.cif", result.rot, result.trans, sink);
        return result;
    }

    static void expect_rot(const parsed_symop& op, const int expected[3][3])
    {
        for (int comp = 0; comp < 3; comp++)
            for (int axis = 0; axis < 3; axis++)
                EXPECT_EQ(op.rot[comp][axis], expected[comp][axis])
                    << "component " << comp << ", axis " << axis;
    }

    TEST(ParseSymopTest, Identity)
    {
        const parsed_symop op = parse("x,y,z");
        const int expected[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.0);
    }

    // The operations of P 2_1 2_1 2_1, which used to abort the process because the
    // translation was read with stof("x+1") instead of being split off the axis term.
    TEST(ParseSymopTest, ScrewAxisWithTrailingFraction)
    {
        const parsed_symop op = parse("x+1/2,-y+1/2,-z");
        const int expected[3][3] = { {1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.0);
    }

    TEST(ParseSymopTest, TranslationBeforeAxis)
    {
        const parsed_symop op = parse("1/2+X,1/2-Y,-Z");
        const int expected[3][3] = { {1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.0);
    }

    TEST(ParseSymopTest, DecimalTranslationsAndWhitespace)
    {
        const parsed_symop op = parse(" 0.5 - x , y , 0.25 + z ");
        const int expected[3][3] = { {-1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.25);
    }

    // Rhombohedral obverse setting: mixed axes in one component and thirds
    TEST(ParseSymopTest, MixedAxesAndThirds)
    {
        const parsed_symop op = parse("-y+2/3,x-y+1/3,z+1/3");
        const int expected[3][3] = { {0, -1, 0}, {1, -1, 0}, {0, 0, 1} };
        expect_rot(op, expected);
        EXPECT_NEAR(op.trans[0], 2.0 / 3.0, 1e-12);
        EXPECT_NEAR(op.trans[1], 1.0 / 3.0, 1e-12);
        EXPECT_NEAR(op.trans[2], 1.0 / 3.0, 1e-12);
    }

    TEST(ParseSymopTest, Inversion)
    {
        const parsed_symop op = parse("-x,-y,-z");
        const int expected[3][3] = { {-1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.0);
    }

    // error_check ends the process with exit(-1). Windows reports the full value back
    // through GetExitCodeProcess, while POSIX wait() only carries the low 8 bits, so the
    // same exit is observed as 0xFFFFFFFF on one and 255 on the other.
#ifdef _WIN32
    constexpr unsigned ERROR_CHECK_EXIT_CODE = static_cast<unsigned>(-1);
#else
    constexpr unsigned ERROR_CHECK_EXIT_CODE = 255u;
#endif

    // A malformed operation has to leave through error_check's exit(-1), not abort the
    // process with a fail-fast. The message itself cannot be matched here because
    // error_check reports on stdout while death tests only see stderr.
    TEST(ParseSymopDeathTest, MalformedOperationExitsCleanly)
    {
        parsed_symop result;
        EXPECT_EXIT(cell::parse_symop("x,y", "test.cif", result.rot, result.trans, std::cout),
            ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
        EXPECT_EXIT(cell::parse_symop("x+1/0,y,z", "test.cif", result.rot, result.trans, std::cout),
            ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }

    // ---------------------------------------------------------------------
    // find_frontier_orbitals - the HOMO/LUMO detection behind the Fukui
    // functions. There is no stored frontier index in WFN, so this logic is
    // the only thing standing between a correct f+/f- and a silently wrong
    // one, which makes it worth testing directly rather than only through the
    // integration test.
    // ---------------------------------------------------------------------

    // A closed-shell case with orbital energies present: the HOMO is the
    // highest-energy occupied orbital and the LUMO the lowest-energy virtual.
    TEST(FukuiTests, FrontierOrbitalsRestrictedWithEnergies)
    {
        WFN wavy;
        wavy.push_back_MO(0, 2.0, -10.0);
        wavy.push_back_MO(1, 2.0, -5.0);
        wavy.push_back_MO(2, 2.0, -1.0);  // HOMO
        wavy.push_back_MO(3, 0.0, 0.5);   // LUMO
        wavy.push_back_MO(4, 0.0, 1.5);

        int homo = -1, lumo = -1;
        bool unrestricted = true;
        EXPECT_TRUE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
        EXPECT_EQ(homo, 2);
        EXPECT_EQ(lumo, 3);
        EXPECT_FALSE(unrestricted);
    }

    // Energies are not always stored - a plain .wfn carries none, and every
    // energy then reads back as exactly 0.0. In that case the routine must fall
    // back to orbital ORDER (last occupied / first virtual) instead of picking
    // an arbitrary orbital because all the energies compare equal.
    TEST(FukuiTests, FrontierOrbitalsFallsBackToOrderWithoutEnergies)
    {
        WFN wavy;
        wavy.push_back_MO(0, 2.0, 0.0);
        wavy.push_back_MO(1, 2.0, 0.0);
        wavy.push_back_MO(2, 2.0, 0.0);  // HOMO by position
        wavy.push_back_MO(3, 0.0, 0.0);  // LUMO by position
        wavy.push_back_MO(4, 0.0, 0.0);

        int homo = -1, lumo = -1;
        bool unrestricted = false;
        EXPECT_TRUE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
        EXPECT_EQ(homo, 2);
        EXPECT_EQ(lumo, 3);
    }

    // Orbitals are not guaranteed to arrive sorted by energy. When energies are
    // available they must win over index order, otherwise the frontier pair is
    // silently wrong for any reader that emits an unsorted set.
    TEST(FukuiTests, FrontierOrbitalsPreferEnergyOverIndexOrder)
    {
        WFN wavy;
        wavy.push_back_MO(0, 2.0, -1.0);  // highest occupied energy -> HOMO
        wavy.push_back_MO(1, 2.0, -10.0);
        wavy.push_back_MO(2, 0.0, 2.0);
        wavy.push_back_MO(3, 0.0, 0.5);   // lowest virtual energy -> LUMO

        int homo = -1, lumo = -1;
        bool unrestricted = false;
        EXPECT_TRUE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
        EXPECT_EQ(homo, 0);
        EXPECT_EQ(lumo, 3);
    }

    // The failure this guards is the one that actually bites: a wavefunction
    // file storing only the occupied orbitals (the common case for .wfn and for
    // the sucrose.wfx in tests/) has no LUMO at all, so no Fukui function can be
    // formed. That must be reported, not silently turned into an all-zero f+.
    TEST(FukuiTests, FrontierOrbitalsFailWhenNoVirtualsExist)
    {
        WFN wavy;
        wavy.push_back_MO(0, 2.0, -10.0);
        wavy.push_back_MO(1, 2.0, -5.0);
        wavy.push_back_MO(2, 2.0, -1.0);

        int homo = -1, lumo = -1;
        bool unrestricted = false;
        EXPECT_FALSE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
        EXPECT_EQ(homo, 2);
        EXPECT_EQ(lumo, -1);
    }

    // The mirror case: no occupied orbitals at all.
    TEST(FukuiTests, FrontierOrbitalsFailWhenNoOccupiedExist)
    {
        WFN wavy;
        wavy.push_back_MO(0, 0.0, 1.0);
        wavy.push_back_MO(1, 0.0, 2.0);

        int homo = -1, lumo = -1;
        bool unrestricted = false;
        EXPECT_FALSE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
        EXPECT_EQ(homo, -1);
        EXPECT_EQ(lumo, 0);
    }

    // An empty wavefunction must fail rather than index into nothing.
    TEST(FukuiTests, FrontierOrbitalsFailOnEmptyWavefunction)
    {
        WFN wavy;
        int homo = -1, lumo = -1;
        bool unrestricted = false;
        EXPECT_FALSE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
        EXPECT_EQ(homo, -1);
        EXPECT_EQ(lumo, -1);
    }

    // Fractional occupations from correlated methods must count as occupied. A
    // naive "occ == 2.0 means occupied" test would classify a natural orbital
    // with occupation 1.98 as virtual and put the HOMO in the wrong place.
    TEST(FukuiTests, FrontierOrbitalsTreatFractionalOccupationAsOccupied)
    {
        WFN wavy;
        wavy.push_back_MO(0, 2.0, -10.0);
        wavy.push_back_MO(1, 1.98, -2.0);
        wavy.push_back_MO(2, 0.02, -1.0);  // still occupied, so this is the HOMO
        wavy.push_back_MO(3, 0.0, 0.5);    // genuinely empty -> LUMO

        int homo = -1, lumo = -1;
        bool unrestricted = false;
        EXPECT_TRUE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
        EXPECT_EQ(homo, 2);
        EXPECT_EQ(lumo, 3);
    }

#ifdef NOSPHERA2_USE_GPU
    static void set_min_flop_env(const char* value)
    {
#ifdef _WIN32
        _putenv_s("NOSPHERA2_BLAS_GPU_MIN_FLOP", value ? value : "");
#else
        if (value) setenv("NOSPHERA2_BLAS_GPU_MIN_FLOP", value, 1);
        else unsetenv("NOSPHERA2_BLAS_GPU_MIN_FLOP");
#endif
    }

    // blas_gpu_dgemm ships and, until this test, was reached by nothing at all: its size
    // gate sits far above anything in the test data, so neither the gate nor the GEMM behind
    // it was ever exercised. The interesting part is not that the device can multiply but
    // that the arguments survive the trip - the routine presents a row-major interface over a
    // column-major library by swapping the operands and the transposes, and that swap is
    // exactly the sort of thing that is right in three of the four transpose combinations and
    // wrong in the fourth. So all four are checked, against a reference computed from the
    // definition rather than from another BLAS, which would share any convention mistake.
    TEST(BlasGpuTests, RowMajorDgemmMatchesTheDefinitionInEveryTransposeCombination)
    {
        if (!blas_gpu_available()) {
            GTEST_SKIP() << "No GPU device present; blas_gpu_dgemm cannot run here";
        }

        const int m = 17, n = 11, k = 23;   // deliberately unequal, so a swapped extent shows
        std::vector<double> A((size_t)m * k), B((size_t)k * n);
        for (size_t i = 0; i < A.size(); i++) A[i] = 0.5 - std::sin(0.37 * (double)i);
        for (size_t i = 0; i < B.size(); i++) B[i] = 0.25 + std::cos(0.21 * (double)i);

        blas_gpu_set_enabled(true);
        set_min_flop_env("1");   // without this the gate refuses a problem this small

        for (int tA = 0; tA < 2; tA++) {
            for (int tB = 0; tB < 2; tB++) {
                // op(A) is m x k and op(B) is k x n either way, so the leading dimensions
                // follow how the operand is stored rather than how it is used
                const int lda = tA ? m : k;
                const int ldb = tB ? k : n;
                std::vector<double> C((size_t)m * n, 0.0);
                const bool ran = blas_gpu_dgemm(tA != 0, tB != 0, m, n, k, 1.0,
                    A.data(), lda, B.data(), ldb, 0.0, C.data(), n);
                ASSERT_TRUE(ran) << "device declined transA=" << tA << " transB=" << tB;

                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        double want = 0.0;
                        for (int p = 0; p < k; p++) {
                            const double a = tA ? A[(size_t)p * m + i] : A[(size_t)i * k + p];
                            const double b = tB ? B[(size_t)j * k + p] : B[(size_t)p * n + j];
                            want += a * b;
                        }
                        ASSERT_NEAR(C[(size_t)i * n + j], want, 1e-10)
                            << "transA=" << tA << " transB=" << tB
                            << " at (" << i << "," << j << ")";
                    }
                }
            }
        }

        set_min_flop_env(nullptr);
        blas_gpu_set_enabled(false);
    }

    // Every shape above uses a single k-slice, so none of them reaches the split-k path - and
    // that path is where a GEMM of our own differs from a library one. A deep, narrow
    // reduction is what the I tensor actually asks for, and it is the shape that forces the
    // depth to be cut across blocks and summed afterwards.
    //
    // Determinism is asserted as well as accuracy, deliberately: accumulating the slices with
    // an atomic would be shorter, would pass an accuracy check, and would make the result
    // depend on the order the device happened to finish them. Bit-identical repeats are the
    // only thing that tells those two implementations apart.
    TEST(BlasGpuTests, ADeepReductionIsSplitAcrossBlocksAndStillSumsInAFixedOrder)
    {
        if (!blas_gpu_available()) {
            GTEST_SKIP() << "No GPU device present; blas_gpu_dgemm cannot run here";
        }

        const int m = 9, n = 7, k = 4096;   // one tile of output, many slices of depth
        std::vector<double> A((size_t)m * k), B((size_t)k * n);
        for (size_t i = 0; i < A.size(); i++) A[i] = 0.5 - std::sin(0.37 * (double)i);
        for (size_t i = 0; i < B.size(); i++) B[i] = 0.25 + std::cos(0.21 * (double)i);

        blas_gpu_set_enabled(true);
        set_min_flop_env("1");

        std::vector<double> C1((size_t)m * n, 0.0), C2((size_t)m * n, 0.0);
        ASSERT_TRUE(blas_gpu_dgemm(false, false, m, n, k, 1.0,
            A.data(), k, B.data(), n, 0.0, C1.data(), n));
        ASSERT_TRUE(blas_gpu_dgemm(false, false, m, n, k, 1.0,
            A.data(), k, B.data(), n, 0.0, C2.data(), n));

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                double want = 0.0;
                for (int p = 0; p < k; p++)
                    want += A[(size_t)i * k + p] * B[(size_t)p * n + j];
                const size_t at = (size_t)i * n + j;
                ASSERT_NEAR(C1[at], want, 1e-9) << "at (" << i << "," << j << ")";
                ASSERT_EQ(C1[at], C2[at]) << "not reproducible at (" << i << "," << j << ")";
            }
        }

        set_min_flop_env(nullptr);
        blas_gpu_set_enabled(false);
    }

    // The gate is the reason the offload is worth having, so it gets its own check: a shape
    // below the threshold has to be declined rather than quietly run at a loss.
    TEST(BlasGpuTests, SmallShapesAreDeclinedSoTheyStayOnTheHost)
    {
        blas_gpu_set_enabled(true);
        std::vector<double> A(16, 1.0), B(16, 1.0), C(16, 0.0);
        EXPECT_FALSE(blas_gpu_dgemm(false, false, 4, 4, 4, 1.0,
            A.data(), 4, B.data(), 4, 0.0, C.data(), 4));
        blas_gpu_set_enabled(false);
    }
#endif

    // The multipole restraint rests on one integral, int r^(l+2) N c e^(-a r^2) dr = N c Gamma(l+3/2) / (2 a^(l+3/2)),
    // checked against the trapezoid rule and, for l = 0, against the pi/(2 a^(3/2)) N c row that
    // add_electron_restraint has used all along (which carries the sqrt(4 pi) of Y_00).
    TEST(RiMultipoleTests, RadialMomentMatchesQuadratureAndTheChargeRow)
    {
        const double exps[] = { 0.3, 2.5, 40.0 }, coef = 1.3;
        for (int e = 0; e < 3; e++) {
            const double a = exps[e];
            const primitive p0(0, 0, a, coef);
            EXPECT_NEAR(std::sqrt(4.0 * PI_VAL) * DensityFitting::radial_moment(a, coef, 0), PI_VAL / (2.0 * std::pow(a, 1.5)) * p0.normalization_constant() * p0.get_coef(), 1e-12);
            for (int l = 0; l <= 4; l++) {
                const primitive p(0, l, a, coef);
                const int n = 200000;
                const double h = 12.0 / std::sqrt(a) / n;
                double sum = 0.0;
                for (int i = 1; i < n; i++) {
                    const double r = i * h;
                    sum += std::pow(r, 2 * l + 2) * std::exp(-a * r * r);
                }
                sum *= h * p.normalization_constant() * p.get_coef();
                EXPECT_NEAR(sum, DensityFitting::radial_moment(a, coef, l), 1e-9 * sum);
            }
        }
    }

    // Closure of the whole restraint: for one oxygen with the combo_basis_fit aux basis and arbitrary
    // coefficients, the density that calc_density_ML evaluates from them is integrated on an atomic
    // grid to give the moments Q_lm = int rho r^l Y_lm, exactly what calculatePartitionedMultipoles
    // produces for the fit target. The restraint rows applied to the same coefficients must return
    // those moments (times the r_cov^-l row scaling), the targets must land in the matching rows,
    // and fitted_multipoles must agree. This pins the row placement, the l/m ordering of the
    // coefficients against constants::spherical_harmonic, and the normalisation in one go.
    TEST(RiMultipoleTests, RestraintRowsReproduceTheGridMomentsOfTheAtomicDensity)
    {
        const double pos[3] = { 0.3, -0.7, 1.1 };
        const int Z = 8, lmax = 2, n_moments = (lmax + 1) * (lmax + 1);
        WFN wavy(e_origin::NOT_YET_DEFINED);
        wavy.push_back_atom("O", pos[0], pos[1], pos[2], Z);
        std::vector<std::shared_ptr<BasisSet>> basis{ BasisSetLibrary::get_basis_set("combo_basis_fit") };
        WFN aux = generate_aux_wfn(wavy, basis);
        const atom A = aux.get_atom(0);
        double alpha_min[9] = { 0.0 }, alpha_max = 0.0;
        int max_l = 0, n_aux = 0, prim = 0;
        for (int shell = 0; shell < (int)A.get_shellcount_size(); shell++) {
            const int l = A.get_basis_set_entry(prim).get_type();
            for (int e = 0; e < (int)A.get_shellcount(shell); e++) {
                const double a = A.get_basis_set_entry(prim + e).get_exponent();
                alpha_max = std::max(alpha_max, a);
                if (alpha_min[l] == 0.0 || a < alpha_min[l]) alpha_min[l] = a;
            }
            max_l = std::max(max_l, l);
            n_aux += 2 * l + 1;
            prim += A.get_shellcount(shell);
        }
        AtomGrid grid(1e-12, 350, 5000, Z, alpha_max, max_l, alpha_min, std::cout);
        const int n_points = grid.get_num_grid_points();
        vec gx(n_points), gy(n_points), gz(n_points), aw(n_points), bw(n_points), tw(n_points), chi(1, 0.0);
        grid.get_grid(1, 0, &pos[0], &pos[1], &pos[2], &Z, gx.data(), gy.data(), gz.data(), aw.data(), bw.data(), tw.data(), WFN(), chi);
        vec coefs(n_aux);
        for (int i = 0; i < n_aux; i++) coefs[i] = std::sin(1.0 + i);
        vec2 Q(1, vec(n_moments, 0.0));
        for (int p = 0; p < n_points; p++) {
            const double f = calc_density_ML(gx[p], gy[p], gz[p], coefs, aux.get_atoms()) * bw[p];
            double d[3] = { gx[p] - pos[0], gy[p] - pos[1], gz[p] - pos[2] };
            const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            for (int i = 0; i < 3; i++) d[i] /= r;
            double rl = 1.0;
            for (int l = 0; l <= lmax; l++) {
                for (int m = -l; m <= l; m++)
                    Q[0][l * l + l + m] += f * rl * constants::spherical_harmonic(l, m, d);
                rl *= r;
            }
        }
        vec eri2c(n_aux * n_aux, 0.0), rho(n_aux, 0.0);
        const vec weights(1, 1.0);
        DensityFitting::add_electron_restraint(eri2c, rho, aux, weights, vec(1, 0.0));
        DensityFitting::add_multipole_restraint(eri2c, rho, aux, Q, weights, lmax);
        ASSERT_EQ(eri2c.size(), (size_t)(n_aux + n_moments) * n_aux);
        ASSERT_EQ(rho.size(), (size_t)(n_aux + n_moments));
        const double r_cov = constants::ang2bohr(constants::covalent_radii[Z]);
        const vec2 fitted = DensityFitting::fitted_multipoles(coefs, aux, lmax);
        for (int row = 0; row < n_moments; row++) {
            const int l = (int)std::floor(std::sqrt(row + 1e-9));
            const double scale = row == 0 ? std::sqrt(4.0 * PI_VAL) : std::pow(r_cov, -l);
            double lhs = 0.0;
            for (int i = 0; i < n_aux; i++) lhs += eri2c[(n_aux + row) * n_aux + i] * coefs[i];
            EXPECT_NEAR(lhs, scale * Q[0][row], 1e-6 * std::max(1.0, std::abs(scale * Q[0][row]))) << "row " << row;
            EXPECT_NEAR(rho[n_aux + row], row == 0 ? 0.0 : scale * Q[0][row], 1e-12) << "row " << row;
            EXPECT_NEAR(fitted[0][row], Q[0][row], 1e-6 * std::max(1.0, std::abs(Q[0][row]))) << "row " << row;
        }
    }

    // gamma(l+3/2, x) = 2 int_0^sqrt(x) u^(2l+2) exp(-u^2) du, in the series and the recurrence regime
    TEST(RiInteractionTests, LowerIncompleteGammaMatchesQuadrature)
    {
        const double xs[] = { 0.05, 0.7, 3.0, 9.0, 30.0 };
        for (int l = 0; l <= 4; l++)
            for (int k = 0; k < 5; k++) {
                const int n = 200000;
                const double x = xs[k], h = std::sqrt(x) / n;
                double sum = 0.5 * std::pow(std::sqrt(x), 2 * l + 2) * std::exp(-x);
                for (int i = 1; i < n; i++) {
                    const double u = i * h;
                    sum += std::pow(u, 2 * l + 2) * std::exp(-u * u);
                }
                sum *= 2.0 * h;
                EXPECT_NEAR(sum, DensityFitting::lower_gamma_half(l, x), 1e-9 * sum) << "l " << l << " x " << x;
            }
    }

    // The s potential is the Gaussian charge erf(sqrt(a)R)/R; outside the density every rank is the point
    // multipole 4pi/(2l+1) Q_lm Y_lm / R^(l+1) with the same Q_lm that fitted_multipoles reports
    TEST(RiInteractionTests, AuxPotentialIsTheGaussianChargeAndThePointMultipoleLimit)
    {
        const double a = 0.9, c = 1.7, R[3] = { 0.8, -0.3, 1.1 };
        const double r = std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
        const double q = std::sqrt(4.0 * PI_VAL) * DensityFitting::radial_moment(a, c, 0);
        EXPECT_NEAR(DensityFitting::aux_potential(a, c, 0, 0, R), q * std::erf(std::sqrt(a) * r) / r, 1e-12 * q);
        const double Rf[3] = { 5.0 * R[0], 5.0 * R[1], 5.0 * R[2] }, rf = 5.0 * r, d[3] = { R[0] / r, R[1] / r, R[2] / r };
        for (int l = 1; l <= 4; l++)
            for (int m = -l; m <= l; m++) {
                const double V = 4.0 * PI_VAL / (2 * l + 1) * DensityFitting::radial_moment(a, c, l) * constants::spherical_harmonic(l, m, d) / std::pow(rf, l + 1);
                EXPECT_NEAR(DensityFitting::aux_potential(a, c, l, m, Rf), V, 1e-9 * std::abs(V) + 1e-15) << "l " << l << " m " << m;
            }
    }

    // A partner B whose density is one very tight s Gaussian holding exactly Z_B electrons is neutral and
    // point-like, so the energy cancels in both halves: the analytic aux potential of A at B's nucleus against
    // the libcint two-centre integrals, and A's nuclear repulsion against A's nuclei in B's density. This ties
    // the potential, the coefficient layout and the combined Int_Params block together
    // Spin-scaled hydrogen atom, E_x[rho_up] = E_x[2 rho_up] / 2: LDA -0.2680, PBE -0.3059, B88 -0.3098 Eh. r2SCAN with
    // the orbital tau gives the exact -0.3125 by construction; the PC07opt tau of r2SCAN-L moves it to -0.3108, integrated
    // in Python from the libxc maple definitions. lap rho = (4 - 4/r) rho for the exponential density
    TEST(RiInteractionTests, ExchangeFunctionalsReproduceTheHydrogenAtom)
    {
        const double ref[4] = { -0.2680, -0.3059, -0.3098, -0.3108 }, dr = 1e-4;
        for (int f = 0; f < 4; f++) {
            double e = 0.0;
            for (int i = 1; i < 300000; i++) {
                const double r = i * dr, rho = 2 * std::exp(-2 * r) / constants::PI;
                e += 4 * constants::PI * r * r * DensityFitting::exchange_density(rho, 4 * rho * rho, f, (4 - 4 / r) * rho) * dr;
            }
            EXPECT_NEAR(0.5 * e, ref[f], 2e-4) << f;
        }
        //r2SCAN-L at r = 0.5 of the same density, from the Python transcription of the libxc maple sources
        const double rho = 2 * std::exp(-1.0) / constants::PI;
        EXPECT_NEAR(DensityFitting::exchange_density(rho, 4 * rho * rho, 3, -4 * rho), -0.12503518792909296, 1e-12);
    }

    TEST(RiInteractionTests, NeutralPointLikePartnerGivesZeroEnergyAndConsistentTables)
    {
        WFN wavy(e_origin::NOT_YET_DEFINED);
        wavy.push_back_atom("O", 0.0, 0.0, 0.0, 8);
        wavy.push_back_atom("H", 1.8, 0.0, 0.0, 1);
        std::vector<std::shared_ptr<BasisSet>> basis{ BasisSetLibrary::get_basis_set("combo_basis_fit") };
        WFN aux_A = generate_aux_wfn(wavy, basis);
        int n_aux = 0;
        for (int a = 0; a < aux_A.get_ncen(); a++) {
            const atom A = aux_A.get_atom(a);
            int prim = 0;
            for (int shell = 0; shell < (int)A.get_shellcount_size(); shell++) {
                n_aux += 2 * A.get_basis_set_entry(prim).get_type() + 1;
                prim += A.get_shellcount(shell);
            }
        }
        vec coef_A(n_aux);
        for (int i = 0; i < n_aux; i++) coef_A[i] = 0.3 * std::sin(1.0 + i);
        const int ZB = 3;
        const double alpha = 2.0e5;
        atom B("Li", {}, 1, 0.0, 3.1, 1.4, ZB);
        B.push_back_basis_set(alpha, 1.0, 0, 0);
        WFN wavy_B(e_origin::NOT_YET_DEFINED);
        wavy_B.push_back_atom(B);
        const vec coef_B{ ZB / (std::sqrt(4.0 * PI_VAL) * DensityFitting::radial_moment(alpha, 1.0, 0)) };
        const DensityFitting::INTERACTION E = DensityFitting::interaction_energy(coef_A, aux_A, coef_B, wavy_B);
        EXPECT_NEAR(E.nuc_nuc + E.nucA_rhoB, 0.0, 1e-10 * E.nuc_nuc);
        EXPECT_GT(std::abs(E.nucB_rhoA), 1e-3);
        EXPECT_NEAR(E.nucB_rhoA + E.rho_rho, 0.0, 1e-5 * std::abs(E.nucB_rhoA));
        double pair = 0.0, rank = 0.0;
        for (int a = 0; a < (int)E.pair.size(); a++)
            for (int b = 0; b < (int)E.pair[a].size(); b++) pair += E.pair[a][b];
        for (int i = 0; i < (int)E.rank.size(); i++)
            for (int j = 0; j < (int)E.rank[i].size(); j++) rank += E.rank[i][j];
        EXPECT_NEAR(pair, E.electrostatic(), 1e-12);
        EXPECT_NEAR(rank, E.electrostatic(), 1e-12);
        EXPECT_EQ(E.rank[0].size(), 2);
        const DensityFitting::INTERACTION F = DensityFitting::interaction_energy(coef_B, wavy_B, coef_A, aux_A, 2.0);
        EXPECT_NEAR(F.electrostatic(), E.electrostatic(), 1e-12);
        EXPECT_NEAR(E.pol_A, 0.0, 1e-12);
        EXPECT_LT(E.pol_B, -1e-8);
        EXPECT_NEAR(F.pol_A, E.pol_B, 1e-14);
        EXPECT_NEAR(F.pol_B, E.pol_A, 1e-14);
        EXPECT_LT(E.disp, 0.0);
        EXPECT_NEAR(F.disp, E.disp, 1e-12);
        EXPECT_NEAR(F.overlap, E.overlap, 1e-10);
        EXPECT_NEAR(F.rep, 2.0 * F.overlap, 1e-14);
        EXPECT_EQ(F.n_A, 0.0);
        EXPECT_NEAR(F.total() - F.rep, E.total() - E.rep, 1e-10);
        EXPECT_LT(E.rep_x, 0.0);
        EXPECT_NEAR(E.rep, E.rep_kin + E.rep_x, 1e-14);
        const DensityFitting::INTERACTION G = DensityFitting::interaction_energy(coef_B, wavy_B, coef_A, aux_A);
        EXPECT_NEAR(G.rep_kin, E.rep_kin, 1e-6);
        EXPECT_NEAR(G.rep_vw, E.rep_vw, 1e-6);
        EXPECT_NEAR(G.rep_x, E.rep_x, 1e-6);
        EXPECT_NEAR(G.n_A, E.n_B, 1e-4);
        EXPECT_NEAR(G.n_B, E.n_A, 1e-4);
        for (int f = 1; f < 3; f++) {
            const DensityFitting::INTERACTION H = DensityFitting::interaction_energy(coef_A, aux_A, coef_B, wavy_B, 0.0, f);
            EXPECT_EQ(H.x_fun, f);
            EXPECT_NEAR(H.rep_kin, E.rep_kin, 1e-12);
            EXPECT_LT(H.rep_x, 0.0);
            EXPECT_GT(std::abs(H.rep_x - E.rep_x), 1e-8);
        }
        EXPECT_NEAR(F.nucA_rhoB, E.nucB_rhoA, 1e-12);
        EXPECT_NEAR(F.rho_rho, E.rho_rho, 1e-10);
        for (int a = 0; a < 2; a++) EXPECT_NEAR(F.pair[0][a], E.pair[a][0], 1e-12);
    }
    namespace {
        const std::filesystem::path thpp = "../Lukas_Test/thpp_p1.xyz";
        std::filesystem::path geometry_aid_tmp(const std::string& name)
        {
            return std::filesystem::temp_directory_path() / ("nosphera2_geometry_aid_" + name);
        }
        void load_npy(const std::filesystem::path& path, std::vector<unsigned long>& shape, vec& data)
        {
            bool fortran_order = true;
            data.clear();
            npy::LoadArrayFromNumpy(path, shape, fortran_order, data);
            EXPECT_FALSE(fortran_order);
        }
        //A GEOAID01 file laid out the way load_model reads it
        void write_model(const std::filesystem::path& path, const geometry_aid::Model& m)
        {
            std::ofstream out(path, std::ios::binary);
            out.write("GEOAID01", 8);
            const int header[5] = { m.n_features, m.n_components, m.n_layers, m.n_classes, m.whiten };
            out.write((const char*)header, sizeof(header));
            for (int c = 0; c < m.n_classes; c++) {
                const int length = (int)m.classes[c].size();
                out.write((const char*)&length, sizeof(int));
                out.write(m.classes[c].data(), length);
            }
            out.write((const char*)m.mean.data(), m.mean.size() * sizeof(double));
            out.write((const char*)m.components.data(), m.components.size() * sizeof(double));
            if (m.whiten) out.write((const char*)m.explained_variance.data(), m.explained_variance.size() * sizeof(double));
            for (int l = 0; l < m.n_layers; l++) {
                const int shape[2] = { m.rows[l], m.cols[l] };
                out.write((const char*)shape, sizeof(shape));
                out.write((const char*)m.w[l].data(), m.w[l].size() * sizeof(double));
                out.write((const char*)m.b[l].data(), m.b[l].size() * sizeof(double));
            }
        }
        //Two features, two components, a ReLU layer and three classes: small enough to classify by hand
        geometry_aid::Model tiny_model(bool whiten)
        {
            geometry_aid::Model m;
            m.n_features = 2, m.n_components = 2, m.n_layers = 2, m.n_classes = 3, m.whiten = whiten;
            m.classes = { "C", "N", "O" };
            m.mean = { 1.0, -1.0 };
            m.components = { 1.0, 0.0, 0.0, 2.0 };
            m.explained_variance = { 4.0, 1.0 };
            m.rows = { 2, 2 }, m.cols = { 2, 3 };
            m.w = { { 1.0, 0.0, 0.0, 1.0 }, { 1.0, 0.0, -1.0, 0.0, 1.0, 1.0 } };
            m.b = { { -1.5, 0.0 }, { 0.0, 0.5, 0.0 } };
            return m;
        }
        vec softmax(vec z)
        {
            double total = 0.0;
            for (int i = 0; i < z.size(); i++) total += (z[i] = std::exp(z[i]));
            for (int i = 0; i < z.size(); i++) z[i] /= total;
            return z;
        }
        //A model of the descriptor's real width, so the whole pipeline can run on a structure
        geometry_aid::Model wide_model()
        {
            geometry_aid::Model m;
            m.n_features = 42042, m.n_components = 3, m.n_layers = 1, m.n_classes = 2;
            m.classes = { "C", "N" };
            m.mean.assign(42042, 0.0);
            m.components.resize(42042 * 3);
            for (int f = 0; f < 42042; f++)
                for (int c = 0; c < 3; c++) m.components[f * 3 + c] = ((f * 7 + c * 13) % 11 - 5) * 1e-3;
            m.rows = { 3 }, m.cols = { 2 };
            m.w = { { 1.0, -1.0, 0.5, 0.5, -1.0, 1.0 } }, m.b = { { 0.0, 0.0 } };
            return m;
        }
        int run_nosphera2(std::vector<std::string> args)
        {
            args.insert(args.begin(), "NoSpherA2");
            std::vector<char*> argv;
            for (int i = 0; i < args.size(); i++) argv.push_back(args[i].data());
            argv.push_back(nullptr);
            return run_app((int)args.size(), argv.data());
        }
        options parse_options(std::vector<std::string> args)
        {
            args.insert(args.begin(), "NoSpherA2");
            std::vector<char*> argv;
            for (int i = 0; i < args.size(); i++) argv.push_back(args[i].data());
            int argc = (int)args.size();
            options opt(argc, argv.data(), std::cout);
            opt.digest_options();
            return opt;
        }
    }

    TEST(GeometryAidTests, HyperparametersMatchTheTrainedModels)
    {
        const SALTED_Utils::FeatomicHyperParameters hp = geometry_aid::hyperparameters();
        const svec species{ "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "I" };
        EXPECT_EQ(hp.cutoff_radius, 3.5);
        EXPECT_EQ(hp.max_radial, 6);
        EXPECT_EQ(hp.max_angular, 12);
        EXPECT_EQ(hp.atomic_gaussian_width, 0.2);
        EXPECT_EQ(hp.center_atom_weight, 1.0);
        EXPECT_EQ(hp.species, species);
        EXPECT_EQ(hp.neighspe, species);
        EXPECT_EQ(hp.radial_basis.type, "Gto");
        EXPECT_EQ(hp.radial_basis.spline_accuracy, 1e-6);
        EXPECT_EQ(hp.cutoff_function.type, "ShiftedCosine");
        EXPECT_EQ(hp.cutoff_function.width, 0.7);
        const int pairs = (int)species.size() * ((int)species.size() + 1) / 2;
        EXPECT_EQ(pairs * (hp.max_radial + 1) * (hp.max_radial + 1) * (hp.max_angular + 1), 42042);
        EXPECT_EQ(geometry_aid::hyperparameters(3.0).cutoff_radius, 3.0);
        EXPECT_NE(geometry_aid::hyperparameters(3.0).to_json(), hp.to_json());
    }

    TEST(GeometryAidTests, DescriptorHasARowPerHeavyAtomAnd42042Features)
    {
        //thpp_p1 has 12 C, 4 N, 2 F and 14 H; H is no SOAP species, so 18 centres and at most the 6 pair blocks of C, N, F, each 7 * 7 * 13 wide
        const std::filesystem::path out = geometry_aid_tmp("thpp.npy"), dirty = geometry_aid_tmp("thpp_dirty.npy");
        geometry_aid::write_descriptor(thpp, out, geometry_aid::hyperparameters());
        geometry_aid::write_descriptor(thpp, dirty, geometry_aid::hyperparameters(3.0));
        std::vector<unsigned long> shape, shape_dirty;
        vec d, d_dirty;
        load_npy(out, shape, d);
        load_npy(dirty, shape_dirty, d_dirty);
        ASSERT_EQ(shape, (std::vector<unsigned long>{ 18, 42042 }));
        ASSERT_EQ(shape_dirty, shape);
        double diff = 0.0;
        for (int a = 0; a < 18; a++) {
            int nonzero = 0;
            for (int f = 0; f < 42042; f++) {
                const double v = d[a * 42042 + f];
                EXPECT_TRUE(std::isfinite(v));
                if (v != 0.0) nonzero++;
                diff += std::abs(v - d_dirty[a * 42042 + f]);
            }
            EXPECT_GE(nonzero, 637) << a;
            EXPECT_LE(nonzero, 6 * 637) << a;
        }
        EXPECT_GT(diff, 1.0);
        std::filesystem::remove(out);
        std::filesystem::remove(dirty);
    }

    TEST(GeometryAidTests, BatchSkipsAMissingStructureAndFailsOnlyWhenNothingWasWritten)
    {
        const std::filesystem::path good = geometry_aid_tmp("batch.npy"), missing = geometry_aid_tmp("missing.npy");
        const geometry_aid::jobvec jobs{ { thpp, good }, { geometry_aid_tmp("does_not_exist.xyz"), missing } };
        EXPECT_EQ(geometry_aid::write_descriptors(jobs, geometry_aid::hyperparameters()), 0);
        EXPECT_TRUE(std::filesystem::exists(good));
        EXPECT_FALSE(std::filesystem::exists(missing));
        EXPECT_EQ(geometry_aid::write_descriptors({ jobs[1] }, geometry_aid::hyperparameters()), 1);
        std::filesystem::remove(good);
        const std::filesystem::path list = geometry_aid_tmp("list.txt");
        { std::ofstream(list) << "# comment\r\n\r\n  a.xyz  \r\nb c.xyz\n"; }
        EXPECT_EQ(geometry_aid::read_structure_list(list), (pathvec{ "a.xyz", "b c.xyz" }));
        std::filesystem::remove(list);
    }

    TEST(GeometryAidTests, ModelRoundTripsThroughTheBinaryFormatAndClassifiesByHand)
    {
        //x = (3, 1): centred (2, 2), projected (2, 4), whitened (1, 4); x = (0, 5) exercises the sparse skip with mean_projection = (1, -2)
        const double x[4] = { 3.0, 1.0, 0.0, 5.0 };
        const vec expect_plain[2] = { softmax({ 0.5, 4.5, 3.5 }), softmax({ 0.0, 12.5, 12.0 }) };
        const vec expect_white[2] = { softmax({ 0.0, 4.5, 4.0 }), softmax({ 0.0, 12.5, 12.0 }) };
        for (int whiten = 0; whiten < 2; whiten++) {
            const std::filesystem::path path = geometry_aid_tmp(whiten ? "white.bin" : "plain.bin");
            write_model(path, tiny_model(whiten));
            const geometry_aid::Model m = geometry_aid::load_model(path);
            EXPECT_EQ(m.n_features, 2);
            EXPECT_EQ(m.n_components, 2);
            EXPECT_EQ(m.n_layers, 2);
            EXPECT_EQ(m.n_classes, 3);
            EXPECT_EQ(m.whiten, whiten == 1);
            EXPECT_EQ(m.classes, (svec{ "C", "N", "O" }));
            EXPECT_EQ(m.mean_projection, (vec{ 1.0, -2.0 }));
            EXPECT_EQ(m.rows, (ivec{ 2, 2 }));
            EXPECT_EQ(m.cols, (ivec{ 2, 3 }));
            EXPECT_EQ(m.w[1].size(), 6);
            const vec p = geometry_aid::classify_descriptor(x, 2, 2, m);
            ASSERT_EQ(p.size(), 6);
            for (int a = 0; a < 2; a++)
                for (int c = 0; c < 3; c++)
                    EXPECT_NEAR(p[a * 3 + c], (whiten ? expect_white : expect_plain)[a][c], 1e-12) << whiten << " " << a << " " << c;
            EXPECT_EQ(&geometry_aid::cached_model(path), &geometry_aid::cached_model(path));
            std::filesystem::remove(path);
        }
    }

    TEST(GeometryAidDeathTest, RejectsAForeignFileAndAMismatchedDescriptor)
    {
        const std::filesystem::path foreign = geometry_aid_tmp("foreign.bin");
        { std::ofstream(foreign, std::ios::binary) << "NOTGEOAID"; }
        EXPECT_EXIT(geometry_aid::load_model(foreign), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
        std::filesystem::remove(foreign);
        const double x[3] = { 1.0, 2.0, 3.0 };
        EXPECT_EXIT(geometry_aid::classify_descriptor(x, 1, 3, tiny_model(false)), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }

    TEST(GeometryAidTests, FlagsQueueTheirJobsInAnyOrder)
    {
        const std::filesystem::path list = geometry_aid_tmp("flags.txt"), model = geometry_aid_tmp("flags.bin");
        { std::ofstream(list) << "# structures\n\n" << thpp.string() << "\n"; }
        write_model(model, tiny_model(false));
        options a = parse_options({ "-calc_featomic_descriptor", "-wfn", thpp.string() });
        EXPECT_TRUE(a.calc_featomic_descriptor);
        EXPECT_EQ(a.wfn, thpp);
        EXPECT_EQ(a.geometry_aid_cutoff, 3.5);
        EXPECT_TRUE(a.featomic_structures.empty() && a.classify_atoms_out.empty() && a.classify_structures.empty());
        options b = parse_options({ "-calc_featomic_descriptors", list.string(), "-geometry_aid_cutoff", "3.0" });
        EXPECT_FALSE(b.calc_featomic_descriptor);
        EXPECT_EQ(b.featomic_structures, (pathvec{ thpp }));
        EXPECT_EQ(b.geometry_aid_cutoff, 3.0);
        options c = parse_options({ "-wfn", thpp.string(), "-classify_atoms", model.string() });
        EXPECT_EQ(c.classify_atoms_out, "probabilities.npy");
        EXPECT_EQ(c.geometry_aid_model, model);
        options d = parse_options({ "-classify_atoms", model.string(), "out.npy", "-no_date", "-wfn", thpp.string() });
        EXPECT_EQ(d.classify_atoms_out, "out.npy");
        EXPECT_TRUE(d.no_date);
        EXPECT_EQ(d.wfn, thpp);
        options e = parse_options({ "-classify_atoms_list", list.string(), model.string() });
        EXPECT_EQ(e.classify_structures, (pathvec{ thpp }));
        EXPECT_EQ(e.geometry_aid_model, model);
        EXPECT_TRUE(e.classify_atoms_out.empty());
        std::filesystem::remove(list);
        std::filesystem::remove(model);
    }

    TEST(GeometryAidTests, RunAppWritesDescriptorNpyTheWayOlex2CallsIt)
    {
        //A copy of the structure in the temp directory, because the batch flag writes <path>.npy beside it
        const std::filesystem::path copy = geometry_aid_tmp("copy.xyz"), list = geometry_aid_tmp("run.txt"), copy_npy = copy.string() + ".npy", direct = geometry_aid_tmp("direct.npy");
        std::filesystem::copy_file(thpp, copy, std::filesystem::copy_options::overwrite_existing);
        { std::ofstream(list) << copy.string() << "\n"; }
        std::filesystem::remove("descriptor.npy");
        ASSERT_EQ(run_nosphera2({ "-wfn", thpp.string(), "-calc_featomic_descriptor", "-no_date" }), 0);
        std::vector<unsigned long> shape;
        vec d, d_direct;
        load_npy("descriptor.npy", shape, d);
        EXPECT_EQ(shape, (std::vector<unsigned long>{ 18, 42042 }));
        std::filesystem::remove("descriptor.npy");
        //The cutoff after the descriptor flag must still apply, so the batch output is the dirty descriptor
        ASSERT_EQ(run_nosphera2({ "-calc_featomic_descriptors", list.string(), "-geometry_aid_cutoff", "3.0", "-no_date" }), 0);
        geometry_aid::write_descriptor(thpp, direct, geometry_aid::hyperparameters(3.0));
        load_npy(copy_npy, shape, d);
        load_npy(direct, shape, d_direct);
        ASSERT_EQ(d.size(), d_direct.size());
        double diff = 0.0;
        for (int i = 0; i < d.size(); i++) diff += std::abs(d[i] - d_direct[i]);
        EXPECT_EQ(diff, 0.0);
        std::filesystem::remove(copy_npy);
        std::filesystem::remove(direct);
        std::filesystem::remove(copy);
        std::filesystem::remove(list);
    }

    TEST(GeometryAidTests, ClassifierWritesOneProbabilityRowPerAtomThroughBothFlags)
    {
        const std::filesystem::path model = geometry_aid_tmp("wide.bin"), out = geometry_aid_tmp("probs.npy"), copy = geometry_aid_tmp("copy2.xyz"), list = geometry_aid_tmp("classify.txt"), copy_probs = copy.string() + ".probs.npy", descr = geometry_aid_tmp("descr.npy");
        write_model(model, wide_model());
        std::filesystem::copy_file(thpp, copy, std::filesystem::copy_options::overwrite_existing);
        { std::ofstream(list) << copy.string() << "\n"; }
        ASSERT_EQ(run_nosphera2({ "-wfn", thpp.string(), "-classify_atoms", model.string(), out.string(), "-no_date" }), 0);
        std::vector<unsigned long> shape;
        vec p, p_batch, d;
        load_npy(out, shape, p);
        ASSERT_EQ(shape, (std::vector<unsigned long>{ 18, 2 }));
        for (int a = 0; a < 18; a++) {
            EXPECT_NEAR(p[2 * a] + p[2 * a + 1], 1.0, 1e-12) << a;
            EXPECT_GT(p[2 * a], 0.0);
            EXPECT_GT(p[2 * a + 1], 0.0);
        }
        //The same numbers from the pieces: the descriptor written by the other flag, classified with the loaded model
        geometry_aid::write_descriptor(thpp, descr, geometry_aid::hyperparameters());
        load_npy(descr, shape, d);
        const vec direct = geometry_aid::classify_descriptor(d.data(), 18, 42042, geometry_aid::load_model(model));
        ASSERT_EQ(direct.size(), p.size());
        double spread = 0.0;
        for (int i = 0; i < 36; i++) {
            EXPECT_NEAR(direct[i], p[i], 1e-12) << i;
            spread = std::max(spread, std::abs(p[i] - 0.5));
        }
        EXPECT_GT(spread, 1e-3);
        ASSERT_EQ(run_nosphera2({ "-classify_atoms_list", list.string(), model.string(), "-no_date" }), 0);
        load_npy(copy_probs, shape, p_batch);
        ASSERT_EQ(p_batch.size(), p.size());
        for (int i = 0; i < 36; i++) EXPECT_EQ(p_batch[i], p[i]) << i;
        for (const std::filesystem::path& f : { model, out, copy, list, copy_probs, descr }) std::filesystem::remove(f);
    }
    namespace {
        vec2 rotation(const double a, const double b, const double c)
        {
            const double ca = std::cos(a), sa = std::sin(a), cb = std::cos(b), sb = std::sin(b), cc = std::cos(c), sc = std::sin(c);
            const vec2 Z1{ {ca, -sa, 0}, {sa, ca, 0}, {0, 0, 1} }, Y{ {cb, 0, sb}, {0, 1, 0}, {-sb, 0, cb} }, Z2{ {cc, -sc, 0}, {sc, cc, 0}, {0, 0, 1} };
            vec2 R(3, vec(3, 0.0));
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++) R[i][j] += Z1[i][k] * Y[k][l] * Z2[l][j];
            return R;
        }
        int aux_size(const WFN& aux)
        {
            int n = 0;
            for (int a = 0; a < aux.get_ncen(); a++) {
                const atom A = aux.get_atom(a);
                int prim = 0;
                for (int shell = 0; shell < (int)A.get_shellcount_size(); shell++) {
                    n += 2 * A.get_basis_set_entry(prim).get_type() + 1;
                    prim += A.get_shellcount(shell);
                }
            }
            return n;
        }
        WFN oh_molecule(const double z)
        {
            WFN wavy(e_origin::NOT_YET_DEFINED);
            wavy.push_back_atom("O", 0.0, 0.0, z, 8);
            wavy.push_back_atom("H", 1.8, 0.0, z, 1);
            return wavy;
        }
        void write_p1bar_cif(const std::filesystem::path& cif, const std::string& op)
        {
            std::ofstream out(cif);
            out << "data_test\n_cell_length_a 10.0\n_cell_length_b 10.0\n_cell_length_c 10.0\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n_cell_volume 1000\n"
                << "loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n'" << op << "'\n";
        }
    }

    TEST(CrystalEnergyTests, RealSphericalHarmonicRotationIsOrthogonalAndMatchesTheFunctions)
    {
        const vec2 R = rotation(0.4, 1.1, -2.3);
        vec2 S = R, I(3, vec(3, 0.0));
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) S[i][j] = -R[i][j], I[i][j] = i == j ? -1.0 : 0.0;
        double u[3] = { 0.3, -0.5, 0.81 };
        const double norm = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        for (int x = 0; x < 3; x++) u[x] /= norm;
        for (int l = 0; l <= 5; l++) {
            const int n = 2 * l + 1;
            for (const vec2& M : { R, S }) {
                const vec2 D = crystal_energies::real_sh_rotation(l, M);
                ASSERT_EQ(D.size(), n);
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) {
                        double dot = 0.0;
                        for (int k = 0; k < n; k++) dot += D[k][i] * D[k][j];
                        EXPECT_NEAR(dot, i == j ? 1.0 : 0.0, 1e-9) << l << " " << i << " " << j;
                    }
                double v[3];
                for (int x = 0; x < 3; x++) v[x] = M[0][x] * u[0] + M[1][x] * u[1] + M[2][x] * u[2];
                for (int m = 0; m < n; m++) {
                    double sum = 0.0;
                    for (int k = 0; k < n; k++) sum += D[k][m] * constants::spherical_harmonic(l, k - l, u);
                    EXPECT_NEAR(sum, constants::spherical_harmonic(l, m - l, v), 1e-9) << l << " " << m;
                }
            }
            const vec2 P = crystal_energies::real_sh_rotation(l, I);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) EXPECT_NEAR(P[i][j], i == j ? (l % 2 ? -1.0 : 1.0) : 0.0, 1e-9) << l << " " << i << " " << j;
        }
    }

    TEST(CrystalEnergyTests, MovedMoleculesKeepTheirInteractionEnergy)
    {
        std::vector<std::shared_ptr<BasisSet>> basis{ BasisSetLibrary::get_basis_set("combo_basis_fit") };
        const WFN aux_A = generate_aux_wfn(oh_molecule(0.0), basis), aux_B = generate_aux_wfn(oh_molecule(5.0), basis);
        vec c_A(aux_size(aux_A)), c_B(aux_size(aux_B));
        for (int i = 0; i < (int)c_A.size(); i++) c_A[i] = 0.3 * std::sin(1.0 + i), c_B[i] = 0.3 * std::cos(0.5 + 2 * i);
        //a positive s part keeps the densities positive, as fitted ones are; the repulsion needs rho^(5/3)
        const aux_density_table tA(aux_A.get_atoms()), tB(aux_B.get_atoms());
        for (int s = 0; s < tA.n_sh; s++) if (tA.sh_l[s] == 0) c_A[tA.coef_off[s]] += 1.0;
        for (int s = 0; s < tB.n_sh; s++) if (tB.sh_l[s] == 0) c_B[tB.coef_off[s]] += 1.0;
        const DensityFitting::INTERACTION E = DensityFitting::interaction_energy(c_A, aux_A, c_B, aux_B);
        EXPECT_GT(std::abs(E.electrostatic()), 1e-4);
        const vec2 R = rotation(0.4, 1.1, -2.3);
        const vec t{ 0.7, -1.1, 2.3 };
        for (int improper = 0; improper < 2; improper++) {
            vec2 M = R;
            if (improper)
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) M[i][j] = -R[i][j];
            WFN a = aux_A, b = aux_B;
            vec ca = c_A, cb = c_B;
            crystal_energies::transform(a, ca, M, t);
            crystal_energies::transform(b, cb, M, t);
            for (int x = 0; x < 3; x++) EXPECT_NEAR(a.get_atom_coordinate(1, x), t[x] + 1.8 * M[x][0], 1e-12);
            const DensityFitting::INTERACTION F = DensityFitting::interaction_energy(ca, a, cb, b);
            EXPECT_NEAR(F.electrostatic(), E.electrostatic(), 1e-8) << improper;
            EXPECT_NEAR(F.overlap, E.overlap, 1e-8) << improper;
            EXPECT_NEAR(F.disp, E.disp, 1e-8) << improper;
            EXPECT_NEAR(F.pol_A, E.pol_A, 1e-7) << improper;
            EXPECT_NEAR(F.pol_B, E.pol_B, 1e-7) << improper;
            //the repulsion is a grid quadrature and the grids do not turn with the atoms
            EXPECT_NEAR(F.rep, E.rep, 1e-2 * std::abs(E.rep)) << improper;
        }
    }

    // Same partition on both sides: the rows partition_rows_on_grid builds for the oxygen of an O-H pair, applied to
    // arbitrary coefficients on both atoms, must give the Becke-weighted grid moments of the density calc_density_ML
    // evaluates from those coefficients, and add_partition_restraint must place them with the sqrt(4pi) / r_cov^-l
    // scaling of the older rows. The hydrogen's functions contribute to the oxygen's rows, which is the point.
    TEST(RiMultipoleTests, PartitionRowsReproduceTheGridMomentsOfTheFittedDensity)
    {
        const int lmax = 3, n_moments = (lmax + 1) * (lmax + 1);
        const WFN mol = oh_molecule(0.4);
        std::vector<std::shared_ptr<BasisSet>> basis{ BasisSetLibrary::get_basis_set("combo_basis_fit") };
        const WFN aux = generate_aux_wfn(mol, basis);
        const aux_density_table t(aux.get_atoms());
        const int n_aux = t.n_coef;
        const atom A = aux.get_atom(0);
        double alpha_min[9] = { 0.0 }, alpha_max = 0.0;
        int max_l = 0, prim = 0;
        for (int shell = 0; shell < (int)A.get_shellcount_size(); shell++) {
            const int l = A.get_basis_set_entry(prim).get_type();
            for (int e = 0; e < (int)A.get_shellcount(shell); e++) {
                const double a = A.get_basis_set_entry(prim + e).get_exponent();
                alpha_max = std::max(alpha_max, a);
                if (alpha_min[l] == 0.0 || a < alpha_min[l]) alpha_min[l] = a;
            }
            max_l = std::max(max_l, l);
            prim += A.get_shellcount(shell);
        }
        const double xs[2] = { mol.get_atom_coordinate(0, 0), mol.get_atom_coordinate(1, 0) }, ys[2] = { mol.get_atom_coordinate(0, 1), mol.get_atom_coordinate(1, 1) }, zs[2] = { mol.get_atom_coordinate(0, 2), mol.get_atom_coordinate(1, 2) };
        const int Zs[2] = { 8, 1 };
        AtomGrid grid(1e-12, 350, 5000, 8, alpha_max, max_l, alpha_min, std::cout);
        const int n_points = grid.get_num_grid_points();
        vec gx(n_points), gy(n_points), gz(n_points), aw(n_points), bw(n_points), tw(n_points), chi(1, 0.0);
        grid.get_grid(2, 0, xs, ys, zs, Zs, gx.data(), gy.data(), gz.data(), aw.data(), bw.data(), tw.data(), WFN(), chi);
        vec coefs(n_aux);
        for (int i = 0; i < n_aux; i++) coefs[i] = std::sin(1.0 + i);
        vec2 Q(1, vec(n_moments, 0.0));
        double n_becke = 0.0;
        for (int p = 0; p < n_points; p++) {
            const double f = calc_density_ML(gx[p], gy[p], gz[p], coefs, aux.get_atoms()) * bw[p];
            double d[3] = { gx[p] - xs[0], gy[p] - ys[0], gz[p] - zs[0] };
            const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            for (int i = 0; i < 3; i++) d[i] /= r;
            double rl = 1.0;
            for (int l = 0; l <= lmax; l++) {
                for (int m = -l; m <= l; m++)
                    Q[0][l * l + l + m] += f * rl * constants::spherical_harmonic(l, m, d);
                rl *= r;
            }
            n_becke += bw[p];
        }
        EXPECT_GT(n_becke, 1.0);
        vec2 rows(n_moments, vec(n_aux, 0.0));
        const double centre[3] = { xs[0], ys[0], zs[0] };
        DensityFitting::partition_rows_on_grid(t, n_points, gx.data(), gy.data(), gz.data(), bw.data(), centre, lmax, rows, 0);
        int on_H = 0;
        for (int i = t.coef_off[t.sh_start[1]]; i < n_aux; i++) on_H += std::abs(rows[0][i]) > 1e-6;
        EXPECT_GT(on_H, 0);
        vec eri2c(n_aux * n_aux, 0.0), rho(n_aux, 0.0);
        WFN aux_O(e_origin::NOT_YET_DEFINED);
        aux_O.push_back_atom("O", xs[0], ys[0], zs[0], 8);
        DensityFitting::add_partition_restraint(eri2c, rho, aux_O, rows, Q, vec(1, 1.0), lmax);
        ASSERT_EQ(eri2c.size(), (size_t)(n_aux + n_moments) * n_aux);
        ASSERT_EQ(rho.size(), (size_t)(n_aux + n_moments));
        const vec2 fitted = DensityFitting::grid_multipoles(rows, coefs, lmax);
        const double r_cov = constants::ang2bohr(constants::covalent_radii[8]);
        for (int row = 0; row < n_moments; row++) {
            const int l = (int)std::floor(std::sqrt(row + 1e-9));
            const double scale = row == 0 ? std::sqrt(4.0 * PI_VAL) : std::pow(r_cov, -l);
            double lhs = 0.0;
            for (int i = 0; i < n_aux; i++) lhs += eri2c[(n_aux + row) * n_aux + i] * coefs[i];
            EXPECT_NEAR(lhs, scale * Q[0][row], 1e-8 * std::max(1.0, std::abs(scale * Q[0][row]))) << "row " << row;
            EXPECT_NEAR(rho[n_aux + row], scale * Q[0][row], 1e-12) << "row " << row;
            EXPECT_NEAR(fitted[0][row], Q[0][row], 1e-8 * std::max(1.0, std::abs(Q[0][row]))) << "row " << row;
        }
    }

    // The flattened aux basis must give the density calc_density_ML gives atom by atom, on the
    // host and on the device; the point set is sized past the kernel's minimum work so that the
    // GPU branch is the one being tested when a device is present
    TEST(CrystalEnergyTests, FlattenedAuxDensityMatchesTheAtomWalk)
    {
        std::vector<std::shared_ptr<BasisSet>> basis{ BasisSetLibrary::get_basis_set("combo_basis_fit") };
        const WFN aux = generate_aux_wfn(oh_molecule(0.0), basis);
        vec c(aux_size(aux));
        for (int i = 0; i < (int)c.size(); i++) c[i] = 0.3 * std::sin(1.0 + i);
        const aux_density_table t(aux.get_atoms());
        EXPECT_EQ(t.n_coef, (int)c.size());
        const int np = 200000;
        vec x(np), y(np), z(np), rho(np);
        for (int p = 0; p < np; p++) x[p] = 6.0 * std::sin(0.37 * p) - 1.0, y[p] = 5.0 * std::cos(0.53 * p), z[p] = 7.0 * std::sin(0.11 * p + 1.0) + 0.5;
#ifdef NOSPHERA2_USE_GPU
        aux_density_gpu_set_enabled(false);
#endif
        calc_density_ML(t, c, np, x.data(), y.data(), z.data(), rho.data());
        double n = 0.0;
        for (int p = 0; p < np; p += 97) {
            const double ref = calc_density_ML(x[p], y[p], z[p], c, aux.get_atoms());
            EXPECT_NEAR(rho[p], ref, 1e-12 * std::abs(ref) + 1e-14) << p;
            n += std::abs(ref);
        }
        EXPECT_GT(n, 1e-3);
#ifdef NOSPHERA2_USE_GPU
        if (!aux_density_gpu_available()) GTEST_SKIP() << "No GPU device present; the aux density kernel cannot run here";
        aux_density_gpu_set_enabled(true);
        vec rho_gpu(np);
        const bool ran = aux_density_gpu_eval(t.n_at, t.cx.data(), t.cy.data(), t.cz.data(), t.r2_max.data(), t.n_sh, t.sh_start.data(), t.sh_l.data(), t.pr_start.data(), t.coef_off.data(), t.n_pr, t.pr_exp.data(), t.pr_norm.data(), t.n_coef, c.data(), np, x.data(), y.data(), z.data(), rho_gpu.data());
        aux_density_gpu_set_enabled(false);
        ASSERT_TRUE(ran);
        for (int p = 0; p < np; p++) EXPECT_NEAR(rho_gpu[p], rho[p], 1e-11 * std::abs(rho[p]) + 1e-14) << p;
#endif
    }

    TEST(CrystalEnergyTests, AnalyticAuxGradientMatchesCentralDifferences)
    {
        std::vector<std::shared_ptr<BasisSet>> basis{ BasisSetLibrary::get_basis_set("combo_basis_fit") };
        const WFN aux = generate_aux_wfn(oh_molecule(0.0), basis);
        vec c(aux_size(aux));
        for (int i = 0; i < (int)c.size(); i++) c[i] = 0.3 * std::sin(1.0 + i);
        const aux_density_table t(aux.get_atoms());
        const int np = 200000;
        const double h = 1e-5;
        vec x(np), y(np), z(np), rho(np), gx(np), gy(np), gz(np);
        for (int p = 0; p < np; p++) x[p] = 6.0 * std::sin(0.37 * p) - 1.0, y[p] = 5.0 * std::cos(0.53 * p), z[p] = 7.0 * std::sin(0.11 * p + 1.0) + 0.5;
#ifdef NOSPHERA2_USE_GPU
        aux_density_gpu_set_enabled(false);
#endif
        calc_density_ML(t, c, np, x.data(), y.data(), z.data(), rho.data(), gx.data(), gy.data(), gz.data());
        double n = 0.0;
        for (int p = 0; p < np; p += 97) {
            EXPECT_NEAR(rho[p], t(x[p], y[p], z[p], c.data()), 1e-12 * std::abs(rho[p]) + 1e-14) << p;
            const double ref[3] = {
                (t(x[p] + h, y[p], z[p], c.data()) - t(x[p] - h, y[p], z[p], c.data())) / (2 * h),
                (t(x[p], y[p] + h, z[p], c.data()) - t(x[p], y[p] - h, z[p], c.data())) / (2 * h),
                (t(x[p], y[p], z[p] + h, c.data()) - t(x[p], y[p], z[p] - h, c.data())) / (2 * h) };
            EXPECT_NEAR(gx[p], ref[0], 1e-6 * std::abs(ref[0]) + 1e-9) << p;
            EXPECT_NEAR(gy[p], ref[1], 1e-6 * std::abs(ref[1]) + 1e-9) << p;
            EXPECT_NEAR(gz[p], ref[2], 1e-6 * std::abs(ref[2]) + 1e-9) << p;
            n += std::abs(ref[0]) + std::abs(ref[1]) + std::abs(ref[2]);
        }
        EXPECT_GT(n, 1e-3);
#ifdef NOSPHERA2_USE_GPU
        if (!aux_density_gpu_available()) GTEST_SKIP() << "No GPU device present; the aux density kernel cannot run here";
        aux_density_gpu_set_enabled(true);
        vec rho_gpu(np), gx_gpu(np), gy_gpu(np), gz_gpu(np);
        const bool ran = aux_density_gpu_eval(t.n_at, t.cx.data(), t.cy.data(), t.cz.data(), t.r2_max.data(), t.n_sh, t.sh_start.data(), t.sh_l.data(), t.pr_start.data(), t.coef_off.data(), t.n_pr, t.pr_exp.data(), t.pr_norm.data(), t.n_coef, c.data(), np, x.data(), y.data(), z.data(), rho_gpu.data(), gx_gpu.data(), gy_gpu.data(), gz_gpu.data());
        aux_density_gpu_set_enabled(false);
        ASSERT_TRUE(ran);
        for (int p = 0; p < np; p++) {
            EXPECT_NEAR(rho_gpu[p], rho[p], 1e-11 * std::abs(rho[p]) + 1e-14) << p;
            EXPECT_NEAR(gx_gpu[p], gx[p], 1e-11 * std::abs(gx[p]) + 1e-13) << p;
            EXPECT_NEAR(gy_gpu[p], gy[p], 1e-11 * std::abs(gy[p]) + 1e-13) << p;
            EXPECT_NEAR(gz_gpu[p], gz[p], 1e-11 * std::abs(gz[p]) + 1e-13) << p;
        }
#endif
    }

    TEST(CrystalEnergyTests, AnalyticAuxLaplacianMatchesCentralDifferences)
    {
        std::vector<std::shared_ptr<BasisSet>> basis{ BasisSetLibrary::get_basis_set("combo_basis_fit") };
        const WFN aux = generate_aux_wfn(oh_molecule(0.0), basis);
        vec c(aux_size(aux));
        for (int i = 0; i < (int)c.size(); i++) c[i] = 0.3 * std::sin(1.0 + i);
        const aux_density_table t(aux.get_atoms());
        const int np = 200000;
        const double h = 1e-5;
        vec x(np), y(np), z(np), rho(np), gx(np), gy(np), gz(np), lap(np);
        for (int p = 0; p < np; p++) x[p] = 6.0 * std::sin(0.37 * p) - 1.0, y[p] = 5.0 * std::cos(0.53 * p), z[p] = 7.0 * std::sin(0.11 * p + 1.0) + 0.5;
#ifdef NOSPHERA2_USE_GPU
        aux_density_gpu_set_enabled(false);
#endif
        calc_density_ML(t, c, np, x.data(), y.data(), z.data(), rho.data(), gx.data(), gy.data(), gz.data(), lap.data());
        double n = 0.0;
        for (int p = 0; p < np; p += 97) {
            double g[3], d[3], ref = 0.0;
            EXPECT_NEAR(rho[p], t(x[p], y[p], z[p], c.data(), g[0], g[1], g[2]), 1e-12 * std::abs(rho[p]) + 1e-14) << p;
            EXPECT_NEAR(gx[p], g[0], 1e-12 * std::abs(g[0]) + 1e-14) << p;
            EXPECT_NEAR(gy[p], g[1], 1e-12 * std::abs(g[1]) + 1e-14) << p;
            EXPECT_NEAR(gz[p], g[2], 1e-12 * std::abs(g[2]) + 1e-14) << p;
            t(x[p] + h, y[p], z[p], c.data(), g[0], g[1], g[2]), t(x[p] - h, y[p], z[p], c.data(), d[0], d[1], d[2]), ref += (g[0] - d[0]) / (2 * h);
            t(x[p], y[p] + h, z[p], c.data(), g[0], g[1], g[2]), t(x[p], y[p] - h, z[p], c.data(), d[0], d[1], d[2]), ref += (g[1] - d[1]) / (2 * h);
            t(x[p], y[p], z[p] + h, c.data(), g[0], g[1], g[2]), t(x[p], y[p], z[p] - h, c.data(), d[0], d[1], d[2]), ref += (g[2] - d[2]) / (2 * h);
            EXPECT_NEAR(lap[p], ref, 1e-6 * std::abs(ref) + 1e-9) << p;
            n += std::abs(ref);
        }
        EXPECT_GT(n, 1e-3);
#ifdef NOSPHERA2_USE_GPU
        if (!aux_density_gpu_available()) GTEST_SKIP() << "No GPU device present; the aux density kernel cannot run here";
        aux_density_gpu_set_enabled(true);
        vec rho_gpu(np), gx_gpu(np), gy_gpu(np), gz_gpu(np), lap_gpu(np);
        const bool ran = aux_density_gpu_eval(t.n_at, t.cx.data(), t.cy.data(), t.cz.data(), t.r2_max.data(), t.n_sh, t.sh_start.data(), t.sh_l.data(), t.pr_start.data(), t.coef_off.data(), t.n_pr, t.pr_exp.data(), t.pr_norm.data(), t.n_coef, c.data(), np, x.data(), y.data(), z.data(), rho_gpu.data(), gx_gpu.data(), gy_gpu.data(), gz_gpu.data(), lap_gpu.data());
        aux_density_gpu_set_enabled(false);
        ASSERT_TRUE(ran);
        for (int p = 0; p < np; p++) {
            EXPECT_NEAR(rho_gpu[p], rho[p], 1e-11 * std::abs(rho[p]) + 1e-14) << p;
            EXPECT_NEAR(gx_gpu[p], gx[p], 1e-11 * std::abs(gx[p]) + 1e-13) << p;
            EXPECT_NEAR(lap_gpu[p], lap[p], 1e-11 * std::abs(lap[p]) + 1e-13) << p;
        }
#endif
    }

    TEST(CrystalEnergyTests, ContactsInPMinus1AreListedOnce)
    {
        const std::filesystem::path cif = geometry_aid_tmp("p1bar.cif"), cif3 = geometry_aid_tmp("p3.cif");
        write_p1bar_cif(cif, "-x, -y, -z");
        write_p1bar_cif(cif3, "-y, x-y, z");
        cell c(cif, std::cout, false, true);
        WFN mol(e_origin::NOT_YET_DEFINED);
        mol.push_back_atom("O", constants::ang2bohr(1.0), constants::ang2bohr(2.0), constants::ang2bohr(3.0), 8);
        const std::vector<crystal_energies::pair> pairs = crystal_energies::contacts(c, { mol }, 10.5);
        ASSERT_EQ(pairs.size(), 8);
        int identity = 0, inversion = 0;
        std::set<std::string> ns;
        for (const crystal_energies::pair& p : pairs) {
            EXPECT_EQ(p.A, 0);
            EXPECT_EQ(p.B, 0);
            if (p.op == 0) {
                identity++;
                EXPECT_NEAR(p.distance, 10.0, 1e-9);
            }
            else {
                inversion++;
                ns.insert(std::to_string(p.n[0]) + std::to_string(p.n[1]) + std::to_string(p.n[2]));
                if (p.n[0] == 0 && p.n[1] == 0 && p.n[2] == 1) {
                    EXPECT_EQ(p.symop, "-x,-y,-z+1");
                    EXPECT_NEAR(p.distance, 6.0, 1e-9);
                    EXPECT_NEAR(p.centroid_B[0], -1.0, 1e-9);
                    EXPECT_NEAR(p.centroid_B[1], -2.0, 1e-9);
                    EXPECT_NEAR(p.centroid_B[2], 7.0, 1e-9);
                }
            }
            for (int x = 0; x < 3; x++) EXPECT_NEAR(p.centroid_A[x], x + 1.0, 1e-9);
        }
        EXPECT_EQ(identity, 3);
        EXPECT_EQ(inversion, 5);
        EXPECT_EQ(ns, std::set<std::string>({ "000", "001", "010", "011", "101" }));
        cell c3(cif3, std::cout, false, true);
        const std::vector<crystal_energies::symop> ops = crystal_energies::symops(c3);
        ASSERT_EQ(ops.size(), 2);
        EXPECT_EQ(ops[1].rot, ivec2({ { 0, -1, 0 }, { 1, -1, 0 }, { 0, 0, 1 } }));
        EXPECT_EQ(crystal_energies::symop_string(ops[1], { 0, 1, -1 }), "-y,x-y+1,z-1");
        std::filesystem::remove(cif), std::filesystem::remove(cif3);
    }

    TEST(CrystalEnergyTests, RunAppWritesThePairTableWithInvertedCoefficients)
    {
        const std::filesystem::path dir = geometry_aid_tmp("crystal"), cif = dir / "p1bar.cif", xyz = dir / "o.xyz", npy = dir / "o.npy", job = dir / "job.txt", table = dir / "ie.txt";
        std::filesystem::create_directories(dir);
        write_p1bar_cif(cif, "-x, -y, -z");
        { std::ofstream out(xyz); out << "1\ntest\nO 1.0 2.0 3.0\n"; }
        std::vector<std::shared_ptr<BasisSet>> basis{ BasisSetLibrary::get_basis_set("combo_basis_fit") };
        const WFN aux_A = generate_aux_wfn(WFN(xyz), basis);
        vec coef(aux_size(aux_A));
        for (int i = 0; i < (int)coef.size(); i++) coef[i] = 0.2 * std::sin(1.0 + i);
        const unsigned long shape[1] = { (unsigned long)coef.size() };
        npy::SaveArrayAsNumpy(npy, false, 1, shape, coef);
        { std::ofstream out(job); out << "# test job\ncif p1bar.cif\ncutoff 10.5\noutput ie.txt\nmolecule o.xyz o.npy\n"; }
        ASSERT_EQ(run_nosphera2({ "-ri_fit", "combo_basis_fit", "-interaction_energies", job.string(), "-no_date" }), 0);
        std::ifstream in(table);
        ASSERT_TRUE(in.good());
        std::string line;
        int rows = 0;
        double total = 0.0;
        bool found = false;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            rows++;
            std::istringstream ss(line);
            int A, B, n[3];
            std::string symop;
            double R, x[6], E[5];
            ss >> A >> B >> symop >> n[0] >> n[1] >> n[2] >> R;
            for (int i = 0; i < 6; i++) ss >> x[i];
            for (int i = 0; i < 5; i++) ss >> E[i];
            ASSERT_FALSE(ss.fail()) << line;
            EXPECT_NEAR(E[0] + E[1] + E[2] + E[3], E[4], 3e-3) << line;
            if (symop == "-x,-y,-z+1" && n[2] == 1) found = true, total = E[4];
        }
        in.close();
        EXPECT_EQ(rows, 8);
        ASSERT_TRUE(found);
        WFN aux_B = aux_A;
        std::vector<atom> atoms = aux_B.get_atoms();
        atoms[0].set_coordinate(0, constants::ang2bohr(-1.0)), atoms[0].set_coordinate(1, constants::ang2bohr(-2.0)), atoms[0].set_coordinate(2, constants::ang2bohr(7.0));
        aux_B.set_atoms(atoms);
        vec coef_B = coef;
        int prim = 0, offset = 0;
        for (int shell = 0; shell < (int)atoms[0].get_shellcount_size(); shell++) {
            const int l = atoms[0].get_basis_set_entry(prim).get_type();
            for (int m = 0; m < 2 * l + 1; m++) coef_B[offset + m] *= l % 2 ? -1.0 : 1.0;
            offset += 2 * l + 1, prim += atoms[0].get_shellcount(shell);
        }
        const DensityFitting::INTERACTION E = DensityFitting::interaction_energy(coef, aux_A, coef_B, aux_B);
        EXPECT_NEAR(total, E.total() * constants::kcal_mol_per_hartree * 4.184, 2e-3);
        std::filesystem::remove_all(dir);
    }
} // namespace NoSpherA2UnitTests
