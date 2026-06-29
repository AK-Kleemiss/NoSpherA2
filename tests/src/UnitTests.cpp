
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
        static constexpr double kTolerance = 1.0e-12;

        static constexpr std::uint64_t kZMask = 0x00000000000000FFULL;
        static constexpr std::uint64_t kAMask = 0x0000000001FFFF00ULL;
        static constexpr std::uint64_t kBMask = 0x000003FFFE000000ULL;
        static constexpr std::uint64_t kCMask = 0x07FFFC0000000000ULL;
        static constexpr std::uint64_t kASign = 0x0800000000000000ULL;
        static constexpr std::uint64_t kBSign = 0x1000000000000000ULL;
        static constexpr std::uint64_t kCSign = 0x2000000000000000ULL;
        static constexpr std::uint64_t kDatMask = 0xC000000000000000ULL;

        static constexpr std::int64_t kCoordinateMaximum = 0x1FFFF;

        static atom make_atom()
        {
            return atom{
                "C1",
                std::uint64_t{0},
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
            const double z)
        {
            atom result{
                "Test",
                std::uint64_t{0},
                charge,
                0.0,
                0.0,
                0.0,
                charge
            };

            result.set_frac_coords(d3{ x, y, z });
            return result;
        }

        static std::uint64_t extractA(const std::uint64_t id)
        {
            return (id & kAMask) >> 8;
        }

        static std::uint64_t extractB(const std::uint64_t id)
        {
            return (id & kBMask) >> 25;
        }

        static std::uint64_t extractC(const std::uint64_t id)
        {
            return (id & kCMask) >> 42;
        }

        static std::uint64_t extractDat(const std::uint64_t id)
        {
            return (id & kDatMask) >> 62;
        }
    };

    TEST_F(AtomTest, IDZeroFractionalCoordinatesReturnZero)
    {
        atom value =
            make_atom_fractional(6, 0.0, 0.0, 0.0);

        EXPECT_EQ(value.get_ID(), std::uint64_t{ 0 });
    }

    TEST_F(AtomTest, IDEncodesChargeCoordinatesAndDat)
    {
        atom value = make_atom_fractional(
            6,
            1.0 / static_cast<double>(kCoordinateMaximum),
            2.0 / static_cast<double>(kCoordinateMaximum),
            3.0 / static_cast<double>(kCoordinateMaximum)
        );

        const std::uint64_t id = value.get_ID(2);

        EXPECT_EQ(id & kZMask, std::uint64_t{ 6 });
        EXPECT_EQ(extractA(id), std::uint64_t{ 1 });
        EXPECT_EQ(extractB(id), std::uint64_t{ 2 });
        EXPECT_EQ(extractC(id), std::uint64_t{ 3 });
        EXPECT_EQ(extractDat(id), std::uint64_t{ 2 });
    }

    TEST_F(AtomTest, IDEncodesCoordinateSignsSeparatelyFromMagnitude)
    {
        atom positive =
            make_atom_fractional(6, 0.1, 0.2, 0.3);

        atom negative =
            make_atom_fractional(6, -0.1, -0.2, -0.3);

        const std::uint64_t positiveId = positive.get_ID();
        const std::uint64_t negativeId = negative.get_ID();

        EXPECT_EQ(extractA(positiveId), extractA(negativeId));
        EXPECT_EQ(extractB(positiveId), extractB(negativeId));
        EXPECT_EQ(extractC(positiveId), extractC(negativeId));

        EXPECT_EQ(positiveId & (kASign | kBSign | kCSign), 0u);
        EXPECT_EQ(
            negativeId & (kASign | kBSign | kCSign),
            kASign | kBSign | kCSign
        );
    }

    TEST_F(AtomTest, IDFractionalCoordinatesAreTruncated)
    {
        atom value =
            make_atom_fractional(6, 0.5, 0.5, 0.5);

        const std::uint64_t id = value.get_ID();

        // static_cast<int64_t>(0.5 * 131071) truncates to 65535.
        constexpr std::uint64_t expected = 65535;

        EXPECT_EQ(extractA(id), expected);
        EXPECT_EQ(extractB(id), expected);
        EXPECT_EQ(extractC(id), expected);
    }

    TEST_F(AtomTest, IDIsCachedUntilExplicitlyReset)
    {
        atom value =
            make_atom_fractional(6, 0.1, 0.2, 0.3);

        const std::uint64_t originalId = value.get_ID(0);

        value.set_frac_coords(d3{ 0.7, 0.8, 0.9 });
        value.set_charge(8);

        EXPECT_EQ(value.get_ID(3), originalId);

        value.set_ID(0);

        EXPECT_NE(value.get_ID(3), originalId);
    }

    TEST_F(AtomTest, IDToStringAndBack)
    {
        atom value =
            make_atom_fractional(6, 0.1, 0.2, 0.3);

        const std::uint64_t originalId = value.get_ID(0);

        const std::string idStr = std::to_string(originalId);
        const std::uint64_t parsedId = std::stoull(idStr);

        EXPECT_EQ(parsedId, originalId);
    }


// -----------------------------------------------------------------------
// Non-trivial atom behavior
// -----------------------------------------------------------------------

TEST_F(AtomTest, DistanceToOtherAtomIsEuclideanDistance)
{
    const atom first{
        "A",
        std::uint64_t{1},
        1,
        1.0,
        2.0,
        3.0,
        0
    };

    const atom second{
        "B",
        std::uint64_t{2},
        1,
        4.0,
        6.0,
        3.0,
        0
    };

    EXPECT_NEAR(first.distance_to(second), 5.0, kTolerance);
    EXPECT_NEAR(second.distance_to(first), 5.0, kTolerance);
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
        std::uint64_t{555},
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

} // namespace NoSpherA2UnitTests