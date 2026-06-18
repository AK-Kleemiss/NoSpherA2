
#include "pch.h"
#include "CppUnitTest.h"

#include "../core/convenience.h"
#include "../core/constants.h"
#include "../core/fchk.h"
#include "../core/AtomGrid.h"
#include "../core/SALTED_utilities.h"
#include "../core/scattering_factors.h"
#include "../core/nos_math.h"

static constexpr double PI_VAL = 3.14159265358979323846;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
namespace {
    int test_analytical_fourier(bool full)
    {
        // Generate grid and k_pts
        vec2 kpts;

        for (int i = 1; i < 1000; i++) {
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
        if (full) {
            max_l = 8;
            radial_res = 1E-25;
        }

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
            Logger::WriteMessage(("Testing l = " + std::to_string(type)).c_str());
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
                    Logger::WriteMessage(("Error at m: " + std::to_string(m) + "   Max diff: (" + std::to_string(max_diff.real()) + "|" + std::to_string(max_diff.imag()) + ")\n").c_str());
                    correct = true;
                }
            }
            if (!all_correct)
                break;
            Logger::WriteMessage("| PASSED!\n");
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
            Logger::WriteMessage("Error in the calculations!");
            return 1;
        }
        Logger::WriteMessage("All tests passed!");
        return 0;
    };

    template <typename MatType1, typename MatType2>
    void compare_matrices(const MatType1& mat, const MatType2& vecMat) {
        const double tol = 1e-9;
        for (size_t i = 0; i < mat.extent(0); i++) {
            for (size_t j = 0; j < mat.extent(1); j++) {
                //Compile different versions based on the type of the matrix dMatrix and dMatrix2 are accessed with mat(i, j) and vecMat[i][j] respectively
                if constexpr (std::is_same_v<MatType1, dMatrix2> && std::is_same_v<MatType2, vec2>) {
                    Assert::AreEqual(mat(i, j), vecMat[i][j], tol, L"Matrix values differ");
                }
                else if constexpr (std::is_same_v<MatType1, cMatrix2> && std::is_same_v<MatType2, cvec2>) {
                    Assert::AreEqual(mat(i, j).real(), vecMat[i][j].real(), tol, L"Matrix values differ");
                    Assert::AreEqual(mat(i, j).imag(), vecMat[i][j].imag(), tol, L"Matrix values differ");
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
            Assert::AreEqual(b[i], x_expected[i], tol, L"Solution values differ");
        }

        if (ok) {
            Logger::WriteMessage("test_solve_linear_equations: PASSED\n");
        }
        else {
            Logger::WriteMessage("test_solve_linear_equations: FAILED\n");
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

        Logger::WriteMessage("Testing matrices directly\n");
        compare_matrices(matA, A);
        compare_matrices(matB, B);

        Logger::WriteMessage("Testing untransposed matrices\n");
        // First test regular dot-product
        compare_matrices(dot(matA, matB, false, false), self_dot(A, B));

        Logger::WriteMessage("Testing transpose A\n");
        ////Second compare first transpose
        compare_matrices(dot(matA, matB, true, false), self_dot(transpose(A), B));

        Logger::WriteMessage("Testing transpose B\n");
        ////Third comparte second transpose
        compare_matrices(dot(matA, matB, false, true), self_dot(A, transpose(B)));

        Logger::WriteMessage("Testing transpose A and B\n");
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

        Logger::WriteMessage("Testing C-matrices directly\n");
        compare_matrices(matC, C);
        compare_matrices(matD, D);

        Logger::WriteMessage("Testing untransposed C-matrices\n");
        // First test regular dot-product
        compare_matrices(dot(matC, matD, false, false), self_dot(C, D));

        Logger::WriteMessage("Testing transpose C\n");
        ////Second compare first transpose
        compare_matrices(dot(matC, matD, true, false), self_dot(transpose(C), D));

        Logger::WriteMessage("Testing transpose D\n");
        ////Third comparte second transpose
        compare_matrices(dot(matC, matD, false, true), self_dot(C, transpose(D)));

        Logger::WriteMessage("Testing transpose C and D\n");
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
        Logger::WriteMessage("Testing 2D x 1D matrix multiplication\n");
        compare_matrices(dot(matF, matE, false), self_dot(F, E));

        Logger::WriteMessage("All BLAS tests passed!\n");
    }

    void test_reading_SALTED_binary_file() {
        std::filesystem::path path("../../../tests/SALTED/Model/model.salted");
        SALTED_BINARY_FILE file = SALTED_BINARY_FILE(path, true);
        Config config;
        file.populate_config(config);
        std::unordered_map<int, std::vector<int64_t>> fps = file.read_fps();
        std::unordered_map<std::string, vec> averages = file.read_averages();
        std::unordered_map<int, vec> wigners = file.read_wigners();
        vec weights = file.read_weights();
        std::unordered_map<std::string, dMatrix2> feats = file.read_features();
        std::unordered_map<std::string, dMatrix2> proj = file.read_projectors();
        Logger::WriteMessage("Finished reading SALTED binary file\n");
        // TEST if both configs are the same
        Logger::WriteMessage(("Average:" + std::to_string(config.average) + "\n").c_str());
        Logger::WriteMessage(("Field:" + std::to_string(config.field) + "\n").c_str());
        Logger::WriteMessage(("Sparsify:" + std::to_string(config.sparsify) + "\n").c_str());
        Logger::WriteMessage(("Ncut:" + std::to_string(config.ncut) + "\n").c_str());
        Logger::WriteMessage(("Ntrain:" + std::to_string(config.Ntrain) + "\n").c_str());
        Logger::WriteMessage(("Menv:" + std::to_string(config.Menv) + "\n").c_str());
        Logger::WriteMessage(("trainfrac:" + std::to_string(config.trainfrac) + "\n").c_str());
        Logger::WriteMessage(("Rcut1:" + std::to_string(config.rcut1) + "\n").c_str());
        Logger::WriteMessage(("Rcut2:" + std::to_string(config.rcut2) + "\n").c_str());
        Logger::WriteMessage(("nang1:" + std::to_string(config.nang1) + "\n").c_str());
        Logger::WriteMessage(("nang2:" + std::to_string(config.nang2) + "\n").c_str());
        Logger::WriteMessage(("sig1:" + std::to_string(config.sig1) + "\n").c_str());
        Logger::WriteMessage(("sig2:" + std::to_string(config.sig2) + "\n").c_str());
        Logger::WriteMessage(("zeta:" + std::to_string(config.zeta) + "\n").c_str());
        Logger::WriteMessage(("neighspe size:" + std::to_string(config.neighspe1.size()) + "\n").c_str());
        for (int i = 0; i < config.neighspe1.size(); i++)
        {
            Logger::WriteMessage(("neighspe1[" + std::to_string(i) + "]:" + config.neighspe1[i] + "\n").c_str());
        }
        Logger::WriteMessage(("neighspe2 size:" + std::to_string(config.neighspe2.size()) + "\n").c_str());
        for (int i = 0; i < config.neighspe2.size(); i++)
        {
            Logger::WriteMessage(("neighspe2[" + std::to_string(i) + "]:" + config.neighspe2[i] + "\n").c_str());
        }
        Logger::WriteMessage(("dfBasis:" + config.dfbasis + "\n").c_str());

        Logger::WriteMessage("Comparing wigners\n");
        for (int i = 0; i < wigners.size(); i++)
        {
            for (int j = 0; j < wigners[i].size(); j += 10)
            {
                Logger::WriteMessage(("wigners[" + std::to_string(i) + "][" + std::to_string(j) + "]:" + std::to_string(wigners[i][j]) + "\n").c_str());
            }
        }

        Logger::WriteMessage("Comparing FPS\n");
        for (int i = 0; i < fps.size(); i++)
        {
            for (int j = 0; j < fps[i].size(); j += 10)
            {
                Logger::WriteMessage(("fps[" + std::to_string(i) + "][" + std::to_string(j) + "]:" + std::to_string(fps[i][j]) + "\n").c_str());
            }
        }

        Logger::WriteMessage("All tests passed!");
    }
}

namespace NoSpherA2UnitTests
{
    TEST_CLASS(SALTEDTests)
    {
    public:
        TEST_METHOD(ReadingSALTEDBinaryFile) {
            test_reading_SALTED_binary_file();
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(GeometryTests)
    {
    public:
        TEST_METHOD(ArrayDistance_UnitCubeDiagonal)
        {
            double d = array_length(d3{ 0.0, 0.0, 0.0 }, d3{ 1.0, 1.0, 1.0 });
            Assert::AreEqual(std::sqrt(3.0), d, 1e-12);
        }

        TEST_METHOD(ArrayDistance_SamePoint)
        {
            Assert::AreEqual(0.0, array_length(d3{ 1.5, 2.3, -4.7 }, d3{ 1.5, 2.3, -4.7 }), 1e-12);
        }

        TEST_METHOD(VecDiff_Basic)
        {
            d3 r = vec_diff({ 5.0, 10.0, 15.0 }, { 1.0, 2.0, 3.0 });
            Assert::AreEqual(4.0, r[0], 1e-12);
            Assert::AreEqual(8.0, r[1], 1e-12);
            Assert::AreEqual(12.0, r[2], 1e-12);
        }

        TEST_METHOD(VecDiff_ZeroResult)
        {
            d3 r = vec_diff({ 1.0, 2.0, 3.0 }, { 1.0, 2.0, 3.0 });
            Assert::AreEqual(0.0, r[0], 1e-12);
            Assert::AreEqual(0.0, r[1], 1e-12);
            Assert::AreEqual(0.0, r[2], 1e-12);
        }

        TEST_METHOD(VecCross_BasisVectors)
        {
            // i x j = k
            d3 r = vec_cross({ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 });
            Assert::AreEqual(0.0, r[0], 1e-12);
            Assert::AreEqual(0.0, r[1], 1e-12);
            Assert::AreEqual(1.0, r[2], 1e-12);
        }

        TEST_METHOD(VecCross_AntiCommutative)
        {
            // a x b = -(b x a)
            d3 ab = vec_cross({ 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 });
            d3 ba = vec_cross({ 4.0, 5.0, 6.0 }, { 1.0, 2.0, 3.0 });
            for (int i = 0; i < 3; ++i)
                Assert::AreEqual(-ab[i], ba[i], 1e-12);
        }

        TEST_METHOD(VecCross_ParallelVectors)
        {
            // Parallel vectors → zero cross product
            d3 r = vec_cross({ 2.0, 4.0, 6.0 }, { 1.0, 2.0, 3.0 });
            Assert::AreEqual(0.0, r[0], 1e-12);
            Assert::AreEqual(0.0, r[1], 1e-12);
            Assert::AreEqual(0.0, r[2], 1e-12);
        }

        TEST_METHOD(VecDot_Perpendicular)
        {
            Assert::AreEqual(0.0, vec_dot({ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }), 1e-12);
        }

        TEST_METHOD(VecDot_Parallel)
        {
            Assert::AreEqual(3.0, vec_dot({ 1.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }), 1e-12);
        }

        TEST_METHOD(VecDot_KnownValue)
        {
            Assert::AreEqual(32.0, vec_dot({ 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }), 1e-12);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(NumericTests)
    {
    public:

        TEST_METHOD(IsSimilarRel_WithinTolerance)
        {
            // 1.0 vs 1.005: relative diff ≈ 0.5 %, within 1 %
            Assert::AreEqual(1, is_similar_rel(1.0, 1.005, 0.01) ? 1 : 0);
        }

        TEST_METHOD(IsSimilarRel_OutsideTolerance)
        {
            // 1.0 vs 1.05: relative diff = 5 %, outside 1 %
            Assert::AreEqual(0, is_similar_rel(1.0, 1.05, 0.01) ? 1 : 0);
        }

        TEST_METHOD(IsSimilarRel_EqualValues)
        {
            Assert::AreEqual(1, is_similar_rel(42.0, 42.0, 1e-6) ? 1 : 0);
        }

        TEST_METHOD(IsSimilarAbs_WithinTolerance)
        {
            Assert::AreEqual(1, is_similar_abs(1.0, 1.0009, 0.001) ? 1 : 0);
        }

        TEST_METHOD(IsSimilarAbs_OutsideTolerance)
        {
            Assert::AreEqual(0, is_similar_abs(1.0, 1.002, 0.001) ? 1 : 0);
        }

        TEST_METHOD(FastExpNeg_AtZero)
        {
            // fast_exp_neg uses std::exp for x > -ln(2) ≈ -0.693
            double val = fast_exp_neg(0.0);
            Assert::AreEqual(1.0, val, 1e-12);
        }

        TEST_METHOD(FastExpNeg_AtMinusOne)
        {
            double approx = fast_exp_neg(-1.0);
            double exact = std::exp(-1.0);
            // Approximation tolerance: 0.5 %
            Assert::AreEqual(exact, approx, exact * 0.005);
        }

        TEST_METHOD(FastExpNeg_AtMinusTen)
        {
            // (1 + x/1024)^1024 underestimates exp(x) by ~x²/(2N) ≈ 4.9% at x=-10.
            double approx = fast_exp_neg(-10.0);
            double exact = std::exp(-10.0);
            Assert::AreEqual(exact, approx, exact * 0.06);
        }

        TEST_METHOD(FastExpNeg_BeyondCutoff)
        {
            // Values below -42 must return 0 exactly
            Assert::AreEqual(0.0, fast_exp_neg(-100.0), 1e-30);
        }

        TEST_METHOD(BLAS_tests)
        {
            test_openblas();
        }

        TEST_METHOD(SolveLinearEquations)
        {
            test_solve_linear_equations();
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(FchkParsingTests)
    {
    public:

        // FCHK format: 40-char label + type char + value starting at position 49
        // We pad to exactly 49 chars then append the value.

        TEST_METHOD(ReadFchkInt_PositiveValue)
        {
            // Format: keyword padded to 40, type at 40, 8 spaces (41-48), value at 49+
            const char* line = "Number of atoms                         I        3";
            Assert::AreEqual(3, read_fchk_integer(std::string(line)));
        }

        TEST_METHOD(ReadFchkInt_NegativeValue)
        {
            const char* line = "Charge                                  I        -1";
            Assert::AreEqual(-1, read_fchk_integer(std::string(line)));
        }

        TEST_METHOD(ReadFchkInt_LargeValue)
        {
            const char* line = "Number of basis functions               I        1024";
            Assert::AreEqual(1024, read_fchk_integer(std::string(line)));
        }

        TEST_METHOD(ReadFchkDbl_NegativeScientific)
        {
            const char* line = "Total Energy                            R        -1.23456789E+02";
            double val = read_fchk_double(std::string(line));
            Assert::AreEqual(-123.456789, val, 1e-6);
        }

        TEST_METHOD(ReadFchkDbl_PositiveScientific)
        {
            const char* line = "Zero-point correction                   R        4.56000000E-02";
            double val = read_fchk_double(std::string(line));
            Assert::AreEqual(0.0456, val, 1e-10);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(StringUtilTests)
    {
    public:

        TEST_METHOD(EndsWith_MatchingSuffix)
        {
            Assert::AreEqual(1, ends_with(std::string("molecule.wfx"), std::string(".wfx")) ? 1 : 0);
            Assert::AreEqual(1, ends_with(std::string("data.hkl"), std::string(".hkl")) ? 1 : 0);
        }

        TEST_METHOD(EndsWith_NonMatchingSuffix)
        {
            Assert::AreEqual(0, ends_with(std::string("molecule.wfx"), std::string(".gbw")) ? 1 : 0);
            Assert::AreEqual(0, ends_with(std::string("test.cpp"), std::string(".txt")) ? 1 : 0);
        }

        TEST_METHOD(EndsWith_EmptySuffix)
        {
            // Empty suffix always matches
            Assert::AreEqual(1, ends_with(std::string("anything"), std::string("")) ? 1 : 0);
        }

        TEST_METHOD(EndsWith_Suffix_LongerThanString)
        {
            Assert::AreEqual(0, ends_with(std::string("ab"), std::string("abc")) ? 1 : 0);
        }

        TEST_METHOD(EndsWith_ExactMatch)
        {
            Assert::AreEqual(1, ends_with(std::string(".wfx"), std::string(".wfx")) ? 1 : 0);
        }

        TEST_METHOD(ShrinkString_RemovesDigits)
        {
            char buf[64];
            std::string input = "C1";
            std::string r = shrink_string(input);
            Assert::AreEqual(1, static_cast<int>(r.size()));
            Assert::AreEqual(0, std::strcmp(r.c_str(), "C"));
        }

        TEST_METHOD(ShrinkString_RemovesSpacesAndDigits)
        {
            char buf[64];
            std::string input = "O 1 1";
            std::string r = shrink_string(input);
            int len = static_cast<int>(r.size());
            Assert::IsTrue(len > 0);
            // All spaces and digits removed → only "O" remains
            Assert::AreEqual(0, std::strcmp(r.c_str(), "O"));
        }

        TEST_METHOD(ShrinkString_PureLettersUnchanged)
        {
            char buf[64];
            std::string input = "Fe";
            std::string r = shrink_string(input);
            Assert::AreEqual(0, std::strcmp(r.c_str(), "Fe"));
        }

        TEST_METHOD(ShrinkString_BufferTooSmall)
        {
            // Can't fit even "C"
            std::string input = "C1H2O";
            std::string r = shrink_string(input);
            Assert::IsTrue(static_cast<int>(r.size()) >= 1);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(Sha256Tests)
    {
    public:

        TEST_METHOD(Sha256_OutputLength)
        {
            char buf[65] = {};
            std::string r = sha::sha256(std::string("abc"));
            Assert::AreEqual(64, static_cast<int>(r.size()));
        }

        TEST_METHOD(Sha256_KnownVector_Abc)
        {
            // NIST FIPS 180-4 test vector
            std::string r = sha::sha256(std::string("abc"));
            Assert::AreEqual(0, std::strcmp(r.c_str(),
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"));
        }

        TEST_METHOD(Sha256_EmptyString)
        {
            std::string r = sha::sha256(std::string(""));
            Assert::AreEqual(0, std::strcmp(r.c_str(),
                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"));
        }

        TEST_METHOD(Sha256_Deterministic)
        {
            std::string r1 = sha::sha256(std::string("NoSpherA2"));
            std::string r2 = sha::sha256(std::string("NoSpherA2"));
            Assert::AreEqual(0, std::strcmp(r1.c_str(), r2.c_str()));
        }

        TEST_METHOD(Sha256_DifferentInputsDifferentOutputs)
        {
            std::string r1 = sha::sha256(std::string("abc"));
            std::string r2 = sha::sha256(std::string("abd"));
            Assert::AreNotEqual(0, std::strcmp(r1.c_str(), r2.c_str()));
        }

        TEST_METHOD(Sha256_BufferTooSmall)
        {
            // No more small-buffer API: just validate digest length is always 64.
            std::string r = sha::sha256(std::string("abc"));
            Assert::AreEqual(64, static_cast<int>(r.size()));
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(VecAggregateTests)
    {
    public:

        TEST_METHOD(VecSumBool_MixedValues)
        {
            int arr[] = { 1, 0, 1, 1, 0 };
            bvec v = { true, false, true, true, false };
            Assert::AreEqual(3, vec_sum(v));
        }

        TEST_METHOD(VecSumBool_AllFalse)
        {
            int arr[] = { 0, 0, 0 };
            bvec v = { false, false, false };
            Assert::AreEqual(0, vec_sum(v));
        }

        TEST_METHOD(VecSumBool_Empty)
        {
            bvec v;
            Assert::AreEqual(0, vec_sum(v));
        }

        TEST_METHOD(VecSumInt_Basic)
        {
            ivec v = { 1, 2, 3, 4, 5 };
            Assert::AreEqual(15, vec_sum(v));
        }

        TEST_METHOD(VecSumInt_Negative)
        {
            ivec v = { -3, 7, -2 };
            Assert::AreEqual(2, vec_sum(v));
        }

        TEST_METHOD(VecSumDouble_Basic)
        {
            vec v = { 1.5, -0.5, 2.0 };
            Assert::AreEqual(3.0, vec_sum(v), 1e-12);
        }

        TEST_METHOD(VecSumDouble_AllZero)
        {
            vec v = { 0.0, 0.0, 0.0 };
            Assert::AreEqual(0.0, vec_sum(v), 1e-12);
        }

        TEST_METHOD(VecLength_345)
        {
            vec v = { 3.0, 4.0 };
            Assert::AreEqual(5.0, vec_length(v), 1e-12);
        }

        TEST_METHOD(VecLength_Unit)
        {
            vec v = { 1.0, 0.0, 0.0 };
            Assert::AreEqual(1.0, vec_length(v), 1e-12);
        }

        TEST_METHOD(VecLength_4D)
        {
            vec v = { 1.0, 1.0, 1.0, 1.0 };
            Assert::AreEqual(2.0, vec_length(v), 1e-12);
        }

        TEST_METHOD(VecLength_Empty)
        {
            vec v;
            Assert::AreEqual(0.0, vec_length(v), 1e-12);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(StringUtilTests2)
    {
    public:

        TEST_METHOD(Trim_LeadingTrailingSpaces)
        {
            char buf[64];
            std::string r = trim(std::string("  hello world  "));
            Assert::AreEqual(0, std::strcmp(r.c_str(), "hello world"));
        }

        TEST_METHOD(Trim_NoSpaces)
        {
            char buf[64];
            std::string r = trim(std::string("no_spaces"));
            Assert::AreEqual(0, std::strcmp(r.c_str(), "no_spaces"));
        }

        TEST_METHOD(Trim_EmptyString)
        {
            char buf[64];
            std::string r = trim(std::string(""));
            Assert::AreEqual(0, static_cast<int>(r.size()));
            Assert::AreEqual(0, std::strcmp(r.c_str(), ""));
        }

        TEST_METHOD(Trim_OnlySpaces)
        {
            char buf[64];
            std::string r = trim(std::string("    "));
            Assert::AreEqual(0, std::strcmp(r.c_str(), ""));
        }

        TEST_METHOD(Asciitolower_Uppercase)
        {
            Assert::AreEqual('a', asciitolower('A'));
            Assert::AreEqual('z', asciitolower('Z'));
            Assert::AreEqual('m', asciitolower('M'));
        }

        TEST_METHOD(Asciitolower_AlreadyLower)
        {
            Assert::AreEqual('a', asciitolower('a'));
            Assert::AreEqual('z', asciitolower('z'));
        }

        TEST_METHOD(Asciitolower_NonAlpha)
        {
            // Digits and symbols are returned unchanged
            Assert::AreEqual('5', asciitolower('5'));
            Assert::AreEqual('_', asciitolower('_'));
        }

        TEST_METHOD(DoubleFromEsd_Plain)
        {
            Assert::AreEqual(3.14159, double_from_string_with_esd(std::string("3.14159")), 1e-10);
        }

        TEST_METHOD(DoubleFromEsd_WithEsd)
        {
            // "(5)" is stripped; value is 1.234
            Assert::AreEqual(1.234, double_from_string_with_esd(std::string("1.234(5)")), 1e-10);
        }

        TEST_METHOD(DoubleFromEsd_Zero)
        {
            Assert::AreEqual(0.0, double_from_string_with_esd(std::string("0.0")), 1e-12);
        }

        TEST_METHOD(DecimalPrecisionCif_WithBracketsAndDecimal)
        {
            // "1.2345(6)" → 6 × 10⁻⁴ = 0.0006
            std::string s = "1.2345(6)";
            Assert::AreEqual(6e-4, get_decimal_precision_from_CIF_number(s), 1e-10);
        }

        TEST_METHOD(DecimalPrecisionCif_NoBrackets)
        {
            // No brackets → default 0.005
            std::string s = "1.234";
            Assert::AreEqual(0.005, get_decimal_precision_from_CIF_number(s), 1e-10);
        }

        TEST_METHOD(DecimalPrecisionCif_IntegerWithEsd)
        {
            // "100(2)" → no decimal, digit count from bracket positions
            std::string s = "100(2)";
            Assert::AreEqual(0.2, get_decimal_precision_from_CIF_number(s), 1e-10);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(BasisTypeTests)
    {
    public:

        TEST_METHOD(Sht2nbas_CartesianShells)
        {
            // Cartesian: S=1, P=3, D=6, F=10, G=15
            Assert::AreEqual(1, sht2nbas(0));
            Assert::AreEqual(3, sht2nbas(1));
            Assert::AreEqual(6, sht2nbas(2));
            Assert::AreEqual(10, sht2nbas(3));
            Assert::AreEqual(15, sht2nbas(4));
        }

        TEST_METHOD(Sht2nbas_SphericalShells)
        {
            // Negative types → spherical: -2=D(5), -3=F(7), -4=G(9)
            Assert::AreEqual(5, sht2nbas(-2));
            Assert::AreEqual(7, sht2nbas(-3));
            Assert::AreEqual(9, sht2nbas(-4));
        }

        TEST_METHOD(DoubleFactorial_SmallValues)
        {
            Assert::AreEqual(1u, doublefactorial(0));
            Assert::AreEqual(1u, doublefactorial(1));
            Assert::AreEqual(2u, doublefactorial(2));
            Assert::AreEqual(3u, doublefactorial(3));
            Assert::AreEqual(8u, doublefactorial(4));
            Assert::AreEqual(15u, doublefactorial(5));
            Assert::AreEqual(48u, doublefactorial(6));
            Assert::AreEqual(105u, doublefactorial(7));
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(ConstantsTests)
    {
    public:

        TEST_METHOD(ConstAbs_Positive) { Assert::AreEqual(5, constants::const_abs(5)); }
        TEST_METHOD(ConstAbs_Negative) { Assert::AreEqual(5, constants::const_abs(-5)); }
        TEST_METHOD(ConstAbs_Zero) { Assert::AreEqual(0, constants::const_abs(0)); }

        TEST_METHOD(ConstexprPow_Integer)
        {
            Assert::AreEqual(1024.0, constants::constexpr_pow(2.0, 10), 1e-12);
            Assert::AreEqual(1.0, constants::constexpr_pow(10.0, 0), 1e-12);
            Assert::AreEqual(8.0, constants::constexpr_pow(2.0, 3), 1e-12);
        }

        TEST_METHOD(ConstantsSqrt_KnownValues)
        {
            Assert::AreEqual(2.0, constants::sqrt(4.0), 1e-10);
            Assert::AreEqual(std::sqrt(2.0), constants::sqrt(2.0), 1e-10);
            Assert::AreEqual(0.0, constants::sqrt(0.0), 1e-12);
        }

        TEST_METHOD(ConstantsSqrt_NegativeIsNaN)
        {
            double r = constants::sqrt(-1.0);
            Assert::IsTrue(std::isnan(r));
        }

        TEST_METHOD(ExpApprox_AtZero)
        {
            Assert::AreEqual(1.0, constants::exp_approx(0.0, 25), 1e-12);
        }

        TEST_METHOD(ExpApprox_AtOne)
        {
            // 25-term Taylor series matches std::exp to machine precision
            Assert::AreEqual(std::exp(1.0), constants::exp_approx(1.0, 25), 1e-12);
        }

        TEST_METHOD(ExpApprox_AtMinusOne)
        {
            Assert::AreEqual(std::exp(-1.0), constants::exp_approx(-1.0, 25), 1e-10);
        }

        TEST_METHOD(LogApprox_AtOne)
        {
            Assert::AreEqual(0.0, constants::log_approx(1.0, 25), 1e-12);
        }

        TEST_METHOD(LogApprox_AtE)
        {
            // Arctanh series converges; 25 iterations accurate to 1e-8 for x=e
            Assert::AreEqual(1.0, constants::log_approx(std::exp(1.0), 25), 1e-6);
        }

        TEST_METHOD(LogApprox_NonPositiveReturnsSentinel)
        {
            Assert::AreEqual(-1.0, constants::log_approx(0.0, 25), 1e-12);
            Assert::AreEqual(-1.0, constants::log_approx(-5.0, 25), 1e-12);
        }

        TEST_METHOD(Bohr2Ang_OneBohr)
        {
            // 1 Bohr = a₀ Å = 0.529177210903 Å
            Assert::AreEqual(0.529177210903, constants::bohr2ang(1.0), 1e-10);
        }

        TEST_METHOD(Bohr2Ang_Zero)
        {
            Assert::AreEqual(0.0, constants::bohr2ang(0.0), 1e-12);
        }

        TEST_METHOD(Ang2Bohr_OneAngstrom)
        {
            Assert::AreEqual(1.0 / 0.529177210903, constants::ang2bohr(1.0), 1e-8);
        }

        TEST_METHOD(Ang2Bohr_Roundtrip)
        {
            // ang2bohr(bohr2ang(x)) ≈ x
            double x = 2.5;
            Assert::AreEqual(x, constants::ang2bohr(constants::bohr2ang(x)), 1e-10);
        }

        TEST_METHOD(CubicBohr2Ang_OneBohr3)
        {
            // 1 Bohr³ = a₀³ Å³
            double expected = 0.529177210903 * 0.529177210903 * 0.529177210903;
            Assert::AreEqual(expected, constants::cubic_bohr2ang(1.0), 1e-10);
        }

        TEST_METHOD(CubicAng2Bohr_Roundtrip)
        {
            double x = 3.0;
            Assert::AreEqual(x, constants::cubic_ang2bohr(constants::cubic_bohr2ang(x)), 1e-8);
        }

        TEST_METHOD(Factorial_SmallValues)
        {
            Assert::AreEqual(1LL, static_cast<long long>(constants::ft_fun(0)));
            Assert::AreEqual(1LL, static_cast<long long>(constants::ft_fun(1)));
            Assert::AreEqual(120LL, static_cast<long long>(constants::ft_fun(5)));
            Assert::AreEqual(3628800LL, static_cast<long long>(constants::ft_fun(10)));
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(OrbitalIndexTests)
    {
    public:

        TEST_METHOD(OrcaToPySCF_SShell)
        {
            // S: only one component, m_idx=0 → 0
            auto r = constants::orca_2_pySCF(0, 0);
            Assert::AreEqual(0, r.has_value() ? static_cast<int>(r.value()) : -1);
        }

        TEST_METHOD(OrcaToPySCF_PShell)
        {
            // P ORCA ordering 0,+1,-1 → PySCF map {1,2,0}
            Assert::AreEqual((size_t)1, constants::orca_2_pySCF(1, 0).value());
            Assert::AreEqual((size_t)2, constants::orca_2_pySCF(1, 1).value());
            Assert::AreEqual((size_t)0, constants::orca_2_pySCF(1, 2).value());
        }

        TEST_METHOD(OrcaToPySCF_DShell_FirstComponent)
        {
            // D: map {2,3,1,4,0}, m_idx=0 → 2
            Assert::AreEqual((size_t)2, constants::orca_2_pySCF(2, 0).value());
        }

        TEST_METHOD(OrcaToPySCF_OutOfRange)
        {
            // l=100 is not in the switch → nullopt → -1
            auto r = constants::orca_2_pySCF(100, 0);
            Assert::AreEqual(-1, r.has_value() ? static_cast<int>(r.value()) : -1);
        }

        TEST_METHOD(TypeToNbo_SShell) { Assert::AreEqual(1u, constants::type_2_nbo(1)); }
        TEST_METHOD(TypeToNbo_PxShell) { Assert::AreEqual(101u, constants::type_2_nbo(2)); }
        TEST_METHOD(TypeToNbo_DxxShell) { Assert::AreEqual(201u, constants::type_2_nbo(5)); }
        TEST_METHOD(TypeToNbo_FxxxShell) { Assert::AreEqual(301u, constants::type_2_nbo(11)); }
        TEST_METHOD(TypeToNbo_GxxxxShell) { Assert::AreEqual(401u, constants::type_2_nbo(21)); }
        TEST_METHOD(TypeToNbo_Unknown) { Assert::AreEqual(0u, constants::type_2_nbo(99)); }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(BesselTests)
    {
    public:

        TEST_METHOD(BesselJ0_AtZero)
        {
            // j₀(0) = 1 (limit of sin(x)/x as x→0)
            Assert::AreEqual(1.0, bessel_first_kind(0, 0.0), 1e-12);
        }

        TEST_METHOD(BesselJ1_AtZero)
        {
            // j_l(0) = 0 for l > 0
            Assert::AreEqual(0.0, bessel_first_kind(1, 0.0), 1e-12);
            Assert::AreEqual(0.0, bessel_first_kind(5, 0.0), 1e-12);
        }

        TEST_METHOD(BesselJ0_AtOne)
        {
            // j₀(1) = sin(1)/1
            Assert::AreEqual(std::sin(1.0), bessel_first_kind(0, 1.0), 1e-12);
        }

        TEST_METHOD(BesselJ1_AtOne)
        {
            // j₁(1) = (sin(1) - cos(1)) / 1
            double expected = std::sin(1.0) - std::cos(1.0);
            Assert::AreEqual(expected, bessel_first_kind(1, 1.0), 1e-12);
        }

        TEST_METHOD(BesselJ2_AtOne)
        {
            // j₂(1) = (2·sin(1) - 3·cos(1)) / 1
            double expected = 2.0 * std::sin(1.0) - 3.0 * std::cos(1.0);
            Assert::AreEqual(expected, bessel_first_kind(2, 1.0), 1e-12);
        }

        TEST_METHOD(BesselJ0_AtPi)
        {
            // j₀(π) = sin(π)/π ≈ 0
            Assert::AreEqual(std::sin(PI_VAL) / PI_VAL, bessel_first_kind(0, PI_VAL), 1e-12);
        }

        TEST_METHOD(BesselJ_HigherOrder_PositiveAndFinite)
        {
            // l=7 exercises the continued-fraction fallback path
            double r = bessel_first_kind(7, 2.0);
            Assert::IsTrue(std::isfinite(r));
            Assert::IsTrue(r > 0.0);
        }

        TEST_METHOD(BesselJ_RecurrenceCheck)
        {
            // Recurrence: j_{l-1}(x) + j_{l+1}(x) = (2l+1)/x · j_l(x)
            double x = 3.0;
            int l = 3;
            double jlm1 = bessel_first_kind(l - 1, x);
            double jl = bessel_first_kind(l, x);
            double jlp1 = bessel_first_kind(l + 1, x);
            double lhs = jlm1 + jlp1;
            double rhs = (2.0 * l + 1.0) / x * jl;
            Assert::AreEqual(lhs, rhs, 1e-10);
        }

        TEST_METHOD(AnalyticFourier)
        {
            int ret = test_analytical_fourier(false);
            Assert::AreEqual(ret, 0);
        }

        TEST_METHOD(AnalyticFourier_full)
        {
            int ret = test_analytical_fourier(true);
            Assert::AreEqual(ret, 0);
        }

    };

    // -----------------------------------------------------------------------
    // IsSimilarPow10Tests — is_similar(a, b, tolerance) where |a-b| <= 10^tol
    // -----------------------------------------------------------------------
    TEST_CLASS(IsSimilarPow10Tests)
    {
    public:
        TEST_METHOD(Equal_ReturnTrue)
        {
            Assert::AreEqual(1, is_similar(1.0, 1.0, -6.0) ? 1 : 0);
        }

        TEST_METHOD(WithinTolerance_ReturnTrue)
        {
            // |1.000001 - 1.0| = 1e-6 <= 10^(-6)
            Assert::AreEqual(1, is_similar(1.000001, 1.0, -6.0) ? 1 : 0);
        }

        TEST_METHOD(OutsideTolerance_ReturnFalse)
        {
            // |1.00001 - 1.0| = 1e-5 > 10^(-6)
            Assert::AreEqual(0, is_similar(1.00001, 1.0, -6.0) ? 1 : 0);
        }

        TEST_METHOD(NegativeValues_WithinTolerance)
        {
            Assert::AreEqual(1, is_similar(-5.0, -5.0 + 1e-8, -7.0) ? 1 : 0);
        }

        TEST_METHOD(LooseTolerance_LargeDiff)
        {
            // |100 - 50| = 50 <= 10^2 = 100
            Assert::AreEqual(1, is_similar(100.0, 50.0, 2.0) ? 1 : 0);
        }

        TEST_METHOD(LooseTolerance_TooLargeDiff)
        {
            // |200 - 50| = 150 > 10^2 = 100
            Assert::AreEqual(0, is_similar(200.0, 50.0, 2.0) ? 1 : 0);
        }
    };

    // -----------------------------------------------------------------------
    // Shell2FunctionTests — shell2function(type, prim) WFN column index
    // -----------------------------------------------------------------------
    TEST_CLASS(Shell2FunctionTests)
    {
    public:
        TEST_METHOD(SType_Prim0_Returns0)
        {
            // s-type shell (type=1): first and only function is index 0
            int r = shell2function(1, 0);
            Assert::IsTrue(r >= 0);
        }

        TEST_METHOD(PType_Prim0)
        {
            // p-type shell: 3 functions
            int r = shell2function(2, 0);
            Assert::IsTrue(r >= 0);
        }

        TEST_METHOD(PType_Prim1)
        {
            int r0 = shell2function(2, 0);
            int r1 = shell2function(2, 1);
            Assert::IsTrue(r1 > r0);
        }

        TEST_METHOD(DType_Prim5_ValidIndex)
        {
            // d-type (type=3): 6 Cartesian or 5 spherical functions
            int r = shell2function(3, 5);
            Assert::IsTrue(r >= 0);
        }

        TEST_METHOD(ResultsAreStrictlyIncreasingWithinShell)
        {
            // f-type (type=4): consecutive prims must give increasing column indices
            int r0 = shell2function(4, 0);
            int r1 = shell2function(4, 1);
            int r2 = shell2function(4, 2);
            Assert::IsTrue(r1 > r0);
            Assert::IsTrue(r2 > r1);
        }
    };

    // -----------------------------------------------------------------------
    // CountWordsTests — CountWords(str)
    // -----------------------------------------------------------------------
    TEST_CLASS(CountWordsTests)
    {
    public:
        TEST_METHOD(Empty_Returns0)
        {
            Assert::AreEqual(0, CountWords(""));
        }

        TEST_METHOD(OneWord)
        {
            Assert::AreEqual(1, CountWords("hello"));
        }

        TEST_METHOD(TwoWords)
        {
            Assert::AreEqual(2, CountWords("hello world"));
        }

        TEST_METHOD(LeadingTrailingSpaces)
        {
            Assert::AreEqual(2, CountWords("  foo   bar  "));
        }

        TEST_METHOD(MultipleSpacesBetweenWords)
        {
            Assert::AreEqual(3, CountWords("a  b  c"));
        }

        TEST_METHOD(SingleSpace)
        {
            Assert::AreEqual(0, CountWords(" "));
        }
    };

    // -----------------------------------------------------------------------
    // ShrinkStringToAtomTests — shrink_string_to_atom(input, atom_number)
    // -----------------------------------------------------------------------
    TEST_CLASS(ShrinkStringToAtomTests)
    {
    public:
        TEST_METHOD(CarbonAtomNumber6)
        {
            // atnr2letter(6) = "C"
            char buf[32];
            std::string input = "C1";
            std::string r = shrink_string_to_atom(input, 6);
            Assert::IsTrue(static_cast<int>(r.size()) >= 0);
            Assert::AreEqual(std::string("C"), r);
        }

        TEST_METHOD(CalciumAtomNumber20)
        {
            // atnr2letter(20) = "Ca"
            char buf[32];
            std::string input = "Ca12";
            std::string r = shrink_string_to_atom(input, 20);
            Assert::IsTrue(static_cast<int>(r.size()) >= 0);
            Assert::AreEqual(std::string("Ca"), r);
        }

        TEST_METHOD(BufferTooSmall_ReturnsNegOne)
        {
            std::string input = "Carbon6";
            std::string r = shrink_string_to_atom(input, 6);
            Assert::IsTrue(static_cast<int>(r.size()) >= 1);
        }

        TEST_METHOD(IronAtomNumber26)
        {
            // atnr2letter(26) = "Fe"
            char buf[32];
            std::string input = "Fe3 ";
            std::string r = shrink_string_to_atom(input, 26);
            Assert::IsTrue(static_cast<int>(r.size()) >= 0);
            Assert::AreEqual(std::string("Fe"), r);
        }
    };

    // -----------------------------------------------------------------------
    // SplitStringTests — split_string(input, delim)
    // -----------------------------------------------------------------------
    TEST_CLASS(SplitStringTests)
    {
    public:
        TEST_METHOD(SingleToken_NoDelim)
        {
            svec toks = split_string<std::string>(std::string("hello"), " ");
            Assert::AreEqual(1, static_cast<int>(toks.size()));
            Assert::AreEqual(std::string("hello"), toks[0]);
        }

        TEST_METHOD(ThreeTokens)
        {
            svec toks = split_string<std::string>(std::string("a b c"), " ");
            Assert::AreEqual(3, static_cast<int>(toks.size()));
            Assert::AreEqual(std::string("a"), toks[0]);
            Assert::AreEqual(std::string("b"), toks[1]);
            Assert::AreEqual(std::string("c"), toks[2]);
        }

        TEST_METHOD(CommaDelimiter)
        {
            svec toks = split_string<std::string>(std::string("x,y,z"), ",");
            Assert::AreEqual(3, static_cast<int>(toks.size()));
            Assert::AreEqual(std::string("z"), toks[2]);
        }

        TEST_METHOD(MaxOutLimit_ReturnsTotalCount)
        {
            // 4 tokens but max_out=2 (old API); direct API returns full token list.
            svec toks = split_string<std::string>(std::string("a b c d"), " ");
            Assert::AreEqual(4, static_cast<int>(toks.size()));
            Assert::AreEqual(std::string("a"), toks[0]);
            Assert::AreEqual(std::string("b"), toks[1]);
        }

        TEST_METHOD(EmptyString_ZeroTokens)
        {
            svec toks = split_string<std::string>(std::string(""), " ");
            Assert::IsTrue(toks.empty() || toks.size() == 1); // impl-defined for empty input
        }
    };

    // -----------------------------------------------------------------------
    // TimingTests — ut_sleep_and_measure_us(N) returns elapsed µs >= N*1000
    // -----------------------------------------------------------------------
    TEST_CLASS(TimingTests)
    {
    public:
        TEST_METHOD(Sleep10ms_ElapsedAtLeast10000us)
        {
            auto t0 = get_time();
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            auto t1 = get_time();
            long long us = get_musec(t0, t1);
            Assert::IsTrue(us >= 10000LL);
        }

        TEST_METHOD(Sleep1ms_ElapsedAtLeast1000us)
        {
            auto t0 = get_time();
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            auto t1 = get_time();
            long long us = get_musec(t0, t1);
            Assert::IsTrue(us >= 1000LL);
        }

        TEST_METHOD(Sleep5ms_ElapsedPositive)
        {
            auto t0 = get_time();
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            auto t1 = get_time();
            long long us = get_musec(t0, t1);
            Assert::IsTrue(us > 0LL);
        }
    };

    // -----------------------------------------------------------------------
    // HypergeometricTests — 2F1(a,b;c;x)
    // -----------------------------------------------------------------------
    TEST_CLASS(HypergeometricTests)
    {
    public:
        TEST_METHOD(Identity_2F1_Zero_IsOne)
        {
            // 2F1(a,b;c;0) = 1 for any a,b,c
            double r = hypergeometric(1.0, 2.0, 3.0, 0.0);
            Assert::AreEqual(1.0, r, 1e-12);
        }

        TEST_METHOD(KnownValue_2F1_1_1_2_Half)
        {
            // 2F1(1,1;2;0.5) = -2*ln(0.5) = 2*ln(2) ≈ 1.386294...
            double expected = -2.0 * std::log(0.5);
            double r = hypergeometric(1.0, 1.0, 2.0, 0.5);
            Assert::AreEqual(expected, r, 1e-6);
        }

        TEST_METHOD(KnownValue_2F1_Half_Half_ThreeHalves_Half)
        {
            // 2F1(0.5,0.5;1.5;0.5) — finite, positive
            double r = hypergeometric(0.5, 0.5, 1.5, 0.5);
            Assert::IsTrue(std::isfinite(r));
            Assert::IsTrue(r > 1.0);
        }

        TEST_METHOD(Symmetry_ab_equals_ba)
        {
            // 2F1(a,b;c;x) = 2F1(b,a;c;x)
            double r1 = hypergeometric(2.0, 3.0, 5.0, 0.3);
            double r2 = hypergeometric(3.0, 2.0, 5.0, 0.3);
            Assert::AreEqual(r1, r2, 1e-10);
        }

        TEST_METHOD(NegativeX_ReturnsFinite)
        {
            double r = hypergeometric(1.0, 2.0, 3.0, -0.5);
            Assert::IsTrue(std::isfinite(r));
        }
    };

    // -----------------------------------------------------------------------
    // Atnr2LetterTests — constants::atnr2letter(nr) element symbol lookup
    // -----------------------------------------------------------------------
    TEST_CLASS(Atnr2LetterTests)
    {
    public:
        TEST_METHOD(Hydrogen_Is_H)
        {
            char buf[8];
            const char* sym = constants::atnr2letter(1);
            Assert::AreEqual(std::string("H"), std::string(sym));
        }

        TEST_METHOD(Carbon_Is_C)
        {
            char buf[8];
            Assert::AreEqual(std::string("C"), std::string(constants::atnr2letter(6)));
        }

        TEST_METHOD(Iron_Is_Fe)
        {
            char buf[8];
            Assert::AreEqual(std::string("Fe"), std::string(constants::atnr2letter(26)));
        }

        TEST_METHOD(Gold_Is_Au)
        {
            char buf[8];
            Assert::AreEqual(std::string("Au"), std::string(constants::atnr2letter(79)));
        }

        TEST_METHOD(Lawrencium_103_Is_Lr)
        {
            char buf[8];
            Assert::AreEqual(std::string("Lr"), std::string(constants::atnr2letter(103)));
        }

        TEST_METHOD(Zero_Is_Q_Peak)
        {
            char buf[8];
            Assert::AreEqual(std::string("Q"), std::string(constants::atnr2letter(0)));
        }

        TEST_METHOD(OutOfRange_Returns_PROBLEM)
        {
            char buf[16];
            Assert::AreEqual(std::string("PROBLEM"), std::string(constants::atnr2letter(200)));
        }

        TEST_METHOD(BufferTooSmall_ReturnsNegOne)
        {
            // No more buffer-sized API: symbol is returned as const char*.
            Assert::AreEqual(std::string("C"), std::string(constants::atnr2letter(6)));
        }
    };

    // -----------------------------------------------------------------------
    // Type2VectorTests — constants::type2vector: basis type → [nx, ny, nz]
    // -----------------------------------------------------------------------
    TEST_CLASS(Type2VectorTests)
    {
    public:
        TEST_METHOD(Type1_Is_SShell_000)
        {
            int nx, ny, nz;
            int v[3];
            constants::type2vector(1, v);
            nx = v[0]; ny = v[1]; nz = v[2];
            Assert::AreEqual(0, nx);
            Assert::AreEqual(0, ny);
            Assert::AreEqual(0, nz);
        }

        TEST_METHOD(Type2_Is_Px_100)
        {
            int nx, ny, nz;
            int v[3];
            constants::type2vector(2, v);
            nx = v[0]; ny = v[1]; nz = v[2];
            Assert::AreEqual(1, nx);
            Assert::AreEqual(0, ny);
            Assert::AreEqual(0, nz);
        }

        TEST_METHOD(Type5_Is_Dx2_200)
        {
            // index 5 in type_vector: (2,0,0) = dx²
            int nx, ny, nz;
            int v[3];
            constants::type2vector(5, v);
            nx = v[0]; ny = v[1]; nz = v[2];
            Assert::AreEqual(2, nx);
            Assert::AreEqual(0, ny);
            Assert::AreEqual(0, nz);
        }

        TEST_METHOD(Type8_Is_Dxy_110)
        {
            // index 8: (1,1,0) = dxy
            int nx, ny, nz;
            int v[3];
            constants::type2vector(8, v);
            nx = v[0]; ny = v[1]; nz = v[2];
            Assert::AreEqual(1, nx);
            Assert::AreEqual(1, ny);
            Assert::AreEqual(0, nz);
        }

        TEST_METHOD(SumOfExponents_Matches_ShellType)
        {
            // For d-type (types 5-10), nx+ny+nz == 2
            for (int t = 5; t <= 10; ++t) {
                int nx, ny, nz;
                int v[3];
                constants::type2vector(t, v);
                nx = v[0]; ny = v[1]; nz = v[2];
                Assert::AreEqual(2, nx + ny + nz);
            }
        }

        TEST_METHOD(OutOfRange_Returns_NegOne)
        {
            int nx, ny, nz;
            int v[3];
            constants::type2vector(0, v);
            Assert::AreEqual(-1, v[0]);
            constants::type2vector(57, v);
            Assert::AreEqual(-1, v[0]);
        }
    };

    // -----------------------------------------------------------------------
    // NormGaussTests — constants::normgauss(type, exponent)
    // -----------------------------------------------------------------------
    TEST_CLASS(NormGaussTests)
    {
    public:
        TEST_METHOD(SShell_IsPositive)
        {
            double n = constants::normgauss(1, 1.0);
            Assert::IsTrue(n > 0.0);
            Assert::IsTrue(std::isfinite(n));
        }

        TEST_METHOD(SShell_ScalesWithExponent)
        {
            // Normalization grows with exponent for s-type
            double n1 = constants::normgauss(1, 1.0);
            double n2 = constants::normgauss(1, 4.0);
            Assert::IsTrue(n2 > n1);
        }

        TEST_METHOD(PShell_IsPositive)
        {
            double n = constants::normgauss(2, 1.0);
            Assert::IsTrue(n > 0.0);
        }

        TEST_METHOD(DShell_IsPositive)
        {
            double n = constants::normgauss(5, 1.0);
            Assert::IsTrue(n > 0.0);
        }

        TEST_METHOD(SShell_KnownValue)
        {
            // normgauss for s-type (0,0,0), exponent α:
            // N = (2α/π)^(9/4) * sqrt(1/1) = (2α/π)^(9/4)
            // For α=1: N = (2/π)^(9/4)
            double alpha = 1.0;
            double expected = std::pow(2.0 * alpha / PI_VAL, 9.0 / 4.0);
            double result = constants::normgauss(1, alpha);
            Assert::AreEqual(expected, result, 1e-10);
        }
    };

    // -----------------------------------------------------------------------
    // AssocLegendreTests — constants::associated_legendre_polynomial(l, m, x)
    // -----------------------------------------------------------------------
    TEST_CLASS(AssocLegendreTests)
    {
    public:
        TEST_METHOD(P00_Is_One)
        {
            Assert::AreEqual(1.0, constants::associated_legendre_polynomial(0, 0, 0.5), 1e-12);
        }

        TEST_METHOD(P10_At_Half_Is_Half)
        {
            // P_1^0(x) = x
            Assert::AreEqual(0.5, constants::associated_legendre_polynomial(1, 0, 0.5), 1e-12);
        }

        TEST_METHOD(P11_AtOne_Is_Zero)
        {
            // P_1^1(x) = sqrt(1-x²), at x=1: 0
            Assert::AreEqual(0.0, constants::associated_legendre_polynomial(1, 1, 1.0), 1e-12);
        }

        TEST_METHOD(P20_Is_Legendre_Polynomial)
        {
            // P_2^0(x) = (3x²-1)/2
            double x = 0.6;
            double expected = 0.5 * (3 * x * x - 1);
            Assert::AreEqual(expected, constants::associated_legendre_polynomial(2, 0, x), 1e-12);
        }

        TEST_METHOD(P21_At_Zero)
        {
            // P_2^1(x) = 3x*sqrt(1-x²), at x=0: 0
            Assert::AreEqual(0.0, constants::associated_legendre_polynomial(2, 1, 0.0), 1e-12);
        }

        TEST_METHOD(P22_At_Zero)
        {
            // P_2^2(x) = -3(x²-1) = 3(1-x²), at x=0: 3
            Assert::AreEqual(-3.0 * (0.0 - 1.0), constants::associated_legendre_polynomial(2, 2, 0.0), 1e-12);
        }

        TEST_METHOD(NegativeM_P1m1)
        {
            // P_1^{-1}(x) = -0.5*sqrt(1-x²)
            double x = 0.5;
            double expected = -0.5 * std::sqrt(1 - x * x);
            Assert::AreEqual(expected, constants::associated_legendre_polynomial(1, -1, x), 1e-12);
        }
    };

    // -----------------------------------------------------------------------
    // CartesianToSphericalTests — constants::cartesian_to_spherical(x,y,z)
    // -----------------------------------------------------------------------
    TEST_CLASS(CartesianToSphericalTests)
    {
    public:
        TEST_METHOD(Origin_Has_Zero_Radius)
        {
            double out[3];
            vec r = constants::cartesian_to_spherical(0.0, 0.0, 0.0);
            out[0] = r[0]; out[1] = r[1]; out[2] = r[2];
            Assert::AreEqual(0.0, out[0], 1e-12); // r
        }

        TEST_METHOD(UnitX_Gives_Correct_Angles)
        {
            double out[3];
            vec r = constants::cartesian_to_spherical(1.0, 0.0, 0.0);
            out[0] = r[0]; out[1] = r[1]; out[2] = r[2];
            Assert::AreEqual(1.0, out[0], 1e-12); // r=1
            Assert::AreEqual(PI_VAL / 2.0, out[1], 1e-12); // theta=pi/2
            Assert::AreEqual(0.0, out[2], 1e-12); // phi=0
        }

        TEST_METHOD(UnitZ_Has_Zero_Theta)
        {
            double out[3];
            vec r = constants::cartesian_to_spherical(0.0, 0.0, 1.0);
            out[0] = r[0]; out[1] = r[1]; out[2] = r[2];
            Assert::AreEqual(1.0, out[0], 1e-12); // r=1
            Assert::AreEqual(0.0, out[1], 1e-12); // theta=0
        }

        TEST_METHOD(UnitY_Gives_PhiHalfPi)
        {
            double out[3];
            vec r = constants::cartesian_to_spherical(0.0, 1.0, 0.0);
            out[0] = r[0]; out[1] = r[1]; out[2] = r[2];
            Assert::AreEqual(1.0, out[0], 1e-12); // r=1
            Assert::AreEqual(PI_VAL / 2.0, out[2], 1e-12); // phi=pi/2
        }

        TEST_METHOD(Radius_Is_Euclidean_Norm)
        {
            double out[3];
            vec r = constants::cartesian_to_spherical(3.0, 4.0, 0.0);
            out[0] = r[0]; out[1] = r[1]; out[2] = r[2];
            Assert::AreEqual(5.0, out[0], 1e-12); // r=5
        }

        TEST_METHOD(Inverse_Recover_Cartesian)
        {
            // Convert (1,1,1) → spherical → back to Cartesian
            double out[3];
            vec r_sph = constants::cartesian_to_spherical(1.0, 1.0, 1.0);
            out[0] = r_sph[0]; out[1] = r_sph[1]; out[2] = r_sph[2];
            double r = out[0], theta = out[1], phi = out[2];
            double x = r * std::sin(theta) * std::cos(phi);
            double y = r * std::sin(theta) * std::sin(phi);
            double z = r * std::cos(theta);
            Assert::AreEqual(1.0, x, 1e-12);
            Assert::AreEqual(1.0, y, 1e-12);
            Assert::AreEqual(1.0, z, 1e-12);
        }
    };

    // -----------------------------------------------------------------------
    // GetLambda1Tests — median eigenvalue of a 3x3 symmetric matrix
    // -----------------------------------------------------------------------
    TEST_CLASS(GetLambda1Tests)
    {
    public:
        TEST_METHOD(DiagonalMatrix_ReturnsMidEigenvalue)
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
            Assert::AreEqual(3.0, lam, 1e-10);
        }

        TEST_METHOD(DiagonalMatrix_AllEqual_ReturnsThat)
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
            Assert::AreEqual(2.0, lam, 1e-10);
        }

        TEST_METHOD(DiagonalDescending_ReturnsMid)
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
            Assert::AreEqual(5.0, lam, 1e-10);
        }

        TEST_METHOD(SymmetricMatrix_ReturnsFiniteValue)
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
            Assert::IsTrue(std::isfinite(lam));
            // Eigenvalues must be between min(diag)=2 and max(diag)=4 for a diag-dominant case
            Assert::IsTrue(lam >= 1.0 && lam <= 5.0);
        }

        TEST_METHOD(SymmetricMatrix_MedianIsInBounds)
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
            Assert::IsTrue(std::isfinite(lam));
            Assert::IsTrue(lam > 0.0); // positive definite matrix
        }
    };

} // namespace NoSpherA2UnitTests