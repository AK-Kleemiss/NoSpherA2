//scf_backend dusk at the boundary to occ: the basis and its function order (the overlap of each library in the other's
//order), and dusk's plain Hartree-Fock against occ's, as ScfSolverTests does for SCF_solver
#include "pch.h"
#include <gtest/gtest.h>

#include "core/xcw/dusk_scf.h"
#include "core/xcw/SCF_solver.h"
#include "core/wfn/basis_set.h"
#include <occ/io/xyz.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#include <filesystem>
#include <fstream>
#include <sstream>

#ifdef NOSPHERA2_USE_DUSK
namespace
{
	const char* water = "3\n\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.469\nH 0.000 -0.757 -0.469\n";
	const char* methyl = "4\n\nC 0.000 0.000 0.000\nH 1.079 0.000 0.000\nH -0.540 0.935 0.000\nH -0.540 -0.935 0.000\n";

	struct pair_run {
		double occ_energy = 0, energy = 0, warm_energy = 0;
		bool converged = false, warm_converged = false;
		int warm_iterations = 0;
		double max_dD = 0;
	};

	//occ's SCF at an energy threshold of 1e-10 and dusk's at its tight preset, then dusk again from its own orbitals
	pair_run solve(const char* xyz, const std::string& basis_name, const occ::qm::SpinorbitalKind kind, const int multiplicity) {
		spdlog::set_level(spdlog::level::err);
		occ::core::Molecule mol = occ::io::molecule_from_xyz_string(xyz);
		mol.set_multiplicity(multiplicity);
		const occ::qm::AOBasis basis = SCF_solver::setup_basis(mol, basis_name);
		occ::qm::HartreeFock hf(basis);
		occ::qm::SCF<occ::qm::HartreeFock> reference(hf, kind);
		reference.set_charge_multiplicity(0, multiplicity);
		reference.convergence_settings.energy_threshold = 1e-10;
		pair_run r;
		r.occ_energy = reference.compute_scf_energy();
		dusk::Input in = dusk_occ::input(mol, *BasisSetLibrary::get_basis_set(basis_name));
		dusk::Options& o = in.options;
		o.reference = kind == occ::qm::SpinorbitalKind::Unrestricted ? dusk::Reference::Unrestricted : dusk::Reference::Restricted;
		o.convergence = dusk::Convergence::tight();
		o.integral_threshold = 1e-14;
		o.verbosity = 0;
		dusk::SCF scf = dusk::make_scf(in);
		const dusk::Result first = scf.run();
		r.converged = first.converged;
		r.energy = first.energy;
		occ::qm::Wavefunction wfn;
		wfn.basis = basis;
		wfn.atoms = mol.atoms();
		dusk_occ::to_occ(first, wfn);
		occ::qm::MolecularOrbitals ref = reference.ctx.mo;
		ref.update_occupied_orbitals();
		ref.update_density_matrix();
		//the total density: occ puts a doublet's unpaired electron in beta (n_alpha 4, n_beta 5 for CH3), dusk in alpha
		r.max_dD = (wfn.mo.total_density_matrix() - ref.total_density_matrix()).cwiseAbs().maxCoeff();
		scf.set_guess(dusk_occ::to_dusk(wfn));
		const dusk::Result warm = scf.run();
		r.warm_converged = warm.converged;
		r.warm_energy = warm.energy;
		r.warm_iterations = warm.iterations;
		return r;
	}

	//dusk's options for a settings file read by loadXCWsettings, and what scf_options logged
	dusk::Options options_of(const std::string& text, std::string& log) {
		const auto p = std::filesystem::temp_directory_path() / ("nosphera2_dusk_settings_" + std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".txt");
		std::ofstream(p) << text;
		options opt;
		opt.xcw_settings_path = p;
		opt.do_XCW = true;
		std::stringstream out;
		std::streambuf* const old = std::cout.rdbuf(out.rdbuf());
		opt.loadXCWsettings();
		std::cout.rdbuf(old);
		std::filesystem::remove(p);
		std::ostringstream l;
		const dusk::Options o = dusk_occ::scf_options(opt.xcw_settings, l);
		log = l.str();
		return o;
	}

	//The coupling of an XCW-like criterion c = (Tr(dm M) - t)^2 for a fixed symmetric M, dm = D / 2 as XCW_wrapper gets it,
	//counting the evaluations dusk let be incremental
	struct mock_coupling : SCF_coupling {
		occ::Mat M;
		double t = 0, x = 0;
		int full = 0, incremental = 0, saves = 0, restores = 0;
		static double trace(const dMatrix2& A, const occ::Mat& B) {
			double s = 0;
			for (int i = 0; i < B.rows(); i++)
				for (int j = 0; j < B.cols(); j++) s += A(i, j) * B(i, j);
			return s;
		}
		void evaluate_cryst_crit(const dMatrix2& dm, const bool inc) override {
			(inc ? incremental : full)++;
			x = trace(dm, M) - t;
		}
		fit_quality cryst_crit() const override { return { x * x, x * x, 0, x * x, 0, 0, 0 }; }
		void perturbation_matrix(occ::Mat& P, const bool, const double, const bool) override { P = x * M; }
		occ::Mat perturbation_matrix_response(const dMatrix2& dD, const double lambda) override { return lambda * trace(dD, M) * M; }
		void save_F_calc() override { saves++; }
		void restore_F_calc() override { restores++; }
		bool grid_on_device() const override { return false; }
	};
}

//The overlap dusk computes over the converted basis equals occ's, element by element: same functions, same order and phases
TEST(DuskBackendTests, OverlapMatchesOccThroughF)
{
	spdlog::set_level(spdlog::level::err);
	const occ::core::Molecule mol = occ::io::molecule_from_xyz_string(water);
	const occ::qm::AOBasis basis = SCF_solver::setup_basis(mol, "def2-tzvp");
	ASSERT_EQ(basis.l_max(), 3);
	occ::qm::HartreeFock hf(basis);
	const occ::Mat S_occ = hf.compute_overlap_matrix();
	const dusk::SCF scf = dusk::make_scf(dusk_occ::input(mol, *BasisSetLibrary::get_basis_set("def2-tzvp")));
	ASSERT_EQ(scf.basis().nbf(), basis.nbf());
	const occ::Mat S = dusk_occ::to_occ(scf.overlap());
	EXPECT_LT((S - S_occ).cwiseAbs().maxCoeff(), 1e-12);
	//and back: the storage change is its own inverse
	EXPECT_EQ(dusk_occ::to_occ(dusk_occ::to_dusk(S_occ)), S_occ);
}

TEST(DuskBackendTests, RestrictedHartreeFockMatchesOcc)
{
	const pair_run r = solve(water, "def2-tzvp", occ::qm::SpinorbitalKind::Restricted, 1);
	ASSERT_TRUE(r.converged);
	EXPECT_NEAR(r.energy, r.occ_energy, 1e-7);
	EXPECT_LT(r.max_dD, 1e-5);
}

TEST(DuskBackendTests, UnrestrictedHartreeFockMatchesOcc)
{
	const pair_run r = solve(methyl, "def2-svp", occ::qm::SpinorbitalKind::Unrestricted, 2);
	ASSERT_TRUE(r.converged);
	EXPECT_NEAR(r.energy, r.occ_energy, 1e-7);
	EXPECT_LT(r.max_dD, 1e-5);
}

//From converged orbitals, handed through occ's wavefunction and back, there is nothing left to do
TEST(DuskBackendTests, WarmRestartStaysAtTheSolution)
{
	const pair_run r = solve(water, "def2-svp", occ::qm::SpinorbitalKind::Restricted, 1);
	ASSERT_TRUE(r.warm_converged);
	EXPECT_NEAR(r.warm_energy, r.energy, 1e-8);
	EXPECT_LE(r.warm_iterations, 3);
}

//A settings file as NoSpherA2 reads it arrives in dusk's options: presets, the ladder switches and the stabilisers
TEST(DuskBackendTests, SettingsFileReachesDusk)
{
	std::string log;
	const dusk::Options o = options_of("basis_set def2-SVP params 1 uhf tight slow_conv stability_check nososcf DIIS_size 6 broken_symmetry 30 check_hessian", log);
	EXPECT_EQ(o.convergence.energy, 5e-7);
	EXPECT_EQ(o.convergence.commutator, 5e-6);
	EXPECT_EQ(o.convergence.gradient, 3e-5);
	EXPECT_EQ(o.convergence.density_rms, 1e-7);
	EXPECT_EQ(o.convergence.density_max, 1e-6);
	EXPECT_TRUE(o.convergence.gradient_norm && o.convergence.commutator_ao);
	EXPECT_EQ(o.reference, dusk::Reference::Unrestricted);
	EXPECT_TRUE(o.stability_check && o.check_hessian && o.confirm_convergence);
	EXPECT_FALSE(o.soscf);
	EXPECT_EQ(o.second_order, dusk::SecondOrder::Auto);
	EXPECT_EQ(o.second_order_start, o.soscf_start);
	EXPECT_EQ(o.diis_size, 6);
	EXPECT_EQ(o.broken_symmetry, 30.0);
	EXPECT_EQ(o.damping, 0.8);
	EXPECT_TRUE(o.damping_eases);
	EXPECT_EQ(o.damping_off, 1e-5);
	EXPECT_EQ(o.level_shift, 1.0);
	EXPECT_EQ(o.level_shift_off, 1e-5);
	EXPECT_EQ(o.rescues, 3);
	EXPECT_EQ(o.full_rebuild_drop, 10.0);
	const dusk::Options r = options_of("basis_set def2-SVP params 1 rohf nodiis notrah fast_conv", log);
	EXPECT_EQ(r.reference, dusk::Reference::RestrictedOpen);
	EXPECT_EQ(r.diis_size, 1);
	EXPECT_EQ(r.second_order, dusk::SecondOrder::Off);
	EXPECT_EQ(r.damping, 0.0);
	EXPECT_EQ(r.level_shift, 0.0);
	EXPECT_NE(log.find("EXPERIMENTAL"), std::string::npos);
}

//A lambda scan through dusk_solver: the criterion evaluated incrementally in DIIS iterations and in full in second-order
//steps and the confirmation, as XCW_solver does; the criterion falls as lambda grows
TEST(DuskBackendTests, LambdaScanThroughTheXcwTerm)
{
	spdlog::set_level(spdlog::level::err);
	const occ::core::Molecule mol = occ::io::molecule_from_xyz_string(water);
	const occ::qm::AOBasis basis = SCF_solver::setup_basis(mol, "def2-svp");
	occ::qm::HartreeFock hf(basis);
	mock_coupling c;
	c.M = hf.compute_kinetic_matrix();
	options::SCF_settings s;
	std::ostringstream log;
	dusk::Input in = dusk_occ::input(mol, *BasisSetLibrary::get_basis_set("def2-svp"));
	in.options = dusk_occ::scf_options(s, log);
	const auto path = std::filesystem::temp_directory_path() / "nosphera2_dusk_scan_SCF.log";
	double last = 0;
	{
		XCW_log_writer writer(path.string());
		dusk_solver solver(std::move(in), mol, basis, writer, c);
		occ::qm::Wavefunction wfn;
		ASSERT_TRUE(solver.run(0.0, wfn, false));
		EXPECT_GT(c.incremental, 0);
		c.t = c.x + c.t - 0.05;
		last = 1e300;
		for (const double lambda : { 0.5, 1.0, 2.0 }) {
			ASSERT_TRUE(solver.run(lambda, wfn, true)) << "lambda " << lambda;
			EXPECT_LE(c.x * c.x, last + 1e-12);
			last = c.x * c.x;
		}
	}
	EXPECT_GT(c.full, 0);
	EXPECT_LE(c.restores, c.saves);
	EXPECT_LT(last, 0.05 * 0.05);
	std::filesystem::remove(path);
}
#endif
