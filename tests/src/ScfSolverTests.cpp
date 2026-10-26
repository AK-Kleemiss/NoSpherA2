//SCF_solver on its own, without a model term: the plain Hartree-Fock the XCW SCF is built on,
//against OCC's own SCF driver, and a warm restart that goes straight to the second-order steps
#include "pch.h"
#include <gtest/gtest.h>

#include "core/SCF_solver.h"
#include <occ/io/xyz.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#include <regex>

namespace
{
	const char* water = "3\n\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.469\nH 0.000 -0.757 -0.469\n";
	const char* methyl = "4\n\nC 0.000 0.000 0.000\nH 1.079 0.000 0.000\nH -0.540 0.935 0.000\nH -0.540 -0.935 0.000\n";

	//The `very_tight` convergence preset and the `normal_conv` damping and level shift of the XCW settings file
	options::SCF_settings tight_settings(const occ::qm::SpinorbitalKind kind, const int multiplicity) {
		options::SCF_settings s;
		s.quant_diff = 1e-7;
		s.max_diis_error = 1e-6;
		s.gradient = 1e-5;
		s.MaxP_diff = 1e-7;
		s.RMSP_diff = 1e-8;
		s.alpha = 0.5;
		s.level_shift = 0.5;
		s.diis_stop_damping = 1e-3;
		s.diis_stop_shift = 1e-2;
		s.max_scf_iterations = 100;
		s.hf_type = kind;
		s.multiplicity = multiplicity;
		return s;
	}

	struct plain_run {
		double energy = 0, warm_energy = 0, reference = 0;
		bool converged = false, warm_converged = false;
		std::string log;
	};

	//Iterations of every "***SCF converged in N iterations***" line of the log, in order
	std::vector<int> iterations(const std::string& log) {
		std::vector<int> n;
		const std::regex re("SCF converged in ([0-9]+) iterations");
		for (std::sregex_iterator it(log.begin(), log.end(), re), end; it != end; ++it) n.push_back(std::stoi((*it)[1]));
		return n;
	}

	//def2-SVP: SCF_solver from a fresh guess, then again from its own orbitals, and OCC's SCF for reference
	plain_run solve(const char* xyz, const occ::qm::SpinorbitalKind kind, const int multiplicity) {
		spdlog::set_level(spdlog::level::err);
		occ::core::Molecule mol = occ::io::molecule_from_xyz_string(xyz);
		mol.set_multiplicity(multiplicity);
		const occ::qm::AOBasis basis = SCF_solver::setup_basis(mol, "def2-svp");
		occ::qm::HartreeFock hf(basis);
		plain_run result;
		{
			occ::qm::SCF<occ::qm::HartreeFock> reference(hf, kind);
			reference.set_charge_multiplicity(0, multiplicity);
			reference.convergence_settings.energy_threshold = 1e-10;
			result.reference = reference.compute_scf_energy();
		}
		//one log per test: ctest runs them side by side
		const std::filesystem::path log_path = std::filesystem::temp_directory_path() /
			("nos_scf_solver_" + std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".log");
		{
			SCF_log_writer writer(log_path.string());
			SCF_solver solver(mol, basis, tight_settings(kind, multiplicity), writer);
			occ::qm::Wavefunction wfn;
			result.converged = solver.run(hf, wfn, false);
			result.energy = wfn.energy.total;
			result.warm_converged = solver.run(hf, wfn, true);
			result.warm_energy = wfn.energy.total;
		}
		std::ifstream in(log_path);
		result.log.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
		in.close();
		std::error_code ec;
		std::filesystem::remove(log_path, ec);
		return result;
	}
}

TEST(ScfSolverTests, PlainRestrictedHartreeFockMatchesOcc)
{
	const plain_run run = solve(water, occ::qm::SpinorbitalKind::Restricted, 1);
	ASSERT_TRUE(run.converged) << run.log;
	EXPECT_NEAR(run.energy, run.reference, 1e-7) << run.log;
	EXPECT_NE(run.log.find("Starting the SCF solver"), std::string::npos);
}

TEST(ScfSolverTests, PlainUnrestrictedHartreeFockMatchesOcc)
{
	const plain_run run = solve(methyl, occ::qm::SpinorbitalKind::Unrestricted, 2);
	ASSERT_TRUE(run.converged) << run.log;
	EXPECT_NEAR(run.energy, run.reference, 1e-7) << run.log;
}

//From converged orbitals the second-order steps start at once and have nothing left to do
TEST(ScfSolverTests, WarmRestartStaysAtTheSolution)
{
	const plain_run run = solve(water, occ::qm::SpinorbitalKind::Restricted, 1);
	ASSERT_TRUE(run.warm_converged) << run.log;
	EXPECT_NEAR(run.warm_energy, run.energy, 1e-8) << run.log;
	EXPECT_NE(run.log.find("Warm start: second-order steps"), std::string::npos) << run.log;
	const std::vector<int> n = iterations(run.log);
	ASSERT_EQ(n.size(), 2u) << run.log;
	EXPECT_LE(n[1], 3) << run.log;
}
