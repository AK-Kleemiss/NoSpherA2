#pragma once
#include "convenience.h"

// Writes what an XCW run reports about its SCF: the line per iteration to SCF.log and the line
// per converged lambda to the console. SCF_wrapper and XCW_solver each fill in what they know in
// `v` and then ask for the line, so neither needs to see the other's data.
class SCF_log_writer {

public:

	SCF_log_writer(const std::string& file = "SCF.log") { log.open(file); }

	struct values {
		int iter = 0;
		double lambda = 0, criterion = 0, GooF2 = 0, R1 = 0, E_total = 0, penalty = 0, quant = 0, criterion_all = 0, R1_all = 0, A2 = 0;
		bool has_A2 = false;
	};

	std::ofstream log;
	values v;

	void WAEmpiRe() { log << WAEmpiRe_message() << std::endl; }

	// The head of the table of one lambda step in SCF.log
	void start_lambda(const double lambda) {
		log << "Starting XCW SCF solver with lambda = " << std::fixed << std::setprecision(5) << lambda << "\n";
		log << "____________________________________________________________________________________\n";
		log << " Iteration\t\tCriterion\tGooF(F^2)\tR1(gt)\t\tTotal Energy\t\tPerturbation\tTarget quantity\n";
		log << "\t\t\t\t\t\t\t\t\t(Eh)\t\t\t(a. u.)\t\t(a. u.)\n";
		log << "____________________________________________________________________________________\n";
	}

	// One SCF iteration: v.iter, criterion, GooF2, R1, E_total, penalty and quant
	void iteration_line() {
		log << "\t" << v.iter << "\t\t" << std::fixed << std::setprecision(4) << v.criterion << "\t\t" << v.GooF2 << "\t\t" << std::setprecision(5) << v.R1 << "\t\t" << std::fixed << std::setprecision(9) << v.E_total << "\t\t" << std::fixed << std::setprecision(3) << v.penalty << "\t\t" << std::fixed << std::setprecision(9) << v.quant << std::endl;
	}

	// The head of the console table, once per run
	static void console_header(const bool with_A2) {
		std::cout << "____________________________________________________________________________________\n";
		std::cout << " Lambda\t\tCriterion\tGooF(F2)\tR1(gt)\t\tTotal Energy\t\tPerturbation\tTarget quantity\t\tCrit(all)\tR1(all)";
		if (with_A2) std::cout << "\t\tA^2 (halt)";
		std::cout << "\n";
		std::cout << "\t\t\t\t\t\t\t\t(Eh)\t\t\t(a. u.)\t\t(a. u.)\n";
		std::cout << "____________________________________________________________________________________\n";
	}

	// The row of a converged lambda: v.lambda, criterion, GooF2, R1, E_total, quant, criterion_all, R1_all and A2 when has_A2
	void lambda_line() const {
		std::cout << std::fixed << std::setprecision(5) << v.lambda << "\t\t" << std::fixed << std::setprecision(4) << v.criterion << "\t\t" << v.GooF2 << "\t\t" << std::setprecision(5) << v.R1 << "\t\t" << std::fixed << std::setprecision(9) << v.E_total << "\t\t" << std::fixed << std::setprecision(3) << v.lambda * v.criterion << "\t\t" << std::fixed << std::setprecision(9) << v.quant
			<< "\t\t" << std::setprecision(4) << v.criterion_all << "\t\t" << std::setprecision(5) << v.R1_all;
		if (v.has_A2) std::cout << "\t\t" << std::setprecision(4) << v.A2;
		std::cout << std::endl;
	}
};
