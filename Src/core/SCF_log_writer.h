#pragma once
#include "convenience.h"

// The SCF.log of a solver: the head of a run, the line of every iteration (its number and the total
// energy) and whatever the solver reports along the way through `log`.
class SCF_log_writer {

public:

	SCF_log_writer(const std::string& file = "SCF.log") { log.open(file); }
	virtual ~SCF_log_writer() = default;

	std::ofstream log;

	void WAEmpiRe() { log << WAEmpiRe_message() << std::endl; }

	// The head of the iteration table of one run of the solver called label
	void start_run(const std::string& label) {
		log << "Starting the " << label << " solver\n";
		log << "____________________________________________________________________________________________________\n";
		log << " Iteration\tTotal Energy (Eh)\n";
	}

	// One SCF iteration
	void iteration_line(const int iter, const double E_total) {
		log << " " << iter << "\t\t\t" << std::fixed << std::setprecision(9) << E_total << std::endl;
	}
};

// The XCW log: the crystallographic indicators (criterion, GooF, R1) next to the energy of every SCF
// iteration in SCF.log, and the row of every converged lambda on the console. XCW_solver and
// XCW_wrapper each fill in what they know in `v` and then ask for the line, so neither needs to see
// the other's data.
class XCW_log_writer : public SCF_log_writer {

public:

	using SCF_log_writer::SCF_log_writer;
	using SCF_log_writer::iteration_line;

	struct values {
		int iter = 0;
		double lambda = 0, criterion = 0, GooF2 = 0, R1 = 0, E_total = 0, penalty = 0, quant = 0, criterion_all = 0, R1_all = 0, A2 = 0;
		bool has_A2 = false;
		double wR2 = 0, wR2_all = 0;
	};

	values v;

	// The head of the table of one lambda step in SCF.log
	void start_lambda(const double lambda) {
		log << "Starting XCW SCF solver with lambda = " << std::fixed << std::setprecision(5) << lambda << "\n";
		log << "____________________________________________________________________________________________________\n";
		log << " Iteration\tCriterion\tGooF(F^2)\tR1(gt)\t\twR2\t\t\tTotal Energy\t\tPerturbation\tTarget quantity\n";
		log << "\t\t\t\t\t\t\t\t\t(%)\t\t\t\t\t\t(Eh)\t\t\t\t(a. u.)\t\t\t(a. u.)\n";
		log << "____________________________________________________________________________________________________\n";
	}

	// One SCF iteration: v.iter, criterion, GooF2, R1, wR2, E_total, penalty and quant
	void iteration_line() {
		log << " " << v.iter << "\t\t\t" << std::fixed << std::setprecision(4) << v.criterion << "\t\t" << v.GooF2 << "\t\t" << std::setprecision(5) << v.R1 << "\t\t" << v.wR2 << "\t\t" << std::fixed << std::setprecision(9) << v.E_total << "\t\t" << std::fixed << std::setprecision(3) << v.penalty << "\t\t\t" << std::fixed << std::setprecision(9) << v.quant << std::endl;
	}

	// The head of the console table, once per run
	static void console_header(const bool with_A2) {
		std::cout << "____________________________________________________________________________________________________________________________\n";
		std::cout << " Lambda\t\tCriterion\tGooF(F2)\tR1(gt)\t\tTotal Energy\t\tPerturbation\tTarget quantity\t\tCrit(all)\tR1(all)\t\twR2\t\twR2(all)";
		if (with_A2) std::cout << "\t\tA^2 (halt)";
		std::cout << "\n";
		std::cout << "\t\t\t\t\t\t\t\t\t(%)\t\t\t(Eh)\t\t\t\t(a. u.)\t\t\t(a. u.)\n";
		std::cout << "____________________________________________________________________________________________________________________________\n";
	}

	// The row of a converged lambda: v.lambda, criterion, GooF2, R1, E_total, quant, criterion_all, R1_all, wR2, wR2_all
	// and A2 when has_A2; to the console and, prefixed "XCW lambda row", to SCF.log
	void lambda_line() {
		std::ostringstream s;
		s << " " << std::fixed << std::setprecision(5) << v.lambda << "\t" << std::fixed << std::setprecision(4) << v.criterion << "\t\t" << v.GooF2 << "\t\t" << std::setprecision(5) << v.R1 << "\t\t" << std::fixed << std::setprecision(9) << v.E_total << "\t\t" << std::fixed << std::setprecision(3) << v.lambda * v.criterion << "\t\t\t" << std::fixed << std::setprecision(9) << v.quant
			<< "\t\t" << std::setprecision(4) << v.criterion_all << "\t\t" << std::setprecision(5) << v.R1_all << "\t\t" << v.wR2 << "\t\t" << v.wR2_all;
		if (v.has_A2) s << "\t\t" << std::setprecision(4) << v.A2;
		std::cout << s.str() << std::endl;
		log << "XCW lambda row:" << s.str() << std::endl;
	}
};
