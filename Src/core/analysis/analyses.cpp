#include "pch.h"
#include "analyses.h"
#include "convenience.h"
#include "cube.h"
#include "properties.h"
#include "eqc.h"
#include "bondwise_analysis.h"
#include "b2c.h"
#include "eli_family.h"
#include "nbo.h"
#include "nbo_run.h"
#include "citations.h"
#include "section_log.h"
#include "GridManager.h"
#include <future>

//Conceptual-DFT reactivity analysis; its table goes to stdout, not the log
static int run_fukui(options& opt, const char*, std::ostream&)
{
	fukui_analysis(opt, std::cout);
	return 0;
}

//EQC energy decomposition; its tables go to stdout and <stem>.eqc_log
static int run_eqc(options& opt, const char* argv0, std::ostream&)
{
	//Only a basis loaded by name needs occ's data directory; the gbw carries its own
	const bool have_data = ensure_occ_data_path(argv0);
	err_checkf(have_data || opt.eqc_basis.empty(), "-eqc_basis needs occ's data directory (set OCC_DATA_PATH)", std::cout);
	return eqc::run(opt);
}

//Full bonding analysis. Each stage rereads the wavefunction: RGBI and the NBO search modify it,
//and the basins must not see that
static int run_fba(options& opt, const char*, std::ostream& log)
{
	using namespace std;
	const auto read_wfn = [&opt]() {
		WFN w(opt.wfn);
		if (opt.ECP) w.set_has_ECPs(true, true, opt.ECP_mode);
		return w;
	};
	log << "\nFull bonding analysis of " << opt.wfn.string() << ": RGBI, bondwise Laplacian, QTAIM and ELI-D, then NBO/NPA with NRT\n" << endl;
	//NBO/NRT on its own thread with its own stream: RGBI and the basins print to std::cout, whose format flags
	//all threads share. RGBI is serial and the NBO search stops scaling past a few threads, so it takes the
	//cores the basins leave.
	std::ostringstream nbo_log;
	citations::cite(citations::Method::NAONPA, nbo_log);
	citations::cite(citations::Method::NBO, nbo_log);
	citations::cite(citations::Method::E2, nbo_log);
	citations::cite(citations::Method::NRT, nbo_log);
	const filesystem::path json = opt.wfn.parent_path() / (opt.wfn.stem().string() + ".native.nbo.json");
	//Refused here, not by write_nbo on the NBO thread, whose exit would tear down RGBI mid-run
	WFN nbo_wfn = read_wfn();
	refuse_unsupported_nbo_source(nbo_wfn, std::cout);
	eli_family::print_noon_header(eli_family::noon_summary(nbo_wfn), log);
	log << endl;
	auto nbo_done = std::async(std::launch::async, [&nbo_log, &opt, &json](WFN w) {
		const auto t0 = std::chrono::steady_clock::now();
		NboOptions nbo;
		nbo.debug = opt.debug;
		nbo.nrt = true;
		if (opt.threads > 0) nbo.threads = opt.threads;
		//The search's many small OpenMP regions would wait at barriers for threads the basin loops hold;
		//NRT's per-candidate loop is coarse and keeps the full count.
		nbo.search_threads = 1;
		NboResults r = native_nbo(w, nbo, nbo_log);
		r.name = opt.wfn.stem().string();
		//Not registered: cout belongs to RGBI and the basins on the main thread
		section_log::section nbo_file(opt.wfn, "nbo", "Natural bond orbitals: NPA, NBO, second-order E(2) and NRT",
			{ citations::Method::NAONPA, citations::Method::NBO, citations::Method::E2, citations::Method::NRT }, opt.no_date, false);
		std::ostringstream tables;
		print_nbo(r, tables);
		nbo_log << tables.str();
		nbo_file.stream() << tables.str() << "\n  Machine-readable copy: " << json.filename().string() << "\n";
		write_nbo_json(r, json);
		nbo_log << "wrote " << json.string() << std::endl;
		//NAO, search and E2 are in the json's timings; NRT is what this leaves over
		if (basin_timing_enabled())
			nbo_log << "  [timing] fba NBO/NPA/NRT thread: " << std::fixed << std::setprecision(2)
			        << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s" << std::endl;
	}, std::move(nbo_wfn));
	//-basin_timing: laps of the main thread, which runs alongside NBO/NRT
	basin_stage_timer fba_timer;
	const section_log::section rgbi_file(opt.wfn, "rgbi", "Roby-Gould Bond Indices (RGBI)", { citations::Method::RGBI }, opt.no_date),
		qtaim_file(opt.wfn, "qtaim", "QTAIM: critical points, atomic basins and delocalization indices",
			{ citations::Method::QTAIM, citations::Method::LIDI }, opt.no_date),
		eli_file(opt.wfn, "eli", "ELI-D: electron localizability basins", { citations::Method::ELID }, opt.no_date);
	{
		WFN w = read_wfn();
		Roby_information Roby(w, opt.rgbi_group_sets, !opt.rgbi_no_sym,
			opt.rgbi_orbital_basis == RGBIOrbitalBasis::ANO, opt.rgbi_EVs, opt.rgbi_theta,
			opt.rgbi_legacy_cutoff);
	}
	fba_timer.lap("fba RGBI");
	bondwise_laplacian_plots(opt.wfn);
	fba_timer.lap("fba bondwise Laplacian");
	ELI_analysis(read_wfn(), opt);
	fba_timer.lap("fba QTAIM and ELI-D");
	{
		//Charges of the grid partitions, on a standard molecular grid: the basin grids carry no density.
		//Labels as in the QTAIM/ELI-D tables, element plus 0-based wfn index.
		WFN w = read_wfn();
		w.delete_unoccupied_MOs();
		GridConfiguration config;
		config.accuracy = opt.accuracy;
		config.all_charges = true;
		GridManager grids(config);
		ivec atoms(w.get_ncen());
		std::iota(atoms.begin(), atoms.end(), 0);
		grids.setup3DGridsForMolecule(w, atoms);
		const PartitionResults charges = grids.calculatePartitionedCharges(w);
		svec labels;
		for (int a : atoms) labels.push_back(w.get_atom_label(a) + std::to_string(a));
		const section_log::section part_file(opt.wfn, "partition", "Atomic charges: Becke, Hirshfeld, TFVC, MBIS and EMBIS partitions",
			{ citations::Method::BeckeGrid, citations::Method::Hirshfeld, citations::Method::TFVC, citations::Method::MBIS, citations::Method::EMBIS }, opt.no_date);
		const section_log::tee part_tee("partition", "Charges by partitioning scheme");
		grids.printChargeTable(labels, w, atoms, std::cout, charges);
	}
	fba_timer.lap("fba partition charges");
	//Rethrows whatever stopped the NBO side
	nbo_done.get();
	fba_timer.lap("fba waiting for NBO/NRT");
	std::cout << "\n" << nbo_log.str() << "\nFull bonding analysis finished." << std::endl;
	return 0;
}

//Charges and atomic dipoles of the selected grid partition, table to the console, atoms to <stem>.xyzq
static int run_partition(options& opt, const char*, std::ostream&)
{
	using namespace std;
	err_checkf(opt.partition_type != PartitionType::RI, "-partition takes a grid partition: -Becke, -TFVC, -mbis, -embis or the Hirshfeld default", cout);
	WFN w(opt.wfn);
	if (opt.ECP) w.set_has_ECPs(true, true, opt.ECP_mode);
	w.delete_unoccupied_MOs();
	GridConfiguration config;
	config.accuracy = opt.accuracy;
	config.partition_type = opt.partition_type;
	config.debug = opt.debug;
	GridManager grids(config);
	ivec atoms(w.get_ncen());
	iota(atoms.begin(), atoms.end(), 0);
	grids.setup3DGridsForMolecule(w, atoms);
	const PartitionResults charges = grids.calculatePartitionedCharges(w);
	const vec2 Q = grids.calculatePartitionedMultipoles(w, 2);
	int s = PartitionResults::S_HIRSH;
	switch (opt.partition_type) {
	case PartitionType::Becke: s = PartitionResults::S_BECKE; break;
	case PartitionType::TFVC:  s = PartitionResults::S_TFVC; break;
	case PartitionType::MBIS:  s = PartitionResults::S_MBIS; break;
	case PartitionType::EMBIS: s = PartitionResults::S_EMBIS; break;
	default: break;
	}
	svec labels(atoms.size());
	for (int a : atoms) {
		labels[a] = w.get_atom_label(a);
		if (labels[a].find_first_of("0123456789") == string::npos) labels[a] += to_string(a + 1);
	}
	grids.printChargeTable(labels, w, atoms, cout, charges);
	ostringstream rows;
	double D[3] = { 0, 0, 0 };
	cout << "\nAtomic dipoles [e bohr] and traceless quadrupoles 1/2 sum q (3 r_i r_j - r^2 delta_ij) [e bohr^2]\n" << setw(10) << "Atom";
	for (const char* h : { "mu_x", "mu_y", "mu_z", "|mu|", "Q_xx", "Q_yy", "Q_zz", "Q_xy", "Q_xz", "Q_yz" }) cout << setw(10) << h;
	cout << "\n";
	for (int a : atoms) {
		//Q_lm = sum rho w r^l Y_lm about the nucleus; Y_1m is c (y, z, x), Y_2m is c (xy, yz, 3z^2-r^2, xz, x^2-y^2); the electrons are negative
		const double q = w.get_atom_charge(a) - charges.atom_charges[s][a];
		const double mu[3] = { -Q[a][3] / constants::c_3_4p, -Q[a][1] / constants::c_3_4p, -Q[a][2] / constants::c_3_4p };
		const double zz = -0.5 * Q[a][6] / constants::c_5_16p, xmy = -1.5 * Q[a][8] / constants::c_15_16p;
		const double th[6] = { 0.5 * (xmy - zz), -0.5 * (xmy + zz), zz, -1.5 * Q[a][4] / constants::c_15_4p, -1.5 * Q[a][7] / constants::c_15_4p, -1.5 * Q[a][5] / constants::c_15_4p };
		rows << left << setw(8) << labels[a] << right << fixed << setprecision(6);
		for (int x = 0; x < 3; x++) {
			rows << setw(14) << constants::bohr2ang(w.get_atom_coordinate(a, x));
			D[x] += q * w.get_atom_coordinate(a, x) + mu[x];
		}
		rows << setw(12) << q;
		for (int x = 0; x < 3; x++) rows << setw(12) << mu[x];
		for (int x = 0; x < 6; x++) rows << setw(12) << th[x];
		rows << "\n";
		cout << setw(10) << labels[a] << fixed << setprecision(4);
		for (int x = 0; x < 3; x++) cout << setw(10) << mu[x];
		cout << setw(10) << sqrt(mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2]);
		for (int x = 0; x < 6; x++) cout << setw(10) << th[x];
		cout << "\n";
	}
	const double Dn = sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]), debye = constants::bohr2ang(1.0) / constants::e_A2Debye;
	cout << "\nMolecular dipole sum(q R + mu) [e bohr]: " << setw(10) << D[0] << setw(10) << D[1] << setw(10) << D[2]
		<< "  |D| = " << Dn << " e bohr = " << Dn * debye << " Debye\n";
	const filesystem::path out = opt.wfn.parent_path() / (opt.wfn.stem().string() + ".xyzq");
	ofstream xyzq(out);
	err_checkf(xyzq.good(), "-partition: cannot write " + out.string(), cout);
	xyzq << atoms.size() << "\n" << config.getPartitionName() << " partition of " << opt.wfn.filename().string()
		<< ": label, x y z [Angstrom], charge [e], dipole x y z [e bohr], quadrupole xx yy zz xy xz yz [e bohr^2]; molecular dipole "
		<< fixed << setprecision(6) << D[0] << " " << D[1] << " " << D[2] << " [e bohr], " << Dn * debye << " Debye\n" << rows.str();
	cout << "wrote " << out.string() << endl;
	return 0;
}

//Basin analysis; the tables stay in the log, which is what the golden test reads
static int run_basins(options& opt, const char*, std::ostream&)
{
	WFN basins(opt.wfn);
	if (opt.ECP) basins.set_has_ECPs(true, true, opt.ECP_mode);
	const section_log::section qtaim_file(opt.wfn, "qtaim", "QTAIM: critical points, atomic basins and delocalization indices",
		{ citations::Method::QTAIM, citations::Method::LIDI }, opt.no_date);
	//-qtaim writes no .eli_log
	std::optional<section_log::section> eli_file;
	if (!opt.qtaim_only)
		eli_file.emplace(opt.wfn, "eli", "ELI-D: electron localizability basins", std::vector<citations::Method>{ citations::Method::ELID }, opt.no_date);
	ELI_analysis(basins, opt);
	return 0;
}

const std::vector<analysis>& analyses()
{
	static const std::vector<analysis> list = {
		{[](const options& o) -> const char* { return o.fukui_analysis_run ? "-fukui_analysis" : nullptr; }, true, run_fukui},
		{[](const options& o) -> const char* { return o.eqc ? "-eqc" : nullptr; }, true, run_eqc},
		{[](const options& o) -> const char* { return o.fba ? "-fba" : nullptr; }, false, run_fba},
		{[](const options& o) -> const char* { return o.partition_run ? "-partition" : nullptr; }, true, run_partition},
		{[](const options& o) -> const char* { return !o.eli_analysis_run ? nullptr : o.qtaim_only ? "-qtaim" : "-eli_analysis"; }, false, run_basins}};
	return list;
}
