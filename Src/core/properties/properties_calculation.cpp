#include "pch.h"
#include "wfn_class.h"
#include "properties.h"
#ifdef NOSPHERA2_USE_GPU
#include "aux_density_gpu.h"
#include "spherical_density_gpu.h"
#endif
#include "convenience.h"
#include "spherical_density.h"
#include "cube.h"
#include "constants.h"
#include "ibo.h"
#include "nbo.h"
#include "GridManager.h"
#include "isosurface.h"
#include "SALTED_predictor.h"
#include "crystal_energies.h"
#include "b2c.h"
#include "density_source.h"
#include "citations.h"
#include "eli_family.h"

void properties_calculation(options &opt)
{
	using namespace std;
	ofstream log2("NoSpherA2_cube.log", ios::out);
	auto _coutbuf = std::cout.rdbuf(log2.rdbuf()); // save and redirect
	log2 << NoSpherA2_message(opt.no_date);
	if (!opt.no_date)
	{
		log2 << build_date;
	}
	log2.flush();

	err_checkf(opt.wfn != "", "Error, no wfn file specified!", log2);
	WFN wavy(opt.wfn);
	//Without MOs (xyz input) rho, ESP, Laplacian and ELI come from the SALTED prediction; with -ri_fit <basis> from the RI fit of the orbitals
	std::unique_ptr<Gaussian_Molecule> ml;
	if (Gaussian_Molecule::requested(wavy, opt)) {
		if (wavy.get_nmo() == 0) log2 << "No orbitals in " << opt.wfn << ", rho, ESP, Laplacian and ELI from the SALTED model " << opt.salted_model_dir << endl;
		else log2 << "rho, ESP, Laplacian and ELI from the RI fit of " << opt.wfn << endl;
		ml = std::make_unique<Gaussian_Molecule>(wavy, opt);
		//Orbital-only properties on the fitted density: say so and drop them, the density properties still run
		auto refuse = [&](bool& flag, const char* what) {
			if (!flag) return;
			log2 << "ERROR: " << what << " need orbitals, the fitted density has none - skipped" << endl;
			flag = false;
		};
		refuse(opt.properties.elf, "ELF");
		if (wavy.get_nmo() == 0) {
			refuse(opt.properties.s_rho, "spin density");
			refuse(opt.properties.fukui, "Fukui functions");
			refuse(opt.properties.all_mos, "MOs");
			bool mos = !opt.properties.MO_numbers.empty();
			refuse(mos, "MOs");
			opt.properties.MO_numbers.clear();
		}
	}
	auto ml_esp = [&](const d3& p) { return ml->esp(p); };
	if (opt.debug)
		log2 << "Starting calculation of properties" << endl;
	if (opt.properties.all_mos)
		for (int mo = 0; mo < wavy.get_nmo(); mo++)
			opt.properties.MO_numbers.push_back(mo);
	if (opt.debug)
		log2 << "Size of MOs: " << opt.properties.MO_numbers.size() << endl;

	vec2 cell_matrix;
	cell_matrix.resize(3);
	for (int i = 0; i < 3; i++)
		cell_matrix[i].resize(3, 0.0);
	if (opt.debug)
	{
		log2 << "cif|resolution|res(bohr)|radius|rad(bohr): " << opt.cif << "|" << opt.properties.resolution << "|" << constants::ang2bohr(opt.properties.resolution) << "|" << opt.properties.radius << "|" << constants::ang2bohr(opt.properties.radius) << endl;
		for (int a = 0; a < wavy.get_ncen(); a++)
			log2 << "Atom " << a << " at " << wavy.get_atom_coordinate(a, 0) << " " << wavy.get_atom_coordinate(a, 1) << " " << wavy.get_atom_coordinate(a, 2) << "\n";
	}
	if (opt.properties.esp_isosurface > 0 && opt.properties.radius < 2.5)
	{
		log2 << "Resetting Radius to at least 2.5 for the isosurface!" << endl;
		opt.properties.radius = 2.5;
	}
	if (opt.cif != "")
		readxyzMinMax_fromCIF(opt.cif, opt.properties, cell_matrix);
	else
	{
		//no even-point rounding here: the step below is the requested resolution and not the span over the
		//point count, so one more point would enlarge the box instead of putting its centre on a grid plane
		readxyzMinMax_fromWFN(wavy, opt.properties, false);
		for (int i = 0; i < 3; i++)
			cell_matrix[i][i] = constants::ang2bohr(opt.properties.resolution);
	}
	if (opt.debug)
	{
		log2 << "MinMax: ";
		for (int i = 0; i < 6; i++)
			log2 << setw(14) << scientific << opt.properties.MinMax[i];
		log2 << endl;
		log2 << "Steps: ";
		for (int i = 0; i < 3; i++)
			log2 << setw(14) << scientific << opt.properties.NbSteps[i];
		log2 << endl;
		log2 << "Cell Matrix:" << endl;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
				log2 << setw(14) << scientific << cell_matrix[i][j];
			log2 << "\n";
		}
	}
	std::vector<cube> cubes;
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), true);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.rdg);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.elf);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.eli);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.lap);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.esp);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), true);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.hdef);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.def);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.hirsh);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.s_rho);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.hdef || opt.properties.hirsh);
	// Fukui_plus / Fukui_minus / Fukui_zero / Dual_Descriptor. All four are
	// allocated together: they come out of a single grid pass over the same two
	// frontier orbitals, so computing one and not the others saves nothing, and
	// an unloaded cube has no storage allocated at all (see cube's grow_values
	// constructor) which would make the shared evaluation path unsafe.
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.fukui);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.fukui);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.fukui);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.fukui);
	//-eli_spin: the members the orbitals carry, aa only for a closed shell (bb identical, triplet a constant multiple)
	std::vector<eli_family::Member> eli_members;
	if (opt.properties.eli_spin)
		eli_members = eli_family::eli_variants_for(wavy);
	auto has_member = [&](const eli_family::Member m) { return std::find(eli_members.begin(), eli_members.end(), m) != eli_members.end(); };
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.eli_spin);
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), has_member(eli_family::Member::eli_d_bb));
	cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), has_member(eli_family::Member::eli_d_triplet));

	for (cube &cube : cubes)
		cube.give_parent_wfn(wavy);


	for (int i = 0; i < 3; i++)
	{
		for (cube &cube : cubes)
			cube.set_origin(i, opt.properties.MinMax[i]);

		for (int j = 0; j < 3; j++)
		{
			for (cube &cube : cubes)
				cube.set_vector(i, j, cell_matrix[i][j]);
		}
	}
	for (cube &cube : cubes)
		cube.calc_dv();
	if (opt.debug)
		log2 << "Origins etc. are set up" << endl;
	cubes[cube_type::Rho].set_comment1("Calculated density using NoSpherA2");
	cubes[cube_type::RDG].set_comment1("Calculated reduced density gradient using NoSpherA2");
	cubes[cube_type::Elf].set_comment1("Calculated electron localization function using NoSpherA2");
	cubes[cube_type::Eli].set_comment1("Calculated same-spin electron localizability indicator using NoSpherA2");
	cubes[cube_type::Lap].set_comment1("Calculated laplacian of electron density using NoSpherA2");
	cubes[cube_type::ESP].set_comment1("Calculated electrostatic potential using NoSpherA2");
	cubes[cube_type::MO_val].set_comment1("Calcualted MO values using NoSpherA2");
	cubes[cube_type::HDEF].set_comment1("Calculated Atomic Hirshfeld deformation density values using NoSpherA2");
	cubes[cube_type::DEF].set_comment1("Calculated static deformation density values using NoSpherA2");
	cubes[cube_type::Hirsh].set_comment1("Calculated Hirshfeld atom density values using NoSpherA2");
	cubes[cube_type::Spin_Density].set_comment1("Calculated spin density using NoSpherA2");
	cubes[cube_type::Fukui_plus].set_comment1("Calculated Fukui f+ (LUMO density, nucleophilic attack) using NoSpherA2");
	cubes[cube_type::Fukui_minus].set_comment1("Calculated Fukui f- (HOMO density, electrophilic attack) using NoSpherA2");
	cubes[cube_type::Fukui_zero].set_comment1("Calculated Fukui f0 (radical attack) using NoSpherA2");
	cubes[cube_type::Dual_Descriptor].set_comment1("Calculated dual descriptor f+ minus f- using NoSpherA2");
	cubes[cube_type::Eli_aa].set_comment1("Calculated alpha-alpha ELI-D using NoSpherA2");
	cubes[cube_type::Eli_bb].set_comment1("Calculated beta-beta ELI-D using NoSpherA2");
	cubes[cube_type::Eli_triplet].set_comment1("Calculated triplet ELI-D using NoSpherA2");
	for (auto cube : cubes)
		cube.set_comment2("from " + wavy.get_path().string());

	cubes[cube_type::Rho].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
	cubes[cube_type::RDG].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rdg.cube");
	cubes[cube_type::Elf].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_elf.cube");
	cubes[cube_type::Eli].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_eli.cube");
	cubes[cube_type::Lap].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_lap.cube");
	cubes[cube_type::ESP].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_esp.cube");
	cubes[cube_type::DEF].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_def.cube");
	cubes[cube_type::Hirsh].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_hirsh.cube");
	cubes[cube_type::Spin_Density].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_s_rho.cube");
	cubes[cube_type::Fukui_plus].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_fukui_plus.cube");
	cubes[cube_type::Fukui_minus].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_fukui_minus.cube");
	cubes[cube_type::Fukui_zero].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_fukui_zero.cube");
	cubes[cube_type::Dual_Descriptor].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_dual_descriptor.cube");
	cubes[cube_type::Eli_aa].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_eli_aa.cube");
	cubes[cube_type::Eli_bb].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_eli_bb.cube");
	cubes[cube_type::Eli_triplet].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_eli_triplet.cube");

	log2 << "\nCalculating:" << endl;
	if (opt.properties.hdef || opt.properties.def || opt.properties.hirsh)
		log2 << "Rho, ";
	if (opt.properties.hdef || opt.properties.def || opt.properties.hirsh)
		log2 << "Spherical Rho, ";
	if (opt.properties.def)
		log2 << "Static deformation density, ";
	if (opt.properties.hdef)
		log2 << "Hirshfeld deformation density, ";
	if (opt.properties.hirsh)
		log2 << "Hirshfeld density, ";
	if (opt.properties.lap)
		log2 << "Laplacian, ";
	if (opt.properties.eli)
		log2 << "ELI, ";
	if (opt.properties.elf)
		log2 << "ELF, ";
	if (opt.properties.rdg)
		log2 << "RDG, ";
	if (opt.properties.esp)
		log2 << "ESP, ";
	if (opt.properties.MO_numbers.size() != 0)
		log2 << "MOs, ";
	if (opt.properties.s_rho)
		log2 << "Spin density, ";
	if (opt.properties.fukui)
		log2 << "Fukui functions and dual descriptor, ";
	if (opt.properties.eli_spin)
		log2 << "spin-resolved ELI-D, ";
	log2 << endl;

	//Which paper each of those grids implements.  Same guards as the list above.
	if (opt.properties.hdef || opt.properties.def || opt.properties.hirsh)
		citations::cite(citations::Method::Hirshfeld, log2);
	if (opt.properties.eli)
		citations::cite(citations::Method::ELID, log2);
	if (opt.properties.elf)
		citations::cite(citations::Method::ELF, log2);
	if (opt.properties.rdg)
		citations::cite(citations::Method::NCI, log2);
	if (opt.properties.esp)
		citations::cite(citations::Method::ESP, log2);
	if (opt.properties.fukui)
		citations::cite(citations::Method::Fukui, log2);
	if (opt.properties.eli_spin) {
		if (!opt.properties.eli) citations::cite(citations::Method::ELID, log2);
		citations::cite(citations::Method::ELIFamily, log2);
	}

	log2 << "Calculating for " << fixed << setprecision(0) << opt.properties.NbSteps[0] * opt.properties.NbSteps[1] * opt.properties.NbSteps[2] << " Gridpoints." << endl;

	// -esp_isosurface with -cif marches its own box below, so a run asking only for that never reads the cell rho
	const properties_options &p = opt.properties;
	if (p.rho || p.lap || p.eli || p.elf || p.rdg || p.esp || p.def || p.hdef || p.hirsh || p.integral_accuracy != -1 || (p.esp_isosurface > 0 && opt.cif == "")) {
		if (ml) Calc_Rho(cubes[cube_type::Rho], *ml, opt.properties.radius, log2, opt.cif != "");
		else Calc_Rho(cubes[cube_type::Rho], wavy, opt.properties.radius, log2, opt.cif != "");
	}

	if (opt.properties.integral_accuracy != -1) {
		log2 << "Refining grid files to integral accuracy of " << opt.properties.integral_accuracy << " ..." << flush;
		vec2 d(16, vec(wavy.get_ncen(), 0.0));
		vec phi(wavy.get_nmo(true), 0.0);
		cubes[cube_type::Rho].adaptive_refine([&wavy](const d3 &pos) { return wavy.compute_dens(pos); }, opt.properties.integral_accuracy, 8, 3);
	}

	if (opt.properties.MO_numbers.size() != 0)
		for (int i = 0; i < opt.properties.MO_numbers.size(); i++)
		{
			log2 << "Calcualting MO: " << opt.properties.MO_numbers[i] << endl;
			cubes[cube_type::MO_val].set_zero();
			cubes[cube_type::MO_val].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_MO_" + to_string(opt.properties.MO_numbers[i]) + ".cube");
			Calc_MO(cubes[cube_type::MO_val], opt.properties.MO_numbers[i], wavy, opt.properties.radius, log2, opt.cif != "");
			cubes[cube_type::MO_val].write_file(true);
		}
	if (!opt.properties.ibo_cube.empty())
	{
		const IBOResult ibo = intrinsic_bond_orbitals(wavy);
		print_ibo(ibo, wavy, log2);
		const WFN local = ibo_wfn(wavy, ibo);
		for (const int k : ibo_selection(ibo, opt.properties.ibo_cube))
		{
			log2 << "Calculating IBO " << k + 1 << endl;
			cubes[cube_type::MO_val].set_zero();
			cubes[cube_type::MO_val].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_IBO_" + to_string(k + 1) + ".cube");
			Calc_MO(cubes[cube_type::MO_val], ibo.mos[k], local, opt.properties.radius, log2, opt.cif != "");
			cubes[cube_type::MO_val].write_file(true);
		}
	}
	if (!opt.properties.nbo_cube.empty() || !opt.properties.nlmo_cube.empty())
	{
		NboOptions nbo_opt;
		nbo_opt.threads = std::max(opt.threads, 0);
		nbo_opt.nrt_parent = opt.properties.nbo_lewis;
		NboOrbitalSpace space;
		//write_nbo renormalises the stored basis; the cubes take their primitives from the untouched wavy
		WFN copy(wavy);
		NboResults r = native_nbo(copy, nbo_opt, log2, &space);
		r.name = wavy.get_path().stem().string();
		print_nbo(r, log2);
		for (const bool nlmo : {false, true})
		{
			const std::string &spec = nlmo ? opt.properties.nlmo_cube : opt.properties.nbo_cube;
			if (spec.empty())
				continue;
			std::vector<std::pair<int, std::string>> picked;
			const WFN local = nbo_cube_wfn(wavy, r, space, spec, nlmo, picked, log2);
			for (int k = 0; k < static_cast<int>(picked.size()); k++)
			{
				log2 << "Calculating " << (nlmo ? "NLMO " : "NBO ") << picked[k].first << " " << picked[k].second << endl;
				cubes[cube_type::MO_val].set_zero();
				cubes[cube_type::MO_val].set_comment1(std::string(nlmo ? "NLMO " : "NBO ") + to_string(picked[k].first) + " " + picked[k].second + " (NoSpherA2 in-house NBO)");
				cubes[cube_type::MO_val].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + (nlmo ? "_NLMO_" : "_NBO_") + to_string(picked[k].first) + ".cube");
				Calc_MO(cubes[cube_type::MO_val], k, local, opt.properties.radius, log2, opt.cif != "");
				cubes[cube_type::MO_val].write_file(true);
			}
		}
	}

	// The Fukui functions need the LUMO, which is an UNOCCUPIED orbital. This
	// block must therefore stay above the delete_unoccupied_MOs() call below -
	// moving it further down leaves the virtual space empty and silently yields
	// an all-zero f+ cube rather than an error that names the cause.
	if (opt.properties.fukui)
	{
		int homo = -1, lumo = -1;
		bool unrestricted = false;
		const bool found = find_frontier_orbitals(wavy, homo, lumo, unrestricted);
		if (!found)
		{
			log2 << "WARNING: Could not identify a HOMO/LUMO pair (found HOMO index " << homo
				 << ", LUMO index " << lumo << "). "
				 << "The wavefunction has no virtual orbitals, or no occupied ones. "
				 << "Skipping the Fukui functions." << endl;
		}
		else
		{
			log2 << "Frontier orbitals: HOMO = MO " << homo
				 << " (energy " << fixed << setprecision(6) << wavy.get_MO_energy(homo)
				 << ", occupation " << setprecision(3) << wavy.get_MO_occ(homo) << "), "
				 << "LUMO = MO " << lumo
				 << " (energy " << setprecision(6) << wavy.get_MO_energy(lumo)
				 << ", occupation " << setprecision(3) << wavy.get_MO_occ(lumo) << ")" << endl;
			if (unrestricted)
				log2 << "WARNING: unrestricted wavefunction. The frontier pair is taken across both "
					 << "spin manifolds, so the result is not a spin-resolved Fukui function." << endl;
			log2 << "Calculating Fukui functions..." << endl;
			Calc_Fukui(cubes, wavy, homo, lumo, opt.properties.radius, log2, opt.cif != "");

			// Written here, next to the MO cubes, rather than in the common
			// write block further down: everything below this point runs after
			// the virtual orbitals have been discarded.
			cubes[cube_type::Fukui_plus].write_file(true);
			cubes[cube_type::Fukui_minus].write_file(true);
			cubes[cube_type::Fukui_zero].write_file(true);
			cubes[cube_type::Dual_Descriptor].write_file(true);

			// Integrated norms. For the exact Fukui function each of f+ and f-
			// integrates to 1; in the frontier-orbital approximation this holds
			// only as far as the orbital is normalised inside the evaluated
			// radius, so a value below 1 means the grid or the radius is
			// clipping the orbital rather than that anything is wrong with the
			// theory. Cheap, and the honest self-check to print.
			// cube::sum() already multiplies by the volume element dv, so it
			// returns the integral rather than a bare sum. Multiplying by
			// get_dv() here as well would apply dv twice and make the reported
			// norm scale with the grid spacing.
			const double int_plus = cubes[cube_type::Fukui_plus].sum();
			const double int_minus = cubes[cube_type::Fukui_minus].sum();
			log2 << "Integrated f+ over the grid: " << fixed << setprecision(4) << int_plus
				 << " (exact value 1.0)" << endl;
			log2 << "Integrated f- over the grid: " << fixed << setprecision(4) << int_minus
				 << " (exact value 1.0)" << endl;

			// Machine-readable summary beside the cubes. This exists separately
			// from the log because the log carries wall-clock timings, which are
			// not reproducible and so cannot be used as golden-file test output.
			const std::filesystem::path summary_path =
				(wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_fukui.dat";
			ofstream summary(summary_path, ios::out);
			summary << "# NoSpherA2 Fukui function summary\n";
			summary << "# Frontier-orbital (frozen-orbital) approximation:\n";
			summary << "#   f+ = |psi_LUMO|^2, f- = |psi_HOMO|^2, f0 = (f+ + f-)/2, df = f+ - f-\n";
			summary << "# The integrals are over the evaluated grid only; each equals 1 for the\n";
			summary << "# exact function, so a smaller value means the grid or radius is clipping\n";
			summary << "# the orbital rather than that anything is wrong.\n";
			summary << fixed;
			summary << "HOMO_index " << homo << "\n";
			summary << "HOMO_energy " << setprecision(6) << wavy.get_MO_energy(homo) << "\n";
			summary << "HOMO_occupation " << setprecision(6) << wavy.get_MO_occ(homo) << "\n";
			summary << "LUMO_index " << lumo << "\n";
			summary << "LUMO_energy " << setprecision(6) << wavy.get_MO_energy(lumo) << "\n";
			summary << "LUMO_occupation " << setprecision(6) << wavy.get_MO_occ(lumo) << "\n";
			summary << "unrestricted " << (unrestricted ? 1 : 0) << "\n";
			summary << "integral_f_plus " << setprecision(6) << int_plus << "\n";
			summary << "integral_f_minus " << setprecision(6) << int_minus << "\n";

			// Condensed (atom-summed) Fukui functions under all five partitions.
			// Uses its own Becke-style atomic integration grid rather than the
			// cube grid: the cube is a display object, clipped by -radius and
			// coarse enough that integrating it per atom would be a much worse
			// number than the one GridManager already knows how to produce.
			log2 << "\nCalculating condensed Fukui functions..." << endl;
			_time_point cf_start = get_time();
			// Default cell with pbc = 0, i.e. the MOLECULAR treatment - the same
			// choice integrator.cpp makes for atomic charges. The condensed
			// values therefore describe the molecule in the wavefunction, even
			// when a CIF was supplied for the cube grid.
			const CondensedFukuiResults condensed =
				Calc_Condensed_Fukui(wavy, homo, lumo, cell(), opt.accuracy, log2);
			_time_point cf_end = get_time();
			if (condensed.valid) {
				print_condensed_fukui(condensed, log2);
				log2 << "Time for condensed Fukui functions: "
					 << fixed << setprecision(1) << get_sec(cf_start, cf_end) << " s" << endl;
				// Same values into the machine-readable summary, Hirshfeld first
				// because it is the scheme the condensed-Fukui literature uses.
				static const int ord[5] = {
					PartitionResults::CHARGE_ORDER::S_HIRSH,
					PartitionResults::CHARGE_ORDER::S_BECKE,
					PartitionResults::CHARGE_ORDER::S_TFVC,
					PartitionResults::CHARGE_ORDER::S_MBIS,
					PartitionResults::CHARGE_ORDER::S_EMBIS
				};
				static const char *ord_name[5] = { "hirshfeld", "becke", "tfvc", "mbis", "embis" };
				summary << "# condensed Fukui: atom f+ f- df, one block per partition\n";
				for (int s = 0; s < 5; s++) {
					summary << "partition " << ord_name[s] << "\n";
					for (size_t a = 0; a < condensed.labels.size(); a++) {
						const double fp = condensed.f_plus[ord[s]][a];
						const double fm = condensed.f_minus[ord[s]][a];
						summary << "  " << condensed.labels[a]
								<< " " << setprecision(6) << fp
								<< " " << setprecision(6) << fm
								<< " " << setprecision(6) << (fp - fm) << "\n";
					}
				}
			}
			else {
				log2 << "Condensed Fukui functions could not be calculated." << endl;
			}

			summary.close();
			log2 << "Wrote Fukui summary to " << summary_path.filename().string() << endl;
		}
	}

	wavy.delete_unoccupied_MOs();
	wavy.delete_Qs();

	for (int i = 1; i < cubes.size(); i++) {
		if (cubes[i].get_loaded())
			cubes[i].resize(cubes[cube_type::Rho].get_sizes());
		cubes[i].set_vectors(cubes[cube_type::Rho].get_vectors());
	}

	if (opt.properties.hdef || opt.properties.def || opt.properties.hirsh)
	{
		cubes[cube_type::spherical_density].resize(cubes[cube_type::Rho].get_sizes());
		for (int i = 0; i < 3; i++)
		{
			cubes[cube_type::spherical_density].set_origin(i, opt.properties.MinMax[i]);
			for (int j = 0; j < 3; j++)
				cubes[cube_type::spherical_density].set_vector(i, j, cell_matrix[i][j]);
		}
		// -def alone used to skip this and hand back rho as the deformation density
		if (opt.properties.hdef || opt.properties.def || opt.properties.hirsh)
		{
			log2 << "Calculating spherical Rho...";
			Calc_Spherical_Dens(cubes[cube_type::spherical_density], wavy, opt.properties.radius, log2, opt.cif != "");
			log2 << " ...done!" << endl;
		}

		if (opt.properties.def)
		{
			log2 << "Calculating static deformation density...";
			Calc_Static_Def(cubes, wavy, opt.properties.radius, log2, opt.cif != "");
			log2 << " ...done!" << endl;
		}

		if (opt.properties.hdef)
		{
			for (int a = 0; a < wavy.get_ncen(); a++)
			{
				log2 << "Calcualting Hirshfeld deformation density for atom: " << a << endl;
				cubes[cube_type::HDEF].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_HDEF_" + to_string(a) + ".cube");
				Calc_Hirshfeld(cubes, wavy, opt.properties.radius, a, log2, opt.cif != "");
				cubes[cube_type::HDEF].write_file(true);
				cubes[cube_type::HDEF].set_zero();
			}
		}

		if (opt.properties.hirsh)
		{
			log2 << "Calcualting Hirshfeld density for atom: " << opt.properties.hirsh_number << endl;
			Calc_Hirshfeld_atom(cubes, wavy, opt.properties.radius, opt.properties.hirsh_number, log2, opt.cif != "");
			log2 << "..done!" << endl;
		}
	}

	if (opt.properties.lap || opt.properties.eli || opt.properties.elf || opt.properties.rdg || opt.properties.esp) {
		if (ml) Calc_Prop(cubes, *ml, opt.properties.radius, log2, opt.no_date, opt.cif != "");
		else Calc_Prop(cubes, wavy, opt.properties.radius, log2, opt.no_date, opt.cif != "");
	}

	if (opt.properties.s_rho)
		Calc_S_Rho(cubes[cube_type::Spin_Density], wavy, log2, opt.no_date);
	if (opt.properties.eli_spin) {
		if (!cubes[cube_type::Eli_bb].get_loaded())
			log2 << "Spin-resolved ELI-D: the wavefunction carries no spin polarisation, so ELI-D(bb) = ELI-D(aa) and the triplet member is a constant multiple of it; writing aa only" << endl;
		log2 << "Calculating spin-resolved ELI-D..." << endl;
		Calc_Eli_Spin(cubes, wavy, log2, opt.no_date);
		cubes[cube_type::Eli_aa].write_file(true);
		if (cubes[cube_type::Eli_bb].get_loaded()) cubes[cube_type::Eli_bb].write_file(true);
		if (cubes[cube_type::Eli_triplet].get_loaded()) cubes[cube_type::Eli_triplet].write_file(true);
	}

	log2 << "Writing cubes to Disk..." << flush;
	if (opt.properties.rdg)
	{
		cubes[cube_type::Rho].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_signed_rho.cube");
		cubes[cube_type::Rho].write_file(true);
		cubes[cube_type::Rho].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
		cubes[cube_type::Rho].write_file(true, true);
	}
	else if (opt.properties.rho || opt.properties.lap || opt.properties.eli || opt.properties.elf || opt.properties.esp)
		cubes[cube_type::Rho].write_file(true);
	if (opt.properties.rdg)
		cubes[cube_type::RDG].write_file(true);
	if (opt.properties.lap)
		cubes[cube_type::Lap].write_file(true);
	if (opt.properties.elf)
		cubes[cube_type::Elf].write_file(true);
	if (opt.properties.eli)
		cubes[cube_type::Eli].write_file(true);
	if (opt.properties.s_rho)
		cubes[cube_type::Spin_Density].write_file(true);
	if (opt.properties.def)
	{
		cubes[cube_type::DEF].write_file(true);
		cubes[cube_type::Rho].write_file(true);
	}
	if (opt.properties.hirsh)
		cubes[cube_type::Hirsh].write_file(true);

	log2 << " done!" << endl;

	if (opt.properties.esp)
	{
		log2 << "Calculating ESP..." << flush;
		WFN temp = wavy;
		temp.delete_unoccupied_MOs();
		temp.delete_Qs();
		_time_point e0 = get_time();
		if (ml) Calc_ESP(cubes[cube_type::ESP], *ml, opt.properties.radius, opt.no_date, log2, opt.cif != "");
		else Calc_ESP(cubes[cube_type::ESP], temp, opt.properties.radius, opt.no_date, log2, opt.cif != "");
		_time_point e1 = get_time();
		log2 << "Writing cube to Disk..." << flush;
		cubes[cube_type::ESP].write_file(true);
		_time_point e2 = get_time();
		log2 << "  done! ESP " << get_msec(e0, e1) << " ms, cube written in " << get_msec(e1, e2) << " ms" << endl;
	}
	if (opt.properties.esp_isosurface > 0)
	{
		log2 << std::defaultfloat << std::setprecision(4) << "Colouring the rho = " << opt.properties.esp_isosurface << " au isosurface by the ESP..." << endl;
		cube box;
		if (opt.cif != "")
		{ // the cell grid clips the molecule, so take rho on its own box
			properties_options box_opts = opt.properties;
			box = box_cube(wavy, box_opts);
			if (ml) Calc_Rho(box, *ml, box_opts.radius, log2, false);
			else Calc_Rho(box, wavy, box_opts.radius, log2, false);
		}
		// the phases' cost varies widely between molecules, so each reports its time
		_time_point t0 = get_time();
		std::vector<Triangle> triangles = marchingCubes(opt.cif != "" ? box : cubes[cube_type::Rho], opt.properties.esp_isosurface);
		_time_point t1 = get_time();
		log2 << "Found " << triangles.size() << " triangles in " << get_msec(t0, t1) << " ms" << endl;
		// the per-face ESP below is the long pole, so the bare isosurface goes out first in its own file: the coloured one
		// must only ever appear complete
		const std::string esp_stem = (wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho_esp";
		writeColourObj(esp_stem + "_shape.obj", triangles);
		{ ofstream stage1(esp_stem + "_shape.obj.stage1"); stage1 << triangles.size() << "\n"; }
		_time_point t2 = get_time();
		if (ml) colour_by_ESP(triangles, surface_ESP(triangles, ml_esp), log2);
		else colour_by_ESP(triangles, wavy, log2);
		_time_point t3 = get_time();
		writeColourObj(esp_stem + ".obj", triangles);
		_time_point t4 = get_time();
		log2 << "shape mesh " << get_msec(t1, t2) << " ms, ESP on the faces " << get_msec(t2, t3)
			 << " ms, coloured mesh " << get_msec(t3, t4) << " ms" << endl;
	}
	// return output tostd::cout
	std::cout.rdbuf(_coutbuf);
	log2.close();
	std::cout << "Properties calculation done!" << std::endl;
}
