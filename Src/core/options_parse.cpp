#include "pch.h"
#include "convenience.h"
#include "cell.h"
#include "tsc_block.h"
#include "tshb_block.h"
#include "scattering_factors.h"
#include "tsc_label_converter.h"
#include "dev_diagnostics.h"
#include "spherical_density.h"
#include "AtomGrid.h"
#include "npy.h"
#include "SALTED_utilities.h"
#include "sphere_lebedev_rule.h"
#include "nos_math.h"
#include "featomic.hpp"
#include "metatensor.h"
#include "GridManager.h"
#include "integrator.h"
#include "properties.h"
#include "wfn_class.h"
#include "atoms.h"
#include "basis_set.h"
#include "fchk.h"
#include "bondwise_analysis.h"
#include "geometry_aid.h"
#include "crystal_energies.h"
#include "SALTED_equicomb.h"
#include "SALTED_predictor.h"
#ifdef NOSPHERA2_USE_GPU
#include "sf_gpu.h"
#endif
#include "nbo_run.h"
#include "nbo.h"
#include "eli_family.h"
#include "topology.h"
#include "citations.h"
#include "b2c.h"
#include "section_log.h"

namespace {
	ivec parse_rgbi_group_indices(const std::string &spec)
	{
		ivec indices;
		for (const auto &raw_part : split_string<std::string>(spec, ","))
		{
			const std::string part = trim(raw_part);
			if (part.empty())
				continue;

			const size_t range_pos = part.find('-', 1);
			if (range_pos == std::string::npos)
			{
				indices.push_back(std::stoi(part));
				continue;
			}

			const std::string first_text = trim(part.substr(0, range_pos));
			const std::string last_text = trim(part.substr(range_pos + 1));
			err_checkf(!first_text.empty() && !last_text.empty(),
				"Invalid RGBI atom index range: " + part, std::cout);

			const int first = std::stoi(first_text);
			const int last = std::stoi(last_text);
			const int step = first <= last ? 1 : -1;
			for (int atom_index = first;; atom_index += step)
			{
				indices.push_back(atom_index);
				if (atom_index == last)
					break;
			}
		}
		return indices;
	}

	void validate_rgbi_group_set(const ivec2 &group_set)
	{
		std::map<int, int> assigned_group;
		for (int group_index = 0; group_index < static_cast<int>(group_set.size()); ++group_index)
		{
			for (int atom_index : group_set[group_index])
			{
				auto [it, inserted] = assigned_group.emplace(atom_index, group_index);
				err_checkf(inserted,
					"Invalid RGBI groups: atom index " + std::to_string(atom_index) +
					" is assigned to both group " + std::to_string(it->second) +
					" and group " + std::to_string(group_index) + ".",
					std::cout);
			}
		}
	}
}

//What each -tune of the current run replaced, in order
static std::vector<std::pair<std::string, std::optional<std::string>>> run_tuned;
void tuning_end_run()
{
	for (auto it = run_tuned.rbegin(); it != run_tuned.rend(); ++it)
		set_tuning(it->first, it->second ? it->second->c_str() : nullptr);
	run_tuned.clear();
}

static std::vector<std::shared_ptr<BasisSet>> get_aux_basis(const int argc, const svec& arguments, const int i) {
	std::vector<std::shared_ptr<BasisSet>> aux_basis;
	int next_basis_set = i + 1;
	// Check if next argument is a valid basis set name or a new argument starting with "-"
	while (next_basis_set < argc && arguments[next_basis_set].find("-") != 0) {
		if (arguments[next_basis_set] == "auto_aux") {
			// an empty BasisSet is filled by gen_auto_aux (ORCA AutoAux, its own beta table);
			// element symbols after it restrict that to those elements, the rest comes from
			// the basis sets that follow: "-ri_fit auto_aux H def2-universal-jkfit"
			auto generated = std::make_shared<BasisSet>();
			ivec elements;
			while (next_basis_set + 1 < argc && constants::get_Z_from_label(arguments[next_basis_set + 1].c_str()) >= 0)
				elements.push_back(constants::get_Z_from_label(arguments[++next_basis_set].c_str()) + 1);
			generated->set_auto_aux_elements(elements);
			aux_basis.push_back(generated);
			next_basis_set++;
			continue;
		}
		err_checkf(BasisSetLibrary::check_basis_set_exists(arguments[next_basis_set]),
			"Basis set " + arguments[next_basis_set] + " not found in the library. Exiting.", std::cout);
		aux_basis.push_back(BasisSetLibrary::get_basis_set(arguments[next_basis_set]));
		next_basis_set++;
	}
	if (aux_basis.size() == 0) {
		std::cout << "No basis set specified. Falling back to automatic generation (auto_aux)!" << std::endl;
		aux_basis.push_back(std::make_shared<BasisSet>());
	}
	return aux_basis;
}

//Options the NBO handlers read from their own tokens, not through a digester. One list for both readers:
//the handlers refuse a -nbo_*/-nrt_* spelling not in it, digest_options does not report these as typos.
const std::set<std::string> &nbo_family_suboptions()
{
	static const std::set<std::string> known = {
		"-nbo", "-nbo_external", "-nbo_parse", "-nbo_native", "-nbo_keywords", "-nbo_exe", "-nbo_json",
		"-nbo_dir", "-nbo_47", "-nbo_ref", "-nbo_e2min", "-nbo_threads", "-nbo_keep47", "-nrt",
		"-nrt_e2", "-nrt_arrows", "-nrt_bond_scale", "-nrt_max", "-nrt_atoms",
		"-nrt_exhaustive", "-nrt_no_symmetry", "-nrt_no_components", "-nrt_no_ion", "-nrt_parent", "-nrt_e2_avg", "-nrt_e2_avg_min",
		"-nbo_centres", "-nbo_tm_pairs", "-nrt_target", "-nrt_screen", "-nrt_rotate", "-nrt_stall", "-nrt_mirror", "-nrt_no_mirror", "-nrt_tie", "-nrt_ref"};
	return known;
}

namespace {
	//No digester sees these tokens, so an unknown one would otherwise be dropped silently.
	void refuse_unknown_nbo_suboption(const std::vector<std::string> &arguments, int from,
									  const std::string &after)
	{
		for (int j = from; j < (int)arguments.size(); j++)
		{
			const std::string &o = arguments[j];
			if ((o.rfind("-nbo", 0) == 0 || o.rfind("-nrt", 0) == 0) &&
				nbo_family_suboptions().count(o) == 0)
				err_checkf(false, "Unknown option " + o + " after " + after +
									  ": it names no NBO/NRT option this build has. Check the "
									  "spelling, -h lists them; running on with it ignored would "
									  "give you a default you did not ask for",
						   std::cout);
		}
	}
}

//file, format and conversion options
//-merge a.tshb b.tshb [-cif x.cif ...]: the CIFs name every scatterer of the run; without them the merged file is partial
static void merge_tshb_args(const pathvec& files, const svec& arguments, const int argc, int at)
{
	pathvec cifs;
	if (at < argc && arguments[at] == "-cif")
		for (at++; at < argc && arguments[at].find("-") != 0; at++)
			cifs.push_back(arguments[at]);
	const std::vector<std::pair<std::string, atomID>> keys = cifs.empty() ? std::vector<std::pair<std::string, atomID>>{} : cif_scatterer_keys(cifs);
	merge_tshbs(files, "combined.tshb", cifs.empty() ? nullptr : &keys);
}

bool options::digest_io_options(const std::string &temp, int &i)
{
	using namespace std;
	const int argc = (int)arguments.size();
	if (temp == "-b")
		basis_set = arguments[i + 1];
	else if (temp == "-coef")
	{
		coef_file = arguments[i + 1];
		aux_basis = get_aux_basis(argc, arguments, i + 1);
		err_checkf(std::filesystem::exists(coef_file), "-coef: coefficient file doesn't exist: " + coef_file.string(), std::cout);
		SALTED = true;
	}
	else if (temp == "-cif")
	{
		cif = arguments[i + 1];
		err_checkf(std::filesystem::exists(cif), "-cif: CIF doesn't exist: " + cif.string(), std::cout);
	}
	else if (temp == "-convert_to_47") {
		err_checkf(argc >= i + 2, "Not enough arguments for -convert_to_47\nPlease provide at least stdout name!", std::cout);
		refuse_unknown_nbo_suboption(arguments, i + 2, "-convert_to_47");
		std::filesystem::path _wfn = arguments[i + 1];
		std::string keys;
		for (int j = i + 2; j + 1 < argc; j++)
			if (arguments[j] == "-nbo_keywords") keys = arguments[j + 1];
		WFN wavy(e_origin::NOT_YET_DEFINED);
		wavy.read_known_wavefunction_format(_wfn, std::cout, debug);
		citations::cite(citations::Method::NBOProgram, std::cout);
		wavy.write_nbo(_wfn.replace_extension(".47"), debug, &std::cout, keys);
		finished = true; return true;
	}
	else if (temp == "-nbo_external") {
		err_checkf(argc >= i + 2, "Not enough arguments for -nbo_external\nPlease provide a wavefunction!", std::cout);
		refuse_unknown_nbo_suboption(arguments, i + 2, "-nbo_external");
		NboRunOptions opt;
		opt.wavefunction = arguments[i + 1];
		opt.debug = debug;
		for (int j = i + 2; j + 1 < argc; j++) {
			if (arguments[j] == "-nbo_keywords") opt.keywords = arguments[j + 1];
			else if (arguments[j] == "-nbo_exe") opt.executable = arguments[j + 1];
			else if (arguments[j] == "-nbo_json") opt.json_out = arguments[j + 1];
			else if (arguments[j] == "-nbo_dir") opt.work_dir = arguments[j + 1];
		}
		finished = true;
		citations::cite(citations::Method::NBOProgram, std::cout);
		//finished is set, so only this check keeps a failed gennbo from exiting 0.
		const int nbo_rc = run_nbo(opt, std::cout);
		err_checkf(nbo_rc == 0, "-nbo_external: gennbo did not finish (exit code " + std::to_string(nbo_rc) +
									"), see the messages above",
				   std::cout);
		return true;
	}
	else if (temp == "-nbo_parse") {
		err_checkf(argc >= i + 2, "Not enough arguments for -nbo_parse\nPlease provide an NBO output file!", std::cout);
		refuse_unknown_nbo_suboption(arguments, i + 2, "-nbo_parse");
		std::filesystem::path out = arguments[i + 1];
		err_checkf(std::filesystem::exists(out), "NBO output doesn't exist: " + out.string(), std::cout);
		std::filesystem::path json = out.parent_path() / (out.stem().string() + ".nbo.json");
		for (int j = i + 2; j + 1 < argc; j++)
			if (arguments[j] == "-nbo_json") json = arguments[j + 1];
		NboResults r = parse_nbo_output(out);
		write_nbo_json(r, json);
		std::cout << "wrote " << json.string() << std::endl;
		finished = true;
		return true;
	}
	else if (temp == "-nbo" || temp == "-nbo_native") {
		//-nbo_native is the old name, kept for scripts and Olex2 builds that still send it.
		//Same JSON as -nbo_external, so tests/nbo_reference/compare_nbo.py compares the two.
		err_checkf(argc >= i + 2, "Not enough arguments for " + temp + "\nPlease provide a wavefunction!", std::cout);
		refuse_unknown_nbo_suboption(arguments, i + 2, temp);
		//-nbo ran gennbo until 29 Oct 2026; its options here mean an old command line, not a native run
		for (int j = i + 2; j < argc; j++)
			err_checkf(arguments[j] != "-nbo_exe" && arguments[j] != "-nbo_keywords" && arguments[j] != "-nbo_dir",
					   arguments[j] + " belongs to the external gennbo run, which is -nbo_external now; " + temp +
						   " is the in-house analysis and would ignore it",
					   std::cout);
		std::filesystem::path _wfn = arguments[i + 1];
		err_checkf(std::filesystem::exists(_wfn), "wavefunction doesn't exist: " + _wfn.string(), std::cout);
		NboOptions opt;
		opt.debug = debug;
		std::filesystem::path json = _wfn.parent_path() / (_wfn.stem().string() + ".native.nbo.json");
		for (int j = i + 2; j + 1 < argc; j++) {
			if (arguments[j] == "-nbo_json") json = arguments[j + 1];
			else if (arguments[j] == "-nbo_47") opt.file47 = arguments[j + 1];
			else if (arguments[j] == "-nbo_ref") opt.reference_scf = arguments[j + 1];
			else if (arguments[j] == "-nbo_e2min") opt.e2_threshold_kcal = std::stod(arguments[j + 1]);
			else if (arguments[j] == "-nbo_threads") opt.threads = std::stoi(arguments[j + 1]);
			//This branch returns before the main digester reads a later -cpus.
			else if (arguments[j] == "-cpus") {
				threads = std::stoi(arguments[j + 1]);
				MKL_Set_Num_Threads(threads);
#ifdef _OPENMP
				omp_set_num_threads(threads);
				omp_set_dynamic(0);
#endif
			}
			else if (arguments[j] == "-nrt") opt.nrt = true;
			else if (arguments[j] == "-nrt_e2") opt.nrt_e2_kcal = std::stod(arguments[j + 1]);
			else if (arguments[j] == "-nrt_e2_avg") {
				opt.nrt_e2_avg_percent = std::stod(arguments[j + 1]);
				opt.nrt = true;
			}
			else if (arguments[j] == "-nrt_e2_avg_min") {
				opt.nrt_e2_avg_min_percent = std::stod(arguments[j + 1]);
				opt.nrt = true;
			}
			else if (arguments[j] == "-nbo_tm_pairs") {
				opt.tm_pairs = std::stoi(arguments[j + 1]);
				err_checkf(opt.tm_pairs >= 0 && opt.tm_pairs <= 9, "-nbo_tm_pairs takes 0 (valence NAO count) to 9", std::cout);
			}
			else if (arguments[j] == "-nbo_centres") {
				opt.max_centres = std::stoi(arguments[j + 1]);
				err_checkf(opt.max_centres >= 2 && opt.max_centres <= 4, "-nbo_centres takes 2, 3 or 4", std::cout);
			}
			else if (arguments[j] == "-nrt_tie") opt.nrt_tie = std::stod(arguments[j + 1]);
			else if (arguments[j] == "-nrt_ref") opt.nrt_ref = std::stod(arguments[j + 1]);
			else if (arguments[j] == "-nrt_parent") {
				opt.nrt_parent = arguments[j + 1];
				opt.nrt = true;
				err_checkf(std::filesystem::exists(opt.nrt_parent),
						   "-nrt_parent: Lewis structure file doesn't exist: " + opt.nrt_parent.string(), std::cout);
			}
			else if (arguments[j] == "-nrt_arrows") opt.nrt_max_arrows = std::stoi(arguments[j + 1]);
			else if (arguments[j] == "-nrt_target") {
				opt.nrt_target = std::stod(arguments[j + 1]);
				opt.nrt = true;
			}
			else if (arguments[j] == "-nrt_screen")
			{
				opt.nrt_screen = arguments[j + 1] == "off" ? -1.0 : std::stod(arguments[j + 1]);
				if (arguments[j + 1] == "off") opt.nrt_screen_min = std::numeric_limits<int>::max();  //-nrt_ref walks unscreened too
			}
			else if (arguments[j] == "-nrt_stall") opt.nrt_stall = std::stoi(arguments[j + 1]);
			else if (arguments[j] == "-nrt_bond_scale")
				opt.nrt_bond_scale = std::stod(arguments[j + 1]);
			else if (arguments[j] == "-nrt_max") {
				opt.nrt_max_candidates = std::stoi(arguments[j + 1]);
				opt.nrt_max_set = true;
			}
			else if (arguments[j] == "-nrt_atoms") {
				//a comma separated 1-based atom list, NBO's "NRT <atoms>" subspace
				std::string s = arguments[j + 1];
				for (char& c : s) if (c == ',') c = ' ';
				std::istringstream is(s);
				int a = 0;
				while (is >> a) opt.nrt_subspace.push_back(a);
			}
		}
		for (int j = i + 2; j < argc; j++)
			if (arguments[j] == "-nrt") opt.nrt = true;
			else if (arguments[j] == "-nrt_exhaustive") opt.nrt_exhaustive = true;
			else if (arguments[j] == "-nrt_no_symmetry") opt.nrt_symmetry = false;
			else if (arguments[j] == "-nrt_no_components") opt.nrt_components = false;
			else if (arguments[j] == "-nrt_no_ion") opt.nrt_ion = false;
			else if (arguments[j] == "-nrt_rotate") opt.nrt_rotate = true;
			else if (arguments[j] == "-nrt_mirror") opt.nrt_mirror = true;
			else if (arguments[j] == "-nrt_no_mirror") opt.nrt_mirror = false;
			else if (arguments[j] == "-nbo_keep47") opt.keep_file47 = true;
		//-nbo_threads wins, then -cpus, then the OpenMP default.
		if (opt.threads <= 0 && threads > 0) opt.threads = threads;
		WFN wavy(e_origin::NOT_YET_DEFINED);
		wavy.read_known_wavefunction_format(_wfn, std::cout, debug);
		citations::cite(citations::Method::NAONPA, std::cout);
		citations::cite(citations::Method::NBO, std::cout);
		citations::cite(citations::Method::E2, std::cout);
		if (opt.nrt)
			citations::cite(citations::Method::NRT, std::cout);
		NboResults r = native_nbo(wavy, opt, std::cout);
		r.name = _wfn.stem().string();
		{
			using citations::Method;
			const section_log::section nbo_file(_wfn, "nbo", std::string("Natural bond orbitals: NPA, NBO, second-order E(2)") + (opt.nrt ? " and NRT" : ""),
				opt.nrt ? std::vector<Method>{ Method::NAONPA, Method::NBO, Method::E2, Method::NRT } : std::vector<Method>{ Method::NAONPA, Method::NBO, Method::E2 }, no_date);
			const section_log::tee nbo_tee("nbo");
			print_nbo(r, std::cout);
		}
		write_nbo_json(r, json);
		std::cout << "wrote " << json.string() << std::endl;
		finished = true;
		return true;
	}
	else if (temp == "-d")
		basis_set_path = arguments[i + 1];
	else if (temp == "-fchk")
		fchk = arguments[i + 1];
	else if (temp == "-gbw2wfn")
		gbw2wfn = true;
	else if (temp == "-hkl")
	{
		hkl = arguments[i + 1];
		err_checkf(std::filesystem::exists(hkl), "-hkl: hkl file doesn't exist: " + hkl.string(), std::cout);
	}
	else if (temp == "-wfn")
	{
		wfn = arguments[i + 1];
		err_checkf(std::filesystem::exists(wfn), "-wfn: wavefunction does not exist: " + wfn.string(), std::cout);
	}
	else if (temp == "-wfn_cif")
	{
		write_CIF = true;
	}
	else if (temp == "-xyz")
	{
		xyz_file = arguments[i + 1];
	}
	else if (temp == "-mask")
	{
		mask_file = arguments[i + 1];
	}
	else
		return false;
	return true;
};

//resources, accuracy and general run control
bool options::digest_run_options(const std::string &temp, int &i)
{
	using namespace std;
	const int argc = (int)arguments.size();
	if (temp == "-acc")
		accuracy = stoi(arguments[i + 1]);
	else if (temp == "-basin_grid")
		basin_grid = std::max(1, stoi(arguments[i + 1]));
	else if (temp == "-basin_cube")
		basin_cube = true;
	else if (temp == "-eli_basin_cube") {
		err_checkf(i + 2 < argc, "-eli_basin_cube needs two atoms, e.g. -eli_basin_cube Pt Zn", std::cout);
		eli_basin_atoms = { arguments[i + 1], arguments[i + 2] };
		i += 2;
	}
	else if (temp == "-eli_mask_cube") {
		err_checkf(i + 1 < argc && std::filesystem::exists(arguments[i + 1]), "-eli_mask_cube needs an existing cube", std::cout);
		eli_mask_cube = arguments[++i];
	}
	else if (temp == "-eli_only")
		eli_only = true;
	else if (temp == "-eli_contributions")
		eli_contributions = true;
	else if (temp == "-no_spin_eli")
		spin_eli = false;
	else if (temp == "-basin_analytic")
		; //accepted and ignored, it is the default
	else if (temp == "-basin_step")
		basin_step_scale_set(stod(arguments[i + 1]));
	else if (temp == "-basin_timing")
		basin_timing_set_enabled(true);
	else if (temp == "-basin_gpu")
		basin_gpu_set_mode(1);
	else if (temp == "-no_basin_gpu")
		basin_gpu_set_mode(0);
	else if (temp == "-no_beta_spheres")
		beta_spheres_set_enabled(false);
	else if (temp == "-no_adaptive_step")
		basin_adaptive_step_set_enabled(false);
	else if (temp == "-adaptive_step")
		basin_adaptive_step_set_enabled(true);
	else if (temp == "-beta_spheres")
		beta_spheres_set_enabled(true);
	else if (temp == "-Anion")
	{
		int n = 1;
		string store;
		if (debug)
			std::cout << "Looking for Anions!" << endl;
		while (i + n < argc && string(arguments[i + n]).find("-") != 0)
		{
			if (i + n - 1 > arguments.size())
				break;
			store = arguments[i + n];
			svec Z = split_string<string>(store, " ");
			for (int r = 0; r < Z.size(); r++)
			{
				if (debug)
					std::cout << Z[r] << endl;
				Anions.push_back(Z[r]);
			}
			n++;
		}
	}
	else if (temp == "-Cation")
	{
		int n = 1;
		string store;
		if (debug)
			std::cout << "Looking for Cations!" << endl;
		while (i + n < argc && string(arguments[i + n]).find("-") != 0)
		{
			if (i + n - 1 > arguments.size())
				break;
			store = arguments[i + n];
			svec Z = split_string<string>(store, " ");
			for (int r = 0; r < Z.size(); r++)
			{
				if (debug)
					std::cout << Z[r] << endl;
				Cations.push_back(Z[r]);
			}
			n++;
		}
	}
	else if (temp == "-charge")
	{
		charge = stoi(arguments[i + 1]);
	}
	else if (temp == "-cpus")
	{
		threads = stoi(arguments[i + 1]);
		MKL_Set_Num_Threads(threads);
#ifdef _OPENMP
		omp_set_num_threads(threads);
		omp_set_dynamic(0);
#endif
	}
	else if (temp == "-dmin")
		dmin = stod(arguments[i + 1]);
	else if (temp == "-e_field")
		efield = stod(arguments[i + 1]);
	else if (temp == "-ECP" || temp == "-ecp" || temp == "-Ecp")
	{
		ECP = true;
		if (argc >= i + 2 && string(arguments[i + 1]).find("-") != 0)
		{
			ECP_mode = stoi(arguments[i + 1]);
		}
	}
	else if (temp == "-group")
	{
		int n = 1;
		while (i + n < argc)
		{
			const string& group_arg = arguments[i + n];
			// Olex2 can emit a bare -group followed by an empty argument.
			// Do not index or parse an empty string: doing so corrupts the
			// CRT state and later manifests as a stack-buffer-overrun.
			if (group_arg.empty() || group_arg.find("-") != string::npos)
				break;
			int group;
			if (group_arg[0] == '+')
				group = -stoi(group_arg);
			else
				group = stoi(group_arg);
			groups[0].push_back(group), n++;
		}
		i += n - 1;
	}
	else if (temp == "-hkl_min_max")
	{
		int h_min(stoi(arguments[i + 1]));
		int h_max(stoi(arguments[i + 2]));
		int k_min(stoi(arguments[i + 3]));
		int k_max(stoi(arguments[i + 4]));
		int l_min(stoi(arguments[i + 5]));
		int l_max(stoi(arguments[i + 6]));
		hkl_min_max = { {h_min, h_max}, {k_min, k_max}, {l_min, l_max} };
	}
	else if (temp == "-hkl_scale")
	{
		hkl_scale = stoi(arguments[i + 1]);
		err_checkf(hkl_scale >= 1, "-hkl_scale needs a positive integer", std::cout);
	}
	else if (temp == "-mem")
	{
		mem = stod(arguments[i + 1]); // In MB
		mem_given = true;
		vec a;
		size_t vec_max_size = a.max_size();
		double doubel_max_size = static_cast<double>(vec_max_size * sizeof(double)) * 1e-6;
		if (mem > doubel_max_size)
		{
			std::cout << "Max memory set to " << mem << " MB, which is larger than the maximum allowed size of " << doubel_max_size << " MB. Setting max memory to " << 50000 << " MB." << endl;
			mem = 50000.0;
		}
	}
	else if (temp == "-method")
		method = arguments[i + 1];
	else if (temp == "-mult")
		mult = stoi(arguments[i + 1]);
	else if (temp == "-no_date")
		no_date = constants::hide_gpu_notes = constants::hide_timings = true;
	else if (temp == "-out")
		; //run_app opened the log before options are digested; listed so it is not "unknown"
	else if (temp == "-no_date_but_gpu")
	{
		no_date = constants::hide_timings = true;
		constants::hide_gpu_notes = false;
	}
	else if (temp == "-pbc")
		pbc = stoi(arguments[i + 1]);
	else if (temp == "-profiling" || temp == "-profile")
	{
		profiling = true;
		if (i + 1 < argc && arguments[i + 1].find("-") != 0)
		{
			profiling_tests_root = arguments[i + 1];
		}
	}
	else if (temp == "-radius")
		properties.radius = stod(arguments[i + 1]);
	else if (temp == "-resolution")
		properties.resolution = stod(arguments[i + 1]);
	else if (temp == "-refine")
		argc > i + 1 ? properties.integral_accuracy = stod(arguments[i + 1]) : properties.integral_accuracy = 0.1;
	else if (temp.find("-rkpts") < 1)
		read_k_pts = true;
	else if (temp == "-skpts")
		save_k_pts = true;
	else if (temp == "-twin")
	{
		twin_law.resize(twin_law.size() + 1);
		twin_law.back().resize(9);
		for (int twl = 0; twl < 9; twl++)
			twin_law.back()[twl] = stod(arguments[i + 1 + twl]);
		if (debug)
		{
			std::cout << "twin_law: ";
			for (int twl = 0; twl < 9; twl++)
				std::cout << setw(7) << setprecision(2) << twin_law.back()[twl];
			std::cout << endl;
		}
		i += 9;
	}
	else
		return false;
	return true;
};

//partitioning schemes, scattering factors and tsc tables
bool options::digest_partition_options(const std::string &temp, int &i)
{
	using namespace std;
	const int argc = (int)arguments.size();
	if (temp == "-atom_sfac")
	{
		wfn = arguments[i + 1];
		wfn2 = arguments[i + 2];
		err_checkf(std::filesystem::exists(wfn), "-atom_sfac: wavefunction doesn't exist: " + wfn.string(), std::cout);
		WFN wavy(e_origin::NOT_YET_DEFINED);
		wavy.read_known_wavefunction_format(wfn, std::cout, debug);
		wavy.delete_unoccupied_MOs();

		WFN wavy2(e_origin::NOT_YET_DEFINED);
		wavy2.read_known_wavefunction_format(wfn2, std::cout, debug);
		wavy2.delete_unoccupied_MOs();

		bvec needs_grid(wavy.get_ncen(), false);
		needs_grid[0] = true;
		GridConfiguration conf;
		conf.partition_type = PartitionType::Hirshfeld;
		conf.accuracy = 4;
		GridManager grid(conf);
		vec2 d1, d2, d3, dens;
		vec2 d1_2, d2_2, d3_2, dens_2;

		cell unit_cell(10.0, 10.0, 10.0, 90, 90, 90);
		ivec asym_atom_list(1, 0);
		// Setup grids for the molecule
		auto grid2 = grid;
		grid.setup3DGridsForMolecule(wavy, asym_atom_list, needs_grid, unit_cell);
		grid.getDensityVectors(wavy, asym_atom_list, d1, d2, d3, dens);
		grid2.setup3DGridsForMolecule(wavy2, asym_atom_list, needs_grid, unit_cell);
		grid2.getDensityVectors(wavy2, asym_atom_list, d1_2, d2_2, d3_2, dens_2);

		for (int j = 0; j < d1.size(); j++)
		{
			for (int p = 0; p < d1[0].size(); p++)
			{
				dens[j][p] -= dens_2[j][p];
			}
		}

		// Calculate partitioned charges
		PartitionResults results = grid.calculatePartitionedCharges(wavy, unit_cell);

		grid.printChargeTable({ "Fe" }, wavy, asym_atom_list, std::cout, results);

		const int points = grid.getTotalGridPoints();
		for (double k = 0.001; k < 10; k += 0.002) {
			cdouble res = calc_spherically_averaged_at_k(d1, d2, d3, dens, k);
			std::cout << "k: " << k << " sfac: " << setprecision(9) << setw(16) << scientific << res.real() << " " << setprecision(9) << setw(16) << scientific << res.imag() << endl;
		}
		finished = true; return true;
	}
	else if (temp == "-becke" || temp == "-BECKE" || temp == "-Becke")
		partition_type = PartitionType::Becke;
	else if (temp == "-cmtc")
	{
		cif_based_combined_tsc_calc = true;
		int n = 1;
		string delimiter = ",";
		groups.pop_back();
		while (i + n < argc && string(arguments[i + n]).find("-") > 0)
		{
			combined_tsc_calc_files.push_back(arguments[i + n]);
			n++;
			combined_tsc_calc_cifs.push_back(arguments[i + n]);
			n++;
			const string _temp = arguments[i + n];
			if (_temp.find(delimiter) == string::npos) {
				if (debug)
					std::cout << "--Delimiter not found, using ." << endl;
				delimiter = ".";
			}
			groups.push_back(split_string<int>(_temp, delimiter));
			if (debug)
			{
				std::cout << "--Group: " << _temp << endl << "--";
				for (int run = 0; run < groups[groups.size() - 1].size(); run++)
					std::cout << groups[groups.size() - 1][run] << " ";
				std::cout << endl;
			}
			n++;
		}
	}
	else if (temp == "-def" || temp == "-DEF")
		properties.def = true;
	else if (temp == "-ED")
		electron_diffraction = true;
	else if (temp == "-embis" || temp == "-EMBIS")
	{
		partition_type = PartitionType::EMBIS;
	}
	else if (temp == "-HDEF" || temp == "-hdef")
		properties.hdef = true;
	else if (temp == "-hirsh")
		properties.hirsh = true, properties.hirsh_number = stoi(arguments[i + 1]);
	else if (temp == "-tsc_block")
	{
		tsc_block_size = static_cast<size_t>(std::stoll(arguments[i + 1]));
		tsc_block_given = true;
	}
	else if (temp == "-IAM")
		iam_switch = true;
	else if (temp == "-mbis" || temp == "-MBIS")
	{
		partition_type = PartitionType::MBIS;
	}
	else if (temp == "-merge")
	{
		pathvec filenames;
		int n = 1;
		while (i + n < argc && string(arguments[i + n]).find("-") > 0)
		{
			filenames.push_back(arguments[i + n]);
			n++;
		}
		//.tshb files carry densities, no reflections: merged by key alone, first wins, as the tables merge their columns
		const auto n_tshb = std::count_if(filenames.begin(), filenames.end(), [](const std::filesystem::path& f) { return f.extension() == ".tshb"; });
		err_checkf(n_tshb == 0 || n_tshb == (long)filenames.size(), "-merge takes either tables or .tshb files, not both", std::cout);
		if (n_tshb > 0)
			merge_tshb_args(filenames, arguments, argc, i + n);
		else
			merge_tscs("combine", filenames, old_tsc);
		finished = true; return true;
	}
	else if (temp == "-merge_nocheck")
	{
		pathvec filenames;
		int n = 1;
		while (i + n < argc && string(arguments[i + n]).find("-") > 0)
		{
			filenames.push_back(arguments[i + n]);
			n++;
		}
		//.tshb files carry densities, no reflections: merged by key alone, first wins, as the tables merge their columns
		const auto n_tshb = std::count_if(filenames.begin(), filenames.end(), [](const std::filesystem::path& f) { return f.extension() == ".tshb"; });
		err_checkf(n_tshb == 0 || n_tshb == (long)filenames.size(), "-merge takes either tables or .tshb files, not both", std::cout);
		if (n_tshb > 0)
			merge_tshb_args(filenames, arguments, argc, i + n);
		else
			merge_tscs_without_checks("combine", filenames, old_tsc);
		finished = true; return true;
	}
	else if (temp == "-salted_fold")
	{
		err_checkf(i + 2 < argc, "-salted_fold needs the model and the file to write", std::cout);
		fold_salted_file(arguments[i + 1], arguments[i + 2]);
		finished = true; return true;
	}
	else if (temp == "-salted_ecp")
	{
		err_checkf(i + 3 < argc, "-salted_ecp needs the model, the file to write and the ECP mode", std::cout);
		write_salted_ecp(arguments[i + 1], arguments[i + 2], stoi(arguments[i + 3]));
		finished = true; return true;
	}
	else if (temp == "-mtc")
	{
		combined_tsc_calc = true;
		int n = 1;
		string delimiter = ",";
		groups.pop_back();
		while (i + n < argc && string(arguments[i + n]).find("-") > 0)
		{
			combined_tsc_calc_files.push_back(arguments[i + n]);
			n++;
			const string _temp = arguments[i + n];
			if (_temp.find(delimiter) == string::npos) {
				if (debug)
					std::cout << "--Delimiter not found, using ." << endl;
				delimiter = ".";
			}
			groups.push_back(split_string<int>(_temp, delimiter));
			if (debug)
			{
				std::cout << "--Group: " << _temp << endl << "--";
				for (int run = 0; run < groups[groups.size() - 1].size(); run++)
					std::cout << groups[groups.size() - 1][run] << " ";
				std::cout << endl;
			}
			n++;
		}
	}
	else if (temp == "-mtc_mult")
	{
		int n = 1;
		while (i + n < argc && string(arguments[i + n]).find("-") > 0)
		{
			combined_tsc_calc_mult.push_back(stoi(arguments[i + n]));
			n++;
		}
	}
	else if (temp == "-mtc_charge")
	{
		int n = 1;
		while (i + n < argc && string(arguments[i + n]).find("-") > 0)
		{
			if (arguments[i + n][0] == 'n')
				combined_tsc_calc_charge.push_back(-stoi(arguments[i + n].substr(1)));
			else
				combined_tsc_calc_charge.push_back(stoi(arguments[i + n]));
			n++;
		}
	}
	else if (temp == "-mtc_ECP")
	{
		int m = 1;
		while (i + m < argc && string(arguments[i + m]).find("-") > 0)
		{
			combined_tsc_calc_ECP.push_back(stoi(arguments[i + m]));
			m++;
		}
	}
	else if (temp == "-SALTED" || temp == "-salted")
	{
		SALTED = true;
		// Several models may follow; the first also sets salted_model_dir for the single-model paths.
		int next = i + 1;
		while (next < argc && arguments[next].find("-") != 0)
			salted_model_dirs.push_back(arguments[next++]);
		err_checkf(!salted_model_dirs.empty(), "No SALTED model given after -SALTED.", std::cout);
		salted_model_dir = salted_model_dirs[0];
	}
	else if (temp == "-SALTED_COEFS" || temp == "-salted_coefs")
	{
		salted_model_dir = arguments[i + 1];

		//Check that wfn is not empty
		err_checkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -SALTED_COEFS to specify a wavefunction.", std::cout);

#ifdef NOSPHERA2_USE_GPU
		//This runs inside the parser, before run_app_impl sets the GPU globals from opt, so the
		//descriptor path was whatever the previous run in the process left on. Needs -no_gpu_salted before -SALTED_COEFS.
		equicomb_set_gpu(use_gpu && gpu_salted);
		//Nor does it get run_app's warm-up: started here, context creation overlaps reading the
		//wavefunction and the model, and predict() joins it
		if (use_gpu && gpu_salted)
			sf_gpu_warmup_start();
#endif
		// The stages around the prediction, which reports its own
		auto mark = std::chrono::steady_clock::now();
		auto lap = [&mark](const char* what) {
			const auto now = std::chrono::steady_clock::now();
			if (ProgressBar::report_counts)
				std::cout << "[stages] " << what << " " << std::chrono::duration<double>(now - mark).count() << " s" << std::endl;
			mark = now;
		};
		WFN wavy(wfn);
		lap("read structure");
		SALTEDPredictor SP(wavy, *this);
		lap("predictor setup");
		filesystem::path salted_model_path = SP.get_salted_filename();
		log_file << "Using " << salted_model_path << " for the prediction" << endl;
		if (!SP.basis_set_loaded()) {
			const string df_basis_name = SP.get_dfbasis_name();
			std::shared_ptr<BasisSet> _aux_basis = BasisSetLibrary::get_basis_set(df_basis_name);
			load_basis_into_WFN(SP.wavy, _aux_basis);
		}
		lap("basis");
		vec coefs = SP.gen_SALTED_densities();
		lap("prediction");

		const aux_density_table t(SP.wavy.get_atoms());

		DensityFitting::analyze_density_fit_quality(
			coefs,
			SP.wavy,
			t,
			vec());
		lap("fit quality");

		npy::npy_data<double> np_coeffs;
		np_coeffs.data = coefs;
		np_coeffs.fortran_order = false;
		np_coeffs.shape = { static_cast<unsigned long>(coefs.size()) };
		npy::write_npy("SALTED_COEFS.npy", np_coeffs);
		lap("write coefficients");
	}
	else if (temp == "-SALTED_Training") {
		err_checkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -SALTED_Training to specify a wavefunction.", std::cout);
		err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEFORE -SALTED_Training to specify an auxiliary basis set.", std::cout);

		WFN wavy(wfn);
		WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);

		create_SALTED_training_data(wavy, wavy_aux, *this);
		finished = true; return true;
	}
	else if (temp == "-sfac_diffuse")
	{
		sfac_diffuse[0] = fromString<double>(arguments[i + 1]);
		sfac_diffuse[1] = fromString<double>(arguments[i + 2]);
		sfac_diffuse[2] = fromString<double>(arguments[i + 3]);
		cif = arguments[i + 4];
		wfn = arguments[i + 5];
		dmin = fromString<double>(arguments[i + 6]);
		calc_sfac_diffuse(*this, std::cout);
		finished = true; return true;
	}
	else if (temp == "-old_tsc")
	{
		old_tsc = true;
	}
	else if (temp == "-tshb")
		write_tshb = true;
	else if (temp == "-tfvc" || temp == "-TFVC")
	{
		partition_type = PartitionType::TFVC;
	}
	else if (temp == "-tscb")
	{
		std::filesystem::path name = arguments[i + 1];
		//the title is the file stem, as the -merge conversion does it; it used to be the placeholder "test"
		if (name.extension() == ".tscb")
		{
			tsc_block<int, cdouble> blocky = read_tsc_table(name);
			name.replace_extension(".tsc");
			blocky.write_tsc_file(name, name);
		}
		else if (name.extension() == ".tsc")
		{
			tsc_block<int, cdouble> blocky = read_tsc_table(name);
			name.replace_extension(".tscb");
			blocky.write_tscb_file(name, name);
		}
		else
			err_checkf(false, "-tscb takes a .tsc or .tscb file, not " + name.string(), std::cout);
		finished = true; return true;
	}
	else if (temp == "-tsc_labels")
	{
		const bool has_table_and_cif =
			i + 2 < argc &&
			arguments[i + 1].find('-') != 0 &&
			arguments[i + 2].find('-') != 0;

		if (has_table_and_cif)
		{
			const std::filesystem::path table = arguments.at(i + 1);
			const std::filesystem::path cif_file = arguments.at(i + 2);
			std::filesystem::path output = table;
			output.replace_extension(".labels.tsc");
			if (i + 3 < argc && arguments[i + 3].find('-') != 0)
				output = arguments[i + 3];
			err_checkf(convert_tsc_ids_to_labels(table, cif_file, output, std::cout), "Could not convert the tsc labels", std::cout);
			finished = true; return true;
		}

		label_tsc_output = true;
	}
	else if (temp == "-occ")
	{
		occ = arguments[i + 1];
		err_checkf(std::filesystem::exists(occ), "-occ: OCC input doesn't exist: " + occ, std::cout);

	}
	else
		return false;
	return true;
};

//cube and property evaluation
bool options::digest_property_options(const std::string &temp, int &i)
{
	using namespace std;
	const int argc = (int)arguments.size();
	if (temp == "-all_charges")
		all_charges = true;
	else if (temp == "-laplacian_bonds")
	{
		wfn = arguments[i + 1];
		bondwise_laplacian_plots(wfn);
		finished = true; return true;
	}
	else if (temp == "-laplacian_bond_paths")
	{
		wfn = arguments[i + 1];
		bondwise_laplacian_plots(wfn, true);
		finished = true; return true;
	}
	else if (temp == "-atom_dens")
	{
		wfn = arguments[i + 1];
		err_checkf(std::filesystem::exists(wfn), "-atom_dens: wavefunction doesn't exist: " + wfn.string(), std::cout);
		ivec val_MOs;
		ivec val_MOs_beta;
		if (argc > i + 3 && arguments[i + 2][0] != '-')
		{
			val_MOs = split_string<int>(arguments[i + 2], ",");
			std::cout << "Alpha MOs to keep: ";
			for (int j = 0; j < val_MOs.size(); j++)
				std::cout << val_MOs[j] << " ";
			std::cout << endl;
			val_MOs_beta = split_string<int>(arguments[i + 3], ",");
			std::cout << "Beta MOs to keep: ";
			for (int j = 0; j < val_MOs_beta.size(); j++)
				std::cout << val_MOs_beta[j] << " ";
			std::cout << endl;
		}
		spherically_averaged_density(*this, val_MOs, val_MOs_beta);
		finished = true; return true;
	}
	else if (temp == "-combine_mos")
	{
		combine_mo.push_back(arguments[i + 1]);
		combine_mo.push_back(arguments[i + 2]);
		do_combine_mo(*this);
		finished = true; return true;
	}
	else if (temp == "-cmos1")
	{
		int j = 1;
		while (i + j < argc && arguments[i + j].find("-") >= 1)
		{
			cmo1.push_back(stoi(arguments[i + j]));
			j++;
		}
	}
	else if (temp == "-cmos2")
	{
		int j = 1;
		while (i + j < argc && arguments[i + j].find("-") >= 1)
		{
			cmo2.push_back(stoi(arguments[i + j]));
			j++;
		}
	}
	else if (temp == "-density_difference")
	{
		wfn2 = arguments[i + 1];
	}
	else if (temp == "-dipole_moments")
	{
		dipole_moments(*this);
		finished = true; return true;
	}
	// Visualize the specified orbital using spherical harmonics.
	// Call as -draw_orbits lambda,m,resolution,radius
	// Where resolution and radius are optional
	else if (temp == "-draw_orbits")
	{
		vec opts = split_string<double>(arguments[i + 1], ",");
		int l = static_cast<int>(opts[0]);
		int m = static_cast<int>(opts[1]);
		properties.resolution = 0.025;
		properties.radius = 3.5;
		if (opts.size() >= 3)
		{
			properties.resolution = opts[2];
		}
		if (opts.size() == 4)
		{
			properties.radius = opts[3];
		}

		draw_orbital(l, m, properties.resolution, properties.radius);
		finished = true; return true;
	}
	else if (temp == "-eli")
		properties.eli = true;
	else if (temp == "-eli_spin")
		properties.eli_spin = true;
	//The basins use the analytic field; resolution and radius only shape the cube built for -basin_cube
	//and a fitted density
	else if (temp == "-eli_analysis" || temp == "-qtaim") {
		err_checkf(argc >= i + 2, "Not enough arguments for " + temp + "\nPlease provide a wfn!", std::cout);
		wfn = arguments[i + 1];
		err_checkf(std::filesystem::exists(wfn), temp + ": wavefunction does not exist: " + wfn.string(), std::cout);
		properties.resolution = 0.05;
		properties.radius = 3.0;
		i += 1;
		if (argc >= i + 3 && arguments[i + 1][0] != '-' && arguments[i + 2][0] != '-') {
			properties.resolution = stod(arguments[i + 1]);
			properties.radius = stod(arguments[i + 2]);
			i += 2;
		}
		eli_analysis_run = true;
		qtaim_only = temp == "-qtaim";
	}
	else if (temp == "-fba") {
		err_checkf(argc >= i + 2, "Not enough arguments for -fba\nPlease provide a wfn!", std::cout);
		wfn = arguments[i + 1];
		err_checkf(std::filesystem::exists(wfn), "-fba: wavefunction does not exist: " + wfn.string(), std::cout);
		properties.resolution = 0.05;
		properties.radius = 3.0;
		i += 1;
		fba = true;
		eli_contributions = true;
	}
	//Runs inline rather than setting a flag: it has no cube or basin output to schedule
	else if (temp == "-eli_family") {
		err_checkf(argc >= i + 2, "Not enough arguments for -eli_family\nPlease provide at least a wfn!", std::cout);
		err_checkf(std::filesystem::exists(arguments[i + 1]), "-eli_family: wavefunction does not exist: " + arguments[i + 1], std::cout);
		const bool has_points = argc >= i + 3 && arguments[i + 2].size() > 0 && arguments[i + 2][0] != '-';
		eli_family::report(arguments[i + 1], has_points ? std::filesystem::path(arguments[i + 2]) : std::filesystem::path());
		finished = true; return true;
	}
	else if (temp == "-topology") {
		err_checkf(argc >= i + 2, "Not enough arguments for -topology\nPlease provide a wfn!", std::cout);
		err_checkf(std::filesystem::exists(arguments[i + 1]), "-topology: wavefunction does not exist: " + arguments[i + 1], std::cout);
		{
			const bool complete = topology::report(arguments[i + 1], std::cout);
			finished = true;
			//An INCOMPLETE set must exit non-zero. Some incomplete sets balance the Poincare-Hopf sum and
			//fail only the bond-graph rank, so the message names both.
			err_checkf(complete, "-topology: the set of critical points printed above is INCOMPLETE "
								 "for " + arguments[i + 1] + " - either the Poincare-Hopf sum does "
								 "not close, or it closes while contradicting the bond graph. See the "
								 "diagnosis printed with the verdict for which class is short or "
								 "spurious",
					   std::cout);
		}
		return true;
	}
	else if (temp == "-qtaim_eli") {
		// Cube-files mode:  -qtaim_eli <rho.cube> <eli.cube> <atoms_csv> [<bg_value>]
		// WFN mode:         -qtaim_eli <wfn_file> <atoms_csv> [<resolution>] [<radius>] [<bg_value>]
		err_checkf(i + 2 < argc,
			"Usage:\n"
			"  -qtaim_eli <rho.cube> <eli.cube> <atoms_csv> [bg_value]\n"
			"  -qtaim_eli <wfn_file> <atoms_csv> [resolution] [radius] [bg_value]\n"
			"atoms_csv: comma-separated 0-based atom indices, e.g. 0,3,7\n"
			"bg_value:  value assigned to non-selected voxels (default 0)",
			std::cout);

		std::filesystem::path arg1 = arguments[i + 1];
		err_checkf(std::filesystem::exists(arg1), "-qtaim_eli: input does not exist: " + arg1.string(), std::cout);
		const bool cube_mode = (arg1.extension() == ".cube");

		std::filesystem::path rho_path, eli_path_arg;
		std::string atoms_csv;
		double bg_val = 0.0;

		if (cube_mode) {
			err_checkf(i + 3 < argc,
				"Cube mode requires: -qtaim_eli <rho.cube> <eli.cube> <atoms_csv>", std::cout);
			rho_path     = arg1;
			eli_path_arg = arguments[i + 2];
			atoms_csv    = arguments[i + 3];
			if (i + 4 < argc) bg_val = stod(arguments[i + 4]);
		} else {
			rho_path  = arg1;  // wfn path
			atoms_csv = arguments[i + 2];
			if (i + 3 < argc && !std::filesystem::path(arguments[i + 3]).has_extension())
				; // next arg looks like atoms_csv already consumed; nothing extra
			if (i + 3 < argc) {
				try { properties.resolution = stod(arguments[i + 3]); }
				catch (...) { /* optional — might be bg_val */ }
			}
			if (i + 4 < argc) {
				try { properties.radius = stod(arguments[i + 4]); }
				catch (...) {}
			}
			if (i + 5 < argc) {
				try { bg_val = stod(arguments[i + 5]); }
				catch (...) {}
			}
		}

		// Parse comma-separated 0-based atom indices
		std::vector<int> indices;
		{
			std::stringstream ss(atoms_csv);
			std::string tok;
			while (std::getline(ss, tok, ',')) {
				std::string t = trim(tok);
				if (!t.empty()) indices.push_back(std::stoi(t));
			}
		}
		err_checkf(!indices.empty(), "No atom indices parsed from: " + atoms_csv, std::cout);

		run_QTAIM_ELI_mask(rho_path, eli_path_arg, indices, bg_val, *this, log_file);
		finished = true; return true;
	}
	else if (temp == "-elf")
		properties.elf = true;
	else if (temp == "-esp")
		properties.esp = true;
	else if (temp == "-esp_isosurface")
	{
		properties.esp_isosurface = 0.002;
		if (argc > i + 1 && !arguments[i + 1].empty() && arguments[i + 1][0] != '-')
			properties.esp_isosurface = stod(arguments[i + 1]);
	}
	else if (temp == "-no_gpu")
		use_gpu = false;
	else if (temp == "-gpu_fp64")
		gpu_fp64 = true;
	else if (temp == "-gpu_fp32")
		gpu_fp32 = true;
	else if (temp == "-gflops")
	{
		track_gflops = true;
		throughput::set_enabled(true);
	}
	else if (temp == "-gpu_itensor")
		gpu_itensor = true;
	else if (temp == "-cpu_itensor_fp32")
		cpu_itensor_fp32 = true;
	else if (temp == "-no_cpu_itensor_fp32")
		cpu_itensor_fp32 = false;
	else if (temp == "-xcw_extrapolate")
		xcw_extrapolate = true;
	else if (temp == "-no_xcw_extrapolate")
		xcw_extrapolate = false;
	else if (temp == "-xcw_incremental")
		xcw_incremental = true;
	else if (temp == "-print_F")
		print_F = true;
	else if (temp == "-no_xcw_incremental")
		xcw_incremental = false;
	else if (temp == "-xcw_int_precision")
		xcw_int_precision = stod(arguments[i + 1]);
	else if (temp == "-xcw_strong_cutoff")
		xcw_settings.xcw_strong_cutoff = stod(arguments[i + 1]);
	else if (temp == "-no_gpu_itensor")
		gpu_itensor = false;
	else if (temp == "-gpu_cublas")
		gpu_cublas = true;
	else if (temp == "-no_gpu_cublas")
		gpu_cublas = false;
	else if (temp == "-gpu_salted")
		gpu_salted = true;
	else if (temp == "-no_gpu_salted")
		gpu_salted = false;
	else if (temp == "-gpu_grid")
		gpu_grid = true;
	else if (temp == "-no_gpu_grid")
		gpu_grid = false;
	else if (temp == "-gpu_density")
		gpu_density = true;
	else if (temp == "-no_gpu_density")
		gpu_density = false;
	else if (temp == "-gpu_blas")
		gpu_blas = true;
	else if (temp == "-fukui" || temp == "-Fukui")
		properties.fukui = true;
	else if (temp == "-fukui_analysis")
	{
		// Only records the request; the analysis itself runs from
		// run_app_impl.
		//
		// Doing the work here instead - the pattern -dipole_moments,
		// -laplacian_bonds and friends use - would put the whole table in
		// NoSpherA2.log and nothing on the terminal, because run_app_impl
		// has already pointed std::cout at that log file by the time
		// digest_options() runs. (Verified: -dipole_moments prints nothing
		// to stdout and 2.9 kB to the log.) That is tolerable for a
		// side-effect command and useless for one whose entire output is
		// meant to be read, so this one is dispatched later, where the
		// console buffer can be restored first.
		fukui_analysis_run = true;
		// Wavefunction may be given inline (-fukui_analysis mol.gbw) or with
		// the usual -wfn. The inline form is the point of the flag, so it
		// wins, but only if the next token is not itself an option.
		if (i + 1 < argc && arguments[i + 1][0] != '-')
		{
			wfn = arguments[i + 1];
			i++;
		}
	}
	else if (temp == "-eqc")
	{
		// Only records the request; eqc::run is called from run_app_impl.
		eqc = true;
		if (i + 1 < argc && arguments[i + 1][0] != '-')
			wfn = arguments[++i];
	}
	else if (temp == "-eqc_frag")
	{
		// Triples <atoms> <charge> <multiplicity> for as long as the next token starts with a
		// digit; a charge may be negative, which is why only the atoms token is tested
		while (i + 3 < argc && std::isdigit(static_cast<unsigned char>(arguments[i + 1][0])))
		{
			eqc_fragment f;
			f.atoms = parse_rgbi_group_indices(arguments[i + 1]);
			f.charge = std::stoi(arguments[i + 2]);
			f.mult = std::stoi(arguments[i + 3]);
			err_checkf(f.mult >= 1, "-eqc_frag: multiplicity must be 1 or more", std::cout);
			eqc_frags.push_back(f);
			i += 3;
		}
		err_checkf(!eqc_frags.empty(), "-eqc_frag wants <atoms> <charge> <multiplicity> triples", std::cout);
	}
	else if (temp == "-eqc_wfn")
	{
		eqc = true;
		while (i + 1 < argc && arguments[i + 1][0] != '-')
			eqc_wfns.emplace_back(arguments[++i]);
	}
	else if (temp == "-eqc_method")
		eqc_method = arguments[++i];
	else if (temp == "-eqc_basis")
		eqc_basis = arguments[++i];
	else if (temp == "-eqc_cold")
		eqc_cold = true;
	else if (temp == "-ewal_sum")
	{
		// bool read, WFN& wave, std::ostream& file,
		WFN *temp_w = new WFN(e_origin::cub);
		cube residual(arguments[i + 1], true, *temp_w, std::cout);
		if (argc >= i + 3)
		{
			int k_max = stoi(arguments[i + 2]);
			if (argc >= i + 4)
				residual.ewald_sum(k_max, stod(arguments[i + 3]));
			else
				residual.ewald_sum(k_max);
		}
		else
			residual.ewald_sum();
		delete (temp_w);
		finished = true; return true;
	}
	else if (temp == "-fractal")
		fract = true, fract_name = arguments[i + 1];
	else if (temp == "-get_g")
		get_g = true;
	else if (temp == "-hirshfeld_surface")
	{
		hirshfeld_surface = arguments[i + 1];
		hirshfeld_surface2 = arguments[i + 2];
	}
	else if (temp == "-lap")
		properties.lap = true;
	else if (temp == "-MO")
	{
		if (string(arguments[i + 1]) != "all")
			properties.MO_numbers.push_back(stoi(arguments[i + 1]));
		else
			properties.all_mos = true;
	}
	else if (temp == "-polarizabilities")
	{
		pol_wfns = { arguments[i + 1],
					arguments[i + 2],
					arguments[i + 3],
					arguments[i + 4],
					arguments[i + 5],
					arguments[i + 6],
					arguments[i + 7] };
	}
	else if (temp == "-QCT" || temp == "-qct")
		qct = true;
	else if (temp == "-promol_nci")
	{
		promol_nci = true;
		//Every following .xyz is a fragment; the optional numeric cutoffs come after the last one
		int n_xyz = 0;
		auto is_xyz = [](const std::string &a) {
			std::string ext = std::filesystem::path(a).extension().string();
			std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });
			return ext == ".xyz";
		};
		while (i + 1 + n_xyz < argc && is_xyz(arguments[i + 1 + n_xyz]))
		{
			promol_nci_xyz.push_back(arguments[i + 1 + n_xyz]);
			err_checkf(std::filesystem::exists(promol_nci_xyz.back()), "XYZ file doesn't exist: " + promol_nci_xyz.back().string(), std::cout);
			n_xyz++;
		}
		err_checkf(n_xyz >= 2,
			"Usage: -promol_nci <frag1.xyz> <frag2.xyz> [frag3.xyz ...] [rcut1=0.95] [rcut2=0.75] [rho_abs_max] [rdg_max] [colour_max=0.015]",
			std::cout);
		i += n_xyz - 2;

		double *optional_values[] = {
			&properties.promol_nci_rcut1,
			&properties.promol_nci_rcut2,
			&properties.promol_nci_rho_abs_max,
			&properties.promol_nci_rdg_max,
			&properties.promol_nci_colour_max
		};
		int optional_index = 0;
		while (optional_index < 5 && i + 3 + optional_index < argc)
		{
			try
			{
				size_t consumed = 0;
				const std::string &candidate = arguments[i + 3 + optional_index];
				const double value = std::stod(candidate, &consumed);
				if (consumed != candidate.size())
					break;
				*optional_values[optional_index] = value;
				optional_index++;
			}
			catch (...)
			{
				break;
			}
		}

		i += 2 + optional_index;
	}
	else if (temp == "-promol_nci_single_thread")
		properties.promol_nci_single_threaded = true;
	else if (temp == "-rdg")
		properties.rdg = true;
	else if (temp == "-rho")
		properties.rho = true; //the help advertised it; until now a density cube was only a by-product of -lap and friends
	else if (temp == "-rho_cube")
	{
		string wfn_name = arguments[i + 1];
		std::cout << "Reading wavefunction: " << wfn_name << endl;
		WFN wavy = WFN(wfn_name);
		std::cout << "Assigning ECPs" << endl;
		if (ECP)
			wavy.set_has_ECPs(true);
		std::cout << "Starting cube calculation" << endl;
		wavy.write_rho_cube(properties.radius, properties.resolution);
		finished = true; return true;
	}
	else if (temp.find("-s_rho") < 1)
		properties.s_rho = true;
	else if (temp == "-atom_dens_diff")
	{
		filesystem::path name_wfn_1 = arguments[i + 1];
		filesystem::path name_wfn_2 = arguments[i + 2];

		subtract_dens_from_gbw(name_wfn_1, name_wfn_2, 2, 0.05);
		finished = true; return true;
	}
	else if (temp == "-spherical_aver_fukui")
	{
		filesystem::path wfn1_name = arguments[i + 1];
		filesystem::path wfn2_name = arguments[i + 2];
		WFN *wavy1 = new WFN(wfn1_name);
		WFN *wavy2 = new WFN(wfn2_name);
		ofstream outputFile("fukui_averaged_density_wfn.dat");
		for (double r = 0.001; r < 10.0; r += 0.001)
		{
			// double dens = calc_grid_averaged_at_r_from_cube(cube_from_file, r, 360, 5800);
			double dens = calc_fukui_averaged_at_r(*wavy1, *wavy2, r, 5810, 5810);
			outputFile << r << " " << dens << "\n";
		}
		outputFile.close();
		std::cout << "Data written to output.dat" << endl;
		delete (wavy1);
		delete (wavy2);
		finished = true; return true;
	}
	else if (temp == "-spherical_aver_hirsh")
	{
		string wfn_name = arguments[i + 1];
		std::cout << "Reading wavefunction: " << wfn_name << endl;
		WFN *wavy = new WFN(wfn_name);
		std::cout << "Assigning ECPs" << endl;
		if (ECP)
			wavy->set_has_ECPs(true);
		std::cout << "Starting spherical averaging" << endl;
		double dens;

		for (int index_atom = 0; index_atom < wavy->get_ncen(); index_atom += 1)
		{
			ofstream outputFile("hirsh_averaged_density_" + std::to_string(index_atom) + ".dat");
			for (double r = 0.001; r < 5.0; r += 0.002)
			{
				dens = calc_hirsh_grid_averaged_at_r(*wavy, index_atom, r, 360, 5800);
				outputFile << r << " " << dens << "\n";
			}
			outputFile.close();
		}
		std::cout << "Data written to output.dat" << endl;
		delete (wavy);
		finished = true; return true;
	}
	else if (temp == "-spherical_harmonic")
	{
		spherical_harmonic_test();
		finished = true; return true;
	}
	else if (temp == "-spherical_atoms") {
		write_spherical_atoms();
		finished = true; return true;
	}
	else if (temp == "-cube_density" || temp == "-cube")
	{
		cube_density = arguments[i + 1];
		err_checkf(std::filesystem::exists(cube_density), temp + ": cube file does not exist: " + cube_density.string(), std::cout);
	}
	else if (temp == "-cubeb")
		cube::binary_output = true;
	else if (temp == "-cube_convert")
	{
		// -cube_convert <in> <out>: text .cube <-> binary .cubeb, the format follows the extension of <out>
		err_checkf(i + 2 < argc, "-cube_convert needs <in> <out>", std::cout);
		WFN wavy(e_origin::cub);
		cube in(arguments[i + 1], true, wavy, std::cout, true);
		err_checkf(in.write_file(std::filesystem::path(arguments[i + 2])), "cube conversion failed", std::cout);
		finished = true; return true;
	}
	else if (temp == "-calc_dens_1D")
	{
		err_checkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -calc_dens_1D to specify a wavefunction.", std::cout);
		err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEFORE -calc_dens_1D to specify an auxiliary basis set.", std::cout);
		int atom_idx_1 = 0;
		int atom_idx_2 = 1;
		int gridpoints = 1000;
		double padding = 2.0;
		if (string(arguments[i + 1]).find("-") != 0) {
			atom_idx_1 = stoi(arguments[i + 1]);
		}
		if (string(arguments[i + 2]).find("-") != 0) {
			atom_idx_2 = stoi(arguments[i + 2]);
		}
		if (string(arguments[i + 3]).find("-") != 0) {
			gridpoints = stoi(arguments[i + 3]);
		}
		if (string(arguments[i + 4]).find("-") != 0) {
			padding = stod(arguments[i + 4]);
		}

		WFN wavy(wfn);
		string out = "Calculating 1D density between atoms " + wavy.get_atom_label(atom_idx_1) + " (" + std::to_string(atom_idx_1) + ") and " + wavy.get_atom_label(atom_idx_2) + " (" + std::to_string(atom_idx_2) + ") with " + std::to_string(gridpoints) + " gridpoints and " + std::to_string(padding) + " Angstrom padding.";
		std::cout << out << std::endl;
		get1DGridData(wavy, aux_basis, atom_idx_1, atom_idx_2, gridpoints, padding);
		finished = true; return true;
	}
	else if (temp == "-rho_at_points")
	{
		//-wfn <file> -rho_at_points <points>: total density at the points (count, then x y z in bohr per line) to <points>.rho
		err_checkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -rho_at_points.", std::cout);
		WFN wavy(wfn);
		ifstream in(arguments[i + 1]);
		int n = 0; in >> n;
		std::vector<d3> pts(n); vec rho(n);
		for (int p = 0; p < n; p++) in >> pts[p][0] >> pts[p][1] >> pts[p][2];
#pragma omp parallel for
		for (int p = 0; p < n; p++) rho[p] = wavy.compute_dens(pts[p]);
		ofstream out(arguments[i + 1] + ".rho"); out << setprecision(12);
		for (int p = 0; p < n; p++) out << rho[p] << "\n";
		finished = true; return true;
	}
	else
		return false;
	return true;
};

//RI fitting, featomic descriptors and atom classification
bool options::digest_ri_options(const std::string &temp, int &i)
{
	using namespace std;
	const int argc = (int)arguments.size();
	if (temp == "-geometry_aid_cutoff") {
		//3.5 selects the c_only model family, 3.0 the dirty one. Both descriptors are 42,042 long and
		//nothing downstream rejects the wrong one, so the value used is echoed.
		err_checkf(i + 1 < argc, "-geometry_aid_cutoff needs a radius in Angstrom", std::cout);
		geometry_aid_cutoff = std::stod(arguments[++i]);
		std::cout << "  geometry-aid SOAP cutoff radius set to " << geometry_aid_cutoff << " A ("
				  << (geometry_aid_cutoff > 3.25 ? "c_only" : "dirty") << " model family)" << std::endl;
	}
	else if (temp == "-geometry_aid_center_weight") {
		err_checkf(i + 1 < argc, "-geometry_aid_center_weight needs a weight", std::cout);
		geometry_aid_center_weight = std::stod(arguments[++i]);
	}
	else if (temp == "-geometry_aid_metals") {
		//Zn joins the species list as the stand-in for every metal: the row-2 specialist sees the heavy atoms
		//labelled and every element outside the eleven as Zn, 12 species give 78 pairs, length 49,686
		geometry_aid_metals = true;
	}
	else if (temp == "-classify_atoms") {
		//-wfn <structure> -classify_atoms <model.bin> [<out.npy>], default probabilities.npy
		err_checkf(i + 1 < argc, "-classify_atoms needs the model file", std::cout);
		geometry_aid_model = arguments[++i];
		err_checkf(std::filesystem::exists(geometry_aid_model), "The geometry-aid model does not exist: " + geometry_aid_model.string(), std::cout);
		classify_atoms_out = "probabilities.npy";
		if (i + 1 < argc && !arguments[i + 1].empty() && arguments[i + 1][0] != '-')
			classify_atoms_out = arguments[++i];
	}
	else if (temp == "-classify_atoms_list") {
		//<list> <model.bin>, one structure path per line, each written as <path>.probs.npy
		err_checkf(i + 2 < argc, "-classify_atoms_list needs the list and the model file", std::cout);
		classify_structures = geometry_aid::read_structure_list(arguments[++i]);
		geometry_aid_model = arguments[++i];
		err_checkf(std::filesystem::exists(geometry_aid_model), "The geometry-aid model does not exist: " + geometry_aid_model.string(), std::cout);
	}
	else if (temp == "-calc_featomic_descriptor")
		//Olex2 calls exactly this and reads descriptor.npy from the working directory
		calc_featomic_descriptor = true;
	else if (temp == "-calc_featomic_descriptors") {
		//One structure path per line, each descriptor written as <path>.npy. One process for all of
		//them because the calculator setup costs 0.7 s against 0.0009 s per atom (cached_calculator).
		err_checkf(i + 1 < argc, "-calc_featomic_descriptors needs the list file", std::cout);
		featomic_structures = geometry_aid::read_structure_list(arguments[++i]);
	}
	else if (temp == "-npa")
		npa = true;
	else if (temp == "-ibo")
		ibo = true;
	else if (temp == "-ibo_cube") {
		err_checkf(i + 1 < argc, "-ibo_cube needs the IBO list", std::cout);
		properties.ibo_cube = arguments[++i];
	}
	else if (temp == "-nbo_cube" || temp == "-nlmo_cube") {
		err_checkf(i + 1 < argc, temp + " needs the NBO list", std::cout);
		(temp == "-nbo_cube" ? properties.nbo_cube : properties.nlmo_cube) = arguments[++i];
	}
	else if (temp == "-nbo_lewis") {
		err_checkf(i + 1 < argc, "-nbo_lewis needs a Lewis structure file", std::cout);
		properties.nbo_lewis = arguments[++i];
		err_checkf(std::filesystem::exists(properties.nbo_lewis),
				   "-nbo_lewis: Lewis structure file doesn't exist: " + properties.nbo_lewis.string(), std::cout);
	}
	else if (temp == "-npa_summary") {
		npa = true;
		npa_orbitals = false;
	}
	//-rgbi <wfn>, like -nbo and -eli_analysis; the bare flag still reads -wfn
	else if (temp == "-rgbi") {
		rgbi = true;
		if (i + 1 < argc && arguments[i + 1][0] != '-') {
			wfn = arguments[i + 1];
			err_checkf(std::filesystem::exists(wfn), "-rgbi: wavefunction does not exist: " + wfn.string(), std::cout);
			i += 1;
		}
	}
	else if (temp == "-rgbi_no_sym") {
		rgbi = true;
		rgbi_no_sym = true;
	}
	else if (temp == "-rgbi_EVs") {
		rgbi = true;
		rgbi_EVs = true;
	}
	else if (temp == "-rgbi_theta") {
		rgbi = true;
		rgbi_theta = true;
	}
	else if (temp == "-rgbi_legacy_cutoff") {
		rgbi = true;
		rgbi_legacy_cutoff = true;
	}
	else if (temp == "-rgbi_basis") {
		err_checkf(i + 1 < argc, "Not enough arguments for -rgbi_basis. Use 'nao' or 'ano'.", std::cout);
		rgbi = true;
		std::string basis = arguments[i + 1];
		std::transform(basis.begin(), basis.end(), basis.begin(),
			[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
		if (basis == "nao")
			rgbi_orbital_basis = RGBIOrbitalBasis::NAO;
		else if (basis == "ano")
			rgbi_orbital_basis = RGBIOrbitalBasis::ANO;
		else
			err_checkf(false, "Invalid -rgbi_basis value '" + basis + "'. Use 'nao' or 'ano'.", std::cout);
	}
	else if (temp == "-rgbi_groups") {
		int n = 1;
		ivec2 group_set;
		while (i + n < argc && string(arguments[i + n]).find("-") > 0) {
			group_set.push_back(parse_rgbi_group_indices(arguments[i + n]));
			n++;
		}
		if (!group_set.empty()) {
			validate_rgbi_group_set(group_set);
			rgbi_group_sets.push_back(group_set);
		}
		i += n - 1;
		rgbi = true;
	}
	else if (temp == "-RI_FIT" || temp == "-ri_fit")
	{
		RI_FIT = true;
		partition_type = PartitionType::RI;
		aux_basis = get_aux_basis(argc, arguments, i);
	}
	else if (temp == "-multipole_moments") {
		err_checkf(i + 1 < argc, "-multipole_moments needs a scheme, e.g. -multipole_moments Hirshfeld 2", std::cout);
		std::string scheme = arguments[++i];
		std::transform(scheme.begin(), scheme.end(), scheme.begin(),
			[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
		if (scheme == "hirshfeld" || scheme == "hirsh")
			multipole_scheme = MultipoleScheme::HIRSHFELD;
		else if (scheme == "tfvc")
			multipole_scheme = MultipoleScheme::TFVC;
		else if (scheme == "mbis")
			multipole_scheme = MultipoleScheme::MBIS;
		else if (scheme == "embis")
			multipole_scheme = MultipoleScheme::EMBIS;
		else if (scheme == "nuclear")
			multipole_scheme = MultipoleScheme::NUCLEAR;
		else if (scheme == "mulliken")
			multipole_scheme = MultipoleScheme::MULLIKEN;
		else if (scheme == "sanderson")
			multipole_scheme = MultipoleScheme::SANDERSON;
		else
			err("Unknown scheme for -multipole_moments: " + arguments[i] + " (Hirshfeld, TFVC, MBIS, EMBIS, Nuclear, Mulliken or Sanderson)", std::cout);
		long order = 0;
		if (i + 1 < argc) {
			char* end = nullptr;
			const long parsed = std::strtol(arguments[i + 1].c_str(), &end, 10);
			if (end != arguments[i + 1].c_str() && *end == '\0') {
				order = parsed;
				++i;
			}
		}
		if (multipole_scheme == MultipoleScheme::NUCLEAR || multipole_scheme == MultipoleScheme::MULLIKEN || multipole_scheme == MultipoleScheme::SANDERSON)
			multipole_lmax = 0;
		else {
			err_checkf(order >= 0 && order <= 8, "-multipole_moments: the order must be between 0 and 8", std::cout);
			multipole_lmax = static_cast<int>(order);
		}
		RI_FIT = true;
		partition_type = PartitionType::RI;
	}
	else if (temp == "-multipole_strength") {
		multipole_strength = std::stod(arguments[++i]);
		err_checkf(multipole_strength > 0.0, "-multipole_strength must be positive", std::cout);
	}
	else if (temp == "-write_ri_coefs") {
		WFN wavy(wfn);
		WFN wavy_aux = generate_aux_wfn(wavy, aux_basis, false);
		DensityFitting::CONFIG config = DensityFitting::config_from_options(*this);
		//config.restrain_type = DensityFitting::RESTRAINT_TYPE::SIMPLE_AND_TIK;
		//config.charge_scheme = DensityFitting::CHARGE_SCHEME::HIRSHFELD;
		vec ri_coefs = DensityFitting::density_fit(wavy, wavy_aux, config);
		npy::npy_data<double> np_coeffs;
		np_coeffs.data = ri_coefs;
		np_coeffs.fortran_order = false;
		np_coeffs.shape = { static_cast<unsigned long>(ri_coefs.size()) };
		npy::write_npy("RI_COEFS.npy", np_coeffs);
		finished = true; return true;
	}
	else if (temp == "-repulsion_overlap") {
		repulsion_overlap = std::stod(arguments[++i]);
		err_checkf(repulsion_overlap >= 0.0, "-repulsion_overlap must not be negative", std::cout);
	}
	else if (temp == "-repulsion_exchange") {
		const std::string f = arguments[++i];
		repulsion_exchange = f == "dirac" ? 0 : f == "pbe" ? 1 : f == "b88" ? 2 : f == "r2scan" ? 3 : -1;
		err_checkf(repulsion_exchange >= 0, "-repulsion_exchange takes dirac, pbe, b88 or r2scan", std::cout);
	}
	else if (temp == "-interaction_energy") {
		//-interaction_energy <A> <B>: electrostatics between two fitted densities. With -SALTED <model-dir> both are predicted
		//from the model, otherwise A and B are wavefunctions and each is RI-fitted with the -ri_fit basis.
		//<A> <A.npy> <B> <B.npy> takes coefficient files in the fitted_multipoles layout instead
		err_checkf(i + 2 < argc, "-interaction_energy needs two structure files", std::cout);
		const bool from_files = i + 4 < argc && std::filesystem::path(arguments[i + 2]).extension() == ".npy";
		if (SALTED && !from_files) std::cout << "Predicting both densities with the SALTED model in " << salted_model_dir << std::endl;
		Gaussian_Molecule A(WFN(arguments[i + 1]), *this, from_files ? arguments[i + 2] : "");
		Gaussian_Molecule B(WFN(arguments[from_files ? i + 3 : i + 2]), *this, from_files ? arguments[i + 4] : "");
		DensityFitting::print_interaction_energy(A.interaction(B, repulsion_overlap, repulsion_exchange), A.basis(), B.basis(), std::cout);
		finished = true; return true;
	}
	else if (temp == "-interaction_energies") {
		//-interaction_energies <job>: every molecule pair in contact in the crystal, run after the options are parsed
		err_checkf(i + 1 < argc, "-interaction_energies needs a job file", std::cout);
		interaction_energies_job = arguments[++i];
		err_checkf(std::filesystem::exists(interaction_energies_job), "The job file does not exist: " + interaction_energies_job.string(), std::cout);
	}
	else if (temp == "-RI_CUBE" || temp == "-ri_cube")
	{
		err_checkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -RI_CUBE to specify a wavefunction.", std::cout);
		err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEFORE -RI_CUBE to specify an auxiliary basis set.", std::cout);
		std::string coef_file = arguments[i + 1];
		std::vector<unsigned long> shape{};
		bool fortran_order;
		vec coefs{};

		npy::LoadArrayFromNumpy(coef_file, shape, fortran_order, coefs);


		WFN wavy(wfn);
		WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);

		int nr_coefs = 0;
		for (const atom &atm : wavy_aux.get_atoms()) {
			int prim = 0;
			for (int shell = 0; shell < atm.get_shellcount_size(); shell++) {
				const int type = atm.get_basis_set_entry(prim).get_type();
				nr_coefs += 2 * type + 1;
				prim += atm.get_shellcount(shell);
			}
		}

		std::cout << coefs.size() << " vs. " << nr_coefs << " ceofficients" << std::endl;

		// First name of coef_file, second name of xyz file
		cube_from_coef_npy(coefs, wavy_aux);

		// std::string aux_basis = arguments[i + 1];
		//gen_CUBE_for_RI(wavy, "def2_qzvppd_rifit", this);
		//gen_CUBE_for_RI(wavy, "def2_universal_jkfit", this);
		//gen_CUBE_for_RI(wavy, "combo-basis-fit", this);
		//gen_CUBE_for_RI(wavy, "cc-pvqz-jkfit", this);

		finished = true; return true;
	}

	else if (temp == "-test_RI")
	{
		err_checkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -test_RI to specify a wavefunction.", std::cout);
		err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEFORE -test_RI to specify an auxiliary basis set.", std::cout);


		WFN wavy(wfn);
		WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);
		DensityFitting::demonstrate_enhanced_density_fitting(wavy, wavy_aux);
		finished = true; return true;

	}
	else if (temp == "-RI_WFN_DIFF") {
		err_checkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -RI_WFN_DIFF to specify a wavefunction.", std::cout);
		err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEFORE -RI_WFN_DIFF to specify an auxiliary basis set.", std::cout);
		WFN wavy(wfn);
		WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);
		DensityFitting::QM_RI_difference_cube(wavy, wavy_aux);
		finished = true; return true;

		//finished = true; return true;
	}
	else
		return false;
	return true;
};

//X-ray constrained wavefunction fitting
//`[robust] sigma | unit | stl <a> | shelx <a> <b> | shelx_opt <a0> <b0> | chebyshev <A0> [A1 ..] | chebyshev_opt <n>`
//and, under robust, Tukey's c (default 6; after chebyshev, whose count is free, the last number and required):
//the settings file's `weight` and the CLI's -xcw_weight. Numbers are taken while the next token starts like one
static void read_xcw_weight(std::istream& is, options::XCW_settings& s)
{
	using ws = options::SF_settings::weight_scheme;
	static const std::map<std::string, std::pair<ws, int>> schemes = { { "sigma", { ws::sigma, 0 } }, { "unit", { ws::unit, 0 } },
		{ "stl", { ws::stl, 1 } }, { "shelx", { ws::shelx, 2 } }, { "shelx_opt", { ws::shelx_opt, 2 } },
		{ "chebyshev", { ws::chebyshev, -1 } }, { "chebyshev_opt", { ws::chebyshev_opt, 1 } } };
	const auto word = [&is]() { std::string w; is >> w; std::transform(w.begin(), w.end(), w.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); }); return w; };
	std::string scheme = word();
	const bool robust = scheme == "robust";
	if (robust) scheme = word();
	const auto it = schemes.find(scheme);
	if (it == schemes.end())
		throw std::runtime_error("Unknown weight '" + scheme + "', expected [robust] sigma, unit, stl, shelx, shelx_opt, chebyshev or chebyshev_opt");
	vec p;
	for (double x; is >> std::ws, (std::isdigit(is.peek()) || is.peek() == '-' || is.peek() == '+' || is.peek() == '.') && (is >> x);) p.push_back(x);
	const int n = it->second.second;
	double c = 6.0;
	if (robust && ((n < 0 && p.size() > 1) || (n >= 0 && p.size() == static_cast<size_t>(n) + 1))) {
		c = p.back();
		p.pop_back();
	}
	else if (robust && n < 0) p.clear();
	if ((n >= 0 && p.size() != static_cast<size_t>(n)) || (n < 0 && p.empty()))
		throw std::runtime_error("Wrong number of values after 'weight " + std::string(robust ? "robust " : "") + scheme + "'");
	if (std::any_of(p.begin(), p.end(), [&](double v) { return (it->second.first == ws::shelx || it->second.first == ws::shelx_opt) && v < 0; })
		|| (it->second.first == ws::chebyshev_opt && (p[0] < 1 || p[0] != std::floor(p[0]))) || !(c > 0))
		throw std::runtime_error("Out-of-range value after 'weight " + scheme + "'");
	s.weight = it->second.first;
	s.weight_params = p;
	s.robust_c = robust ? c : 0.0;
}

bool options::digest_xcw_options(const std::string &temp, int &i)
{
	using namespace std;
	const int argc = (int)arguments.size();
	if (temp == "-convert_XCW")
	{
		err_checkf(argc >= i + 3, "Not enough arguments for -convert_XCW\nPlease provide at least stdout name and lambda step!", std::cout);
		std::string stdo = arguments[i + 1];
		std::string step = arguments[i + 2];
		convert_tonto_XCW_lambda_steps(stdo, step, debug, *this);
		finished = true; return true;
	}
	else if (temp == "-do_XCW") {
		do_XCW = true;
	}
	else if (temp == "-XCW_settings") {
		xcw_settings_path = arguments[i + 1];
	}
	else if (temp == "-xcw_weight") {
		err_checkf(i + 1 < argc, "-xcw_weight takes [robust] sigma, unit, stl, shelx, shelx_opt, chebyshev or chebyshev_opt", std::cout);
		//every following argument up to the next option (a negative number is not one)
		std::string words = arguments[i + 1];
		for (int k = i + 2; k < argc && (arguments[k][0] != '-' || (arguments[k].size() > 1 && (std::isdigit(static_cast<unsigned char>(arguments[k][1])) || arguments[k][1] == '.'))); k++)
			words += " " + arguments[k];
		std::istringstream is(words);
		read_xcw_weight(is, xcw_settings);
	}
	else if (temp == "-xcw_weight_rescale") {
		err_checkf(i + 1 < argc && (arguments[i + 1] == "0" || arguments[i + 1] == "1"), "-xcw_weight_rescale takes 0 or 1", std::cout);
		xcw_settings.weight_rescale = arguments[++i] == "1";
	}
	else if (temp == "-xcw_partition") {
		err_checkf(i + 1 < argc, "-xcw_partition takes grid or ri", std::cout);
		xcw_partition = arguments[++i];
		err_checkf(xcw_partition == "grid" || xcw_partition == "ri", "-xcw_partition takes grid or ri, not " + xcw_partition, std::cout);
	}
	else
		return false;
	return true;
};

//development and test-only switches
bool options::digest_dev_options(const std::string &temp, int &i)
{
	using namespace std;
	if (temp == "-tune") //applied before the main loop; this only checks there is a value
	{
		(void)arguments[++i];
		return true;
	}
	if (temp == "-lahvatest")
	{
		//_test_lahva();
		finished = true; return true;
	}
	else if (temp == "-NNLS_TEST")
	{
		test_NNLS();
		finished = true; return true;
	}
	else if (temp == "-test")
		std::cout << "Running in test mode!" << endl, test = true;
	else if (temp == "-partitioning_test")
	{
		calc_partition_densities();
	}
	else if (temp == "-lukas_test")
	{
		//Check that wfn is not empty
		err_checkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -lukas_test to specify a wavefunction.", std::cout);
		err_checkf(!salted_model_dir.empty(), "No SALTED model directory specified! Use -SALTED option BEFORE -lukas_test to specify a model directory.", std::cout);

		WFN wavy(wfn);
		SALTEDPredictor SP(wavy, *this);
		string df_basis_name = SP.get_dfbasis_name();
		filesystem::path salted_model_path = SP.get_salted_filename();
		log_file << "Using " << salted_model_path << " for the prediction" << endl;
		if (!SP.basis_set_loaded()) {
			std::shared_ptr<BasisSet> aux_basis = BasisSetLibrary::get_basis_set(df_basis_name);
			load_basis_into_WFN(SP.wavy, aux_basis);
		}
		vec coefs = SP.gen_SALTED_densities();

		cube atom_cube = calc_cube_ML(coefs, SP.wavy);
		atom_cube.write_file("DBA_total.cube");

		for (int atm_idx = 0; atm_idx < wavy.get_ncen(); atm_idx++) {
			atom_cube = calc_cube_ML(coefs, SP.wavy, atm_idx);
			atom_cube.write_file("DBA_atom_" + std::to_string(atm_idx) + ".cube");
		}
	}
	else
		return false;
	return true;
};

//The analysis an option's prefix names: an unclaimed flag with one of these prefixes is a misspelt
//option of that analysis and is refused rather than silently ignored.
const char *owning_analysis(const std::string &flag)
{
	static const std::pair<const char *, const char *> families[] = {
		{"-rgbi", "RGBI"}, {"-npa", "NPA"}, {"-ibo", "IBO"}, {"-nbo", "NBO"}, {"-nrt", "NRT"},
		{"-nao", "NAO"}, {"-eli", "ELI-D"}, {"-elf", "ELF"}, {"-qtaim", "QTAIM"},
		{"-basin", "QTAIM basin"}, {"-topology", "topology"}, {"-eqc", "EQC"}};
	for (const auto &f : families)
		if (flag.rfind(f.first, 0) == 0)
			return f.second;
	return nullptr;
}

namespace {
	//Mirrors run_app_impl's chain of analyses that return before the RGBI/NPA block is reached;
	//an analysis added to that chain needs a line here.
	const char *quit_early_analysis(const options &o)
	{
		if (o.qct) return "-qct";
		if (o.fukui_analysis_run) return "-fukui_analysis";
		if (o.eqc) return "-eqc";
		if (o.fba) return "-fba";
		if (o.eli_analysis_run) return o.qtaim_only ? "-qtaim" : "-eli_analysis";
		if (o.fract) return "-fractal";
		if (o.promol_nci) return "-promol_nci";
		if (!o.interaction_energies_job.empty()) return "-interaction_energy";
		if (!o.hirshfeld_surface.empty()) return "-hirshfeld_surface";
		if (!o.wfn2.empty()) return "-wfn2";
		if (!o.pol_wfns.empty()) return "-polarizabilities";
		if (o.combined_tsc_calc || o.cif_based_combined_tsc_calc) return "-merge";
		if (o.iam_switch) return "-IAM";
		if (!o.cube_density.empty()) return "-cube_density";
		//the conditions that skip the wavefunction block (d_sfac_scan has no option that sets it)
		if (o.properties.calc()) return "a property cube (-rho, -lap, -esp, ...)";
		if (o.do_XCW) return "-do_XCW";
		if (o.gbw2wfn) return "-gbw2wfn";
		return nullptr;
	}
}

void options::refuse_unread_bonding_options()
{
	const char *early = quit_early_analysis(*this);
	//-fba runs RGBI (reading the -rgbi_* modifiers, each of which sets rgbi) and its own NPA, which
	//reads neither -npa nor -npa_summary; -do_XCW runs RGBI on each refined wavefunction unless a
	//property cube takes the run first
	const bool rgbi_read_early = fba || (do_XCW && !properties.calc());
	const char *bonding = (rgbi && !rgbi_read_early) ? "-rgbi" : (npa ? "-npa" : (ibo ? "-ibo" : nullptr));
	if (early != nullptr && bonding != nullptr)
		err_checkf(false, std::string("Cannot do both ") + early + " and " + bonding + " in one run: " +
							  early + " ends the run before " + bonding + " would be reached, so " +
							  bonding + " would be dropped without a word. Run them one at a time",
				   log_file);
	//The NBO handlers finish the command line where they run, so an -nbo_*/-nrt_* option reaching
	//here is read by no analysis.
	const char *running = early != nullptr ? early : bonding;
	if (running == nullptr)
		return;
	for (const std::string &raw : arguments)
	{
		std::string o = raw;
		if (o.size() > 1 && o[0] == '-' && isalpha(static_cast<unsigned char>(o[1])))
			std::replace(o.begin() + 1, o.end(), '-', '_');
		if (nbo_family_suboptions().count(o) == 0)
			continue;
		err_checkf(false, "Option " + o + " is an NBO/NRT option and this run is " + running +
							  (fba ? ", which runs NBO and NRT at their defaults and does not read it: the "
									 "option would be ignored and the run would still exit 0. Use "
									 "-nbo <wavefunction> -nrt for NRT options, or drop " + o
								   : ", which does not read it: no NBO analysis was asked for, so the "
									 "option would be ignored and the run would still exit 0. Add "
									 "-nbo <wavefunction>, or drop " + o),
				   log_file);
	}
}

void options::digest_options()
{
	using namespace std;
	// Lets print what was the command line, for debugging
	if (debug)
	{
		std::cout << " Recap of input:\nsize: " << arguments.size() << endl;
	}
	//Flags with process-wide state start from their defaults, not from the previous run in this process
	tuning_end_run();
	basin_flags_reset();
	cube::binary_output = false;
	throughput::set_enabled(false);
	throughput::reset();
	//-tune first: some options run their job inside the loop below
	for (size_t i = 0; i + 1 < arguments.size(); i++)
		if (arguments[i] == "-tune")
		{
			const std::string &knob = arguments[++i];
			const size_t eq = knob.find('=');
			const std::string name = knob.substr(0, eq);
			const char *old = tuning(name.c_str());
			run_tuned.emplace_back(name, old ? std::optional<std::string>(old) : std::nullopt);
			set_tuning(name, eq == std::string::npos ? "1" : knob.c_str() + eq + 1);
		}
	//the NBO handlers compare raw tokens, so -nbo-json is spelled -nbo_json for them here
	for (int i = 0; i < arguments.size(); i++)
		if (arguments[i].size() > 1 && arguments[i][0] == '-' && isalpha(static_cast<unsigned char>(arguments[i][1])))
		{
			string o = arguments[i];
			replace(o.begin() + 1, o.end(), '-', '_');
			if (nbo_family_suboptions().count(o))
				arguments[i] = o;
		}
	// This loop figures out command line options
	for (int i = 0; i < arguments.size(); i++)
	{
		if (debug)
			std::cout << arguments[i] << endl;
		string temp = arguments[i];
		if (temp.find("-") > 0)
			continue;
		if (finished)
			return;
		//-rgbi-groups is -rgbi_groups; a negative number keeps its dashes (no letter after the first)
		if (temp.size() > 1 && isalpha(static_cast<unsigned char>(temp[1])))
			replace(temp.begin() + 1, temp.end(), '-', '_');
		//The digesters index arguments[i + n] and call stoi/stod directly; a flag that is
		//last on the line or followed by a non-number used to die as a bare "invalid stod
		//argument" with no hint which option it was
		bool claimed = false;
		try
		{
			claimed = digest_io_options(temp, i) || digest_run_options(temp, i) ||
					  digest_partition_options(temp, i) || digest_property_options(temp, i) ||
					  digest_ri_options(temp, i) || digest_xcw_options(temp, i) ||
					  digest_dev_options(temp, i);
		}
		catch (const missing_argument &)
		{
			err_checkf(false, "Option " + temp + " needs more values than were given after it", log_file);
		}
		//one-shot options run their whole calculation in here; running out of memory is not a bad value
		catch (const std::bad_alloc &)
		{
			err_checkf(false, "Option " + temp + " ran out of memory", log_file);
		}
		catch (const std::exception &e)
		{
			err_checkf(false, "Option " + temp + ": " + e.what() + " - check the value(s) that follow it", log_file);
		}
		if (claimed)
			continue;
		//An inline handler that ran and failed returns false with finished set; not an unknown option.
		if (finished)
			return;
		//No digester wanted it. A value starting with a dash (a negative number) is not an option.
		if (!(temp.size() > 1 && isalpha(static_cast<unsigned char>(temp[1]))))
			continue;
		//look_for_debug() handles these before any digester sees them
		if (temp == "-v" || temp == "-v2" || temp == "-debug" || temp == "-h" || temp == "--h" ||
			temp == "-help" || temp == "--help")
			continue;
		//the -nbo/-nbo_external handlers finish the line and read only the tokens after their wavefunction,
		//so one reaching here is read by nothing
		if (nbo_family_suboptions().count(temp))
		{
			std::cerr << "WARNING: ignoring " << temp << ": NBO/NRT options are read only after -nbo <wavefunction>, -nbo_external, -nbo_parse or -convert_to_47" << endl;
			continue;
		}
		if (const char *family = owning_analysis(temp))
			err_checkf(false, "Unknown option " + temp + ": it names no " + family +
								  " option this build has. Check the spelling, -h lists them; "
								  "running on with it ignored would give you a default you did "
								  "not ask for",
					   log_file);
		log_file << "WARNING: ignoring unknown option " << temp << endl;
		//run_app_impl has pointed cout at the log, so stderr is the console
		std::cerr << "WARNING: ignoring unknown option " << temp << endl;
	}

	// SALTED predicts a density from atom positions.  Historically its
	// calculation path read those positions through -wfn, although the CLI
	// also documents -xyz for SALTED input.  Accept the documented form and
	// preserve an explicitly supplied -wfn when both options are present.
	if (SALTED && wfn.empty() && !xyz_file.empty())
	{
		wfn = xyz_file;
		if (debug)
			log_file << "Using -xyz input as the SALTED structure: " << wfn << endl;
	}
	//-multipole_moments without -ri_fit: auto_aux, whichever order the two came in
	if (multipole_lmax >= 0 && aux_basis.empty())
		aux_basis.push_back(std::make_shared<BasisSet>());
	err_checkf(!(qtaim_only && (eli_only || eli_contributions)),
		"-qtaim computes no ELI-D basins: use -eli_analysis for -eli_only or -eli_contributions", log_file);
	//last, so it sees every flag whichever order they were written in
	refuse_unread_bonding_options();
};

void options::loadXCWsettings() {
	auto lowercase = [](std::string s) {
		std::transform(s.begin(), s.end(), s.begin(),
			[](unsigned char c) { return std::tolower(c); });
		return s;
		};

	std::string conv_preset = "normal";
	bool check_hessian = false;
	std::string speed_preset = "normal_conv";

	/* 1: Classical Jayatilaka XWR
	   2: Ewald sum XWR */
	xcw_settings.XWR_type = 1;
	double quant_diff = 32768, diis_stop_damping = 32768, diis_stop_shift = 32768, max_diis_error = 32768, gradient = 32768, MaxP_diff = 32768, RMSP_diff = 32768, alpha = 32768, level_shift = 32768, start = 32768, end = 32768, step_size = 32768, xcw_strong_cutoff = xcw_settings.xcw_strong_cutoff;
	int max_scf_iterations = 32768, charge = 32768, multiplicity = 32768, n_params = 32768, refine_against = 32768;
	std::string basis_set_name = "Undefined";
	std::string df_basis_name, hessian_df_basis_name, guess_basis_name;
	bool grown = false, read_first_guess = false, nbo_output = false, xcw_gaussian_halt = false, mask = false;
	extinction::model extinction_model = extinction::model::none;
	bool extinction_aniso = false, extinction_refine = true;
	vec extinction_start = { 1e-4 };
	double wavelength = 0.0;
	bool grid_double = false;
	std::filesystem::path eri_read_path, eri_save_path, basis_overrides;
	occ::qm::SpinorbitalKind hf_type = occ::qm::SpinorbitalKind::Restricted;
	if (!std::filesystem::exists(xcw_settings_path)) {
		throw std::runtime_error("Settings file not found! Aborting run!");
	}
	else {

		std::ifstream input(xcw_settings_path);
		using Handler = std::function<void(std::istream&)>;
		std::unordered_map<std::string, Handler> handlers;
		handlers["conv"] = [&](std::istream& is) {
			if (!(is >> quant_diff))
				throw std::runtime_error("Expected value after 'conv'");
			};

		handlers["extinction"] = [&](std::istream& is) {
			//`extinction <model> [iso|aniso] [fixed] [start value]` in any order: a word names the model
			//or a flag, a number is the start value, and the next settings keyword ends it and is read
			//as usual, so `extinction shelx` need not be the last keyword; six numbers are the aniso tensor
			std::string token, next;
			vec starts;
			while (is >> token) {
				const std::string low = lowercase(token);
				if (low == "aniso" || low == "anisotropic") extinction_aniso = true;
				else if (low == "iso" || low == "isotropic") extinction_aniso = false;
				else if (low == "fixed" || low == "fix") extinction_refine = false;
				else if (extinction::from_string(low) != extinction::model::none) extinction_model = extinction::from_string(low);
				else if (handlers.count(low)) { next = low; break; }
				else {
					try { starts.push_back(std::stod(token)); }
					catch (const std::exception&) { throw std::runtime_error("Could not read '" + token + "' after 'extinction'"); }
				}
			}
			if (extinction_model == extinction::model::none)
				throw std::runtime_error("Expected a model (shelx, bc_gaussian or bc_lorentzian) after 'extinction'");
			if (starts.size() == 6 && !extinction_aniso)
				throw std::runtime_error("Six extinction start values are the anisotropic tensor and need 'aniso'");
			if (!starts.empty() && starts.size() != 1 && starts.size() != 6)
				throw std::runtime_error("Expected 1 or 6 extinction start values, got " + std::to_string(starts.size()));
			if (!starts.empty()) extinction_start = starts;
			if (!next.empty()) handlers.at(next)(is);
			};

		handlers["wavelength"] = [&](std::istream& is) {
			if (!(is >> wavelength))
				throw std::runtime_error("Expected value after 'wavelength'");
			};

		handlers["diis_damping"] = [&](std::istream& is) {
			if (!(is >> diis_stop_damping))
				throw std::runtime_error("Expected value after 'diis_damping'");
			};

		handlers["diis_shift"] = [&](std::istream& is) {
			if (!(is >> diis_stop_shift))
				throw std::runtime_error("Expected value after 'diis_shift'");
			};

		handlers["conv_diis"] = [&](std::istream& is) {
			if (!(is >> max_diis_error))
				throw std::runtime_error("Expected value after 'conv_diis'");
			};

		handlers["gradient"] = [&](std::istream& is) {
			if (!(is >> gradient))
				throw std::runtime_error("Expected value after 'gradient'");
			};

		handlers["maxp_diff"] = [&](std::istream& is) {
			if (!(is >> MaxP_diff))
				throw std::runtime_error("Expected value after 'MaxP_diff'");
			};

		handlers["rmsp_diff"] = [&](std::istream& is) {
			if (!(is >> RMSP_diff))
				throw std::runtime_error("Expected value after 'RMSP_diff'");
			};

		handlers["params"] = [&](std::istream& is) {
			if (!(is >> n_params))
				throw std::runtime_error("Expected value after 'params'");
			};

		handlers["damp"] = [&](std::istream& is) {
			if (!(is >> alpha))
				throw std::runtime_error("Expected value after 'damp'");
			};

		handlers["shift"] = [&](std::istream& is) {
			if (!(is >> level_shift))
				throw std::runtime_error("Expected value after 'shift'");
			};

		handlers["max_iter"] = [&](std::istream& is) {
			if (!(is >> max_scf_iterations))
				throw std::runtime_error("Expected value after 'max_iter'");
			};

		handlers["charge"] = [&](std::istream& is) {
			if (!(is >> charge))
				throw std::runtime_error("Expected value after 'charge'");
			};

		handlers["mult"] = [&](std::istream& is) {
			if (!(is >> multiplicity))
				throw std::runtime_error("Expected value after 'mult'");
			};

		handlers["f"] = [&](std::istream&) {
			refine_against = 1;
			};

		handlers["f2"] = [&](std::istream&) {
			refine_against = 2;
			};

		handlers["i_sigma"] = [&](std::istream& is) {
			if (!(is >> xcw_settings.i_sigma_cutoff))
				throw std::runtime_error("Expected value after 'i_sigma'");
			};

		handlers["weighted"] = [&](std::istream&) {
			xcw_settings.XWR_type = 2;
			};

		handlers["weight"] = [&](std::istream& is) {
			read_xcw_weight(is, xcw_settings);
			};

		handlers["weight_rescale"] = [&](std::istream& is) {
			std::string v;
			is >> v;
			v = lowercase(v);
			if (v == "off" || v == "0") xcw_settings.weight_rescale = false;
			else if (v == "on" || v == "1") xcw_settings.weight_rescale = true;
			else throw std::runtime_error("Expected on or off after 'weight_rescale'");
			};

		handlers["basis_set"] = [&](std::istream& is) {
			if (!(is >> basis_set_name))
				throw std::runtime_error("Expected basis set name");
			};
		handlers["basis_overrides"] = [&](std::istream& is) {
			std::string path;
			if (!(is >> std::quoted(path, '"', '\0')) || path.empty())
				throw std::runtime_error("Expected a JSON file after 'basis_overrides'");
			basis_overrides = path;
			if (basis_overrides.is_relative())
				basis_overrides = std::filesystem::absolute(xcw_settings_path).parent_path() / basis_overrides;
			};
		handlers["df_basis"] = [&](std::istream& is) {
			if (!(is >> df_basis_name))
				throw std::runtime_error("Expected a fitting basis name after 'df_basis'");
			};
		handlers["hessian_df_basis"] = [&](std::istream& is) {
			if (!(is >> hessian_df_basis_name))
				throw std::runtime_error("Expected a fitting basis name after 'hessian_df_basis'");
			};
		handlers["guess_basis"] = [&](std::istream& is) {
			if (!(is >> guess_basis_name))
				throw std::runtime_error("Expected a basis name after 'guess_basis'");
			};

		handlers["start"] = [&](std::istream& is) {
			if (!(is >> start))
				throw std::runtime_error("Expected value after 'start'");
			};

		handlers["end"] = [&](std::istream& is) {
			if (!(is >> end))
				throw std::runtime_error("Expected value after 'end'");
			};

		handlers["step_size"] = [&](std::istream& is) {
			if (!(is >> step_size))
				throw std::runtime_error("Expected value after 'step_size'");
			};

		handlers["grown"] = [&](std::istream&) {
			grown = true;
			};

		handlers["mask"] = [&](std::istream&) {
			mask = true;
			};

		handlers["rhf"] = [&](std::istream&) {
			hf_type = occ::qm::SpinorbitalKind::Restricted;
			};

		handlers["uhf"] = [&](std::istream&) {
			hf_type = occ::qm::SpinorbitalKind::Unrestricted;
			};

		handlers["sloppy"] = [&](std::istream&) {
			conv_preset = "sloppy";
			};

		handlers["normal"] = [&](std::istream&) {
			conv_preset = "normal";
			};

		handlers["tight"] = [&](std::istream&) {
			conv_preset = "tight";
			};

		handlers["very_tight"] = [&](std::istream&) {
			conv_preset = "very_tight";
			};

		handlers["slow_conv"] = [&](std::istream&) {
			speed_preset = "slow_conv";
			};

		handlers["normal_conv"] = [&](std::istream&) {
			speed_preset = "normal_conv";
			};

		//`soscf` is accepted and ignored: cold DIIS + warm TRAH is fixed (ccb566ce), and the key
		//turned the cold start second-order, which costs Fock builds (P1 lambda 0: 14 -> 26)
		handlers["soscf"] = [&](std::istream&) {};

		handlers["check_hessian"] = [&](std::istream&) {
			check_hessian = true;
			};

		handlers["fast_conv"] = [&](std::istream&) {
			speed_preset = "fast_conv";
			};

		handlers["load_wfn"] = [&](std::istream&) {
			read_first_guess = true;
			};
		handlers["nbo"] = [&](std::istream&) {
			nbo_output = true;
			};
		//The stored two-electron integrals depend on geometry and basis only; read_eri falls back to building
		//them when the file holds another basis or geometry
		handlers["save_eri"] = [&](std::istream& is) {
			std::string path;
			if (!(is >> path))
				throw std::runtime_error("Expected a path after 'save_eri'");
			eri_save_path = path;
			};
		handlers["read_eri"] = [&](std::istream& is) {
			std::string path;
			if (!(is >> path))
				throw std::runtime_error("Expected a path after 'read_eri'");
			eri_read_path = path;
			};

		handlers["i_double"] = [&](std::istream&) {
			grid_double = true;
			};


		handlers["gaussian_halt"] = [&](std::istream&) {
			xcw_gaussian_halt = true;
			};

		handlers["strong_cutoff"] = [&](std::istream& is) {
			if (!(is >> xcw_strong_cutoff))
				throw std::runtime_error("Expected value after 'strong_cutoff'");
			};

		std::string keyword;

		while (input >> keyword)
		{
			keyword = lowercase(keyword);

			auto it = handlers.find(keyword);
			if (it == handlers.end()) {
				throw std::runtime_error("Unknown keyword '" + keyword + "'");
			}

			it->second(input);
		}
	}

	if (conv_preset == "sloppy") {
		xcw_settings.quant_diff = 3e-5;
		xcw_settings.max_diis_error = 1e-4;
		xcw_settings.gradient = 5e-4;
		xcw_settings.MaxP_diff = 1e-4;
		xcw_settings.RMSP_diff = 1e-5;
		xcw_settings.max_scf_iterations = 100;
	}
	else if (conv_preset == "normal") {
		xcw_settings.quant_diff = 1e-6;
		xcw_settings.max_diis_error = 1e-5;
		xcw_settings.gradient = 7e-5;
		xcw_settings.MaxP_diff = 1e-5;
		xcw_settings.RMSP_diff = 1e-6;
		xcw_settings.max_scf_iterations = 100;
	}
	else if (conv_preset == "tight") {
		xcw_settings.quant_diff = 5e-7;
		xcw_settings.max_diis_error = 5e-6;
		xcw_settings.gradient = 3e-5;
		xcw_settings.MaxP_diff = 1e-6;
		xcw_settings.RMSP_diff = 1e-7;
		xcw_settings.max_scf_iterations = 100;
	}
	else if (conv_preset == "very_tight") {
		xcw_settings.quant_diff = 1e-7;
		xcw_settings.max_diis_error = 1e-6;
		xcw_settings.gradient = 1e-5;
		xcw_settings.MaxP_diff = 1e-7;
		xcw_settings.RMSP_diff = 1e-8;
		xcw_settings.max_scf_iterations = 100;
	}

	if (speed_preset == "slow_conv") {
		xcw_settings.alpha_damp = 0.8;
		xcw_settings.level_shift = 1;
		xcw_settings.diis_stop_damping = 1e-5;
		xcw_settings.diis_stop_shift = 1e-5;
	}
	else if (speed_preset == "normal_conv") {
		xcw_settings.alpha_damp = 0.5;
		xcw_settings.level_shift = 0.5;
		xcw_settings.diis_stop_damping = 1e-3;
		xcw_settings.diis_stop_shift = 1e-2;
	}
	else if (speed_preset == "fast_conv") {
		xcw_settings.method_apply_damping = false;
		xcw_settings.method_apply_shift = false;
	}

	if (basis_set_name == "Undefined") {
		throw std::runtime_error("Basis set name not specified in settings file! Aborting run!");
	}
	xcw_settings.basis_set_name = basis_set_name;
	xcw_settings.basis_overrides = basis_overrides;
	if (quant_diff != 32768) xcw_settings.quant_diff = quant_diff;
	if (diis_stop_damping != 32768) xcw_settings.diis_stop_damping = diis_stop_damping;
	if (diis_stop_shift != 32768) xcw_settings.diis_stop_shift = diis_stop_shift;
	if (max_diis_error != 32768) xcw_settings.max_diis_error = max_diis_error;
	if (gradient != 32768) xcw_settings.gradient = gradient;
	if (MaxP_diff != 32768) xcw_settings.MaxP_diff = MaxP_diff;
	if (RMSP_diff != 32768) xcw_settings.RMSP_diff = RMSP_diff;
	if (alpha != 32768) xcw_settings.alpha_damp = alpha;
	if (level_shift != 32768) xcw_settings.level_shift = level_shift;
	if (max_scf_iterations != 32768) xcw_settings.max_scf_iterations = max_scf_iterations;
	xcw_settings.charge = (charge != 32768) ? charge : xcw_settings.charge;
	xcw_settings.multiplicity = (multiplicity != 32768) ? multiplicity : xcw_settings.multiplicity;
	if (n_params != 32768) xcw_settings.n_params = n_params;
	else {
		xcw_settings.n_params = 1;
		std::cout << "Warning: Number of crystallographic parameters not set (use the params keyword)! Continuing with a default value of 1..." << std::endl;
	}
	xcw_settings.extinction_model = extinction_model;
	xcw_settings.extinction_aniso = extinction_aniso;
	xcw_settings.extinction_refine = extinction_refine;
	xcw_settings.extinction_start = extinction_start;
	xcw_settings.wavelength = wavelength;
	if (refine_against != 32768) xcw_settings.refine_against = refine_against;
	xcw_settings.xcw_start_value = (start != 32768) ? start : 0;
	xcw_settings.xcw_step_size = (step_size != 32768) ? step_size : 0.01;
	if (end != 32768) {
		xcw_settings.num_xcw_steps = static_cast<int>(std::round((end - xcw_settings.xcw_start_value) / xcw_settings.xcw_step_size)) + 1;
	}
	else {
		xcw_settings.num_xcw_steps = static_cast<int>(std::round((1 - xcw_settings.xcw_start_value) / xcw_settings.xcw_step_size)) + 1;
	}
	xcw_settings.grown = grown;
	xcw_settings.mask = mask;
	xcw_settings.hf_type = hf_type;
	xcw_settings.read_first_guess = read_first_guess;
	xcw_settings.grid_double = grid_double;
	xcw_settings.eri_read_path = eri_read_path;
	xcw_settings.eri_save_path = eri_save_path;
	xcw_settings.nbo_output = nbo_output;
	xcw_settings.df_basis_name = df_basis_name;
	xcw_settings.hessian_df_basis_name = hessian_df_basis_name;
	xcw_settings.guess_basis_name = guess_basis_name;
	xcw_settings.check_hessian = check_hessian;
	xcw_settings.xcw_gaussian_halt = xcw_gaussian_halt;
	xcw_settings.xcw_strong_cutoff = xcw_strong_cutoff;
}

std::string options::unrunnable_analysis() const
{
	if (!wfn.empty() || !occ.empty())
		return "";
	//These run only in run_app_impl's wavefunction branch, which is skipped silently without a wavefunction.
	const std::pair<bool, const char *> needs_a_wavefunction[] = {
		{rgbi, "-rgbi (Roby-Gould bond indices)"},
		{npa, "-npa (natural population analysis)"},
		{ibo, "-ibo (intrinsic bond orbitals)"},
		{!fchk.empty(), "-fchk (conversion to a formatted checkpoint)"},
		{!cif.empty() || !hkl.empty(), "-cif/-hkl (scattering factor table)"}};
	for (const auto &a : needs_a_wavefunction)
		if (a.first)
			return std::string(a.second) + " was asked for, but no wavefunction was given.\n" +
				   "       Name one with -wfn <file.gbw|.wfx|.fchk|.molden|.wfn> or with -occ;\n" +
				   "       this option does not take the wavefunction as a bare argument.";
	return "";
}

void convert_tonto_XCW_lambda_steps(const std::string &str, const std::string &lambda_step, bool debug, options &opt) {
	double lambda = 0.0;
	const double ls = stod(lambda_step);
	std::string jobname, line;
	std::filesystem::path stdout_file = str;
	err_checkf(std::filesystem::exists(stdout_file), "couldn't open or find " + stdout_file.string() + ", leaving", std::cout);
	std::ifstream rf(stdout_file.string().c_str(), std::ios::in);
	rf.seekg(0);
	while (rf.good() && line.find("Name ...") == std::string::npos) {
		getline_universal(rf, line);
	}
	jobname = split_string<std::string>(line, " ")[2];
	std::cout << "Conervting XCW wavefunctions with lambda step " + std::to_string(ls) + " and jobname: " + jobname << std::endl;

	//iterate over files in the same folder looking for .orbital_energies,restricted or .MO_energies,r and .molecular_orbitals,restricted or .MOs,r at the given lambda value
	//example filename: NiP3.molecular_orbitals,lambda=0.000000,restricted
	//make lambda value into a string with 6 decimal places
	std::stringstream ss;
	ss << std::fixed << std::setprecision(6) << lambda;
	std::string formatted_lambda = ss.str();
	std::filesystem::path energies_file = stdout_file.parent_path() / (jobname + ".orbital_energies,lambda=" + formatted_lambda + ",restricted");
	if (!std::filesystem::exists(energies_file))
		energies_file = stdout_file.parent_path() / (jobname + ".MO_energies,lambda=" + formatted_lambda + ",r");
	if (!std::filesystem::exists(energies_file))
		energies_file = stdout_file.parent_path() / (jobname + ".orbital_energies,lambda=" + formatted_lambda + ",alpha");
	if (!std::filesystem::exists(energies_file))
		energies_file = stdout_file.parent_path() / (jobname + ".MO_energies,lambda=" + formatted_lambda + ",a");

	std::filesystem::path orbitals_file = stdout_file.parent_path() / (jobname + ".molecular_orbitals,lambda=" + formatted_lambda + ",restricted");
	if (!std::filesystem::exists(orbitals_file))
		orbitals_file = stdout_file.parent_path() / (jobname + ".MOs,lambda=" + formatted_lambda + ",r");
	if (!std::filesystem::exists(orbitals_file))
		orbitals_file = stdout_file.parent_path() / (jobname + ".molecular_orbitals,lambda=" + formatted_lambda + ",alpha");
	if (!std::filesystem::exists(orbitals_file))
		orbitals_file = stdout_file.parent_path() / (jobname + ".MOs,lambda=" + formatted_lambda + ",a");

	while (std::filesystem::exists(energies_file) && std::filesystem::exists(orbitals_file)) {
		const std::filesystem::path of = orbitals_file;
		const std::filesystem::path ef = energies_file;

		err_checkf(of.string() != "", "Orbitals file name is empty?!", std::cout);
		err_checkf(ef.string() != "", "Energy file name is empty?!", std::cout);
		err_checkf(std::filesystem::exists(of), "couldn't open or find " + of.string() + ", leaving", std::cout);
		err_checkf(std::filesystem::exists(ef), "couldn't open or find " + ef.string() + ", leaving", std::cout);
		std::cout << "lambda = " + std::to_string(lambda) + "..." << std::flush;

		std::vector<WFN> wavy;
		wavy.emplace_back(e_origin::tonto);
		wavy.back().read_tonto(stdout_file, std::cout, debug, ef, of);
		std::filesystem::path basename = stdout_file.parent_path() / (jobname + "_l_" + formatted_lambda);
		wavy.back().write_wfn(basename.string() + ".wfn", debug, false);
		free_fchk(std::cout, basename.string() + ".fchk", "", wavy.back(), debug, true);

		if (opt.cif != "") {
			svec ka;
			int nr = 0;
			opt.groups[0] = { 0 };
			tsc_block<int, cdouble> result = calculate_scattering_factors<itsc_block, std::vector<WFN> &>(opt, wavy, std::cout, ka, nr);
			result.write_tscb_file(opt.cif, basename.string() + ".tscb");
		}

		lambda += ls;
		std::stringstream ss_l;
		ss_l << std::fixed << std::setprecision(6) << lambda;
		formatted_lambda = ss_l.str();
		energies_file = stdout_file.parent_path() / (jobname + ".orbital_energies,lambda=" + formatted_lambda + ",restricted");
		if (!std::filesystem::exists(energies_file))
			energies_file = stdout_file.parent_path() / (jobname + ".MO_energies,lambda=" + formatted_lambda + ",r");
		if (!std::filesystem::exists(energies_file))
			energies_file = stdout_file.parent_path() / (jobname + ".orbital_energies,lambda=" + formatted_lambda + ",alpha");
		if (!std::filesystem::exists(energies_file))
			energies_file = stdout_file.parent_path() / (jobname + ".MO_energies,lambda=" + formatted_lambda + ",a");

		orbitals_file = stdout_file.parent_path() / (jobname + ".molecular_orbitals,lambda=" + formatted_lambda + ",restricted");
		if (!std::filesystem::exists(orbitals_file))
			orbitals_file = stdout_file.parent_path() / (jobname + ".MOs,lambda=" + formatted_lambda + ",r");
		if (!std::filesystem::exists(orbitals_file))
			orbitals_file = stdout_file.parent_path() / (jobname + ".molecular_orbitals,lambda=" + formatted_lambda + ",alpha");
		if (!std::filesystem::exists(orbitals_file))
			orbitals_file = stdout_file.parent_path() / (jobname + ".MOs,lambda=" + formatted_lambda + ",a");
		std::cout << " .. done!" << std::endl;
	}
};
