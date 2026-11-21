/**
 * @file scattering_factors.cpp
 * @brief The tsc/tshb drivers (calculate_scattering_factors and its streamed parts), the ECP and electron-diffraction
 * conversions and the diffuse scattering factors. hkl and k-point I/O: hkl_io.cpp; CIF atoms: cif_atoms.cpp; calc_SF and
 * its kernels: sf_kernels.cpp; SALTED form factors: sf_salted.cpp.
 *
 */
#include "pch.h"
#include "tsc_block.h"
#include "tsc_stream.h"
#include "scattering_factors.h"
#include <map>
#include "convenience.h"
#include "cell.h"
#include "wfn_class.h"
#include "spherical_density.h"
#include "AtomGrid.h"
#include "integrator.h"
#include "basis_set.h"
#include "citations.h"
#ifdef NOSPHERA2_USE_GPU
#include "SALTED_equicomb.h"
#include "grid_gpu.h"
#include "blas_gpu.h"
#endif
#include "SALTED_utilities.h"
#include "gaussian_atom.h"
#include "GridManager.h"
#include "cube.h"
#include "tuning.h"
#ifdef NOSPHERA2_USE_GPU
#include "sf_gpu.h"
#include "itensor_gpu.h"
#endif


#ifdef __APPLE__
#include "TargetConditionals.h"
#endif

//Declared in the header so the XCW grid screens on the same ladder: what counts as
//negligible is the run's accuracy setting, and there should be one answer to that.
double cutoff(const int& accuracy)
{
	if (accuracy < 3)
		return 1E-10;
	else if (accuracy == 3)
		return 1E-14;
	else
		return 1E-30;
}

// This function yields the fourier bessel transform of the radial integral of a gaussian density function (compare equation 1.2.7.9 in 10.1107/97809553602060000759), assuming that H = 2 \pi S
double fourier_bessel_integral(
	const primitive& p,
	const double& H,
	const int& l)
{
	const double b = p.get_exp();
	return (pow(H, l) * exp(-H * H / (4.0 * b))) / (constants::pow_2[l] * p.get_exp_l_plus_3_2());
}

cdouble sfac_bessel(
	const primitive& p,
	const double* k_point,
	const double* coefs)
{
	const int l = p.get_type();
	switch (l % 4) {
	case 0:
		return cdouble(constants::PI3_2 * fourier_bessel_integral(p, k_point[3], l) * p.get_normalized_coefficient() * constants::spherical_harmonic(l, k_point, coefs), 0);
	case 1:
		return cdouble(0, constants::PI3_2 * fourier_bessel_integral(p, k_point[3], l) * p.get_normalized_coefficient() * constants::spherical_harmonic(l, k_point, coefs));
	case 2:
		return cdouble(-constants::PI3_2 * fourier_bessel_integral(p, k_point[3], l) * p.get_normalized_coefficient() * constants::spherical_harmonic(l, k_point, coefs), 0);
	case 3:
		return cdouble(0, -constants::PI3_2 * fourier_bessel_integral(p, k_point[3], l) * p.get_normalized_coefficient() * constants::spherical_harmonic(l, k_point, coefs));
	default:
		//should normally never happen
		return constants::cnull;
	}
}

/**
 * Adds the ECP (Effective Core Potential) contribution to the scattering factors.
 *
 * @param asym_atom_list The list of asymmetric atoms.
 * @param wave The WFN (Wave Function) object.
 * @param sf The scattering factors.
 * @param cell The cell object.
 * @param hkl The hkl_list object.
 * @param file The output stream for writing the results.
 * @param mode The mode of operation. 0 = Gaussian tight core function, 1,2,3 = Thakkar core density based on the ECP type used
 * @param debug Flag indicating whether to enable debug mode.
 */
//column s of sf belongs to stl[s], so a streamed block gets its core as readily as a whole table
static void add_ECP_contribution(const ivec& asym_atom_list,
	const WFN& wave,
	cvec2& sf,
	const vec& stl,
	std::ostream& file,
	const int& mode,
	const bool debug)
{
	using namespace std;
	double k = 1.0;
	err_checkf(mode >= 0, "Invalid mode for ECP contribution!", file);
	if (mode == 0)
	{ // Using a gaussian tight core function
		if (debug)
		{
			file << "Using a gaussian tight core function" << endl;
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				if (wave.get_atom_ECP_electrons(asym_atom_list[i]) != 0)
					file << "Atom nr: " << wave.get_atom_charge(asym_atom_list[i]) << " core f000: "
					<< scientific << setw(14) << setprecision(8)
					<< wave.get_atom_ECP_electrons(asym_atom_list[i])
					<< " and at 1 angstrom: " << exp(-pow(constants::bohr2ang(k), 2) / 16.0 / constants::PI) * wave.get_atom_ECP_electrons(asym_atom_list[i]) << "\n";
			}
		}
#pragma omp parallel for private(k)
		for (int s = 0; s < (int)stl.size(); s++)
		{
			k = constants::FOUR_PI * stl[s];
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				sf[i][s] += wave.get_atom_ECP_electrons(asym_atom_list[i]) * exp(-k / 16.0 / constants::PI);
			}
		}
	}
	else if (mode == 1 || mode == 2 || mode == 3)
	{ // Using a Thakkar core density
		if (debug)
			file << "Using a Thakkar core density" << endl;
		//vector<Thakkar> temp;
		map<int, Thakkar> temp;
		map<int, Spherical_Gaussian_Density> temp_G;
		if (debug) {
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				const int charge = wave.get_atom_charge(asym_atom_list[i]);
				if (temp.find(charge) == temp.end()) {
					temp.emplace(charge, Thakkar{ charge, mode });
					temp_G.emplace(charge, Spherical_Gaussian_Density{ charge, mode });
					if (wave.get_atom_ECP_electrons(asym_atom_list[i]) != 0)
					{
						double k_0001 = temp[charge].get_core_form_factor(0, wave.get_atom_ECP_electrons(asym_atom_list[i]));
						double k_1 = temp[charge].get_core_form_factor(constants::FOUR_PI * constants::bohr2ang(1.0), wave.get_atom_ECP_electrons(asym_atom_list[i]));
						file << "Atom nr: " << charge << " number of ECP electrons: " << wave.get_atom_ECP_electrons(asym_atom_list[i]) << " core f(0) : "
							<< scientific << setw(14) << setprecision(8) << k_0001 << " and at 1 Ang: " << k_1 << "\n";
					}
				}
			}
		}
		else {
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				const int charge = wave.get_atom_charge(asym_atom_list[i]);
				if (temp.find(charge) == temp.end()) {
					temp.emplace(charge, Thakkar{ charge, mode });
					temp_G.emplace(charge, Spherical_Gaussian_Density{ charge, mode });
				}
			}
		}

#pragma omp parallel for private(k)
		for (int s = 0; s < (int)stl.size(); s++)
		{
			int n_el = 0, charge = 0;
			k = constants::FOUR_PI * constants::bohr2ang(stl[s]);
			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				n_el = wave.get_atom_ECP_electrons(asym_atom_list[i]);
				charge = wave.get_atom_charge(asym_atom_list[i]);
				if (n_el != 0)
				{
					sf[i][s] += temp_G.at(charge).get_form_factor(k); // This bit will correct for the error of the valence denisty of ECP atoms
					sf[i][s] += temp.at(charge).get_core_form_factor(k, n_el);
				}
			}
		}
	}
	else
	{
		err_not_impl_f("No higher ECP mode than 3 implemented!", file);
	}
}

static void add_ECP_contribution(const ivec& asym_atom_list,
	const WFN& wave,
	cvec2& sf,
	const cell& cell,
	hkl_list& hkl,
	std::ostream& file,
	const int& mode,
	const bool debug)
{
	vec stl;
	stl.reserve(hkl.size());
	for (const auto& h : hkl)
		stl.push_back(cell.get_stl_of_hkl(h));
	add_ECP_contribution(asym_atom_list, wave, sf, stl, file, mode, debug);
}

/**
 * Converts the given asymmetric atom list, wavefunction, unit cell, and hkl list
 * to electron scattering factors using the given X-ray scattering factors.
 * This will be performed in place
 *
 * @param asym_atom_list The list of asymmetric atoms.
 * @param wave The wavefunction.
 * @param sf The scattering factors.
 * @param unit_cell The unit cell.
 * @param hkl The hkl list.
 */
//sf[i][s] depends on stl[s] and atom i alone, so this applies to a block as readily as to a whole table
void convert_to_ED(const ivec& asym_atom_list,
	const WFN& wave,
	cvec2& sf,
	const vec& stl)
{
	const int n = (int)stl.size();
#pragma omp parallel for
	for (int s = 0; s < n; s++)
	{
		const double h2 = pow(stl[s], 2);
		for (int i = 0; i < asym_atom_list.size(); i++)
			sf[i][s] = cdouble(constants::ED_fact * (wave.get_atom_charge(asym_atom_list[i]) - sf[i][s].real()) / h2, -constants::ED_fact * sf[i][s].imag() / h2);
	}
}

void convert_to_ED(const ivec& asym_atom_list,
	const WFN& wave,
	cvec2& sf,
	const cell& unit_cell,
	const std::vector<i3>& hkl_vector)
{
	vec stl(hkl_vector.size());
	for (size_t s = 0; s < hkl_vector.size(); s++)
		stl[s] = unit_cell.get_stl_of_hkl(hkl_vector[s]);
	convert_to_ED(asym_atom_list, wave, sf, stl);
}

void convert_to_ED(const ivec& asym_atom_list,
	const WFN& wave,
	cvec2& sf,
	const cell& unit_cell,
	const hkl_list& hkl)
{
	convert_to_ED(asym_atom_list, wave, sf, unit_cell,
		std::vector<i3>(hkl.begin(), hkl.end()));
}

//shared machinery of the streamed tsc writes; the fill of one reflection block is all that differs
namespace {

//a table must hold either all labels or all atomIDs, mixed ones are rejected on write
void append_scatterer_ids(ScattererLabels& ids,
	const options& opt,
	const svec& labels,
	WFN& wave,                 // non-const: get_id_for_atom() is not a const method
	const ivec& asym_atom_list)
{
	for (size_t a = 0; a < asym_atom_list.size(); a++)
	{
		if (opt.label_tsc_output)
			ids.emplace_back(labels[a]);
		else
			ids.emplace_back(wave.get_id_for_atom(asym_atom_list[a]));
	}
}

//reflections [lo, hi) transposed to [dimension][reflection], the layout write_tscb_reflection_block expects
ivec2 slice_hkl(const std::vector<i3>& hkl_v,
	const size_t lo, const size_t hi)
{
	ivec2 idx(3, ivec(hi - lo));
	for (size_t r = lo; r < hi; r++)
		for (int dm = 0; dm < 3; dm++)
			idx[dm][r - lo] = hkl_v[r][dm];
	return idx;
}

//calc_SF_SALTED sizes its output from the k-points it is handed, so a slice needs no change on its side
vec2 slice_k_points(const vec2& k_pt, const size_t lo, const size_t hi)
{
	vec2 k_slice(3, vec(hi - lo));
	for (size_t r = lo; r < hi; r++)
		for (int dm = 0; dm < 3; dm++)
			k_slice[dm][r - lo] = k_pt[dm][r];
	return k_slice;
}

//one Thakkar evaluator per element; depends on the element list alone, so it is built outside the block loop
std::vector<Thakkar> make_spherical_evaluators(const ivec& atom_type_list)
{
	std::vector<Thakkar> spheres;
	spheres.reserve(atom_type_list.size());
	for (size_t t = 0; t < atom_type_list.size(); t++)
		spheres.emplace_back(atom_type_list[t]);
	return spheres;
}

//the file list a spherical fill runs over, and the index into it; nr also selects opt.groups[nr] and the part's cif,
//so it must be the whole file list: a one-element vector takes part 0's atoms whatever part is filled, or reads tempy[1] if nr is passed;
//only tempy[nr] is read, the other parts stay empty
int build_fill_wavefunctions(const options& opt, const int nr, std::vector<WFN>& tempy)
{
	if (!opt.wfn.empty())
	{
		tempy.emplace_back(opt.wfn);
		return 0;
	}
	for (int i = 0; i < (int)opt.combined_tsc_calc_files.size(); i++)
		if (i == nr)
			tempy.emplace_back(opt.combined_tsc_calc_files[i]);
		else
			tempy.emplace_back();
	return nr;
}

//One evaluator per atom: same-element atoms carry different EEQ charges. Charges are matched by position, as the
//fill rebuilds its wavefunction from the original file; no recorded charge or no tabulated ion falls back to neutral.
std::vector<HE_Spherical_Atom> make_he_evaluators(const salted_part_prep& sph,
	const WFN& fill_wavy, const options& opt, std::ostream& file)
{
	std::vector<HE_Spherical_Atom> out;
	out.reserve(sph.asym_atom_list.size());
	int n_charged = 0, n_extrap = 0, n_nofit = 0, n_delta = 0;
	double worst_neg = 0.0, q_min = 1e300, q_max = -1e300;
	for (size_t a = 0; a < sph.asym_atom_list.size(); a++)
	{
		const int idx = sph.asym_atom_list[a];
		const int Z = fill_wavy.get_atom_charge(idx);
		double q = 0.0;
		double best = 1e-3;   // squared tolerance in the wavefunction's units
		for (const auto& e : opt.spherical_fill_charges)
		{
			double d = 0.0;
			for (int ax = 0; ax < 3; ax++)
			{
				const double dx = fill_wavy.get_atom_coordinate(idx, ax) - e[ax];
				d += dx * dx;
			}
			if (d < best) { best = d; q = e[3]; }
		}
		if (std::abs(q) > 1e-6) { n_charged++; q_min = std::min(q_min, q); q_max = std::max(q_max, q); }
		else if (!opt.spherical_fill_charges.empty()) n_nofit++;
		out.emplace_back(Z, q);
		if (out.back().is_extrapolating()) n_extrap++;
		if (out.back().uses_delta_series()) n_delta++;
		worst_neg = std::min(worst_neg, out.back().most_negative_density());
	}
	if (n_charged)
	{
		file << "Spherical fill: " << n_charged << " of " << out.size()
			 << " atom(s) given a fractional charge (" << q_min << " to " << q_max << ")";
		if (n_delta)
			file << ", " << n_delta << " past +1 via the delta_k series (interpolated between bound states)";
		if (n_extrap)
			file << ", " << n_extrap << " beyond the tabulated ions (shape extrapolated)";
		file << "." << std::endl;
		if (worst_neg < -1e-8)
			file << "      most negative density from the blend: " << worst_neg
				 << " e/bohr^3 - an extrapolation artefact, not clipped." << std::endl;
		if (n_nofit)
			file << "      " << n_nofit << " atom(s) had no charge recorded and stay neutral." << std::endl;
	}
	return out;
}

//the filled atoms of a part into the run's .tshb: the whole atom (every occupied orbital) as Slater records in the charge
//state make_he_evaluators gives the streamed table, Z on each, and the fitted delta_k or GTO differences as Gaussians
void fill_tshb(const salted_part_prep& sph, WFN& fill_wavy, options& opt, std::ostream& file)
{
	if (sph.asym_atom_list.empty())
		return;
	const std::vector<HE_Spherical_Atom> spheres = make_he_evaluators(sph, fill_wavy, opt, file);
	ScattererLabels keys;
	append_scatterer_ids(keys, opt, sph.labels, fill_wavy, sph.asym_atom_list);
	std::vector<tshb_record> records(sph.asym_atom_list.size());
	std::vector<std::array<double, 3>> slaters;
	std::vector<std::array<double, 2>> gaussians;
	const double b = constants::bohr2ang(1.0);
	for (size_t a = 0; a < records.size(); a++)
	{
		const int idx = sph.asym_atom_list[a];
		tshb_record& r = records[a];
		r.label = sph.labels[a];
		r.part = fill_wavy.get_atom_group_nr(idx);
		r.Z = fill_wavy.get_atom_charge(idx);
		spheres[a].terms(slaters, gaussians);
		for (const std::array<double, 3>& t : slaters)
			r.slaters.push_back(tshb_slater(t));
		for (const std::array<double, 2>& g : gaussians)
			r.gaussians.push_back({ g[0] / (b * b * b), g[1] / (b * b) });
	}
	opt.tshb->append(tshb_block(std::move(records), keys, opt.electron_diffraction), file);
	file << "Spherical fill: " << sph.asym_atom_list.size() << " atom(s) go into the .tshb as Slater records" << std::endl;
}

//rows for atoms the SALTED model could not predict (unknown species, or nothing inside the descriptor cutoff)
//a Thakkar factor depends only on element and reflection, so these chunk like everything else here
void append_spherical_rows(cvec2& chunk,
	const salted_part_prep& sph,
	const std::vector<HE_Spherical_Atom>& spheres,
	const bool electron_diffraction,
	const size_t lo, const size_t hi)
{
	const size_t first = chunk.size();
	const size_t n_atoms = sph.asym_atom_list.size();
	chunk.resize(first + n_atoms, cvec(hi - lo));
	//every value is independent, so the threads split the reflections; this ran serially and outlasted the prediction
#pragma omp parallel for
	for (int s = 0; s < (int)(hi - lo); s++)
		for (size_t a = 0; a < n_atoms; a++)
		{
			const int t = sph.asym_atom_to_type_list[a];
			const double f = spheres[a].get_form_factor(sph.k_of_reflection[lo + s]);
			//IAM form of Mott-Bethe: tabulated charge, no imaginary part
			chunk[first + a][s] = electron_diffraction
				? cdouble(constants::ED_fact * (sph.atom_type_list[t] - f) /
					pow(sph.stl_of_reflection[lo + s], 2), 0.0)
				: cdouble(f, 0.0);
		}
}

//writes a tsc as a sequence of reflection blocks; peak memory is queue depth * scatterers * block size * 16 bytes
//fill_block(lo, hi, progress) returns the rows for reflections [lo, hi) in the order `ids` declares them
//progress_items is the bar total: reflections, or reflections * parts when a block is filled part by part
template <typename FillBlock>
void stream_blocks(options& opt,
	std::ostream& file,
	const std::filesystem::path& name,
	const ScattererLabels& ids,
	const std::vector<i3>& hkl_v,
	const size_t progress_items,
	FillBlock fill_block)
{
	const size_t n_refl = hkl_v.size();
	const size_t derived = opt.tsc_block_for(n_refl, ids.size());
	const size_t bs = std::min(derived ? derived : n_refl, n_refl ? n_refl : 1);
	file << "Streaming tsc in blocks of " << bs << " reflections" << std::endl;

	tsc_stream_writer<int, cdouble> writer(name, ids, std::string(), n_refl, 2);
	//declared after the writer so it is destroyed first: the bar rewinds to its own line and must finish first
	ProgressBar progress(progress_items, 60, "#", " ", "Generating scattering factors...");
	size_t block_id = 0;
	for (size_t lo = 0; lo < n_refl; lo += bs)
	{
		const size_t hi = std::min(lo + bs, n_refl);
		writer.submit(block_id++, slice_hkl(hkl_v, lo, hi), fill_block(lo, hi, progress));
	}
	writer.finish();
	opt.tsc_written_by_stream = true;
}

//runs a nested calculate_scattering_factors() as a spherical (Thakkar) fill over exactly the reflections written
//restores every flag on the way out, exceptions included: these flags are read all over this file
struct spherical_fill_scope
{
	options& opt;
	const bool saved_iam;
	const hkl_list saved_hkl;

	spherical_fill_scope(options& o, const hkl_list& reflections)
		: opt(o), saved_iam(o.iam_switch), saved_hkl(o.m_hkl_list)
	{
		//pin the fill to our reflections, else it builds its own list of a different length and is read off the end
		opt.m_hkl_list = reflections;
		opt.iam_switch = true;
		opt.allow_empty_asym = true;
		opt.spherical_fill = true;
	}
	~spherical_fill_scope()
	{
		opt.m_hkl_list = saved_hkl;
		opt.iam_switch = saved_iam;
		opt.allow_empty_asym = false;
		opt.spherical_fill = false;
	}

	spherical_fill_scope(const spherical_fill_scope&) = delete;
	spherical_fill_scope& operator=(const spherical_fill_scope&) = delete;
};

}


int make_atomic_grids_wrapper(
	const WFN& wave, const bvec& needs_grid, const ivec& asym_atom_list, const cell& unit_cell, const svec& labels, //
	std::vector<_time_point>& time_points, svec& time_descriptions, vec2& d1, vec2& d2, vec2& d3, vec2& dens,
	const options& opt, std::ostream& file = std::cout, sf_shell_layout* layout = nullptr) {

	citations::cite(citations::Method::BeckeGrid, file);
	if (opt.partition_type == PartitionType::Hirshfeld)
		citations::cite(citations::Method::Hirshfeld, file);
	else if (opt.partition_type == PartitionType::TFVC)
		citations::cite(citations::Method::TFVC, file);
	else if (opt.partition_type == PartitionType::MBIS || opt.partition_type == PartitionType::EMBIS) {
		citations::cite(citations::Method::MBIS, file);
		//EMBIS is MBIS with an ellipsoidal sigma: both papers
		if (opt.partition_type == PartitionType::EMBIS)
			citations::cite(citations::Method::EMBIS, file);
	}

	const int atoms_with_grids = vec_sum(needs_grid);
	err_checkf(atoms_with_grids > 0, "No atoms with grids to generate!", file);
	err_checkf(atoms_with_grids <= wave.get_ncen(), "More atoms with grids than in the wavefunction! Aborting!", file);
	err_checkf(atoms_with_grids == asym_atom_list.size(), "Number of atoms with grids does not match the number of atoms in the CIF file!", file);
	file << "There are:\n"
		<< std::setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
		//<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
		<< std::setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << std::endl;


	file << "\nSelected accuracy: " << opt.accuracy << "\nMaking Integration Grids..." << std::endl;

	GridConfiguration config;
#ifdef NOSPHERA2_USE_GPU
	grid_gpu_set_enabled(opt.use_gpu && opt.gpu_grid);
	blas_gpu_set_enabled(opt.gpu_blas);
#endif
	config.accuracy = opt.accuracy;
	config.partition_type = opt.partition_type;
	config.pbc = opt.pbc;
	config.debug = opt.debug;
	config.all_charges = opt.all_charges;

	GridManager local_manager(config);
	GridManager& grid_manager = opt.grid_cache ? *opt.grid_cache : local_manager;
	if (opt.grid_cache)
		grid_manager.setConfiguration(config);

	WFN temp = wave;
	temp.delete_unoccupied_MOs();

	// Setup grids for the molecule
	grid_manager.setup3DGridsForMolecule(temp, asym_atom_list, needs_grid, unit_cell, opt.get_g, file);
	grid_manager.addTimingInfoToVecs(time_points, time_descriptions);


	// Calculate partitioned charges
	PartitionResults results = grid_manager.calculatePartitionedCharges(temp, unit_cell);
	grid_manager.printChargeTable(labels, temp, asym_atom_list, file, results);
	time_points.push_back(get_time());
	time_descriptions.push_back("calculate charges");

	grid_manager.getDensityVectors(temp, asym_atom_list, d1, d2, d3, dens, opt.get_g, layout);
	time_points.push_back(get_time());
	time_descriptions.push_back("combined density vectors");

	return grid_manager.getTotalGridPoints();
}

itsc_block calculate_scattering_factors_from_cube(
	options& opt,
	WFN& wave,
	const cube& density_cube,
	std::ostream& file)
{
	using namespace std;
	err_checkf(opt.cif != "", "Cube-density SF calculation requires -cif.", file);
	err_checkf(filesystem::exists(opt.cif), "CIF " + opt.cif.string() + " does not exists!", file);
	err_checkf(!opt.groups.empty() && !opt.groups[0].empty(), "No CIF disorder group selected. Use -group 0 for the default group.", file);
	err_checkf(wave.get_ncen() != 0, "Cube file does not contain atom definitions in the header.", file);

	vector<_time_point> time_points;
	svec time_descriptions;
	time_points.push_back(get_time());
#ifdef NOSPHERA2_USE_GPU
	if (opt.use_gpu)
		sf_gpu_warmup_start();
#endif

	cell unit_cell(opt.cif, file, opt.debug);
	unit_cell.set_hkl_scale(opt.hkl_scale);
	ifstream cif_input(opt.cif.c_str(), ios::in);
	ivec atom_type_list;
	ivec asym_atom_to_type_list;
	ivec asym_atom_list;
	bvec needs_grid(wave.get_ncen(), false);
	svec known_atoms;

	auto labels = read_atoms_from_CIF(cif_input,
		opt.groups[0],
		unit_cell,
		wave,
		known_atoms,
		atom_type_list,
		asym_atom_to_type_list,
		asym_atom_list,
		needs_grid,
		file,
		opt.debug);

	cif_input.close();
	time_points.push_back(get_time());
	time_descriptions.push_back("cif reading");

	vec2 k_pt;
	hkl_list hkl;
	if (opt.m_hkl_list.size() != 0)
	{
		hkl = opt.m_hkl_list;
	}
	else if (opt.read_k_pts == false)
	{
		generate_hkl_from_options(opt, hkl, unit_cell, file);
		opt.m_hkl_list = hkl;
	}

	make_k_pts(
		hkl.size() == 0,
		opt.save_k_pts,
		unit_cell,
		hkl,
		k_pt,
		file,
		opt.debug);

	time_points.push_back(get_time());
	time_descriptions.push_back("k-points preparation");

	GridConfiguration config;
#ifdef NOSPHERA2_USE_GPU
	grid_gpu_set_enabled(opt.use_gpu && opt.gpu_grid);
	blas_gpu_set_enabled(opt.gpu_blas);
#endif
	config.accuracy = opt.accuracy;
	config.partition_type = opt.partition_type;
	config.pbc = opt.pbc;
	config.debug = opt.debug;
	config.all_charges = opt.all_charges;
	config.no_density_eval = true;

	GridManager grid_manager(config);
	grid_manager.setup3DGridsForMolecule(wave, asym_atom_list, needs_grid, unit_cell, false);
	grid_manager.addTimingInfoToVecs(time_points, time_descriptions);

	vec2 d1, d2, d3, dens;
	vec atom_electrons;
	grid_manager.getDensityVectorsFromCube(wave, asym_atom_list, density_cube, d1, d2, d3, dens, atom_electrons);
	time_points.push_back(get_time());
	time_descriptions.push_back("combined density vectors");

	file << "Table of Charges in electrons\n"
		<< setw(10) << "Atom" << setw(12) << "CubePart" << endl;
	for (int i = 0; i < labels.size(); i++)
	{
		const int atom_index = asym_atom_list[i];
		file << setw(10) << labels[i]
			<< fixed << setw(12) << setprecision(3) << wave.get_atom_charge(atom_index) - atom_electrons[i] << "\n";
	}
	const double electron_sum = reduce(atom_electrons.begin(), atom_electrons.end(), 0.0);
	file << setprecision(4) << "Total number of partitioned electrons from cube: " << electron_sum << endl;
	time_points.push_back(get_time());
	time_descriptions.push_back("calculate charges");

	cvec2 sf;
	_time_point end1;
	calc_SF(grid_manager.getTotalGridPoints(),
		k_pt,
		d1, d2, d3, dens,
		sf,
		file,
		time_points.front(),
		end1,
		opt.debug,
		opt.no_date,
		false,
		opt.use_gpu,
		opt.gpu_fp64,
		opt.gpu_fp32);
	time_points.push_back(get_time());
	time_descriptions.push_back("Fourier transform");

	if (opt.electron_diffraction)
	{
		convert_to_ED(asym_atom_list, wave, sf, unit_cell, hkl);
	}

	itsc_block blocky;
	if (opt.label_tsc_output)
	{
		blocky = itsc_block(sf, labels, hkl);
	}
	else
	{
		vector<atomID> ids(asym_atom_list.size());
		for (int atm_idx = 0; atm_idx < asym_atom_list.size(); atm_idx++)
		{
			ids[atm_idx] = wave.get_id_for_atom(asym_atom_list[atm_idx]);
		}
		blocky = itsc_block(sf, ids, hkl);
	}
	time_points.push_back(get_time());
	time_descriptions.push_back("tsc calculation");

	if (!opt.no_date)
	{
		write_timing_to_file(file, time_points, time_descriptions);
	}

	return blocky;
}

/**
 * @brief Calculates scattering factors using the provided calculator and options.
 *
 * This templated function computes scattering factors for a set of atoms and k-points.
 *
 * @tparam tsc_block_type The type used to store the resulting scattering factors block.
 *         Must support assignment and access as used in the function.
 * @tparam calculator_type The type of the calculator object.
 *         Must provide the required interface for calculating scattering factors.
 * @param opt Options for the calculation.
 * @param calculator The calculator object used to compute scattering factors.
 * @param file Output stream for logging messages.
 * @param known_atoms List of known atom symbols.
 * @param nr Number of atoms.
 * @param kpts Pointer to the k-points vector.
 * @return A block of scattering factors of type tsc_block_type.
 */
template <typename tsc_block_type, typename calculator_type>
tsc_block_type calculate_scattering_factors(
	options& opt,
	calculator_type calculator,
	std::ostream& file,
	svec& known_atoms,
	const int& nr,
	vec2* kpts,
	salted_part_prep* prep_out
) {
	using namespace std;
#ifdef NOSPHERA2_USE_GPU
	if (opt.use_gpu)
		sf_gpu_warmup_start();
#endif
	int nat = 0;
	WFN* wavy = NULL;
	if constexpr (std::is_same_v<calculator_type, std::vector<WFN> &>) {
		wavy = &calculator[nr];
		wavy->delete_Qs();
		err_checkf(wavy->get_ncen() != 0, "No Atoms in the wavefunction, this will not work!!ABORTING!!", file);
		if (!opt.cif_based_combined_tsc_calc)
		{
			err_checkf(filesystem::exists(opt.cif), "CIF " + opt.cif.string() + " does not exists!", file);
		}
		else
		{
			for (int i = 0; i < opt.combined_tsc_calc_cifs.size(); i++)
			{
				err_checkf(filesystem::exists(opt.combined_tsc_calc_cifs[i]), "CIF " + opt.combined_tsc_calc_cifs[i].string() + " does not exists!", file);
			}
		}
		err_checkf(opt.groups[nr].size() >= 1, "Not enough groups specified to work with!", file);
		citations::cite(citations::Method::NoSpherA2, file);
		if (opt.iam_switch)
			citations::cite(citations::Method::IAM, file);
		else
			citations::cite(citations::Method::HAR, file);
		file << "Number of protons: " << wavy->get_nr_electrons() << endl
			<< "Number of electrons: " << fixed << wavy->count_nr_electrons() << endl;
		if (wavy->get_has_ECPs()) {
			file << "Number of ECP electrons: " << wavy->get_nr_ECP_electrons() << endl;
			citations::cite(citations::Method::ECP, file);
		}
		// err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);
		if (opt.debug)
			file << "Working with: " << wavy->get_path() << endl;
		nat = wavy->get_ncen();
	}
	else if constexpr (std::is_same_v<calculator_type, SALTEDPredictor&>) {
		wavy = &calculator.wavy;
		nat = wavy->get_ncen();
	}

	vector<_time_point> time_points;
	svec time_descriptions;
	time_points.push_back(get_time());

	filesystem::path cif;
	if (opt.cif_based_combined_tsc_calc)
		cif = opt.combined_tsc_calc_cifs[nr];
	else
		cif = opt.cif;


	cell unit_cell(cif, file, opt.debug);
	unit_cell.set_hkl_scale(opt.hkl_scale);
	ifstream cif_input(cif.c_str(), std::ios::in);
	ivec atom_type_list;
	ivec asym_atom_to_type_list;
	ivec asym_atom_list;
	bvec needs_grid(nat, false);
	if (opt.debug)
		file << "Reading atoms!!!!" << endl;

	auto labels = read_atoms_from_CIF(cif_input,
		opt.groups[nr],
		unit_cell,
		*wavy,
		known_atoms,
		atom_type_list,
		asym_atom_to_type_list,
		asym_atom_list,
		needs_grid,
		file,
		opt.debug,
		//a -mtc part whose new atoms the model all dropped (a lone metal or isolated ion in its own PART) predicts
		//nothing; stream_mtc_salted gives it only its spherical rows
		opt.allow_empty_asym || (prep_out && opt.needs_Thakkar_fill && std::is_same_v<calculator_type, SALTEDPredictor&>));

	cif_input.close();

	vec2 k_pt;
	hkl_list hkl;
	if (opt.m_hkl_list.size() != 0)
	{
		hkl = opt.m_hkl_list;
	}
	else if (nr == 0 && opt.read_k_pts == false)
	{
		generate_hkl_from_options(opt, hkl, unit_cell, file);
		opt.m_hkl_list = hkl;
	}

	//empty only means a broken CIF unless the caller allowed it: a spherical fill of an already covered part finds nothing.
	//After the hkl list, which the later parts take from part 0 even when part 0 is the empty one (COD 2018514)
	if (asym_atom_list.empty())
	{
		if (prep_out) *prep_out = salted_part_prep();
		return tsc_block_type();
	}

	if (opt.debug)
		file << "There are " << atom_type_list.size() << " Types of atoms and " << asym_atom_to_type_list.size() << " atoms in total" << endl;

	time_points.push_back(get_time());
	time_descriptions.push_back("cif reading");

	//a part whose atoms go into the .tshb, and the fill's preparation for it, write records and evaluate no reflection
	const bool records_only = opt.tshb != nullptr && (opt.spherical_fill ? prep_out != NULL : prep_out == NULL);
	if (records_only)
		;
	else if (kpts == NULL || kpts->size() == 0)
	{
		make_k_pts(
			nr != 0 && hkl.size() == 0,
			opt.save_k_pts,
			unit_cell,
			hkl,
			k_pt,
			file,
			opt.debug);
		if (kpts != NULL)
		{
			*kpts = k_pt;
		}
	}
	else
	{
		k_pt = *kpts;
	}

	time_points.push_back(get_time());
	time_descriptions.push_back("k-points preparation");
	//streaming emits the table in reflection blocks instead of holding scatterers * reflections * 16 bytes at once
	//excluded: prep_out and spherical_fill must hand their rows back rather than write a file, and a combined
	//table needs every part present - stream_mtc_salted() at the bottom inverts those loops instead
	//tsc_block_for returns 0 when -mem says the table fits whole; not named tsc_block, that is the class template below
	const size_t block_reflections = opt.tsc_block_for(hkl.size(), asym_atom_list.size());
	//the stream writer only produces the binary format, so a text table (-old_tsc or no tscb) keeps the whole table in memory
	const bool stream_tsc = block_reflections > 0
		&& prep_out == NULL
		&& !opt.spherical_fill
		&& opt.combined_tsc_calc_files.size() <= 1
		&& opt.binary_tsc && !opt.old_tsc;
	if (block_reflections > 0 && !stream_tsc && (!opt.binary_tsc || opt.old_tsc))
		file << "Text tsc requested: the table is kept in memory instead of streamed in blocks" << std::endl;
	//-tshb: the fitted densities go into the run's .tshb (opt.tshb, set by the caller once the run qualifies) and no
	//table is made; olex2.refine evaluates them at any reflection. The fill call and the -mtc preparation never decide it
	const bool tshb_only = opt.tshb != nullptr && !opt.iam_switch
		&& (std::is_same_v<calculator_type, SALTEDPredictor &> || opt.partition_type == PartitionType::RI)
		&& prep_out == NULL && !opt.spherical_fill;
	cvec2 sf;
	if (!stream_tsc && !tshb_only)
	{
		sf.resize(asym_atom_list.size());
#pragma omp parallel for
		for (int i = 0; i < asym_atom_list.size(); i++)
			sf[i].resize(hkl.size());
	}

	if (opt.iam_switch) {
		if (prep_out != NULL)
		{
			// hand back what a per-block spherical calculation needs and stop
			const std::vector<i3> hv(records_only ? std::vector<i3>() : std::vector<i3>(hkl.begin(), hkl.end()));
			prep_out->asym_atom_list = asym_atom_list;
			prep_out->labels = labels;
			prep_out->atom_type_list = atom_type_list;
			prep_out->asym_atom_to_type_list = asym_atom_to_type_list;
			prep_out->hkl_v = hv;
			prep_out->k_of_reflection.resize(hv.size());
			prep_out->stl_of_reflection.resize(hv.size());
			for (size_t s = 0; s < hv.size(); s++)
			{
				const double stl = unit_cell.get_stl_of_hkl(hv[s]);
				prep_out->stl_of_reflection[s] = stl;
				prep_out->k_of_reflection[s] = constants::bohr2ang(constants::FOUR_PI * stl);
			}
			return tsc_block_type();
		}
		//a part without orbitals in a -tshb run: its atoms as the fill's spherical records, no table
		if (opt.tshb != nullptr && !opt.spherical_fill)
		{
			salted_part_prep sph;
			sph.asym_atom_list = asym_atom_list;
			sph.labels = labels;
			sph.atom_type_list = atom_type_list;
			sph.asym_atom_to_type_list = asym_atom_to_type_list;
			fill_tshb(sph, *wavy, opt, file);
			return tsc_block_type();
		}
		std::vector<Thakkar> spherical_atoms = make_spherical_evaluators(atom_type_list);
		const int imax = (int)asym_atom_list.size();
		const std::vector<i3> hkl_vector(hkl.begin(), hkl.end());
		const int hkl_max = hkl.size();

		if (stream_tsc)
		{
			//a Thakkar factor depends on element and reflection alone, so no state crosses a block boundary
			ScattererLabels stream_ids;
			append_scatterer_ids(stream_ids, opt, labels, *wavy, asym_atom_list);

			stream_blocks(opt, file, "experimental.tscb",
				stream_ids, hkl_vector, hkl_vector.size(),
				[&](const size_t lo, const size_t hi, ProgressBar& progress)
				{
					cvec2 chunk(imax, cvec(hi - lo));
#pragma omp parallel for
					for (int s = 0; s < (int)(hi - lo); s++)
					{
						const double stl = unit_cell.get_stl_of_hkl(hkl_vector[lo + s]);
						const double k = constants::bohr2ang(constants::FOUR_PI * stl);
						const double h2 = pow(stl, 2);
						//once per element, not once per atom
						vec f_type(spherical_atoms.size());
						for (size_t t = 0; t < spherical_atoms.size(); t++)
							f_type[t] = spherical_atoms[t].get_form_factor(k);
						for (int i = 0; i < imax; i++)
						{
							const int type = asym_atom_to_type_list[i];
							const double f = f_type[type];
							//IAM form of Mott-Bethe: tabulated charge, no imaginary part
							chunk[i][s] = opt.electron_diffraction
								? cdouble(constants::ED_fact * (atom_type_list[type] - f) / h2, 0.0)
								: cdouble(f, 0.0);
						}
					}
					progress.update(hi - lo);
					return chunk;
				});
		}
		else
		{
#pragma omp parallel for shared(hkl_vector)
			for (int s = 0; s < hkl_max; s++)
			{
				const double stl = unit_cell.get_stl_of_hkl(hkl_vector[s]);
				const double k = constants::bohr2ang(constants::FOUR_PI * stl);
				const double h2 = pow(stl, 2);
				vec f_type(spherical_atoms.size());
				for (size_t t = 0; t < spherical_atoms.size(); t++)
					f_type[t] = spherical_atoms[t].get_form_factor(k);
				for (int i = 0; i < imax; i++)
				{
					const int type = asym_atom_to_type_list[i];
					const double f = f_type[type];
					//IAM form of Mott-Bethe: tabulated charge, no imaginary part
					sf[i][s] = opt.electron_diffraction
						? cdouble(constants::ED_fact * (atom_type_list[type] - f) / h2, 0.0)
						: cdouble(f, 0.0);
				}
			}
		}
	}
	else if constexpr (std::is_same_v<calculator_type, SALTEDPredictor &>)
	{
		// Generation of SALTED density coefficients
		file << "\nGenerating densities... " << endl;
#ifdef NOSPHERA2_USE_GPU
		equicomb_set_gpu(opt.use_gpu && opt.gpu_salted);
#endif
		auto ml = std::make_shared<Gaussian_Molecule>(calculator.wavy, calculator.gen_SALTED_densities());
		file << setw(13 * 4) << "... done!" << endl;
		time_points.push_back(get_time());
		time_descriptions.push_back("SALTED prediction");
		calculator.shrink_intermediate_vectors();

		err_checkf(labels.size() == asym_atom_list.size(),
			"Inconsistent SALTED atom bookkeeping after disorder filtering!", file);

		vec atom_elecs = ml->populations();
		file << "Table of Charges in electrons\n"
			<< "       Atom      ML" << endl;

		for (int i = 0; i < labels.size(); i++)
		{
			const int atom_index = asym_atom_list[i];
			file << setw(10) << labels[i]
				<< fixed << setw(10) << setprecision(3) << wavy->get_atom_charge(atom_index) - atom_elecs[atom_index];
			if (opt.debug)
				file << " " << setw(4) << wavy->get_atom_charge(atom_index) << " " << fixed << setw(10) << setprecision(3) << atom_elecs[atom_index];
			file << "\n";
		}
		auto el_sum = reduce(atom_elecs.begin(), atom_elecs.end(), 0.0);
		file << setprecision(4) << "Total number of analytical Electrons: " << el_sum << endl;
		time_points.push_back(get_time());
		time_descriptions.push_back("Calculation of Charges");
		if (tshb_only) {
			//keyed as the table's rows, the core as the table adds it: the model's ECP mode
			const auto t_rec = std::chrono::steady_clock::now();
			ScattererLabels keys;
			append_scatterer_ids(keys, opt, labels, *wavy, asym_atom_list);
			opt.tshb->append(ml->tshb(labels, keys, asym_atom_list, *wavy, calculator.wavy.get_ECP_mode(), opt.electron_diffraction), file);
			file << labels.size() << " atom(s) go into the .tshb, no table" << endl;
			if (ProgressBar::report_counts)
				std::cout << "[stages] tshb records " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_rec).count() << " s" << endl;
		}

		if (prep_out != NULL)
		{
			//-mtc streaming: the reflection loop lives outside, so hand back the reflection-independent part and stop
			prep_out->mol = ml;
			prep_out->asym_atom_list = asym_atom_list;
			prep_out->labels = labels;
			prep_out->k_pt = k_pt;
			prep_out->hkl_v.assign(hkl.begin(), hkl.end());
			//carried so the -mtc loop can convert to ED per block without a unit cell of its own
			prep_out->stl_of_reflection.resize(prep_out->hkl_v.size());
			for (size_t s = 0; s < prep_out->hkl_v.size(); s++)
				prep_out->stl_of_reflection[s] = unit_cell.get_stl_of_hkl(prep_out->hkl_v[s]);
			return tsc_block_type();
		}

		if (tshb_only)
			;
		else if (stream_tsc)
		{
			ScattererLabels stream_ids;
			append_scatterer_ids(stream_ids, opt, labels, *wavy, asym_atom_list);

			const std::vector<i3> hkl_v(hkl.begin(), hkl.end());
			const size_t n_refl = hkl_v.size();

			//atoms the model cannot predict were erased before the prediction; append_spherical_rows() emits them per block
			salted_part_prep spherical;
			std::vector<HE_Spherical_Atom> spheres;
			if (opt.needs_Thakkar_fill)
			{
				std::vector<WFN> tempy;
				const int fill_nr = build_fill_wavefunctions(opt, nr, tempy);
				{
					const spherical_fill_scope fill(opt, hkl);
					calculate_scattering_factors<itsc_block, std::vector<WFN>&>(
						opt, tempy, file, labels, fill_nr, kpts, &spherical);
				}
				if (!spherical.asym_atom_list.empty())
				{
					err_checkf(spherical.k_of_reflection.size() == n_refl,
						"Spherical remainder covers " + std::to_string(spherical.k_of_reflection.size()) +
						" reflections, the table covers " + std::to_string(n_refl), file);
					append_scatterer_ids(stream_ids, opt, spherical.labels,
						tempy[fill_nr], spherical.asym_atom_list);
					spheres = make_he_evaluators(spherical, tempy[fill_nr], opt, file);
					file << "Spherical remainder: " << spherical.asym_atom_list.size()
						 << " atom(s) the model cannot predict" << endl;
				}
				opt.needs_Thakkar_fill = false;   // emitted here, not appended later
			}

			stream_blocks(opt, file, "experimental.tscb", stream_ids, hkl_v, n_refl,
				[&](const size_t lo, const size_t hi, ProgressBar& progress)
				{
					cvec2 chunk = ml->scattering_factors(slice_k_points(k_pt, lo, hi), asym_atom_list, &progress);
					//the core before Mott-Bethe, which wants all electrons
					if (calculator.wavy.get_has_ECPs())
					{
						vec stl(hi - lo);
						for (size_t s = lo; s < hi; s++)
							stl[s - lo] = unit_cell.get_stl_of_hkl(hkl_v[s]);
						add_ECP_contribution(asym_atom_list, calculator.wavy, chunk, stl, file,
							calculator.wavy.get_ECP_mode(), false);
					}
					//Mott-Bethe never looks outside one reflection, so a block is as valid a unit as a table
					if (opt.electron_diffraction)
						convert_to_ED(asym_atom_list, *wavy, chunk, unit_cell,
							std::vector<i3>(hkl_v.begin() + lo, hkl_v.begin() + hi));
					if (!spheres.empty())
						append_spherical_rows(chunk, spherical, spheres,
							opt.electron_diffraction, lo, hi);
					return chunk;
				});
		}
		else
			sf = ml->scattering_factors(k_pt, asym_atom_list);
		file << setw(13 * 4) << "... done!\n"
			<< flush;
		time_points.push_back(get_time());
		time_descriptions.push_back("Fourier transform");

		//streamed blocks got theirs as they were made; the mode is the model's (-salted_ecp) or -ECP's
		if (calculator.wavy.get_has_ECPs() && !sf.empty())
		{
			add_ECP_contribution(
				asym_atom_list,
				calculator.wavy,
				sf,
				unit_cell,
				hkl,
				file,
				calculator.wavy.get_ECP_mode(),
				opt.debug);
		}
	}
	else if constexpr (std::is_same_v<calculator_type, std::vector<WFN> &>) {
		if (opt.partition_type == PartitionType::Hirshfeld ||
			opt.partition_type == PartitionType::Becke ||
			opt.partition_type == PartitionType::TFVC ||
			opt.partition_type == PartitionType::MBIS ||
			opt.partition_type == PartitionType::EMBIS)
		{
			vec2 d1, d2, d3, dens;
			sf_shell_layout layout;
			const int points = make_atomic_grids_wrapper(
				*wavy,
				needs_grid,
				asym_atom_list,
				unit_cell,
				labels,
				time_points,
				time_descriptions,
				d1, d2, d3, dens,
				opt, file, &layout);

			_time_point end1;
			calc_SF(points,
				k_pt,
				d1, d2, d3, dens,
				sf,
				file,
				time_points.front(),
				end1,
				opt.debug,
				opt.no_date,
				false,
				opt.use_gpu,
				opt.gpu_fp64,
				opt.gpu_fp32,
				&layout);

			time_points.push_back(get_time());
			time_descriptions.push_back("Fourier transform");
		}
		else if (opt.partition_type == PartitionType::RI)
		{
			file << "\nGenerating densities... " << endl;
			WFN wavy_aux = generate_aux_wfn(*wavy, opt.aux_basis, false);

			//TODO: only compute coefs for atoms that are actually in the symmetric unit!
			DensityFitting::CONFIG config = DensityFitting::config_from_options(opt);
			config.asym_atm_list = asym_atom_list;
			//config.restrain_type = DensityFitting::RESTRAINT_TYPE::SIMPLE_AND_TIK;
			//config.charge_scheme = DensityFitting::CHARGE_SCHEME::HIRSHFELD;
			//if (wavy->get_origin() == e_origin::ptb)
			//    config.restraint_strength = 1.0e-4;

			const Gaussian_Molecule ml(wavy_aux, DensityFitting::density_fit(*wavy, wavy_aux, config));
			file << setw(12 * 4 + 2) << "... done!\n"
				<< flush;
			time_points.push_back(get_time());
			time_descriptions.push_back("RI-Fit");

			vec atom_elecs = ml.populations();
			file << "Table of Charges in electrons\n"
				<< "       Atom  Charge_RI" << endl;

			for (int i = 0; i < asym_atom_list.size(); i++)
			{
				int a = asym_atom_list[i];
				file << setw(10) << labels[i]
					<< fixed << setw(10) << setprecision(3) << wavy_aux.get_atom_charge(a) - atom_elecs[a];
				if (opt.debug)
					file << " " << setw(4) << wavy_aux.get_atom_charge(a) << " " << fixed << setw(10) << setprecision(3) << atom_elecs[a];
				file << "\n";
			}

			auto el_sum = reduce(atom_elecs.begin(), atom_elecs.end(), 0.0);
			file << setprecision(4) << "Total number of analytical Electrons: " << el_sum << endl;
			time_points.push_back(get_time());
			time_descriptions.push_back("Calculation of Charges");
			if (tshb_only) {
				ScattererLabels keys;
				append_scatterer_ids(keys, opt, labels, *wavy, asym_atom_list);
				opt.tshb->append(ml.tshb(labels, keys, asym_atom_list, *wavy, wavy->get_ECP_mode(), opt.electron_diffraction), file);
				file << labels.size() << " atom(s) go into the .tshb, no table" << endl;
			}
			else if (stream_tsc)
			{
				//as the SALTED case above, but the atoms come from the auxiliary wavefunction, so no spherical remainder
				ScattererLabels stream_ids;
				append_scatterer_ids(stream_ids, opt, labels, *wavy, asym_atom_list);

				const std::vector<i3> hkl_v(hkl.begin(), hkl.end());
				stream_blocks(opt, file, "experimental.tscb",
					stream_ids, hkl_v, hkl_v.size(),
					[&](const size_t lo, const size_t hi, ProgressBar& progress)
					{
						cvec2 chunk = ml.scattering_factors(slice_k_points(k_pt, lo, hi), asym_atom_list, &progress);
						if (wavy->get_has_ECPs())
						{
							vec stl(hi - lo);
							for (size_t s = lo; s < hi; s++)
								stl[s - lo] = unit_cell.get_stl_of_hkl(hkl_v[s]);
							add_ECP_contribution(asym_atom_list, *wavy, chunk, stl, file, wavy->get_ECP_mode(), false);
						}
						if (opt.electron_diffraction)
							convert_to_ED(asym_atom_list, *wavy, chunk, unit_cell,
								std::vector<i3>(hkl_v.begin() + lo, hkl_v.begin() + hi));
						return chunk;
					});
			}
			else
				sf = ml.scattering_factors(k_pt, asym_atom_list);
			file << setw(12 * 4 + 2) << "... done!" << endl;
			time_points.push_back(get_time());
			time_descriptions.push_back("Fourier transform");
		}
		else {
			std::cout << "Unknown Partition type, stopping here!" << std::endl;
		}
		//a streamed RI table is empty here, its blocks got the core as they were made; Hirshfeld fills sf even when streaming
		if (wavy->get_has_ECPs() && !sf.empty())
		{
			add_ECP_contribution(
				asym_atom_list,
				*wavy,
				sf,
				unit_cell,
				hkl,
				file,
				wavy->get_ECP_mode(),
				opt.debug);
		}
	}

	//not when streamed: sf is empty there and each block was converted as it was produced (calc_SF fills sf whole
	//even under stream_tsc, so a Hirshfeld table still needs it)
	if (opt.electron_diffraction && !opt.iam_switch && !sf.empty())
	{
		convert_to_ED(asym_atom_list,
			*wavy,
			sf,
			unit_cell,
			hkl);
	}

	tsc_block_type blocky;
	if (opt.label_tsc_output)
	{
		blocky = tsc_block_type(sf, labels, hkl);
	}
	else
	{
		std::vector<atomID> IDs(asym_atom_list.size());
		for (int atm_idx = 0; atm_idx < asym_atom_list.size(); atm_idx++) {
			IDs[atm_idx] = wavy->get_id_for_atom(asym_atom_list[atm_idx]);
		}
		blocky = tsc_block_type(sf, IDs, hkl);
	}


	if (opt.needs_Thakkar_fill)
	{
		file << "Performing the remaining calculation of spherical atoms..." << std::endl;
		opt.needs_Thakkar_fill = false;
		const auto t_fill = std::chrono::steady_clock::now();
		std::vector<WFN> tempy;
		const int fill_nr = build_fill_wavefunctions(opt, nr, tempy);
		const auto t_wfn = std::chrono::steady_clock::now();
		//no_date is not part of spherical_fill_scope: only this caller suppresses the banner
		const bool no_date_was = opt.no_date;
		opt.no_date = true;
		//the atoms of earlier parts too: the table and the .tshb keep the first copy of a key, so a fill of them was discarded
		svec fill_known(known_atoms);
		fill_known.insert(fill_known.end(), labels.begin(), labels.end());
		if (tshb_only)
		{
			//the atoms alone: their records need no reflections
			salted_part_prep spherical;
			{
				const spherical_fill_scope fill(opt, hkl);
				calculate_scattering_factors<itsc_block, std::vector<WFN> &>(opt, tempy, file, fill_known, fill_nr, kpts, &spherical);
			}
			opt.no_date = no_date_was;
			const auto t_atoms = std::chrono::steady_clock::now();
			fill_tshb(spherical, tempy[fill_nr], opt, file);
			if (ProgressBar::report_counts)
				std::cout << "[stages] fill: wavefunctions " << std::chrono::duration<double>(t_wfn - t_fill).count()
						  << " s, atoms " << std::chrono::duration<double>(t_atoms - t_wfn).count()
						  << " s, records " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_atoms).count() << " s" << endl;
		}
		else
		{
			tsc_block<int, cdouble> blocky_thakkar;
			{
				const spherical_fill_scope fill(opt, hkl);
				blocky_thakkar = calculate_scattering_factors<itsc_block, std::vector<WFN> &>(opt, tempy, file, fill_known, fill_nr);
			}
			opt.no_date = no_date_was;
			if (!blocky_thakkar.is_empty())
				blocky.append(std::move(blocky_thakkar), file);
		}
		time_points.push_back(get_time());
		time_descriptions.push_back("Spherical Atoms");
	}

	time_points.push_back(get_time());
	time_descriptions.push_back("tsc calculation");
	if (!opt.no_date)
	{
		write_timing_to_file(file,
			time_points,
			time_descriptions);
	}
	return blocky;
}
//streams a combined (-mtc) table one block of reflections at a time, inverting the sequential loops:
//prepare every disorder part once, then per block compute every part and hand the assembled block to the writer
//the inversion works because the expensive parts (CIF read, SALTED prediction) do not depend on reflections,
//which is what salted_part_prep holds
//row order must equal the merged path's exactly, or values are attributed to the wrong atom without any error
//it is a plain concatenation because each part is prepared with the identifiers of the parts before it (the
//growing `known` list) and read_atoms_from_CIF skips atoms already in it
//that list must hold atomID hex strings, as result.get_scatterers_string() gives: labels are not unique across
//parts, so a later part skips atoms it should keep and the table silently loses rows - check any change against
//a multi-part structure byte for byte, not by size
//species the model does not know are erased from every part and return as extra Thakkar rows of each block
bool stream_mtc_salted(options& opt, std::vector<WFN>& wavy, std::ostream& file, vec2* known_kpts)
{
	const size_t n_parts = opt.combined_tsc_calc_files.size();
	//electron diffraction is not excluded: the prep carries the stl, so no unit cell is needed down here
	//the stream writer is binary only, a text table has to go through the in-memory path
	//a .tshb run makes no table, so the sequential loop (equal row for row) collects the parts instead
	if (opt.tsc_block_size == 0 || !opt.SALTED || n_parts < 2 || opt.iam_switch || !opt.binary_tsc || opt.old_tsc || opt.tshb)
		return false;

	std::vector<std::shared_ptr<SALTEDPredictor>> preds;
	std::vector<salted_part_prep> preps;
	//each predictor replaces opt.spherical_fill_charges with its part's EEQ charges; the fill below needs every part's
	std::vector<std::vector<std::array<double, 4>>> charges(n_parts);
	svec known;
	auto mark = std::chrono::steady_clock::now();
	auto lap = [&mark](const std::string& what) {
		const auto now = std::chrono::steady_clock::now();
		if (ProgressBar::report_counts)
			std::cout << "[stages] " << what << " " << std::chrono::duration<double>(now - mark).count() << " s" << std::endl;
		mark = now;
	};
	std::shared_ptr<SALTED_model_store> model_store = std::make_shared<SALTED_model_store>();
	for (size_t i = 0; i < n_parts; i++)
	{
		auto pred = std::make_shared<SALTEDPredictor>(wavy[i], opt);
		pred->share_model(model_store);
		lap("part " + std::to_string(i + 1) + ": predictor");
		charges[i] = opt.spherical_fill_charges;
		if (!pred->basis_set_loaded())
			load_basis_into_WFN(pred->wavy, BasisSetLibrary::get_basis_set(pred->get_dfbasis_name()));
		salted_part_prep prep;
		calculate_scattering_factors<itsc_block, SALTEDPredictor&>(
			opt, *pred, file, known, static_cast<int>(i), known_kpts, &prep);
		//feed the next part what the sequential path feeds it: get_scatterers_string(), the atomID hex string
		for (size_t a = 0; a < prep.asym_atom_list.size(); a++)
		{
			if (opt.label_tsc_output)
				known.push_back(prep.labels[a]);
			else
				known.push_back(pred->wavy.get_id_for_atom(prep.asym_atom_list[a]).to_hex_string());
		}
		preds.push_back(pred);
		preps.push_back(std::move(prep));
		lap("part " + std::to_string(i + 1) + ": prediction and setup");
	}
	model_store.reset();

	//atoms the model cannot predict (unknown species, or nothing inside the descriptor cutoff) become extra rows here
	//every part must be walked: a part file only yields atoms of the parts it covers, so part 0 alone drops the rest
	//feed `known` as we go so no atom is produced twice
	std::vector<salted_part_prep> spherical(n_parts);
	std::vector<char> have_spherical(n_parts, 0);
	//a table must hold either all atomIDs or all labels; salted_part_prep::labels holds hex strings while the
	//predicted parts contribute atomID objects, so the spherical rows are converted here rather than passed through
	std::vector<ScattererLabels> spherical_ids(n_parts);
	//built in the block below, where the fill wavefunctions the charges are matched against exist
	std::vector<std::vector<HE_Spherical_Atom>> spheres(n_parts);
	//every part reads the same reflections, but a part the model predicts nothing of returns none
	const auto ref = std::find_if(preps.begin(), preps.end(), [](const salted_part_prep& p) { return !p.hkl_v.empty(); });
	err_checkf(ref != preps.end(), "The SALTED model predicts no atom of any part", file);
	const std::vector<i3>& hkl_v = ref->hkl_v;
	if (opt.needs_Thakkar_fill)
	{
		hkl_list fill_reflections;
		for (const auto& h : hkl_v)
			fill_reflections.emplace(h);
		const spherical_fill_scope fill(opt, fill_reflections);

		//the whole file list, not one file: the index also selects opt.groups[nr], see the sequential path above
		std::vector<WFN> tempy;
		for (const auto& part_file : opt.combined_tsc_calc_files)
			tempy.emplace_back(part_file);

		size_t n_filled = 0;
		for (size_t i = 0; i < n_parts; i++)
		{
			calculate_scattering_factors<itsc_block, std::vector<WFN>&>(
				opt, tempy, file, known, static_cast<int>(i), known_kpts, &spherical[i]);
			have_spherical[i] = spherical[i].asym_atom_list.empty() ? 0 : 1;
			if (!have_spherical[i]) continue;
			n_filled += spherical[i].asym_atom_list.size();
			err_checkf(spherical[i].k_of_reflection.size() == hkl_v.size(),
				"Spherical remainder of part " + std::to_string(i + 1) + " covers " +
				std::to_string(spherical[i].k_of_reflection.size()) +
				" reflections, the parts cover " + std::to_string(hkl_v.size()), file);
			// keep feeding `known`, so a later part cannot claim these atoms again
			for (size_t a = 0; a < spherical[i].labels.size(); a++)
				known.push_back(spherical[i].labels[a]);
			append_scatterer_ids(spherical_ids[i], opt, spherical[i].labels,
				tempy[i], spherical[i].asym_atom_list);
			opt.spherical_fill_charges = charges[i];
			spheres[i] = make_he_evaluators(spherical[i], tempy[i], opt, file);
		}
		file << "Spherical remainder: " << n_filled
			 << " atom(s) the model cannot predict" << std::endl;
		lap("spherical fill");
	}

	//each spherical remainder straight after its own predicted atoms, as the sequential path appends it per part;
	//collecting them at the end gives the same scatterers in a different order and breaks a byte-for-byte comparison
	ScattererLabels ids;
	for (size_t p = 0; p < preps.size(); p++)
	{
		append_scatterer_ids(ids, opt, preps[p].labels, preds[p]->wavy, preps[p].asym_atom_list);
		for (const auto& sid : spherical_ids[p])
			ids.emplace_back(sid);
	}

	const size_t n_refl = hkl_v.size();
	file << "Combined tsc: " << ids.size() << " scatterers from "
		<< n_parts << " parts" << std::endl;
	//the bar counts reflections * parts, since every part is evaluated for every block
	stream_blocks(opt, file, "experimental.tscb", ids, hkl_v,
		n_refl * preps.size(),
		[&](const size_t lo, const size_t hi, ProgressBar& progress)
		{
			cvec2 combined;
			combined.reserve(ids.size());
			for (size_t p = 0; p < preps.size(); p++)
			{
				if (!preps[p].mol)  // the model predicts nothing in this part, only its spherical rows
				{
					if (have_spherical[p])
						append_spherical_rows(combined, spherical[p], spheres[p], opt.electron_diffraction, lo, hi);
					continue;
				}
				cvec2 chunk = preps[p].mol->scattering_factors(slice_k_points(preps[p].k_pt, lo, hi), preps[p].asym_atom_list, &progress);
				const vec stl(preps[p].stl_of_reflection.begin() + lo, preps[p].stl_of_reflection.begin() + hi);
				//the core before Mott-Bethe, as in the single-table stream
				if (preds[p]->wavy.get_has_ECPs())
					add_ECP_contribution(preps[p].asym_atom_list, preds[p]->wavy, chunk, stl, file,
						preds[p]->wavy.get_ECP_mode(), false);
				if (opt.electron_diffraction)
					convert_to_ED(preps[p].asym_atom_list, preds[p]->wavy, chunk, stl);
				for (auto& row : chunk)
					combined.push_back(std::move(row));
				// this part's spherical remainder, where the id list put it
				if (have_spherical[p])
					append_spherical_rows(combined, spherical[p], spheres[p],
						opt.electron_diffraction, lo, hi);
			}
			return combined;
		});
	lap("factors and table, streamed");
	return true;
}

template itsc_block calculate_scattering_factors(options& opt,
	std::vector<WFN>& calculator,
	std::ostream& file,
	svec& known_atoms,
	const int& nr,
	vec2* kpts,
	salted_part_prep* prep_out);
template itsc_block calculate_scattering_factors(options& opt,
	SALTEDPredictor& calculator,
	std::ostream& file,
	svec& known_atoms,
	const int& nr,
	vec2* kpts,
	salted_part_prep* prep_out);


/**
 * Calculates the diffuse (that is non integer hkl) scattering factors based on the given options.
 *
 * @param opt The options for calculating the diffuse scattering factors.
 * @param log_file The output stream for logging the calculation process.
 */
void calc_sfac_diffuse(const options& opt, std::ostream& log_file)
{
	using namespace std;
	std::vector<WFN> wavy;
	wavy.emplace_back(e_origin::wfn);
	wavy[0].read_known_wavefunction_format(opt.wfn, std::cout, opt.debug);
	// set number of threads
#ifdef _OPENMP
	if (opt.threads > 0)
		omp_set_num_threads(opt.threads);
#endif
	err_checkf(wavy[0].get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
	err_checkf(exists(opt.cif), "CIF does not exists!", std::cout);
	// err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);

	// time_point start = get_time();
	// time_point end_becke, end_prototypes, end_spherical, end_prune, end_aspherical;
	vector<_time_point> time_points;
	svec time_descriptions;
	time_points.push_back(get_time());

	cell unit_cell(opt.cif, std::cout, opt.debug);
	ifstream cif_input(opt.cif.c_str(), std::ios::in);
	ivec atom_type_list;
	ivec asym_atom_to_type_list;
	ivec asym_atom_list;
	bvec needs_grid(wavy[0].get_ncen(), false);
	svec known_atoms;

	auto labels = read_atoms_from_CIF(cif_input,
		opt.groups[0],
		unit_cell,
		wavy[0],
		known_atoms,
		atom_type_list,
		asym_atom_to_type_list,
		asym_atom_list,
		needs_grid,
		std::cout,
		opt.debug);

	cif_input.close();
	vec2 d1_, d2_, d3_, dens;

	make_atomic_grids_wrapper(
		wavy[0],
		needs_grid,
		asym_atom_list,
		unit_cell,
		labels,
		time_points,
		time_descriptions,
		d1_, d2_, d3_, dens, opt);

	hkl_list_d hkl;
	generate_fractional_hkl(opt.dmin, hkl, opt.twin_law, unit_cell, log_file, opt.sfac_diffuse, opt.debug);

	const long long int size = static_cast<long long int>(hkl.size());
	//k_pt, sf and the tsc_block's copies of sf and the indices are alive together when the block is built
	const double need = double(size) * (48.0 + 32.0 * dens.size());
	const size_t avail = available_memory_bytes();
	err_checkf(avail == 0 || need < avail, "-sfac_diffuse: " + std::to_string(size) + " fractional reflections x " +
		std::to_string(dens.size()) + " atoms need ~" + std::to_string(int(need / 1073741824.0) + 1) + " GB, " +
		std::to_string(avail >> 30) + " GB are free; use a larger step or dmin", std::cout);
	vec2 k_pt;
	k_pt.resize(3);
#pragma omp parallel for
	for (int i = 0; i < 3; i++)
		k_pt[i].resize(size, 0.0);

	if (opt.debug)
	{
		log_file << "K_point_vector is here! size: " << k_pt[0].size() << endl;
	}
	int i_ = 0;
	for (const d3& hkl_ : hkl)
	{
		for (int x = 0; x < 3; x++)
		{
			for (int j = 0; j < 3; j++)
			{
				k_pt[x][i_] += unit_cell.get_rcm(x, j) * hkl_[j];
			}
		}
		i_++;
	}

	// below is a strip of Calc_SF without the file IO or progress bar
	cvec2 sf;

	const long long int imax = static_cast<long long int>(dens.size());
	const long long int smax = static_cast<long long int>(k_pt[0].size());
	long long int pmax = static_cast<long long int>(dens[0].size());
	//const long long int step = max(static_cast<long long int>(floor(smax / 20)), 1LL);
	std::cout << "Done with making k_pt " << smax << " " << imax << " " << pmax << endl;
	sf.resize(imax);
#pragma omp parallel for
	for (int i = 0; i < imax; i++)
		sf[i].resize(smax);
	double* dens_local, * d1_local, * d2_local, * d3_local;
	complex<double>* sf_local;
	const double* k1_local = k_pt[0].data();
	const double* k2_local = k_pt[1].data();
	const double* k3_local = k_pt[2].data();
	double work, rho;
	ProgressBar* progress = new ProgressBar(imax * smax, 60, "=", " ", "Calculating Scattering Factors");
	for (long long int i = 0; i < imax; i++)
	{
		pmax = static_cast<long long int>(dens[i].size());
		dens_local = dens[i].data();
		d1_local = d1_[i].data();
		d2_local = d2_[i].data();
		d3_local = d3_[i].data();
		sf_local = sf[i].data();
#pragma omp parallel for private(work, rho)
		for (long long int s = 0; s < smax; s++)
		{
			double re = 0.0, im = 0.0, si, c;
			const double& _k1_local = k1_local[s];
			const double& _k2_local = k2_local[s];
			const double& _k3_local = k3_local[s];
			for (long long int p = pmax - 1; p >= 0; p--)
			{
				rho = dens_local[p];
				work = _k1_local * d1_local[p] + _k2_local * d2_local[p] + _k3_local * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
				if (rho < 0)
				{
					rho = -rho;
					work += M_PI;
				}
#endif
#endif
				c = cos(work);
				si = sin(work);
				re += rho * c;
				im += rho * si;
			}
			sf_local[s].real(re);
			sf_local[s].imag(im);
			progress->update();
		}
	}
	delete (progress);
	tsc_block<double, cdouble> result(sf, labels, hkl);
	result.write_tsc_file_non_integer(opt.cif);
}
