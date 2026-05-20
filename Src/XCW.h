#pragma once
#include "convenience.h"
#include "GridManager.h"
#include "integration_params.h"
#include "scattering_factors.h"
#include "cell.h"
#include <occ/qm/hf.h>

class XCW {
public:
	XCW(const options& opt_in)
	{
		opt = &opt_in;
		wave.read_known_wavefunction_format(opt_in.wfn, std::cout, opt_in.debug);
		wave.set_method(opt_in.method);
		wave.set_charge(opt_in.charge);
		wave.set_multi(opt_in.mult);
		params = Int_Params(wave);
		std::filesystem::path hkl_filename = opt_in.hkl;
		std::filesystem::path cif = opt_in.cif;
		std::ifstream cif_input(cif.c_str(), std::ios::in);
		unit_cell = cell(cif, std::cout, opt_in.debug, opt_in.do_XCW);
		obs.resize(2);
		hkl_enlarged = read_hkl_full(hkl_filename, hkl, opt_in.twin_law, unit_cell, std::cout, obs, opt_in.debug);
		nmo = wave.get_nmo();
		svec empty({});
		ivec atom_type_list;
		ivec asym_atom_to_type_list;
		bvec constant_atoms;
		needs_grid = bvec(wave.get_ncen(), false);
		std::ofstream log3("log3.txt", std::ios::out);
		labels = read_atoms_from_CIF(cif_input, opt_in.groups[0], unit_cell, wave, empty, atom_type_list, asym_atom_to_type_list, asym_atom_list, needs_grid, std::cout, constant_atoms, opt_in.SALTED, opt_in.debug);
		U_iso = read_U_iso_from_CIF(cif, wave, unit_cell, log3, opt_in.debug);
		asym_atoms = unit_cell.get_asym_atoms(wave, labels, atom_type_list, asym_atom_to_type_list, asym_atom_list);
		ncen = asym_atoms.size();
		nr = hkl_enlarged.size();
		nr_small = hkl.size();
		make_k_pts(nr != 0 && hkl.size() == 0, opt_in.save_k_pts, unit_cell, hkl_enlarged, k_pt, std::cout, opt_in.debug);
		//Read U_iso, U_ij, C_ijk and Dijkl from CIF file
		read_fracs_ADPs_from_CIF(cif, wave, unit_cell, log3, opt_in.debug);
		unit_cell.eval_symm(asym_atoms);
		n_params = nmo * nmo / 2;
		DW_fact.resize(ncen, vec(nr, 1.0));

	};
	void eval_DW();
	void eval_phase();
	void eval_translation_phase();
	cvec3 eval_I();
	cvec eval_anom_disp();
	void calc_F_calc(const cvec3& I, const cvec& corr, cvec& F_calc);
	void calc_perturb(vec2& perturb, const cvec& F_calc, const double& scale, cvec3& I);
	void do_SCF();
	void calc_F_calc_fast(const cvec& corr);
	void run_XCW_fitting();

private:
	struct ao_data {
		std::vector<primitive> prims;
		d3 pos;
		int m;
	};
	struct anom_atom {
		std::string identifier;
		cdouble dispersion;
	};
	void set_phase(cvec2& phase_in) { phase_fact = phase_in; }
	void set_DW(vec2& DW_in) { DW_fact = DW_in; }
	void set_translation_phase(cvec2& phase_in) { translation_phase = phase_in; }
	void U_frac2U_rec();
	void U_rec2U_cart();
	ivec generate_asym_lookup(const int r);
	cvec2 eval_integrals(const ao_data& mu, const ao_data& nu, GridManager& grid_manager, const vec2& k_pt);
	cvec2 calculateXCWintegral(GridManager& grid_manager, const ao_data& mu_data, const ao_data& nu_data, const ivec& asym_atom_list, vec2& k_pt, bool& equal, vec2& mu_vals);
	void parse_anom_atoms(std::vector<anom_atom>& anom_atoms);
	void SCF_iteration(occ::qm::SCF<occ::qm::HartreeFock>& scf, double& ehf_last, occ::Mat& D_last, occ::Mat& D_diff, bool& incremental, occ::Mat& F_diis);
	void SCF_convergence_check(occ::qm::SCF<occ::qm::HartreeFock>& scf, bool& converged, const double& energy_conv, const double& commutator_conv, const double& ehf_last);
	const options* opt;
	WFN wave;
	Int_Params params;
	cell unit_cell;
	hkl_list hkl;
	hkl_list hkl_enlarged;
	vec2 obs;
	vec2 k_pt;
	cvec2 phase_fact;
	cvec2 translation_phase;
	bvec needs_grid;
	vec U_iso;
	ivec asym_atom_list;
	std::vector<cell::asym_atom> asym_atoms;
	int nmo;
	int ncen;
	int nr;
	int nr_small;
	int n_params;
	svec labels;
	vec2 DW_fact;
};