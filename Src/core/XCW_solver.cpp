#include "pch.h"
#include "XCW_solver.h"
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
#include "itensor_gpu.h"
#endif
#include "scattering_factors.h"
#include "basis_set.h"
#include "bondwise_analysis.h"
#include <limits>

XCW_solver::XCW_solver(structure_factors& sf_in) : scf_solver(*this, sf_in) {
	sf = &sf_in;
	opt = sf->opt;
	std::cout << "XCW orbital basis set: " << opt->xcw_settings.basis_set_name << std::endl;
	std::cout << "XCW: I/sigma(I) >= " << opt->xcw_settings.i_sigma_cutoff << " (F/sigma(F) >= " << 2 * opt->xcw_settings.i_sigma_cutoff << "): " << sf->nr_fit << " of " << sf->model_data.nr << " reflections in the fit; R1 and Criterion are over these, R1(all) and Crit(all) over all" << std::endl;
}

//z_h = (|F_obs,h| - |F_calc,h|) / sigma_h for the converged F_calc/F_scale over
//strong reflections (|F_obs|/sigma >= opt->xcw_strong_cutoff), tested against
//N(0,1) after a global shape/scale rescale. Full reflection set only; the
//free/working-set cross-validation variant is not implemented.
//Background: tests/P1_test/XCW_plan.md, Src/core/xcw_halting.h.
void XCW_solver::evaluate_gaussian_halting(const double lambda) {
	sf->ensure_hkl_ordered();

	const int n_total = sf->model_data.nr;
	vec z_raw;
	vec resolution;
	vec abs_F;
	z_raw.reserve(n_total);
	resolution.reserve(n_total);
	abs_F.reserve(n_total);

	for (int i = 0; i < n_total; i++) {
		if (sf->scatter_data.sigma_obs[i] <= 0.0) {
			continue;
		}
		const double f_over_sigma = sf->scatter_data.abs_F_obs[i] / sf->scatter_data.sigma_obs[i];
		if (f_over_sigma < opt->xcw_settings.xcw_strong_cutoff) {
			continue;
		}
		const double scaled_F_calc = sf->scatter_data.scale * std::abs(sf->scatter_data.F_calc[i]);
		const double diff = scaled_F_calc - sf->scatter_data.abs_F_obs[i];
		z_raw.push_back(diff / sf->scatter_data.sigma_obs[i]);
		resolution.push_back(sf->unit_cell.get_stl_of_hkl(sf->hkl_ordered_[i]));
		abs_F.push_back(sf->scatter_data.abs_F_obs[i]);
	}

	GaussianHaltEntry entry;
	entry.lambda = lambda;
	entry.n_total = n_total;
	entry.n_used = static_cast<int>(z_raw.size());

	if (entry.n_used < 8) {
		sf->XCW_log << "Gaussian halting criterion: only " << entry.n_used
			<< " strong reflections (|F|/sigma >= " << opt->xcw_settings.xcw_strong_cutoff
			<< ") at lambda=" << lambda << ", skipping (need >= 8)." << std::endl;
		gaussian_halt_history_.push_back(entry);
		return;
	}

	//Decouple shape from scale: one global factor rescaling z so <z^2> ~ 1, not the
	//resolution-uniform weighting of XCW_plan.md 4.1
	double mean_z2 = 0.0;
	for (const double v : z_raw) {
		mean_z2 += v * v;
	}
	mean_z2 /= entry.n_used;
	entry.sigma_scale = (mean_z2 > 0.0) ? std::sqrt(mean_z2) : 1.0;

	vec z(entry.n_used);
	for (int i = 0; i < entry.n_used; i++) {
		z[i] = z_raw[i] / entry.sigma_scale;
	}

	entry.A2 = anderson_darling_statistic(z);
	entry.ad_reject_5pct = entry.A2 > ANDERSON_DARLING_CRITICAL_5PCT;

	const ProbabilityPlotFit pp = normal_probability_plot_fit(z);
	entry.pp_slope = pp.slope;
	entry.pp_intercept = pp.intercept;

	entry.skewness = sample_skewness(z);
	entry.excess_kurtosis = sample_excess_kurtosis(z);
	entry.jarque_bera = jarque_bera_statistic(z);

	const int n_bins = std::clamp(entry.n_used / 20, 2, 10);
	const BinnedTrend res_trend = binned_z_squared_trend(z, resolution, n_bins);
	entry.resolution_trend_slope = res_trend.slope;
	entry.resolution_trend_r = res_trend.spearman_r;
	entry.resolution_trend_flagged = res_trend.flagged;

	const BinnedTrend int_trend = binned_z_squared_trend(z, abs_F, n_bins);
	entry.intensity_trend_slope = int_trend.slope;
	entry.intensity_trend_r = int_trend.spearman_r;
	entry.intensity_trend_flagged = int_trend.flagged;

	sf->XCW_log << "Gaussian halting criterion at lambda=" << std::fixed << std::setprecision(5) << lambda << ":\n"
		<< "  n_used=" << entry.n_used << "/" << entry.n_total << " (|F|/sigma >= " << opt->xcw_settings.xcw_strong_cutoff
		<< "), sigma_scale=" << entry.sigma_scale << "\n"
		<< "  A^2=" << entry.A2 << (entry.ad_reject_5pct ? " (rejects N(0,1) at 5%)" : " (consistent with N(0,1) at 5%)") << "\n"
		<< "  probability-plot slope=" << entry.pp_slope << " intercept=" << entry.pp_intercept << "\n"
		<< "  skewness=" << entry.skewness << " excess_kurtosis=" << entry.excess_kurtosis
		<< " Jarque-Bera=" << entry.jarque_bera << "\n"
		<< "  resolution-binned <z^2> trend: slope=" << entry.resolution_trend_slope
		<< " spearman_r=" << entry.resolution_trend_r << (entry.resolution_trend_flagged ? " [FLAGGED]" : "") << "\n"
		<< "  |F|-binned <z^2> trend: slope=" << entry.intensity_trend_slope
		<< " spearman_r=" << entry.intensity_trend_r << (entry.intensity_trend_flagged ? " [FLAGGED]" : "") << std::endl;

	gaussian_halt_history_.push_back(entry);
}

//Full per-lambda table to XCW_log, then the final recommendation
void XCW_solver::report_gaussian_halting_summary() {
	if (gaussian_halt_history_.empty()) {
		return;
	}

	sf->XCW_log << "\n____________________________________________________________________________\n"
		<< "Gaussian halting criterion summary (tests/P1_test/XCW_plan.md)\n"
		<< " Lambda\t\tA^2\treject5%\tpp_slope\tpp_intercept\tskew\tkurt\tres_trend_r\tint_trend_r\tn_used\n";
	for (const GaussianHaltEntry& e : gaussian_halt_history_) {
		sf->XCW_log << "\t" << std::fixed << std::setprecision(5) << e.lambda
			<< "\t" << std::setprecision(4) << e.A2
			<< "\t" << (e.ad_reject_5pct ? "yes" : "no")
			<< "\t\t" << e.pp_slope << "\t\t" << e.pp_intercept
			<< "\t" << e.skewness << "\t" << e.excess_kurtosis
			<< "\t" << e.resolution_trend_r << "\t\t" << e.intensity_trend_r
			<< "\t" << e.n_used << "\n";
	}

	report_halting_progress_estimate(true);
}

//lambda* is the argmin of A^2 over steps with enough strong reflections. A flagged
//binned-trend test means spatially correlated residuals that a marginal normality
//test would miss, so it is warned about. See XCW.h for the full behaviour.
void XCW_solver::report_halting_progress_estimate(bool is_final) {
	const GaussianHaltEntry* best = nullptr;
	double max_valid_lambda = 0.0;
	vec fit_lambda, fit_A2;
	for (const GaussianHaltEntry& e : gaussian_halt_history_) {
		if (e.n_used < 8) {
			continue;
		}
		if (!best || e.A2 < best->A2) {
			best = &e;
		}
		max_valid_lambda = std::max(max_valid_lambda, e.lambda);
		fit_lambda.push_back(e.lambda);
		fit_A2.push_back(e.A2);
	}
	if (!best) {
		return;
	}

	//A minimum at the last evaluated lambda is a scan-boundary artifact, not a
	//reached lambda*, so it is flagged rather than reported as the answer
	const bool at_boundary = (best->lambda >= max_valid_lambda - 1e-12) && (fit_lambda.size() > 1);
	const bool trend_ok = !best->resolution_trend_flagged && !best->intensity_trend_flagged;

	//AIC picks between candidate forms of the A^2(lambda) trend; quartic only enters
	//once the scan has enough steps, via fit_polynomial's degree + 3 minimum
	std::vector<PolynomialFit> candidates;
	const PolynomialFit fit = choose_best_polynomial_fit(fit_lambda, fit_A2, { 2, 4 }, &candidates);

	std::ostream* streams[2] = { &sf->XCW_log, &std::cout };
	for (std::ostream* s : streams) {
		*s << "____________________________________________________________________________\n";
		if (!is_final) {
			*s << "Gaussian halting criterion: progress update after " << fit_lambda.size() << " lambda steps\n";
		}
		*s << "Recommended halting lambda* = " << std::fixed << std::setprecision(5) << best->lambda
			<< " (A^2=" << std::setprecision(4) << best->A2 << ")";
		if (!trend_ok) {
			*s << " -- WARNING: binned <z^2> trend test flagged at this lambda; "
				<< "residuals may be spatially correlated, inspect before trusting lambda*.";
		}
		*s << "\n";

		for (const PolynomialFit& c : candidates) {
			*s << "  candidate fit: degree=" << c.degree;
			if (!c.valid) {
				*s << " -- not enough points yet (need >= degree+3 evaluated lambda steps)\n";
				continue;
			}
			*s << " RSS=" << std::setprecision(4) << c.rss << " R^2=" << c.r_squared
				<< " AIC=" << c.aic << ((fit.valid && c.degree == fit.degree) ? " [chosen]" : "")
				<< "\n";
		}

		if (at_boundary) {
			*s << "WARNING: A^2 is still falling at the last evaluated lambda (" << std::setprecision(5)
				<< max_valid_lambda << "); this is a scan-boundary value, not a confirmed interior minimum. "
				<< "Extend the scan (larger max_value in -do_XCW stepsize max_value) to find the true optimum.";
			if (fit.valid && fit.has_minimum && fit.vertex_x > max_valid_lambda) {
				*s << " Degree-" << fit.degree << " polynomial extrapolation (best fit by AIC, R^2="
					<< std::setprecision(4) << fit.r_squared << ") of the A^2(lambda) trend so far estimates "
					<< "the minimum near lambda ~= " << std::setprecision(5) << fit.vertex_x
					<< " (extrapolated, not yet observed -- treat as a rough guide to where to extend the scan, not a final answer).";
			}
			else if (!fit.valid) {
				*s << " Not enough points yet to extrapolate an estimate (need >= 5 evaluated lambda steps).";
			}
			else if (!fit.has_minimum) {
				*s << " The trend so far is not curving upward yet within the search window; "
					<< "cannot extrapolate a stopping estimate, extend the scan further.";
			}
			else {
				//fit.vertex_x <= max_valid_lambda: the fitted curve turns before the last step
				//while the observed A^2 is still falling there - a flat, noisy trend, so no
				//estimate is offered
				*s << " The degree-" << fit.degree << " fit places its minimum at lambda ~= "
					<< std::setprecision(5) << fit.vertex_x << ", inside the range already scanned, "
					<< "which contradicts A^2 still falling at the last step; the trend is too flat "
					<< "here to extrapolate from, only further steps can decide it.";
			}
			*s << "\n";
		}
	}
}

void XCW_solver::calc_perturb(occ::Mat& perturb, const occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	sf->ensure_inv_H2_weights();

	//The four (XWR_type, refine_against) combinations differ only in the per-reflection
	//scalar and the prefactor, so one walk over I serves all of them
	const int xwr = static_cast<int>(opt->xcw_settings.XWR_type), ref = static_cast<int>(opt->xcw_settings.refine_against);
	const bool against_F2 = (ref == 2);
	const bool weighted = (xwr == 2);
	const bool valid = (xwr == 1 || xwr == 2) && (ref == 1 || ref == 2);
	if (!valid) sf->XCW_log << "Invalid refinement option" << std::endl;
	const double scale_sq = sf->scatter_data.scale * sf->scatter_data.scale;
	const double prefactor = against_F2
		? 4.0 * scale_sq / (sf->nr_fit - sf->n_params())
		: 2.0 * sf->scatter_data.scale / (sf->nr_fit - sf->n_params());

	cvec pre(sf->model_data.nr);
#pragma omp parallel for
	for (int r = 0; r < sf->model_data.nr; r++) {
		if (!valid || !sf->scatter_data.hkl_mask[r]) continue;
		cdouble precompute;
		//with extinction the model is I = y |Fc|^2, so the residual carries sqrt(y)|Fc| and the
		//carrier d/dD picks up dI/d|Fc|^2 (F^2) or d(sqrt(y)|Fc|)/d|Fc| (F); both are 1 without
		//a model, and the expressions below are then exactly the ones this always used
		if (against_F2) {
			const double I = sf->ext_y(r) * std::pow(std::abs(sf->scatter_data.F_calc[r]), 2);
			precompute = sf->ext_g(r) * std::conj(sf->scatter_data.F_calc[r]) * (scale_sq * I - sf->scatter_data.F_obs2[r]) / (sf->scatter_data.sigma_obs2[r] * sf->scatter_data.sigma_obs2[r]);
		}
		else {
			const double F_calc_abs = std::abs(sf->scatter_data.F_calc[r]);
			precompute = sf->ext_m(r) * std::conj(sf->scatter_data.F_calc[r]) * (sf->scatter_data.scale * sf->ext_sqrt_y(r) * F_calc_abs - sf->scatter_data.abs_F_obs[r]) / (sf->scatter_data.sigma_obs[r] * sf->scatter_data.sigma_obs[r] * F_calc_abs);
		}
		if (weighted) precompute *= sf->inv_H2_[r];
		pre[r] = precompute;
	}
	contract_I(perturb, pre);
	perturb *= prefactor;
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted) {
		perturb.conservativeResize(2 * sf->model_data.nmo, Eigen::NoChange);
		perturb.bottomRows(sf->model_data.nmo) = perturb.topRows(sf->model_data.nmo);
	}
	//closing function
}

void XCW_solver::contract_I(occ::Mat& out, const cvec& pre) {
	out.setZero(sf->model_data.nmo, sf->model_data.nmo);
	const int step = std::max(1, I_tens->i_streamed_ ? I_tens->i_window_ : sf->model_data.nr);
	bool on_device = false;
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (I_tens->i_on_device_) {
		vec dev(I_tens->i_compact_);
		err_checkf(itensor_gpu_cols(pre.data(), dev.data()), "I tensor walk on the device failed", std::cout);
		for (size_t k = 0; k < I_tens->i_compact_; k++) out(I_tens->i_pair_mu_[k], I_tens->i_pair_nu_[k]) = dev[k];
		on_device = true;
	}
#endif
	// See calc_F_calc: an exception must not leave an OpenMP structured block.
	std::string io_error;
	//One region, one accumulator per thread, one reduction, however many windows the
	//budget implies: an nmo x nmo matrix per window is what made narrow windows dear
	std::vector<occ::Mat> parts(omp_get_max_threads());
#pragma omp parallel if (!on_device)
	{
		occ::Mat local = occ::Mat::Zero(sf->model_data.nmo, sf->model_data.nmo);
		double* local_ptr = local.data();
		for (int r0 = 0; !on_device && r0 < sf->model_data.nr; r0 += step) {
			const int r1 = std::min(r0 + step, sf->model_data.nr);
			if (I_tens->i_streamed_) {
#pragma omp single
				{
					try { I_tens->i_file_.load(r0, r1); }
					catch (const std::exception& e) { io_error = e.what(); }
				}
			}
			//No nowait: the next window's read must not start until every thread has read this one
#pragma omp for
			for (int r = r0; r < r1; r++) {
				if (!io_error.empty()) continue;
				const cdouble precompute = pre[r];
				if (precompute == cdouble(0.0, 0.0)) continue;
				//As in calc_F_calc: one walk over whichever element type is resident, with
				//the accumulation in double either way.
				auto accumulate = [&](const auto* I_r) {
					const int* pmu = I_tens->i_pair_mu_.data();
					const int* pnu = I_tens->i_pair_nu_.data();
					for (size_t k = 0; k < I_tens->i_compact_; k++) {
						const double vr = static_cast<double>(I_r[k].real());
						const double vi = static_cast<double>(I_r[k].imag());
						local_ptr[pnu[k] * sf->model_data.nmo + pmu[k]] += precompute.real() * vr - precompute.imag() * vi;
					}
				};
				if (I_tens->i_float_) accumulate(I_tens->i_block32(r)); else accumulate(I_tens->i_block(r));
			}
		}
		parts[omp_get_thread_num()].swap(local);
	}
	if (!io_error.empty()) throw std::runtime_error(io_error);
	//The partials outweigh the matrix many times over, so their sum is parallel too
	if (!on_device) {
#pragma omp parallel for schedule(static)
		for (int mu = 0; mu < sf->model_data.nmo; mu++) {
			for (int nu = mu; nu < sf->model_data.nmo; nu++) {
				double sum = 0.0;
				for (int t = 0; t < static_cast<int>(parts.size()); t++)
					if (parts[t].size() != 0) sum += parts[t](mu, nu);
				out(mu, nu) = sum;
			}
		}
	}
	for (int mu = 0; mu < sf->model_data.nmo; mu++) {
		for (int nu = mu + 1; nu < sf->model_data.nmo; nu++) {
			out(nu, mu) = out(mu, nu);
		}
	}
}

//The sign of f(+-3), g(+-3) and g(+-4) in every orbital, in OCC's own m = -l..l order,
//and the density matrix rebuilt from them. Applied once on the way out and once on the
//way back in, it is its own inverse.
void XCW_solver::flip_high_m_phases(occ::qm::Wavefunction& w) {
	int row = 0;
	const int spins = w.mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
	for (const auto& shell : w.basis.shells()) {
		const int nsph = 2 * shell.l + 1;
		if (shell.l >= 3)
			for (int spin = 0; spin < spins; spin++) {
				const int base = spin * w.nbf + row;
				w.mo.C.row(base).array() *= -1.0;
				w.mo.C.row(base + nsph - 1).array() *= -1.0;
				if (shell.l >= 4) {
					w.mo.C.row(base + 1).array() *= -1.0;
					w.mo.C.row(base + nsph - 2).array() *= -1.0;
				}
			}
		row += nsph;
	}
	w.mo.update_occupied_orbitals();
	w.mo.update_density_matrix();
}

void XCW_solver::create_tscb(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda) {
	sf->XCW_log << "Creating .tscb file from converged SCF calculation..." << std::endl;
	std::vector<WFN> sf_wave_vec(1, { scf.wavefunction(), false });
	//The constructor marks anything taken from OCC as OCC-origin. What this refinement holds
	//is an OCC result over a basis this program loaded, so say that: Int_Params then reads the
	//shells with the convention they actually have.
	sf_wave_vec[0].set_origin(e_origin::XCW_fit);
	svec known_atoms_;
	tsc_block<int, cdouble> result;
	vec2 known_kpts_;
	opt->m_hkl_list = sf->scatter_data.hkl_enlarged;
	opt->grid_cache = &tsc_grids;
	result.append(calculate_scattering_factors<itsc_block, std::vector<WFN>&>(
		*opt,
		sf_wave_vec,
		sf->XCW_log,
		known_atoms_,
		0,
		&sf->k_pt),
		sf->XCW_log);
	std::string value = std::to_string(lambda);
	value.erase(std::remove(value.begin(), value.end(), '.'), value.end());
	while (value.length() < 7) {
		value += '0';
	}
	std::ostringstream oss;
	oss << "NA2_" << value << ".tscb";
	result.write_tscb_file("test.cif", oss.str());
	std::ostringstream oss2;
	oss2 << "NA2_" << value << ".wfn";
	sf_wave_vec[0].write_wfn(oss2.str(), false, true);
	//the structure factors this step was scored on, for an R-factor against another route
	//(an ORCA tsc through Olex2, or a Laue check between symmetry-equivalent rows)
	{
		sf->ensure_hkl_ordered();
		std::ofstream fc("NA2_" + value + "_Fcalc.txt");
		if (!sf->ext_p_.empty()) {
			std::cout << "XCW lambda " << lambda << ": " << sf->extinction_report() << std::endl;
			sf->XCW_log << "XCW lambda " << lambda << ": " << sf->extinction_report() << std::endl;
		}
		fc << "#    h    k    l          F_obs        sig(F)   scale*|F_calc|     phase(deg)   R1(gt) = " << std::setprecision(5) << sf->quality_criteria.R1 << " R1(all) = " << sf->quality_criteria.R1_all << " scale = " << std::setprecision(10) << sf->scatter_data.scale << "\n";
		for (int r = 0; r < sf->model_data.nr; r++) {
			const cdouble& f = sf->scatter_data.F_calc[r];
			fc << std::setw(5) << sf->hkl_ordered_[r][0] << std::setw(5) << sf->hkl_ordered_[r][1] << std::setw(5) << sf->hkl_ordered_[r][2]
				<< std::fixed << std::setprecision(4) << std::setw(15) << sf->scatter_data.F_obs[r] << std::setw(14) << sf->scatter_data.sigma_obs[r]
				<< std::setw(17) << sf->scatter_data.scale * sf->ext_sqrt_y(r) * std::abs(f) << std::setw(15) << std::arg(f) * 180.0 / constants::PI << "\n";
		}
	}
	std::ostringstream oss3;
	oss3 << "NA2_" << value << ".fchk";
	//OCC's fchk writer reorders to Gaussian's basis functions but keeps libcint's phases,
	//and Gaussian's f(+-3), g(+-3), g(+-4) are the opposite sign; the file has to carry
	//Gaussian's so that anything reading an fchk gets the density right
	{
		occ::qm::Wavefunction w = scf.wavefunction();
		flip_high_m_phases(w);
		w.save(oss3.str());
	}
	if (opt->xcw_settings.nbo_output) {
		std::ostringstream oss4;
		oss4 << "NA2_" << value << ".47";
		sf_wave_vec[0].write_nbo(oss4.str(), opt->debug, &sf->XCW_log);
	}
	//Neither file written above can carry this analysis - a .wfn has bare primitives and
	//the fchk reader keeps no shells - so -rgbi runs it here, on the refined wavefunction,
	//and its report goes to a file of its own per lambda
	if (opt->rgbi) {
		std::ostringstream oss5;
		oss5 << "NA2_" << value << "_RGBI.txt";
		std::ofstream rgbi_out(oss5.str());
		std::streambuf* const cout_buf = std::cout.rdbuf(rgbi_out.rdbuf());
		Roby_information Roby(sf_wave_vec[0], opt->rgbi_group_sets, !opt->rgbi_no_sym,
			opt->rgbi_orbital_basis == RGBIOrbitalBasis::ANO, opt->rgbi_EVs, opt->rgbi_theta,
			opt->rgbi_legacy_cutoff);
		std::cout.rdbuf(cout_buf);
		sf->XCW_log << "RGBI analysis written to " << oss5.str() << std::endl;
	}
}

occ::qm::HartreeFock XCW_solver::setup_XCW_procedure() {
	std::vector<structure_factors::ao_data> ao_data_shells;
	occ::core::Molecule mol;
	scf_solver.setup_SCF_mol(mol);
	occ::qm::AOBasis occ_basis_set;
	std::string temp_basis = opt->xcw_settings.basis_set_name;
	scf_solver.setup_basis(mol, temp_basis, occ_basis_set);
	occ::qm::HartreeFock hf(occ_basis_set);
	if (!opt->xcw_settings.df_basis_name.empty()) {
		//OCC loads a fitting basis by name from a data directory this build does not ship,
		//but reads a .json path as given: the library's set is written out once, for the
		//elements present, and handed over that way
		std::shared_ptr<BasisSet> aux = BasisSetLibrary::get_basis_set(opt->xcw_settings.df_basis_name);
		//a Coulomb-only set fits J and leaves the exchange to a basis never meant for it
		if (aux->get_name().find("jkfit") == std::string::npos)
			std::cout << "WARNING: " << aux->get_name() << " is not a JK-fitting basis; Hartree-Fock exchange is fitted with it all the same. "
				"def2-universal-jkfit serves the def2 family, cc-pvXz-jkfit the cc-pVXZ family." << std::endl;
		ivec elements;
		for (int i = 0; i < static_cast<int>(mol.atoms().size()); i++)
			if (std::find(elements.begin(), elements.end(), mol.atoms()[i].atomic_number) == elements.end())
				elements.push_back(mol.atoms()[i].atomic_number);
		const std::string file = aux->get_name() + "_df.json";
		aux->write_occ_json(file, elements);
		hf.set_density_fitting_basis(file);
		//OCC keeps the three-index integrals only under a 512 MB limit and otherwise recomputes
		//them every iteration, which costs twice a direct build here. Held whenever they fit in
		//half of what the process can have; they are computed once for the whole lambda scan.
		const size_t naux = aux->to_AOBasis(mol.atoms()).nbf();
		const size_t nbf = occ_basis_set.nbf();
		const size_t store = naux * nbf * (nbf + 1) / 2 * sizeof(double);
		const size_t avail = available_memory_bytes();
		const bool stored = avail == 0 || store < avail / 2;
		hf.set_density_fitting_policy(stored ? occ::qm::IntegralEngineDF::Policy::Stored : occ::qm::IntegralEngineDF::Policy::Direct);
		std::cout << "XCW density fitting with " << aux->get_name() << ": " << naux << " functions, "
			<< (store / 1048576.0) << " MB of three-index integrals " << (stored ? "held in memory" : "recomputed every iteration") << std::endl;
	}
	if (opt->xcw_int_precision > 0.0) hf.set_precision(opt->xcw_int_precision);
	sf->create_prims(ao_data_shells, occ_basis_set);
	I_tens = &sf->eval_I_anom_disp(ao_data_shells);
	return hf;
	// closing function
}

void XCW_solver::run() {
	//OCC parallelises through TBB, which does not read OMP_NUM_THREADS. Not a speedup - it
	//already used every core - but it makes -cpus bind the 82% of a run that OCC owns.
	occ::parallel::set_num_threads(opt->threads > 0 ? opt->threads : omp_get_max_threads());
	occ::qm::HartreeFock hf = setup_XCW_procedure();
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	//The two walks of an iteration read the whole tensor and the host is bound by its memory
	//bandwidth doing so; the device reads it several times faster. The host copy stays for
	//the background writer.
	if (opt->gpu_itensor && opt->use_gpu && !I_tens->i_streamed_) {
		I_tens->i_on_device_ = I_tens->i_float_ ? itensor_gpu_hold(I_tens->I32.data(), sf->model_data.nr, static_cast<int>(I_tens->i_compact_))
			: itensor_gpu_hold(I_tens->I.data(), sf->model_data.nr, static_cast<int>(I_tens->i_compact_));
		if (!(opt->no_date))
			std::cerr << "GPU in use: XCW structure factors and perturbation on "
			<< (I_tens->i_on_device_ ? "the device" : "the CPU - device unavailable or the tensor too large") << std::endl;
	}
#endif
	//Four fifths of the free memory, the I tensor's budget, for the stored integrals
	scf_solver.eri_.clear();
	if (opt->xcw_settings.hf_type != occ::qm::SpinorbitalKind::General) {
		const _time_point eri_t0 = get_time();
		const size_t avail = available_memory_bytes();
		if (scf_solver.eri_.build(hf, avail ? avail / 5 * 4 : 0, sf->XCW_log))
			throughput::record_time("XCW two-electron integrals", false, get_msec(eri_t0, get_time()));
	}
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (scf_solver.eri_ && opt->gpu_itensor && opt->use_gpu) {
		const auto up_t0 = get_time();
		scf_solver.eri_on_device_ = eri_gpu_hold(scf_solver.eri_.data(), scf_solver.eri_.nbf(), scf_solver.eri_.npairs(), scf_solver.eri_.pair_a().data(), scf_solver.eri_.pair_b().data(),
			scf_solver.eri_.first_pair().data(), scf_solver.eri_.pair_index().data());
		if (scf_solver.eri_on_device_) throughput::record_time("XCW two-electron integrals upload", true, get_msec(up_t0, get_time()));
		if (!(opt->no_date))
			std::cerr << "GPU in use: XCW Fock build from the stored integrals on "
			<< (scf_solver.eri_on_device_ ? "the device" : "the CPU - device unavailable or the integrals too large") << std::endl;
	}
#endif
	bool has_guess = false;
	occ::qm::Wavefunction last_wfn, prev_wfn;
	const occ::Mat S_ao = hf.compute_overlap_matrix();
	if (opt->xcw_settings.read_first_guess) {
		std::ostringstream oss2;
		std::string start_value_str = std::to_string(opt->xcw_settings.xcw_start_value);
		start_value_str.erase(std::remove(start_value_str.begin(), start_value_str.end(), '.'), start_value_str.end());
		if (start_value_str.length() > 7) {
			start_value_str = start_value_str.substr(0, 7);
		}
		else if (start_value_str.length() < 7) {
			start_value_str.append(7 - start_value_str.length(), '0');
		}
		oss2 << "NA2_" << start_value_str << ".fchk";
		last_wfn = occ::qm::Wavefunction::load(oss2.str());
		//Written in Gaussian's phases by create_tscb; OCC's loader does not undo that
		flip_high_m_phases(last_wfn);
		has_guess = true;
	}

	std::cout << "More detailed output in XCW.log file..." << std::endl;
	if (opt->xcw_settings.XWR_type == 2) {
		std::cout << "XCW: fitting against the 1/|H|^2-weighted residual self-energy criterion "
			<< "Criterion below is this weighted quantity, not the classical GoF." << std::endl;
		sf->XCW_log << "XCW: fitting against the 1/|H|^2-weighted residual self-energy criterion "
			<< "Criterion below are this weighted quantity, not the classical GoF." << std::endl;
	}
	std::cout << "____________________________________________________________________________________\n";
	std::cout << " Lambda\t\tCriterion\tGooF(F2)\tR1(gt)\t\tTotal Energy\t\tPerturbation\tTarget quantity\t\tCrit(all)\tR1(all)";
	if (opt->xcw_settings.xcw_gaussian_halt) {
		std::cout << "\t\tA^2 (halt)";
	}
	std::cout << "\n";
	std::cout << "\t\t\t\t\t\t\t\t(Eh)\t\t\t(a. u.)\t\t(a. u.)\n";
	std::cout << "____________________________________________________________________________________\n";

	// Runs the lambda steps for XCW fitting
	auto run_lambda = [&](const double lambda, occ::qm::Wavefunction guess, const bool use_guess, const bool write_result) {
		occ::qm::SCF scf(hf, opt->xcw_settings.hf_type);
		double alpha = opt->xcw_settings.alpha;
		bool has_local_guess = use_guess;
		scf.set_charge_multiplicity(opt->xcw_settings.charge, opt->xcw_settings.multiplicity);
		scf.maxiter = opt->xcw_settings.max_scf_iterations;
		scf.convergence_settings.level_shift = opt->xcw_settings.level_shift;
		scf.convergence_settings.level_shift_threshold = 0;
		scf.update_occupied_orbital_count();
		const bool converged = scf_solver.do_SCF(lambda, alpha, scf, guess, has_local_guess, write_result);
		return std::make_pair(converged, scf.wavefunction());
	};
	const double min_lambda_step = opt->xcw_settings.xcw_step_size / 128.0;
	double last_lambda = opt->xcw_settings.xcw_start_value;
	//slow_conv's damping (alpha 0.8, shift 1) is for the perturbed steps, which start from
	//converged orbitals; from the core guess it runs out of iterations (Fe(phen)2(SCN)2 UHF:
	//200 iterations, RMSD still 4e-5). The first step without a guess takes the normal
	//schedule, every later one the slow settings from the file
	const bool slow_start = opt->xcw_settings.slow_conv && !has_guess;
	const double slow_alpha = opt->xcw_settings.alpha, slow_shift = opt->xcw_settings.level_shift, slow_stop_damping = opt->xcw_settings.diis_stop_damping, slow_stop_shift = opt->xcw_settings.diis_stop_shift;
	if (slow_start) {
		opt->xcw_settings.alpha = 0.5; opt->xcw_settings.level_shift = 0.5; opt->xcw_settings.diis_stop_damping = 1e-3; opt->xcw_settings.diis_stop_shift = 1e-2;
		std::cout << "XCW: slow_conv - the unperturbed first step runs the normal schedule, slow damping from the second step on" << std::endl;
		sf->XCW_log << "XCW: slow_conv - the unperturbed first step runs the normal schedule, slow damping from the second step on" << std::endl;
	}
	//The scan used to stop at max_value even when A^2 was still falling there, i.e. on the
	//scan boundary instead of on lambda*, and only printed "extend the scan" afterwards.
	//The perturbed steps start from converged orbitals and TRAH takes the ones that DIIS
	//cannot, so running past the requested range is stable enough to just do it. Capped at
	//the requested number of steps, so a trend that never curves upward at most doubles the
	//scan rather than running away
	int planned_steps = opt->xcw_settings.num_xcw_steps;
	const int max_extra_steps = opt->xcw_settings.num_xcw_steps;
	int extra_steps = 0;
	for (int step = 0; step < planned_steps; step++) {
		const double lambda = step * opt->xcw_settings.xcw_step_size + opt->xcw_settings.xcw_start_value;
		if (slow_start && step == 1) {
			opt->xcw_settings.alpha = slow_alpha; opt->xcw_settings.level_shift = slow_shift; opt->xcw_settings.diis_stop_damping = slow_stop_damping; opt->xcw_settings.diis_stop_shift = slow_stop_shift;
		}
		const occ::qm::Wavefunction previous_wfn = last_wfn;
		occ::qm::Wavefunction guess = last_wfn;
		if (opt->xcw_extrapolate && step >= 2 && opt->xcw_settings.hf_type == occ::qm::SpinorbitalKind::Restricted) {
			//The density extrapolated through the two previous steps, pulled back to
			//idempotency by two McWeeny steps D <- 3DSD - 2DSDSD; the orbitals stay those of
			//the last step, they only seed the level shift and the gradient
			occ::Mat D = 2.0 * last_wfn.mo.D - prev_wfn.mo.D;
			for (int k = 0; k < 2; k++) {
				const occ::Mat DS = D * S_ao;
				D = 3.0 * DS * D - 2.0 * DS * DS * D;
			}
			guess.mo.D = D;
		}
		auto result = run_lambda(lambda, guess, has_guess, true);
		if (!result.first && step > 0) {
			double lambda_step = lambda - last_lambda;
			while (!result.first && lambda_step > min_lambda_step) {
				lambda_step *= 0.5;
				const double trial_lambda = std::min(lambda, last_lambda + lambda_step);
				sf->XCW_log << "XCW: retrying lambda " << std::fixed << std::setprecision(8) << trial_lambda
					<< " from converged lambda " << last_lambda << " with step " << lambda_step << std::endl;
				std::cout << "XCW: retrying lambda " << std::fixed << std::setprecision(8) << trial_lambda
					<< " with step " << lambda_step << std::endl;
				auto trial = run_lambda(trial_lambda, last_wfn, true, trial_lambda == lambda);
				if (trial.first) {
					if (trial_lambda == lambda) {
						result = std::move(trial);
						break;
					}
					last_lambda = trial_lambda;
					last_wfn = trial.second;
					lambda_step = std::min(2.0 * lambda_step, lambda - last_lambda);
					result = run_lambda(lambda, last_wfn, true, true);
				}
			}
		}
		if (!result.first) {
			//step 0 has no converged neighbour to continue from, so the halving loop never ran
			std::ostringstream why;
			why << "XCW: unable to converge lambda " << std::fixed << std::setprecision(8) << lambda;
			if (step == 0)
				why << " in " << opt->xcw_settings.max_scf_iterations << " SCF iterations (raise max_iter or loosen the criteria); stopping scan.";
			else
				why << " with a continuation step above " << min_lambda_step << "; stopping scan.";
			sf->XCW_log << why.str() << std::endl;
			std::cout << why.str() << std::endl;
			break;
		}
		prev_wfn = previous_wfn;
		last_wfn = result.second;
		last_lambda = lambda;
		has_guess = true;

		//At the end of the planned scan, keep going while the minimum is still ahead
		if (opt->xcw_settings.xcw_gaussian_halt && step + 1 == planned_steps) {
			double estimated_minimum = 0.0;
			if (halting_minimum_beyond_scan(gaussian_halt_history_, estimated_minimum)) {
				std::ostringstream msg;
				msg << std::fixed << std::setprecision(5)
					<< "XCW: A^2 is still falling at lambda " << lambda << ", so the minimum is outside the requested range - ";
				if (extra_steps < max_extra_steps) {
					planned_steps++;
					extra_steps++;
					msg << "extending the scan to lambda " << (lambda + opt->xcw_settings.xcw_step_size)
						<< " (extra step " << extra_steps << " of at most " << max_extra_steps << ")";
					if (estimated_minimum > 0.0) {
						msg << ", extrapolated minimum near lambda ~= " << estimated_minimum;
					}
				}
				else {
					msg << "stopping anyway, the scan has already been extended by its limit of "
						<< max_extra_steps << " steps; raise max_value in -do_XCW to continue";
				}
				sf->XCW_log << msg.str() << std::endl;
				std::cout << msg.str() << std::endl;
			}
		}

		//Progress estimate every 5 lambda steps; the last step is skipped because the
		//summary below always prints a final one
		if (opt->xcw_settings.xcw_gaussian_halt && (step + 1) % 5 == 0 && step + 1 < planned_steps) {
			report_halting_progress_estimate(false);
		}
	}

	if (opt->xcw_settings.xcw_gaussian_halt) {
		report_gaussian_halting_summary();
	}

	//Before the run ends: the writer holds a file handle and reads the resident tensor, and
	//by now it has usually been finished for a long while - the refinement takes far longer
	//than the write. Joining is what keeps it from outliving the process.
	sf->finish_i_save();
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	itensor_gpu_release();
	eri_gpu_release();
	I_tens->i_on_device_ = scf_solver.eri_on_device_ = false;
#endif

	std::cout << "Finished XCW fitting procedure." << std::endl;
}
