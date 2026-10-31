#include "pch.h"
#include "XCW_wrapper.h"
#ifdef NOSPHERA2_USE_GPU
#include "xcw_grid_gpu.h"
#endif
#include "scattering_factors.h"
#include "basis_set.h"
#include "gaussian_atom.h"
#include "SALTED_utilities.h"
#include "bondwise_analysis.h"
#include "tuning.h"
#include <limits>

XCW_wrapper::XCW_wrapper(options& opt_in) : opt(&opt_in), sf(opt_in) {
	std::cout << "XCW orbital basis set: " << sf.orbital_basis->get_name() << std::endl;
	writer.WAEmpiRe();
	writer.log << "XCW orbital basis set: " << sf.orbital_basis->get_name() << std::endl;
	writer.log << sf.model_summary << std::endl;
}

//z_h = (|F_obs,h| - |F_calc,h|) / sigma_h for the converged F_calc/F_scale over
//strong reflections (|F_obs|/sigma >= xcw_strong_cutoff), tested against
//N(0,1) after a global shape/scale rescale. Full reflection set only; the
//free/working-set cross-validation variant is not implemented.
//Background: tests/P1_test/XCW_plan.md, Src/core/xcw_halting.h.
void XCW_wrapper::evaluate_gaussian_halting(const double lambda) {
	sf.ensure_hkl_ordered();

	const int n_total = sf.model_data.nr;
	vec z_raw;
	vec resolution;
	vec abs_F;
	z_raw.reserve(n_total);
	resolution.reserve(n_total);
	abs_F.reserve(n_total);

	for (int i = 0; i < n_total; i++) {
		if (sf.scatter_data.sigma_obs[i] <= 0.0) {
			continue;
		}
		const double f_over_sigma = sf.scatter_data.abs_F_obs[i] / sf.scatter_data.sigma_obs[i];
		if (f_over_sigma < opt->xcw_settings.xcw_strong_cutoff) {
			continue;
		}
		const double scaled_F_calc = sf.scatter_data.scale * std::abs(sf.scatter_data.F_calc[i]);
		const double diff = scaled_F_calc - sf.scatter_data.abs_F_obs[i];
		z_raw.push_back(diff / sf.scatter_data.sigma_obs[i]);
		resolution.push_back(sf.unit_cell.get_stl_of_hkl(sf.hkl_ordered_[i]));
		abs_F.push_back(sf.scatter_data.abs_F_obs[i]);
	}

	GaussianHaltEntry entry;
	entry.lambda = lambda;
	entry.n_total = n_total;
	entry.n_used = static_cast<int>(z_raw.size());

	if (entry.n_used < 8) {
		writer.log << "Gaussian halting criterion: only " << entry.n_used
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

	writer.log << "Gaussian halting criterion at lambda=" << std::fixed << std::setprecision(5) << lambda << ":\n"
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

//Full per-lambda table to the SCF log, then the final recommendation
void XCW_wrapper::report_gaussian_halting_summary() {
	if (gaussian_halt_history_.empty()) {
		return;
	}

	writer.log << "\n____________________________________________________________________________\n"
		<< "Gaussian halting criterion summary (tests/P1_test/XCW_plan.md)\n"
		<< " Lambda\t\tA^2\treject5%\tpp_slope\tpp_intercept\tskew\tkurt\tres_trend_r\tint_trend_r\tn_used\n";
	for (const GaussianHaltEntry& e : gaussian_halt_history_) {
		writer.log << "\t" << std::fixed << std::setprecision(5) << e.lambda
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
//test would miss, so it is warned about. See XCW_wrapper.h for the full behaviour.
void XCW_wrapper::report_halting_progress_estimate(bool is_final) {
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

	std::ostream* streams[2] = { &writer.log, &std::cout };
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

//lambda as the result files name it: its digits without the point, padded to seven (0.01 -> 0010000)
static std::string lambda_tag(const double lambda) {
	std::string value = std::to_string(lambda);
	value.erase(std::remove(value.begin(), value.end(), '.'), value.end());
	while (value.length() < 7) value += '0';
	return value;
}

//Guess for the next lambda step: the last orbitals rotated on by t times the rotation that led
//to them from the step before. In the last step's MO frame the previous occupied projector has
//the occupied-virtual block P(a,i) ~ kappa(a,i) of that rotation, so C exp(-t kappa) by the
//Cayley transform, per spin. C and D stay consistent, which TRAH needs: it steps from the orbitals
static void extrapolate_orbitals(occ::qm::MolecularOrbitals& mo, const occ::qm::MolecularOrbitals& prev, const occ::Mat& S, const double t) {
	const Eigen::Index n = static_cast<Eigen::Index>(mo.n_ao), m = mo.C.cols();
	const int nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
	for (int b = 0; b < nb; b++) {
		const Eigen::Index nocc = static_cast<Eigen::Index>(b == 0 ? mo.n_alpha : mo.n_beta);
		const occ::Mat Cb = mo.C.middleRows(b * n, n);
		const occ::Mat X = Cb.transpose() * S * prev.C.middleRows(b * n, n).leftCols(nocc);
		const occ::Mat Pvo = X.bottomRows(m - nocc) * X.topRows(nocc).transpose();
		occ::Mat K = occ::Mat::Zero(m, m);
		K.bottomLeftCorner(m - nocc, nocc) = -0.5 * t * Pvo;
		K.topRightCorner(nocc, m - nocc) = 0.5 * t * Pvo.transpose();
		const occ::Mat I = occ::Mat::Identity(m, m);
		mo.C.middleRows(b * n, n) = Cb * (I - K).partialPivLu().solve(I + K);
	}
	mo.update_occupied_orbitals();
	mo.update_density_matrix();
}

//The spin blocks of a density matrix summed, as the structure factors take it
static dMatrix2 spin_summed(const occ::Mat& Dw) {
	const int n = static_cast<int>(Dw.cols());
	const bool unrestricted = Dw.rows() == 2 * Dw.cols();
	dMatrix2 D(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) D(i, j) = Dw(i, j) + (unrestricted ? Dw(i + n, j) : 0.0);
	return D;
}

void XCW_wrapper::evaluate_cryst_crit(const dMatrix2& dm_eff) {
	//density-following shares refit at this evaluation only when NOSPHERA2_XCW_WEIGHT_UPDATE asks for it
	sf.XG.refresh = sf.XG.ws && scf_solver.lambda != 0.0 && (weight_every_ > 0 ? weight_count_ % weight_every_ == 0 : weight_every_ == 0 && weight_count_ == 0);
	weight_count_++;
	//between solves refits also follow the density inside one while it moves by more than NOSPHERA2_XCW_WEIGHT_DTOL, and hold
	//still for the rest of the solve once it moves less, so that TRAH near convergence sees one objective
	if (sf.XG.ws && scf_solver.lambda != 0.0 && weight_every_ < 0 && weight_dtol_ > 0.0 && !weight_hold_) {
		const int n = static_cast<int>(dm_eff.extent(0));
		double d = static_cast<int>(weight_D_fit_.extent(0)) == n ? 0.0 : 1e300;
		for (int i = 0; i < n && d <= weight_dtol_; i++)
			for (int j = 0; j < n; j++) d = std::max(d, std::abs(dm_eff(i, j) - weight_D_fit_(i, j)));
		sf.XG.refresh = d > weight_dtol_;
		weight_hold_ = !sf.XG.refresh && weight_moved_;
		weight_moved_ = weight_moved_ || sf.XG.refresh;
		//a refit mid-solve needs no more than the density it follows has: 1e-2 of that move, between the tight 1e-6 and 1e-3
		sf.XG.tol = std::min(1e-3, 1e-2 * d);
	}
	if (sf.XG.refresh) weight_D_fit_ = dm_eff;
	sf.calc_F_calc(dm_eff);
	sf.XG.refresh = false;
	sf.XG.tol = 0.0;
	sf.eval_scale();
	sf.calc_criteria();
}

void XCW_wrapper::perturbation_matrix(occ::Mat& perturb, const bool unrestricted) {
	sf.perturbation(perturb, unrestricted, writer.log);
}

occ::Mat XCW_wrapper::perturbation_matrix_response(const dMatrix2& dD_eff, const double lambda) {
	return sf.perturbation_response(dD_eff, lambda);
}

fit_quality XCW_wrapper::cryst_crit() const {
	return { sf.criterion(false), sf.quality_criteria.GooF2, sf.quality_criteria.R1, sf.criterion(true), sf.quality_criteria.R1_all };
}

void XCW_wrapper::save_F_calc() {
	sf.save_F_calc();
}

void XCW_wrapper::restore_F_calc() {
	sf.restore_F_calc();
}

bool XCW_wrapper::grid_on_device() const {
	return sf.XG.on_device;
}

vec XCW_wrapper::ri_tshb_coefficients() {
	const occ::qm::AOBasis& A = sf.XR.aux;
	if (tshb_map_.empty()) {
		const aux_density_table t(*tshb_aux_.get_atoms_ptr());
		if (t.n_sh != static_cast<int>(A.size()) || t.n_coef != sf.XR.naux) return vec();
		std::mt19937 gen(5);
		std::normal_distribution<double> nd;
		vec unit(t.n_coef, 0.0);
		const int np = 48;
		for (int s = 0; s < t.n_sh; s++) {
			const int l = t.sh_l[s], m = 2 * l + 1, o = t.coef_off[s];
			if (static_cast<int>(A[s].l) != l || static_cast<int>(A.first_bf()[s]) != o) { tshb_map_.clear(); return vec(); }
			//points where the shell is not negligible: the spread of its most diffuse primitive
			double amin = A[s].exponents(0);
			for (int q = 1; q < static_cast<int>(A[s].num_primitives()); q++) amin = std::min(amin, A[s].exponents(q));
			occ::Mat3N p(3, np);
			for (int i = 0; i < np; i++)
				for (int d = 0; d < 3; d++) p(d, i) = A[s].origin(d) + nd(gen) / std::sqrt(2.0 * amin);
			occ::Mat N(np, m);
			for (int k = 0; k < m; k++) {
				unit[o + k] = 1.0;
				for (int i = 0; i < np; i++) N(i, k) = t(p(0, i), p(1, i), p(2, i), unit.data());
				unit[o + k] = 0.0;
			}
			const occ::Mat O = occ::gto::evaluate_basis(A, p, 0).phi.middleCols(o, m);
			//phi_occ = phi_nos X, so rho = phi_occ c = phi_nos (X c)
			const occ::Mat X = N.colPivHouseholderQr().solve(O);
			if ((N * X - O).cwiseAbs().maxCoeff() > 1e-9 * O.cwiseAbs().maxCoeff()) { tshb_map_.clear(); return vec(); }
			tshb_map_.push_back(X);
		}
	}
	const aux_density_table t(*tshb_aux_.get_atoms_ptr());
	vec c(t.n_coef, 0.0);
	for (int s = 0; s < t.n_sh; s++) {
		const int m = 2 * t.sh_l[s] + 1, o = t.coef_off[s];
		for (int i = 0; i < m; i++)
			for (int k = 0; k < m; k++) c[static_cast<size_t>(o) + i] += tshb_map_[s](i, k) * sf.XR.c[static_cast<size_t>(o) + k];
	}
	return c;
}

//The fit is another partition than the table's, so the two differ at the 1e-2 level in F; -ri_fit names
//the aux basis, AutoAux otherwise
void XCW_wrapper::write_tshb(WFN& wave, const tsc_block<int, cdouble>& table, const std::string& name) {
	if (tshb_aux_.get_ncen() == 0) {
		std::vector<std::shared_ptr<BasisSet>> aux = opt->aux_basis;
		if (aux.empty())
			aux.push_back(std::make_shared<BasisSet>(wave));
		//on the RI partition the basis it fits in, so its coefficients carry over
		if (sf.XR.on) aux = { sf.XR.aux_set };
		tshb_aux_ = generate_aux_wfn(wave, aux, false);
	}
	DensityFitting::CONFIG config = DensityFitting::config_from_options(*opt);
	config.partition_cache = &tshb_cache_;
	//the RI partition's own fit when it runs, so the .tshb and the scored F share one fit; NOSPHERA2_XCW_TSHB_FIT fits anew
	vec c = sf.XR.on && !tuning("NOSPHERA2_XCW_TSHB_FIT") ? ri_tshb_coefficients() : vec();
	const bool own = !c.empty();
	if (!own) c = DensityFitting::density_fit(wave, tshb_aux_, config);
	const Gaussian_Molecule ml(tshb_aux_, c);
	const ScattererLabels& keys = table.get_scatterers();
	const svec names = table.get_scatterers_string();
	//the scatterers keep their order over the steps, and only the first step's wavefunction went through
	//the CIF reader that gives its atoms their ids
	if (tshb_atoms_.empty()) {
		tshb_atoms_.assign(keys.size(), -1);
		for (int scatterer = 0; scatterer < keys.size(); scatterer++) {
			for (int a = 0; a < wave.get_ncen() && tshb_atoms_[scatterer] < 0; a++)
				if (wave.get_id_for_atom(a).to_hex_string() == names[scatterer] || wave.get_atom_label(a) == names[scatterer])
					tshb_atoms_[scatterer] = a;
			err_checkf(tshb_atoms_[scatterer] >= 0, "-tshb: no atom of the XCW wavefunction for scatterer " + names[scatterer], writer.log);
		}
	}
	svec labels(keys.size());
	for (int scatterer = 0; scatterer < keys.size(); scatterer++)
		labels[scatterer] = sf.asym_atoms[tshb_atoms_[scatterer]].label;
	ml.tshb(labels, keys, tshb_atoms_, wave, wave.get_ECP_mode(), opt->electron_diffraction).write(name);
	writer.log << "Wrote " << name << ": " << keys.size() << " atom(s), " << (own ? "the RI partition's fit" : "RI fit of the converged density") << std::endl;
}

//The files of a converged lambda step, all named NA2_<lambda>: the scattering factor table (.tscb, or the
//.tshb fit with -tshb), the wavefunction (.wfn, .fchk, and .47 for NBO when asked), with -print_F the F_calc it was
//scored on and, with -rgbi, the RGBI analysis of the refined wavefunction
void XCW_wrapper::write_step_results(const occ::qm::Wavefunction& wfn, const double lambda) {
	//what the step's writers print (tscb check, AutoAux, the .tshb fit) goes to SCF.log, so the console keeps
	//only the lambda table; fatal errors still reach the console, error_check restores cout first
	struct cout_restore { std::streambuf* b; ~cout_restore() { std::cout.rdbuf(b); } } console{ std::cout.rdbuf(writer.log.rdbuf()) };
	writer.log << "Creating .tscb file from converged SCF calculation..." << std::endl;
	std::vector<WFN> sf_wave_vec(1, { wfn, false });
	//The constructor marks anything taken from OCC as OCC-origin. What this refinement holds
	//is an OCC result over a basis this program loaded, so say that: Int_Params then reads the
	//shells with the convention they actually have.
	sf_wave_vec[0].set_origin(e_origin::XCW_fit);
	const _time_point fwd_t0 = get_time();
	//the scattering factor table of the converged density: by the tsc route at the first step, from the grid after it
	const itsc_block result = sf.result_table(wfn.mo.D, sf_wave_vec, writer.log, weight_every_ > 0 && lambda != 0.0);
	//the atomic charges of density-following shares, from the partitioned rho of the table
	if (sf.XG.ws && (!sf.XR.on || sf.XR.lmax >= 0)) {
		writer.log << "XCW grid charges at lambda " << std::fixed << std::setprecision(5) << lambda << ":" << std::setprecision(6);
		for (int a = 0; a < sf.model_data.ncen; a++) {
			double e = 0;
			for (int p = 0; p < static_cast<int>(sf.XG.rho_s[a].size()); p++) e += sf.XG.rho_s[a][p];
			writer.log << " " << sf.asym_atoms[a].label << " " << sf.asym_atoms[a].type - e;
		}
		writer.log << std::endl;
	}
	//NOSPHERA2_XCW_GRID_FCHECK: F_calc from the scattering factors of this table against the one the SCF
	//used at the same density
	if (tuning("NOSPHERA2_XCW_GRID_FCHECK")) {
		const double fwd_ms = get_msec(fwd_t0, get_time());
		const cvec F_I = sf.scatter_data.F_calc;
		sf.calc_F_calc(result);
		double dmax = 0.0, fmax = 0.0, num = 0.0, den = 0.0;
		for (size_t r = 0; r < F_I.size(); r++) {
			const double d = std::abs(sf.scatter_data.F_calc[r] - F_I[r]);
			dmax = std::max(dmax, d);
			fmax = std::max(fmax, std::abs(F_I[r]));
			num += std::abs(std::abs(sf.scatter_data.F_calc[r]) - std::abs(F_I[r]));
			den += std::abs(F_I[r]);
		}
		sf.scatter_data.F_calc = F_I;
		std::cout << "grid F check lambda " << lambda << ": max|dF|/max|F| " << std::scientific << std::setprecision(2) << dmax / fmax
			<< " sum||F_grid|-|F_I||/sum|F_I| " << num / den << std::fixed << " forward " << fwd_ms << " ms" << std::endl;
	}
	const std::string base = "NA2_" + lambda_tag(lambda);
	if (opt->write_tshb)
		write_tshb(sf_wave_vec[0], result, base + ".tshb");
	//NOSPHERA2_XCW_TSHB_TABLE: the table beside the .tshb, to compare the two
	if (!opt->write_tshb || tuning("NOSPHERA2_XCW_TSHB_TABLE"))
		result.write_tscb_file("test.cif", base + ".tscb");
	sf_wave_vec[0].write_wfn(base + ".wfn", false, true);
	if (!sf.ext_p_.empty())
		writer.log << "XCW lambda " << lambda << ": " << sf.extinction_report() << std::endl;
	if (opt->print_F)
		sf.write_F_calc(base + "_Fcalc.txt");
	//OCC's fchk is already in Gaussian's phases, which are libcint's (FormatConsistencyTests.OccFchkMatchesBridge);
	//flipping f(+-3), g(+-3), g(+-4) here wrote ORCA's and every fchk reader lost ~0.5 % of the electrons
	occ::qm::Wavefunction(wfn).save(base + ".fchk"); //save() is not const
	if (opt->xcw_settings.nbo_output)
		sf_wave_vec[0].write_nbo(base + ".47", opt->debug, &writer.log);
	//Neither file written above can carry this analysis - a .wfn has bare primitives and
	//the fchk reader keeps no shells - so -rgbi runs it here, on the refined wavefunction,
	//and its report goes to a file of its own per lambda
	if (opt->rgbi) {
		std::ofstream rgbi_out(base + "_RGBI.txt");
		std::streambuf* const cout_buf = std::cout.rdbuf(rgbi_out.rdbuf());
		Roby_information Roby(sf_wave_vec[0], opt->rgbi_group_sets, !opt->rgbi_no_sym,
			opt->rgbi_orbital_basis == RGBIOrbitalBasis::ANO, opt->rgbi_EVs, opt->rgbi_theta,
			opt->rgbi_legacy_cutoff);
		std::cout.rdbuf(cout_buf);
		writer.log << "RGBI analysis written to " << base << "_RGBI.txt" << std::endl;
	}
}

//The OCC molecule, the orbital basis (and the density fitting basis when one is set) and the SCF solver
occ::qm::HartreeFock XCW_wrapper::setup_system() {

	// Sets up the OCC molecule
	const occ::core::Molecule mol = SCF_solver::setup_SCF_mol(sf.asym_atoms, opt->xcw_settings.charge, opt->xcw_settings.multiplicity);
	// Sets up the OCC basis set
	const occ::qm::AOBasis occ_basis_set = sf.orbital_basis->to_AOBasis(mol.atoms());
	occ::qm::HartreeFock hf(occ_basis_set);

	// Handles density fitting basis
	sf.df_reserved_ = 0;
	if (!opt->xcw_settings.df_basis_name.empty()) density_fitting(hf, opt->xcw_settings.df_basis_name, "of the Fock build", sf.free_memory());

	if (opt->xcw_int_precision > 0.0) hf.set_precision(opt->xcw_int_precision);

	// The SCF gets its slice of the settings and no more
	options::SCF_settings scf_settings = opt->xcw_settings;
	scf_settings.incremental = opt->xcw_incremental;
	scf_settings.gpu_eri = opt->gpu_itensor && opt->use_gpu;
	scf_settings.no_date = opt->no_date;
	scf_solver = XCW_solver(mol, occ_basis_set, scf_settings, writer, *this);
	return hf;
	// closing function
}

void XCW_wrapper::density_fitting(occ::qm::HartreeFock& target, const std::string& name, const char* what, const size_t avail) {
	const occ::qm::AOBasis& basis = target.aobasis();
	//a copy: the library's set is shared, and AutoAux adds the elements it lacks (U in def2-universal-jkfit)
	std::shared_ptr<BasisSet> aux = std::make_shared<BasisSet>(*BasisSetLibrary::get_basis_set(name));
	if (!std::all_of(basis.atoms().begin(), basis.atoms().end(), [&](const occ::core::Atom& a) { return aux->has_element(a.atomic_number); })) {
		//gen_auto_aux reads the WFN-file convention, types from 1; the library stores l
		WFN atoms_wfn(e_origin::NOT_YET_DEFINED);
		for (const auto& a : basis.atoms()) {
			const int i = atoms_wfn.get_ncen();
			atoms_wfn.push_back_atom(constants::atnr2letter(a.atomic_number), a.x, a.y, a.z, a.atomic_number);
			for (const SimplePrimitive& p : (*sf.orbital_basis)[a.atomic_number - 1])
				atoms_wfn.push_back_atom_basis_set(i, p.exp, p.coefficient, p.type + 1, p.shell);
		}
		gen_missing_basis_auto_aux(atoms_wfn, aux);
	}
	if (aux->get_name().find("jkfit") == std::string::npos)
		std::cout << "WARNING: " << aux->get_name() << " is not a JK-fitting basis; Hartree-Fock exchange is fitted with it all the same. "
			"def2-universal-jkfit serves the def2 family, cc-pvXz-jkfit the cc-pVXZ family." << std::endl;
	// Create a list storing all elements of the molecule
	ivec elements;
	for (const auto& a : basis.atoms())
		if (std::find(elements.begin(), elements.end(), a.atomic_number) == elements.end())
			elements.push_back(a.atomic_number);
	const std::string file = aux->get_name() + "_df.json";
	aux->write_occ_json(file, elements);
	target.set_density_fitting_basis(file);
	//OCC keeps the three-index integrals only under a 512 MB limit and otherwise recomputes
	//them every iteration, which costs twice a direct build here. Held whenever they fit in
	//half of what the process can have; they are computed once for the whole lambda scan.
	const size_t naux = aux->to_AOBasis(basis.atoms()).nbf();
	const size_t nbf = basis.nbf();
	//OCC holds the full square nbf x nbf per auxiliary function, not the triangle
	const size_t store = naux * nbf * nbf * sizeof(double);
	const bool stored = avail == 0 || store < avail / 2;
	if (stored) sf.df_reserved_ += store;
	target.set_density_fitting_policy(stored ? occ::qm::IntegralEngineDF::Policy::Stored : occ::qm::IntegralEngineDF::Policy::Direct);
	if (opt->xcw_int_precision > 0.0) target.set_precision(opt->xcw_int_precision);
	std::cout << "XCW density fitting " << what << " with " << aux->get_name() << ": " << naux << " functions, "
		<< (store / 1048576.0) << " MB of three-index integrals " << (stored ? "held in memory" : "recomputed every build") << std::endl;
}

//Everything a converged step puts on screen, in the order it always has: the halting statistic
//first so its A^2 can be the last column of the row
void XCW_wrapper::report_lambda(const double lambda, const occ::qm::Wavefunction& wfn, const bool write_result) {
	//lambda 0 does not feel F: density-following shares are fit once, to the converged density, and the row is of them
	if (lambda == 0.0 && sf.XG.ws) {
		weight_D_fit_ = spin_summed(wfn.mo.D);
		sf.XG.refresh = true;
		sf.calc_F_calc(weight_D_fit_);
		sf.XG.refresh = false;
		sf.eval_scale();
		sf.calc_criteria();
	}
	if (opt->xcw_settings.xcw_gaussian_halt)
		evaluate_gaussian_halting(lambda);
	const fit_quality fit = cryst_crit();
	writer.v.criterion = fit.criterion;
	writer.v.GooF2 = fit.GooF2;
	writer.v.R1 = fit.R1;
	writer.v.criterion_all = fit.criterion_all;
	writer.v.R1_all = fit.R1_all;
	writer.v.has_A2 = opt->xcw_settings.xcw_gaussian_halt && !gaussian_halt_history_.empty();
	if (writer.v.has_A2) writer.v.A2 = gaussian_halt_history_.back().A2;
	writer.lambda_line();
	if (write_result) {
		const _time_point tscb_t0 = get_time();
		write_step_results(wfn, lambda);
		throughput::record_time("XCW tscb", false, get_msec(tscb_t0, get_time()));
	}
	if (tuning("NOSPHERA2_XCW_WEIGHT_CHECK") && sf.XG.ws)
		weight_check(wfn, lambda);
	const std::string shares = sf.grid_weight_report();
	if (!shares.empty()) {
		std::ostringstream s;
		s << "XCW lambda " << std::fixed << std::setprecision(5) << lambda << ": " << shares;
		if (weight_every_ < 0 && lambda != 0.0) s << ", " << weight_rounds_ << " more solves";
		std::cout << s.str() << std::endl;
		writer.log << s.str() << std::endl;
	}
	weight_rounds_ = 0;
}

//F_calc, scale and criteria are left of the refit shares
bool XCW_wrapper::refit_moves(const occ::qm::Wavefunction& wfn, const double lambda) {
	const dMatrix2 D = spin_summed(wfn.mo.D);
	sf.calc_F_calc(D);
	const cvec F0 = sf.scatter_data.F_calc;
	weight_D_fit_ = D;
	sf.XG.refresh = true;
	sf.calc_F_calc(D);
	sf.XG.refresh = false;
	sf.eval_scale();
	sf.calc_criteria();
	double d = 0, m = 0;
	for (int r = 0; r < sf.model_data.nr; r++) {
		d = std::max(d, std::abs(sf.scatter_data.F_calc[r] - F0[r]));
		m = std::max(m, std::abs(F0[r]));
	}
	const char* ft = tuning("NOSPHERA2_XCW_WEIGHT_FTOL");
	const bool moved = d > (ft ? std::atof(ft) : 1e-5) * m;
	writer.log << "XCW lambda " << std::fixed << std::setprecision(5) << lambda << ": shares refit to the converged density, max |dF| "
		<< std::scientific << std::setprecision(3) << d << " of max |F| " << m << (moved ? ", solved again" : "") << std::defaultfloat << std::endl;
	weight_rounds_ += moved;
	return moved;
}

//F(D1, s1) against F(D1, s0), what the shares' move over the step does, and F(D1, s0) against F(D0, s0), what the density's does,
//D0 and s0 the step before's; the first is the derivative term the perturbation leaves out, taken over a whole step
void XCW_wrapper::weight_check(const occ::qm::Wavefunction& wfn, const double lambda) {
	const dMatrix2 D1 = spin_summed(wfn.mo.D);
	const vec2 w1 = sf.XG.w;
	if (weight_D_.extent(0) == 0) {
		weight_D_ = D1;
		weight_w_ = w1;
		return;
	}
	const cvec F_keep = sf.scatter_data.F_calc;
	const double scale_keep = sf.scatter_data.scale;
	const auto q_keep = sf.quality_criteria;
	double R1[3];
	cvec F[3];
	const dMatrix2* Ds[3] = { &D1, &D1, &weight_D_ };
	const vec2* ws[3] = { &w1, &weight_w_, &weight_w_ };
	for (int c = 0; c < 3; c++) {
		sf.XG.w = *ws[c];
		sf.grid_F(*Ds[c]);
		sf.solve_scale();
		sf.calc_criteria();
		F[c] = sf.scatter_data.F_calc;
		R1[c] = sf.quality_criteria.R1;
	}
	double dw = 0, dd = 0, fm = 0;
	for (int r = 0; r < sf.model_data.nr; r++) {
		dw = std::max(dw, std::abs(F[0][r] - F[1][r]));
		dd = std::max(dd, std::abs(F[1][r] - F[2][r]));
		fm = std::max(fm, std::abs(F[0][r]));
	}
	std::ostringstream s;
	s << "XCW share check lambda " << std::fixed << std::setprecision(5) << lambda << ": max |dF| over the step from the shares "
		<< std::scientific << std::setprecision(3) << dw << ", from the density " << dd << " (max |F| " << fm << "); R1(gt) "
		<< std::fixed << std::setprecision(5) << R1[0] << " with this step's shares, " << R1[1] << " with the last step's";
	std::cout << s.str() << std::endl;
	writer.log << s.str() << std::endl;
	sf.XG.w = w1;
	sf.scatter_data.F_calc = F_keep;
	sf.scatter_data.scale = scale_keep;
	sf.quality_criteria = q_keep;
	weight_D_ = D1;
	weight_w_ = w1;
}

//The structure factor side of the run: Debye-Waller and phase factors, the anomalous dispersion
//corrections, XCW's integration grid and the RI partition when asked for, on the device when asked to
void XCW_wrapper::setup_structure_factors(const occ::qm::HartreeFock& hf) {
	sf.set_DW();
	sf.set_phases();
	sf.set_anom();
	//the symmetry links of the reflections, the AO pair screen, the atom grids and the AO values on them in blocks
	sf.setup_grid(hf.aobasis());
	//(mu nu|P) and the factorised fit metric of the RI partition
	if (opt->xcw_partition == "ri")
		sf.setup_ri(hf.aobasis());
	const char* wu = tuning("NOSPHERA2_XCW_WEIGHT_UPDATE");
	weight_every_ = wu ? std::max(0, std::atoi(wu)) : -1;
	const char* wd = tuning("NOSPHERA2_XCW_WEIGHT_DTOL");
	weight_dtol_ = wd ? std::atof(wd) : 1e-3;
#ifdef NOSPHERA2_USE_GPU
	//GPU: the grid's two block GEMMs on the device, the transform between them stays on the host.
	//Single precision there, so not under i_double
	if (opt->gpu_itensor && opt->use_gpu && !opt->xcw_settings.grid_double) {
		const bool held = sf.hold_grid_gpu();
		if (!(opt->no_date))
			std::cerr << "GPU in use: XCW grid blocks on " << (held ? "the device" : "the CPU - device unavailable or the blocks too large") << std::endl;
	}
#endif
}

//The two-electron integrals of the SCF: stored when they fit in memory, and a density-fitted Fock
//build for the TRAH Hessian products
void XCW_wrapper::setup_scf_integrals(const occ::qm::HartreeFock& hf) {
	const size_t avail = sf.free_memory();
	//read (read_eri) or built when they fit what is left; NOSPHERA2_XCW_ERI_NO_BUDGET holds them whatever -mem leaves
	scf_solver.setup_eri(hf, tuning("NOSPHERA2_XCW_ERI_NO_BUDGET") ? 0 : avail);
	//Fitted Hessian products only change the TRAH steps, never the converged result. Without stored
	//integrals a product costs a direct build and two fitted builds are cheaper, so they are fitted
	//then by default; `hessian_df_basis none` keeps them exact. Set up after the stored integrals so
	//it never squeezes them out
	const std::string& hdf = opt->xcw_settings.hessian_df_basis_name;
	std::string hdf_name = hdf == "none" ? "" : hdf;
	//elements def2-universal-jkfit lacks get AutoAux in density_fitting
	if (hdf.empty() && !scf_solver.has_eri() && opt->xcw_settings.df_basis_name.empty())
		hdf_name = "def2-universal-jkfit";
	if (!hdf_name.empty()) {
		occ::qm::HartreeFock hf_df(hf.aobasis());
		density_fitting(hf_df, hdf_name, "of the Hessian products", scf_solver.has_eri() ? avail / 5 : avail);
		scf_solver.set_hessian_df(std::move(hf_df));
	}
}

//The orbitals a previous run wrote for the start lambda, when the settings ask for them
bool XCW_wrapper::read_start_guess(occ::qm::Wavefunction& wfn) {
	if (!opt->xcw_settings.read_first_guess) return false;
	wfn = occ::qm::Wavefunction::load("NA2_" + lambda_tag(opt->xcw_settings.xcw_start_value) + ".fchk");
	return true;
}

//One SCF at lambda from the orbitals in wfn (a fresh guess unless use_guess), reported when it
//converged; wfn holds the result
bool XCW_wrapper::run_lambda(occ::qm::HartreeFock& hf, const double lambda, occ::qm::Wavefunction& wfn, const bool use_guess, const bool write_result) {
	scf_solver.lambda = lambda;
	weight_count_ = 0;
	weight_hold_ = weight_moved_ = false;
	//DIIS iterations with damping, level shift and rescue, then TRAH steps where DIIS stalls
	bool converged = scf_solver.run(hf, wfn, use_guess);
	//shares held within a solve and refit to its density until that no longer moves F: the fixed point refitting at every
	//evaluation goes for, with every solve one objective for TRAH to accept steps on
	for (int round = 0; converged && weight_every_ < 0 && sf.XG.ws && lambda != 0.0 && round < 20 && refit_moves(wfn, lambda); round++)
		converged = (weight_hold_ = true) && scf_solver.run(hf, wfn, true);
	if (converged) report_lambda(lambda, wfn, write_result);
	return converged;
}

//After lambda failed from last_lambda: halve the step until one converges, then try lambda again
//from there, until lambda converges or the step falls below min_step. last_lambda and last_wfn
//follow every intermediate step that converged
bool XCW_wrapper::approach_lambda(occ::qm::HartreeFock& hf, const double lambda, const double min_step, double& last_lambda, occ::qm::Wavefunction& last_wfn, occ::qm::Wavefunction& wfn) {
	double lambda_step = lambda - last_lambda;
	bool converged = false;
	while (!converged && lambda_step > min_step) {
		lambda_step *= 0.5;
		const double trial_lambda = std::min(lambda, last_lambda + lambda_step);
		writer.log << "XCW: retrying lambda " << std::fixed << std::setprecision(8) << trial_lambda
			<< " from converged lambda " << last_lambda << " with step " << lambda_step << std::endl;
		std::cout << "XCW: retrying lambda " << std::fixed << std::setprecision(8) << trial_lambda
			<< " with step " << lambda_step << std::endl;
		occ::qm::Wavefunction trial = last_wfn;
		if (!run_lambda(hf, trial_lambda, trial, true, trial_lambda == lambda)) continue;
		if (trial_lambda == lambda) {
			wfn = std::move(trial);
			return true;
		}
		last_lambda = trial_lambda;
		last_wfn = trial;
		lambda_step = std::min(2.0 * lambda_step, lambda - last_lambda);
		wfn = last_wfn;
		converged = run_lambda(hf, lambda, wfn, true, true);
	}
	return converged;
}

//At the end of the planned scan with A^2 still falling, one more step, up to max_extra_steps
//beyond the requested range, so a trend that never turns at most doubles the scan
void XCW_wrapper::extend_scan(const double lambda, int& planned_steps, int& extra_steps, const int max_extra_steps) {
	double estimated_minimum = 0.0;
	if (!halting_minimum_beyond_scan(gaussian_halt_history_, estimated_minimum)) return;
	std::ostringstream msg;
	msg << std::fixed << std::setprecision(5)
		<< "XCW: A^2 is still falling at lambda " << lambda << ", so the minimum is outside the requested range - ";
	if (extra_steps < max_extra_steps) {
		planned_steps++;
		extra_steps++;
		msg << "extending the scan to lambda " << (lambda + opt->xcw_settings.xcw_step_size)
			<< " (extra step " << extra_steps << " of at most " << max_extra_steps << ")";
		if (estimated_minimum > 0.0)
			msg << ", extrapolated minimum near lambda ~= " << estimated_minimum;
	}
	else
		msg << "stopping anyway, the scan has already been extended by its limit of "
			<< max_extra_steps << " steps; raise max_value in -do_XCW to continue";
	writer.log << msg.str() << std::endl;
	std::cout << msg.str() << std::endl;
}

void XCW_wrapper::run() {
	//OCC parallelises through TBB, which does not read OMP_NUM_THREADS; this makes -cpus bind OCC's part of the run too
	occ::parallel::set_num_threads(opt->threads > 0 ? opt->threads : omp_get_max_threads());
	//the molecule, the orbital basis and the SCF solver
	occ::qm::HartreeFock hf = setup_system();
	//DW, phase and anomalous dispersion factors, the integration grid with its AO values and the RI partition
	setup_structure_factors(hf);
	//NOSPHERA2_XCW_RCHECK=<fchk>: scale, extinction and R1 of that converged wavefunction against this model and
	//data, its F_calc in rcheck_Fcalc.txt, then stop; no integrals, no SCF
	if (const char* rc = tuning("NOSPHERA2_XCW_RCHECK")) {
		const occ::qm::Wavefunction w = occ::qm::Wavefunction::load(rc);
		weight_D_fit_ = spin_summed(w.mo.D);
		sf.XG.refresh = true;
		sf.calc_F_calc(weight_D_fit_);
		sf.XG.refresh = false;
		//the SCF calls eval_scale every iteration, so scale and extinction have settled by convergence
		for (int i = 0; i < 10; i++) sf.eval_scale();
		sf.calc_criteria();
		std::ostringstream s;
		s << "XCW RCHECK " << rc << ": R1 " << sf.quality_criteria.R1 << " R1(all) " << sf.quality_criteria.R1_all
			<< " scale " << sf.scatter_data.scale << " " << sf.extinction_report();
		std::cout << s.str() << std::endl;
		writer.log << s.str() << std::endl;
		sf.write_F_calc("rcheck_Fcalc.txt");
		return;
	}
	//the stored two-electron integrals and the density-fitted Hessian products of the SCF
	setup_scf_integrals(hf);
	occ::qm::Wavefunction last_wfn, prev_wfn;
	bool has_guess = read_start_guess(last_wfn);
	const occ::Mat S_ao = hf.compute_overlap_matrix();
	std::cout << "More detailed output in SCF.log file..." << std::endl;
	if (opt->xcw_settings.XWR_type == 2) {
		std::cout << "XCW: fitting against the 1/|H|^2-weighted residual self-energy criterion "
			<< "Criterion below is the square of this weighted quantity, not the classical GooF^2." << std::endl;
		writer.log << "XCW: fitting against the 1/|H|^2-weighted residual self-energy criterion "
			<< "Criterion below is the square of this weighted quantity, not the classical GooF^2." << std::endl;
	}
	XCW_log_writer::console_header(opt->xcw_settings.xcw_gaussian_halt);
	//The lambda scan: each step starts from the orbitals of the last converged one, extrapolated
	//along the path of the two before when asked to (-xcw_extrapolate)
	const double min_lambda_step = opt->xcw_settings.xcw_step_size / 128.0;
	double last_lambda = opt->xcw_settings.xcw_start_value, prev_lambda = last_lambda;
	bool have_prev = false;
	int planned_steps = opt->xcw_settings.num_xcw_steps, extra_steps = 0;
	const int max_extra_steps = opt->xcw_settings.num_xcw_steps;
	for (int step = 0; step < planned_steps; step++) {
		const double lambda = step * opt->xcw_settings.xcw_step_size + opt->xcw_settings.xcw_start_value;
		occ::qm::Wavefunction wfn = last_wfn;
		if (opt->xcw_extrapolate && have_prev)
			extrapolate_orbitals(wfn.mo, prev_wfn.mo, S_ao, (lambda - last_lambda) / (last_lambda - prev_lambda));
		//one SCF at this lambda and, when it converged, its row and its files
		bool converged = run_lambda(hf, lambda, wfn, has_guess, true);
		//after the first step, so that density-following shares are the scheme's
		if (step == 0 && converged && tuning("NOSPHERA2_XCW_ADJOINT_CHECK"))
			sf.check_adjoint();
		if (!converged && step > 0)
			converged = approach_lambda(hf, lambda, min_lambda_step, last_lambda, last_wfn, wfn);
		if (!converged) {
			//step 0 has no converged neighbour to continue from, so it is not retried
			std::ostringstream why;
			why << "XCW: unable to converge lambda " << std::fixed << std::setprecision(8) << lambda;
			if (step == 0)
				why << " in " << opt->xcw_settings.max_scf_iterations << " SCF iterations (raise max_iter or loosen the criteria); stopping scan.";
			else
				why << " with a continuation step above " << min_lambda_step << "; stopping scan.";
			writer.log << why.str() << std::endl;
			std::cout << why.str() << std::endl;
			break;
		}
		//the jump off the unconstrained lambda 0 is far from linear, so it is never extrapolated
		have_prev = has_guess && lambda > last_lambda && last_lambda != 0.0;
		prev_wfn = std::move(last_wfn);
		prev_lambda = last_lambda;
		last_wfn = std::move(wfn);
		last_lambda = lambda;
		has_guess = true;
		if (opt->xcw_settings.xcw_gaussian_halt && step + 1 == planned_steps)
			extend_scan(lambda, planned_steps, extra_steps, max_extra_steps);
		//Progress estimate every 5 lambda steps; not at the last one, the summary below prints a final one
		if (opt->xcw_settings.xcw_gaussian_halt && (step + 1) % 5 == 0 && step + 1 < planned_steps)
			report_halting_progress_estimate(false);
	}
	if (opt->xcw_settings.xcw_gaussian_halt)
		report_gaussian_halting_summary();
#ifdef NOSPHERA2_USE_GPU
	//GPU: gives the grid blocks back
	xcw_grid_gpu_release();
	sf.XG.on_device = false;
#endif
	scf_solver.release_device();
	std::cout << "Finished XCW fitting procedure." << std::endl;
}
