#include "pch.h"
#include "ccsd_density.h"
#include "convenience.h"

#include <occ/core/util.h>
#include <occ/driver/method_parser.h>
#include <occ/main/occ_scf.h>
#include <toml++/toml.hpp>

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace
{
	//occ throws on bad input or non-convergence; NoSpherA2 ends the run through err_checkf
	template <typename F>
	CcsdDensityResult guarded(F &&f, std::ostream &log)
	{
		try
		{
			return f();
		}
		catch (const std::exception &e)
		{
			err_checkf(false, e.what(), log);
		}
		return {};
	}

	std::string integrals_label(const CcsdDensityResult &r)
	{
		return r.backend == "exact" ? std::string("exact integrals") : "DF integrals, " + r.aux_basis;
	}

	void occupations_summary(const CcsdDensityResult &r, std::ostream &log)
	{
		const int nmo = static_cast<int>(r.occupations.size());
		for (int k = 0; k < std::min(nmo, 8); k++) log << " " << std::setprecision(6) << r.occupations(k);
		log << (nmo > 8 ? " ..." : "") << "\n  sum of occupations: " << std::setprecision(8) << r.occupations.sum()
			<< ", smallest " << std::scientific << std::setprecision(3) << r.occupations(nmo - 1) << std::defaultfloat << std::endl;
	}
}

CcsdDensityResult ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log)
{
	const CcsdDensityResult r = guarded([&] { return occ::qm::cc::ccsd_density(scf, opts); }, log);
	log << "\n" << (opts.triples ? "CCSD(T)" : "CCSD") << " density (unrelaxed, Lambda equations), " << integrals_label(r) << "\n"
		<< "  correlated: " << r.n_occ << " occupied, " << r.n_vir << " virtual; frozen core: " << r.n_frozen << " orbitals\n"
		<< std::fixed << std::setprecision(10) << "  E(CCSD correlation) = " << r.e_corr << "  (" << r.cc_iterations << " iterations)\n";
	if (opts.triples) log << "  E(T) = " << r.e_triples << "\n";
	log << "  Lambda equations converged in " << r.lambda_iterations << " iterations\n"
		<< std::setprecision(6) << "  |t1| = " << r.t1_norm << "  |t2| = " << r.t2_norm << "  |l1| = " << r.l1_norm << "  |l2| = " << r.l2_norm << "\n"
		<< "  natural occupations:";
	occupations_summary(r, log);
	return r;
}

CcsdDensityResult mp2_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log)
{
	const CcsdDensityResult r = guarded([&] { return occ::qm::cc::mp2_relaxed_density(scf, opts); }, log);
	log << "\nMP2 density (relaxed, Z-vector), " << integrals_label(r) << "\n"
		<< "  correlated: " << r.n_occ << " occupied, " << r.n_vir << " virtual; frozen core: " << r.n_frozen << " orbitals\n"
		<< std::fixed << std::setprecision(10) << "  E(MP2 correlation) = " << r.e_corr << "\n"
		<< "  Z-vector equations converged in " << r.lambda_iterations << " iterations\n" << "  natural occupations:";
	occupations_summary(r, log);
	return r;
}

CcsdDensityResult local_ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log)
{
	const CcsdDensityResult r = guarded([&] { return occ::qm::cc::local_ccsd_density(scf, opts); }, log);
	log << "\n" << (opts.triples ? "local CCSD(T)" : "local CCSD") << " density (cluster-in-molecule, unrelaxed Lambda densities), "
		<< (r.backend == "df" ? "DF integrals, " + r.aux_basis : std::string("exact integrals")) << "\n"
		<< "  " << r.n_lmo << " Pipek-Mezey LMOs, frozen core: " << r.n_frozen << " orbitals; candidate radius " << opts.cluster_radius
		<< " A, domain threshold " << opts.domain_threshold << ", pair threshold " << opts.pair_threshold << " Eh, LNO threshold " << opts.lno_threshold << "\n"
		<< "  candidates " << std::fixed << std::setprecision(1) << r.mean_candidates << " LMOs, PAO domains " << r.mean_domain << " virtuals on average\n"
		<< "  " << r.n_clusters << " clusters: occupied " << r.mean_cluster_occ << " on average (max " << r.max_cluster_occ << "), virtual "
		<< r.mean_cluster_vir << " (max " << r.max_cluster_vir << ")\n"
		<< std::setprecision(10) << "  E(CCSD correlation, sum over clusters) = " << r.e_corr << "\n"
		<< "  natural occupations:";
	occupations_summary(r, log);
	return r;
}

CcsdDensityResult ccsd_density_from_occ_input(const std::string &toml_path, occ::io::OccInput config, std::ostream &log)
{
	const auto spec = occ::driver::parse_method_string(config.method.name);
	const bool mp2 = spec.kind == occ::driver::MethodKind::MP2;
	err_checkf(mp2 || spec.kind == occ::driver::MethodKind::CCSD || spec.kind == occ::driver::MethodKind::CCSD_T,
		"Method '" + config.method.name + "' has no density here; MP2, CCSD and CCSD(T) do (method = \"mp2\", \"ccsd\" or \"ccsd(t)\")", log);
	err_checkf(!mp2 || occ::util::to_lower_copy(config.method.mp2_spin_scaling) == "none",
		"SCS/SOS-MP2 has no density here (mp2_spin_scaling = \"" + config.method.mp2_spin_scaling + "\"); the plain MP2 density is the relaxed one", log);
	const auto input = toml::parse_file(toml_path);  //the view below points into it, so it must outlive the reads
	const auto scf_table = input["scf"];
	CcsdDensityOptions opts;
	opts.triples = spec.kind == occ::driver::MethodKind::CCSD_T;
	if (mp2 && config.method.mp2_backend != "auto") opts.backend = config.method.mp2_backend;
	if (!spec.backend.empty()) opts.backend = spec.backend;  //"df-ccsd", as occ reads it
	opts.n_frozen = scf_table["ccsd_frozen_core"].value_or(opts.n_frozen);
	opts.backend = scf_table["ccsd_backend"].value_or(opts.backend);
	opts.max_memory_gb = scf_table["ccsd_max_memory_gb"].value_or(opts.max_memory_gb);
	opts.tol = scf_table["ccsd_tol"].value_or(opts.tol);
	opts.aux_basis = config.basis.ri_basis;
	const bool local = scf_table["ccsd_local"].value_or(false);
	opts.cluster_radius = scf_table["ccsd_cluster_radius"].value_or(opts.cluster_radius);
	opts.domain_threshold = scf_table["ccsd_domain_threshold"].value_or(opts.domain_threshold);
	opts.pair_threshold = scf_table["ccsd_pair_threshold"].value_or(opts.pair_threshold);
	opts.lno_threshold = scf_table["ccsd_lno_threshold"].value_or(opts.lno_threshold);
	config.method.name = "rhf";
	//exact CCSD on an exact-J/K reference, as ORCA's NoRI
	const std::string backend = occ::util::to_lower_copy(opts.backend);
	if (backend != "df" && backend != "ri") config.method.ri_policy = occ::io::RIPolicy::None;
	err_checkf(!mp2 || !local, "ccsd_local applies to CCSD and CCSD(T), not to MP2", log);
	const occ::qm::Wavefunction scf = occ::main::run_scf_external(config, false);
	if (mp2) return mp2_natural_orbitals(scf, opts, log);
	return local ? local_ccsd_natural_orbitals(scf, opts, log) : ccsd_natural_orbitals(scf, opts, log);
}

occ::qm::Wavefunction wavefunction_from_occ_input(const std::string &toml_path, const occ::io::OccInput &config, std::ostream &log)
{
	const auto kind = occ::driver::parse_method_string(config.method.name).kind;
	if (kind != occ::driver::MethodKind::MP2 && kind != occ::driver::MethodKind::CCSD && kind != occ::driver::MethodKind::CCSD_T)
		return occ::main::run_scf_external(config, true);
	return ccsd_density_from_occ_input(toml_path, config, log).natural;
}
