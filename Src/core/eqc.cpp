#include "pch.h"
#include "eqc.h"
#include "convenience.h"
#include "wfn_class.h"
#include "section_log.h"
#include "citations.h"
#include <occ/dft/dft.h>
#include <occ/qm/hf.h>

namespace eqc
{
	double covalency(const double x, const double v)
	{
		if (x < 0.0) return 100.0 * std::abs(x) / (std::abs(v) + std::abs(x));
		if (x > 0.0) return 100.0 * v / (v - x);
		return std::numeric_limits<double>::quiet_NaN();
	}

	reaction react(const std::vector<terms> &reactants, const terms &product)
	{
		double E = 0.0, nX = 0.0, Vnn = 0.0, Eee = 0.0;
		for (const terms &t : reactants)
		{
			E += t.E;
			nX += t.nX;
			Vnn += t.Vnn;
			Eee += t.Eee();
		}
		reaction r;
		r.dE = product.E - E;
		r.dnX = product.nX - nX;
		r.dVnn = product.Vnn - Vnn;
		r.dEee = product.Eee() - Eee;
		r.Q = 2.0 * r.dnX / r.dE - 1.0;
		r.covalency = covalency(r.dnX, r.dV());
		r.dEeeE = 100.0 * (std::abs(product.Eee() / product.E) - std::abs(Eee / E));
		r.n = product.n;
		return r;
	}

	double nuclear_repulsion(const ivec &Z, const std::vector<d3> &pos)
	{
		double V = 0.0;
		for (size_t i = 0; i < Z.size(); i++)
			for (size_t j = 0; j < i; j++)
			{
				const double dx = pos[i][0] - pos[j][0], dy = pos[i][1] - pos[j][1], dz = pos[i][2] - pos[j][2];
				V += Z[i] * Z[j] / std::sqrt(dx * dx + dy * dy + dz * dz);
			}
		return V;
	}

	double orca_energy(const std::string &out_file)
	{
		std::ifstream f(out_file);
		double E = std::numeric_limits<double>::quiet_NaN();
		const std::string key = "FINAL SINGLE POINT ENERGY";
		for (std::string line; std::getline(f, line);)
		{
			const size_t p = line.find(key);
			if (p != std::string::npos)
				E = std::stod(line.substr(p + key.size()));
		}
		return E;
	}

	namespace
	{
		using occ::Mat;
		using occ::Vec;
		using occ::qm::MolecularOrbitals;
		using occ::qm::SpinorbitalKind;

		//One wavefunction of the analysis: its terms and where they came from
		struct entry
		{
			std::string name, program;
			terms t;
			int iterations = -1, cold_iterations = -1;  //-1: no SCF ran
			bool converged = true;
		};

		ivec atomic_numbers(const std::vector<occ::core::Atom> &atoms)
		{
			ivec Z;
			for (const auto &a : atoms) Z.push_back(a.atomic_number);
			return Z;
		}

		std::vector<d3> positions(const std::vector<occ::core::Atom> &atoms)
		{
			std::vector<d3> p;
			for (const auto &a : atoms) p.push_back({ a.x, a.y, a.z });
			return p;
		}

		//Twice the orbital energy sum over the doubly occupied orbitals, or alpha plus beta
		double orbital_energy_sum(const MolecularOrbitals &mo)
		{
			double s = 0.0;
			if (mo.kind == SpinorbitalKind::Restricted)
				for (size_t i = 0; i < mo.n_alpha; i++) s += 2.0 * mo.energies(i);
			else
			{
				for (size_t i = 0; i < mo.n_alpha; i++) s += mo.energies(i);
				for (size_t i = 0; i < mo.n_beta; i++) s += mo.energies(mo.n_ao + i);
			}
			return s;
		}

		//Half the spin-summed density, the one matrix a closed-shell Fock build takes
		Mat half_density(const MolecularOrbitals &mo)
		{
			const size_t nbf = mo.n_ao;
			if (mo.kind == SpinorbitalKind::Restricted)
			{
				const auto Co = mo.C.leftCols(mo.n_alpha);
				return Co * Co.transpose();
			}
			const auto Ca = mo.C.topRows(nbf).leftCols(mo.n_alpha);
			const auto Cb = mo.C.bottomRows(nbf).leftCols(mo.n_beta);
			return 0.5 * (Ca * Ca.transpose() + Cb * Cb.transpose());
		}

		terms occ_terms(const occ::qm::Wavefunction &w, const int charge)
		{
			terms t;
			t.E = w.energy.total;
			t.nX = orbital_energy_sum(w.mo);
			t.Vnn = nuclear_repulsion(atomic_numbers(w.atoms), positions(w.atoms));
			for (const auto &a : w.atoms) t.n += a.atomic_number;
			t.n -= charge;
			return t;
		}

		struct scf_result
		{
			occ::qm::Wavefunction wfn;
			int iterations = 0;
			bool converged = false;
		};

		template <class P>
		void tighten(occ::qm::SCF<P> &scf)
		{
			//Orbital energies are linear in the density error, E only quadratic: nX needs the tight commutator
			scf.convergence_settings.energy_threshold = 1e-10;
			scf.convergence_settings.commutator_threshold = 1e-9;
		}

		//occ's add_guess_density from a given density: H + G[D_half], diagonalised
		template <class P>
		void seed_density(occ::qm::SCF<P> &scf, P &proc, const Mat &D_half)
		{
			scf.update_occupied_orbital_count();
			scf.set_core_matrices();
			scf.ctx.F = scf.ctx.H;
			scf.set_conditioning_orthogonalizer();
			scf.ctx.K = proc.compute_schwarz_ints();
			MolecularOrbitals g;
			g.kind = SpinorbitalKind::Restricted;
			g.n_ao = D_half.rows();
			g.n_alpha = scf.n_alpha();
			g.n_beta = scf.n_beta();
			g.D = D_half;
			const Mat G = proc.compute_fock_from_density(g, scf.ctx.K);
			if (scf.ctx.mo.kind == SpinorbitalKind::Unrestricted)
			{
				occ::qm::block::a(scf.ctx.F) += G;
				occ::qm::block::b(scf.ctx.F) += G;
			}
			else
				scf.ctx.F += G;
			scf.ctx.orthogonalizer.orthogonalize_molecular_orbitals(scf.ctx.mo, scf.ctx.F);
			scf.m_have_initial_guess = true;
		}

		//Seeds from orbitals, a density block or occ's own guess. e_first is iteration 1, i.e. E of the seed itself.
		template <class P>
		scf_result converge(P &proc, const int charge, const int mult, const SpinorbitalKind kind,
			const MolecularOrbitals *orbitals, const Mat *D_half, double *e_first = nullptr)
		{
			const auto seeded = [&](occ::qm::SCF<P> &scf) {
				tighten(scf);
				scf.set_charge_multiplicity(charge, mult);
				if (orbitals)
				{
					occ::qm::Wavefunction w;
					w.mo = *orbitals;
					scf.set_initial_guess_from_wfn(w);
				}
				else if (D_half)
					seed_density(scf, proc, *D_half);
			};
			if (e_first)
			{
				occ::qm::SCF<P> one(proc, kind);
				seeded(one);
				one.maxiter = 1;
				*e_first = one.compute_scf_energy();
			}
			occ::qm::SCF<P> scf(proc, kind);
			seeded(scf);
			scf.compute_scf_energy();
			return { scf.wavefunction(), scf.iter, scf.ctx.converged };
		}

		//The gbw's basis. ORCA stores normalised contractions, occ wants them bare and normalises itself.
		occ::gto::AOBasis gbw_basis(const WFN &w, const std::vector<occ::core::Atom> &atoms)
		{
			std::vector<occ::gto::Shell> shells;
			for (int a = 0; a < w.get_ncen(); a++)
			{
				const atom at = w.get_atom(a);
				const std::vector<basis_set_entry> basis = at.get_basis_set();
				for (size_t s = 0; s < basis.size();)
				{
					const int l = static_cast<int>(basis[s].get_type()) - 1;
					vec expo, con;
					const unsigned int shell = basis[s].get_shell();
					for (; s < basis.size() && basis[s].get_shell() == shell; s++)
					{
						expo.push_back(basis[s].get_exponent());
						con.push_back(basis[s].get_coefficient() / occ::gto::gto_norm(l, expo.back()));
					}
					occ::gto::Shell sh(l, expo, { con }, { atoms[a].x, atoms[a].y, atoms[a].z });
					sh.incorporate_shell_norm();
					shells.push_back(sh);
				}
			}
			occ::gto::AOBasis b(atoms, shells, "gbw");
			b.set_pure(true);
			return b;
		}

		//Same shells, order and exponents (1e-6 relative): the gbw coefficients then belong to this basis
		std::string basis_mismatch(const occ::gto::AOBasis &a, const occ::gto::AOBasis &b)
		{
			if (a.size() != b.size())
				return std::to_string(a.size()) + " shells against " + std::to_string(b.size());
			for (size_t s = 0; s < a.size(); s++)
			{
				const auto &x = a[s], &y = b[s];
				bool same = x.l == y.l && x.num_primitives() == y.num_primitives()
					&& a.shell_to_atom()[s] == b.shell_to_atom()[s];
				for (size_t p = 0; same && p < x.num_primitives(); p++)
					same = std::abs(x.exponents(p) - y.exponents(p)) <= 1e-6 * x.exponents(p);
				if (!same) return "shell " + std::to_string(s + 1) + " differs";
			}
			return "";
		}

		//The reader already orders m = -l..l but keeps ORCA's phases (orca_pure_sign_flips)
		MolecularOrbitals gbw_orbitals(const WFN &w, const occ::gto::AOBasis &basis)
		{
			const int nbf = static_cast<int>(basis.nbf());
			const int ops = w.get_is_unrestricted() ? 2 : 1;
			err_checkf(static_cast<int>(w.get_MO_sph().extent(0)) == ops * nbf && static_cast<int>(w.get_MO_sph().extent(1)) == nbf,
				"-eqc needs the gbw's spherical coefficients for " + std::to_string(nbf) + " basis functions", std::cout);
			MolecularOrbitals mo;
			mo.kind = ops == 2 ? SpinorbitalKind::Unrestricted : SpinorbitalKind::Restricted;
			mo.n_ao = nbf;
			mo.C = Mat(ops * nbf, nbf);
			mo.energies = Vec(ops * nbf);
			for (int r = 0; r < ops * nbf; r++)
				for (int c = 0; c < nbf; c++)
					mo.C(r, c) = w.get_MO_sph()(r, c);
			size_t n_occ[2] = { 0, 0 };
			for (int s = 0; s < ops; s++)
				for (int j = 0; j < nbf; j++)
				{
					const int k = s * nbf + j;
					mo.energies(k) = w.get_MO_energy(k);
					const double o = w.get_MO_occ(k);
					err_checkf(ops == 2 || o < 0.5 || o > 1.5, "-eqc: the gbw is restricted open-shell, which occ does not converge; run UHF/UKS", std::cout);
					if (o > 0.5)
					{
						err_checkf(static_cast<size_t>(j) == n_occ[s], "-eqc: the gbw's occupied orbitals are not the lowest ones", std::cout);
						n_occ[s]++;
					}
				}
			mo.n_alpha = n_occ[0];
			mo.n_beta = ops == 2 ? n_occ[1] : n_occ[0];
			const auto &first = basis.first_bf();
			for (size_t s = 0; s < basis.size(); s++)
			{
				const int l = basis[s].l;
				for (int m = 3; m <= l; m++)
					if (orca_pure_sign_flips(m))
						for (int spin = 0; spin < ops; spin++)
						{
							mo.C.row(spin * nbf + first[s] + l + m) *= -1.0;
							mo.C.row(spin * nbf + first[s] + l - m) *= -1.0;
						}
			}
			return mo;
		}

		//max |C_occ^T S C_occ - 1| over both spins: zero when the coefficients belong to this basis
		double orthonormality_error(const MolecularOrbitals &mo, const Mat &S)
		{
			const size_t nbf = mo.n_ao;
			double err = 0.0;
			const auto check = [&](const Mat &C) {
				const Mat M = C.transpose() * S * C - Mat::Identity(C.cols(), C.cols());
				err = std::max(err, M.cwiseAbs().maxCoeff());
			};
			if (mo.kind == SpinorbitalKind::Restricted)
				check(mo.C.leftCols(mo.n_alpha));
			else
			{
				check(mo.C.topRows(nbf).leftCols(mo.n_alpha));
				if (mo.n_beta > 0) check(mo.C.bottomRows(nbf).leftCols(mo.n_beta));
			}
			return err;
		}

		//occ keeps the larger spin count in beta, ORCA in alpha
		void match_spins(MolecularOrbitals &mo, const size_t n_alpha)
		{
			if (mo.kind != SpinorbitalKind::Unrestricted || mo.n_alpha == n_alpha) return;
			const size_t nbf = mo.n_ao;
			Mat a = mo.C.topRows(nbf);
			mo.C.topRows(nbf) = mo.C.bottomRows(nbf);
			mo.C.bottomRows(nbf) = a;
			Vec e = mo.energies.head(nbf);
			mo.energies.head(nbf) = mo.energies.tail(nbf);
			mo.energies.tail(nbf) = e;
			std::swap(mo.n_alpha, mo.n_beta);
		}

		struct fragment_basis
		{
			std::vector<occ::core::Atom> atoms;
			occ::gto::AOBasis basis;
			ivec ao;  //parent AO index of every fragment AO
		};

		fragment_basis cut(const occ::gto::AOBasis &parent, ivec atoms)
		{
			std::sort(atoms.begin(), atoms.end());
			fragment_basis f;
			std::vector<occ::gto::Shell> shells, ecp;
			ivec ecp_electrons;
			for (const int a : atoms)
			{
				f.atoms.push_back(parent.atoms()[a]);
				ecp_electrons.push_back(parent.ecp_electrons()[a]);
			}
			for (size_t s = 0; s < parent.size(); s++)
				if (std::binary_search(atoms.begin(), atoms.end(), parent.shell_to_atom()[s]))
				{
					shells.push_back(parent[s]);
					for (size_t k = 0; k < parent[s].size(); k++) f.ao.push_back(static_cast<int>(parent.first_bf()[s] + k));
				}
			for (size_t s = 0; s < parent.ecp_shells().size(); s++)
				if (std::binary_search(atoms.begin(), atoms.end(), parent.ecp_shell_to_atom()[s]))
					ecp.push_back(parent.ecp_shells()[s]);
			f.basis = occ::gto::AOBasis(f.atoms, shells, parent.name(), ecp);
			f.basis.set_ecp_electrons(ecp_electrons);
			f.basis.set_pure(true);
			return f;
		}

		struct mode1_out
		{
			std::vector<entry> rows;  //parent first
			double e_orca = 0.0, e_seed = 0.0, ortho = 0.0;
			double nx_orca = 0.0;  //nX from the gbw's own orbital energies, the mixed-program check
		};

		template <class P>
		std::unique_ptr<P> procedure(const std::string &method, const occ::gto::AOBasis &basis);
		template <>
		std::unique_ptr<occ::qm::HartreeFock> procedure(const std::string &, const occ::gto::AOBasis &basis)
		{
			return std::make_unique<occ::qm::HartreeFock>(basis);
		}
		template <>
		std::unique_ptr<occ::dft::DFT> procedure(const std::string &method, const occ::gto::AOBasis &basis)
		{
			return std::make_unique<occ::dft::DFT>(method, basis);
		}

		template <class P>
		mode1_out mode1(const options &opt, const WFN &w, const occ::gto::AOBasis &basis, const MolecularOrbitals &gbw_mo,
			const int charge, const int mult)
		{
			mode1_out out;
			const std::string label = opt.eqc_method == "hf" ? "HF" : opt.eqc_method;
			auto proc = procedure<P>(opt.eqc_method, basis);
			out.ortho = orthonormality_error(gbw_mo, proc->compute_overlap_matrix());
			const bool usable = out.ortho < 1e-6;
			if (!usable)
				std::cout << "The gbw orbitals are not orthonormal in the rebuilt basis (max |C^T S C - 1| = "
				<< std::scientific << std::setprecision(2) << out.ortho << std::fixed
				<< "); the parent starts from occ's own guess instead" << std::endl;
			const SpinorbitalKind kind = gbw_mo.kind;
			out.nx_orca = orbital_energy_sum(gbw_mo);
			MolecularOrbitals seed = gbw_mo;
			{
				occ::qm::SCF<P> probe(*proc, kind);
				probe.set_charge_multiplicity(charge, mult);
				match_spins(seed, probe.n_alpha());
			}
			seed.update_occupied_orbitals();
			seed.update_density_matrix();
			const scf_result parent = converge(*proc, charge, mult, kind, usable ? &seed : nullptr, nullptr,
				usable ? &out.e_seed : nullptr);
			if (!usable) out.e_seed = std::numeric_limits<double>::quiet_NaN();
			out.e_orca = orca_energy((w.get_path().parent_path() / (w.get_path().stem().string() + ".out")).string());
			if (std::isnan(out.e_orca))
				out.e_orca = orca_energy((w.get_path().parent_path() / (w.get_path().stem().string() + ".log")).string());
			entry p{ "parent", "ORCA gbw, re-converged in occ " + label, occ_terms(parent.wfn, charge), parent.iterations, -1, parent.converged };
			if (opt.eqc_cold)
				p.cold_iterations = converge(*proc, charge, mult, kind, nullptr, nullptr).iterations;
			out.rows.push_back(p);
			const Mat D_parent = half_density(parent.wfn.mo);
			for (size_t i = 0; i < opt.eqc_frags.size(); i++)
			{
				const auto &fr = opt.eqc_frags[i];
				const fragment_basis f = cut(basis, fr.atoms);
				Mat D(f.ao.size(), f.ao.size());
				for (size_t r = 0; r < f.ao.size(); r++)
					for (size_t c = 0; c < f.ao.size(); c++)
						D(r, c) = D_parent(f.ao[r], f.ao[c]);
				auto fproc = procedure<P>(opt.eqc_method, f.basis);
				const SpinorbitalKind fkind = fr.mult > 1 ? SpinorbitalKind::Unrestricted : SpinorbitalKind::Restricted;
				const scf_result res = converge(*fproc, fr.charge, fr.mult, fkind, nullptr, &D);
				entry e{ "fragment " + std::to_string(i + 1), "occ " + label + ", seeded from the parent", occ_terms(res.wfn, fr.charge),
					res.iterations, -1, res.converged };
				if (opt.eqc_cold)
					e.cold_iterations = converge(*fproc, fr.charge, fr.mult, fkind, nullptr, nullptr).iterations;
				out.rows.push_back(e);
			}
			return out;
		}

		//Mode 2: all from the file; wfn, gbw and molden take E from the ORCA output beside them
		entry file_entry(const std::filesystem::path &path, const std::string &name)
		{
			WFN w(path);
			entry e;
			e.name = name;
			e.program = path.filename().string();
			double E = w.get_total_energy();
			if (E == 0.0)
			{
				E = orca_energy((path.parent_path() / (path.stem().string() + ".out")).string());
				e.program += " + " + path.stem().string() + ".out";
			}
			err_checkf(!std::isnan(E), "-eqc_wfn: " + path.string() + " carries no total energy, and no ORCA .out with one sits next to it", std::cout);
			e.t.E = E;
			double n = 0.0;
			for (int k = 0; k < w.get_nmo(); k++)
			{
				e.t.nX += w.get_MO_occ(k) * w.get_MO_energy(k);
				n += w.get_MO_occ(k);
			}
			ivec Z;
			std::vector<d3> pos;
			int ecp = 0;
			for (int a = 0; a < w.get_ncen(); a++)
			{
				const atom at = w.get_atom(a);
				Z.push_back(at.get_charge());
				const d3 p = at.get_pos();
				pos.push_back({ p[0], p[1], p[2] });
				ecp += at.get_ECP_electrons();
			}
			e.t.Vnn = nuclear_repulsion(Z, pos);
			e.t.n = static_cast<int>(std::lround(n)) + ecp;
			return e;
		}

		void print(const std::vector<entry> &rows, const std::string &mode_line)
		{
			using std::setw;
			std::cout << "\n" << mode_line << "\n\n"
				<< "Terms per wavefunction (eV; Eee = Vnn - (E - nX))\n"
				<< std::left << setw(12) << "" << std::right << setw(5) << "n" << setw(18) << "E" << setw(18) << "nX"
				<< setw(18) << "Vnn" << setw(18) << "Eee" << setw(7) << "iter" << setw(7) << "cold" << "  source\n";
			for (const entry &e : rows)
			{
				std::cout << std::left << setw(12) << e.name << std::right << setw(5) << e.t.n << std::fixed << std::setprecision(6)
					<< setw(18) << e.t.E * hartree2eV << setw(18) << e.t.nX * hartree2eV << setw(18) << e.t.Vnn * hartree2eV
					<< setw(18) << e.t.Eee() * hartree2eV;
				const auto it = [](const int i) { return i < 0 ? std::string("-") : std::to_string(i); };
				std::cout << setw(7) << it(e.iterations) << setw(7) << it(e.cold_iterations) << "  " << e.program
					<< (e.converged ? "" : " (NOT CONVERGED)") << "\n";
			}
			std::vector<terms> reactants;
			for (size_t i = 1; i < rows.size(); i++) reactants.push_back(rows[i].t);
			const reaction r = react(reactants, rows[0].t);
			const double n = r.n;
			std::cout << "\nFragments -> parent, products minus reactants (eV), as X-analysis -m 2 prints it\n"
				<< std::setprecision(6)
				<< "  Delta E          " << setw(14) << r.dE * hartree2eV << "\n"
				<< "  Delta(nX-bar)    " << setw(14) << r.dnX * hartree2eV << "\n"
				<< "  Delta(Vnn-Eee)   " << setw(14) << r.dV() * hartree2eV << "\n"
				<< "  Delta(E/n)       " << setw(14) << r.dE * hartree2eV / n << "\n"
				<< "  Delta(X-bar)     " << setw(14) << r.dnX * hartree2eV / n << "\n"
				<< "  Delta(Vnn/n)     " << setw(14) << r.dVnn * hartree2eV / n << "\n"
				<< "  Delta(Eee/n)     " << setw(14) << r.dEee * hartree2eV / n << "\n"
				<< "  Delta(Vnn-Eee)/n " << setw(14) << r.dV() * hartree2eV / n << "\n"
				<< "  Q                " << setw(14) << r.Q << "\n"
				<< "  Delta|Eee/E| (%) " << setw(14) << r.dEeeE << "\n"
				<< "  Covalency (%)    " << setw(14) << r.covalency << "\n"
				<< "  Ionicity (%)     " << setw(14) << 100.0 - r.covalency << std::endl;
		}
	}

	int run(const options &opt)
	{
		const section_log::section file(opt.eqc_wfns.empty() ? opt.wfn : opt.eqc_wfns[0], "eqc", "EQC energy decomposition analysis", { citations::Method::EQC }, opt.no_date);
		const section_log::tee to_file("eqc");
		citations::cite(citations::Method::EQC);
		std::vector<entry> rows;
		std::string mode_line;
		if (!opt.eqc_wfns.empty())
		{
			err_checkf(opt.eqc_wfns.size() >= 3, "-eqc_wfn wants the parent first and then at least two fragments", std::cout);
			std::cout << "EQC analysis, Mode 2: parent and fragments read from files" << std::endl;
			for (size_t i = 0; i < opt.eqc_wfns.size(); i++)
				rows.push_back(file_entry(opt.eqc_wfns[i], i == 0 ? "parent" : "fragment " + std::to_string(i)));
			mode_line = "Mode 2 (external wavefunctions): every term read from its file, no SCF run";
		}
		else
		{
			err_checkf(opt.eqc_frags.size() >= 2, "-eqc needs at least two -eqc_frag <atoms> <charge> <multiplicity>, or -eqc_wfn files", std::cout);
			WFN w(opt.wfn);
			err_checkf(w.get_origin() == e_origin::gbw, "-eqc with -eqc_frag starts from an ORCA .gbw; for other files use -eqc_wfn", std::cout);
			ivec seen(w.get_ncen(), 0);
			for (const auto &f : opt.eqc_frags)
				for (const int a : f.atoms)
				{
					err_checkf(a >= 0 && a < w.get_ncen(), "-eqc_frag: atom " + std::to_string(a) + " does not exist (indices are 0-based)", std::cout);
					seen[a]++;
				}
			for (int a = 0; a < w.get_ncen(); a++)
				err_checkf(seen[a] == 1, "-eqc_frag: atom " + std::to_string(a) + " is in " + std::to_string(seen[a]) + " fragments, it must be in exactly one", std::cout);
			std::vector<occ::core::Atom> atoms;
			int Z = 0, ecp = 0;
			for (int a = 0; a < w.get_ncen(); a++)
			{
				const atom at = w.get_atom(a);
				const d3 p = at.get_pos();
				atoms.push_back({ at.get_charge(), p[0], p[1], p[2] });
				Z += at.get_charge();
				ecp += at.get_ECP_electrons();
			}
			double n_explicit = 0.0, n_a = 0.0, n_b = 0.0;
			for (int k = 0; k < w.get_nmo(); k++)
			{
				const double o = w.get_MO_occ(k);
				n_explicit += o;
				if (w.get_MO_op(k) == 1) n_b += o;
				else if (w.get_is_unrestricted()) n_a += o;
			}
			const int charge = Z - ecp - static_cast<int>(std::lround(n_explicit));
			const int mult = w.get_is_unrestricted() ? static_cast<int>(std::lround(std::abs(n_a - n_b))) + 1 : 1;
			occ::gto::AOBasis basis;
			if (!opt.eqc_basis.empty())
			{
				basis = occ::gto::AOBasis::load(atoms, opt.eqc_basis);
				basis.set_pure(true);
				const std::string why = basis_mismatch(basis, gbw_basis(w, atoms));
				err_checkf(why.empty(), "-eqc_basis " + opt.eqc_basis + " is not the basis of the gbw: " + why, std::cout);
			}
			else
			{
				err_checkf(ecp == 0, "The gbw uses ECPs, whose parameters it does not store in a form occ reads; name the basis with -eqc_basis (e.g. -eqc_basis def2-svp)", std::cout);
				basis = gbw_basis(w, atoms);
			}
			const MolecularOrbitals gbw_mo = gbw_orbitals(w, basis);
			std::cout << "EQC analysis, Mode 1: " << opt.wfn.filename().string() << ", charge " << charge << ", multiplicity " << mult
				<< ", " << basis.nbf() << " basis functions (" << (opt.eqc_basis.empty() ? "from the gbw" : opt.eqc_basis) << "), "
				<< opt.eqc_frags.size() << " fragments, method " << opt.eqc_method << std::endl;
			occ::log::set_log_file((opt.wfn.parent_path() / (opt.wfn.stem().string() + "_eqc_occ.log")).string());
			occ::parallel::set_num_threads(opt.threads > 0 ? opt.threads : omp_get_max_threads());
			const mode1_out m = opt.eqc_method == "hf" ? mode1<occ::qm::HartreeFock>(opt, w, basis, gbw_mo, charge, mult)
				: mode1<occ::dft::DFT>(opt, w, basis, gbw_mo, charge, mult);
			rows = m.rows;
			const double e_occ = rows[0].t.E;
			std::ostringstream s;
			//E(ORCA) comes from the .out beside the gbw, which need not be there
			const bool have_orca = !std::isnan(m.e_orca);
			s << std::fixed << std::setprecision(9) << "Parent: E(ORCA) ";
			if (have_orca) s << m.e_orca;
			else s << "n/a (no " << opt.wfn.stem().string() << ".out)";
			s << ", E(occ) at the gbw orbitals " << m.e_seed << ", E(occ) converged " << e_occ << " hartree\n        "
				<< std::scientific << std::setprecision(3);
			if (have_orca) s << "E(occ) - E(ORCA) " << e_occ - m.e_orca << ", ";
			s << "E(occ, gbw orbitals) - E(occ) " << m.e_seed - e_occ << ", nX(occ) - nX(ORCA) " << rows[0].t.nX - m.nx_orca
				<< "\n        gbw max|C^T S C - 1| " << m.ortho << ", " << rows[0].iterations << " iterations to re-converge\n"
				<< "Mode 1 (seeded): parent from the ORCA gbw re-converged in occ, fragments converged in occ from the parent's density block; SCF iterations";
			for (const entry &e : rows) s << " " << e.name << " " << e.iterations << (e.cold_iterations >= 0 ? " (cold " + std::to_string(e.cold_iterations) + ")" : "");
			mode_line = s.str();
		}
		print(rows, mode_line);
		return 0;
	}
}
