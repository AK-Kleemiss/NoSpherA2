#include "pch.h"
#include "eli_family.h"
#include "constants.h"
#include "citations.h"

namespace eli_family
{
	SpinSplit spin_split(const WFN& wave)
	{
		const int nmo = wave.get_nmo();
		int singles = 0;
		for (int mo = 0; mo < nmo; mo++)
		{
			if (wave.get_MO_op(mo) == 1) return SpinSplit::unrestricted;
			const double occ = wave.get_MO_occ(mo);
			if (occ == 1.0) singles++;
			else if (occ != 0.0 && occ != 2.0) return SpinSplit::halves;
		}
		const int stated = (int)wave.get_multi();
		return singles > 0 && (stated == 0 || stated == singles + 1) ? SpinSplit::restricted_open : SpinSplit::halves;
	}

	bool alpha_beta_orbitals_differ(const WFN& wave)
	{
		std::vector<int> occ[2];
		for (int mo = 0; mo < wave.get_nmo(); mo++)
			if (wave.get_MO_occ(mo) != 0.0) occ[wave.get_MO_op(mo) == 1 ? 1 : 0].push_back(mo);
		if (occ[0].size() != occ[1].size()) return true;
		const int nex = wave.get_nex();
		bool differ = false;
		for (size_t k = 0; k < occ[0].size() && !differ; k++)
		{
			const int a = occ[0][k], b = occ[1][k];
			if (wave.get_MO_occ(a) != wave.get_MO_occ(b)) return true;
			double big = 0.0, same = 0.0, flip = 0.0;
			for (int j = 0; j < nex; j++)
			{
				const double ca = wave.get_MO_coef(a, j), cb = wave.get_MO_coef(b, j);
				big = std::max(big, std::abs(ca));
				same = std::max(same, std::abs(ca - cb));
				flip = std::max(flip, std::abs(ca + cb));
			}
			differ = std::min(same, flip) > 1e-4 * big;
		}
		if (!differ) return false;
		//Differing orbitals can span the same space, so the spin density decides: 0.3/1/2 bohr from each
		//nucleus along two skew directions.
		//ponytail: 6 points per atom, a spin density that vanishes at all of them is missed
		const std::vector<atom>& atoms = *wave.get_atoms_ptr();
		const double dir[2][3] = { { 0.48, 0.64, 0.6 }, { -0.6, 0.48, -0.64 } }, rad[3] = { 0.3, 1.0, 2.0 };
		bool polarised = false;
		wave.get_coef_primitive_major(); //build the lazy cache before the threads read it
#pragma omp parallel for reduction(||:polarised)
		for (int i = 0; i < (int)atoms.size(); i++)
			for (int s = 0; s < 6 && !polarised; s++)
			{
				d3 p;
				for (int c = 0; c < 3; c++) p[c] = atoms[i].get_coordinate(c) + rad[s % 3] * dir[s / 3][c];
				SpinFields f;
				spin_fields(wave, p, f);
				polarised = std::abs(f.rho[0] - f.rho[1]) > 1e-10 + 1e-6 * (f.rho[0] + f.rho[1]);
			}
		return polarised;
	}

	void spin_fields(const WFN& wave, const d3& p, SpinFields& f)
	{
		const int _nmo = wave.get_nmo();
		const int _nex = wave.get_nex();
		const std::vector<atom>& atoms = *wave.get_atoms_ptr();
		//Primitive-major: every MO for one primitive is contiguous
		const double* const coefs = wave.get_coef_primitive_major();

		vec phi(4 * (size_t)_nmo, 0.0);
		double chi[4]{ 0, 0, 0, 0 };
		double d[3]{ 0, 0, 0 };
		int l[3]{ 0, 0, 0 };
		double xl[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

		for (int j = 0; j < _nex; j++)
		{
			const int iat = wave.get_center(j) - 1;
			constants::type2vector(wave.get_type(j), l);
			d[0] = p[0] - atoms[iat].get_coordinate(0);
			d[1] = p[1] - atoms[iat].get_coordinate(1);
			d[2] = p[2] - atoms[iat].get_coordinate(2);
			const double temp = -wave.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
			if (temp < constants::exp_cutoff) continue;
			const double ex = exp(temp);
			for (int k = 0; k < 3; k++)
			{
				const double dk = d[k], d2 = dk * dk;
				switch (l[k])
				{
				case 0: xl[k][0] = 1.0;          xl[k][1] = 0.0;           break;
				case 1: xl[k][0] = dk;           xl[k][1] = 1.0;           break;
				case 2: xl[k][0] = d2;           xl[k][1] = 2 * dk;        break;
				case 3: xl[k][0] = d2 * dk;      xl[k][1] = 3 * d2;        break;
				case 4: xl[k][0] = d2 * d2;      xl[k][1] = 4 * d2 * dk;   break;
				case 5: xl[k][0] = d2 * d2 * dk; xl[k][1] = 5 * d2 * d2;   break;
				case 6: xl[k][0] = d2 * d2 * d2; xl[k][1] = 6 * d2 * d2 * dk; break;
				default: return;
				}
			}
			const double ex2 = 2 * wave.get_exponent(j);
			chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
			chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
			chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
			chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;

			const double* c_row = coefs + (size_t)j * _nmo;
			double* phi_ptr = phi.data();
			for (int mo = 0; mo < _nmo; ++mo, phi_ptr += 4)
			{
				const double c = c_row[mo];
				for (int k = 0; k < 4; k++) phi_ptr[k] += c * chi[k];
			}
		}

		const SpinSplit how = spin_split(wave);
		f = SpinFields{};
		for (int mo = 0; mo < _nmo; mo++)
		{
			const double occ = wave.get_MO_occ(mo);
			if (occ == 0.0) continue;
			const double* ph = &phi[(size_t)mo * 4];
			double nc[2];
			mo_spin_occupations(occ, wave.get_MO_op(mo), how, nc);
			for (int c = 0; c < 2; c++)
			{
				const double n = nc[c];
				if (n == 0.0) continue;
				f.rho[c] += n * ph[0] * ph[0];
				f.grad[c][0] += 2 * n * ph[0] * ph[1];
				f.grad[c][1] += 2 * n * ph[0] * ph[2];
				f.grad[c][2] += 2 * n * ph[0] * ph[3];
				f.T[c] += n * (ph[1] * ph[1] + ph[2] * ph[2] + ph[3] * ph[3]);
			}
		}
	}

	//rho^(t) / rho = 1 - N_beta / (2(N-1)), DGrid 5.2's convention
	double triplet_density_factor(const WFN& wave)
	{
		const SpinSplit how = spin_split(wave);
		double N = 0.0, Nb = 0.0;
		for (int mo = 0; mo < wave.get_nmo(); mo++)
		{
			const double occ = wave.get_MO_occ(mo);
			double n[2];
			mo_spin_occupations(occ, wave.get_MO_op(mo), how, n);
			N += occ;
			Nb += n[1];
		}
		return N > 1.0 ? 1.0 - Nb / (2.0 * (N - 1.0)) : 0.0;
	}

	const char* member_name(const Member m)
	{
		switch (m)
		{
		case Member::eli_d_aa:      return "eli_d_aa";
		case Member::eli_d_bb:      return "eli_d_bb";
		case Member::eli_d_triplet: return "eli_d_triplet";
		case Member::eli_q_aa:      return "eli_q_aa";
		case Member::eli_q_bb:      return "eli_q_bb";
		case Member::elia_singlet:  return "elia_singlet";
		}
		return "?";
	}

	const char* member_column(const Member m)
	{
		switch (m)
		{
		case Member::eli_d_aa:      return "ELI-D_aa";
		case Member::eli_d_bb:      return "ELI-D_bb";
		case Member::eli_d_triplet: return "ELI-D_t";
		case Member::eli_q_aa:      return "ELI-q_aa";
		case Member::eli_q_bb:      return "ELI-q_bb";
		case Member::elia_singlet:  return "ELIA_s";
		}
		return "?";
	}

	//ELI-q = ELI-D^(-8/3) ascends into ELI-D's minima and shatters into tail basins; 1 - zeta^2 has no pair topology
	bool basins_independent(const Member m)
	{
		return m == Member::eli_d_aa || m == Member::eli_d_bb || m == Member::eli_d_triplet;
	}

	std::vector<Member> eli_variants_for(const WFN& wave, std::string* warning)
	{
		//From the orbitals, never the label; a broken-symmetry singlet counts as spin-polarised
		double Na = 0.0, Nb = 0.0;
		const SpinSplit how = spin_split(wave);
		const bool has_beta_set = how == SpinSplit::unrestricted;
		for (int mo = 0; mo < wave.get_nmo(); mo++)
		{
			const double occ = wave.get_MO_occ(mo);
			if (occ == 0.0) continue;
			if (how == SpinSplit::restricted_open) { double n[2]; mo_spin_occupations(occ, 0, how, n); Na += n[0]; Nb += n[1]; }
			else if (has_beta_set && wave.get_MO_op(mo) == 1) Nb += occ; else Na += occ;
		}
		const bool spin_polarised = how == SpinSplit::restricted_open
			|| (has_beta_set && (std::abs(Na - Nb) > 1e-8 || alpha_beta_orbitals_differ(wave)));

		const bool has_alpha_electrons = Na > 1e-8;
		const bool has_beta_electrons = Nb > 1e-8;
		const bool has_a_pair = Na + Nb > 1.0 + 1e-8;
		std::vector<Member> out;
		if (has_alpha_electrons) out.push_back(Member::eli_d_aa);
		if (spin_polarised && has_beta_electrons) out.push_back(Member::eli_d_bb);
		if (spin_polarised && has_beta_electrons && has_a_pair) out.push_back(Member::eli_d_triplet);
		if (has_alpha_electrons) out.push_back(Member::eli_q_aa);
		if (spin_polarised && has_beta_electrons) out.push_back(Member::eli_q_bb);
		if (spin_polarised && has_beta_electrons) out.push_back(Member::elia_singlet);
		if (warning)
		{
			warning->clear();
			if (has_beta_set && !has_beta_electrons)
				*warning = "the beta orbital set of this wavefunction holds no electrons (N_alpha = "
				+ std::to_string(Na) + ", N_beta = 0), so ELI-D(bb), ELI-q(bb), the triplet member"
				" and the singlet ELI-q are identically zero everywhere and are not listed; the"
				" alpha-alpha members are the whole family this input has. ";
			if (!has_alpha_electrons)
				*warning += "the alpha orbital set holds no electrons either: there is no ELI"
				" family to compute from this wavefunction at all. ";
			const int stated = (int)wave.get_multi();
			const int actual = (int)std::lround(std::abs(Na - Nb)) + 1;
			if (stated > 0 && stated != actual)
				*warning += "multiplicity label says " + std::to_string(stated)
				+ " but the orbital occupations give " + std::to_string(actual)
				+ (spin_polarised ? "" : " (no beta orbital set: this file is restricted)")
				+ "; the ELI member list follows the occupations.";
		}
		return out;
	}

	const char* elia_status()
	{
		return "ELIA (antiparallel-spin ELI, Kohout et al., Theor. Chem. Acc. 113 (2005) 287) is not\n"
			"computable from this wavefunction. It samples the opposite-spin pair density, whose\n"
			"leading term is the on-top density rho_2^ab(r,r). For a single determinant that is\n"
			"exactly rho_alpha*rho_beta, so the antiparallel pair population per micro-cell is\n"
			"uniform and ELIA is constant throughout space. The same holds for the singlet ELI-q of\n"
			"Eq. 53 of Part III, which is identically 1 for a closed-shell determinant. Both need the\n"
			"alpha-beta block of a correlated 2-matrix, which no wfn/wfx/molden/fchk file carries.\n"
			"DGrid 5.2 agrees: fed a single-determinant molden it writes an identically zero ELIA\n"
			"field and reports 'ALPHA-BETA 2-matrix element missing'.";
	}

	void report(const std::filesystem::path& wfn_path, const std::filesystem::path& points_file)
	{
		WFN wave(wfn_path);
		const SpinSplit how = spin_split(wave);
		const bool unrestricted = how == SpinSplit::unrestricted;
		double N = 0.0, Nb = 0.0;
		for (int mo = 0; mo < wave.get_nmo(); mo++)
		{
			const double occ = wave.get_MO_occ(mo);
			double n[2];
			mo_spin_occupations(occ, wave.get_MO_op(mo), how, n);
			N += occ;
			Nb += n[1];
		}
		const double factor = triplet_density_factor(wave);
		std::string warning;
		const std::vector<Member> members = eli_variants_for(wave, &warning);
		citations::cite(citations::Method::ELIFamily, std::cout);
		std::cout << "ELI family for " << wfn_path.string() << "\n"
			<< "  " << (unrestricted ? "unrestricted" : how == SpinSplit::restricted_open ? "restricted open shell" : "restricted") << ", N = " << N
			<< ", N_alpha = " << N - Nb << ", N_beta = " << Nb << "\n"
			<< "  members worth computing:";
		for (const Member m : members)
			std::cout << " " << member_name(m) << (basins_independent(m) ? "" : "(no basins)");
		std::cout << "\n";
		if (!warning.empty()) std::cout << "  WARNING: " << warning << "\n";
		std::cout << "  rho^(t)/rho = " << factor << "  (DGrid 5.2 convention, 1 - N_beta/(2(N-1)))\n"
			<< "  ELI-D(UEG) = " << ueg_eli_d() << ", ELI-q(UEG) = " << ueg_eli_q() << "\n"
			<< "  " << elia_status() << std::endl;
		if (points_file.empty()) return;
		std::ifstream in(points_file);
		err_checkf(in.good(), "Could not open " + points_file.string(), std::cout);
		//All ingredients, T_s = sum n_i |grad phi_i|^2, then the members, so the dump compares against DGrid
		std::cout << "#         x           y           z         rho_a         rho_b"
			<< "           T_a           T_b     |grad_a|^2     |grad_b|^2      grad_a.b"
			<< "      ELI-D_aa      ELI-D_bb      ELI-q_aa      ELI-q_bb       ELI-D_t        ELIA_s\n";
		std::string line;
		while (std::getline(in, line))
		{
			if (line.empty() || line[0] == '#') continue;
			std::istringstream s(line);
			d3 p{ 0, 0, 0 };
			if (!(s >> p[0] >> p[1] >> p[2])) continue;
			SpinFields f;
			spin_fields(wave, p, f);
			const double ga = g_same_spin(f, 0), gb = g_same_spin(f, 1);
			double dd[3]{ 0, 0, 0 };
			for (int k = 0; k < 3; k++) {
				dd[0] += f.grad[0][k] * f.grad[0][k];
				dd[1] += f.grad[1][k] * f.grad[1][k];
				dd[2] += f.grad[0][k] * f.grad[1][k];
			}
			std::cout << std::fixed << std::setprecision(6) << std::setw(11) << p[0] << std::setw(12) << p[1] << std::setw(12) << p[2]
				<< std::scientific << std::setprecision(6)
				<< std::setw(14) << f.rho[0] << std::setw(14) << f.rho[1]
				<< std::setw(14) << f.T[0] << std::setw(14) << f.T[1]
				<< std::setw(14) << dd[0] << std::setw(14) << dd[1] << std::setw(14) << dd[2]
				<< std::setw(14) << eli_d(f.rho[0], ga) << std::setw(14) << eli_d(f.rho[1], gb)
				<< std::setw(14) << eli_q(f.rho[0], ga) << std::setw(14) << eli_q(f.rho[1], gb)
				<< std::setw(14) << eli_d(factor * (f.rho[0] + f.rho[1]), g_triplet(f))
				<< std::setw(14) << elia_singlet_eli_q(f) << "\n";
		}
		std::cout << std::flush;
	}
}
