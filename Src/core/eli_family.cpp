#include "pch.h"
#include "eli_family.h"
#include "constants.h"

//The spin-resolved evaluator.  Structurally WFN::computeELIELF (wfn_density.cpp), i.e. the same
//4-component primitive loop (value + three Cartesian derivatives), but the MO accumulation is split
//into the alpha and beta channels instead of being spin-summed.  Free function on the public WFN
//accessors, so nothing in wfn_class.h / wfn_density.cpp has to change.
namespace eli_family
{
	void spin_fields(const WFN& wave, const d3& p, SpinFields& f)
	{
		const int _nmo = wave.get_nmo();
		const int _nex = wave.get_nex();
		const std::vector<atom>& atoms = *wave.get_atoms_ptr();
		//Primitive-major coefficients: every MO for one primitive is contiguous (see computeELIELF).
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

		//A restricted wavefunction carries no beta MOs; its doubly occupied orbitals contribute
		//occ/2 to each channel, giving rho_alpha = rho_beta = rho/2 exactly.
		const bool unrestricted = wave.get_MO_op_count(1) > 0;
		f = SpinFields{};
		for (int mo = 0; mo < _nmo; mo++)
		{
			const double occ = wave.get_MO_occ(mo);
			if (occ == 0.0) continue;
			const double* ph = &phi[(size_t)mo * 4];
			const double n = unrestricted ? occ : 0.5 * occ;
			const int s = unrestricted ? wave.get_MO_op(mo) : 0;
			const int lo = unrestricted ? s : 0, hi = unrestricted ? s : 1;
			for (int c = lo; c <= hi; c++)
			{
				f.rho[c] += n * ph[0] * ph[0];
				f.grad[c][0] += 2 * n * ph[0] * ph[1];
				f.grad[c][1] += 2 * n * ph[0] * ph[2];
				f.grad[c][2] += 2 * n * ph[0] * ph[3];
				f.T[c] += n * (ph[1] * ph[1] + ph[2] * ph[2] + ph[3] * ph[3]);
			}
		}
	}

	//rho^(t) / rho = 1 - N_beta / (2(N-1)), DGrid 5.2's convention - see the header.
	double triplet_density_factor(const WFN& wave)
	{
		const bool unrestricted = wave.get_MO_op_count(1) > 0;
		double N = 0.0, Nb = 0.0;
		for (int mo = 0; mo < wave.get_nmo(); mo++)
		{
			const double occ = wave.get_MO_occ(mo);
			N += occ;
			if (unrestricted) { if (wave.get_MO_op(mo) == 1) Nb += occ; }
			else Nb += 0.5 * occ;
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

	//Only ELI-D has maxima to walk up to.  ELI-q is ELI-D^(-8/3), so an ascent on it runs into ELI-D's
	//minima and shatters into tail basins; the singlet member is 1 - zeta^2 and has no pair topology at
	//all.  The measurements behind this are in the header, above eli_variants_for's declaration.
	bool basins_independent(const Member m)
	{
		return m == Member::eli_d_aa || m == Member::eli_d_bb || m == Member::eli_d_triplet;
	}

	std::vector<Member> eli_variants_for(const WFN& wave, std::string* warning)
	{
		//Read off the wavefunction, never off a label: a beta MO set must actually be present, and the
		//two channels must actually hold different numbers of electrons.
		double Na = 0.0, Nb = 0.0;
		const bool has_beta_set = wave.get_MO_op_count(1) > 0;
		for (int mo = 0; mo < wave.get_nmo(); mo++)
		{
			const double occ = wave.get_MO_occ(mo);
			if (occ == 0.0) continue;
			if (has_beta_set && wave.get_MO_op(mo) == 1) Nb += occ; else Na += occ;
		}
		const bool spin_polarised = has_beta_set && std::abs(Na - Nb) > 1e-8;

		std::vector<Member> out{ Member::eli_d_aa, Member::eli_q_aa };
		if (spin_polarised)
			out = { Member::eli_d_aa, Member::eli_d_bb, Member::eli_d_triplet,
					Member::eli_q_aa, Member::eli_q_bb, Member::elia_singlet };
		if (warning)
		{
			warning->clear();
			//The stated multiplicity is a cross-check only; 2S+1 from the occupations is the truth.
			const int stated = (int)wave.get_multi();
			const int actual = (int)std::lround(std::abs(Na - Nb)) + 1;
			if (stated > 0 && stated != actual)
				*warning = "multiplicity label says " + std::to_string(stated)
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
		const bool unrestricted = wave.get_MO_op_count(1) > 0;
		double N = 0.0, Nb = 0.0;
		for (int mo = 0; mo < wave.get_nmo(); mo++)
		{
			const double occ = wave.get_MO_occ(mo);
			N += occ;
			Nb += unrestricted ? (wave.get_MO_op(mo) == 1 ? occ : 0.0) : 0.5 * occ;
		}
		const double factor = triplet_density_factor(wave);
		std::string warning;
		const std::vector<Member> members = eli_variants_for(wave, &warning);
		std::cout << "ELI family for " << wfn_path.string() << "\n"
			<< "  " << (unrestricted ? "unrestricted" : "restricted") << ", N = " << N
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
		//Everything the family is built from, so a point dump doubles as the DGrid comparison input:
		//the two densities, the two T_s = sum n_i |grad phi_i|^2, the three gradient invariants, then
		//the members themselves.
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
