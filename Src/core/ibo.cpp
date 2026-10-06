#include "pch.h"
#include "ibo.h"
#include "wfn_class.h"
#include "integration_params.h"
#include "libCintMain.h"
#include "basis_set.h"
#include "constants.h"
#include "citations.h"

#include <Eigen/Dense>

using Eigen::MatrixXd;

namespace
{
	MatrixXd inverse_sqrt(const MatrixXd &M)
	{
		Eigen::SelfAdjointEigenSolver<MatrixXd> es(M);
		return es.operatorInverseSqrt();
	}

	//atom of every spherical function of an Int_Params basis, libcint order
	ivec function_atoms(const Int_Params &p)
	{
		const ivec bas = p.get_bas();
		ivec at;
		for (int s = 0; s < static_cast<int>(p.get_nbas()); s++)
			at.insert(at.end(), (2 * bas[8 * s + 1] + 1) * bas[8 * s + 3], bas[8 * s]);
		return at;
	}

	//Electrons in the closed shells under the valence of element Z: the noble-gas core plus filled d/f shells
	int core_electrons(const int Z)
	{
		static const int lim[] = {2, 10, 18, 30, 36, 48, 54, 71, 80, 86}, core[] = {0, 2, 10, 18, 28, 36, 46, 54, 68, 78};
		for (int i = 0; i < 10; i++)
			if (Z <= lim[i]) return core[i];
		return 86;
	}

	std::string atom_label(const WFN &wavy, const int a)
	{
		return std::string(constants::atnr2letter(wavy.get_atom_charge(a))) + std::to_string(a + 1);
	}

	int parse_atom(const std::string &s)
	{
		size_t p = 0;
		while (p < s.size() && isalpha(static_cast<unsigned char>(s[p]))) p++;
		err_checkf(p < s.size(), "-ibo_cube: no atom number in '" + s + "'", std::cout);
		return std::stoi(s.substr(p)) - 1;
	}
}

IBOResult intrinsic_bond_orbitals(const WFN &wavy)
{
	err_checkf(!wavy.get_is_unrestricted(), "IBO: only closed-shell (RHF/RKS) wavefunctions are supported; this one is unrestricted. "
				"Localise alpha and beta separately in the program that wrote it.", std::cout);
	err_checkf(origin_has_orca_pure_phases(wavy.get_origin()), "IBO needs the contracted spherical basis and its MO coefficients: "
				"read a .gbw or a .molden (spherical), not a .wfn/.wfx/.fchk.", std::cout);
	IBOResult r;
	for (int i = 0; i < wavy.get_nmo(); i++) {
		const double occ = wavy.get_MO_occ(i);
		if (occ < 1e-6) continue;
		err_checkf(std::abs(occ - 2.0) < 1e-6, "IBO: MO " + std::to_string(i + 1) + " has occupation " + std::to_string(occ) +
					"; only closed-shell wavefunctions with occupations 0 or 2 are supported.", std::cout);
		r.mos.push_back(i);
	}
	const int nocc = static_cast<int>(r.mos.size()), ncen = wavy.get_ncen();
	err_checkf(nocc > 0, "IBO: no occupied orbitals", std::cout);
	//molecular basis + MINAO, one combined overlap: S1, S12 and S2 are its blocks
	WFN minao(e_origin::NOT_YET_DEFINED);
	minao.set_atoms(wavy.get_atoms());
	minao.set_ncen(ncen);
	minao.delete_basis_set();
	const std::shared_ptr<BasisSet> mb = BasisSetLibrary::get_basis_set("minao");
	for (int a = 0; a < ncen; a++)
		err_checkf(mb->has_element(wavy.get_atom_charge(a)), std::string("IBO: the MINAO reference basis has no ") +
					constants::atnr2letter(wavy.get_atom_charge(a)) + " (it lacks K, Rb, Sr, Cs, Ba, the lanthanides and Z > 86)", std::cout);
	load_basis_into_WFN(minao, mb);
	Int_Params pw(wavy), pm(minao);
	Int_Params pc(pw, pm);
	vec flat;
	compute2C<Overlap2C_SPH>(pc, flat);
	const int n1 = static_cast<int>(pw.get_nao()), n2 = static_cast<int>(pm.get_nao()), n = n1 + n2;
	err_checkf(static_cast<int>(flat.size()) == n * n, "IBO: combined overlap has the wrong size", std::cout);
	const ivec ao_atom = function_atoms(pw), iao_atom = function_atoms(pm);
	err_checkf(static_cast<int>(ao_atom.size()) == n1 && static_cast<int>(iao_atom.size()) == n2, "IBO: shell layout does not match the AO count", std::cout);
	//ORCA phases on the molecular rows only, as ao_overlap does it
	bvec flip(n, false);
	const ivec bas = pw.get_bas();
	for (int s = 0, k = 0; s < static_cast<int>(pw.get_nbas()); s++)
		for (int m = -bas[8 * s + 1]; m <= bas[8 * s + 1]; m++, k++)
			flip[k] = orca_pure_sign_flips(m);
	MatrixXd S(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			S(i, j) = flip[i] != flip[j] ? -flat[i * n + j] : flat[i * n + j];
	const MatrixXd S1 = S.topLeftCorner(n1, n1), S12 = S.topRightCorner(n1, n2), S2 = S.bottomRightCorner(n2, n2);
	//MO_sph rows are in file shell order (libcint m order); Int_Params groups an atom's shells by l
	const dMatrix2 &sph = wavy.get_MO_sph();
	err_checkf(static_cast<int>(sph.extent(0)) == n1 && static_cast<int>(sph.extent(1)) > r.mos.back(),
				"IBO: this reader kept no spherical MO coefficients matching the basis (" + std::to_string(sph.extent(0)) + " rows, " +
				std::to_string(n1) + " AOs); use a .gbw or a spherical .molden", std::cout);
	ivec row(n1, -1);
	const std::vector<atom> atoms = wavy.get_atoms();
	for (int a = 0, file = 0, internal = 0; a < ncen; a++) {
		const atom &at = atoms[a];
		ivec shell_l, shell_off;
		for (int s = 0, p = 0; s < static_cast<int>(at.get_shellcount_size()); p += at.get_shellcount(s), s++) {
			shell_l.push_back(at.get_basis_set_entry(p).get_type() - 1);
			shell_off.push_back(file);
			file += 2 * shell_l.back() + 1;
		}
		for (int l = 0; l <= 20; l++)
			for (int s = 0; s < static_cast<int>(shell_l.size()); s++)
				if (shell_l[s] == l)
					for (int m = 0; m <= 2 * l; m++)
						row.at(internal++) = shell_off[s] + m;
	}
	err_checkf(*std::min_element(row.begin(), row.end()) >= 0 && *std::max_element(row.begin(), row.end()) < n1,
				"IBO: the atoms' shells do not add up to the AO count", std::cout);
	MatrixXd C(n1, nocc);
	for (int i = 0; i < n1; i++)
		for (int k = 0; k < nocc; k++)
			C(i, k) = sph(row[i], r.mos[k]);
	err_checkf((C.transpose() * S1 * C - MatrixXd::Identity(nocc, nocc)).cwiseAbs().maxCoeff() < 1e-5,
				"IBO: the occupied MOs are not orthonormal over the AO overlap - basis and coefficients do not match", std::cout);
	const dMatrix2 P = wavy.get_dm();
	if (static_cast<int>(P.extent(0)) == n1) {
		const MatrixXd D = 2.0 * C * C.transpose();
		double dev = 0.0;
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < n1; j++)
				dev = std::max(dev, std::abs(D(i, j) - P(i, j)));
		err_checkf(dev < 1e-5, "IBO: the MO coefficients do not reproduce the density matrix (max dev " + std::to_string(dev) + ")", std::cout);
	}
	//IAOs, Knizia eq. 1 in the form of PySCF lo/iao.py
	Eigen::LLT<MatrixXd> l1(S1), l2(S2);
	const MatrixXd P12 = l1.solve(S12);
	MatrixXd Ct = l1.solve(S12 * l2.solve(S12.transpose() * C));
	Ct = Ct * inverse_sqrt(Ct.transpose() * S1 * Ct);
	const MatrixXd CCS = C * (C.transpose() * S1), CtS = Ct * (Ct.transpose() * S1);
	MatrixXd A = P12 + 2.0 * CCS * (CtS * P12) - CCS * P12 - CtS * P12;
	A = A * inverse_sqrt(A.transpose() * S1 * A);
	MatrixXd B = A.transpose() * S1 * C;
	const double span = (A * B - C).cwiseAbs().maxCoeff();
	err_checkf(n2 >= nocc && span < 1e-5, "IBO: the IAOs do not span the occupied space (max dev " + std::to_string(span) +
				"); MINAO and the wavefunction disagree on the electrons per atom - an ECP calculation needs the PP MINAO "
				"(Z > 38), an all-electron one Z <= 36", std::cout);
	r.charge.assign(ncen, 0.0);
	for (int a = 0; a < ncen; a++)
		r.charge[a] = wavy.get_atom_charge(a) - wavy.get_atom_ECP_electrons(a);
	for (int mu = 0; mu < n2; mu++)
		r.charge[iao_atom[mu]] -= 2.0 * B.row(mu).squaredNorm();
	//Jacobi sweeps, exponent 4 (Knizia's ibo-ref)
	ivec first(ncen + 1, n2);
	for (int mu = n2 - 1; mu >= 0; mu--) first[iao_atom[mu]] = mu;
	for (int a = ncen - 1; a >= 0; a--) first[a] = std::min(first[a], first[a + 1]);
	auto pop = [&](const int a, const int i, const int j) { return B.block(first[a], i, first[a + 1] - first[a], 1).cwiseProduct(B.block(first[a], j, first[a + 1] - first[a], 1)).sum(); };
	auto functional = [&]() {
		double f = 0.0;
		for (int k = 0; k < nocc; k++)
			for (int a = 0; a < ncen; a++) f += std::pow(pop(a, k, k), 4);
		return f;
	};
	r.functional_start = functional();
	MatrixXd U = MatrixXd::Identity(nocc, nocc);
	for (r.sweeps = 1; r.sweeps <= 1000; r.sweeps++) {
		double grad = 0.0;
		for (int i = 0; i < nocc; i++)
			for (int j = 0; j < i; j++) {
				double Aij = 0.0, Bij = 0.0;
				for (int a = 0; a < ncen; a++) {
					const double qii = pop(a, i, i), qjj = pop(a, j, j), qij = pop(a, i, j);
					Aij += -std::pow(qii, 4) - std::pow(qjj, 4) + 6.0 * (qii * qii + qjj * qjj) * qij * qij + qii * qii * qii * qjj + qii * qjj * qjj * qjj;
					Bij += 4.0 * qij * (qii * qii * qii - qjj * qjj * qjj);
				}
				grad += Bij * Bij;
				const double phi = 0.25 * std::atan2(Bij, -Aij), c = std::cos(phi), s = std::sin(phi);
				const Eigen::VectorXd bi = B.col(i), ui = U.col(i);
				B.col(i) = c * bi + s * B.col(j);
				B.col(j) = -s * bi + c * B.col(j);
				U.col(i) = c * ui + s * U.col(j);
				U.col(j) = -s * ui + c * U.col(j);
			}
		if (std::sqrt(grad) < 1e-10) break;
	}
	err_checkf(r.sweeps <= 1000, "IBO: localisation did not converge in 1000 sweeps", std::cout);
	r.functional = functional();
	//populations, energies, classification
	const MatrixXd L = C * U, SL = S1 * L;
	MatrixXd iao = MatrixXd::Zero(ncen, nocc), mul = MatrixXd::Zero(ncen, nocc);
	for (int k = 0; k < nocc; k++) {
		for (int a = 0; a < ncen; a++) iao(a, k) = pop(a, k, k);
		for (int mu = 0; mu < n1; mu++) mul(ao_atom[mu], k) += L(mu, k) * SL(mu, k);
	}
	vec e(nocc, 0.0);
	std::vector<IBOKind> kind(nocc, IBOKind::Delocalised);
	ivec2 cen(nocc);
	for (int k = 0; k < nocc; k++) {
		for (int i = 0; i < nocc; i++) e[k] += U(i, k) * U(i, k) * wavy.get_MO_energy(r.mos[i]);
		ivec order(ncen);
		std::iota(order.begin(), order.end(), 0);
		std::sort(order.begin(), order.end(), [&](int x, int y) { return iao(x, k) > iao(y, k); });
		const double top = iao(order[0], k), two = ncen > 1 ? top + iao(order[1], k) : top;
		if (top >= 0.95) {
			kind[k] = IBOKind::LonePair;
			cen[k] = {order[0]};
		}
		else if (two >= 0.85) {
			kind[k] = IBOKind::Bond;
			cen[k] = {std::min(order[0], order[1]), std::max(order[0], order[1])};
		}
		else
			cen[k] = {order[0], order[1]};
	}
	for (int a = 0; a < ncen; a++) {
		ivec one;
		for (int k = 0; k < nocc; k++)
			if (kind[k] == IBOKind::LonePair && cen[k][0] == a) one.push_back(k);
		std::sort(one.begin(), one.end(), [&](int x, int y) { return e[x] < e[y]; });
		const int ncore = std::max(0, core_electrons(wavy.get_atom_charge(a)) - wavy.get_atom_ECP_electrons(a)) / 2;
		for (int c = 0; c < std::min(ncore, static_cast<int>(one.size())); c++) kind[one[c]] = IBOKind::Core;
	}
	ivec order(nocc);
	std::iota(order.begin(), order.end(), 0);
	std::sort(order.begin(), order.end(), [&](int x, int y) {
		if (kind[x] != kind[y]) return kind[x] < kind[y];
		if (cen[x] != cen[y]) return cen[x] < cen[y];
		return e[x] < e[y];
	});
	r.U = dMatrix2(nocc, nocc);
	r.iao_pop = dMatrix2(ncen, nocc);
	r.mulliken = dMatrix2(ncen, nocc);
	for (int k = 0; k < nocc; k++) {
		const int o = order[k];
		//sign convention: the largest AO coefficient is positive
		int big = 0;
		for (int mu = 0; mu < n1; mu++)
			if (std::abs(L(mu, o)) > std::abs(L(big, o))) big = mu;
		const double sg = L(big, o) < 0 ? -1.0 : 1.0;
		for (int i = 0; i < nocc; i++) r.U(i, k) = sg * U(i, o);
		for (int a = 0; a < ncen; a++) r.iao_pop(a, k) = iao(a, o), r.mulliken(a, k) = mul(a, o);
		r.kind.push_back(kind[o]);
		r.centres.push_back(cen[o]);
		r.energy.push_back(e[o]);
	}
	return r;
}

void print_ibo(const IBOResult &r, const WFN &wavy, std::ostream &file)
{
	using namespace std;
	const int ncen = wavy.get_ncen(), nocc = static_cast<int>(r.mos.size());
	file << "\nIntrinsic atomic orbitals and intrinsic bond orbitals (IAO/IBO, MINAO reference, exponent 4)\n";
	citations::cite(citations::Method::IBO, file);
	file << "\nIAO partial charges\n" << fixed << setprecision(5);
	double total = 0.0;
	for (int a = 0; a < ncen; a++) {
		file << "  " << left << setw(6) << atom_label(wavy, a) << right << setw(11) << r.charge[a] << "\n";
		total += r.charge[a];
	}
	file << "  " << left << setw(6) << "Total" << right << setw(11) << total << "\n";
	file << "\nLocalisation sum_k sum_A (q_A^k)^4: " << setprecision(10) << r.functional_start << " (canonical) -> " << r.functional
		 << " after " << r.sweeps << " sweeps, " << nocc << " doubly occupied orbitals\n" << setprecision(5);
	file << "IBO, orbital energy (Eh, Fock diagonal), then each atom with IAO population >= 0.01 and its Mulliken population [ ]\n";
	const char *title[] = {"Core", "Lone pairs", "Bond-like (two-centre, top two >= 0.85)", "Delocalised"};
	for (int t = 0; t < 4; t++) {
		bool head = false;
		for (int k = 0; k < nocc; k++) {
			if (static_cast<int>(r.kind[k]) != t) continue;
			if (!head) file << "\n" << title[t] << "\n", head = true;
			ivec at(ncen);
			iota(at.begin(), at.end(), 0);
			sort(at.begin(), at.end(), [&](int x, int y) { return r.iao_pop(x, k) > r.iao_pop(y, k); });
			file << setw(5) << k + 1 << setw(12) << r.energy[k];
			for (int a : at) {
				if (r.iao_pop(a, k) < 0.01) break;
				file << "   " << left << setw(5) << atom_label(wavy, a) << right << setw(8) << r.iao_pop(a, k) << " [" << r.mulliken(a, k) << "]";
			}
			file << "\n";
		}
	}
	file << endl;
}

ivec ibo_selection(const IBOResult &r, const std::string &spec)
{
	ivec sel;
	const int nibo = static_cast<int>(r.kind.size());
	std::stringstream ss(spec);
	std::string tok;
	while (std::getline(ss, tok, ',')) {
		if (tok.empty()) continue;
		if (tok == "all") {
			for (int k = 0; k < nibo; k++) sel.push_back(k);
			continue;
		}
		const size_t colon = tok.find(':');
		if (colon == std::string::npos) {
			const int k = std::stoi(tok) - 1;
			err_checkf(k >= 0 && k < nibo, "-ibo_cube: IBO " + tok + " does not exist, there are " + std::to_string(nibo), std::cout);
			sel.push_back(k);
			continue;
		}
		const int a = parse_atom(tok.substr(0, colon)), b = parse_atom(tok.substr(colon + 1));
		const ivec pair = {std::min(a, b), std::max(a, b)};
		const size_t before = sel.size();
		for (int k = 0; k < nibo; k++)
			if (r.kind[k] == IBOKind::Bond && r.centres[k] == pair) sel.push_back(k);
		err_checkf(sel.size() > before, "-ibo_cube: no bond-like IBO between atoms " + tok, std::cout);
	}
	return sel;
}

WFN ibo_wfn(const WFN &wavy, const IBOResult &r)
{
	WFN w(wavy);
	const int nocc = static_cast<int>(r.mos.size());
	for (int p = 0; p < wavy.get_nex(); p++)
		for (int k = 0; k < nocc; k++) {
			double v = 0.0;
			for (int i = 0; i < nocc; i++) v += r.U(i, k) * wavy.get_MO_coef(r.mos[i], p);
			w.set_MO_coef(r.mos[k], p, v);
		}
	return w;
}
