#include "pch.h"
#include "nbo.h"

namespace {
	//"BD O1-H2" from the hybrids, which a parsed gennbo result carries too; the description is the fallback
	std::string nbo_label(const NboOrbital& o)
	{
		if (o.hybrids.empty()) return o.description;
		std::string s = o.type + " ";
		for (size_t k = 0; k < o.hybrids.size(); k++)
			s += (k ? "-" : "") + o.hybrids[k].element + std::to_string(o.hybrids[k].center);
		return s;
	}

	//The spins a table carries, in the order they first appear; one empty entry for a closed shell.
	template <class Row>
	svec spins_of(const std::vector<Row>& rows)
	{
		svec spins;
		for (const Row& x : rows)
			if (std::find(spins.begin(), spins.end(), x.spin) == spins.end()) spins.push_back(x.spin);
		return spins;
	}

	//"Resonance weights (alpha spin):" - the spin goes into the title, not into a column of every row
	std::string title(const std::string& what, const std::string& spin, std::string note = "")
	{
		if (!spin.empty())
			note += (note.empty() ? "" : ", ") + spin + (spin == "alpha" || spin == "beta" ? " spin" : "");
		return "\n" + what + (note.empty() ? "" : " (" + note + ")") + ":\n";
	}

	//a nearly empty orbital comes out at either sign, and "-0.00000" reads as a negative occupancy
	double printable(const double x)
	{
		return std::abs(x) < 5e-6 ? 0.0 : x;
	}
}

void print_nrt(const NboResults& r, std::ostream& out)
{
	using namespace std;
	const NboNrt& n = r.nrt;
	if (!n.present) return;
	const ostream_format_guard restore_format(out);
	//bond orders carry atom indices only; the populations table is where the elements are
	std::map<int, std::string> element;
	for (const NboAtomPopulation& p : r.npa) element[p.index] = p.element;
	for (const NboValency& v : n.valencies) element[v.atom] = v.element;
	const auto label = [&element](const int a) {
		const auto it = element.find(a);
		return (it == element.end() ? std::string("?") : it->second) + std::to_string(a);
	};
	//the parent the search starts from: dense in column blocks up to 40 atoms, one row per atom beyond
	for (const NboTopo& t : n.parent_topo) {
		const int na = static_cast<int>(t.matrix.size());
		out << title("Parent Lewis structure", t.spin);
		if (na <= 40) {
			out << "  valence lone pairs on the diagonal, bond orders off it\n";
			for (int c0 = 0; c0 < na; c0 += 12) {
				out << "\n        ";
				for (int b = c0; b < std::min(na, c0 + 12); b++) out << setw(6) << label(b + 1);
				out << "\n";
				for (int a = 0; a < na; a++) {
					out << "  " << left << setw(6) << label(a + 1) << right;
					for (int b = c0; b < std::min(na, c0 + 12); b++)
						out << setw(6) << (t.matrix[a][b] ? std::to_string(t.matrix[a][b]) : ".");
					out << "\n";
				}
			}
		}
		else {
			out << "  atom   lone pairs   partner:bond order\n";
			for (int a = 0; a < na; a++) {
				out << "  " << left << setw(7) << label(a + 1) << right << setw(10) << t.matrix[a][a] << "  ";
				for (int b = 0; b < na; b++)
					if (b != a && t.matrix[a][b]) out << " " << label(b + 1) << ":" << t.matrix[a][b];
				out << "\n";
			}
		}
	}
	if (!n.parent_topo.empty())
		out << "\nParent Lewis structure as NRT input (JSON):\n" << parent_lewis_json(r) << "\n";

	out << "\nNatural resonance theory (in house):\n"
		<< "  " << n.structures_used << " of " << n.structures_found
		<< " resonance structures carry the fit, D(0) = " << fixed << setprecision(5) << n.d_0
		<< ", D(w) = " << n.d_w << "\n";
	for (const NboNrtD& d : n.spin_d)
		out << "  in NBO's units" << (d.spin.empty() ? "" : " (" + d.spin + ")") << ": D(0) = " << d.d_0 << ", D(w) = " << d.d_w << "\n";
	//an open shell runs the search once per spin, and both runs push the same budget line
	svec seen;
	for (const std::string& s : n.notes)
		if (std::find(seen.begin(), seen.end(), s) == seen.end()) {
			seen.push_back(s);
			out << "  " << s << "\n";
		}

	for (const std::string& spin : spins_of(n.weights)) {
		out << title("Resonance weights", spin) << "  structure   weight %   change from the parent\n";
		//structure 0 is the parent, listed only when the fit drops it (a retained one is its own row)
		const bool retained = std::any_of(n.weights.begin(), n.weights.end(), [&](const NboResonanceWeight& w) {
			return w.spin == spin && w.changes == "parent";
		});
		for (const NboNrtD& d : n.spin_d)
			if (!retained && d.spin == spin && d.parent_percent >= 0.0)
				out << setw(11) << 0 << setprecision(2) << setw(11) << d.parent_percent << "   parent\n";
		for (const NboResonanceWeight& w : n.weights)
			if (w.spin == spin)
				out << setw(11) << w.structure << setprecision(2) << setw(11) << w.weight_percent << "   "
					<< w.changes << "\n";
	}
	//the bond-order diagonal is the atom's lone-pair count, so it goes into the valency table
	std::map<std::pair<std::string, int>, double> lone_pairs;
	for (const NboBondOrder& b : n.bond_orders)
		if (b.diagonal) lone_pairs[{ b.spin, b.atom1 }] = b.total;
	for (const std::string& spin : spins_of(n.bond_orders)) {
		out << title("Natural bond orders", spin) << "  atom   atom        total  covalent     ionic   range at D(w) + 0.1 %\n"
			<< setprecision(4);
		for (const NboBondOrder& b : n.bond_orders)
			if (b.spin == spin && !b.diagonal) {
				out << "  " << left << setw(7) << label(b.atom1) << setw(7) << label(b.atom2) << right
					<< setw(10) << b.total << setw(10) << b.covalent << setw(10) << b.ionic;
				if (b.range_lo >= 0.0) out << "   " << b.range_lo << " - " << b.range_hi;
				out << "\n";
			}
	}
	for (const std::string& spin : spins_of(n.valencies)) {
		out << title("Natural atomic valencies", spin)
			<< "  atom      valency  covalency  electrovalency  electrons  lone pairs\n" << setprecision(4);
		for (const NboValency& v : n.valencies) {
			if (v.spin != spin) continue;
			out << "  " << left << setw(7) << (v.element + std::to_string(v.atom)) << right << setw(10)
				<< v.valency << setw(11) << v.covalency << setw(16) << v.electrovalency << setw(11)
				<< v.electron_count;
			const auto lp = lone_pairs.find({ v.spin, v.atom });
			if (lp != lone_pairs.end()) out << setw(12) << lp->second;
			out << "\n";
		}
	}
	for (const std::string& s : n.symmetry_forms) out << "  " << s << "\n";
}

void print_nbo(const NboResults& r, std::ostream& out)
{
	using namespace std;
	//fixed/setprecision would otherwise outlive this table
	const ostream_format_guard restore_format(out);
	if (!r.npa.empty()) {
		const bool spin = r.npa.front().has_spin_density;
		out << "\nNatural population analysis (in house):\n"
			<< "  atom       charge      core   valence   Rydberg     total" << (spin ? "  spin density\n" : "\n")
			<< fixed << setprecision(5);
		double charge_sum = 0.0, total_sum = 0.0;
		for (const NboAtomPopulation& p : r.npa) {
			out << "  " << left << setw(7) << (p.element + std::to_string(p.index)) << right << setw(10)
				<< p.charge << setw(10) << p.core << setw(10) << p.valence << setw(10) << p.rydberg
				<< setw(10) << p.total;
			if (p.has_spin_density) out << setw(14) << p.spin_density;
			out << "\n";
			charge_sum += p.charge;
			total_sum += p.total;
		}
		//charges sum to the molecular charge, populations to the electron count
		out << "  " << left << setw(7) << "total" << right << setw(10) << printable(charge_sum) << setw(40) << total_sum
			<< "\n";
	}


	//diffuse or polarisation functions leave dozens of near-empty Rydberg NBOs per atom: counted, not
	//listed; the JSON keeps every one
	constexpr double rydberg_floor = 1e-4;
	const auto hybrid = [&out](const NboHybrid& h) {
		out << "   " << left << setw(6) << (h.element + std::to_string(h.center)) << right << setprecision(2)
			<< setw(9) << h.weight_percent << setw(8) << h.s << setw(8) << h.p << setw(8) << h.d << setw(8)
			<< h.f;
		if (h.sp_exponent() > 0.0) out << setw(8) << h.sp_exponent();
		out << "\n";
	};
	for (const std::string& spin : spins_of(r.orbitals)) {
		out << title("Natural bond orbitals", spin, "in house")
			<< "    nbo  type  centres       occupancy     energy   hybrid weight %     s %     p %     d %     f %    sp^x\n";
		int hidden = 0;
		double hidden_occupancy = 0.0;
		for (const NboOrbital& o : r.orbitals) {
			if (o.spin != spin) continue;
			if (o.type.rfind("RY", 0) == 0 && o.occupancy < rydberg_floor) {
				hidden++;
				hidden_occupancy += o.occupancy;
				continue;
			}
			const std::string centres = nbo_label(o).substr(o.hybrids.empty() ? 0 : o.type.size() + 1);
			out << setw(7) << o.index << "  " << left << setw(6) << o.type << setw(12) << centres << right
				<< fixed << setprecision(5) << setw(11) << printable(o.occupancy) << setw(11) << o.energy;
			if (o.hybrids.empty()) out << "\n";
			for (size_t k = 0; k < o.hybrids.size(); k++) {
				if (k) out << string(49, ' ');
				hybrid(o.hybrids[k]);
			}
		}
		if (hidden)
			out << "  " << hidden << " further Rydberg NBOs below " << scientific << setprecision(0)
				<< rydberg_floor << " e, together " << fixed << setprecision(5) << printable(hidden_occupancy)
				<< " e\n";
	}

	//"4 BD O1-H2": the number finds the orbital in the table above, the label says what it is
	std::map<std::pair<std::string, int>, const NboOrbital*> by_index;
	for (const NboOrbital& o : r.orbitals) by_index[{ o.spin, o.index }] = &o;
	const auto name = [&by_index](const std::string& spin, const int index, const std::string& fallback) {
		const auto it = by_index.find({ spin, index });
		return it == by_index.end() ? fallback : std::to_string(index) + " " + nbo_label(*it->second);
	};
	if (!r.e2_refused.empty()) out << "\n" << r.e2_refused << "\n";
	if (!r.e2.empty()) {
		for (const std::string& spin : spins_of(r.e2)) {
			out << title("Second-order donor-acceptor energies in the NBO basis", spin, "in house") << "  "
				<< left << setw(21) << "donor" << setw(21) << "acceptor" << right << setw(14) << "E(2) kcal/mol"
				<< setw(15) << "E(j)-E(i) Eh" << setw(12) << "F(i,j) Eh" << "\n";
			bool flagged = false;
			for (const NboE2Entry& e : r.e2)
				if (e.spin == spin) {
					out << "  " << left << setw(21) << name(spin, e.donor_index, e.donor) << setw(21)
						<< name(spin, e.acceptor_index, e.acceptor) << right << fixed << setprecision(2)
						<< setw(14) << e.energy_kcal << setprecision(3) << setw(15) << e.e_diff << setw(12)
						<< e.fij << (e.beyond_second_order ? " !" : "") << "\n";
					flagged |= e.beyond_second_order;
				}
			if (flagged)
				out << "  ! |F(i,j)| >= E(j)-E(i): beyond second order, the exact two-level mixing gains at most q|F|\n";
		}
	}
	for (const std::string& spin : spins_of(r.budget)) {
		const double c = spin.empty() ? 2.0 : 1.0;
		const bool priced = !r.transfer_fock.empty();
		out << title(priced ? "Donor -> acceptor occupancies of the density, density-derived E2(D)"
							: "Donor -> acceptor occupancies of the density, no Fock matrix", spin, "in house")
			<< "  q(k->j) = D(j,k)^2/(" << c << " - D(j,j)) in the NBO basis, largest first, q >= 1e-4 e";
		if (priced)
			out << ";\n  E2(D) = q (F(j,j) - F(k,k)) on the NBO-basis Fock diagonal of " << r.transfer_fock
				<< ", E(2) to second order for one determinant";
		out << ":\n  " << left << setw(21) << "donor" << setw(21) << "acceptor" << right << setw(11) << "q / e";
		if (priced) out << setw(16) << "E2(D) kcal/mol";
		out << "\n";
		std::vector<NboTransfer> sorted;
		for (const NboTransfer& t : r.transfer)
			if (t.spin == spin) sorted.push_back(t);
		std::stable_sort(sorted.begin(), sorted.end(), [](const NboTransfer& a, const NboTransfer& b) { return a.q > b.q; });
		for (const NboTransfer& t : sorted) {
			out << "  " << left << setw(21) << name(spin, t.donor_index, t.donor) << setw(21)
				<< name(spin, t.acceptor_index, t.acceptor) << right << fixed << setprecision(5) << setw(11) << t.q;
			if (priced) out << setprecision(2) << setw(16) << t.e2d_kcal;
			out << "\n";
		}
		//D(j,j) = q(Lewis) + q(non-Lewis) + s; s, the part no donor explains, is 0 for one determinant
		out << "  Acceptor budget, occupancy = q from Lewis + q from non-Lewis + s:\n  " << left << setw(21) << "acceptor"
			<< right << setw(11) << "occupancy" << setw(11) << "q(L)" << setw(11) << "q(NL)" << setw(11) << "s" << "\n";
		int hidden = 0;
		double sum[4] = { 0.0, 0.0, 0.0, 0.0 };
		for (const NboAcceptorBudget& b : r.budget) {
			if (b.spin != spin) continue;
			sum[0] += b.occupancy, sum[1] += b.q_lewis, sum[2] += b.q_nonlewis, sum[3] += b.s;
			if (b.occupancy < rydberg_floor) {
				hidden++;
				continue;
			}
			out << "  " << left << setw(21) << name(spin, b.index, b.label) << right << fixed << setprecision(5)
				<< setw(11) << b.occupancy << setw(11) << b.q_lewis << setw(11) << b.q_nonlewis << setw(11) << printable(b.s) << "\n";
		}
		out << "  " << left << setw(21) << ("total" + (hidden ? " (" + std::to_string(hidden) + " < 1e-4 e unlisted)" : std::string())) << right
			<< setw(11) << sum[0] << setw(11) << sum[1] << setw(11) << sum[2] << setw(11) << printable(sum[3]) << "\n";
		bool any = false;
		for (const NboLocalNoon& l : r.local_noon) {
			if (l.spin != spin) continue;
			if (!any) out << "  Local natural occupations of each BD/BD* pair, b_eff = (n+ - n-)/2:\n  " << left << setw(21) << "bond"
						  << right << setw(11) << "n+" << setw(11) << "n-" << setw(11) << "b_eff" << "\n";
			any = true;
			out << "  " << left << setw(21) << name(spin, l.bond_index, l.label) << right << fixed << setprecision(5)
				<< setw(11) << l.n_plus << setw(11) << l.n_minus << setw(11) << l.b_eff << "\n";
		}
	}
	if (!r.e2_avg.structures.empty()) {
		const NboE2Average& a = r.e2_avg;
		out << "\nSecond-order E(2) averaged over NRT structures " << e2_avg_window(a)
			<< " (in house; channels sum orbitals of one type on the same centres, weights renormalised):\n";
		for (size_t i = 0; i < a.structures.size(); i++)
			out << "  structure " << setw(3) << a.structures[i] << "  weight " << setprecision(4) << setw(7) << a.weights[i]
				<< "  rho(NL) " << setprecision(5) << a.rho_nl[i] << "\n";
		out << "  " << left << setw(21) << "donor" << setw(21) << "acceptor" << right << setw(10) << "<E(2)>"
			<< setw(9) << "min" << setw(9) << "max" << setw(9) << "present" << setw(9) << "<q>/e" << "\n";
		for (const NboE2Channel& c : a.channels)
			out << "  " << left << setw(21) << c.donor << setw(21) << c.acceptor << right << setprecision(2)
				<< setw(10) << c.mean_kcal << setw(9) << c.min_kcal << setw(9) << c.max_kcal << setw(9) << c.present
				<< setprecision(4) << setw(9) << c.mean_e << "\n";
	}
	print_nrt(r, out);
}
