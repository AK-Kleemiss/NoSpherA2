#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "mo_class.h"
#include "cube.h"
#include "constants.h"
#include "fchk.h"
#include "basis_set.h"
#include "nos_math.h"
#include "libCintMain.h"
#include "integrator.h"
#include "cell.h"

bool WFN::write_wfn(const std::filesystem::path &fileName, const bool &debug, const bool occupied) const
{
	using namespace std;
	ofstream rf(fileName, ios::out);
	if (!rf.is_open())
	{
		std::cout << "Sorry, can't open the file...\n";
		return false;
	}
	if (debug)
		std::cout << "Writing " << fileName << ": ncen " << ncen << " nex " << nex << " nmo " << nmo << endl;
	rf << comment << '\n' << hdr(occupied);
	for (int i = 0; i < ncen; i++)
	{
		rf << setw(3) << atoms[i].get_label() << ' ';
		rf << left << setw(8) << i + 1 << right << "(CENTRE";
		rf << setw(3) << i + 1 << ") ";
		rf << fixed << showpoint << setprecision(8);
		rf << setw(12) << get_atom_coordinate(i, 0);
		rf << setw(12) << get_atom_coordinate(i, 1);
		rf << setw(12) << get_atom_coordinate(i, 2);
		rf << "  CHARGE = " << setw(2) << get_atom_charge(i) << ".0\n";
	}
	//Fortran fixed-width blocks: every line starts with tag and holds per values, value(i) prints the i-th
	char buf[32];
	auto write_block = [&](const string &tag, const int per, auto value)
	{
		for (int i = 0; i < nex; i++)
		{
			if (i % per == 0)
				rf << (i ? "\n" : "") << tag;
			value(i);
		}
		if (nex > 0)
			rf << '\n';
	};
	write_block("CENTRE ASSIGNMENTS  ", 20, [&](const int i) { rf << setw(3) << centers[i]; });
	write_block("TYPE ASSIGNMENTS    ", 20, [&](const int i) { rf << setw(3) << types[i]; });
	write_block("EXPONENTS ", 5, [&](const int i) { snprintf(buf, sizeof(buf), "%14.7E", exponents[i]); rf << buf; });
	int mo_run = 1;
	for (int m = 0; m < nmo; m++)
	{
		if (occupied && MOs[m].get_occ() == 0)
			continue;
		rf << "MO" << setw(3) << mo_run++ << setw(29) << "OCC NO =" << setw(13) << fixed << setprecision(8) << MOs[m].get_occ()
		   << setw(14) << "ORB. ENERGY =" << setw(13) << fixed << setprecision(8) << MOs[m].get_energy() << '\n';
		write_block("", 5, [&](const int i) { snprintf(buf, sizeof(buf), "%16.8E", MOs[m].get_coefficient(i)); rf << buf; });
	}
	rf << "END DATA\n";
	rf << " THE SCF ENERGY =" << setw(20) << fixed << setprecision(12) << total_energy << " THE VIRIAL(-V/T)=   0.00000000" << endl;
	return rf.good();
};

//AIM wfx: every quantity in its own <Tag> ... </Tag> block, the ORCA/AIMAll layout that read_wfx expects
bool WFN::write_wfx(const std::filesystem::path &fileName, const bool occupied) const
{
	using namespace std;
	ofstream rf(fileName, ios::out);
	if (!rf.is_open())
	{
		std::cout << "Sorry, can't open the file...\n";
		return false;
	}
	char buf[32];
	auto block = [&](const string &tag, auto body) { rf << '<' << tag << ">\n"; body(); rf << "</" << tag << ">\n"; };
	//per values per line; value(i) prints the i-th of n
	auto numbers = [&](const string &tag, const int n, const int per, auto value)
	{
		block(tag, [&]() { for (int i = 0; i < n; i++) { value(i); rf << ((i + 1) % per == 0 || i + 1 == n ? "\n" : " "); } });
	};
	auto sci = [&](const double x) { snprintf(buf, sizeof(buf), "%16.8E", x); rf << buf; };
	ivec sel;
	double nel = 0, nalpha = 0;
	//The spin structure follows the orbitals, not only is_unrestricted: a wavefunction assembled MO by MO can hold beta
	//orbitals without the flag set.
	bool unrestricted = is_unrestricted;
	for (int m = 0; m < nmo; m++)
		if (!(occupied && MOs[m].get_occ() == 0) && MOs[m].get_op() != 0)
			unrestricted = true;
	for (int m = 0; m < nmo; m++)
	{
		if (occupied && MOs[m].get_occ() == 0)
			continue;
		sel.push_back(m);
		nel += MOs[m].get_occ();
		nalpha += unrestricted ? (MOs[m].get_op() == 0 ? MOs[m].get_occ() : 0) : MOs[m].get_occ() / 2;
	}
	const int n = static_cast<int>(sel.size()), i_nel = static_cast<int>(round(nel)), i_nalpha = static_cast<int>(round(nalpha));
	block("Title", [&]() { rf << comment << '\n'; });
	block("Keywords", [&]() { rf << "GTO\n"; });
	block("Number of Nuclei", [&]() { rf << ncen << '\n'; });
	block("Number of Primitives", [&]() { rf << nex << '\n'; });
	block("Number of Occupied Molecular Orbitals", [&]() { rf << n << '\n'; });
	block("Number of Perturbations", [&]() { rf << "0\n"; });
	numbers("Nuclear Names", ncen, 1, [&](const int i) { rf << (atoms[i].get_label().empty() ? constants::atnr2letter(get_atom_charge(i)) + to_string(i + 1) : atoms[i].get_label()); });
	numbers("Atomic Numbers", ncen, 1, [&](const int i) { rf << get_atom_charge(i); });
	numbers("Nuclear Charges", ncen, 1, [&](const int i) { sci(get_atom_charge(i)); });
	numbers("Nuclear Cartesian Coordinates", 3 * ncen, 3, [&](const int i) { sci(get_atom_coordinate(i / 3, i % 3)); });
	block("Net Charge", [&]() { rf << charge << '\n'; });
	block("Number of Electrons", [&]() { rf << i_nel << '\n'; });
	block("Number of Alpha Electrons", [&]() { rf << i_nalpha << '\n'; });
	block("Number of Beta Electrons", [&]() { rf << i_nel - i_nalpha << '\n'; });
	block("Electronic Spin Multiplicity", [&]() { rf << (multi > 0 ? static_cast<int>(multi) : 2 * i_nalpha - i_nel + 1) << '\n'; });
	numbers("Primitive Centers", nex, 20, [&](const int i) { rf << centers[i]; });
	numbers("Primitive Types", nex, 20, [&](const int i) { rf << types[i]; });
	numbers("Primitive Exponents", nex, 5, [&](const int i) { sci(exponents[i]); });
	numbers("Molecular Orbital Occupation Numbers", n, 1, [&](const int i) { sci(MOs[sel[i]].get_occ()); });
	numbers("Molecular Orbital Energies", n, 1, [&](const int i) { sci(MOs[sel[i]].get_energy()); });
	numbers("Molecular Orbital Spin Types", n, 1, [&](const int i) { rf << (!unrestricted ? "Alpha and Beta" : MOs[sel[i]].get_op() == 0 ? "Alpha" : "Beta"); });
	block("Molecular Orbital Primitive Coefficients", [&]()
	{
		for (int i = 0; i < n; i++)
		{
			block("MO Number", [&]() { rf << i + 1 << '\n'; });
			for (int j = 0; j < nex; j++)
			{
				sci(MOs[sel[i]].get_coefficient(j));
				rf << ((j + 1) % 5 == 0 || j + 1 == nex ? "\n" : " ");
			}
		}
	});
	block("Energy = T + Vne + Vee + Vnn", [&]() { snprintf(buf, sizeof(buf), "%22.14E", total_energy); rf << buf << '\n'; });
	block("Virial Ratio (-V/T)", [&]() { sci(virial_ratio); rf << '\n'; });
	return rf.good();
};

bool WFN::write_xyz(const std::filesystem::path &fileName)
{
	using namespace std;
	ofstream f(fileName, ios::out);
	if (!f.is_open())
	{
		err("Error writing the xyz file! Aborting!", std::cout);
		return false;
	}
	f << ncen << '\n' << "XYZ File written by NoSpherA2 based on " << path << '\n';
	for (int i = 0; i < ncen; i++)
	{
		f << (atoms[i].get_label().empty() ? constants::atnr2letter(get_atom_charge(i)) : atoms[i].get_label());
		for (int c = 0; c < 3; c++)
			f << setw(14) << setprecision(8) << (isBohr ? constants::bohr2ang(get_atom_coordinate(i, c)) : get_atom_coordinate(i, c));
		f << '\n';
	}
	return f.good();
};
