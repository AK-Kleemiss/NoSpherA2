#include "pch.h"

#include "core/base/convenience.h"
#include "core/base/nos_math.h"
#include "core/xcw/cell.h"
#include <cstring>
#include <random>
#include <ranges>

namespace NoSpherA2UnitTests
{
	namespace
	{
		std::filesystem::path cleanup_tmp(const std::string& name)
		{
			const std::string test = ::testing::UnitTest::GetInstance()->current_test_info()->name();
			return std::filesystem::temp_directory_path() / ("nosphera2_xcwcleanup_" + test + "_" + name);
		}

		cell cell_from_ops(const std::string& name, const svec& ops)
		{
			const std::filesystem::path p = cleanup_tmp(name);
			{
				std::ofstream f(p);
				f << "data_test\n_cell_length_a 8.0\n_cell_length_b 8.5\n_cell_length_c 6.0\n_cell_angle_alpha 90.0\n"
					"_cell_angle_beta 97.0\n_cell_angle_gamma 90.0\nloop_\n_space_group_symop_operation_xyz\n";
				for (int i = 0; i < ops.size(); i++) f << "'" << ops[i] << "'\n";
			}
			std::ostringstream log;
			cell cl(p, log, false, true);
			std::filesystem::remove(p);
			return cl;
		}

		const svec p21c_ops = { "x, y, z", "-x, y+1/2, -z+1/2", "-x, -y, -z", "x, -y+1/2, z+1/2" };

		asym_atom atom_at(const cell& c, const std::string& label, const int type, const d3& frac)
		{
			asym_atom a;
			a.label = label;
			a.type = type;
			const vec cart = c.get_coords_cartesian(frac[0], frac[1], frac[2], true);
			a.pos = { cart[0], cart[1], cart[2] };
			a.frac_pos = frac;
			a.cryst_mult = 1;
			a.anom = cdouble(0.0, 0.0);
			return a;
		}

		bool same_bits(const vec& a, const vec& b)
		{
			return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
		}

		// master's cell::keep_operations with the cell's sym and trans as arguments
		void keep_operations_master_ref(ivec3& sym, vec2& trans, const ivec& kept, ivec& applied_symmetry, hkl_list& hkl_enlarged, const hkl_list& hkl) {
			const int num_sym_ops = static_cast<int>(sym[0][0].size());
			applied_symmetry.clear();
			for (int sym_op = 0; sym_op < num_sym_ops; sym_op++)
				if (std::find(kept.begin(), kept.end(), sym_op) == kept.end()) {
					applied_symmetry.push_back(sym_op);
				}
			for (const int sym_op : std::views::reverse(applied_symmetry)) {
				for (ivec2& middle : sym) {
					for (ivec& inner : middle) {
						inner.erase(inner.begin() + sym_op);
					}
				}
				for (vec& inner : trans) {
					inner.erase(inner.begin() + sym_op);
				}
			}

			vec3 rotations;
			for (int sym_op = 0; sym_op < sym[0][0].size(); ++sym_op) {
				rotations.push_back({ {
						static_cast<double>(sym[0][0][sym_op]),
						static_cast<double>(sym[0][1][sym_op]),
						static_cast<double>(sym[0][2][sym_op])},{
						static_cast<double>(sym[1][0][sym_op]),
						static_cast<double>(sym[1][1][sym_op]),
						static_cast<double>(sym[1][2][sym_op])},{
						static_cast<double>(sym[2][0][sym_op]),
						static_cast<double>(sym[2][1][sym_op]),
						static_cast<double>(sym[2][2][sym_op])} });
			}

			const int nr = hkl.size();
			std::vector<i3> hkl_vec(hkl.begin(), hkl.end());
			hkl_list new_enlarged;
			for (const auto& h : hkl) {
				vec hkl_temp = { static_cast<double>(h[0]), static_cast<double>(h[1]), static_cast<double>(h[2]) };
				for (const vec2& rot : rotations) {
					vec new_hkl = self_dot(rot, hkl_temp, false);
					i3 new_hkl_int = { (int)std::round(new_hkl[0]), (int)std::round(new_hkl[1]), (int)std::round(new_hkl[2]) };
					new_enlarged.insert(new_hkl_int);
				}
			}
			hkl_enlarged = new_enlarged;
		}

		// Grows asym by xyz in cl as structure_factors does, then checks apply_grown's deletion and enlargement
		// against master's keep_operations run on the operations before it with the same kept set
		void expect_keep_operations_as_master(const svec& ops, std::vector<asym_atom> asym, const std::vector<d3>& xyz_frac, const hkl_list& hkl)
		{
			cell cl = cell_from_ops("keep.cif", ops);
			std::vector<asym_atom> xyz;
			for (int i = 0; i < xyz_frac.size(); i++) xyz.push_back(atom_at(cl, "x" + std::to_string(i), asym[0].type, xyz_frac[i]));
			cl.asym_unit_cryst_mult(asym);
			cl.grow_asym_atoms(asym, xyz);
			cl.resolve_internal_symmetry(asym);
			ivec3 sym = cl.get_sym();
			vec2 trans = cl.get_trans();
			const int nsym = static_cast<int>(sym[0][0].size());
			hkl_list enlarged;
			const ivec applied = cl.apply_grown(hkl, enlarged, asym);
			ivec kept;
			for (int s = 0; s < nsym; s++)
				if (std::find(applied.begin(), applied.end(), s) == applied.end()) kept.push_back(s);
			ivec applied_ref;
			hkl_list enlarged_ref;
			keep_operations_master_ref(sym, trans, kept, applied_ref, enlarged_ref, hkl);
			EXPECT_EQ(applied, applied_ref);
			EXPECT_EQ(enlarged, enlarged_ref);
			EXPECT_EQ(cl.get_sym(), sym);
			const vec2 t = cl.get_trans();
			ASSERT_EQ(t.size(), trans.size());
			for (int i = 0; i < t.size(); i++) EXPECT_TRUE(same_bits(t[i], trans[i]));
		}

		hkl_list test_reflections()
		{
			hkl_list hkl = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 1, 1 }, { -2, 3, -1 }, { 0, -4, 2 }, { 5, 0, -3 } };
			std::mt19937 gen(4711);
			std::uniform_int_distribution<int> d(-12, 12);
			for (int i = 0; i < 200; i++) hkl.insert({ d(gen), d(gen), d(gen) });
			return hkl;
		}
	}

	// apply_grown's kept operations, deleted operations and enlarged reflections against master's keep_operations
	// in P21/c and P-1: nothing grown, a screw image, an inversion pair about a centre, a whole unit cell
	TEST(XcwCleanupTests, KeepOperationsMatchesMaster)
	{
		const hkl_list hkl = test_reflections();
		const cell c = cell_from_ops("ref.cif", p21c_ops);
		auto at = [&](const std::string& l, const int t, const d3& f) { return atom_at(c, l, t, f); };
		expect_keep_operations_as_master(p21c_ops, { at("C1", 6, { 0.1, 0.2, 0.3 }) }, {}, hkl);
		expect_keep_operations_as_master(p21c_ops, { at("C1", 6, { 0.1, 0.2, 0.3 }) }, { { 0.1, 0.2, 0.3 }, { -0.1, 0.7, 0.2 } }, hkl);
		expect_keep_operations_as_master(p21c_ops, { at("C1", 6, { 0.1, 0.2, 0.3 }) },
			{ { 0.1, 0.2, 0.3 }, { -0.1, 0.7, 0.2 }, { -0.1, -0.2, -0.3 }, { 0.1, 0.3, 0.8 } }, hkl);
		expect_keep_operations_as_master(p21c_ops, { at("C1", 6, { 0.1, 0.2, 0.3 }) }, { { -0.1, -0.2, -0.3 } }, hkl);
		expect_keep_operations_as_master(p21c_ops, { at("Fe1", 26, { 0.0, 0.0, 0.0 }) }, { { 0.0, 0.5, 0.5 } }, hkl);
		expect_keep_operations_as_master({ "x, y, z", "-x, -y, -z" }, { at("C1", 6, { 0.1, 0.2, 0.3 }) }, { { -0.1, -0.2, -0.3 } }, hkl);
		expect_keep_operations_as_master({ "x, y, z", "-x, -y, -z" }, { at("C1", 6, { 0.1, 0.2, 0.3 }) }, {}, hkl);
	}

	// lookup[r][s] is the position of h_r R_s in the sorted enlarged list, R_s as the cell stores it (transposed)
	TEST(XcwCleanupTests, AsymLookupFindsEveryImage)
	{
		const cell cl = cell_from_ops("lookup.cif", p21c_ops);
		const ivec3 sym = cl.get_sym();
		const int nsym = static_cast<int>(cl.get_num_syms());
		const hkl_list hkl = test_reflections();
		hkl_list enlarged;
		for (const i3& h : hkl)
			for (int s = 0; s < nsym; s++) {
				i3 img = { 0, 0, 0 };
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++) img[j] += h[i] * sym[j][i][s];
				enlarged.insert(img);
			}
		const std::vector<i3> sorted(enlarged.begin(), enlarged.end());
		const ivec2 lookup = cl.generate_asym_lookup(hkl, enlarged);
		ASSERT_EQ(lookup.size(), hkl.size());
		int r = 0;
		for (const i3& h : hkl) {
			ASSERT_EQ(lookup[r].size(), static_cast<size_t>(nsym));
			for (int s = 0; s < nsym; s++) {
				const i3 img = sorted[lookup[r][s]];
				for (int j = 0; j < 3; j++)
					EXPECT_EQ(img[j], h[0] * sym[j][0][s] + h[1] * sym[j][1][s] + h[2] * sym[j][2][s]);
			}
			r++;
		}
	}
}
