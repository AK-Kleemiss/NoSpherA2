#include "pch.h"
#include <cstring>

#include "core/base/convenience.h"
#include "core/base/constants.h"
#include "core/wfn/basis_set.h"
#include "core/integration/libCintMain.h"
#include "core/salted/SALTED_utilities.h"

double NOS_CINTcommon_fac_sp(int l);

namespace
{
	// verbatim copies of master d624d8c4, frozen: the refactored functions must reproduce them bit for bit
	double NOS_CINTcommon_fac_sp_master_ref(int l)
	{
		switch (l)
		{
		case 0:
			return 0.282094791773878143;
		case 1:
			return 0.488602511902919921;
		default:
			return 1;
		}
	}

	dMatrix2 cart2sph_master_ref(const int l, const bool normalized) {
		int n_cart = (l + 1) * (l + 2) / 2;

		dMatrix2 c_tensor(n_cart, n_cart);
		for (int i = 0; i < n_cart; i++) {
			c_tensor(i, i) = 1.0;
		}

		if (l == 0 || l == 1) { //For s and p functions, the transformation is trivial
			if (l == 1) { //libcint is built with PYPZPX: spherical p comes as py, pz, px
				std::fill(c_tensor.container().begin(), c_tensor.container().end(), 0.0);
				c_tensor(1, 0) = c_tensor(2, 1) = c_tensor(0, 2) = 1.0;
			}
			if (normalized) {
				return c_tensor;
			}
			else {
				double norm_factor = (l == 0) ? 0.282094791773878143 : 0.488602511902919921;
				for (auto &val : c_tensor.container()) {
					val *= norm_factor;
				}
				return c_tensor;
			}
		}

		err_checkf(l <= 15, "cart2sph_matrix: l must be <= 15", std::cout);

		int n_sph = 2 * l + 1;
		vec c_sph(static_cast<size_t>(n_sph) * n_cart, 0.0);

		libcint::CINTc2s_ket_sph(c_sph.data(), n_cart, c_tensor.data(), l);
		//Transform back to row-major order
		dMatrix2 c_sph_RM(n_cart, n_sph);
		for (int i = 0; i < n_cart; i++) {
			for (int j = 0; j < n_sph; j++) {
				c_sph_RM(i, j) = c_sph[j * n_cart + i];
			}
		}
		return c_sph_RM;
	}
}

TEST(SaltedRiCleanupTests, CommonFactorSpIsBitwiseMaster)
{
	for (int l = -5; l <= 64; l++) {
		const double got = NOS_CINTcommon_fac_sp(l), want = NOS_CINTcommon_fac_sp_master_ref(l);
		EXPECT_EQ(std::memcmp(&got, &want, sizeof(double)), 0) << "l " << l;
	}
}

TEST(SaltedRiCleanupTests, Cart2SphIsBitwiseMaster)
{
	for (int l = 0; l <= 15; l++)
		for (int n = 0; n < 2; n++) {
			const dMatrix2 got = cart2sph(l, n == 1), want = cart2sph_master_ref(l, n == 1);
			ASSERT_EQ(got.extent(0), want.extent(0)) << "l " << l;
			ASSERT_EQ(got.extent(1), want.extent(1)) << "l " << l;
			EXPECT_EQ(std::memcmp(got.data(), want.data(), want.size() * sizeof(double)), 0) << "l " << l << " normalized " << n;
		}
}

// lmax is the highest primitive type of the species, nmax[<label><l>] its number of distinct shells of that l
TEST(SaltedRiCleanupTests, LmaxNmaxCountDistinctShellsPerL)
{
	BasisSet b;
	const SimplePrimitive h[] = { { 0, 0, 3.0, 0.5, 0 }, { 0, 0, 1.0, 0.5, 0 }, { 0, 0, 0.3, 1.0, 1 }, { 0, 1, 0.8, 1.0, 2 }, { 0, 2, 0.5, 1.0, 3 }, { 0, 2, 0.2, 1.0, 4 } };
	for (int i = 0; i < 6; i++) b.add_owned_primitive(h[i]);
	b.set_range_for_element(0, 0, 6);
	const SimplePrimitive c[] = { { 0, 0, 9.0, 1.0, 0 }, { 0, 3, 1.1, 1.0, 1 } };
	for (int i = 0; i < 2; i++) b.add_owned_primitive(c[i]);
	b.set_range_for_element(5, 6, 2);
	std::unordered_map<std::string, int> lmax, nmax;
	SALTED_Utils::set_lmax_nmax(lmax, nmax, b, { "H", "C" });
	EXPECT_EQ(lmax["H"], 2);
	EXPECT_EQ(nmax["H0"], 2);
	EXPECT_EQ(nmax["H1"], 1);
	EXPECT_EQ(nmax["H2"], 2);
	EXPECT_EQ(lmax["C"], 3);
	EXPECT_EQ(nmax["C0"], 1);
	EXPECT_EQ(nmax["C1"], 0);
	EXPECT_EQ(nmax["C2"], 0);
	EXPECT_EQ(nmax["C3"], 1);
	EXPECT_EQ(SALTED_Utils::get_lmax_max(lmax), 3);
}
