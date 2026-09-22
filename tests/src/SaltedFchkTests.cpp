#include "pch.h"
#include "core/convenience.h"
#include "core/constants.h"
#include "core/fchk.h"
#include "core/SALTED_io.h"
#include "core/SALTED_utilities.h"
#include "core/SALTED_equicomb.h"
#include "core/SALTED_predictor.h"
#include "core/wfn_class.h"
#include "core/basis_set.h"
#include "core/cube.h"
#include "core/npy.h"
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#undef I

namespace
{
	std::filesystem::path tmp_path(const std::string& name)
	{
		return std::filesystem::temp_directory_path() / ("nosphera2_saltedfchk_" + name);
	}

	// byte builder for the .salted container: every reader walks these exact byte layouts
	struct salted_writer
	{
		std::string buf;
		template <class T>
		void raw(const T& v)
		{
			buf.append(reinterpret_cast<const char*>(&v), sizeof(T));
		}
		void tag(std::string s)
		{
			s.resize(5, ' ');
			buf += s;
		}
		template <class T>
		void dataset(const std::vector<T>& d, const std::vector<uint32_t>& dims)
		{
			raw(static_cast<int32_t>(dims.size()));
			for (const uint32_t x : dims)
				raw(x);
			buf.append(reinterpret_cast<const char*>(d.data()), d.size() * sizeof(T));
		}
		// config fields carry a 9-byte header the reader steps over
		template <class T>
		void field(const T& v)
		{
			buf.append(9, 'x');
			raw(v);
		}
		void str(const std::string& s)
		{
			buf.append(9, 'x');
			raw(static_cast<int32_t>(s.size()));
			buf += s;
		}
		void block_head(const int n)
		{
			raw(static_cast<int32_t>(1));
			raw(static_cast<int32_t>(n));
		}
	};

	// the synthetic model: two species, two lambdas for H and one for O
	const vec feats_h0{ 1, 2, 3, 4, 5, 6 }, feats_h1{ 7, 8 }, feats_o0{ 9, 10, 11, 12, 13, 14 };
	const vec weights{ 0.1, 0.2, 0.3, 0.4 };

	void lambda_block(salted_writer& w, const double scale)
	{
		w.block_head(2);
		w.tag("H");
		w.raw(static_cast<int32_t>(2));
		vec a(feats_h0), b(feats_h1), c(feats_o0);
		for (auto* v : { &a, &b, &c })
			for (double& x : *v)
				x *= scale;
		w.dataset(a, { 2, 3 });
		w.dataset(b, { 1, 2 });
		w.tag("O");
		w.raw(static_cast<int32_t>(1));
		w.dataset(c, { 3, 2 });
	}

	void write_synthetic_model(const std::filesystem::path& path, const int version, const bool with_basis, const bool with_normc)
	{
		std::vector<std::pair<std::string, std::string>> blocks;
		{
			salted_writer w;
			w.field(static_cast<char>(1));
			w.field(static_cast<char>(0));
			w.field(static_cast<char>(1));
			for (const int32_t v : { 4, 2, 2, 1, 1, 3, 5 })
				w.field(v);
			for (const double v : { 3.5, 3.5, 0.3, 0.4, 2.0, 0.8 })
				w.field(v);
			w.str("H O");
			w.str("H O");
			w.str("H O");
			w.str("tiny-jkfit");
			blocks.emplace_back("CONFG", w.buf);
		}
		{
			salted_writer w;
			w.block_head(2);
			w.tag("H");
			w.dataset(vec{ 0.5 }, { 1 });
			w.tag("O");
			w.dataset(vec{ 1.0, 2.0 }, { 2 });
			blocks.emplace_back("AVERG", w.buf);
		}
		{
			salted_writer w;
			w.block_head(2);
			w.dataset(vec{ 1.0, 2.0, 3.0 }, { 3 });
			w.dataset(vec{ 4.0 }, { 1 });
			blocks.emplace_back("WIG", w.buf);
		}
		{
			salted_writer w;
			w.block_head(2);
			w.dataset(std::vector<int64_t>{ 0, 3 }, { 2 });
			w.dataset(std::vector<int64_t>{ 1 }, { 1 });
			blocks.emplace_back("FPS", w.buf);
		}
		{
			salted_writer w;
			w.block_head(1);
			w.dataset(weights, { 4 });
			blocks.emplace_back("WEIGH", w.buf);
		}
		{
			salted_writer w;
			lambda_block(w, 1.0);
			blocks.emplace_back("FEATS", w.buf);
		}
		{
			salted_writer w;
			lambda_block(w, 10.0);
			blocks.emplace_back("PROJ", w.buf);
		}
		if (with_normc)
		{
			salted_writer w;
			w.block_head(3);
			w.tag("MODE");
			w.dataset(vec{ 1.0 }, { 1 });
			w.tag("DEFCT");
			w.dataset(vec{ -0.00235 }, { 1 });
			w.tag("NCAL");
			w.dataset(vec{ 600.0 }, { 1 });
			blocks.emplace_back("NORMC", w.buf);
		}
		if (with_basis)
		{
			salted_writer w;
			w.block_head(3);
			w.raw(static_cast<int32_t>(1));
			w.dataset(std::vector<int32_t>{ 1 }, { 1 });
			w.dataset(std::vector<int32_t>{ 0 }, { 1 });
			w.dataset(vec{ 1.5 }, { 1 });
			w.dataset(vec{ 1.0 }, { 1 });
			w.raw(static_cast<int32_t>(8));
			w.dataset(std::vector<int32_t>{ 1, 1 }, { 2 });
			w.dataset(std::vector<int32_t>{ 0, 1 }, { 2 });
			w.dataset(vec{ 2.0, 0.7 }, { 2 });
			w.dataset(vec{ 1.0, 1.0 }, { 2 });
			// carbon carries one CONTRACTED shell: two primitives sharing one angular momentum
			w.raw(static_cast<int32_t>(6));
			w.dataset(std::vector<int32_t>{ 2 }, { 1 });
			w.dataset(std::vector<int32_t>{ 1 }, { 1 });
			w.dataset(vec{ 3.0, 0.9 }, { 2 });
			w.dataset(vec{ 0.6, 0.4 }, { 2 });
			blocks.emplace_back("BASIS", w.buf);
		}
		salted_writer h;
		h.buf = "SALTD";
		h.raw(static_cast<int32_t>(version));
		h.raw(static_cast<int32_t>(blocks.size()));
		int32_t loc = static_cast<int32_t>(5 + 4 + 4 + 9 * blocks.size());
		for (const auto& b : blocks)
		{
			h.tag(b.first);
			h.raw(loc);
			loc += static_cast<int32_t>(b.second.size());
		}
		std::ofstream out(path, std::ios::binary);
		out << h.buf;
		for (const auto& b : blocks)
			out << b.second;
	}

	// an atom carrying a decontracted auxiliary basis, 0-based types as the RI code stores them
	atom aux_atom(const std::string& label, const int Z, const double x, const double y, const double z, const std::vector<std::pair<double, int>>& shells)
	{
		atom A(label, {}, 1, x, y, z, Z);
		std::vector<unsigned int> counts;
		for (size_t s = 0; s < shells.size(); s++)
		{
			A.push_back_basis_set(shells[s].first, 1.0, shells[s].second, static_cast<int>(s));
			counts.push_back(1u);
		}
		A.set_shellcount(counts);
		return A;
	}

	// the aux function of exponent a and angular momentum l is N r^l exp(-a r^2) Y_lm with the radial norm
	// N = 1 / sqrt(int r^(2l+2) exp(-2 a r^2) dr) = sqrt(2 (2a)^(l+3/2) / Gamma(l+3/2)), hand-derived so the
	// checks below do not lean on Int_Params::normalize_gto
	double radial_norm(const double a, const int l)
	{
		return std::sqrt(2.0 * std::pow(2.0 * a, l + 1.5) / std::tgamma(l + 1.5));
	}

	// electrons one unit of s coefficient carries: N Y00 int exp(-a r^2) d3r = N pi / (2 a^(3/2)) = (2 pi / a)^(3/4)
	double electrons_per_unit(const double a)
	{
		return std::pow(2.0 * constants::PI / a, 0.75);
	}

	// H2+ doublet from occ with one s and one p primitive per centre: the pure and the cartesian counts agree,
	// which free_fchk's beta branch relies on
	WFN occ_h2_plus()
	{
		spdlog::set_level(spdlog::level::err);
		const std::vector<occ::core::Atom> atoms{ { 1, 0.0, 0.0, -0.7 }, { 1, 0.0, 0.0, 0.7 } };
		std::vector<occ::gto::Shell> shells;
		for (const auto& at : atoms)
			for (int l = 0; l <= 1; l++)
			{
				shells.emplace_back(l, vec{ l == 0 ? 1.2 : 0.9 }, vec2{ { 1.0 } }, std::array<double, 3>{ at.x, at.y, at.z });
				shells.back().kind = occ::gto::Shell::Kind::Spherical;
				shells.back().incorporate_shell_norm();
			}
		occ::gto::AOBasis basis(atoms, shells, "sp");
		basis.set_pure(true);
		occ::qm::HartreeFock hf(basis);
		occ::qm::SCF<occ::qm::HartreeFock> scf(hf, occ::qm::SpinorbitalKind::Unrestricted);
		scf.set_charge_multiplicity(1, 2);
		scf.compute_initial_guess();
		scf.compute_scf_energy();
		WFN w(scf.wavefunction(), false);
		w.set_multi(2);
		return w;
	}

	void expect_same_density(const WFN& a, const WFN& b, const double rtol, const std::string& what)
	{
		ASSERT_EQ(a.get_ncen(), b.get_ncen()) << what;
		for (int i = 0; i < a.get_ncen(); i++)
			for (const double r : { 0.3, 1.1, 2.5 })
			{
				const d3 pos{ a.get_atom_coordinate(i, 0) + r, a.get_atom_coordinate(i, 1) + 0.5 * r, a.get_atom_coordinate(i, 2) - 0.7 * r };
				const double da = a.compute_dens(pos);
				EXPECT_NEAR(da, b.compute_dens(pos), rtol * std::max(1.0, da)) << what << " atom " << i << " r " << r;
			}
	}

	std::string second_line(const std::filesystem::path& p)
	{
		std::ifstream in(p);
		std::string line;
		std::getline(in, line);
		std::getline(in, line);
		return line;
	}

	bool full_tests_enabled()
	{
		const char* env = std::getenv("RUN_FULL_TEST");
		return env && std::string(env) != "0" && std::string(env) != "false";
	}
}

// ---------------------------------------------------------------- SALTED_io

// the text config parser: booleans, numbers, quoted species lists and the derived nspe counts
TEST(SaltedFchkIoTests, ConfigTextFileRoundTrip)
{
	const auto p = tmp_path("config.txt");
	{
		std::ofstream out(p);
		out << "average = True\nfield = False\nsparsify = True\nncut = 12\n"
			<< "species = [ 'H', \"C\", O ]\nrcut1 = 4.0\nrcut2 = 3.5\nnang1 = 6\nnang2 = 5\nnrad1 = 7\nnrad2 = 8\n"
			<< "sig1 = 0.3\nsig2 = 0.25\nneighspe1 = ['H','C']\nneighspe2 = [O]\nzeta = 2.0\nMenv = 100\nNtrain = 50\n"
			<< "trainfrac = 0.8\ndfbasis = cc-pvqz-jkfit\nunknown = ignored\n";
	}
	SALTEDConfig c{};
	c.populateFromFile(p);
	std::filesystem::remove(p);
	EXPECT_TRUE(c.average);
	EXPECT_FALSE(c.field);
	EXPECT_TRUE(c.sparsify);
	EXPECT_EQ(c.ncut, 12);
	ASSERT_EQ(c.species.size(), 3u);
	EXPECT_EQ(c.species[0], "H");
	EXPECT_EQ(c.species[1], "C");
	EXPECT_EQ(c.species[2], "O");
	EXPECT_NEAR(c.rcut1, 4.0, 1e-12);
	EXPECT_NEAR(c.rcut2, 3.5, 1e-12);
	EXPECT_EQ(c.nang1, 6);
	EXPECT_EQ(c.nang2, 5);
	EXPECT_EQ(c.nrad1, 7);
	EXPECT_EQ(c.nrad2, 8);
	EXPECT_NEAR(c.sig1, 0.3, 1e-12);
	EXPECT_NEAR(c.sig2, 0.25, 1e-12);
	EXPECT_EQ(c.nspe1, 2);
	EXPECT_EQ(c.nspe2, 1);
	EXPECT_EQ(c.neighspe2[0], "O");
	EXPECT_NEAR(c.zeta, 2.0, 1e-12);
	EXPECT_EQ(c.Menv, 100);
	EXPECT_EQ(c.Ntrain, 50);
	EXPECT_NEAR(c.trainfrac, 0.8, 1e-12);
	EXPECT_EQ(c.dfbasis, "cc-pvqz-jkfit");
	EXPECT_FALSE(c.from_binary);
}

// one value per line, a line that does not parse is dropped rather than aborting
TEST(SaltedFchkIoTests, ReadVectorFromFileSkipsUnparsableLines)
{
	const auto p = tmp_path("vector.txt");
	{
		std::ofstream out(p);
		out << "1.5\nabc\n-2.25\n\n7\n";
	}
	const vec d = readVectorFromFile<double>(p);
	const ivec i = readVectorFromFile<int>(p);
	std::filesystem::remove(p);
	ASSERT_EQ(d.size(), 3u);
	EXPECT_NEAR(d[0], 1.5, 1e-15);
	EXPECT_NEAR(d[1], -2.25, 1e-15);
	EXPECT_NEAR(d[2], 7.0, 1e-15);
	ASSERT_EQ(i.size(), 3u);
	EXPECT_EQ(i[0], 1);
	EXPECT_EQ(i[1], -2);
	EXPECT_EQ(i[2], 7);
}

// a missing vector file is an err_checkf exit, not an empty vector
TEST(SaltedFchkIoTests, ReadVectorFromFileMissingExits)
{
	EXPECT_EXIT(readVectorFromFile<double>(tmp_path("does_not_exist.txt")), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

// read_fps appends the lambda and .npy to the prefix, one int64 array per lambda; read_npy is the double path
TEST(SaltedFchkIoTests, ReadFpsFromNpyPerLambda)
{
	const std::string prefix = tmp_path("fps").string();
	for (int lam = 0; lam < 3; lam++)
	{
		npy::npy_data<int64_t> d;
		d.data = { static_cast<int64_t>(lam), static_cast<int64_t>(10 * lam + 1) };
		d.shape = { 2 };
		npy::write_npy(prefix + std::to_string(lam) + ".npy", d);
	}
	std::filesystem::path p(prefix);
	const auto fps = read_fps<int64_t>(p, 2);
	ASSERT_EQ(fps.size(), 3u);
	for (int lam = 0; lam < 3; lam++)
	{
		ASSERT_EQ(fps.at(lam).size(), 2u);
		EXPECT_EQ(fps.at(lam)[0], lam);
		EXPECT_EQ(fps.at(lam)[1], 10 * lam + 1);
		std::filesystem::remove(prefix + std::to_string(lam) + ".npy");
	}
	npy::npy_data<double> dd;
	dd.data = { 0.25, -1.0, 3.5 };
	dd.shape = { 3 };
	std::filesystem::path dp = tmp_path("doubles.npy");
	npy::write_npy(dp.string(), dd);
	vec back;
	read_npy<double>(dp, back);
	std::filesystem::remove(dp);
	ASSERT_EQ(back.size(), 3u);
	EXPECT_NEAR(back[1], -1.0, 1e-15);
	EXPECT_NEAR(back[2], 3.5, 1e-15);
}

// the model directory scan returns the bare file name, an empty path without a model and survives a bad directory
TEST(SaltedFchkIoTests, FindFirstSaltedFile)
{
	const auto dir = tmp_path("modeldir");
	std::filesystem::create_directories(dir);
	{
		std::ofstream(dir / "notes.txt") << "x";
	}
	EXPECT_TRUE(find_first_salted_file(dir).empty());
	{
		std::ofstream(dir / "m.salted") << "x";
	}
	EXPECT_EQ(find_first_salted_file(dir).string(), "m.salted");
	std::filesystem::remove_all(dir);
	EXPECT_TRUE(find_first_salted_file(dir).empty());
}

// the CONFG block of the synthetic model reads back field for field
TEST(SaltedFchkIoTests, SyntheticModelConfig)
{
	const auto p = tmp_path("cfg.salted");
	write_synthetic_model(p, 3, false, false);
	SALTEDConfig c{};
	{
		SALTED_BINARY_FILE f(p);
		f.populate_config(c);
		EXPECT_FALSE(f.basis_set_defined());
		EXPECT_FALSE(f.charge_constraint_defined());
	}
	std::filesystem::remove(p);
	EXPECT_TRUE(c.average);
	EXPECT_FALSE(c.field);
	EXPECT_TRUE(c.sparsify);
	EXPECT_EQ(c.ncut, 4);
	EXPECT_EQ(c.nang1, 2);
	EXPECT_EQ(c.nang2, 2);
	EXPECT_EQ(c.nrad1, 1);
	EXPECT_EQ(c.nrad2, 1);
	EXPECT_EQ(c.Menv, 3);
	EXPECT_EQ(c.Ntrain, 5);
	EXPECT_NEAR(c.rcut1, 3.5, 1e-12);
	EXPECT_NEAR(c.sig2, 0.4, 1e-12);
	EXPECT_NEAR(c.zeta, 2.0, 1e-12);
	EXPECT_NEAR(c.trainfrac, 0.8, 1e-12);
	ASSERT_EQ(c.species.size(), 2u);
	EXPECT_EQ(c.species[1], "O");
	EXPECT_EQ(c.nspe1, 2);
	EXPECT_EQ(c.nspe2, 2);
	EXPECT_EQ(c.dfbasis, "tiny-jkfit");
}

// averages, wigners, fps and weights come back with the written keys, sizes and values
TEST(SaltedFchkIoTests, SyntheticModelSimpleBlocks)
{
	const auto p = tmp_path("simple.salted");
	write_synthetic_model(p, 3, false, false);
	{
		SALTED_BINARY_FILE f(p);
		const auto av = f.read_averages();
		ASSERT_EQ(av.size(), 2u);
		ASSERT_EQ(av.at("O").size(), 2u);
		EXPECT_NEAR(av.at("H")[0], 0.5, 1e-15);
		EXPECT_NEAR(av.at("O")[1], 2.0, 1e-15);
		const auto wig = f.read_wigners();
		ASSERT_EQ(wig.size(), 2u);
		ASSERT_EQ(wig.at(0).size(), 3u);
		EXPECT_NEAR(wig.at(0)[2], 3.0, 1e-15);
		EXPECT_NEAR(wig.at(1)[0], 4.0, 1e-15);
		const auto fps = f.read_fps();
		ASSERT_EQ(fps.size(), 2u);
		ASSERT_EQ(fps.at(0).size(), 2u);
		EXPECT_EQ(fps.at(0)[1], 3);
		EXPECT_EQ(fps.at(1)[0], 1);
		const vec w = f.read_weights();
		ASSERT_EQ(w.size(), 4u);
		EXPECT_NEAR(w[3], 0.4, 1e-15);
	}
	std::filesystem::remove(p);
}

// the VERSION 3 NORMC block: presence flag and the three keyed entries
TEST(SaltedFchkIoTests, SyntheticModelChargeConstraint)
{
	const auto p = tmp_path("normc.salted");
	write_synthetic_model(p, 3, false, true);
	{
		SALTED_BINARY_FILE f(p);
		ASSERT_TRUE(f.charge_constraint_defined());
		const auto e = f.read_charge_constraint();
		ASSERT_EQ(e.size(), 3u);
		EXPECT_EQ(std::lround(e.at("MODE")[0]), 1);
		EXPECT_NEAR(e.at("DEFCT")[0], -0.00235, 1e-15);
		EXPECT_NEAR(e.at("NCAL")[0], 600.0, 1e-15);
	}
	std::filesystem::remove(p);
}

// the BASIS block becomes owned primitives with per-element ranges; a species not in the block is absent
TEST(SaltedFchkIoTests, SyntheticModelBasisSet)
{
	const auto p = tmp_path("basis.salted");
	write_synthetic_model(p, 3, true, false);
	std::shared_ptr<BasisSet> b;
	{
		SALTED_BINARY_FILE f(p);
		ASSERT_TRUE(f.basis_set_defined());
		b = f.read_basis_set();
	}
	std::filesystem::remove(p);
	ASSERT_TRUE(b);
	EXPECT_EQ(b->get_owned_primitive_count(), 5u);
	EXPECT_EQ(b->get_primitive_count(), 5u);
	EXPECT_TRUE(b->has_element(1));
	EXPECT_TRUE(b->has_element(8));
	EXPECT_TRUE(b->has_element(6));
	EXPECT_FALSE(b->has_element(7));
	const auto h = (*b)[0];
	ASSERT_EQ(h.size(), 1u);
	EXPECT_NEAR(h[0].exp, 1.5, 1e-15);
	EXPECT_EQ(h[0].type, 0);
	const auto o = (*b)[7];
	ASSERT_EQ(o.size(), 2u);
	EXPECT_NEAR(o[1].exp, 0.7, 1e-15);
	EXPECT_EQ(o[1].type, 1);
	EXPECT_EQ(o[1].shell, 1);
	EXPECT_EQ(o[0].shell, 0);
	// both primitives of the contracted shell keep the shell's angular momentum
	const auto c = (*b)[5];
	ASSERT_EQ(c.size(), 2u);
	EXPECT_EQ(c[0].type, 1);
	EXPECT_EQ(c[1].type, 1);
	EXPECT_EQ(c[0].shell, 0);
	EXPECT_EQ(c[1].shell, 0);
	EXPECT_NEAR(c[1].exp, 0.9, 1e-15);
	EXPECT_NEAR(c[1].coefficient, 0.4, 1e-15);
}

// wanted species are loaded, the rest contribute only their shape; features load everything and keep row-major order
TEST(SaltedFchkIoTests, SyntheticModelProjectorsWanted)
{
	const auto p = tmp_path("proj.salted");
	write_synthetic_model(p, 3, false, false);
	{
		SALTED_BINARY_FILE f(p);
		const std::unordered_set<std::string> wanted{ "H" };
		std::unordered_map<std::string, std::array<size_t, 2>> dims;
		const auto proj = f.read_projectors(&wanted, &dims);
		ASSERT_EQ(proj.size(), 2u);
		EXPECT_TRUE(proj.count("H0"));
		EXPECT_TRUE(proj.count("H1"));
		EXPECT_FALSE(proj.count("O0"));
		ASSERT_EQ(dims.size(), 3u);
		EXPECT_EQ(dims.at("O0")[0], 3u);
		EXPECT_EQ(dims.at("O0")[1], 2u);
		EXPECT_EQ(dims.at("H1")[1], 2u);
		EXPECT_NEAR(proj.at("H0")(1, 2), 60.0, 1e-15);
		EXPECT_NEAR(proj.at("H1")(0, 1), 80.0, 1e-15);
		const auto feats = f.read_features();
		ASSERT_EQ(feats.size(), 3u);
		EXPECT_EQ(feats.at("O0").extent(0), 3u);
		EXPECT_NEAR(feats.at("O0")(2, 1), 14.0, 1e-15);
		EXPECT_NEAR(feats.at("H0")(0, 0), 1.0, 1e-15);
	}
	std::filesystem::remove(p);
}

// the lazy index records offsets and shapes only; loading a block through it equals the eager read
TEST(SaltedFchkIoTests, SyntheticModelIndexAndLoadBlock)
{
	const auto p = tmp_path("index.salted");
	write_synthetic_model(p, 3, false, false);
	{
		SALTED_BINARY_FILE f(p);
		const auto idx = f.index_lambda_based_data("FEATS");
		ASSERT_EQ(idx.size(), 3u);
		EXPECT_EQ(idx.at("H0").rows, 2u);
		EXPECT_EQ(idx.at("H0").cols, 3u);
		EXPECT_EQ(idx.at("O0").rows, 3u);
		const auto eager = f.read_features();
		for (const auto& [key, ref] : idx)
		{
			const dMatrix2 lazy = f.load_block(ref);
			const dMatrix2& e = eager.at(key);
			ASSERT_EQ(lazy.extent(0), e.extent(0)) << key;
			ASSERT_EQ(lazy.extent(1), e.extent(1)) << key;
			for (size_t i = 0; i < e.extent(0); i++)
				for (size_t j = 0; j < e.extent(1); j++)
					EXPECT_NEAR(lazy(i, j), e(i, j), 1e-15) << key;
		}
		const SALTED_BINARY_FILE::block_ref empty{};
		EXPECT_EQ(f.load_block(empty).extent(0), 0u);
	}
	std::filesystem::remove(p);
}

// a file from the future warns but still reads through its table of contents
TEST(SaltedFchkIoTests, NewerVersionStillReads)
{
	const auto p = tmp_path("future.salted");
	write_synthetic_model(p, 4, false, true);
	{
		SALTED_BINARY_FILE f(p);
		EXPECT_TRUE(f.charge_constraint_defined());
		EXPECT_EQ(f.read_weights().size(), 4u);
	}
	std::filesystem::remove(p);
}

// a wrong magic number, a negative block count and a truncated table of contents all abort in the constructor
TEST(SaltedFchkIoTests, CorruptHeaderExits)
{
	const auto bad_magic = tmp_path("badmagic.salted"), neg = tmp_path("negblocks.salted"), trunc = tmp_path("trunc.salted");
	{
		salted_writer w;
		w.buf = "SALTX";
		w.raw(static_cast<int32_t>(3));
		w.raw(static_cast<int32_t>(0));
		std::ofstream(bad_magic, std::ios::binary) << w.buf;
	}
	{
		salted_writer w;
		w.buf = "SALTD";
		w.raw(static_cast<int32_t>(3));
		w.raw(static_cast<int32_t>(-1));
		std::ofstream(neg, std::ios::binary) << w.buf;
	}
	{
		salted_writer w;
		w.buf = "SALTD";
		w.raw(static_cast<int32_t>(3));
		w.raw(static_cast<int32_t>(2));
		w.tag("CONFG");
		std::ofstream(trunc, std::ios::binary) << w.buf;
	}
	EXPECT_EXIT(SALTED_BINARY_FILE f(bad_magic), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(SALTED_BINARY_FILE f(neg), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(SALTED_BINARY_FILE f(trunc), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(SALTED_BINARY_FILE f(tmp_path("absent.salted")), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	std::filesystem::remove(bad_magic);
	std::filesystem::remove(neg);
	std::filesystem::remove(trunc);
}

// a model that stopped copying part-way keeps a valid header listing blocks that are
// no longer in the file; it has to say so instead of failing inside the first block read
TEST(SaltedFchkIoTests, TruncatedFileExits)
{
	const auto p = tmp_path("truncated.salted");
	write_synthetic_model(p, 3, true, true);
	const auto full = std::filesystem::file_size(p);
	std::filesystem::resize_file(p, full / 2);
	EXPECT_EXIT(SALTED_BINARY_FILE f(p), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), "is incomplete");
	std::filesystem::remove(p);
}

// asking for a block the table of contents does not list is fatal
TEST(SaltedFchkIoTests, MissingBlockExits)
{
	const auto p = tmp_path("nonormc.salted");
	write_synthetic_model(p, 3, false, false);
	EXPECT_EXIT({ SALTED_BINARY_FILE f(p); f.read_charge_constraint(); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT({ SALTED_BINARY_FILE f(p); f.read_basis_set(); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	std::filesystem::remove(p);
}

// the shipped model fixture: config and block sizes as decoded from the file
TEST(SaltedFchkIoTests, ShippedModelHeader)
{
	const auto p = nos_test_repo_root() / "tests" / "SALTED" / "Model" / "model.salted";
	SALTED_BINARY_FILE f(p);
	SALTEDConfig c{};
	f.populate_config(c);
	EXPECT_EQ(c.dfbasis, "cc-pvqz-jkfit");
	ASSERT_EQ(c.species.size(), 5u);
	EXPECT_EQ(c.species[4], "S");
	EXPECT_EQ(c.ncut, 500);
	EXPECT_EQ(c.Menv, 150);
	EXPECT_EQ(c.nrad1, 6);
	EXPECT_EQ(c.nang1, 7);
	EXPECT_NEAR(c.rcut1, 4.0, 1e-12);
	EXPECT_FALSE(f.basis_set_defined());
	EXPECT_FALSE(f.charge_constraint_defined());
	EXPECT_EQ(f.read_weights().size(), 10176u);
	EXPECT_EQ(f.read_fps().at(0).size(), 500u);
	EXPECT_EQ(f.read_averages().at("S").size(), 13u);
	EXPECT_EQ(f.read_wigners().at(0).size(), 64u);
	const auto idx = f.index_lambda_based_data("FEATS");
	EXPECT_EQ(idx.size(), 29u); //H lambda 0..4, C N O S lambda 0..5
	EXPECT_EQ(idx.at("C0").rows, 17u);
	EXPECT_EQ(idx.at("C0").cols, 500u);
}

// ---------------------------------------------------------- SALTED_utilities

// the complex-to-real matrices are unitary and the l = 1 one has the documented entries
TEST(SaltedFchkUtilTests, ComplexToRealIsUnitary)
{
	const auto mats = SALTED_Utils::complex_to_real_transformation({ 1, 3, 5, 7 });
	ASSERT_EQ(mats.size(), 4u);
	for (const auto& m : mats)
	{
		const size_t n = m.size();
		for (size_t i = 0; i < n; i++)
			for (size_t j = 0; j < n; j++)
			{
				cdouble s = 0.0;
				for (size_t k = 0; k < n; k++)
					s += m[i][k] * std::conj(m[j][k]);
				EXPECT_NEAR(s.real(), i == j ? 1.0 : 0.0, 1e-14) << n << " " << i << " " << j;
				EXPECT_NEAR(s.imag(), 0.0, 1e-14) << n << " " << i << " " << j;
			}
	}
	const double r = 1.0 / std::sqrt(2.0);
	EXPECT_NEAR(mats[1][0][0].imag(), r, 1e-15);
	EXPECT_NEAR(mats[1][0][2].imag(), r, 1e-15);
	EXPECT_NEAR(mats[1][1][1].real(), 1.0, 1e-15);
	EXPECT_NEAR(mats[1][2][0].real(), r, 1e-15);
	EXPECT_NEAR(mats[1][2][2].real(), -r, 1e-15);
	EXPECT_NEAR(std::abs(mats[1][0][1]), 0.0, 1e-15);
	std::unordered_map<std::string, int> lmax{ { "H", 2 }, { "C", 5 }, { "O", 4 } };
	EXPECT_EQ(SALTED_Utils::get_lmax_max(lmax), 5);
}

// every atom known to the model and inside the cutoff: nothing is removed and no fill is requested
TEST(SaltedFchkUtilTests, FilterInputKeepsCompleteKnownStructure)
{
	WFN w;
	w.push_back_atom("O", 0.0, 0.0, 0.0, 8);
	w.push_back_atom("H", 1.0, 0.0, 0.0, 1);
	w.push_back_atom("H", 0.0, 1.0, 0.0, 1);
	options opt;
	SALTEDConfig cfg{};
	cfg.species = { "H", "O" };
	cfg.rcut1 = cfg.rcut2 = 4.0;
	EXPECT_TRUE(SALTED_Utils::filter_input(w, opt, cfg).empty());
	EXPECT_EQ(w.get_ncen(), 3);
	EXPECT_FALSE(opt.needs_Thakkar_fill);
}

// a species the model never saw is erased and flagged for the spherical fill
TEST(SaltedFchkUtilTests, FilterInputRemovesUnknownSpecies)
{
	WFN w;
	w.push_back_atom("O", 0.0, 0.0, 0.0, 8);
	w.push_back_atom("C", 1.0, 0.0, 0.0, 6);
	w.push_back_atom("H", 0.0, 1.0, 0.0, 1);
	options opt;
	SALTEDConfig cfg{};
	cfg.species = { "H", "O" };
	cfg.rcut1 = 4.0;
	cfg.rcut2 = 3.0;
	const std::vector<char> removed = SALTED_Utils::filter_input(w, opt, cfg);
	ASSERT_EQ(removed.size(), 3u);
	EXPECT_EQ(removed[0], 0);
	EXPECT_EQ(removed[1], 1);
	EXPECT_EQ(removed[2], 0);
	EXPECT_EQ(w.get_ncen(), 2);
	EXPECT_EQ(w.get_atom_charge(1), 1);
	EXPECT_TRUE(opt.needs_Thakkar_fill);
}

// a known atom with nothing inside min(rcut1, rcut2) has no environment and goes to the fill too
TEST(SaltedFchkUtilTests, FilterInputRemovesIsolatedAtom)
{
	WFN w;
	w.push_back_atom("O", 0.0, 0.0, 0.0, 8);
	w.push_back_atom("H", 1.0, 0.0, 0.0, 1);
	w.push_back_atom("O", 60.0, 0.0, 0.0, 8);
	options opt;
	SALTEDConfig cfg{};
	cfg.species = { "H", "O" };
	cfg.rcut1 = cfg.rcut2 = 4.0;
	const std::vector<char> removed = SALTED_Utils::filter_input(w, opt, cfg);
	ASSERT_EQ(removed.size(), 3u);
	EXPECT_EQ(removed[2], 1);
	EXPECT_EQ(removed[0] + removed[1], 0);
	EXPECT_EQ(w.get_ncen(), 2);
	EXPECT_TRUE(opt.needs_Thakkar_fill);
}

// a 2 % surplus is scaled out of the s coefficient only, the p coefficients stay, and the electron count becomes exact
TEST(SaltedFchkUtilTests, ChargeConstraintScalesOnlyS)
{
	const atom A = aux_atom("H", 1, 0.0, 0.0, 0.0, { { 1.0, 0 }, { 0.8, 1 } });
	const double per_unit = electrons_per_unit(1.0);
	vec coefs{ 1.02 / per_unit, 0.3, 0.4, 0.5 };
	std::ostringstream log;
	const double f = apply_charge_constraint({ A }, coefs, 0, false, 0, 0.0, 0.0, log);
	EXPECT_NEAR(f, 1.0 / 1.02, 1e-12);
	EXPECT_NEAR(calc_atomic_density({ A }, coefs)[0], 1.0, 1e-10);
	EXPECT_NEAR(coefs[1], 0.3, 1e-15);
	EXPECT_NEAR(coefs[3], 0.5, 1e-15);
	EXPECT_NE(log.str().find("Charge constraint applied"), std::string::npos);
	EXPECT_NE(log.str().find("NOTE: correction"), std::string::npos);
	// 0.2 % is inside the training scatter: applied, but without the caution note
	vec tiny_surplus{ 1.002 / per_unit };
	std::ostringstream quiet;
	const atom S = aux_atom("H", 1, 0.0, 0.0, 0.0, { { 1.0, 0 } });
	EXPECT_NEAR(apply_charge_constraint({ S }, tiny_surplus, 0, false, 0, 0.0, 0.0, quiet), 1.0 / 1.002, 1e-12);
	EXPECT_EQ(quiet.str().find("NOTE: correction"), std::string::npos);
}

// more than 5 % off means something else is wrong: refused, coefficients untouched
TEST(SaltedFchkUtilTests, ChargeConstraintRefusesLargeFactor)
{
	const atom A = aux_atom("H", 1, 0.0, 0.0, 0.0, { { 1.0, 0 } });
	const double per_unit = electrons_per_unit(1.0);
	vec coefs{ 1.10 / per_unit };
	const vec before = coefs;
	std::ostringstream log;
	EXPECT_NEAR(apply_charge_constraint({ A }, coefs, 0, false, 0, 0.0, 0.0, log), 1.0, 1e-15);
	EXPECT_NEAR(coefs[0], before[0], 1e-15);
	EXPECT_NE(log.str().find("SKIPPED: factor"), std::string::npos);
}

// a net charge moves the target: an anion of two Z = 1 atoms holds 3 electrons
TEST(SaltedFchkUtilTests, ChargeConstraintHonoursNetCharge)
{
	const atom A = aux_atom("H", 1, 0.0, 0.0, 0.0, { { 1.0, 0 } });
	const atom B = aux_atom("H", 1, 0.0, 0.0, 1.5, { { 1.0, 0 } });
	const double per_unit = electrons_per_unit(1.0);
	vec coefs{ 1.48 / per_unit, 1.48 / per_unit };
	std::ostringstream log;
	const double f = apply_charge_constraint({ A, B }, coefs, -1, false, 0, 0.0, 0.0, log);
	EXPECT_NEAR(f, 3.0 / 2.96, 1e-12);
	const vec e = calc_atomic_density({ A, B }, coefs);
	EXPECT_NEAR(e[0] + e[1], 3.0, 1e-10);
	EXPECT_NE(log.str().find("net charge -1 taken into account"), std::string::npos);
	// with a spherical fill the split of that charge is undefined: refused
	vec again{ 1.48 / per_unit, 1.48 / per_unit };
	std::ostringstream log2;
	EXPECT_NEAR(apply_charge_constraint({ A, B }, again, -1, true, 1, 0.0, 0.0, log2), 1.0, 1e-15);
	EXPECT_NEAR(again[0], 1.48 / per_unit, 1e-15);
	EXPECT_NE(log2.str().find("SKIPPED: net charge"), std::string::npos);
}

// the fill notes: charge moved onto filled ions shifts the target, an unknown eeq charge is said so, an unappliable one too
TEST(SaltedFchkUtilTests, ChargeConstraintFillNotes)
{
	const atom A = aux_atom("H", 1, 0.0, 0.0, 0.0, { { 1.0, 0 } });
	const double per_unit = electrons_per_unit(1.0);
	{
		// +0.02 e went onto a filled cation, so the predicted region must hold 1.02
		vec coefs{ 1.0 / per_unit };
		std::ostringstream log;
		const double f = apply_charge_constraint({ A }, coefs, 0, true, 1, 0.02, 0.02, log);
		EXPECT_NEAR(f, 1.02, 1e-12);
		EXPECT_NEAR(calc_atomic_density({ A }, coefs)[0], 1.02, 1e-10);
		EXPECT_NE(log.str().find("moved to the spherically filled"), std::string::npos);
		EXPECT_NE(log.str().find("EEQ puts"), std::string::npos);
		EXPECT_EQ(log.str().find("could not be applied"), std::string::npos);
		EXPECT_NE(log.str().find("ML-predicted"), std::string::npos);
	}
	{
		vec coefs{ 1.0 / per_unit };
		std::ostringstream log;
		apply_charge_constraint({ A }, coefs, 0, true, 2, std::numeric_limits<double>::quiet_NaN(), 0.0, log);
		EXPECT_NE(log.str().find("could not be estimated"), std::string::npos);
	}
	{
		vec coefs{ 1.0 / per_unit };
		std::ostringstream log;
		apply_charge_constraint({ A }, coefs, 0, true, 1, 0.3, 0.0, log);
		EXPECT_NE(log.str().find("could not be applied"), std::string::npos);
	}
}

// a negative coefficient gives a non-positive electron count: skipped, factor one
TEST(SaltedFchkUtilTests, ChargeConstraintSkipsNonPositive)
{
	const atom A = aux_atom("H", 1, 0.0, 0.0, 0.0, { { 1.0, 0 } });
	vec coefs{ -0.5 };
	std::ostringstream log;
	EXPECT_NEAR(apply_charge_constraint({ A }, coefs, 0, false, 0, 0.0, 0.0, log), 1.0, 1e-15);
	EXPECT_NEAR(coefs[0], -0.5, 1e-15);
	EXPECT_NE(log.str().find("non-positive"), std::string::npos);
}

// the cube overload evaluates the table on the grid; atom_nr slices out that atom's coefficients
TEST(SaltedFchkUtilTests, CubeMLMatchesTable)
{
	WFN w;
	w.push_back_atom(aux_atom("H", 1, 0.0, 0.0, 0.0, { { 1.0, 0 }, { 0.6, 1 } }));
	w.push_back_atom(aux_atom("H", 1, 0.0, 0.0, 1.4, { { 0.8, 0 } }));
	const vec coefs{ 0.7, 0.1, -0.2, 0.3, 0.9 };
	const std::vector<atom> atoms = w.get_atoms();
	const aux_density_table full(atoms), second({ atoms[1] });
	ASSERT_EQ(full.n_coef, 5);
	cube c(std::array<int, 3>{ 3, 3, 3 }, 2, true);
	for (int i = 0; i < 3; i++)
	{
		c.set_origin(i, -0.5);
		for (int j = 0; j < 3; j++)
			c.set_vector(i, j, i == j ? 0.5 : 0.0);
	}
	calc_cube_ML(coefs, w, c);
	double biggest = 0.0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
			{
				const double x = -0.5 + 0.5 * i, y = -0.5 + 0.5 * j, z = -0.5 + 0.5 * k;
				EXPECT_NEAR(c.get_value(i, j, k), full(x, y, z, coefs.data()), 1e-14) << i << j << k;
				biggest = std::max(biggest, std::abs(c.get_value(i, j, k)));
			}
	EXPECT_GT(biggest, 1e-3);
	EXPECT_NEAR(c.get_dv(), 0.125, 1e-15);
	// hand-derived value at grid point (2, 0, 1) = (0.5, -0.5, 0): atom 0 sees r^2 = 0.5 and its p shell, in the
	// m = -1, 0, +1 order of the coefficients, is sqrt(3/4pi) (c0 y + c1 z + c2 x); atom 1 sees dz = -1.4, r^2 = 2.46.
	// dx != dy != dz here, so a swapped m order or a lost normalisation changes the number
	const double sq14 = 1.0 / std::sqrt(4.0 * constants::PI), sq34 = std::sqrt(3.0 / (4.0 * constants::PI));
	const double want = 0.7 * radial_norm(1.0, 0) * std::exp(-0.5) * sq14
		+ radial_norm(0.6, 1) * std::exp(-0.6 * 0.5) * sq34 * (0.1 * (-0.5) + (-0.2) * 0.0 + 0.3 * 0.5)
		+ 0.9 * radial_norm(0.8, 0) * std::exp(-0.8 * 2.46) * sq14;
	EXPECT_NEAR(c.get_value(2, 0, 1), want, 1e-13);
	cube one(std::array<int, 3>{ 3, 3, 3 }, 2, true);
	for (int i = 0; i < 3; i++)
	{
		one.set_origin(i, -0.5);
		for (int j = 0; j < 3; j++)
			one.set_vector(i, j, i == j ? 0.5 : 0.0);
	}
	calc_cube_ML(coefs, w, one, 1);
	const double slice[1]{ 0.9 };
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				EXPECT_NEAR(one.get_value(i, j, k), second(-0.5 + 0.5 * i, -0.5 + 0.5 * j, -0.5 + 0.5 * k, slice), 1e-14) << i << j << k;
}

// the WFN overload sizes the box itself; the integrated s density is the analytic shell population, and both
// analytic routes (table and calc_atomic_density share normalize_gto) agree with the hand-derived (2 pi / a)^(3/4)
TEST(SaltedFchkUtilTests, CubeMLIntegratesToPopulation)
{
	WFN w;
	w.push_back_atom(aux_atom("H", 1, 0.0, 0.0, 0.0, { { 1.0, 0 } }));
	const aux_density_table t(w.get_atoms());
	const vec coefs{ 1.3 };
	const double want = 1.3 * electrons_per_unit(1.0);
	EXPECT_NEAR(t.shell_population_integral(0), electrons_per_unit(1.0), 1e-12);
	EXPECT_NEAR(calc_atomic_density(w.get_atoms(), coefs)[0], want, 1e-12);
	const cube c = calc_cube_ML(coefs, w);
	EXPECT_GE(c.get_size(0), 30);
	EXPECT_NEAR(c.sum(), want, 1e-3 * want);
}

// ------------------------------------------------------------- equicomb

// the gpu switch is a plain flag; the sparse path checks dimensions before touching any buffer
TEST(SaltedFchkEquicombTests, GpuSwitchAndArgumentChecks)
{
	const bool before = equicomb_gpu_enabled();
	equicomb_set_gpu(true);
	EXPECT_TRUE(equicomb_gpu_enabled());
	equicomb_set_gpu(false);
	EXPECT_FALSE(equicomb_gpu_enabled());
	equicomb_set_gpu(before);

	SALTEDDescriptors v(1, 1, 0);
	const vec w3j{ 1.0 };
	const ivec2 llvec{ { 0 }, { 0 } };
	const cvec2 c2r{ { 1.0 } };
	const std::vector<int64_t> vfps{ 0 };
	vec p(1, 42.0);
	EXPECT_THROW(equicomb(-1, 1, 1, v, v, w3j, llvec, 0, c2r, 1, 1, vfps, p), std::invalid_argument);
	EXPECT_THROW(equicomb(1, 1, 1, v, v, w3j, llvec, -2, c2r, 1, 1, vfps, p), std::invalid_argument);
	vec tiny;
	EXPECT_THROW(equicomb(1, 1, 1, v, v, w3j, llvec, 0, c2r, 1, 1, vfps, tiny), std::out_of_range);
	// zero-sized problems return before writing anything
	equicomb(0, 1, 1, v, v, w3j, llvec, 0, c2r, 1, 1, vfps, p);
	equicomb(1, 0, 1, v, v, w3j, llvec, 0, c2r, 1, 1, vfps, p);
	EXPECT_NEAR(p[0], 42.0, 1e-15);

	vec w3j_d{ 1.0 };
	ivec2 llvec_d{ { 0 }, { 0 } };
	cvec2 c2r_d{ { 1.0 } };
	EXPECT_THROW(equicomb(1, -1, 1, v, v, w3j_d, 1, llvec_d, 0, c2r_d, 1, p), std::invalid_argument);
	EXPECT_THROW(equicomb(1, 1, 1, v, v, w3j_d, 1, llvec_d, -1, c2r_d, 1, p), std::invalid_argument);
	EXPECT_THROW(equicomb(1, 1, 1, v, v, w3j_d, 1, llvec_d, 0, c2r_d, 1, tiny), std::out_of_range);
	equicomb(1, 1, 0, v, v, w3j_d, 1, llvec_d, 0, c2r_d, 1, p);
	EXPECT_NEAR(p[0], 42.0, 1e-15);
}

// lam = 0 with two radial channels: the feature vector is (a, b) c w normalised, so the answer is (a, b)/|(a, b)|;
// the sparse path reorders through vfps and the conj flag conjugates v1 in place of v2
TEST(SaltedFchkEquicombTests, Lam0AnalyticNormalisation)
{
	equicomb_set_gpu(false);
	SALTEDDescriptors v1(1, 2, 0), v2(1, 1, 0);
	v1.block(0, 0, 0)[0] = 3.0;
	v1.block(0, 1, 0)[0] = 4.0;
	v2.block(0, 0, 0)[0] = 2.0;
	vec w3j{ 0.5 };
	ivec2 llvec{ { 0 }, { 0 } };
	cvec2 c2r{ { 1.0 } };
	vec dense(2, 0.0);
	equicomb(1, 2, 1, v1, v2, w3j, 1, llvec, 0, c2r, 2, dense);
	EXPECT_NEAR(dense[0], 0.6, 1e-14);
	EXPECT_NEAR(dense[1], 0.8, 1e-14);
	vec sparse(2, 0.0);
	equicomb(1, 2, 1, v1, v2, w3j, llvec, 0, c2r, 2, 2, std::vector<int64_t>{ 1, 0 }, sparse);
	EXPECT_NEAR(sparse[0], 0.8, 1e-14);
	EXPECT_NEAR(sparse[1], 0.6, 1e-14);
	// complex a: with conjugation a conj(a) is real and b conj(a) imaginary, without it a a is real and negative
	v1.block(0, 0, 0)[0] = cdouble(0.0, 3.0);
	vec conj_p(2, 0.0), plain(2, 0.0);
	equicomb(1, 2, 1, v1, v1, w3j, llvec, 0, c2r, 2, 2, std::vector<int64_t>{ 0, 1 }, conj_p, true);
	EXPECT_NEAR(conj_p[0], 1.0, 1e-14);
	EXPECT_NEAR(conj_p[1], 0.0, 1e-14);
	equicomb(1, 2, 1, v1, v1, w3j, 1, llvec, 0, c2r, 2, plain, false);
	EXPECT_NEAR(plain[0], -1.0, 1e-14);
	EXPECT_NEAR(plain[1], 0.0, 1e-14);
	vec plain_sparse(2, 0.0);
	equicomb(1, 2, 1, v1, v1, w3j, llvec, 0, c2r, 2, 2, std::vector<int64_t>{ 0, 1 }, plain_sparse, false);
	EXPECT_NEAR(plain_sparse[0], -1.0, 1e-14);
}

// lam = 1 over (l1, l2) = (1, 0) and (0, 1) with complex descriptors: the run-based sparse loop must reproduce
// the m-by-m dense loop, in both feature orders and with the conj shortcut
TEST(SaltedFchkEquicombTests, SparseMatchesDenseLam1)
{
	equicomb_set_gpu(false);
	const int natoms = 2, nrad1 = 2, nrad2 = 2;
	SALTEDDescriptors v1(natoms, nrad1, 1), v2(natoms, nrad2, 1);
	double seed = 0.37;
	auto next = [&seed]()
	{
		seed = std::fmod(seed * 7.31 + 0.113, 1.0);
		return seed - 0.5;
	};
	for (auto& x : v1.values())
		x = cdouble(next(), next());
	for (auto& x : v2.values())
		x = cdouble(next(), next());
	vec w3j{ 0.3, -0.5, 0.7, 0.2, 0.9, -0.4 };
	ivec2 llvec{ { 1, 0 }, { 0, 1 } };
	cvec2 c2r = SALTED_Utils::complex_to_real_transformation({ 3 })[0];
	const int llmax = 2, featsize = nrad1 * nrad2 * llmax, l21 = 3;
	vec dense(static_cast<size_t>(natoms) * l21 * featsize, 0.0);
	equicomb(natoms, nrad1, nrad2, v1, v2, w3j, llmax, llvec, 1, c2r, featsize, dense);
	double norm = 0.0;
	for (int f = 0; f < l21 * featsize; f++)
		norm += dense[f] * dense[f];
	EXPECT_NEAR(norm, 1.0, 1e-12);
	std::vector<int64_t> identity(featsize);
	for (int f = 0; f < featsize; f++)
		identity[f] = f;
	vec sparse(dense.size(), 0.0);
	equicomb(natoms, nrad1, nrad2, v1, v2, w3j, llvec, 1, c2r, featsize, featsize, identity, sparse);
	for (size_t i = 0; i < dense.size(); i++)
		EXPECT_NEAR(sparse[i], dense[i], 1e-13) << i;
	// a subset in another order: p[iat*l21*nfps + imu*nfps + i] = dense[iat*l21*featsize + imu*featsize + vfps[i]]
	const std::vector<int64_t> vfps{ 5, 0, 7 };
	vec sub(static_cast<size_t>(natoms) * l21 * 3, 0.0);
	equicomb(natoms, nrad1, nrad2, v1, v2, w3j, llvec, 1, c2r, featsize, 3, vfps, sub);
	for (int iat = 0; iat < natoms; iat++)
		for (int imu = 0; imu < l21; imu++)
			for (int i = 0; i < 3; i++)
				EXPECT_NEAR(sub[iat * l21 * 3 + imu * 3 + i], dense[iat * l21 * featsize + imu * featsize + vfps[i]], 1e-13) << iat << imu << i;
	// v2 = conj(v1) handed explicitly must equal the flag
	SALTEDDescriptors v1c = v1;
	for (auto& x : v1c.values())
		x = std::conj(x);
	vec explicit_conj(dense.size(), 0.0), flagged(dense.size(), 0.0), flagged_dense(dense.size(), 0.0);
	equicomb(natoms, nrad1, nrad1, v1, v1c, w3j, llvec, 1, c2r, featsize, featsize, identity, explicit_conj);
	equicomb(natoms, nrad1, nrad1, v1, v1, w3j, llvec, 1, c2r, featsize, featsize, identity, flagged, true);
	equicomb(natoms, nrad1, nrad1, v1, v1, w3j, llmax, llvec, 1, c2r, featsize, flagged_dense, true);
	for (size_t i = 0; i < dense.size(); i++)
	{
		EXPECT_NEAR(flagged[i], explicit_conj[i], 1e-13) << i;
		EXPECT_NEAR(flagged_dense[i], explicit_conj[i], 1e-13) << i;
	}
}

// lam = 1 from (l1, l2) = (1, 0) with one radial channel, hand-derived: m2 = m1 - mu has to be 0, so only
// m1 = mu survives, p_mu = w3j[mu + 1] v1[mu] conj(v0) in (il, imu, im1) order, and the real form applies the
// l = 1 rows (i, 0, i)/sqrt2, (0, 1, 0), (1, 0, -1)/sqrt2 before the unit normalisation. Sparse and dense
// share none of this arithmetic with the test, so a wrong m selection, w3j order, conj sign or c2r row fails here
TEST(SaltedFchkEquicombTests, Lam1HandDerived)
{
	equicomb_set_gpu(false);
	SALTEDDescriptors v(1, 1, 1);
	const cdouble a(0.3, -0.7), b(-0.4, 0.2), c(0.9, 0.5), d(0.6, -0.8);
	v.block(0, 0, 1)[0] = a;
	v.block(0, 0, 1)[1] = b;
	v.block(0, 0, 1)[2] = c;
	v.block(0, 0, 0)[0] = d;
	vec w3j{ 0.5, -0.25, 0.75 };
	ivec2 llvec{ { 1 }, { 0 } };
	cvec2 c2r = SALTED_Utils::complex_to_real_transformation({ 3 })[0];
	const cdouble p0 = w3j[0] * a * std::conj(d), p1 = w3j[1] * b * std::conj(d), p2 = w3j[2] * c * std::conj(d);
	const double r = 1.0 / std::sqrt(2.0);
	vec want{ -r * (p0 + p2).imag(), p1.real(), r * (p0 - p2).real() };
	const double norm = std::sqrt(want[0] * want[0] + want[1] * want[1] + want[2] * want[2]);
	for (double& x : want)
		x /= norm;
	vec sparse(3, 0.0), dense(3, 0.0);
	equicomb(1, 1, 1, v, v, w3j, llvec, 1, c2r, 1, 1, std::vector<int64_t>{ 0 }, sparse, true);
	equicomb(1, 1, 1, v, v, w3j, 1, llvec, 1, c2r, 1, dense, true);
	for (int imu = 0; imu < 3; imu++)
	{
		EXPECT_NEAR(sparse[imu], want[imu], 1e-14) << imu;
		EXPECT_NEAR(dense[imu], want[imu], 1e-14) << imu;
	}
}

// an all-zero descriptor is an empty environment: zeros, never NaN, in both paths
TEST(SaltedFchkEquicombTests, EmptyEnvironmentGivesZeros)
{
	equicomb_set_gpu(false);
	SALTEDDescriptors v(1, 1, 0);
	vec w3j{ 1.0 };
	ivec2 llvec{ { 0 }, { 0 } };
	cvec2 c2r{ { 1.0 } };
	vec sparse(1, 5.0), dense(1, 5.0);
	equicomb(1, 1, 1, v, v, w3j, llvec, 0, c2r, 1, 1, std::vector<int64_t>{ 0 }, sparse, true);
	equicomb(1, 1, 1, v, v, w3j, 1, llvec, 0, c2r, 1, dense, true);
	EXPECT_NEAR(sparse[0], 0.0, 1e-300);
	EXPECT_NEAR(dense[0], 0.0, 1e-300);
	EXPECT_FALSE(std::isnan(sparse[0]));
}

// --------------------------------------------------------- SALTED_predictor

// a coefficient file with an auxiliary basis needs no model: the aux wavefunction is built and the npy returned as is
TEST(SaltedFchkPredictorTests, CoefficientFilePath)
{
	const auto root = nos_test_repo_root() / "tests" / "SALTED";
	WFN w(root / "test_cysteine.xyz", false);
	ASSERT_EQ(w.get_ncen(), 14);
	options opt;
	opt.coef_file = root / "test_cysteine.npy";
	opt.aux_basis = { BasisSetLibrary::get_basis_set("cc-pvqz-jkfit") };
	SALTEDPredictor SP(w, opt);
	EXPECT_TRUE(SP.basis_set_loaded());
	EXPECT_EQ(SP.get_dfbasis_name(), opt.aux_basis[0]->get_name());
	EXPECT_EQ(SP.get_salted_filename().string(), "coefficient file");
	EXPECT_EQ(SP.wavy.get_ncen(), 14);
	const vec coefs = SP.gen_SALTED_densities();
	EXPECT_EQ(coefs.size(), 1141u);
	EXPECT_EQ(coefs.size(), static_cast<size_t>(aux_density_table(SP.wavy.get_atoms()).n_coef));
	vec direct;
	std::filesystem::path np = opt.coef_file;
	read_npy<double>(np, direct);
	ASSERT_EQ(direct.size(), coefs.size());
	EXPECT_NEAR(coefs[0], direct[0], 1e-15);
	EXPECT_NEAR(coefs.back(), direct.back(), 1e-15);
}

// a species the shipped model never saw is handed to the spherical fill with an eeq charge estimate
TEST(SaltedFchkPredictorTests, UnknownSpeciesGoesToThakkarFill)
{
	const auto root = nos_test_repo_root() / "tests";
	WFN w(root / "reading_SALTED" / "water_monomer.xyz", false);
	ASSERT_EQ(w.get_ncen(), 3);
	w.push_back_atom("Cl", 5.0, 0.0, 0.0, 17);
	options opt;
	opt.salted_model_dir = root / "SALTED" / "Model";
	SALTEDPredictor SP(w, opt);
	EXPECT_EQ(SP.get_salted_filename().string(), "model.salted");
	EXPECT_EQ(SP.get_dfbasis_name(), "cc-pvqz-jkfit");
	EXPECT_FALSE(SP.basis_set_loaded());
	EXPECT_TRUE(opt.needs_Thakkar_fill);
	EXPECT_EQ(SP.wavy.get_ncen(), 3);
	EXPECT_EQ(SP.wavy.get_nmo(), 0);
	ASSERT_EQ(opt.spherical_fill_charges.size(), 1u);
	for (int ax = 0; ax < 3; ax++)
		EXPECT_NEAR(opt.spherical_fill_charges[0][ax], w.get_atom_coordinate(3, ax), 1e-12);
	EXPECT_TRUE(std::isfinite(opt.spherical_fill_charges[0][3]));
	// the value itself is occ's EEQ, not ours, and depends on its parameter table: only its plausibility is checked
	EXPECT_LT(std::abs(opt.spherical_fill_charges[0][3]), 2.0);
}

// a model that ships its own BASIS block loads it into the predictor's wavefunction
TEST(SaltedFchkPredictorTests, SyntheticModelLoadsBasis)
{
	const auto dir = tmp_path("basismodel");
	std::filesystem::create_directories(dir);
	write_synthetic_model(dir / "tiny.salted", 3, true, true);
	WFN w;
	w.push_back_atom("O", 0.0, 0.0, 0.0, 8);
	w.push_back_atom("H", 1.8, 0.0, 0.0, 1);
	w.push_back_atom("H", 0.0, 1.8, 0.0, 1);
	w.push_back_MO(1, 2.0, -1.0);
	options opt;
	opt.salted_model_dir = dir;
	{
		SALTEDPredictor SP(w, opt);
		EXPECT_EQ(SP.get_salted_filename().string(), "tiny.salted");
		EXPECT_EQ(SP.get_dfbasis_name(), "tiny-jkfit");
		EXPECT_TRUE(SP.basis_set_loaded());
		EXPECT_FALSE(opt.needs_Thakkar_fill);
		EXPECT_EQ(SP.wavy.get_nmo(), 0);
		ASSERT_EQ(SP.wavy.get_ncen(), 3);
		const aux_density_table t(SP.wavy.get_atoms());
		EXPECT_EQ(t.n_pr, 4);
		EXPECT_EQ(t.n_coef, 6);
		EXPECT_NEAR(SP.wavy.get_atoms()[0].get_basis_set_entry(1).get_exponent(), 0.7, 1e-15);
		EXPECT_NEAR(SP.wavy.get_atoms()[1].get_basis_set_entry(0).get_exponent(), 1.5, 1e-15);
	}
	std::filesystem::remove_all(dir);
}

// a model directory without a .salted file ends the run with exit code 1
TEST(SaltedFchkPredictorTests, EmptyModelDirExits)
{
	const auto dir = tmp_path("emptymodel");
	std::filesystem::create_directories(dir);
	WFN w;
	w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
	options opt;
	opt.salted_model_dir = dir;
	EXPECT_EXIT(SALTEDPredictor SP(w, opt), ::testing::ExitedWithCode(1), ".*");
	std::filesystem::remove_all(dir);
}

// the full prediction on the water monomer with the shipped model: the fitted density holds ten electrons
TEST(SaltedFchkPredictorTests, PredictWaterMonomer_full)
{
	if (!full_tests_enabled())
		GTEST_SKIP() << "Set RUN_FULL_TEST=1 to run the SALTED prediction";
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	WFN w(root / "reading_SALTED" / "water_monomer.xyz", false);
	options opt;
	opt.salted_model_dir = root / "SALTED" / "Model";
	SALTEDPredictor SP(w, opt);
	load_basis_into_WFN(SP.wavy, BasisSetLibrary::get_basis_set(SP.get_dfbasis_name()));
	const vec coefs = SP.gen_SALTED_densities();
	const aux_density_table t(SP.wavy.get_atoms());
	ASSERT_EQ(coefs.size(), static_cast<size_t>(t.n_coef));
	const vec e = calc_atomic_density(SP.wavy.get_atoms(), coefs);
	// the fit is ~0.24 % short on training-like systems and 0.016 % on the cysteine golden; 1 % is four times
	// the worst of those and still fails on any lost term, since the species average alone is not 10
	EXPECT_NEAR(e[0] + e[1] + e[2], 10.0, 0.1);
}

// ------------------------------------------------------------------ fchk

// an unrestricted doublet written through the beta branch and read back keeps its header and coefficient blocks;
// occ puts the unpaired electron of a 1-electron doublet in the beta channel, so the file carries the UHF header,
// zero alpha and one beta electron (the counts follow the occupations per MO operator) and a full beta block
TEST(SaltedFchkTests, UnrestrictedDoubletRoundTrip)
{
	const auto tmp = tmp_path("h2p.fchk");
	std::filesystem::remove(tmp);
	WFN w = occ_h2_plus();
	ASSERT_TRUE(w.get_is_unrestricted());
	ASSERT_EQ(w.get_nmo(), 16);
	w.set_origin(e_origin::wfn);
	w.set_method("rhf");
	std::ostringstream log;
	ASSERT_TRUE(free_fchk(log, tmp, "", w, false, true));
	EXPECT_EQ(second_line(tmp).substr(0, 13), "SP        UHF");
	{
		std::ifstream in(tmp);
		EXPECT_EQ(read_fchk_integer(in, "Number of alpha electrons"), 0);
		EXPECT_EQ(read_fchk_integer(in, "Number of beta electrons"), 1);
		EXPECT_EQ(read_fchk_integer(in, "Multiplicity"), 2);
		EXPECT_EQ(read_fchk_integer(in, "Number of basis functions"), 8);
		vec beta, alpha, be;
		ASSERT_TRUE(read_fchk_double_block(in, "Beta MO coefficients", beta));
		ASSERT_TRUE(read_fchk_double_block(in, "Alpha MO coefficients", alpha));
		ASSERT_TRUE(read_fchk_double_block(in, "Beta Orbital Energies", be));
		EXPECT_EQ(beta.size(), 64u);
		EXPECT_EQ(alpha.size(), 64u);
		ASSERT_EQ(be.size(), 8u);
		EXPECT_NEAR(be[0], w.get_MO_energy(8), 1e-7);
		double diff = 0.0;
		for (size_t i = 0; i < 64; i++)
			diff += std::abs(beta[i] - alpha[i]);
		EXPECT_GT(diff, 1e-3);
	}
	WFN back(tmp, false);
	EXPECT_TRUE(back.get_is_unrestricted());
	EXPECT_EQ(back.get_multi(), 2u);
	EXPECT_EQ(back.get_charge(), 1);
	EXPECT_EQ(back.get_nmo(), 16);
	std::filesystem::remove(tmp);
}

// the occ-built H2+ doublet written to fchk and read back keeps its density (the occupied orbital sits in the beta block)
TEST(SaltedFchkTests, UnrestrictedDoubletRoundTripKeepsDensity)
{
	const auto tmp = tmp_path("h2p_density.fchk");
	std::filesystem::remove(tmp);
	WFN w = occ_h2_plus();
	w.set_origin(e_origin::wfn);
	w.set_method("rhf");
	std::ostringstream log;
	ASSERT_TRUE(free_fchk(log, tmp, "", w, false, true));
	WFN back(tmp, false);
	expect_same_density(w, back, 1e-6, "H2+ via fchk");
	std::filesystem::remove(tmp);
}

// the method line: rks/rhf with a singlet or higher multiplicity pick the four Gaussian labels, anything else B3LYP
TEST(SaltedFchkTests, MethodHeaderStrings)
{
	const auto tmp = tmp_path("method.fchk");
	std::filesystem::remove(tmp);
	WFN w = occ_h2_plus();
	w.set_origin(e_origin::wfn);
	std::ostringstream log;
	w.set_method("rks");
	ASSERT_TRUE(free_fchk(log, tmp, "", w, false, true));
	EXPECT_EQ(second_line(tmp).substr(0, 16), "SP        UB3LYP");
	w.set_method("");
	ASSERT_TRUE(free_fchk(log, tmp, "", w, false, true));
	EXPECT_EQ(second_line(tmp).substr(0, 15), "SP        B3LYP");
	std::filesystem::remove(tmp);
}

// without force_overwrite an existing fchk is left alone and the call reports it
TEST(SaltedFchkTests, ExistingFileNotOverwritten)
{
	const auto tmp = tmp_path("exists.fchk");
	{
		std::ofstream(tmp) << "placeholder\n";
	}
	WFN w = occ_h2_plus();
	w.set_origin(e_origin::wfn);
	w.set_method("rhf");
	std::ostringstream log;
	EXPECT_FALSE(free_fchk(log, tmp, "", w, false, false));
	EXPECT_NE(log.str().find("already exists"), std::string::npos);
	EXPECT_EQ(second_line(tmp), "");
	std::filesystem::remove(tmp);
}

// a wavefunction whose origin free_fchk does not know is refused with an err_checkf exit
TEST(SaltedFchkTests, UnsupportedOriginExits)
{
	const auto tmp = tmp_path("occ_origin.fchk");
	WFN w = occ_h2_plus();
	ASSERT_EQ(w.get_origin(), e_origin::OCC);
	std::ostringstream log;
	EXPECT_EXIT(free_fchk(log, tmp, "", w, false, true), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	std::filesystem::remove(tmp);
}
