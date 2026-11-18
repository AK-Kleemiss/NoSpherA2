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
#include "core/nos_math.h"
#include "core/NoSpherA2.h"
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#include <random>
#include <regex>
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
			// carbon carries one contracted shell: two primitives sharing one angular momentum
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

// a model folder: the lead is the first non-element name, element models join only for elements present
TEST(SaltedFchkIoTests, SaltedFolderModels)
{
	const auto dir = tmp_path("modelfolder");
	std::filesystem::create_directories(dir);
	for (const char* f : {"Co.salted", "Ag.salted", "v8.salted", "z_old.salted", "Cl.salted", "xyz.salted", "Fe.txt"})
		std::ofstream(dir / f) << "x";
	EXPECT_EQ(salted_element_file_Z("Co.salted"), 27);
	EXPECT_EQ(salted_element_file_Z("CO.salted"), 0);
	EXPECT_EQ(salted_element_file_Z("D.salted"), 0);
	EXPECT_EQ(find_first_salted_file(dir).string(), "v8.salted");
	const pathvec m = salted_folder_models(dir, {1, 6, 17, 27});
	ASSERT_EQ(m.size(), 3u);
	EXPECT_EQ(m[0].filename().string(), "v8.salted");
	EXPECT_EQ(m[1].filename().string(), "Cl.salted");
	EXPECT_EQ(m[2].filename().string(), "Co.salted");
	EXPECT_EQ(salted_folder_models(dir, {1, 6}).size(), 1u);
	for (const char* f : {"v8.salted", "z_old.salted", "xyz.salted"})
		std::filesystem::remove(dir / f);
	EXPECT_EQ(find_first_salted_file(dir).string(), "Ag.salted"); // element models only: the first file
	EXPECT_EQ(salted_folder_models(dir, {27}).size(), 1u);
	std::filesystem::remove_all(dir);
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

// load_blocks reads every block of the shipped model as load_block does, byte for byte, with one reader and with several
// whose pieces (an odd size, so they cut doubles) cross the blocks, more of them than blocks; a block past the end of the file still stops the run
TEST(SaltedFchkIoTests, ParallelReadMatchesLoadBlock)
{
	const auto p = nos_test_repo_root() / "tests" / "SALTED" / "Model" / "model.salted";
	SALTED_BINARY_FILE f(p);
	std::vector<SALTED_BINARY_FILE::block_ref> refs;
	for (const std::string key : { "FEATS", "PROJ" })
		for (const auto& [k, r] : f.index_lambda_based_data(key)) refs.push_back(r);
	refs.insert(refs.begin() + 3, SALTED_BINARY_FILE::block_ref{});
	for (const std::size_t readers : { 1, 2, 7, 64 })
		for (const std::size_t piece : { std::size_t(4u << 20), std::size_t(1001) })
		{
			f.readers = readers;
			f.piece = piece;
			const std::vector<dMatrix2> got = f.load_blocks(refs);
			ASSERT_EQ(got.size(), refs.size());
			for (size_t i = 0; i < refs.size(); i++)
			{
				const dMatrix2 want = f.load_block(refs[i]);
				ASSERT_EQ(got[i].size(), want.size()) << i;
				EXPECT_EQ(std::memcmp(got[i].data(), want.data(), want.size() * sizeof(double)), 0) << i << " readers " << readers << " piece " << piece;
			}
		}
	const SALTED_BINARY_FILE::block_ref past{ static_cast<std::streamoff>(std::filesystem::file_size(p)) - 800, 10, 20 };
	EXPECT_EXIT({ f.readers = 4; f.piece = 1001; f.load_blocks({ refs[0], past }); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

// -salted_fold's writer: PROJ becomes PROJW, the blocks behind it move and still read back; a block
// dropped without a replacement is gone
TEST(SaltedFchkIoTests, BlockReplacedCopy)
{
	const auto p = tmp_path("fold_in.salted"), q = tmp_path("fold_out.salted"), r = tmp_path("fold_drop.salted");
	write_synthetic_model(p, 3, true, true);
	{
		salted_writer w;
		w.block_head(1);
		w.tag("O");
		w.raw(static_cast<int32_t>(1));
		w.dataset(vec{ 1.5, 2.5, 3.5 }, { 3, 1 });
		SALTED_BINARY_FILE f(p);
		f.write_with_blocks(q, 4, { "PROJ" }, { { "PROJW", w.buf } });
		f.write_with_blocks(r, 4, { "WEIGH", "FEATS" }, {});
	}
	{
		SALTED_BINARY_FILE g(q);
		EXPECT_FALSE(g.has_block("PROJ"));
		const auto idx = g.index_lambda_based_data("PROJW");
		ASSERT_EQ(idx.size(), 1u);
		const dMatrix2 o = g.load_block(idx.at("O0"));
		ASSERT_EQ(o.extent(0), 3u);
		EXPECT_EQ(o(2, 0), 3.5);
		EXPECT_EQ(g.read_weights(), weights);
		EXPECT_EQ(g.read_features().at("O0")(2, 1), 14.0);
		ASSERT_TRUE(g.basis_set_defined());
		EXPECT_EQ(g.read_basis_set()->get_owned_primitive_count(), 5u);
	}
	{
		SALTED_BINARY_FILE g(r);
		EXPECT_FALSE(g.has_block("WEIGH"));
		EXPECT_FALSE(g.has_block("FEATS"));
		EXPECT_NEAR(g.read_projectors().at("H1")(0, 1), 80.0, 1e-15);
		EXPECT_EQ(g.read_basis_set()->get_owned_primitive_count(), 5u);
	}
	std::filesystem::remove(p);
	std::filesystem::remove(q);
	std::filesystem::remove(r);
}

// -salted_ecp: the ECP block reads back as the mode, a flagged copy is VERSION 6 and still reads, an unflagged
// model is all-electron (0), and re-flagging replaces the block instead of adding a second one
TEST(SaltedFchkIoTests, EcpBlockRoundTrip)
{
	const auto p = tmp_path("ecp_in.salted"), q = tmp_path("ecp_out.salted"), r = tmp_path("ecp_again.salted");
	write_synthetic_model(p, 3, true, true);
	write_salted_ecp(p, q, 1);
	write_salted_ecp(q, r, 3);
	SALTEDConfig c{}, cq{}, cr{};
	{
		SALTED_BINARY_FILE f(p), g(q), h(r);
		f.populate_config(c);
		g.populate_config(cq);
		h.populate_config(cr);
		EXPECT_EQ(g.read_weights(), weights);
		EXPECT_EQ(g.read_basis_set()->get_owned_primitive_count(), 5u);
	}
	EXPECT_EQ(c.ecp_mode, 0);
	EXPECT_EQ(cq.ecp_mode, 1);
	EXPECT_EQ(cr.ecp_mode, 3);
	EXPECT_EQ(std::filesystem::file_size(r), std::filesystem::file_size(q));
	EXPECT_EXIT(write_salted_ecp(p, q, 4), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	std::filesystem::remove(p);
	std::filesystem::remove(q);
	std::filesystem::remove(r);
}

// a file from the future warns but still reads through its table of contents
TEST(SaltedFchkIoTests, NewerVersionStillReads)
{
	const auto p = tmp_path("future.salted");
	write_synthetic_model(p, 7, false, true);
	{
		SALTED_BINARY_FILE f(p);
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

// the cell list finds exactly the atoms an all-pairs scan finds: a sparse cloud across many
// cells and negative coordinates, so pairs straddle every kind of cell boundary
TEST(SaltedFchkUtilTests, FilterInputCellListMatchesAllPairs)
{
	WFN w;
	std::mt19937 gen(7);
	std::uniform_real_distribution<double> u(-25.0, 25.0);
	for (int a = 0; a < 300; a++)
		w.push_back_atom(a % 2 ? "H" : "O", u(gen), u(gen), u(gen), a % 2 ? 1 : 8);
	const double rcut = 3.0;
	const double r = w.get_isBohr() ? constants::ang2bohr(rcut) : rcut;
	std::vector<char> expect(w.get_ncen(), 1);
	for (int a = 0; a < w.get_ncen(); a++)
		for (int b = 0; b < w.get_ncen(); b++)
		{
			double d_sq = 0.0;
			for (int ax = 0; ax < 3; ax++)
				d_sq += std::pow(w.get_atom_coordinate(a, ax) - w.get_atom_coordinate(b, ax), 2);
			if (b != a && d_sq < r * r) expect[a] = 0;
		}
	const int n_lonely = (int)std::count(expect.begin(), expect.end(), (char)1);
	ASSERT_GT(n_lonely, 10);
	ASSERT_LT(n_lonely, 290);
	options opt;
	SALTEDConfig cfg{};
	cfg.species = { "H", "O" };
	cfg.rcut1 = cfg.rcut2 = rcut;
	EXPECT_EQ(SALTED_Utils::filter_input(w, opt, cfg), expect);
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

namespace
{
	// the expansion of a real density, as featomic's are: v[n,l,-m] = (-1)^m conj(v[n,l,m]), so v[n,l,0] is real
	void make_real_density(SALTEDDescriptors& v, const int natoms)
	{
		for (int l = 0; l < static_cast<int>(v.offsets().size()); l++)
			for (int a = 0; a < natoms; a++)
				for (int n = 0; n < v.nchannels(); n++)
				{
					cdouble* x = v.block(a, n, l);
					x[l] = x[l].real();
					for (int m = 1; m <= l; m++)
						x[l + m] = (m % 2 ? -1.0 : 1.0) * std::conj(x[l - m]);
				}
	}

	// Racah's formula for (j1 j2 j3; m1 m2 m3) at integer arguments
	double wigner_3j(const int j1, const int j2, const int j3, const int m1, const int m2, const int m3)
	{
		if (m1 + m2 + m3 != 0 || j3 < std::abs(j1 - j2) || j3 > j1 + j2 || std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(m3) > j3)
			return 0.0;
		auto f = [](const int n) { return std::tgamma(n + 1.0); };
		double s = 0.0;
		for (int k = std::max({ 0, j2 - j3 - m1, j1 - j3 + m2 }); k <= std::min({ j1 + j2 - j3, j1 - m1, j2 + m2 }); k++)
			s += (k % 2 ? -1.0 : 1.0) / (f(k) * f(j3 - j2 + k + m1) * f(j3 - j1 + k - m2) * f(j1 + j2 - j3 - k) * f(j1 - k - m1) * f(j2 - k + m2));
		return ((j1 - j2 - m3) % 2 ? -1.0 : 1.0) * s * std::sqrt(f(j1 + j2 - j3) * f(j1 - j2 + j3) * f(j2 + j3 - j1) / f(j1 + j2 + j3 + 1)
			* f(j1 + m1) * f(j1 - m1) * f(j2 + m2) * f(j2 - m2) * f(j3 + m3) * f(j3 - m3));
	}

	// SALTED's coupling list: per (l1, l2), mu, m1 with m2 = m1 - mu, (lam l2 l1; mu m2 -m1) (-1)^m1
	vec salted_w3j(const int lam, const ivec2& llvec)
	{
		vec w;
		for (size_t il = 0; il < llvec[0].size(); il++)
			for (int mu = -lam; mu <= lam; mu++)
				for (int m1 = -llvec[0][il]; m1 <= llvec[0][il]; m1++)
					if (std::abs(m1 - mu) <= llvec[1][il])
						w.push_back(wigner_3j(lam, llvec[1][il], llvec[0][il], mu, m1 - mu, -m1) * (m1 % 2 ? -1.0 : 1.0));
		return w;
	}
}

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
	vec conj_p(2, 0.0);
	equicomb(1, 2, 1, v1, v1, w3j, llvec, 0, c2r, 2, 2, std::vector<int64_t>{ 0, 1 }, conj_p, true);
	EXPECT_NEAR(conj_p[0], 0.6, 1e-14);
	EXPECT_NEAR(conj_p[1], 0.8, 1e-14);
	// complex a, which no real density has, in the dense path only: without conjugation a a is real and negative
	v1.block(0, 0, 0)[0] = cdouble(0.0, 3.0);
	vec plain(2, 0.0);
	equicomb(1, 2, 1, v1, v1, w3j, 1, llvec, 0, c2r, 2, plain, false);
	EXPECT_NEAR(plain[0], -1.0, 1e-14);
	EXPECT_NEAR(plain[1], 0.0, 1e-14);
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
	make_real_density(v1, natoms);
	make_real_density(v2, natoms);
	ivec2 llvec{ { 1, 0 }, { 0, 1 } };
	vec w3j = salted_w3j(1, llvec);
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
	// a real density: c = -conj(a), b and d real
	const cdouble a(0.3, -0.7), b(-0.4, 0.0), c(-0.3, -0.7), d(0.6, 0.0);
	v.block(0, 0, 1)[0] = a;
	v.block(0, 0, 1)[1] = b;
	v.block(0, 0, 1)[2] = c;
	v.block(0, 0, 0)[0] = d;
	ivec2 llvec{ { 1 }, { 0 } };
	vec w3j = salted_w3j(1, llvec);
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
	++constants::run_id;
	testing::internal::CaptureStdout();
	const vec2 norms = equicomb_norms(1, 1, 1, v, v, { &w3j }, { llvec }, { c2r }, true);
	EXPECT_NE(testing::internal::GetCapturedStdout().find("1 atom(s) have no neighbour"), std::string::npos);
	EXPECT_EQ(norms[0][0], 0.0);
	EXPECT_TRUE(equicomb_norms(0, 1, 1, v, v, { &w3j }, { llvec }, { c2r }, true)[0].empty());
}

// lam = 2 over every (l1, l2) up to 3 that couples to it with l1 + l2 + lam even, real densities, unequal channel counts: the sparse path builds only the
// selected features and takes the norm from per-l density matrices, so a subset in shuffled order, with slots past
// nrad1*nrad2*llmax, must still equal the dense path that builds and normalises all of them. Both conj forms, and
// the device when one is present (it falls back to the CPU otherwise)
TEST(SaltedFchkEquicombTests, SelectedFeaturesMatchDenseLam2)
{
	const int natoms = 3, nrad1 = 3, nrad2 = 2, lam = 2, l21 = 5, lmax = 3;
	SALTEDDescriptors v1(natoms, nrad1, lmax), v2(natoms, nrad2, lmax);
	double seed = 0.61;
	auto next = [&seed]()
	{
		seed = std::fmod(seed * 7.31 + 0.113, 1.0);
		return seed - 0.5;
	};
	for (auto& x : v1.values())
		x = cdouble(next(), next());
	for (auto& x : v2.values())
		x = cdouble(next(), next());
	make_real_density(v1, natoms);
	make_real_density(v2, natoms);
	// every lambda up to lam, for equicomb_norms' one pass over all of them
	std::vector<ivec2> llvec_all(lam + 1, ivec2(2));
	std::vector<vec> w3j_k(lam + 1);
	std::vector<const vec*> w3j_all(lam + 1);
	for (int k = 0; k <= lam; k++)
	{
		for (int l1 = 0; l1 <= lmax; l1++)
			for (int l2 = 0; l2 <= lmax; l2++)
				if (std::abs(l1 - l2) <= k && k <= l1 + l2 && (l1 + l2 + k) % 2 == 0)
				{
					llvec_all[k][0].push_back(l1);
					llvec_all[k][1].push_back(l2);
				}
		w3j_k[k] = salted_w3j(k, llvec_all[k]);
		w3j_all[k] = &w3j_k[k];
	}
	ivec2& llvec = llvec_all[lam];
	const int llmax = static_cast<int>(llvec[0].size());
	vec& w3j = w3j_k[lam];
	const std::vector<cvec2> c2r_all = SALTED_Utils::complex_to_real_transformation({ 1, 3, l21 });
	cvec2 c2r = c2r_all[lam];
	const bool gpu_before = equicomb_gpu_enabled();
	for (const bool conj : { false, true })
	{
		const int n2 = conj ? nrad1 : nrad2;
		const SALTEDDescriptors& u = conj ? v1 : v2;
		const int shells = nrad1 * n2 * llmax, featsize = shells + 4;
		vec dense(static_cast<size_t>(natoms) * l21 * featsize, 0.0);
		equicomb_set_gpu(false);
		equicomb(natoms, nrad1, n2, v1, u, w3j, llmax, llvec, lam, c2r, featsize, dense, conj);
		std::vector<int64_t> vfps;
		for (int f = shells + 3; f >= 0; f -= 3)
			vfps.push_back(f);
		const int nfps = static_cast<int>(vfps.size());
		// 0: CPU, norm computed inside; 1: CPU, norm from equicomb_norms; 2: GPU
		const vec2 norms = equicomb_norms(natoms, nrad1, n2, v1, u, w3j_all, llvec_all, c2r_all, conj);
		for (const int mode : { 0, 1, 2 })
		{
			const bool gpu = mode == 2;
			equicomb_set_gpu(gpu);
			vec sub(static_cast<size_t>(natoms) * l21 * nfps, 7.0);
			equicomb(natoms, nrad1, n2, v1, u, w3j, llvec, lam, c2r, featsize, nfps, vfps, sub, conj,
				mode == 1 ? norms[lam].data() : nullptr);
			for (int iat = 0; iat < natoms; iat++)
				for (int imu = 0; imu < l21; imu++)
					for (int i = 0; i < nfps; i++)
					{
						const double want = vfps[i] < shells ? dense[iat * l21 * featsize + imu * featsize + vfps[i]] : 0.0;
						EXPECT_NEAR(sub[iat * l21 * nfps + imu * nfps + i], want, 1e-13) << conj << mode << iat << imu << i;
					}
		}
	}
	equicomb_set_gpu(gpu_before);
}

// NaCl rocksalt atoms sit on cubic inversion centres, where the lam = 1, 2, 3 features vanish. Displaced by d their
// squared norm is a small fraction of the bound sum_il tr A1[l1] tr A2[l2] / (2 lam + 1), while the density-matrix norm
// is a sum of terms of the size of the bound, so it cancels. For every d > 0 the norm from equicomb_norms must be the
// sum of the squared features to 1e-12, and the sparse features, norm made inside or by equicomb_norms, must stay on the
// dense path to the rounding of features that small
TEST(SaltedFchkEquicombTests, NearSymmetricSitesKeepTheNorm)
{
	EXPECT_NEAR(wigner_3j(1, 1, 0, 0, 0, 0), -1.0 / std::sqrt(3.0), 1e-15);
	EXPECT_NEAR(wigner_3j(2, 2, 0, 0, 0, 0), 1.0 / std::sqrt(5.0), 1e-15);
	EXPECT_NEAR(wigner_3j(1, 1, 2, 0, 0, 0), std::sqrt(2.0 / 15.0), 1e-15);
	EXPECT_NEAR(wigner_3j(1, 1, 1, 1, -1, 0), 1.0 / std::sqrt(6.0), 1e-15);
	const bool gpu_before = equicomb_gpu_enabled();
	equicomb_set_gpu(false);
	const double a = 5.64;
	const int natoms = 8, nrad = 4, nang = 6, nchan = 2 * nrad, lam_max = 4;
	const double frac[natoms][3] = { { 0, 0, 0 }, { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 }, { 0.5, 0.5, 0 },
		{ 0.5, 0.5, 0.5 }, { 0.5, 0, 0 }, { 0, 0.5, 0 }, { 0, 0, 0.5 } };
	// unit rms per component, as ase's rattle(stdev = d)
	std::mt19937 gen(5);
	double u[natoms][3];
	for (int i = 0; i < natoms; i++)
		for (int x = 0; x < 3; x++)
			u[i][x] = (gen() / 4294967296.0 - 0.5) * std::sqrt(12.0);
	SALTED_Utils::FeatomicHyperParameters hp;
	hp.cutoff_radius = 4.0;
	hp.max_radial = nrad;
	hp.max_angular = nang;
	hp.atomic_gaussian_width = 0.3;
	hp.species = { "Na", "Cl" };
	hp.neighspe = { "Na", "Cl" };
	std::vector<ivec2> llvec_all(lam_max + 1);
	std::vector<cvec2> c2r_all(lam_max + 1);
	std::vector<vec> w3j_all(lam_max + 1);
	std::vector<const vec*> w3j_ptr(lam_max + 1);
	for (int lam = 0; lam <= lam_max; lam++)
	{
		ivec2 llvec;
		for (int l1 = 0; l1 <= nang; l1++)
			for (int l2 = 0; l2 <= nang; l2++)
				if ((lam + l1 + l2) % 2 == 0 && std::abs(l2 - lam) <= l1 && l1 <= l2 + lam)
					llvec.push_back({ l1, l2 });
		llvec_all[lam] = transpose<int>(llvec);
		c2r_all[lam] = SALTED_Utils::complex_to_real_transformation({ 2 * lam + 1 })[0];
		w3j_all[lam] = salted_w3j(lam, llvec_all[lam]);
		w3j_ptr[lam] = &w3j_all[lam];
	}
	std::printf("     d/A  lam  dev(norm inside)  dev(equicomb_norms)  |norm2 err|/norm2  min norm2/bound\n");
	for (const double d : { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0.0 })
	{
		featomic::SimpleSystem sys({ { { a, 0, 0 }, { 0, a, 0 }, { 0, 0, a } } });
		for (int i = 0; i < natoms; i++)
			sys.add_atom(i < 4 ? 11 : 17, { (frac[i][0] + 0.125) * a + d * u[i][0], (frac[i][1] + 0.125) * a + d * u[i][1], (frac[i][2] + 0.125) * a + d * u[i][2] });
		const SALTEDDescriptors v = SALTED_Utils::calculate_SALTED_descriptors(sys, hp);
		ASSERT_EQ(v.nchannels(), nchan);
		ASSERT_EQ(v.values().size(), static_cast<size_t>(natoms) * nchan * (nang + 1) * (nang + 1));
		vec2 tr(natoms, vec(nang + 1, 0.0));
		for (int i = 0; i < natoms; i++)
			for (int l = 0; l <= nang; l++)
				for (int n = 0; n < nchan; n++)
					for (int m = 0; m < 2 * l + 1; m++)
						tr[i][l] += std::norm(v.block(i, n, l)[m]);
		const vec2 norms = equicomb_norms(natoms, nchan, nchan, v, v, w3j_ptr, llvec_all, c2r_all, true);
		const vec ones(natoms, 1.0);
		for (int lam = 0; lam <= lam_max; lam++)
		{
			const int l21 = 2 * lam + 1, llmax = static_cast<int>(llvec_all[lam][0].size()), featsize = nchan * nchan * llmax, block = l21 * featsize;
			std::vector<int64_t> identity(featsize);
			for (int f = 0; f < featsize; f++)
				identity[f] = f;
			const size_t size = static_cast<size_t>(natoms) * block;
			vec dense(size, 0.0), inside(size, 0.0), given(size, 0.0), raw(size, 0.0);
			equicomb(natoms, nchan, nchan, v, v, w3j_all[lam], llmax, llvec_all[lam], lam, c2r_all[lam], featsize, dense, true);
			equicomb(natoms, nchan, nchan, v, v, w3j_all[lam], llvec_all[lam], lam, c2r_all[lam], featsize, featsize, identity, inside, true);
			equicomb(natoms, nchan, nchan, v, v, w3j_all[lam], llvec_all[lam], lam, c2r_all[lam], featsize, featsize, identity, given, true, norms[lam].data());
			equicomb(natoms, nchan, nchan, v, v, w3j_all[lam], llvec_all[lam], lam, c2r_all[lam], featsize, featsize, identity, raw, true, ones.data());
			// the device when one is present, the CPU otherwise
			vec device(size, 0.0);
			equicomb_set_gpu(true);
			equicomb(natoms, nchan, nchan, v, v, w3j_all[lam], llvec_all[lam], lam, c2r_all[lam], featsize, featsize, identity, device, true, norms[lam].data());
			equicomb_set_gpu(false);
			double dev_inside = 0.0, dev_given = 0.0, norm_err = 0.0, ratio = 1e300;
			for (int i = 0; i < natoms; i++)
			{
				double top = 0.0, di = 0.0, dg = 0.0, dd = 0.0, direct = 0.0, bound = 0.0;
				for (int f = i * block; f < (i + 1) * block; f++)
				{
					top = std::max(top, std::abs(dense[f]));
					di = std::max(di, std::abs(inside[f] - dense[f]));
					dg = std::max(dg, std::abs(given[f] - dense[f]));
					dd = std::max(dd, std::abs(device[f] - dense[f]));
					direct += raw[f] * raw[f];
				}
				for (int il = 0; il < llmax; il++)
					bound += tr[i][llvec_all[lam][0][il]] * tr[i][llvec_all[lam][1][il]] / l21;
				dev_inside = std::max(dev_inside, di / top);
				dev_given = std::max(dev_given, dg / top);
				const double nf = norms[lam][i];
				const double err = nf > 0.0 ? std::abs(1.0 / (nf * nf) - direct) / direct : 1.0;
				norm_err = std::max(norm_err, err);
				ratio = std::min(ratio, direct / bound);
				// features a fraction r of the size of their terms carry rounding of about 1e-16 / sqrt(r) in either path
				if (d > 0.0)
				{
					EXPECT_LE(err, 1e-12) << "d " << d << " lam " << lam << " atom " << i;
					EXPECT_LE(di / top, 1e-14 / std::sqrt(direct / bound)) << "d " << d << " lam " << lam << " atom " << i;
					EXPECT_LE(dg / top, 1e-14 / std::sqrt(direct / bound)) << "d " << d << " lam " << lam << " atom " << i;
					EXPECT_LE(dd / top, 1e-14 / std::sqrt(direct / bound)) << "device, d " << d << " lam " << lam << " atom " << i;
				}
			}
			std::printf("%8.0e  %3d  %16.2e  %19.2e  %17.2e  %15.2e\n", d, lam, dev_inside, dev_given, norm_err, ratio);
		}
	}
	equicomb_set_gpu(gpu_before);
}

// CH4 (Td) and SF6 (Oh) at exact symmetry, at the origin and shifted by 100 A: the central atom's lambda 1, 2, 5 (Td) and
// 1, 2, 3, 5 (Oh) features vanish, and what is computed for them is rounding noise of the positions, which equicomb_norms
// with noise2 = (100 eps L)^2 zeroes. Every other atom and lambda keeps its norm, and so does the centre when one ligand moves
// by 1e-8 A (moving the centre alone changes its even-l coefficients only at second order, the shift adding u.grad rho)
TEST(SaltedFchkEquicombTests, ExactSymmetryZeroesTheVanishingFeatures)
{
	equicomb_set_gpu(false);
	const int nrad = 7, nang = 7, lam_max = 5;
	const double a = 1.09 / std::sqrt(3.0), s = 1.56;
	struct molecule { std::vector<std::string> el; std::vector<d3> xyz; std::vector<int> vanish; };
	const std::vector<molecule> mols{
		{ { "C", "H", "H", "H", "H" }, { { 0, 0, 0 }, { a, a, a }, { a, -a, -a }, { -a, a, -a }, { -a, -a, a } }, { 1, 2, 5 } },
		{ { "S", "F", "F", "F", "F", "F", "F" }, { { 0, 0, 0 }, { s, 0, 0 }, { -s, 0, 0 }, { 0, s, 0 }, { 0, -s, 0 }, { 0, 0, s }, { 0, 0, -s } }, { 1, 2, 3, 5 } } };
	for (const molecule& m : mols)
		for (const double shift : { 0.0, 100.0 })
			for (const double d : { 0.0, 1e-8 })
			{
				const int n = static_cast<int>(m.el.size());
				std::vector<std::string> species(m.el.begin(), m.el.begin() + 2);
				SALTED_Utils::FeatomicHyperParameters hp;
				hp.cutoff_radius = 4.0;
				hp.max_radial = nrad;
				hp.max_angular = nang;
				hp.atomic_gaussian_width = 0.3;
				hp.species = species;
				hp.neighspe = species;
				featomic::SimpleSystem sys;
				double L = 1.0;
				for (int i = 0; i < n; i++)
				{
					const double di = i == 1 ? d : 0.0;
					const d3 r{ m.xyz[i][0] + (shift + 0.8 * di) * 0.6, m.xyz[i][1] - (shift - 0.3 * di) * 0.45, m.xyz[i][2] + (shift + 0.5 * di) * 0.8 };
					sys.add_atom(constants::get_Z_from_label(m.el[i].c_str()) + 1, r);
					for (const double x : r)
						L = std::max(L, std::abs(x));
				}
				const SALTEDDescriptors v = SALTED_Utils::calculate_SALTED_descriptors(sys, hp);
				const int nchan = v.nchannels();
				std::vector<ivec2> llvec(lam_max + 1);
				std::vector<cvec2> c2r(lam_max + 1);
				std::vector<vec> w3j(lam_max + 1);
				std::vector<const vec*> w3j_ptr(lam_max + 1);
				for (int lam = 0; lam <= lam_max; lam++)
				{
					ivec2 ll;
					for (int l1 = 0; l1 <= nang; l1++)
						for (int l2 = 0; l2 <= nang; l2++)
							if ((lam + l1 + l2) % 2 == 0 && std::abs(l2 - lam) <= l1 && l1 <= l2 + lam)
								ll.push_back({ l1, l2 });
					llvec[lam] = transpose<int>(ll);
					c2r[lam] = SALTED_Utils::complex_to_real_transformation({ 2 * lam + 1 })[0];
					w3j[lam] = salted_w3j(lam, llvec[lam]);
					w3j_ptr[lam] = &w3j[lam];
				}
				const double noise = 100.0 * std::numeric_limits<double>::epsilon() * L;
				++constants::run_id;
				testing::internal::CaptureStdout();
				const vec2 norms = equicomb_norms(n, nchan, nchan, v, v, w3j_ptr, llvec, c2r, true, noise * noise);
				EXPECT_EQ(testing::internal::GetCapturedStdout().find("no neighbour"), std::string::npos) << m.el[0] << " shift " << shift << " d " << d;
				for (int lam = 0; lam <= lam_max; lam++)
					for (int i = 0; i < n; i++)
					{
						const bool zero = i == 0 && d == 0.0 && std::find(m.vanish.begin(), m.vanish.end(), lam) != m.vanish.end();
						EXPECT_EQ(norms[lam][i] == 0.0, zero) << m.el[0] << " shift " << shift << " d " << d << " lam " << lam << " atom " << i;
					}
				// with the centre's norm factor 0 its features are zero, so the kernel gives no coefficient for that lambda
				const int lam = 2, llmax = static_cast<int>(llvec[lam][0].size()), featsize = nchan * nchan * llmax;
				std::vector<int64_t> identity(featsize);
				for (int f = 0; f < featsize; f++)
					identity[f] = f;
				vec p(static_cast<size_t>(n) * 5 * featsize, 1.0);
				equicomb(n, nchan, nchan, v, v, w3j[lam], llvec[lam], lam, c2r[lam], featsize, featsize, identity, p, true, norms[lam].data());
				double centre = 0.0;
				for (int f = 0; f < 5 * featsize; f++)
					centre = std::max(centre, std::abs(p[f]));
				if (d == 0.0)
					EXPECT_EQ(centre, 0.0) << m.el[0] << " shift " << shift;
				else
					EXPECT_GT(centre, 1e-3) << m.el[0] << " shift " << shift;
			}
}

// rel_power against SALTED's own kernels (branch soap-norm-symmetric-sites, 61df639): CH4 with H1 moved 0.01 A along x,
// rcut 4, sigma 0.3, nrad 4, nang 4, no sparsification; carbon lambda 0-2 divided by its lambda = 0 length
TEST(SaltedFchkEquicombTests, RelPowerMatchesSalted)
{
	equicomb_set_gpu(false);
	const auto ref = nos_test_repo_root() / "tests" / "SALTED" / "relpower_ch4";
	const int nrad = 4, nang = 4;
	SALTED_Utils::FeatomicHyperParameters hp;
	hp.cutoff_radius = 4.0;
	hp.max_radial = nrad;
	hp.max_angular = nang;
	hp.atomic_gaussian_width = 0.3;
	hp.species = { "C", "H" };
	hp.neighspe = { "C", "H" };
	featomic::SimpleSystem sys;
	const double a = 0.629118;
	const std::vector<d3> xyz{ { 0, 0, 0 }, { a + 0.01, a, a }, { -a, -a, a }, { a, -a, -a }, { -a, a, -a } };
	for (int i = 0; i < 5; i++)
		sys.add_atom(constants::get_Z_from_label(i == 0 ? "C" : "H") + 1, xyz[i]);
	const SALTEDDescriptors v = SALTED_Utils::calculate_SALTED_descriptors(sys, hp);
	const int nchan = v.nchannels();
	ASSERT_EQ(nchan, 8);
	std::vector<ivec2> llvec(3);
	std::vector<cvec2> c2r(3);
	std::vector<vec> w3j(3);
	for (int lam = 0; lam <= 2; lam++)
	{
		ivec2 ll;
		for (int l1 = 0; l1 <= nang; l1++)
			for (int l2 = 0; l2 <= nang; l2++)
				if ((lam + l1 + l2) % 2 == 0 && std::abs(l2 - lam) <= l1 && l1 <= l2 + lam)
					ll.push_back({ l1, l2 });
		llvec[lam] = transpose<int>(ll);
		c2r[lam] = SALTED_Utils::complex_to_real_transformation({ 2 * lam + 1 })[0];
		w3j[lam] = salted_w3j(lam, llvec[lam]);
	}
	const vec norm0 = equicomb_norms(5, nchan, nchan, v, v, { &w3j[0] }, { llvec[0] }, { c2r[0] }, true)[0];
	vec raw0;
	std::filesystem::path rp = ref / "C_raw_lam0.npy";
	read_npy<double>(rp, raw0);
	double q0 = 0.0;
	for (const double x : raw0) q0 += x * x;
	EXPECT_NEAR(norm0[0] * std::sqrt(q0), 1.0, 1e-13);
	for (int lam = 0; lam <= 2; lam++)
	{
		const int llmax = static_cast<int>(llvec[lam][0].size()), featsize = nchan * nchan * llmax, l21 = 2 * lam + 1;
		vec p(5 * static_cast<size_t>(l21) * featsize);
		equicomb(5, nchan, nchan, v, v, w3j[lam], llmax, llvec[lam], lam, c2r[lam], featsize, p, true, norm0.data());
		vec want;
		std::filesystem::path wp = ref / ("C_relpower_lam" + std::to_string(lam) + ".npy");
		read_npy<double>(wp, want);
		ASSERT_EQ(want.size(), static_cast<size_t>(l21) * featsize);
		double top = 0.0, dev = 0.0;
		for (size_t i = 0; i < want.size(); i++)
		{
			top = std::max(top, std::abs(want[i]));
			dev = std::max(dev, std::abs(p[i] - want[i]));
		}
		EXPECT_LT(dev, 1e-12 * top) << "lambda " << lam;
	}
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

// the full prediction on the water monomer with the shipped model, which does not describe isolated water: the
// count is pinned, not compared to ten
TEST(SaltedFchkPredictorTests, PredictWaterMonomer)
{
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
	// each H gets 0.460 e against 0.807 e in the reference fit of this geometry (reading_SALTED/coefficients_conf0.npy,
	// same basis, 10.000 e); the Aug and Sep 2026 builds predict the same coefficients to 1e-16
	EXPECT_NEAR(e[0] + e[1] + e[2], 9.1937, 1e-3);
}

// -v accounts for the prediction: the payload line names exactly the bytes of the H and O blocks of the shipped model
// and the eager line those of the datasets read whole; without -v none of the stage lines is written
TEST(SaltedFchkPredictorTests, VerboseStagesCountThePayload)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	auto run = [&root](const bool verbose) {
		const bool before = ProgressBar::report_counts;
		ProgressBar::report_counts = verbose;
		WFN w(root / "reading_SALTED" / "water_monomer.xyz", false);
		options opt;
		opt.salted_model_dir = root / "SALTED" / "Model";
		++constants::run_id;
		testing::internal::CaptureStdout();
		SALTEDPredictor SP(w, opt);
		load_basis_into_WFN(SP.wavy, BasisSetLibrary::get_basis_set(SP.get_dfbasis_name()));
		SP.gen_SALTED_densities();
		const std::string out = testing::internal::GetCapturedStdout();
		ProgressBar::report_counts = before;
		return out;
	};
	SALTED_BINARY_FILE f(root / "SALTED" / "Model" / "model.salted");
	SALTEDConfig c{};
	f.populate_config(c);
	ASSERT_FALSE(f.has_block("FEATL") || f.has_block("ENVW") || f.has_block("GENV"));
	size_t bytes = 0, blocks = 0;
	for (const std::string key : { "FEATS", "PROJ" })
		for (const auto& [k, r] : f.index_lambda_based_data(key))
		{
			const std::string spe = k.substr(0, k.size() - 1);
			if (spe != "H" && spe != "O") continue;
			bytes += r.rows * r.cols * sizeof(double);
			blocks += r.rows * r.cols != 0;
		}
	f.read_weights();
	f.read_wigners();
	if (c.average) f.read_averages();
	if (c.sparsify) f.read_fps();
	const std::string out = run(true);
	for (const char* line : { "[stages] featomic calculator ", "[stages] featomic descriptors ", "[stages] predict ", " + model install ", " + assemble " })
		EXPECT_NE(out.find(line), std::string::npos) << line << " missing in\n" << out;
	std::smatch m;
	ASSERT_TRUE(std::regex_search(out, m, std::regex("\\[model\\] payload read \\S+ s, (\\d+) bytes in (\\d+) blocks; eager blocks \\S+ s, (\\d+) bytes"))) << out;
	EXPECT_EQ(std::stoull(m[1]), bytes);
	EXPECT_EQ(std::stoull(m[2]), blocks);
	EXPECT_EQ(std::stoull(m[3]), f.dataset_bytes);
	const std::string quiet = run(false);
	EXPECT_EQ(quiet.find("[stages]"), std::string::npos) << quiet;
	EXPECT_EQ(quiet.find("[model]"), std::string::npos) << quiet;
}

// Two parts sharing one model store (as -mtc does): water, then cysteine. The second reads only the blocks of the species
// water lacks, each block of the run comes from the file once, and both predict what a predictor of their own predicts
TEST(SaltedFchkPredictorTests, PartsShareOneModelLoad)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const std::filesystem::path xyz[2] = { root / "reading_SALTED" / "water_monomer.xyz", root / "SALTED" / "test_cysteine.xyz" };
	auto predict = [&root](const std::filesystem::path& p, const std::shared_ptr<SALTED_model_store>& s, std::string& out) {
		const bool before = ProgressBar::report_counts;
		ProgressBar::report_counts = true;
		WFN w(p, false);
		options opt;
		opt.salted_model_dir = root / "SALTED" / "Model";
		++constants::run_id;
		testing::internal::CaptureStdout();
		SALTEDPredictor SP(w, opt);
		SP.share_model(s);
		load_basis_into_WFN(SP.wavy, BasisSetLibrary::get_basis_set(SP.get_dfbasis_name()));
		const vec c = SP.gen_SALTED_densities();
		out = testing::internal::GetCapturedStdout();
		ProgressBar::report_counts = before;
		return c;
	};
	auto payload = [](const std::string& out) {
		std::smatch m;
		EXPECT_TRUE(std::regex_search(out, m, std::regex("\\[model\\] payload read \\S+ s, (\\d+) bytes in (\\d+) blocks"))) << out;
		return m.empty() ? std::array<size_t, 2>{} : std::array<size_t, 2>{ std::stoull(m[1]), std::stoull(m[2]) };
	};
	SALTED_BINARY_FILE f(root / "SALTED" / "Model" / "model.salted");
	std::map<std::string, std::array<size_t, 2>> per_species;
	for (const std::string key : { "FEATS", "PROJ" })
		for (const auto& [k, r] : f.index_lambda_based_data(key))
		{
			auto& b = per_species[k.substr(0, k.size() - 1)];
			b[0] += r.rows * r.cols * sizeof(double);
			b[1] += r.rows * r.cols != 0;
		}
	ASSERT_TRUE(per_species.count("C") && per_species.count("N"));
	const std::array<size_t, 2> water = { per_species["H"][0] + per_species["O"][0], per_species["H"][1] + per_species["O"][1] };
	std::array<size_t, 2> new_in_cys{};
	for (const std::string spe : { "C", "N", "S" })
		if (per_species.count(spe)) new_in_cys[0] += per_species[spe][0], new_in_cys[1] += per_species[spe][1];
	std::string out;
	const vec alone0 = predict(xyz[0], nullptr, out), alone1 = predict(xyz[1], nullptr, out);
	auto store = std::make_shared<SALTED_model_store>();
	const vec shared0 = predict(xyz[0], store, out);
	EXPECT_EQ(payload(out), water);
	const vec shared1 = predict(xyz[1], store, out);
	EXPECT_EQ(payload(out), new_in_cys);
	EXPECT_EQ(store->files.size(), 1u);
	EXPECT_EQ(store->files.begin()->second->lazy_blocks, water[1] + new_in_cys[1]);
	ASSERT_EQ(shared0.size(), alone0.size());
	ASSERT_EQ(shared1.size(), alone1.size());
	for (size_t i = 0; i < alone0.size(); i++) ASSERT_EQ(shared0[i], alone0[i]) << i;
	for (size_t i = 0; i < alone1.size(); i++) ASSERT_EQ(shared1[i], alone1[i]) << i;
}

// Every lambda's blocks are read by an OpenMP team of at most -cpus threads: reader threads of their own (a team of one
// each) run while the team's workers still spin, beyond -cpus
TEST(SaltedFchkPredictorTests, ModelIsReadInLineWithCpusReaders)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const int before = omp_get_max_threads();
	for (const int n : { 1, 4 })
	{
		omp_set_num_threads(n);
		WFN w(root / "reading_SALTED" / "water_monomer.xyz", false);
		options opt;
		opt.salted_model_dir = root / "SALTED" / "Model";
		++constants::run_id;
		auto store = std::make_shared<SALTED_model_store>();
		SALTEDPredictor SP(w, opt);
		SP.share_model(store);
		EXPECT_FALSE(SP.gen_SALTED_densities().empty());
		ASSERT_EQ(store->files.size(), 1u);
		const SALTED_BINARY_FILE& f = *store->files.begin()->second;
		EXPECT_GT(f.lazy_blocks, 0u);
		EXPECT_EQ(f.readers, (size_t)n);
		if (n == 1) EXPECT_EQ(f.team, 1);
		else EXPECT_TRUE(f.team > 1 && f.team <= n) << f.team;
	}
	omp_set_num_threads(before);
}

// A merged folder (the model and a copy of it as an alternative in a subfolder): the l = 0 blocks the SOAP score reads
// are kept for the prediction, so in a run without a store the predictions read no more blocks than in a run whose
// store holds every block (the score passes print no payload line), and predict the same
TEST(SaltedFchkPredictorTests, MergedScoreReadsLambda0Once)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const auto dir = tmp_path("merged");
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir / "alt");
	std::filesystem::copy_file(root / "SALTED" / "Model" / "model.salted", dir / "model.salted");
	std::filesystem::copy_file(root / "SALTED" / "Model" / "model.salted", dir / "alt" / "model.salted");
	auto predict = [&](const std::shared_ptr<SALTED_model_store>& s, size_t& blocks) {
		const bool before = ProgressBar::report_counts;
		ProgressBar::report_counts = true;
		WFN w(root / "reading_SALTED" / "water_monomer.xyz", false);
		options opt;
		opt.salted_model_dir = dir;
		++constants::run_id;
		testing::internal::CaptureStdout();
		SALTEDPredictor SP(w, opt);
		if (s) SP.share_model(s);
		const vec c = SP.gen_SALTED_densities();
		const std::string out = testing::internal::GetCapturedStdout();
		ProgressBar::report_counts = before;
		EXPECT_NE(out.find("Combining 2 SALTED models"), std::string::npos) << out;
		blocks = 0;
		const std::regex re("\\[model\\] payload read \\S+ s, \\d+ bytes in (\\d+) blocks");
		for (std::sregex_iterator m(out.begin(), out.end(), re), e; m != e; ++m) blocks += std::stoull((*m)[1]);
		return c;
	};
	size_t alone_blocks = 0, held_blocks = 0;
	const vec alone = predict(nullptr, alone_blocks);
	auto store = std::make_shared<SALTED_model_store>();
	const vec held = predict(store, held_blocks);
	EXPECT_EQ(store->files.size(), 2u);
	EXPECT_GT(alone_blocks, 0u);
	EXPECT_EQ(alone_blocks, held_blocks);
	ASSERT_EQ(alone.size(), held.size());
	for (size_t i = 0; i < alone.size(); i++) ASSERT_EQ(alone[i], held[i]) << i;
	store.reset();   // the store keeps both files open, and Windows deletes no open file
	std::filesystem::remove_all(dir);
}

// rel_power (NORM block): the lambda > 0 features are divided by the lambda = 0 length, so the coefficients of the
// carbon of CH4 move with one H by about as much at 1e-10 A as the displacement; unit-normalised its lambda 1 and 2
// blocks jump to full size at any displacement
TEST(SaltedFchkPredictorTests, RelPowerIsContinuousAtASymmetricSite)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const auto dir = tmp_path("relpower");
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir / "unit");
	std::filesystem::create_directories(dir / "rel");
	std::filesystem::copy_file(root / "SALTED" / "Model" / "model.salted", dir / "unit" / "model.salted");
	{
		const int32_t head[2] = { 4, 9 };
		std::string block(reinterpret_cast<const char*>(head), sizeof(head));
		block += "rel_power";
		SALTED_BINARY_FILE f(dir / "unit" / "model.salted");
		f.write_with_blocks(dir / "rel" / "model.salted", 7, {}, { { "NORM", block } });
	}
	{
		SALTEDConfig c{}, u{};
		SALTED_BINARY_FILE f(dir / "rel" / "model.salted"), g(dir / "unit" / "model.salted");
		f.populate_config(c);
		g.populate_config(u);
		EXPECT_EQ(c.normalisation, "rel_power");
		EXPECT_EQ(u.normalisation, "unit");
	}
	const double a = 1.09 / std::sqrt(3.0);
	auto predict = [&](const std::string& model, const double d) {
		const auto xyz = dir / "ch4.xyz";
		{
			std::ofstream o(xyz);
			o << std::setprecision(17) << "5\nCH4\nC 0 0 0\nH " << a + 0.3 * d << " " << a - 0.5 * d << " " << a + 0.8 * d
			  << "\nH " << a << " " << -a << " " << -a << "\nH " << -a << " " << a << " " << -a << "\nH " << -a << " " << -a << " " << a << "\n";
		}
		WFN w(xyz, false);
		options opt;
		opt.salted_model_dir = dir / model;
		SALTEDPredictor SP(w, opt);
		load_basis_into_WFN(SP.wavy, BasisSetLibrary::get_basis_set(SP.get_dfbasis_name()));
		return SP.gen_SALTED_densities();
	};
	auto dist = [](const vec& x, const vec& y) {
		double s = 0.0;
		for (size_t i = 0; i < x.size(); i++) s += (x[i] - y[i]) * (x[i] - y[i]);
		return std::sqrt(s);
	};
	for (const std::string model : { "unit", "rel" })
	{
		const vec c0 = predict(model, 0.0), c1 = predict(model, 1e-10), c2 = predict(model, 0.05);
		ASSERT_EQ(c0.size(), c2.size());
		const double step = dist(c1, c0) / dist(c2, c0);
		if (model == "unit")
			EXPECT_GT(step, 0.1);
		else
			EXPECT_LT(step, 1e-6);
	}
	std::filesystem::remove_all(dir);
}

// a .wfx names its atoms "C1", "H2": the species come from the nuclear charge, so it predicts what the same atoms
// named by element do. Read from the label, no atom matched a species and the model predicted nothing
TEST(SaltedFchkPredictorTests, WfxLabelsPredictLikeElements)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	auto predict = [&root](const WFN& w) {
		options opt;
		opt.salted_model_dir = root / "SALTED" / "Model";
		SALTEDPredictor SP(w, opt);
		load_basis_into_WFN(SP.wavy, BasisSetLibrary::get_basis_set(SP.get_dfbasis_name()));
		const vec coefs = SP.gen_SALTED_densities();
		EXPECT_EQ(coefs.size(), static_cast<size_t>(aux_density_table(SP.wavy.get_atoms()).n_coef));
		return coefs;
	};
	WFN wfx(root / "eqc_ch3f_wfx" / "ch3p.wfx", false);
	ASSERT_EQ(wfx.get_ncen(), 4);
	ASSERT_EQ(wfx.get_atom_label(0), "C1");
	WFN named = wfx;
	for (int a = 0; a < named.get_ncen(); a++)
		named.set_atom_label(a, constants::atnr2letter(named.get_atom_charge(a)));
	const vec a = predict(wfx), b = predict(named);
	ASSERT_FALSE(a.empty());
	ASSERT_EQ(a.size(), b.size());
	double big = 0.0;
	for (const double x : b) big = std::max(big, std::abs(x));
	for (size_t i = 0; i < a.size(); i++)
		EXPECT_NEAR(a[i], b[i], 1e-12 * big) << "coefficient " << i;
}

// -salted_fold: the folded model predicts what the unfolded one does, through ENVW for the shipped
// zeta = 1 and, in a copy with zeta set to 2, through GENV from the lambda on where it is no larger
TEST(SaltedFchkPredictorTests, FoldedModelPredictsTheSame)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	auto predict = [&root](const std::filesystem::path& dir) {
		WFN w(root / "reading_SALTED" / "water_monomer.xyz", false);
		options opt;
		opt.salted_model_dir = dir;
		SALTEDPredictor SP(w, opt);
		load_basis_into_WFN(SP.wavy, BasisSetLibrary::get_basis_set(SP.get_dfbasis_name()));
		return SP.gen_SALTED_densities();
	};
	const auto dir = tmp_path("fold");
	for (const double zeta : { 1.0, 2.0 })
	{
		std::filesystem::remove_all(dir);
		std::filesystem::create_directories(dir / "plain");
		std::filesystem::create_directories(dir / "folded");
		const auto plain = dir / "plain" / "model.salted", folded = dir / "folded" / "model.salted";
		std::filesystem::copy_file(root / "SALTED" / "Model" / "model.salted", plain);
		if (zeta != 1.0)
		{
			// CONFG: the value follows its 5-byte tag and a 4-byte type word
			std::fstream f(plain, std::ios::in | std::ios::out | std::ios::binary);
			std::string head(4096, '\0');
			f.read(head.data(), head.size());
			const size_t at = head.find("zeta");
			ASSERT_NE(at, std::string::npos);
			f.clear();
			f.seekp(at + 9);
			f.write(reinterpret_cast<const char*>(&zeta), sizeof(zeta));
		}
		const int threads = omp_get_max_threads();
		fold_salted_file(plain, folded);
		MKL_Set_Num_Threads(threads);   // the converter leaves BLAS on one thread
		{
			SALTED_BINARY_FILE f(folded);
			EXPECT_FALSE(f.has_block("FEATS"));
			EXPECT_EQ(f.has_block("ENVW"), zeta == 1.0);
			EXPECT_EQ(f.has_block("GENV"), zeta != 1.0);
		}
		const vec a = predict(dir / "plain"), b = predict(dir / "folded");
		ASSERT_EQ(a.size(), b.size());
		double big = 0.0;
		for (const double x : a) big = std::max(big, std::abs(x));
		for (size_t i = 0; i < a.size(); i++)
			EXPECT_NEAR(b[i], a[i], 1e-10 * big) << "zeta " << zeta << ", coefficient " << i;
	}
	std::filesystem::remove_all(dir);
}

// weights left over mean the BASIS block is not the set the model was trained with (a contracted
// cc-pvtz-jkfit packed for an uncontracted model): stop instead of misreading them
TEST(SaltedFchkPredictorTests, BasisNotFittingWeightsExits)
{
	const auto extra = tmp_path("extra_weight.salted");
	{
		SALTED_BINARY_FILE f(nos_test_repo_root() / "tests" / "SALTED" / "Model" / "model.salted");
		vec w = f.read_weights();
		w.push_back(0.0);
		salted_writer b;
		b.block_head(1);
		b.dataset(w, { static_cast<uint32_t>(w.size()) });
		f.write_with_blocks(extra, 2, { "WEIGH" }, { { "WEIGH", b.buf } });
	}
	EXPECT_EXIT(fold_salted_file(extra, tmp_path("extra_weight_folded.salted")), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	std::filesystem::remove(extra);
}

namespace
{
	// CONFG: the value follows its 5-byte tag and a 4-byte type word
	void set_salted_zeta(const std::filesystem::path& file, const double zeta)
	{
		std::fstream f(file, std::ios::in | std::ios::out | std::ios::binary);
		std::string head(4096, '\0');
		f.read(head.data(), head.size());
		const size_t at = head.find("zeta");
		ASSERT_NE(at, std::string::npos);
		f.clear();
		f.seekp(at + 9);
		f.write(reinterpret_cast<const char*>(&zeta), sizeof(zeta));
	}
	std::vector<salted_novelty> novelty_of(const std::filesystem::path& xyz, const std::filesystem::path& model, const double eigcut = 1e-10)
	{
		WFN w(xyz, false);
		options opt;
		opt.salted_model_dir = model;
		SALTEDPredictor SP(w, opt);
		return SP.kernel_novelty(w, eigcut);
	}
	// kernel_novelty_salted.py (SALTED fork-minimise-speed in salted_new.sif, Slurm job 1827715) on the shipped model and on
	// its copy with zeta = 2, lambda 0 read from FEATS where the script reads FEATL, in input order: kmax0, row of kmax0, nov0;
	// nov0 at eigcut 1e-2 from job 1828598, where K_MM keeps 39 of 98 H eigenvalues at zeta 1 and 55 at zeta 2
	struct novelty_ref
	{
		std::string xyz;
		vec kmax0;
		ivec nearest;
		vec nov_z1, nov_z2, cut_z1, cut_z2;
	};
	const std::vector<novelty_ref> novelty_refs{
		{ "SALTED/test_cysteine.xyz",
			{ 0.754045478178509, 0.774495929367651, 0.739953488699972, 0.645008575763813, 0.642696179320275, 0.655477621548586, 0.659872654593108,
			  0.677864282273147, 0.726285576508273, 0.734167059646307, 0.777846061933896, 0.732884044430286, 0.706836571086824, 0.717942810068487 },
			{ 3, 3, 3, 0, 0, 3, 4, 45, 30, 8, 20, 21, 52, 21 },
			{ 0.0117797375535706, 0.00453695737563131, 0.00870231169675484, 0.0140974756828659, 0.0136450478036562, 0.0112595798480108, 0.00898310975134875,
			  0.00227153525551615, 0.00177424163193673, 0.00192304699029833, 0.00187489131554963, 0.00136885396671171, 0.00195198563145882, 0.00100435685484646 },
			{ 0.0246668153831703, 0.0107230858551686, 0.0204508598398003, 0.0290194162101628, 0.0284509692653839, 0.0256681637554955, 0.0190661102894034,
			  0.00883340846900049, 0.00604142975824773, 0.00730866100164373, 0.0069412133254575, 0.00457980990759987, 0.00617739440626297, 0.00407103576622481 },
			{ 0.011940626883210181, 0.006138993432647011, 0.01058855057455077, 0.01425811370291219, 0.015607642989628356, 0.011403145399314973, 0.009412558448220865,
			  0.00449702981002964, 0.004750807783221789, 0.005324366268608993, 0.003394392392558787, 0.0034507106346224337, 0.0037741959876437248, 0.002621675021036163 },
			{ 0.024977741801844977, 0.012888844661895216, 0.022581311214535438, 0.029019416210162796, 0.028523416539094182, 0.025928586711795876, 0.019067724367184735,
			  0.011997258072159789, 0.009472759536513609, 0.01150040212299408, 0.009575878553352069, 0.007763792639124967, 0.009372900525734917, 0.006304757733920807 } },
		{ "reading_SALTED/water_monomer.xyz",
			{ 0.646376601897029, 0.792180180330442, 0.792180180330442 },
			{ 14, 65, 65 },
			{ 0.136930450574212, 0.016164218956451, 0.016164218956451 },
			{ 0.299485602667512, 0.0615321490271709, 0.0615321490271709 },
			{ 0.14806827511795106, 0.026896128672088127, 0.026896128672088127 },
			{ 0.30654672016102014, 0.0683310039578382, 0.0683310039578382 } } };
}

// -SALTED_distance: kmax0, its reference row and nov0 of every atom as SALTED's own descriptors give them, for zeta 1 and 2;
// the zeta = 2 copy folded (FEATL) gives the same numbers as unfolded (FEATS)
TEST(SaltedFchkPredictorTests, KernelNoveltyMatchesSaltedScript)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const auto dir = tmp_path("novelty");
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir);
	std::filesystem::copy_file(root / "SALTED" / "Model" / "model.salted", dir / "z2.salted");
	set_salted_zeta(dir / "z2.salted", 2.0);
	const int threads = omp_get_max_threads();
	fold_salted_file(dir / "z2.salted", dir / "z2_folded.salted");
	MKL_Set_Num_Threads(threads);
	{
		SALTED_BINARY_FILE f(dir / "z2_folded.salted");
		ASSERT_TRUE(f.has_block("FEATL") && !f.has_block("FEATS"));
	}
	for (const auto& ref : novelty_refs)
	{
		const auto z1 = novelty_of(root / ref.xyz, root / "SALTED" / "Model" / "model.salted");
		const auto z2 = novelty_of(root / ref.xyz, dir / "z2.salted"), z2f = novelty_of(root / ref.xyz, dir / "z2_folded.salted");
		ASSERT_EQ(z1.size(), ref.kmax0.size());
		ASSERT_EQ(z2.size(), ref.kmax0.size());
		ASSERT_EQ(z2f.size(), ref.kmax0.size());
		for (int a = 0; a < (int)z1.size(); a++)
		{
			for (const auto* r : { &z1[a], &z2[a], &z2f[a] })
			{
				EXPECT_EQ(r->atom, a);
				EXPECT_EQ(r->status, "ok");
				EXPECT_EQ(r->role, "assigned");
				EXPECT_NEAR(r->kmax0, ref.kmax0[a], 1e-8) << ref.xyz << " atom " << a;
				EXPECT_EQ(r->nearest, ref.nearest[a]) << ref.xyz << " atom " << a;
				EXPECT_GT(r->nov0, -1e-12);
				EXPECT_LT(r->nov0, 1.0 + 1e-12);
			}
			EXPECT_NEAR(z1[a].nov0, ref.nov_z1[a], 1e-8) << ref.xyz << " atom " << a;
			EXPECT_NEAR(z2[a].nov0, ref.nov_z2[a], 1e-8) << ref.xyz << " atom " << a;
			EXPECT_NEAR(z2f[a].nov0, z2[a].nov0, 1e-12) << ref.xyz << " atom " << a;
		}
	}
	std::filesystem::remove_all(dir);
}

// eigcut 1e-2 drops eigenvalues of K_MM for every species but N at zeta 2 (1e-4 drops none), nov0 grows by up to 1.1e-2
TEST(SaltedFchkPredictorTests, KernelNoveltyDropsEigenvaluesBelowEigcut)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const auto dir = tmp_path("novelty_cut");
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir);
	std::filesystem::copy_file(root / "SALTED" / "Model" / "model.salted", dir / "z2.salted");
	set_salted_zeta(dir / "z2.salted", 2.0);
	const int threads = omp_get_max_threads();
	fold_salted_file(dir / "z2.salted", dir / "z2_folded.salted");
	MKL_Set_Num_Threads(threads);
	for (const auto& ref : novelty_refs)
	{
		const auto z1 = novelty_of(root / ref.xyz, root / "SALTED" / "Model" / "model.salted", 1e-2);
		const auto z2 = novelty_of(root / ref.xyz, dir / "z2.salted", 1e-2), z2f = novelty_of(root / ref.xyz, dir / "z2_folded.salted", 1e-2);
		ASSERT_EQ(z1.size(), ref.cut_z1.size());
		ASSERT_EQ(z2.size(), ref.cut_z2.size());
		ASSERT_EQ(z2f.size(), ref.cut_z2.size());
		for (int a = 0; a < (int)z1.size(); a++)
		{
			EXPECT_NEAR(z1[a].nov0, ref.cut_z1[a], 1e-8) << ref.xyz << " atom " << a;
			EXPECT_NEAR(z2[a].nov0, ref.cut_z2[a], 1e-8) << ref.xyz << " atom " << a;
			EXPECT_NEAR(z2f[a].nov0, ref.cut_z2[a], 1e-8) << ref.xyz << " atom " << a;
		}
	}
	std::filesystem::remove_all(dir);
}

// the eigcut after the xyz reaches the novelty, and one that is no number stops
TEST(SaltedFchkPredictorTests, KernelNoveltyTakesTheEigcutFromTheFlag)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	auto args = [&root](const std::string& cut) {
		return svec{ "NoSpherA2", "-no_gpu_salted", "-SALTED_distance", (root / "SALTED" / "Model" / "model.salted").string(),
			(root / "reading_SALTED" / "water_monomer.xyz").string(), cut };
	};
	//the parser stops before anything is read or written, so no working directory of its own
	auto parse = [&args](const std::string& cut) {
		svec a = args(cut);
		std::vector<char*> argv;
		for (std::string& s : a) argv.push_back(s.data());
		int argc = static_cast<int>(argv.size());
		std::ostringstream log;
		options opt(argc, argv.data(), log);
		opt.digest_options();
	};
	EXPECT_EXIT(parse("1e-2x"), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	const auto dir = tmp_path("novelty_flag"), old = std::filesystem::current_path();
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir);
	std::filesystem::current_path(dir);
	svec a = args("1e-2");
	std::vector<char*> argv;
	for (std::string& s : a) argv.push_back(s.data());
	EXPECT_EQ(run_app(static_cast<int>(argv.size()), argv.data()), 0);
	std::filesystem::current_path(old);
	svec lines;
	{
		std::ifstream in(dir / "water_monomer_salted_distance.csv");
		for (std::string line; std::getline(in, line);) lines.push_back(line);
	}
	ASSERT_EQ(lines.size(), 4u);
	for (int a = 0; a < 3; a++)
	{
		svec f;
		std::stringstream s(lines[a + 1]);
		for (std::string x; std::getline(s, x, ',');) f.push_back(x);
		ASSERT_GE(f.size(), 8u) << lines[a + 1];
		EXPECT_NEAR(std::stod(f[7]), novelty_refs[1].cut_z1[a], 1e-8) << lines[a + 1];
	}
	std::filesystem::remove_all(dir);
}

// ENVW (folded at zeta = 1) keeps only V W^T F, no reference features: the novelty stops and says so
TEST(SaltedFchkPredictorTests, KernelNoveltyNeedsTheFeatures)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const auto folded = tmp_path("novelty_envw.salted");
	const int threads = omp_get_max_threads();
	fold_salted_file(root / "SALTED" / "Model" / "model.salted", folded);
	MKL_Set_Num_Threads(threads);
	EXPECT_EXIT(novelty_of(root / "reading_SALTED" / "water_monomer.xyz", folded), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	std::filesystem::remove(folded);
}

// an element the model lacks goes to the spherical fill, and so does an atom of a model element with nothing within rcut:
// their rows say filled and empty environment at their place in the input, the others keep their numbers; salted_distance
// writes the same rows to <xyz stem>_salted_distance.csv, the model quoted with forward slashes
TEST(SaltedFchkPredictorTests, KernelNoveltyListsFilledAndEmptyAtomsInPlace)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const auto dir = tmp_path("novelty_filled");
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir);
	{
		std::ofstream o(dir / "water_cl.xyz");
		o << "5\nwater, a far Cl and a far O\nO 0.066334 0.0 0.003701\nCl 20.0 0.0 0.0\nH -0.526383 -0.769327 -0.029366\nH -0.526383 0.769327 -0.029366\n"
		  << "O 0.0 20.0 0.0\n";
	}
	const auto old = std::filesystem::current_path();
	std::filesystem::current_path(dir);
	options opt;
	const auto rows = salted_distance(root / "SALTED" / "Model" / "model.salted", dir / "water_cl.xyz", 1e-10, opt);
	std::filesystem::current_path(old);
	ASSERT_EQ(rows.size(), 5u);
	const novelty_ref& w = novelty_refs[1];
	for (const auto& [a, el, status] : { std::make_tuple(1, "Cl", "filled"), std::make_tuple(4, "O", "empty environment") })
	{
		EXPECT_EQ(rows[a].atom, a);
		EXPECT_EQ(rows[a].element, el);
		EXPECT_EQ(rows[a].status, status);
		EXPECT_TRUE(rows[a].model.empty());
		EXPECT_TRUE(std::isnan(rows[a].nov0));
	}
	for (const int a : { 0, 2, 3 })
	{
		const int r = a == 0 ? 0 : a - 1;
		EXPECT_EQ(rows[a].atom, a);
		EXPECT_EQ(rows[a].status, "ok");
		EXPECT_NEAR(rows[a].kmax0, w.kmax0[r], 1e-8) << a;
		EXPECT_NEAR(rows[a].nov0, w.nov_z1[r], 1e-8) << a;
	}
	svec lines;
	{
		std::ifstream in(dir / "water_cl_salted_distance.csv");
		std::string line;
		while (std::getline(in, line)) lines.push_back(line);
	}
	ASSERT_EQ(lines.size(), 6u);
	EXPECT_EQ(lines[0], "atom,label,element,model,role,kmax0,cosine,nov0,status,nearest,p0p0");
	EXPECT_EQ(lines[2], "2,Cl,Cl,\"\",,,,,filled,,");
	EXPECT_EQ(lines[5], "5,O,O,\"\",,,,,empty environment,,");
	EXPECT_EQ(lines[1].rfind("1,O,O,\"", 0), 0u) << lines[1];
	EXPECT_NE(lines[1].find("/model.salted\",assigned,"), std::string::npos) << lines[1];
	EXPECT_EQ(lines[1].find('\\'), std::string::npos) << lines[1];
	EXPECT_NE(lines[1].find(",ok,14,"), std::string::npos) << lines[1];
	std::filesystem::remove_all(dir);
}

// a folder with a lead and an O.salted, and an O.salted in the subgroup alt: H belongs to the lead and O to O.salted, and
// by the first shell each atom goes to the model of its bonded partner's element, so O to the lead and H to O.salted, with
// alt/O.salted the alternative listed after it; all the same file, so all with the numbers of the model itself
TEST(SaltedFchkPredictorTests, KernelNoveltyReportsTheAlternatives)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const auto dir = tmp_path("novelty_folder");
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir / "alt");
	for (const auto& f : { dir / "model.salted", dir / "O.salted", dir / "alt" / "O.salted" })
		std::filesystem::copy_file(root / "SALTED" / "Model" / "model.salted", f);
	const auto rows = novelty_of(root / "reading_SALTED" / "water_monomer.xyz", dir);
	const novelty_ref& w = novelty_refs[1];
	//pins build_merged as it is: O going to the lead because H is the lead's element is a known bug, a fixed rule changes this row
	const std::vector<std::tuple<int, std::string, std::string, std::string>> want{
		{ 0, dir.filename().string(), "model.salted", "assigned" },
		{ 1, dir.filename().string(), "O.salted", "assigned" }, { 1, "alt", "O.salted", "alternative" },
		{ 2, dir.filename().string(), "O.salted", "assigned" }, { 2, "alt", "O.salted", "alternative" } };
	ASSERT_EQ(rows.size(), want.size());
	for (int i = 0; i < (int)rows.size(); i++)
	{
		const auto& [a, folder, file, role] = want[i];
		const std::filesystem::path m(rows[i].model);
		EXPECT_EQ(rows[i].atom, a);
		EXPECT_EQ(rows[i].role, role) << i;
		EXPECT_EQ(m.filename().string(), file) << i;
		EXPECT_EQ(m.parent_path().filename().string(), folder) << i;
		EXPECT_NEAR(rows[i].nov0, w.nov_z1[a], 1e-8) << i;
		EXPECT_NEAR(rows[i].kmax0, w.kmax0[a], 1e-8) << i;
	}
	std::filesystem::remove_all(dir);
}

// a lead, S.salted and H.salted on cysteine: S is bonded to C (lead) and H, Cb to S and H, so both are first-shell means of
// two models, the first bonded partner's model assigned and the other shared; every row has the numbers of the model itself
TEST(SaltedFchkPredictorTests, KernelNoveltyReportsSharedModels)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	const auto dir = tmp_path("novelty_shared");
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir);
	for (const auto& f : { dir / "model.salted", dir / "S.salted", dir / "H.salted" })
		std::filesystem::copy_file(root / "SALTED" / "Model" / "model.salted", f);
	const auto rows = novelty_of(root / "SALTED" / "test_cysteine.xyz", dir);
	const novelty_ref& c = novelty_refs[0];
	std::map<int, svec> roles;
	for (const auto& r : rows)
	{
		roles[r.atom].push_back(std::filesystem::path(r.model).filename().string() + " " + r.role);
		EXPECT_EQ(r.status, "ok") << r.atom;
		EXPECT_NEAR(r.nov0, c.nov_z1[r.atom], 1e-8) << r.atom;
		EXPECT_NEAR(r.kmax0, c.kmax0[r.atom], 1e-8) << r.atom;
	}
	ASSERT_EQ(roles.size(), c.kmax0.size());
	EXPECT_EQ(roles[0], svec({ "model.salted assigned", "H.salted shared" }));
	EXPECT_EQ(roles[6], svec({ "S.salted assigned", "H.salted shared" }));
	for (const auto& [a, v] : roles)
		if (a != 0 && a != 6) EXPECT_EQ(v.size(), 1u) << a;
	std::filesystem::remove_all(dir);
}

// the input's cores (-ECP, -mtc_ECP, an ECP gbw) belong to its basis: an all-electron model declares none, a model
// with an ECP block its own mode's table and nothing of the input's, a merged folder the cores of the models it uses
TEST(SaltedFchkPredictorTests, OnlyTheModelDeclaresACore)
{
	const auto root = nos_test_repo_root() / "tests" / "SALTED";
	const auto model = root / "Model" / "model.salted", dir = tmp_path("ecp_core");
	std::filesystem::remove_all(dir);
	std::filesystem::create_directories(dir);
	write_salted_ecp(model, dir / "def2.salted", 1);
	write_salted_ecp(model, dir / "xtb.salted", 2);
	WFN plain(root / "test_cysteine.xyz", false);
	WFN flagged = plain;
	flagged.set_has_ECPs(true, true, 2);
	auto cores = [](const WFN& w, const pathvec& models) {
		options opt;
		opt.salted_model_dir = models[0];
		if (models.size() > 1) opt.salted_model_dirs = models;
		SALTEDPredictor SP(w, opt);
		ivec c(SP.wavy.get_ncen());
		for (int a = 0; a < SP.wavy.get_ncen(); a++) c[a] = SP.wavy.get_atom_ECP_electrons(a);
		return std::make_pair(SP.wavy.get_has_ECPs() ? SP.wavy.get_ECP_mode() : 0, c);
	};
	const ivec none(plain.get_ncen(), 0);
	ivec xtb(plain.get_ncen());
	for (int a = 0; a < plain.get_ncen(); a++) xtb[a] = constants::ECP_core_electrons(constants::ECP_electrons_xTB, plain.get_atom_charge(a));
	ASSERT_EQ(vec_sum(xtb), 22);
	EXPECT_EQ(cores(flagged, { model }), std::make_pair(0, none));
	EXPECT_EQ(cores(flagged, { dir / "def2.salted" }), std::make_pair(1, none));
	EXPECT_EQ(cores(plain, { dir / "xtb.salted" }), std::make_pair(2, xtb));
	EXPECT_EQ(cores(flagged, { model, dir / "def2.salted" }), std::make_pair(0, none));
	EXPECT_EQ(cores(plain, { dir / "xtb.salted", model }), std::make_pair(2, xtb));
	std::filesystem::remove_all(dir);
}

namespace
{
	// run_app() in a directory of its own, so these runs share no log with the toml cases
	std::filesystem::path run_salted_cli(const std::string& name, svec args)
	{
		const auto dir = tmp_path(name), previous = std::filesystem::current_path();
		std::filesystem::remove_all(dir);
		std::filesystem::create_directories(dir);
		std::filesystem::current_path(dir);
		for (const char* d : { "-all_charges", "-no_date", "-no_gpu_grid", "-no_gpu_salted" }) args.push_back(d);
		args.insert(args.begin(), "NoSpherA2");
		std::vector<char*> argv;
		for (std::string& a : args) argv.push_back(a.data());
		EXPECT_EQ(run_app(static_cast<int>(argv.size()), argv.data()), 0) << name;
		std::filesystem::current_path(previous);
		return dir;
	}

	// the numbers of a log's charge tables and electron totals
	vec charge_numbers(const std::filesystem::path& log)
	{
		std::ifstream f(log);
		vec v;
		bool in = false;
		for (std::string l; std::getline(f, l);)
		{
			if (l.rfind("Table of Charges", 0) == 0) in = true;
			if (!in) continue;
			std::istringstream s(l);
			for (std::string t; s >> t;)
				if (std::isdigit(static_cast<unsigned char>(t.back())) && (std::isdigit(static_cast<unsigned char>(t[0])) || t[0] == '-'))
					v.push_back(std::stod(t));
			if (l.rfind("Total number of analytical Electrons", 0) == 0) in = false;
		}
		return v;
	}

	void expect_charges_of(const std::filesystem::path& good, const std::filesystem::path& log)
	{
		const vec g = charge_numbers(good), r = charge_numbers(log);
		ASSERT_FALSE(g.empty());
		ASSERT_EQ(g.size(), r.size()) << log;
		for (size_t i = 0; i < g.size(); i++) EXPECT_NEAR(r[i], g[i], 0.01) << log << " number " << i;
	}

	// a text .tsc's rows: h, k, l and the real and imaginary part of every scatterer
	vec2 tsc_rows(const std::filesystem::path& p)
	{
		std::ifstream f(p);
		vec2 rows;
		bool data = false;
		for (std::string l; std::getline(f, l);)
		{
			if (data)
			{
				std::replace(l.begin(), l.end(), ',', ' ');
				std::istringstream s(l);
				rows.emplace_back();
				for (double x; s >> x;) rows.back().push_back(x);
			}
			data = data || l.rfind("DATA:", 0) == 0;
		}
		return rows;
	}
}

// the toml cases SALTED_ECP2 and SALTED_ECP1: the all-electron model with -ECP gives the golden and the table of the
// run without it; the same model flagged as xTB valence adds the 22 core electrons of cysteine
TEST(SaltedFchkCliTests, EcpFlagAddsNoCoreToAnAllElectronModel)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests" / "SALTED";
	const svec base = { "-SALTED", (root / "Model").string(), "-cif", (root / "test_cysteine.cif").string(),
		"-wfn", (root / "test_cysteine.xyz").string(), "-dmin", "0.73", "-old_tsc" };
	auto with = [&base](const svec& extra) { svec a = base; a.insert(a.end(), extra.begin(), extra.end()); return a; };
	const auto none = run_salted_cli("cli_none", base), ecp2 = run_salted_cli("cli_ecp2", with({ "-ECP", "2" })),
		ecp1 = run_salted_cli("cli_ecp1", with({ "-ECP", "1" }));
	const auto xtb_model = tmp_path("cli_xtb_model");
	std::filesystem::remove_all(xtb_model);
	std::filesystem::create_directories(xtb_model);
	write_salted_ecp(root / "Model" / "model.salted", xtb_model / "model.salted", 2);
	svec xa = base;
	xa[1] = xtb_model.string();
	const auto xtb = run_salted_cli("cli_xtb", xa);
	const vec2 t0 = tsc_rows(none / "experimental.tsc"), tx = tsc_rows(xtb / "experimental.tsc");
	ASSERT_FALSE(t0.empty());
	for (const auto& d : { ecp2, ecp1 })
	{
		expect_charges_of(root / "SALTED.good", d / "NoSpherA2.log");
		const vec2 t = tsc_rows(d / "experimental.tsc");
		ASSERT_EQ(t.size(), t0.size()) << d;
		for (size_t r = 0; r < t.size(); r++)
		{
			ASSERT_EQ(t[r].size(), t0[r].size()) << d;
			for (size_t i = 0; i < t[r].size(); i++) EXPECT_NEAR(t[r][i], t0[r][i], 1e-8 * (1.0 + std::abs(t0[r][i]))) << d << " row " << r;
		}
	}
	const vec c0 = charge_numbers(none / "NoSpherA2.log"), cx = charge_numbers(xtb / "NoSpherA2.log");
	ASSERT_EQ(c0.size(), cx.size());
	EXPECT_NEAR(cx.back() - c0.back(), 22.0, 1e-3);
	// at the lowest reflection the core adds to the real part of every scatterer, nearly its electron count
	ASSERT_EQ(tx.size(), t0.size());
	auto hsum = [&t0](const size_t r) { return std::abs(t0[r][0]) + std::abs(t0[r][1]) + std::abs(t0[r][2]); };
	size_t low = 0;
	for (size_t r = 0; r < t0.size(); r++)
		if (hsum(r) < hsum(low)) low = r;
	double core = 0.0;
	for (size_t i = 3; i + 1 < t0[low].size(); i += 2)
	{
		EXPECT_GT(tx[low][i] - t0[low][i], -1e-8) << "scatterer " << (i - 3) / 2;
		EXPECT_NEAR(tx[low][i + 1], t0[low][i + 1], 1e-8);
		core += tx[low][i] - t0[low][i];
	}
	EXPECT_GT(core, 0.95 * 22.0);
	EXPECT_LT(core, 22.0 + 1e-6);
	for (const auto& d : { none, ecp2, ecp1, xtb, xtb_model }) std::filesystem::remove_all(d);
}

// the toml case SALTED_disorder_tshb_ECP2: -mtc_ECP 2 2 with the all-electron model gives the golden of -mtc_ECP 0 0
// and the same .tshb
TEST(SaltedFchkCliTests, MtcEcpAddsNoCoreToAnAllElectronModel)
{
	equicomb_set_gpu(false);
	const auto root = nos_test_repo_root() / "tests";
	auto args = [&root](const std::string& ecp) {
		return svec{ "-SALTED", (root / "SALTED" / "Model").string(), "-cif", (root / "Hybrid" / "ZP2.cif").string(),
			"-dmin", "0.9", "-tshb", "-mtc", (root / "Hybrid" / "ZP2_part1.gbw").string(), "0.1",
			(root / "Hybrid" / "ZP2_part2.gbw").string(), "0.2", "-mtc_mult", "1", "1", "-mtc_charge", "0", "0",
			"-mtc_ECP", ecp, ecp };
	};
	const auto zero = run_salted_cli("cli_mtc0", args("0")), two = run_salted_cli("cli_mtc2", args("2"));
	expect_charges_of(root / "Hybrid" / "SALTED_disorder_tshb.good", two / "NoSpherA2.log");
	auto bytes = [](const std::filesystem::path& p) { std::stringstream s; s << std::ifstream(p, std::ios::binary).rdbuf(); return s.str(); };
	const std::string a = bytes(zero / "experimental.tshb");
	EXPECT_FALSE(a.empty());
	EXPECT_TRUE(a == bytes(two / "experimental.tshb"));
	for (const auto& d : { zero, two }) std::filesystem::remove_all(d);
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

// the occ-built H2+ doublet written to fchk and read back keeps its density and spin density (the occupied orbital sits in the beta block)
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
	EXPECT_EQ(back.get_MO_op_count(1), 8);
	const d3 pos{ w.get_atom_coordinate(0, 0) + 0.3, w.get_atom_coordinate(0, 1), w.get_atom_coordinate(0, 2) };
	EXPECT_LT(w.compute_spin_dens(pos), 0.0);
	EXPECT_NEAR(back.compute_spin_dens(pos), w.compute_spin_dens(pos), 1e-6);
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
