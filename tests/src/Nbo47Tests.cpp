#include "pch.h"

#include "core/wfn_class.h"
#include "core/nbo.h"
#include "core/nbo_run.h"

#include <stdexcept>

namespace {

std::filesystem::path repo_root()
{
	return nos_test_repo_root();
}

std::string read_file(const std::filesystem::path& path)
{
	std::ifstream in(path);
	std::ostringstream out;
	out << in.rdbuf();
	return out.str();
}

vec extract_section_numbers(const std::string& text, const std::string& section)
{
	const auto start = text.find(section);
	if (start == std::string::npos) {
		return {};
	}
	const auto content_start = text.find('\n', start);
	const auto end = text.find("$END", content_start);
	if (content_start == std::string::npos || end == std::string::npos) {
		return {};
	}
	static const std::regex number_pattern(R"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)");
	const std::string content = text.substr(content_start, end - content_start);
	vec values;
	for (auto it = std::sregex_iterator(content.begin(), content.end(), number_pattern);
		 it != std::sregex_iterator();
		 ++it) {
		values.push_back(std::stod((*it).str()));
	}
	return values;
}

std::optional<int> parse_key_int(const std::string& text, const std::string& key)
{
	const std::regex pattern(key + R"(\s*=\s*(\d+))");
	std::smatch match;
	if (std::regex_search(text, match, pattern)) {
		return std::stoi(match[1].str());
	}
	return std::nullopt;
}

std::filesystem::path make_temp_dir()
{
	const auto parent = std::filesystem::temp_directory_path();
	const auto timestamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();

	// Each NBO test is registered as a separate CTest process.  They can run in
	// parallel, so a shared directory (and remove_all) lets one test delete the
	// .47 file another test has just produced.
	for (unsigned int attempt = 0; attempt != 100; ++attempt) {
		const auto directory = parent / ("nos_nbo47_test_" + std::to_string(timestamp) + "_" +
										 std::to_string(attempt));
		if (std::filesystem::create_directory(directory)) {
			return directory;
		}
	}

	throw std::runtime_error("Could not create a unique NBO test directory");
}

std::string windows_path_to_wsl(std::filesystem::path path)
{
	path = std::filesystem::absolute(path);
	std::string s = path.string();
	std::replace(s.begin(), s.end(), '\\', '/');
	if (s.size() >= 3 && s[1] == ':' && s[2] == '/') {
		const char drive = static_cast<char>(std::tolower(static_cast<unsigned char>(s[0])));
		s = std::string("/mnt/") + drive + s.substr(2);
	}
	return s;
}

bool wsl_gennbo_available()
{
	return std::system("wsl bash -lc \"test -x ~/nbo7/gennbo\"") == 0;
}

std::optional<double> parse_total_electrons(const std::filesystem::path& nbo_path)
{
	std::ifstream in(nbo_path);
	std::string line;
	const std::regex total_line(R"(\*\s+Total\s+\*\s+[-+]?\d*\.?\d+\s+[-+]?\d*\.?\d+\s+[-+]?\d*\.?\d+\s+[-+]?\d*\.?\d+\s+([-+]?\d*\.?\d+))");
	while (std::getline(in, line)) {
		std::smatch match;
		if (std::regex_search(line, match, total_line)) {
			return std::stod(match[1].str());
		}
	}
	return std::nullopt;
}

double packed_trace_product(const vec& density, const vec& overlap, int nbasis)
{
	double trace = 0.0;
	for (int i = 0; i < nbasis; i++) {
		for (int j = 0; j <= i; j++) {
			const int ij = i * (i + 1) / 2 + j;
			trace += density[ij] * overlap[ij] * (i == j ? 1.0 : 2.0);
		}
	}
	return trace;
}

} // namespace

TEST(Nbo47, EpoxideGbwWritesValidFile47)
{
	const auto root = repo_root();
	const auto fixture_dir = root / "tests" / "epoxide_gbw" / "NBO";
	const auto input_gbw = fixture_dir / "epoxide.gbw";
	ASSERT_TRUE(std::filesystem::exists(input_gbw));

	const auto temp_dir = make_temp_dir();
	const auto generated_47 = temp_dir / "epoxide.47";

	WFN wave(input_gbw, false);
	ASSERT_TRUE(wave.write_nbo(generated_47, false));
	ASSERT_TRUE(std::filesystem::exists(generated_47));

	const std::string text = read_file(generated_47);
	ASSERT_NE(text.find("$GENNBO"), std::string::npos);
	ASSERT_NE(text.find("$COORD"), std::string::npos);
	ASSERT_NE(text.find("$BASIS"), std::string::npos);
	ASSERT_NE(text.find("$CONTRACT"), std::string::npos);
	ASSERT_NE(text.find("$OVERLAP"), std::string::npos);
	ASSERT_NE(text.find("$DENSITY"), std::string::npos);
	ASSERT_NE(text.find("$FOCK"), std::string::npos);
	ASSERT_NE(text.find("$LCAOMO"), std::string::npos);

	EXPECT_EQ(parse_key_int(text, "NATOMS").value_or(-1), 7);
	EXPECT_EQ(parse_key_int(text, "NBAS").value_or(-1), 62);
	EXPECT_EQ(parse_key_int(text, "NSHELL").value_or(-1), 30);
	EXPECT_EQ(parse_key_int(text, "NEXP").value_or(-1), 56);

	const auto overlap = extract_section_numbers(text, "$OVERLAP");
	const auto density = extract_section_numbers(text, "$DENSITY");
	const auto fock = extract_section_numbers(text, "$FOCK");
	const auto lcaomo = extract_section_numbers(text, "$LCAOMO");
	const int nbasis = 62;
	EXPECT_EQ(overlap.size(), static_cast<size_t>(nbasis * (nbasis + 1) / 2));
	EXPECT_EQ(density.size(), static_cast<size_t>(nbasis * (nbasis + 1) / 2));
	EXPECT_EQ(fock.size(), static_cast<size_t>(nbasis * (nbasis + 1) / 2));
	EXPECT_EQ(lcaomo.size(), static_cast<size_t>(nbasis * nbasis));
	EXPECT_NEAR(packed_trace_product(density, overlap, nbasis), 24.0, 1.0e-5);
}

TEST(Nbo47, WriteNboReportsProgressWhenRequested)
{
	const auto root = repo_root();
	const auto fixture_dir = root / "tests" / "epoxide_gbw" / "NBO";
	const auto input_gbw = fixture_dir / "epoxide.gbw";
	ASSERT_TRUE(std::filesystem::exists(input_gbw));

	const auto temp_dir = make_temp_dir();
	const auto generated_47 = temp_dir / "epoxide.47";

	WFN wave(input_gbw, false);
	std::ostringstream progress_log;
	ASSERT_TRUE(wave.write_nbo(generated_47, false, &progress_log));

	const std::string log_text = progress_log.str();
	EXPECT_NE(log_text.find("[FILE47] Starting .47 conversion"), std::string::npos);
	EXPECT_NE(log_text.find("[FILE47] Computing spherical AO overlap integrals"), std::string::npos);
	EXPECT_NE(log_text.find("[FILE47] Writing FILE47 sections"), std::string::npos);
	EXPECT_NE(log_text.find("[FILE47] Finished .47 conversion"), std::string::npos);
}

TEST(Nbo47, OpenShellNh3LiGbwWritesValidFile47)
{
	const auto root = repo_root();
	const auto fixture_dir = root / "tests" / "RGBI_groups";
	const auto input_gbw = fixture_dir / "nh3li.gbw";
	ASSERT_TRUE(std::filesystem::exists(input_gbw));

	const auto temp_dir = make_temp_dir();
	const auto generated_47 = temp_dir / "nh3li.47";

	WFN wave(input_gbw, false);
	ASSERT_TRUE(wave.get_is_unrestricted());
	ASSERT_TRUE(wave.write_nbo(generated_47, false));
	ASSERT_TRUE(std::filesystem::exists(generated_47));

	const std::string text = read_file(generated_47);
	ASSERT_NE(text.find("$GENNBO"), std::string::npos);
	ASSERT_NE(text.find("$COORD"), std::string::npos);
	ASSERT_NE(text.find("$BASIS"), std::string::npos);
	ASSERT_NE(text.find("$CONTRACT"), std::string::npos);
	ASSERT_NE(text.find("$OVERLAP"), std::string::npos);
	ASSERT_NE(text.find("$DENSITY"), std::string::npos);
	ASSERT_NE(text.find("$FOCK"), std::string::npos);
	ASSERT_NE(text.find("$LCAOMO"), std::string::npos);
	//Without OPEN, NBO reads the archive as restricted and silently halves the electron
	//count it finds; the doubled blocks below are only meaningful together with it.
	ASSERT_NE(text.find(" OPEN "), std::string::npos);

	EXPECT_EQ(parse_key_int(text, "NATOMS").value_or(-1), 5);
	EXPECT_EQ(parse_key_int(text, "NBAS").value_or(-1), 63);
	EXPECT_EQ(parse_key_int(text, "NSHELL").value_or(-1), 31);
	EXPECT_EQ(parse_key_int(text, "NEXP").value_or(-1), 52);

	const auto overlap = extract_section_numbers(text, "$OVERLAP");
	const auto density = extract_section_numbers(text, "$DENSITY");
	const auto fock = extract_section_numbers(text, "$FOCK");
	const auto lcaomo = extract_section_numbers(text, "$LCAOMO");
	const int nbasis = 63;
	const size_t ntri = static_cast<size_t>(nbasis) * (nbasis + 1) / 2;
	//$OVERLAP stays single; $DENSITY, $FOCK and $LCAOMO carry an alpha block then a beta one.
	EXPECT_EQ(overlap.size(), ntri);
	ASSERT_EQ(density.size(), 2 * ntri);
	EXPECT_EQ(fock.size(), 2 * ntri);
	EXPECT_EQ(lcaomo.size(), 2 * static_cast<size_t>(nbasis) * nbasis);

	const vec alpha_density(density.begin(), density.begin() + ntri);
	const vec beta_density(density.begin() + ntri, density.end());
	//13 electrons in a doublet: 7 alpha, 6 beta. A spin-summed archive would give 13 here
	//twice, and a swapped one 6 then 7.
	EXPECT_NEAR(packed_trace_product(alpha_density, overlap, nbasis), 7.0, 1.0e-5);
	EXPECT_NEAR(packed_trace_product(beta_density, overlap, nbasis), 6.0, 1.0e-5);
}

TEST(Nbo47, OpenShellNh3LiGennboProducesEnergyAnalysisWhenAvailable)
{
	if (!wsl_gennbo_available()) {
		GTEST_SKIP() << "WSL ~/nbo7/gennbo is not available";
	}

	const auto root = repo_root();
	const auto fixture_dir = root / "tests" / "RGBI_groups";
	const auto input_gbw = fixture_dir / "nh3li.gbw";
	ASSERT_TRUE(std::filesystem::exists(input_gbw));

	const auto temp_dir = make_temp_dir();
	const auto generated_47 = temp_dir / "nh3li.47";
	const auto generated_nbo = temp_dir / "nh3li.nbo";

	WFN wave(input_gbw, false);
	ASSERT_TRUE(wave.get_is_unrestricted());
	ASSERT_TRUE(wave.write_nbo(generated_47, false));

	const std::string wsl_dir = windows_path_to_wsl(temp_dir);
	const std::string command = "wsl bash -lc \"cd '" + wsl_dir + "' && ~/nbo7/gennbo nh3li\"";
	ASSERT_EQ(std::system(command.c_str()), 0);
	ASSERT_TRUE(std::filesystem::exists(generated_nbo));

	const std::string generated_nbo_text = read_file(generated_nbo);
	EXPECT_EQ(generated_nbo_text.find("Basis functions are not in expected form"), std::string::npos);
	EXPECT_NE(generated_nbo_text.find("NAO Atom No lang   Type(AO)    Occupancy      Energy"), std::string::npos);
	EXPECT_NE(generated_nbo_text.find("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX"), std::string::npos);
	EXPECT_NE(generated_nbo_text.find("NHO DIRECTIONALITY AND BOND BENDING"), std::string::npos);

	const auto actual_electrons = parse_total_electrons(generated_nbo);
	ASSERT_TRUE(actual_electrons.has_value());
	EXPECT_NEAR(*actual_electrons, 13.0, 1.0e-5);
}

TEST(NboRun, OpenShellNh3LiSpinResolvedNpaMatchesOrcaSpinPopulations)
{
	if (!wsl_gennbo_available()) {
		GTEST_SKIP() << "WSL ~/nbo7/gennbo is not available";
	}

	const auto input_gbw = repo_root() / "tests" / "RGBI_groups" / "nh3li.gbw";
	ASSERT_TRUE(std::filesystem::exists(input_gbw));

	const auto temp_dir = make_temp_dir();
	const auto generated_47 = temp_dir / "nh3li.47";
	const auto generated_nbo = temp_dir / "nh3li.nbo";

	WFN wave(input_gbw, false);
	ASSERT_TRUE(wave.write_nbo(generated_47, false));
	const std::string command = "wsl bash -lc \"cd '" + windows_path_to_wsl(temp_dir) + "' && ~/nbo7/gennbo nh3li\"";
	ASSERT_EQ(std::system(command.c_str()), 0);
	ASSERT_TRUE(std::filesystem::exists(generated_nbo));

	const NboResults r = parse_nbo_output(generated_nbo);
	EXPECT_TRUE(r.open_shell);
	ASSERT_EQ(r.npa.size(), 5u);

	double spin_sum = 0.0, charge_sum = 0.0;
	for (const auto& a : r.npa) { spin_sum += a.spin_density; charge_sum += a.charge; ASSERT_TRUE(a.has_spin_density); }
	EXPECT_NEAR(spin_sum, 1.0, 1.0e-4);
	EXPECT_NEAR(charge_sum, 0.0, 1.0e-4);

	//Independent reference: ORCA 6.1.1 on the same wavefunction puts 0.99 (Mulliken) / 0.90
	//(Loewdin) of the unpaired electron on Li and leaves N slightly negative. NPA is a third
	//partitioning, so only the pattern is compared - but a spin-summed or spin-swapped
	//archive gets the pattern wrong, which is what this pins down.
	const NboAtomPopulation* li = nullptr;
	const NboAtomPopulation* n = nullptr;
	for (const auto& a : r.npa) { if (a.element == "Li") li = &a; if (a.element == "N") n = &a; }
	ASSERT_NE(li, nullptr);
	ASSERT_NE(n, nullptr);
	EXPECT_GT(li->spin_density, 0.80);
	EXPECT_LT(std::abs(n->spin_density), 0.15);
	EXPECT_LT(n->charge, 0.0);

	//The unrestricted analysis has to reach the spin-resolved NBO sections as well.
	bool alpha = false, beta = false;
	for (const auto& o : r.orbitals) { alpha |= o.spin == "alpha"; beta |= o.spin == "beta"; }
	EXPECT_TRUE(alpha);
	EXPECT_TRUE(beta);
	bool spin_e2 = false;
	for (const auto& e : r.e2) spin_e2 |= !e.spin.empty();
	EXPECT_TRUE(spin_e2);
}

TEST(Nbo47, GShellWavefunctionWritesFile47WithCorrectElectronCount)
{
	//A basis with g functions is what broke get_shell_start_in_primitives: its switch covered
	//s, p, d and f and added nothing for g, so every primitive index behind the first g shell
	//was short by 15 per g shell.  Fe.gbw's atom 2 then asked for its s shell and was handed a
	//g primitive, wrote past the end of a one-component buffer and aborted in the heap later.
	//Tr(P S) is what says the coefficients that came back are the right ones, not merely that
	//nothing crashed.
	const auto root = repo_root();
	const auto input_gbw = root / "tests" / "Fe_gbw" / "Fe.gbw";
	ASSERT_TRUE(std::filesystem::exists(input_gbw));

	const auto temp_dir = make_temp_dir();
	const auto generated_47 = temp_dir / "fe.47";

	WFN wave(input_gbw, false);
	int highest_shell = 0;
	for (int a = 0; a < wave.get_ncen(); a++)
		for (int s = 0; s < wave.get_atom_shell_count(a); s++)
			highest_shell = std::max(highest_shell, wave.get_shell_type(a, s));
	ASSERT_GE(highest_shell, 5) << "fixture no longer carries g functions";

	ASSERT_TRUE(wave.write_nbo(generated_47, false));
	ASSERT_TRUE(std::filesystem::exists(generated_47));

	const std::string text = read_file(generated_47);
	const int nbasis = parse_key_int(text, "NBAS").value_or(-1);
	ASSERT_GT(nbasis, 0);
	const size_t ntri = static_cast<size_t>(nbasis) * (nbasis + 1) / 2;
	const auto overlap = extract_section_numbers(text, "$OVERLAP");
	const auto density = extract_section_numbers(text, "$DENSITY");
	ASSERT_EQ(overlap.size(), ntri);
	const bool open_shell = density.size() == 2 * ntri;
	ASSERT_TRUE(open_shell || density.size() == ntri);

	double electrons = 0.0;
	for (size_t block = 0; block < density.size() / ntri; block++)
		electrons += packed_trace_product(
			vec(density.begin() + block * ntri, density.begin() + (block + 1) * ntri), overlap, nbasis);
	//Against the wavefunction's own occupations, not against the nuclear charges: this fixture
	//integrates to 128 while Z - charge gives 126, which is a question about the gbw charge field
	//and not about whether the archive reproduces the wavefunction it was written from.
	double occupied = 0.0;
	for (int m = 0; m < wave.get_nmo(); m++)
		occupied += wave.get_MO_occ(m);
	EXPECT_NEAR(electrons, occupied, 1.0e-3);
}
