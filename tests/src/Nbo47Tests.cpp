#include "pch.h"

#include "core/wfn_class.h"
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

vec parse_natural_charges(const std::filesystem::path& nbo_path)
{
	std::ifstream in(nbo_path);
	std::string line;
	bool in_summary = false;
	vec charges;
	const std::regex charge_line(R"(^\s+[A-Z][a-z]?\s+\d+\s+([-+]?\d*\.?\d+))");
	while (std::getline(in, line)) {
		if (line.find("Summary of Natural Population Analysis") != std::string::npos) {
			in_summary = true;
			continue;
		}
		if (!in_summary) {
			continue;
		}
		if (line.find("* Total *") != std::string::npos) {
			break;
		}
		std::smatch match;
		if (std::regex_search(line, match, charge_line)) {
			charges.push_back(std::stod(match[1].str()));
		}
	}
	return charges;
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

TEST(Nbo47, EpoxideGennboMatchesReferenceWhenAvailable)
{
	if (!wsl_gennbo_available()) {
		GTEST_SKIP() << "WSL ~/nbo7/gennbo is not available";
	}

	const auto root = repo_root();
	const auto fixture_dir = root / "tests" / "epoxide_gbw" / "NBO";
	const auto input_gbw = fixture_dir / "epoxide.gbw";
	const auto reference_nbo = fixture_dir / "reference.nbo";
	ASSERT_TRUE(std::filesystem::exists(input_gbw));
	if (!std::filesystem::exists(reference_nbo)) {
		GTEST_SKIP() << "NBO reference fixture is not available";
	}

	const auto temp_dir = make_temp_dir();
	const auto generated_47 = temp_dir / "epoxide.47";
	const auto generated_nbo = temp_dir / "epoxide.nbo";

	WFN wave(input_gbw, false);
	ASSERT_TRUE(wave.write_nbo(generated_47, false));

	const std::string wsl_dir = windows_path_to_wsl(temp_dir);
	const std::string command = "wsl bash -lc \"cd '" + wsl_dir + "' && ~/nbo7/gennbo epoxide\"";
	ASSERT_EQ(std::system(command.c_str()), 0);
	ASSERT_TRUE(std::filesystem::exists(generated_nbo));

	const std::string generated_nbo_text = read_file(generated_nbo);
	EXPECT_EQ(generated_nbo_text.find("Basis functions are not in expected form"), std::string::npos);
	EXPECT_NE(generated_nbo_text.find("NAO Atom No lang   Type(AO)    Occupancy      Energy"), std::string::npos);
	EXPECT_NE(generated_nbo_text.find("SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX"), std::string::npos);
	EXPECT_NE(generated_nbo_text.find("NHO DIRECTIONALITY AND BOND BENDING"), std::string::npos);

	const auto expected_charges = parse_natural_charges(reference_nbo);
	const auto actual_charges = parse_natural_charges(generated_nbo);
	ASSERT_EQ(actual_charges.size(), expected_charges.size());
	for (size_t i = 0; i < expected_charges.size(); i++) {
		EXPECT_NEAR(actual_charges[i], expected_charges[i], 2.0e-3);
	}

	const auto expected_electrons = parse_total_electrons(reference_nbo);
	const auto actual_electrons = parse_total_electrons(generated_nbo);
	ASSERT_TRUE(expected_electrons.has_value());
	ASSERT_TRUE(actual_electrons.has_value());
	EXPECT_NEAR(*actual_electrons, *expected_electrons, 1.0e-5);
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

TEST(NboRun, ParsesReferenceOutputOfTheEpoxideFixture)
{
	const auto reference_nbo = repo_root() / "tests" / "epoxide_gbw" / "NBO" / "reference.nbo";
	if (!std::filesystem::exists(reference_nbo)) {
		GTEST_SKIP() << "NBO reference fixture is not available";
	}

	const NboResults r = parse_nbo_output(reference_nbo);
	EXPECT_FALSE(r.open_shell);
	ASSERT_EQ(r.npa.size(), 7u);
	EXPECT_EQ(r.npa[0].element, "O");
	EXPECT_NEAR(r.npa[0].charge, -0.56088, 1.0e-5);
	EXPECT_NEAR(r.npa[0].total, 8.56088, 1.0e-5);
	double charge_sum = 0.0;
	for (const auto& a : r.npa) charge_sum += a.charge;
	EXPECT_NEAR(charge_sum, 0.0, 1.0e-4);

	ASSERT_FALSE(r.nao.empty());
	EXPECT_EQ(r.nao.front().type, "Cor");
	EXPECT_NEAR(r.nao.front().occupancy, 1.99997, 1.0e-5);

	//An NBO with a bond has two hybrids that add up to the whole orbital, and each hybrid's
	//s/p/d percentages add up to 100 - the two things a wrong parse gets wrong first.
	ASSERT_FALSE(r.orbitals.empty());
	const NboOrbital* bond = nullptr;
	for (const auto& o : r.orbitals) if (o.type == "BD" && o.centers.size() == 2) { bond = &o; break; }
	ASSERT_NE(bond, nullptr);
	ASSERT_EQ(bond->hybrids.size(), 2u);
	EXPECT_NEAR(bond->hybrids[0].weight_percent + bond->hybrids[1].weight_percent, 100.0, 0.05);
	for (const auto& h : bond->hybrids) EXPECT_NEAR(h.s + h.p + h.d + h.f, 100.0, 0.05);
	EXPECT_LT(bond->energy, 0.0);

	ASSERT_FALSE(r.e2.empty());
	EXPECT_NE(r.e2.front().donor.find("LP"), std::string::npos);
	for (const auto& e : r.e2) {
		EXPECT_GT(e.energy_kcal, 0.0);
		EXPECT_NE(e.donor_index, e.acceptor_index);
	}
}

TEST(NboRun, ComparisonPassesAgainstItselfAndCatchesAShiftedCharge)
{
	const auto reference_nbo = repo_root() / "tests" / "epoxide_gbw" / "NBO" / "reference.nbo";
	if (!std::filesystem::exists(reference_nbo)) {
		GTEST_SKIP() << "NBO reference fixture is not available";
	}

	const NboResults reference = parse_nbo_output(reference_nbo);
	const NboComparison same = compare_nbo_results(reference, reference);
	EXPECT_TRUE(same.ok) << same.report();
	for (const auto& q : same.quantities) EXPECT_EQ(q.missing, 0) << q.quantity;

	NboResults shifted = reference;
	shifted.npa.front().charge += 0.01;
	shifted.e2.front().energy_kcal += 1.0;
	const NboComparison differs = compare_nbo_results(reference, shifted);
	EXPECT_FALSE(differs.ok);

	NboResults truncated = reference;
	truncated.orbitals.pop_back();
	const NboComparison missing = compare_nbo_results(reference, truncated);
	EXPECT_FALSE(missing.ok);
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

/*
 * The NRT capture, against the acetylene output of the reference set (NRT E2PERT NRTLST=0.1
 * NRTDTL). Everything asserted here is a number NBO 7.0.9 printed, so a parser that starts
 * dropping rows - the zero-weight tail, the diagonal of the bond-order matrix, a resonance
 * structure whose Added(Removed) column wrapped onto a second line - fails here rather than
 * silently shipping a short reference.
 */
TEST(NboRun, CapturesTheFullNrtSectionOfTheAcetyleneReference)
{
	const auto nbo = repo_root() / "tests" / "nbo_reference" / "acetylene_nrtdtl.nbo";
	if (!std::filesystem::exists(nbo)) GTEST_SKIP() << "NRT reference fixture is not available";

	const NboResults r = parse_nbo_output(nbo);
	//The keylist as NBO echoed it back. Re-deriving a stored reference after a parser change goes
	//through -nbo_parse, and then this is the only record of what the run was asked for:
	//r.keywords stays empty because nothing passed a keylist in.
	EXPECT_EQ(r.keywords_reported, "NRT NRTLST NRTDTL E2PERT");
	ASSERT_TRUE(r.nrt.present);
	EXPECT_EQ(r.nrt.structures_used, 7);
	EXPECT_EQ(r.nrt.structures_found, 15);
	EXPECT_NEAR(r.nrt.d_w, 0.01830453, 1.0e-8);
	EXPECT_NEAR(r.nrt.d_0, 0.01884235, 1.0e-8);
	EXPECT_EQ(r.nrt.max_search_cycles, 3);
	EXPECT_EQ(r.nrt.initial_topo, 1);
	EXPECT_NE(r.nrt.symmetry.find("symmetry operator"), std::string::npos);
	EXPECT_GT(r.nbo_cpu_seconds, 0.0);

	//The search table: two cycles, the second one generating nothing new.
	ASSERT_EQ(r.nrt.cycles.size(), 2u);
	EXPECT_EQ(r.nrt.cycles[0].structures_found, 1);
	EXPECT_EQ(r.nrt.cycles[1].structures_used, 7);
	EXPECT_EQ(r.nrt.cycles[1].structures_found, 15);
	EXPECT_EQ(r.nrt.cycles[1].e2, 0);

	//All 15 candidates, not only the 7 with weight: the ratio is what a screening scheme has
	//to beat, so the zero-weight tail has to survive the parse.
	ASSERT_EQ(r.nrt.weights.size(), 15u);
	EXPECT_NEAR(r.nrt.weights[0].weight_percent, 95.70, 1.0e-6);
	EXPECT_NEAR(r.nrt.weights[0].weight_fraction, 0.95696, 1.0e-6);
	EXPECT_NE(r.nrt.weights[1].changes.find("C 1- C 2"), std::string::npos);
	int zero_weight = 0;
	for (const auto& w : r.nrt.weights) if (w.weight_fraction == 0.0) zero_weight++;
	EXPECT_EQ(zero_weight, 8);
	ASSERT_EQ(r.nrt.candidates.size(), 15u);
	EXPECT_NEAR(r.nrt.candidates[0].rho_nl, 0.02527, 1.0e-6);
	ASSERT_EQ(r.nrt.candidates[0].topo.size(), 4u);
	EXPECT_EQ(r.nrt.candidates[0].topo[0][1], 3);   //the C-C triple bond of the leading structure

	//The QP path, so a candidate implementation can be compared step by step and not only at
	//the converged answer.
	ASSERT_GE(r.nrt.qp_iterations.size(), 8u);
	EXPECT_NEAR(r.nrt.qp_iterations.back().d_w, 0.01830453, 1.0e-8);
	//Two "Perform ARROWS on structures of weight > X%" lines (the parent threshold as NBO applied
	//it, once per cycle) and the one line naming the parent structure and the E2 depth used.
	ASSERT_EQ(r.nrt.arrows.size(), 3u);
	EXPECT_NE(r.nrt.arrows[1].find("generates 6 new structures from structure 1"), std::string::npos);
	EXPECT_NE(r.nrt.arrows[1].find("E(2)=1.0 kcal/mol"), std::string::npos);

	//The bond-order matrix as printed: upper triangle plus diagonal of a 4-atom system.
	ASSERT_EQ(r.nrt.bond_orders.size(), 10u);
	const NboBondOrder* cc = nullptr;
	const NboBondOrder* diag = nullptr;
	for (const auto& b : r.nrt.bond_orders) {
		if (b.atom1 == 1 && b.atom2 == 2) cc = &b;
		if (b.atom1 == 1 && b.atom2 == 1) diag = &b;
	}
	ASSERT_NE(cc, nullptr);
	ASSERT_NE(diag, nullptr);
	EXPECT_NEAR(cc->total, 2.9938, 1.0e-6);
	EXPECT_NEAR(cc->covalent, 2.9938, 1.0e-6);
	EXPECT_NEAR(cc->ionic, 0.0, 1.0e-6);
	EXPECT_TRUE(diag->diagonal);
	EXPECT_NEAR(diag->total, 0.0198, 1.0e-6);

	//Valencies, the atom-by-atom sum of those bond orders.
	ASSERT_EQ(r.nrt.valencies.size(), 4u);
	EXPECT_EQ(r.nrt.valencies[0].element, "C");
	EXPECT_NEAR(r.nrt.valencies[0].valency, 3.9631, 1.0e-6);
	EXPECT_NEAR(r.nrt.valencies[0].covalency, 3.7369, 1.0e-6);
	EXPECT_NEAR(r.nrt.valencies[0].electron_count, 7.9657, 1.0e-6);

	ASSERT_EQ(r.nrt.leading_topo.size(), 1u);
	EXPECT_EQ(r.nrt.leading_topo[0].matrix[1][0], 3);
	EXPECT_NE(r.nrt.nrtstr_keylist.find("STR"), std::string::npos);
}
