#include "pch.h"

#include "core/wfn_class.h"

namespace {

std::filesystem::path repo_root()
{
    if (const char* env = std::getenv("NOS_REPO_ROOT")) {
        return std::filesystem::path(env);
    }

    auto p = std::filesystem::current_path();
    for (int i = 0; i < 8; ++i) {
        if (std::filesystem::exists(p / "tests" / "tests.toml")) {
            return p;
        }
        if (!p.has_parent_path()) {
            break;
        }
        p = p.parent_path();
    }
    return std::filesystem::current_path();
}

std::string read_file(const std::filesystem::path& path)
{
    std::ifstream in(path);
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

std::vector<double> extract_section_numbers(const std::string& text, const std::string& section)
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
    std::vector<double> values;
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
    const auto base = std::filesystem::temp_directory_path() / "nos_nbo47_test";
    std::filesystem::remove_all(base);
    std::filesystem::create_directories(base);
    return base;
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

std::vector<double> parse_natural_charges(const std::filesystem::path& nbo_path)
{
    std::ifstream in(nbo_path);
    std::string line;
    bool in_summary = false;
    std::vector<double> charges;
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
    ASSERT_TRUE(std::filesystem::exists(reference_nbo));

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

    EXPECT_EQ(parse_key_int(text, "NATOMS").value_or(-1), 5);
    EXPECT_EQ(parse_key_int(text, "NBAS").value_or(-1), 63);
    EXPECT_EQ(parse_key_int(text, "NSHELL").value_or(-1), 31);
    EXPECT_EQ(parse_key_int(text, "NEXP").value_or(-1), 52);

    const auto overlap = extract_section_numbers(text, "$OVERLAP");
    const auto density = extract_section_numbers(text, "$DENSITY");
    const auto fock = extract_section_numbers(text, "$FOCK");
    const auto lcaomo = extract_section_numbers(text, "$LCAOMO");
    const int nbasis = 63;
    EXPECT_EQ(overlap.size(), static_cast<size_t>(nbasis * (nbasis + 1) / 2));
    EXPECT_EQ(density.size(), static_cast<size_t>(nbasis * (nbasis + 1) / 2));
    EXPECT_EQ(fock.size(), static_cast<size_t>(nbasis * (nbasis + 1) / 2));
    EXPECT_EQ(lcaomo.size(), static_cast<size_t>(nbasis * nbasis));
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
