#include "pch.h"
#include "CppUnitTest.h"

#include "../core/NoSpherA2.h"

#include <Windows.h>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <optional>
#include <string>
#include <vector>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
namespace {

static const std::regex kNumberPattern(R"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)");

static std::vector<std::string> read_lines_stripped(const std::filesystem::path& path)
{
    std::ifstream f(path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    return lines;
}

static std::vector<double> extract_floats(const std::string& line)
{
    std::vector<double> result;
    for (auto it = std::sregex_iterator(line.begin(), line.end(), kNumberPattern);
         it != std::sregex_iterator();
         ++it) {
        result.push_back(std::stod((*it).str()));
    }
    return result;
}

static std::string normalize_whitespace(const std::string& s)
{
    std::istringstream iss(s);
    std::ostringstream oss;
    std::string token;
    bool first = true;
    while (iss >> token) {
        if (!first) {
            oss << ' ';
        }
        oss << token;
        first = false;
    }
    return oss.str();
}

static std::string make_skeleton(const std::string& line)
{
    return normalize_whitespace(std::regex_replace(line, kNumberPattern, "<num>"));
}

static bool numeric_line_matches(const std::string& expected,
    const std::string& actual,
    double rtol)
{
    const auto exp_nums = extract_floats(expected);
    const auto act_nums = extract_floats(actual);
    if (exp_nums.empty() || act_nums.empty() || exp_nums.size() != act_nums.size()) {
        return false;
    }
    if (make_skeleton(expected) != make_skeleton(actual)) {
        return false;
    }
    for (size_t i = 0; i < exp_nums.size(); ++i) {
        const double diff = std::abs(exp_nums[i] - act_nums[i]);
        const double allowed = std::abs(act_nums[i]) * rtol;
        if (diff > allowed) {
            return false;
        }
    }
    return true;
}

static std::string compare_files_numeric(const std::filesystem::path& good_path,
    const std::filesystem::path& actual_path,
    double rtol)
{
    const auto expected = read_lines_stripped(good_path);
    const auto actual = read_lines_stripped(actual_path);

    std::ostringstream failures;
    bool bad = false;
    const size_t n = (expected.size() > actual.size()) ? expected.size() : actual.size();
    for (size_t i = 0; i < n; ++i) {
        if (i >= expected.size()) {
            failures << "  Extra line " << (i + 1) << " in actual: " << actual[i] << "\n";
            bad = true;
        }
        else if (i >= actual.size()) {
            failures << "  Missing line " << (i + 1) << ": " << expected[i] << "\n";
            bad = true;
        }
        else if (expected[i] == actual[i]) {
            continue;
        }
        else if (numeric_line_matches(expected[i], actual[i], rtol)) {
            continue;
        }
        else {
            failures << "  Line " << (i + 1) << ":\n"
                     << "    expected: " << expected[i] << "\n"
                     << "    actual:   " << actual[i] << "\n";
            bad = true;
        }
    }
    return bad ? failures.str() : std::string{};
}

static std::filesystem::path get_repo_root_from_module()
{
    HMODULE hMod = nullptr;
    GetModuleHandleExA(
        GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
        reinterpret_cast<LPCSTR>(&get_repo_root_from_module), &hMod);
    char buf[MAX_PATH] = {};
    GetModuleFileNameA(hMod, buf, MAX_PATH);
    auto dll_dir = std::filesystem::path(buf).parent_path();
    // Windows/x64/<Config>/ → repo root
    return dll_dir.parent_path().parent_path().parent_path();
}

struct TomlTestDef {
    std::string directory;
    std::string goodFile;
    std::string actualFile;
    bool full = false;
    // Each entry is a flat list of CLI tokens to append after "-key".
    // For scalar args this is a single value; for array args this is N values.
    std::vector<std::pair<std::string, std::vector<std::string>>> args;
};

static std::string trim(std::string s)
{
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    while (!s.empty() && is_space(static_cast<unsigned char>(s.front()))) {
        s.erase(s.begin());
    }
    while (!s.empty() && is_space(static_cast<unsigned char>(s.back()))) {
        s.pop_back();
    }
    return s;
}

static std::string unquote_if_needed(std::string v)
{
    v = trim(v);
    if (v.size() >= 2 && v.front() == '"' && v.back() == '"') {
        return v.substr(1, v.size() - 2);
    }
    return v;
}

static std::string strip_comment(const std::string& line)
{
    // TOML comments start with #, but ignore if inside a quoted string.
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '"') {
            in_quotes = !in_quotes;
        }
        if (!in_quotes && c == '#') {
            return line.substr(0, i);
        }
    }
    return line;
}

static std::optional<std::string> parse_section_name(const std::string& line)
{
    // [section]
    auto t = trim(line);
    if (t.size() < 3 || t.front() != '[' || t.back() != ']') {
        return std::nullopt;
    }
    auto inner = t.substr(1, t.size() - 2);
    inner = trim(inner);
    if (inner.empty()) {
        return std::nullopt;
    }
    return inner;
}

static std::optional<std::pair<std::string, std::string>> parse_kv(const std::string& line)
{
    // key = value
    auto pos = line.find('=');
    if (pos == std::string::npos) {
        return std::nullopt;
    }

    std::string key = trim(line.substr(0, pos));
    std::string value = trim(line.substr(pos + 1));
    if (key.empty() || value.empty()) {
        return std::nullopt;
    }
    return std::make_pair(key, value);
}

static bool is_array_start(const std::string& v)
{
    auto t = trim(v);
    return !t.empty() && t.front() == '[';
}

static bool is_array_end(const std::string& v)
{
    auto t = trim(v);
    return !t.empty() && t.back() == ']';
}

static std::vector<std::string> split_array_items(const std::string& inner)
{
    // Minimal TOML array splitting: handles quoted strings and numbers/bools.
    // No nested arrays.
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (size_t i = 0; i < inner.size(); ++i) {
        char c = inner[i];
        if (c == '"') {
            in_quotes = !in_quotes;
            cur.push_back(c);
            continue;
        }
        if (!in_quotes && c == ',') {
            auto tok = trim(cur);
            if (!tok.empty()) {
                out.push_back(unquote_if_needed(tok));
            }
            cur.clear();
            continue;
        }
        cur.push_back(c);
    }
    auto tok = trim(cur);
    if (!tok.empty()) {
        out.push_back(unquote_if_needed(tok));
    }
    return out;
}

static std::optional<std::vector<std::string>> parse_array_value(std::string value,
    std::istream& in)
{
    // Supports:
    //   key = [1, 2, 3]
    //   key = [
    //     "a",
    //     1,
    //   ]
    value = strip_comment(value);
    value = trim(value);
    if (!is_array_start(value)) {
        return std::nullopt;
    }

    std::string buf = value;
    while (!is_array_end(buf)) {
        std::string more;
        if (!std::getline(in, more)) {
            break;
        }
        more = strip_comment(more);
        more = trim(more);
        if (more.empty()) {
            continue;
        }
        buf += " " + more;
    }

    auto t = trim(buf);
    if (t.size() < 2 || t.front() != '[' || t.back() != ']') {
        return std::nullopt;
    }
    auto inner = t.substr(1, t.size() - 2);
    inner = trim(inner);
    if (inner.empty()) {
        return std::vector<std::string>{};
    }
    return split_array_items(inner);
}

static std::optional<TomlTestDef> load_test_def_from_toml(const std::filesystem::path& toml_path,
    const std::string& test_name)
{
    std::ifstream in(toml_path);
    if (!in) {
        return std::nullopt;
    }

    const std::string section_main = test_name;
    const std::string section_args = test_name + ".args";

    enum class Mode { None, InMain, InArgs };
    Mode mode = Mode::None;

    TomlTestDef def;
    // Defaults from [defaults] in tests/tests.toml are currently two flags.
    def.args.emplace_back("all_charges", std::vector<std::string>{"true"});
    def.args.emplace_back("no_date", std::vector<std::string>{"true"});
    std::string line;
    while (std::getline(in, line)) {
        line = strip_comment(line);
        line = trim(line);
        if (line.empty()) {
            continue;
        }

        if (auto sec = parse_section_name(line); sec.has_value()) {
            if (sec.value() == section_main) {
                mode = Mode::InMain;
            } else if (sec.value() == section_args) {
                mode = Mode::InArgs;
            } else {
                mode = Mode::None;
            }
            continue;
        }

        auto kv = parse_kv(line);
        if (!kv.has_value()) {
            continue;
        }

        const std::string key = kv->first;
        std::string value = kv->second;

        if (mode == Mode::InMain) {
            if (key == "directory") {
                def.directory = unquote_if_needed(value);
            } else if (key == "good") {
                def.goodFile = unquote_if_needed(value);
            } else if (key == "actual") {
                def.actualFile = unquote_if_needed(value);
            } else if (key == "full") {
                auto v = unquote_if_needed(value);
                def.full = (v == "true" || v == "True" || v == "1");
            }
        } else if (mode == Mode::InArgs) {
            // Scalars or arrays.
            if (auto arr = parse_array_value(value, in); arr.has_value()) {
                std::vector<std::string> vals;
                vals.reserve(arr->size());
                for (const auto& v : *arr) {
                    vals.push_back(unquote_if_needed(v));
                }
                def.args.emplace_back(key, vals);
            } else {
                def.args.emplace_back(key, std::vector<std::string>{unquote_if_needed(value)});
            }
        }
    }

    if (def.directory.empty()) {
        return std::nullopt;
    }

    return def;
}

static int run_inprocess_test(const std::filesystem::path& repo_root,
    const std::string& test_name)
{
    const auto toml_path = repo_root / "tests" / "tests.toml";
    auto def_opt = load_test_def_from_toml(toml_path, test_name);
    if (!def_opt.has_value()) {
        return -1;
    }

    const auto test_src_dir = repo_root / "tests" / def_opt->directory;
    if (!std::filesystem::exists(test_src_dir)) {
        return -2;
    }

    // Inject OCC_DATA_PATH for in-process runs.
    const std::string occDataPath = (repo_root / "Lib" / "occ" / "share" / "occ").string();
    SetEnvironmentVariableA("OCC_DATA_PATH", occDataPath.c_str());

    // Build argv: NoSpherA2 -key value...
    std::vector<std::string> argv_storage;
    argv_storage.emplace_back("NoSpherA2");
    for (const auto& kv : def_opt->args) {
        const std::string flag = "-" + kv.first;
        argv_storage.emplace_back(flag);
        for (const auto& val : kv.second) {
            if (val == "true" || val == "True" || val == "") {
                // flag-only semantics: keep just -key
                continue;
            }
            if (val == "false" || val == "False") {
                // omit false values
                continue;
            }
            argv_storage.emplace_back(val);
        }
    }

    std::vector<char*> argv;
    argv.reserve(argv_storage.size() + 1);
    for (auto& s : argv_storage) {
        argv.push_back(s.data());
    }
    argv.push_back(nullptr);

    const auto old_cwd = std::filesystem::current_path();
    // Move to the directory of the test:
    std::filesystem::current_path(test_src_dir);
    int rc = run_app(static_cast<int>(argv_storage.size()), argv.data());
    std::filesystem::current_path(old_cwd);
    if (rc != 0) {
        return rc;
    }

    std::filesystem::path out_file = test_src_dir / "NoSpherA2.log";
    if (!def_opt->actualFile.empty()) {
        out_file = test_src_dir / def_opt->actualFile;
    }

    if (!std::filesystem::exists(out_file)) {
        return -2;
    }
    Logger::WriteMessage(std::wstring(L"Log file: " + out_file.wstring()).c_str());

    const std::string good_name = def_opt->goodFile.empty() ? (test_name + ".good") : def_opt->goodFile;
    const auto good_path = test_src_dir / good_name;

    if (!std::filesystem::exists(good_path)) {
        return -2;
    }

    const auto diff = compare_files_numeric(good_path, out_file, 0.01);
    if (!diff.empty()) {
        Logger::WriteMessage(std::wstring(diff.begin(), diff.end()).c_str());
        return -3;
    }

    return 0;
}

static std::filesystem::path guess_repo_root()
{
    // Prefer a path derived from the test module location (stable under vstest).
    auto root = get_repo_root_from_module();
    if (std::filesystem::exists(root / "tests" / "tests.toml")) {
        return root;
    }

    // Fallback: walk up from CWD.
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

} // namespace

namespace NoSpherA2Integration {

TEST_CLASS(TomlIntegrationTests)
{
public:
    static std::filesystem::path repo_root;

    TEST_CLASS_INITIALIZE(ClassInitialize)
    {
        repo_root = guess_repo_root();
    }

    TEST_METHOD(AlanineOcc)
    {
        const int rc = run_inprocess_test(repo_root, "alanine_occ");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(AlanineIntegratedOcc)
    {
        const int rc = run_inprocess_test(repo_root, "alanine_integrated_occ");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(DisorderTHPP)
    {
        const int rc = run_inprocess_test(repo_root, "disorder_THPP");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(Fractal)
    {
        const int rc = run_inprocess_test(repo_root, "fractal");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(GrownWater)
    {
        const int rc = run_inprocess_test(repo_root, "grown_water");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(HybridMode)
    {
        const int rc = run_inprocess_test(repo_root, "Hybrid_mode");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(MalbacSfEcp)
    {
        const int rc = run_inprocess_test(repo_root, "malbac_SF_ECP");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(Properties)
    {
        const int rc = run_inprocess_test(repo_root, "properties");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(RiFit)
    {
        const int rc = run_inprocess_test(repo_root, "ri_fit");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(RubredoxinCmtc)
    {
        const int rc = run_inprocess_test(repo_root, "rubredoxin_cmtc");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(SALTED)
    {
        const int rc = run_inprocess_test(repo_root, "SALTED");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(SucroseIAM)
    {
        const int rc = run_inprocess_test(repo_root, "sucrose_IAM");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(SucrosePtb)
    {
        const int rc = run_inprocess_test(repo_root, "sucrose_ptb");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(SucroseSF)
    {
        const int rc = run_inprocess_test(repo_root, "sucrose_SF");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(SucroseTwin)
    {
        const int rc = run_inprocess_test(repo_root, "sucrose_twin");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(WfnReading)
    {
        const int rc = run_inprocess_test(repo_root, "wfn_reading");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(FchkConversion)
    {
        const int rc = run_inprocess_test(repo_root, "fchk_conversion");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(TFVC)
    {
        const int rc = run_inprocess_test(repo_root, "TFVC");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(TFVCEcp)
    {
        const int rc = run_inprocess_test(repo_root, "TFVC_ECP");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(RGBI_Groups_NH3BH3)
    {
        const int rc = run_inprocess_test(repo_root, "RGBI_Groups_NH3BH3");
        Assert::AreEqual(0, rc);
    }

    TEST_METHOD(RGBI_Groups_NH3Li)
    {
        const int rc = run_inprocess_test(repo_root, "RGBI_Groups_NH3Li");
        Assert::AreEqual(0, rc);
    }
};

std::filesystem::path TomlIntegrationTests::repo_root;

}
