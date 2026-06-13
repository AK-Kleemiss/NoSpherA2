#pragma once
#include <Windows.h>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <regex>
#include <optional>
#include <algorithm>
#include <cstdlib>
#include "CppUnitTest.h"

// In-process entry point exported by NoSpherA2_DLL.dll.
// Linked via the import library; the DLL must be in the same output directory.
extern "C" __declspec(dllimport) int __cdecl NoSpherA2_run(
    const char* workDir,
    int         argc,
    char**      argv);

namespace NosTestFramework
{
    using namespace Microsoft::VisualStudio::CppUnitTestFramework;

    // ---------------------------------------------------------------------------
    // Path helpers
    // ---------------------------------------------------------------------------

    inline std::filesystem::path GetDllDirectory()
    {
        HMODULE hMod = nullptr;
        GetModuleHandleExA(
            GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            reinterpret_cast<LPCSTR>(&GetDllDirectory), &hMod);
        char buf[MAX_PATH] = {};
        GetModuleFileNameA(hMod, buf, MAX_PATH);
        return std::filesystem::path(buf).parent_path();
    }

    // DLL lives at  Windows/x64/<Config>/  →  three levels up = repo root
    inline std::filesystem::path GetRepoRoot()
    {
        return GetDllDirectory().parent_path().parent_path().parent_path();
    }

    // ---------------------------------------------------------------------------
    // Numeric comparison helpers  (mirrors run_test.py logic)
    // ---------------------------------------------------------------------------

    static const std::regex kNumberPattern(R"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)");

    inline std::vector<double> ExtractFloats(const std::string &line)
    {
        std::vector<double> result;
        auto begin = std::sregex_iterator(line.begin(), line.end(), kNumberPattern);
        for (auto it = begin; it != std::sregex_iterator(); ++it)
            result.push_back(std::stod((*it).str()));
        return result;
    }

    inline std::string NormalizeWhitespace(const std::string &s)
    {
        std::istringstream iss(s);
        std::ostringstream oss;
        std::string token;
        bool first = true;
        while (iss >> token) {
            if (!first) oss << ' ';
            oss << token;
            first = false;
        }
        return oss.str();
    }

    inline std::string MakeSkeleton(const std::string &line)
    {
        return NormalizeWhitespace(std::regex_replace(line, kNumberPattern, "<num>"));
    }

    // Returns true if two lines agree numerically within rtol (default 1 %)
    inline bool NumericLineMatches(const std::string &expected, const std::string &actual,
        double rtol = 0.01)
    {
        auto exp_nums = ExtractFloats(expected);
        auto act_nums = ExtractFloats(actual);
        if (exp_nums.empty() || act_nums.empty() || exp_nums.size() != act_nums.size())
            return false;
        if (MakeSkeleton(expected) != MakeSkeleton(actual))
            return false;
        for (size_t i = 0; i < exp_nums.size(); ++i) {
            double diff = std::abs(exp_nums[i] - act_nums[i]);
            double allowed = std::abs(act_nums[i]) * rtol;
            if (diff > allowed) return false;
        }
        return true;
    }

    // Returns an empty string on success; a human-readable diff on failure.
    inline std::string CompareFiles(const std::filesystem::path &goodPath,
        const std::filesystem::path &actualPath,
        double rtol = 0.01)
    {
        auto readLines = [](const std::filesystem::path &p) -> std::vector<std::string> {
            std::ifstream f(p);
            std::vector<std::string> lines;
            std::string line;
            while (std::getline(f, line)) {
                // strip CR, skip blank lines (mirrors Python strip + filter)
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (!line.empty()) lines.push_back(line);
            }
            return lines;
            };

        auto expected = readLines(goodPath);
        auto actual = readLines(actualPath);

        std::ostringstream failures;
        bool bad = false;

        size_t n = (((expected.size()) > (actual.size())) ? (expected.size()) : (actual.size()));
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
                continue; // exact match
            }
            else if (NumericLineMatches(expected[i], actual[i], rtol)) {
                continue; // within numeric tolerance (tolerated difference)
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

    // ---------------------------------------------------------------------------
    // Test definition
    // ---------------------------------------------------------------------------

    struct TestDef
    {
        std::string name;           // test name (for temp dir naming)
        std::string directory;      // subdirectory under tests/
        std::string goodFile;       // e.g. "alanine_occ.good"
        // actual output: if empty, use goodFile with .good→.log replacement
        std::string actualFile;
        std::vector<std::string> args; // flat list of CLI tokens
        bool full = false;             // skip when RUN_FULL_TEST not set
    };

    // ---------------------------------------------------------------------------
    // Run a test
    // ---------------------------------------------------------------------------

    inline void RunTest(const TestDef &test)
    {
        // Skip full tests unless opted in
        if (test.full) {
            const char *env = std::getenv("RUN_FULL_TEST");
            if (!env || std::string(env).empty()) {
                Logger::WriteMessage(L"Skipping: RUN_FULL_TEST not set");
                return;
            }
        }

        auto repoRoot = GetRepoRoot();
        auto originalCwd = std::filesystem::current_path();
        auto testSrcDir = repoRoot / "tests" / test.directory;
        if (!std::filesystem::exists(testSrcDir)) {
            auto msg = std::wstring(L"Test directory not found: ") +
                testSrcDir.wstring();
            Assert::Fail(msg.c_str());
        }

        // Create isolated temp work directory. Include process/run identity so
        // a stale Visual Studio test host cannot lock the next run's folder.
        auto tempRoot = std::filesystem::temp_directory_path() / "NosTests";
        auto tempBase = tempRoot /
            (test.name + "_" + std::to_string(GetCurrentProcessId()) + "_" +
             std::to_string(GetTickCount64()));
        auto cleanupTemp = [&]() {
            std::error_code ignored;
            std::filesystem::current_path(originalCwd, ignored);
            SetCurrentDirectoryW(originalCwd.wstring().c_str());

            std::error_code cleanupEc;
            std::filesystem::remove_all(tempBase, cleanupEc);
            if (cleanupEc) {
                auto cleanupMessage = cleanupEc.message();
                auto msg = std::wstring(L"Warning: could not clean temp test directory: ") +
                    tempBase.wstring() + L" (" +
                    std::wstring(cleanupMessage.begin(), cleanupMessage.end()) +
                    L")";
                Logger::WriteMessage(msg.c_str());
            }
        };
        {
            std::error_code ignored;
            std::filesystem::remove_all(tempRoot / test.name, ignored);
        }
        std::filesystem::create_directories(tempBase);
        std::filesystem::copy(testSrcDir, tempBase,
            std::filesystem::copy_options::recursive |
            std::filesystem::copy_options::overwrite_existing);

        // OCC_DATA_PATH — injected for in-process runs.
        // The prebuilt OCC data lives at Lib/occ/share/occ (not occ/share,
        // which is the gh-pages documentation checkout).
        std::string occDataPath = (repoRoot / "Lib" / "occ" / "share" / "occ").string();
        SetEnvironmentVariableA("OCC_DATA_PATH", occDataPath.c_str());

        // ---- in-process path -----------------------------------------------
        // All tests run in-process through NoSpherA2_DLL.dll so Visual Studio
        // can debug Src/ code directly in the test host.
        std::vector<std::string> argStrings;
        // argv[0] = DLL path (used by ensure_occ_data_path fallback logic)
        argStrings.push_back((GetDllDirectory() / "NoSpherA2_DLL.dll").string());
        for (const auto &a : test.args)
            argStrings.push_back(a);
        std::vector<char*> argPtrs;
        for (auto &s : argStrings)
            argPtrs.push_back(const_cast<char*>(s.c_str()));

        std::string workDirStr = tempBase.string();
        int exitCode = NoSpherA2_run(
            workDirStr.c_str(),
            static_cast<int>(argPtrs.size()),
            argPtrs.data());

        if (exitCode != 0) {
            cleanupTemp();
            Assert::Fail(
                (std::wstring(L"NoSpherA2 exited with code ") +
                 std::to_wstring(exitCode) + L" for test " +
                 std::wstring(test.name.begin(), test.name.end())).c_str());
        }

        // Determine actual output file path
        std::filesystem::path actualPath;
        if (!test.actualFile.empty()) {
            // Custom actual file (e.g. fractal test)
            actualPath = tempBase / test.actualFile;
        }
        else {
            // Default: exe writes NoSpherA2.log; compare under good-file name with .log ext
            std::string logName = test.goodFile;
            actualPath = tempBase / logName;
            // If the renamed file doesn't exist, fall back to NoSpherA2.log
            if (!std::filesystem::exists(actualPath)) {
                auto fallback = tempBase / "NoSpherA2.log";
                if (std::filesystem::exists(fallback))
                    actualPath = fallback;
            }
        }

        if (!std::filesystem::exists(actualPath)) {
            cleanupTemp();
            auto msg = std::wstring(L"Output file not found: ") + actualPath.wstring();
            Assert::Fail(msg.c_str());
        }

        auto goodPath = testSrcDir / test.goodFile;
        if (!std::filesystem::exists(goodPath)) {
            cleanupTemp();
            auto msg = std::wstring(L"Good file not found: ") + goodPath.wstring();
            Assert::Fail(msg.c_str());
        }

        auto diffResult = CompareFiles(goodPath, actualPath);
        cleanupTemp();

        if (!diffResult.empty()) {
            auto wmsg = std::wstring(diffResult.begin(), diffResult.end());
            Assert::Fail((std::wstring(L"Output mismatch:\n") + wmsg).c_str());
        }
    }

} // namespace NosTestFramework
