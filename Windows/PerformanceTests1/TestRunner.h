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

        // Find NoSpherA2.exe (same output dir as this DLL)
        auto dllDir = GetDllDirectory();
        auto exePath = dllDir / "NoSpherA2.exe";
        if (!std::filesystem::exists(exePath)) {
            // Also try NOS_EXE env var
            const char *envExe = std::getenv("NOS_EXE");
            if (envExe && std::filesystem::exists(envExe))
                exePath = envExe;
            else
                Assert::Fail(L"NoSpherA2.exe not found. Build it first or set NOS_EXE.");
        }

        auto repoRoot = GetRepoRoot();
        auto testSrcDir = repoRoot / "tests" / test.directory;
        if (!std::filesystem::exists(testSrcDir)) {
            auto msg = std::wstring(L"Test directory not found: ") +
                testSrcDir.wstring();
            Assert::Fail(msg.c_str());
        }

        // Create isolated temp work directory
        auto tempBase = std::filesystem::temp_directory_path() / "NosTests" / test.name;
        std::error_code ec;
        std::filesystem::remove_all(tempBase, ec);
        std::filesystem::create_directories(tempBase);
        std::filesystem::copy(testSrcDir, tempBase,
            std::filesystem::copy_options::recursive |
            std::filesystem::copy_options::overwrite_existing);

        // Build command line string
        std::string cmdLine = "\"" + exePath.string() + "\"";
        for (const auto &a : test.args)
            cmdLine += " " + a;

        // Set OCC_DATA_PATH in the child environment
        std::string occDataPath = (repoRoot / "occ" / "share").string();
        std::wstring envBlock;
        // Collect current environment
        LPWCH envStrings = GetEnvironmentStrings();
        for (LPWCH p = envStrings; *p; ) {
            envBlock += p;
            envBlock += wchar_t('\0');
            p += wcslen(p) + 1;
        }
        FreeEnvironmentStrings(envStrings);
        // Override/add OCC_DATA_PATH
        envBlock += L"OCC_DATA_PATH=" + std::wstring(occDataPath.begin(), occDataPath.end()) + L'\0';
        // When VS debugger is attached to this test host, tell the child to wait
        // so you can attach to it too (Debug → Attach to Process, use the PID
        // printed in the Test Explorer output window).
        bool debuggerPresent = IsDebuggerPresent() != 0;
        if (debuggerPresent) {
            envBlock += L"DEBUG_WAIT=1";
            envBlock += L'\0';
        }
        envBlock += L'\0'; // double-null terminator

        // Launch process
        STARTUPINFOA si = {};
        si.cb = sizeof(si);
        si.dwFlags = STARTF_USESTDHANDLES;
        si.hStdInput = INVALID_HANDLE_VALUE;
        si.hStdOutput = INVALID_HANDLE_VALUE;
        si.hStdError = INVALID_HANDLE_VALUE;

        PROCESS_INFORMATION pi = {};
        std::string workDir = tempBase.string();

        BOOL ok = CreateProcessA(
            nullptr,
            const_cast<LPSTR>(cmdLine.c_str()),
            nullptr, nullptr, FALSE,
            CREATE_NO_WINDOW,
            reinterpret_cast<LPVOID>(const_cast<wchar_t *>(envBlock.data())),
            workDir.c_str(),
            &si, &pi);

        if (!ok) {
            auto err = std::to_wstring(GetLastError());
            Assert::Fail((std::wstring(L"CreateProcess failed: ") + err).c_str());
        }

        if (debuggerPresent) {
            // Child is spinning in wait_for_debugger(). Attach now:
            //   Debug → Attach to Process (Ctrl+Alt+P) → PID shown below
            auto pidMsg = std::wstring(L"[DEBUG] NoSpherA2.exe PID: ") +
                std::to_wstring(pi.dwProcessId) +
                L" — attach with Debug → Attach to Process (Ctrl+Alt+P)";
            Logger::WriteMessage(pidMsg.c_str());
        }

        WaitForSingleObject(pi.hProcess, INFINITE);
        DWORD exitCode = 0;
        GetExitCodeProcess(pi.hProcess, &exitCode);
        CloseHandle(pi.hProcess);
        CloseHandle(pi.hThread);

        if (exitCode != 0) {
            std::filesystem::remove_all(tempBase);
            auto msg = std::wstring(L"NoSpherA2 exited with code ") +
                std::to_wstring(exitCode) + L" for test " +
                std::wstring(test.name.begin(), test.name.end());
            Assert::Fail(msg.c_str());
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
            auto pos = logName.rfind(".good");
            if (pos != std::string::npos)
                logName.replace(pos, 5, ".log");
            actualPath = tempBase / logName;
            // If the renamed file doesn't exist, fall back to NoSpherA2.log
            if (!std::filesystem::exists(actualPath)) {
                auto fallback = tempBase / "NoSpherA2.log";
                if (std::filesystem::exists(fallback))
                    actualPath = fallback;
            }
        }

        if (!std::filesystem::exists(actualPath)) {
            std::filesystem::remove_all(tempBase);
            auto msg = std::wstring(L"Output file not found: ") + actualPath.wstring();
            Assert::Fail(msg.c_str());
        }

        auto goodPath = testSrcDir / test.goodFile;
        if (!std::filesystem::exists(goodPath)) {
            std::filesystem::remove_all(tempBase);
            auto msg = std::wstring(L"Good file not found: ") + goodPath.wstring();
            Assert::Fail(msg.c_str());
        }

        auto diffResult = CompareFiles(goodPath, actualPath);
        std::filesystem::remove_all(tempBase);

        if (!diffResult.empty()) {
            auto wmsg = std::wstring(diffResult.begin(), diffResult.end());
            Assert::Fail((std::wstring(L"Output mismatch:\n") + wmsg).c_str());
        }
    }

} // namespace NosTestFramework
