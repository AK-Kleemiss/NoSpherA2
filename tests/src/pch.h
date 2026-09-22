// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.

#ifndef PCH_H
#define PCH_H

#include <cmath>
#include <limits>
#include <string>
#include <chrono>
#include <thread>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <optional>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cctype>

#include "gtest/gtest.h"

//The repository root: NOS_REPO_ROOT, else the nearest parent of the working directory that holds
//tests/tests.toml. Empty when neither is found.
inline std::filesystem::path nos_test_find_repo_root()
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
	return {};
}

//The root, or the working directory when there is none, so a failing test names the path it tried
inline std::filesystem::path nos_test_repo_root()
{
	const auto root = nos_test_find_repo_root();
	return root.empty() ? std::filesystem::current_path() : root;
}

//The working directory the tests are written for: tests/src, which ctest sets
//(tests/src/CMakeLists.txt). The geometry-aid and SALTED-model cases open their inputs relative to it.
inline std::filesystem::path nos_test_expected_cwd()
{
	return nos_test_find_repo_root() / "tests" / "src";
}

// error_check ends the process with exit(-1). Windows reports the full value back
// through GetExitCodeProcess, while POSIX wait() only carries the low 8 bits, so the
// same exit is observed as 0xFFFFFFFF on one and 255 on the other.
#ifdef _WIN32
constexpr unsigned ERROR_CHECK_EXIT_CODE = static_cast<unsigned>(-1);
#else
constexpr unsigned ERROR_CHECK_EXIT_CODE = 255u;
#endif

#endif //PCH_H
