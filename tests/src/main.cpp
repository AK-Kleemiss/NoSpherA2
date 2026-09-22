#include "pch.h"

//Say where the executable should run before any test fails on a relative path
static void warn_about_working_directory()
{
	const auto cwd = std::filesystem::current_path();
	const auto root = nos_test_find_repo_root();
	if (root.empty()) {
		std::cerr << "NoSpherA2_Tests: no repository found, neither " << cwd.string()
				  << " nor its parents hold tests/tests.toml.\n"
				  << "  Run the tests through ctest (ctest --preset <preset>), start the executable in "
				  << "<repository>/tests/src, or set NOS_REPO_ROOT=<repository>.\n";
		return;
	}
	const auto expected = nos_test_expected_cwd();
	if (!std::filesystem::is_directory(expected) || !std::filesystem::equivalent(cwd, expected)) {
		std::cerr << "NoSpherA2_Tests: started in " << cwd.string() << ", but the tests expect "
				  << expected.string() << " as working directory (ctest sets it); "
				  << "cases that open their inputs relative to it will fail here.\n";
	}
}

//Production code leaves std::fixed / setprecision on std::cout (run_app restores it, a direct call does not); a
//test that then checks printed numbers would read the previous test's format. Reset before every test.
struct CoutFormatReset : ::testing::EmptyTestEventListener
{
	void OnTestStart(const ::testing::TestInfo&) override
	{
		std::cout.flags(std::ios::dec | std::ios::skipws);
		std::cout.precision(6);
		std::cout.width(0);
		std::cout.clear();
	}
};

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	::testing::UnitTest::GetInstance()->listeners().Append(new CoutFormatReset);
	warn_about_working_directory();
	return RUN_ALL_TESTS();
}
