#include "pch.h"
#include "core/tuning.h"

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

//Report an exit that interrupts a running test.
static void report_exit_during_test()
{
	const ::testing::TestInfo* info = ::testing::UnitTest::GetInstance()->current_test_info();
	if (info == nullptr)
		return;
	std::cout.flush();
	std::cerr << "\nNoSpherA2_Tests: the process exited while " << info->test_suite_name() << "."
			  << info->name() << " was still running. Production code called exit() (error_check "
			  << "does) outside a death test, so no verdict was printed and the tests after it "
			  << "never ran.\n";
	std::cerr.flush();
}

int main(int argc, char** argv)
{
	//A death-test child is expected to exit.
	const bool death_test_child = std::any_of(argv, argv + argc, [](const char* a) {
		return std::string(a).rfind("--gtest_internal_run_death_test", 0) == 0;
	});
	if (!death_test_child)
		std::atexit(report_exit_during_test);

	//The "fast" death test style forks: the child owns copies of every thread object but none of the
	//threads, so its exit path can throw (pthread_detach) and abort with signal 6 instead of the expected
	//exit code.  "threadsafe" re-executes the binary instead.  Set before InitGoogleTest so
	//--gtest_death_test_style on the command line still wins.
	GTEST_FLAG_SET(death_test_style, "threadsafe");
	//The RGBI tests count free-atom SCFs and compare their digits run against run; a density read back
	//from a previous run's on-disk cache would answer them without running anything. A test that wants
	//the disk cache sets NOS_FREEATOM_CACHE_DIR itself.
	set_tuning("NOS_FREEATOM_CACHE_DIR", "off");
	::testing::InitGoogleTest(&argc, argv);
	::testing::UnitTest::GetInstance()->listeners().Append(new CoutFormatReset);
	warn_about_working_directory();
	return RUN_ALL_TESTS();
}
