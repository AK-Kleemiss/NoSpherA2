#include "pch.h"

#include "core/convenience.h"
#include "core/NoSpherA2.h"

//What the command line tells a user who typed it wrong; a valid line must print exactly what it did before
namespace {

struct Parsed
{
	std::string log, err;
};

//digest_options with stderr captured; the exe points cout at its log, so stderr is what reaches the console
Parsed digest(const std::vector<std::string>& args)
{
	std::vector<std::string> owned{"NoSpherA2"};
	owned.insert(owned.end(), args.begin(), args.end());
	std::vector<char*> argv;
	for (auto& a : owned) argv.push_back(a.data());
	int argc = static_cast<int>(argv.size());
	std::ostringstream log, err;
	std::streambuf* previous = std::cerr.rdbuf(err.rdbuf());
	options opt(argc, argv.data(), log);
	opt.digest_options();
	std::cerr.rdbuf(previous);
	return {log.str(), err.str()};
}

//err_checkf writes to cout and exits; a death test matches stderr only
void die_on(const std::vector<std::string>& args)
{
	std::cout.rdbuf(std::cerr.rdbuf());
	digest(args);
}

} // namespace

//The refusal names the flag and the file it could not find
TEST(CliUsability, MissingInputFileIsNamed)
{
	const auto code = ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE);
	EXPECT_EXIT(die_on({"-cif", "no_such_file.cif"}), code, "-cif: CIF doesn't exist: no_such_file.cif");
	EXPECT_EXIT(die_on({"-hkl", "no_such_file.hkl"}), code, "-hkl: hkl file doesn't exist: no_such_file.hkl");
	EXPECT_EXIT(die_on({"-wfn", "no_such_file.gbw"}), code, "-wfn: wavefunction does not exist: no_such_file.gbw");
	EXPECT_EXIT(die_on({"-coef", "no_such_file.npy"}), code, "-coef: coefficient file doesn't exist: no_such_file.npy");
	EXPECT_EXIT(die_on({"-cube", "no_such_file.cube"}), code, "-cube: cube file does not exist: no_such_file.cube");
	EXPECT_EXIT(die_on({"-tscb", "table.txt"}), code, "-tscb takes a .tsc or .tscb file, not table.txt");
}

//The hint names the option that asked, not the one its text was copied from
TEST(CliUsability, MissingPrerequisiteNamesTheOptionThatNeedsIt)
{
	const auto code = ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE);
	EXPECT_EXIT(die_on({"-RI_CUBE", "c.npy"}), code, "BEFORE -RI_CUBE");
	EXPECT_EXIT(die_on({"-SALTED_Training"}), code, "BEFORE -SALTED_Training");
	EXPECT_EXIT(die_on({"-lukas_test"}), code, "BEFORE -lukas_test");
}

TEST(CliUsability, UnknownOptionReachesTheConsoleAndTheLog)
{
	const Parsed p = digest({"-hkll", "x.hkl"});
	EXPECT_NE(p.log.find("WARNING: ignoring unknown option -hkll"), std::string::npos) << p.log;
	EXPECT_NE(p.err.find("WARNING: ignoring unknown option -hkll"), std::string::npos) << p.err;
}

TEST(CliUsability, KnownOptionsPrintNothingToTheConsole)
{
	const Parsed p = digest({"-acc", "1", "-no_date", "-all_charges"});
	EXPECT_TRUE(p.err.empty()) << p.err;
}
