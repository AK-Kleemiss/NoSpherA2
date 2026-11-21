#include "pch.h"

#include "core/base/convenience.h"
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

//err_checkf writes to cout or the log and exits; a death test matches stderr only
void die_on(const std::vector<std::string>& args)
{
	std::vector<std::string> owned{"NoSpherA2"};
	owned.insert(owned.end(), args.begin(), args.end());
	std::vector<char*> argv;
	for (auto& a : owned) argv.push_back(a.data());
	int argc = static_cast<int>(argv.size());
	std::cout.rdbuf(std::cerr.rdbuf());
	options opt(argc, argv.data(), std::cerr);
	opt.digest_options();
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

//-nbo-json is -nbo_json, as -no-date is -no_date: the run gets past the sub-option check to the missing file
TEST(CliUsability, NboSubOptionsTakeDashesLikeEveryFlag)
{
	EXPECT_EXIT(die_on({"-nbo", "no_such_file.gbw", "-nrt-max", "5", "-nbo-json", "x.json"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), "wavefunction doesn't exist: no_such_file.gbw");
}

//the -nbo handler reads only what follows its wavefunction; one written before it is dropped, and says so
TEST(CliUsability, NboSubOptionBeforeItsParentIsReported)
{
	const Parsed p = digest({"-nrt", "-nbo_threads", "2"});
	EXPECT_NE(p.err.find("WARNING: ignoring -nrt:"), std::string::npos) << p.err;
	EXPECT_NE(p.err.find("WARNING: ignoring -nbo_threads:"), std::string::npos) << p.err;
	EXPECT_EQ(p.log.find("WARNING"), std::string::npos) << p.log;
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

//the nearest spelling -h lists; the log line stays as it was
TEST(CliUsability, UnknownOptionGetsTheNearestSpelling)
{
	Parsed p = digest({"-hkll", "x.hkl"});
	EXPECT_NE(p.err.find("-hkll (did you mean -hkl?)"), std::string::npos) << p.err;
	EXPECT_EQ(p.log.find("did you mean"), std::string::npos) << p.log;
	p = digest({"-Cif", "x.cif"});
	EXPECT_NE(p.err.find("(did you mean -cif?)"), std::string::npos) << p.err;
	p = digest({"-completely_unrelated"});
	EXPECT_EQ(p.err.find("did you mean"), std::string::npos) << p.err;
	EXPECT_EXIT(die_on({"-rgbi_gruops", "0,1"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), "did you mean -rgbi_groups");
}

TEST(CliUsability, KnownOptionsPrintNothingToTheConsole)
{
	const Parsed p = digest({"-acc", "1", "-no_date", "-all_charges"});
	EXPECT_TRUE(p.err.empty()) << p.err;
}

//-h <word> prints the sections whose title holds the word, and names the sections when none does
TEST(CliUsability, HelpSectionPrintsOnlyThatSection)
{
	const std::string gpu = help_section("gpu");
	EXPECT_NE(gpu.find("GPU ACCELERATION"), std::string::npos) << gpu;
	EXPECT_NE(gpu.find("-gpu_cublas"), std::string::npos) << gpu;
	EXPECT_EQ(gpu.find("BONDING ANALYSIS"), std::string::npos) << gpu;
	EXPECT_LT(gpu.size(), help_message.size() / 4);
	const std::string none = help_section("no_such_section");
	EXPECT_NE(none.find("No help section matches"), std::string::npos) << none;
	EXPECT_NE(none.find("BONDING ANALYSIS"), std::string::npos) << none;
}

//every flag a digester compares against is in -h; case variants and old names count as listed
TEST(CliUsability, EveryParsedFlagIsInTheHelp)
{
	const auto src = nos_test_repo_root() / "Src" / "core" / "options_parse.cpp";
	std::ifstream in(src);
	ASSERT_TRUE(in.good()) << src;
	const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	auto lower = [](std::string s) { for (char& c : s) c = (char)tolower((unsigned char)c); return s; };
	const std::string help = lower(help_message);
	const std::set<std::string> aliases = {"-cube", "-nbo_native", "-lahvatest"};
	std::vector<std::string> missing;
	int found = 0;
	for (size_t p = text.find("temp == \"-"); p != std::string::npos; p = text.find("temp == \"-", p + 1))
	{
		const size_t b = p + 9, e = text.find('"', b);
		const std::string flag = lower(text.substr(b, e - b));
		found++;
		bool listed = false;
		for (size_t q = help.find(flag); q != std::string::npos && !listed; q = help.find(flag, q + 1))
			listed = q + flag.size() == help.size() || !(isalnum((unsigned char)help[q + flag.size()]) || help[q + flag.size()] == '_');
		if (!listed && !aliases.count(text.substr(b, e - b)))
			missing.push_back(text.substr(b, e - b));
	}
	ASSERT_GT(found, 100) << "the scan found nothing - did the parser move?";
	EXPECT_TRUE(missing.empty()) << [&] { std::string s; for (const auto& m : missing) s += m + " "; return s; }();
}

//no job: the help goes to the log as before, and the console is told where it went
TEST(CliUsability, NoJobTellsTheConsole)
{
	const std::filesystem::path log = std::filesystem::temp_directory_path() / "CliUsability_nojob.log";
	std::vector<std::string> args = {"NoSpherA2", "-out", log.string(), "-no_date"};
	std::vector<char*> argv;
	for (std::string& a : args) argv.push_back(a.data());
	std::ostringstream err;
	std::streambuf* previous = std::cerr.rdbuf(err.rdbuf());
	const int rc = run_app((int)args.size(), argv.data());
	std::cerr.rdbuf(previous);
	std::filesystem::remove(log);
	EXPECT_EQ(rc, 0);
	EXPECT_NE(err.str().find("Did not understand the task to perform"), std::string::npos) << err.str();
	EXPECT_NE(err.str().find(log.string()), std::string::npos) << err.str();
}

//-tune NAME is set whatever NAME is, as before; one no source documents is reported
TEST(CliUsability, UnknownTuningKnobIsReportedAndStillSet)
{
	EXPECT_TRUE(tuning_known("NOSPHERA2_XCW_RCHECK"));
	EXPECT_TRUE(tuning_known("NOSPHERA2_ERI_DYNAMIC"));
	EXPECT_FALSE(tuning_known("NOSPHERA2_NO_SUCH_KNOB"));
	Parsed p = digest({"-tune", "NOSPHERA2_NO_SUCH_KNOB=2"});
	ASSERT_NE(tuning("NOSPHERA2_NO_SUCH_KNOB"), nullptr);
	EXPECT_STREQ(tuning("NOSPHERA2_NO_SUCH_KNOB"), "2");
	tuning_end_run();
	EXPECT_NE(p.err.find("WARNING: -tune NOSPHERA2_NO_SUCH_KNOB"), std::string::npos) << p.err;
	EXPECT_EQ(p.log.find("NOSPHERA2_NO_SUCH_KNOB"), std::string::npos) << p.log;
	p = digest({"-tune", "NOSPHERA2_ERI_DYNAMIC"});
	tuning_end_run();
	EXPECT_TRUE(p.err.empty()) << p.err;
}

//-tune list prints the documented knobs and stops before any input is read
TEST(CliUsability, TuneListPrintsTheKnobsAndStops)
{
	std::ostringstream out;
	tuning_list(out);
	EXPECT_NE(out.str().find("NOSPHERA2_ERI_DYNAMIC"), std::string::npos) << out.str();
	EXPECT_NE(out.str().find("NOSPHERA2_GPU_BACKEND=cuda|hip"), std::string::npos) << out.str();
	EXPECT_EXIT(digest({"-tune", "list", "-cif", "no_such_file.cif"}), ::testing::ExitedWithCode(0), "");
}

//a one-shot option ends the command line: what follows is named, except what is read wherever it stands
TEST(CliUsability, FlagsAfterAOneShotOptionAreReported)
{
	const Parsed p = digest({"-acc", "1", "-lahvatest", "-cif", "x.cif", "-tune", "NOSPHERA2_ERI_DYNAMIC", "-no_date", "-cpus", "2"});
	tuning_end_run();
	const size_t at = p.err.find("WARNING: ignored after -lahvatest");
	ASSERT_NE(at, std::string::npos) << p.err;
	const std::string line = p.err.substr(at, p.err.find('\n', at) - at);
	EXPECT_NE(line.find(": -cif -cpus"), std::string::npos) << line;
	EXPECT_EQ(line.find("-tune"), std::string::npos) << line;
	EXPECT_EQ(line.find("-no_date"), std::string::npos) << line;
	EXPECT_EQ(line.find("-acc"), std::string::npos) << line;
	EXPECT_EQ(p.log.find("ignored after"), std::string::npos) << p.log;
	EXPECT_TRUE(digest({"-acc", "1", "-lahvatest"}).err.empty());
}
