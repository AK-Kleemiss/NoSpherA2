#include "pch.h"

#include "core/convenience.h"
#include "core/NoSpherA2.h"

//An analysis that cannot run, or an option nothing reads, must exit non-zero: exit 0 looks like success
namespace {

std::filesystem::path scratch(const std::string& name)
{
	const auto dir = std::filesystem::temp_directory_path() /
					 ("nos_cli_refusal_" + name + "_" +
					  std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
	std::filesystem::create_directories(dir);
	return dir;
}

//run_app() logs std::cout to NoSpherA2.log in the cwd and restores the buffer it found, so the
//capture sees the final message and the log lands in the scratch
struct Cli
{
	std::filesystem::path dir;
	std::filesystem::path previous_cwd;
	std::ostringstream console;
	std::streambuf* previous_cout;

	explicit Cli(const std::string& name)
		: dir(scratch(name)), previous_cwd(std::filesystem::current_path()),
		  previous_cout(std::cout.rdbuf(console.rdbuf()))
	{
		std::filesystem::current_path(dir);
	}
	~Cli()
	{
		std::cout.rdbuf(previous_cout);
		std::error_code ec;
		std::filesystem::current_path(previous_cwd, ec);
		std::filesystem::remove_all(dir, ec);
	}
	int run(const std::vector<std::string>& args)
	{
		std::vector<std::string> owned{"NoSpherA2"};
		owned.insert(owned.end(), args.begin(), args.end());
		std::vector<char*> argv;
		for (auto& a : owned) argv.push_back(a.data());
		int argc = static_cast<int>(argv.size());
		return run_app(argc, argv.data());
	}
	std::string output() const { return console.str(); }
};

//The parser exits through err_checkf, hence death tests; an unknown option is refused before any input is read
void parse(const std::vector<std::string>& args)
{
	std::vector<std::string> owned{"NoSpherA2"};
	owned.insert(owned.end(), args.begin(), args.end());
	std::vector<char*> argv;
	for (auto& a : owned) argv.push_back(a.data());
	int argc = static_cast<int>(argv.size());
	std::ostringstream log;
	options opt(argc, argv.data(), log);
	opt.digest_options();
}

std::filesystem::path fixture(const std::string& rel)
{
	return nos_test_repo_root() / "tests" / rel;
}

} // namespace

//RGBI has no positional form: `-rgbi water.gbw` sets the flag and leaves opt.wfn empty
TEST(CliRefusal, RgbiWithoutWavefunctionExitsNonZeroAndNamesTheInput)
{
	Cli cli("rgbi");
	const int rc = cli.run({"-rgbi", fixture("RGBI_groups/nh3li.gbw").string()});
	EXPECT_NE(rc, 0);
	const std::string out = cli.output();
	EXPECT_NE(out.find("-rgbi"), std::string::npos) << out;
	EXPECT_NE(out.find("-wfn"), std::string::npos) << out;
}

//NPA shares the RGBI branch
TEST(CliRefusal, NpaWithoutWavefunctionExitsNonZero)
{
	Cli cli("npa");
	EXPECT_NE(cli.run({"-npa", fixture("RGBI_groups/nh3li.gbw").string()}), 0);
	EXPECT_NE(cli.output().find("-npa"), std::string::npos) << cli.output();
}

//The message is a pure function of the options, so its wording is checked without running anything
TEST(CliRefusal, UnrunnableAnalysisNamesTheAnalysisAndTheOptionItWanted)
{
	options opt;
	EXPECT_TRUE(opt.unrunnable_analysis().empty());
	opt.rgbi = true;
	const std::string msg = opt.unrunnable_analysis();
	EXPECT_NE(msg.find("-rgbi"), std::string::npos) << msg;
	EXPECT_NE(msg.find("-wfn"), std::string::npos) << msg;
	opt.wfn = "some.gbw";
	EXPECT_TRUE(opt.unrunnable_analysis().empty());
	opt.wfn.clear();
	opt.rgbi = false;
	opt.npa = true;
	EXPECT_NE(opt.unrunnable_analysis().find("-npa"), std::string::npos);
	opt.npa = false;
	opt.hkl = "some.hkl";
	EXPECT_FALSE(opt.unrunnable_analysis().empty());
}

//one misspelling per analysis family
TEST(CliRefusal, MisspelledOptionInAnAnalysisFamilyIsFatal)
{
	EXPECT_EXIT(parse({"-rgbi_gruops", "0,1"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-npa_orbitalz"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-nrt_maxx", "5"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-nbo_e2mn", "0.5"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-eli_familly", "x.gbw"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-basin_gridd"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-topologyy", "x.gbw"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//-nbo/-nrt options are read by the -nbo_native handler from the tokens after its wavefunction, out of
//sight of the check above; a real fixture, since the run must die on the option before computing
TEST(CliRefusal, MisspelledNrtOptionAfterNboNativeIsFatal)
{
	const auto wfn = fixture("RGBI_groups/nh3li.gbw");
	ASSERT_TRUE(std::filesystem::exists(wfn));
	EXPECT_EXIT(parse({"-nbo_native", wfn.string(), "-nrt_maxx", "5"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//A lone hydrogen's empty spin channel gives a rank-zero candidate (a 0x0 eigensolve); it must be
//infeasible so the parent-structure check refuses the run
TEST(CliRefusal, NrtOnAOneElectronWavefunctionRefusesInsteadOfCrashing)
{
	const auto wfn = fixture("ptb_H_file/H.gbw");
	if (!std::filesystem::exists(wfn))
		GTEST_SKIP() << wfn.string() << " not found";
	EXPECT_EXIT(parse({"-nbo_native", wfn.string(), "-nrt"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//the positional wavefunction of the ELI and topology analyses is checked like -wfn
TEST(CliRefusal, PositionalWavefunctionThatDoesNotExistIsFatal)
{
	EXPECT_EXIT(parse({"-topology", "no_such_wavefunction.gbw"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-eli_family", "no_such_wavefunction.gbw"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-eli_analysis", "no_such_wavefunction.gbw", "0.2", "2.0"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-qtaim_eli", "no_such_wavefunction.gbw", "0,1"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//An unclaimed family flag is fatal, and the -nbo/-nbo_native handlers read their -nbo_*/-nrt_* options
//where no digester sees them, so those must be listed in nbo_family_suboptions().  Scanning the parser
//source catches a new flag that is neither digested nor listed before it aborts real runs
TEST(CliRefusal, AnalysisFlagsAreEitherDigestedOrListed)
{
	const auto src = nos_test_repo_root() / "Src" / "core" / "convenience.cpp";
	std::ifstream in(src);
	ASSERT_TRUE(in.good()) << src;
	const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

	//every "-flag" string literal in the file, and separately the ones a digester compares against
	std::set<std::string> mentioned, digested;
	for (size_t p = text.find("\"-"); p != std::string::npos; p = text.find("\"-", p + 1))
	{
		const size_t end = text.find('"', p + 1);
		if (end == std::string::npos) break;
		const std::string flag = text.substr(p + 1, end - p - 1);
		if (flag.size() < 2 || !isalpha(static_cast<unsigned char>(flag[1]))) continue;
		if (flag.find_first_of(" \t<>") != std::string::npos) continue; // help text, not a flag
		if (!owning_analysis(flag)) continue;
		size_t before = p;
		while (before > 0 && (text[before - 1] == ' ' || text[before - 1] == '\t')) before--;
		if (before > 0 && text[before - 1] == '{')
			continue; // {"-basin", "QTAIM basin"} - owning_analysis' own table of prefixes
		mentioned.insert(flag);
		//"temp == \"-flag\"" on the same line - the digesters' only way of claiming one
		const size_t line_start = text.rfind('\n', p) + 1;
		if (text.substr(line_start, p - line_start).find("temp ==") != std::string::npos)
			digested.insert(flag);
	}
	ASSERT_GT(mentioned.size(), 10u) << "the scan found nothing - did the parser move?";

	std::vector<std::string> orphans;
	for (const auto& flag : mentioned)
		if (digested.count(flag) == 0 && nbo_family_suboptions().count(flag) == 0)
			orphans.push_back(flag);
	EXPECT_TRUE(orphans.empty())
		<< "these flags name an analysis family but no digester compares against them and they are "
		   "not in nbo_family_suboptions(), so digest_options will refuse them as misspellings: "
		<< [&] { std::string s; for (const auto& o : orphans) s += o + " "; return s; }();
}

//A real option the running analysis never reads: -nbo_* belong to the -nbo_native handler, and
//run_app_impl's early-exit analyses return before the RGBI/NPA block, so the unknown-option check misses both
TEST(CliRefusal, RealOptionNoAnalysisReadsIsFatal)
{
	//-nrt and -nbo_threads are read by the NBO handlers, and this line runs RGBI
	EXPECT_EXIT(parse({"-rgbi", "-nrt"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-rgbi", "-nbo_threads", "2"}), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	//an analysis that ends the run, plus one that sits after it in the chain
	const auto wfn = fixture("ptb_H_file/H.gbw");
	if (!std::filesystem::exists(wfn))
		GTEST_SKIP() << wfn.string() << " not found";
	EXPECT_EXIT(parse({"-eli_analysis", wfn.string(), "0.3", "2.0", "-rgbi_basis", "nao"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-eli_analysis", wfn.string(), "0.3", "2.0", "-npa"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	//a property cube skips the wavefunction block too
	EXPECT_EXIT(parse({"-wfn", wfn.string(), "-lap", "-rgbi"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	//-fba runs its own NPA; -calc_F stops -do_XCW before the fit that runs RGBI
	EXPECT_EXIT(parse({"-fba", wfn.string(), "-npa_summary"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({"-do_XCW", "-calc_F", "-rgbi"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	//-fba reads the -rgbi_* modifiers, each of which also sets rgbi
	parse({"-fba", wfn.string(), "-rgbi_no_sym", "-rgbi_basis", "nao"});
}

//-nrt and -nbo_json belong to the -nbo_native handler, so no top-level digester claims them; they
//must still parse in any order
TEST(CliRefusal, RealOptionsStillParse)
{
	parse({"-rgbi", "-rgbi_no_sym", "-rgbi_EVs"});
	parse({"-nrt", "-nrt_exhaustive", "-nbo_json", "out.json", "-nbo_threads", "4"});
	parse({"-basin_grid", "3", "-basin_cube"});
	parse({"-rgbi", "-rgbi_theta", "-rgbi_legacy_cutoff"});
	SUCCEED();
}
