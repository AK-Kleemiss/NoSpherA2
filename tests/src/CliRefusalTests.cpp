#include "pch.h"

#include "core/convenience.h"
#include "core/NoSpherA2.h"

//A command line that names an analysis and cannot run it must fail, loudly. Every case here used
//to end in "Did not understand the task to perform!" written into NoSpherA2.log with exit code 0,
//or in the option being dropped and the analysis running with a default nobody asked for - both
//read, from the outside, exactly like a successful run.
namespace {

std::filesystem::path scratch(const std::string& name)
{
	const auto dir = std::filesystem::temp_directory_path() /
					 ("nos_cli_refusal_" + name + "_" +
					  std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
	std::filesystem::create_directories(dir);
	return dir;
}

//run_app() redirects std::cout into NoSpherA2.log in the working directory and restores the buffer
//it found on entry, so the caller sees the final message and the log lands in the scratch.
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

//The parser dies through err_checkf, so the unknown-option cases are death tests. They need no
//input file: a misspelled option is refused before anything is read.
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

//RGBI has no positional form: `-rgbi water.gbw` sets the flag, leaves opt.wfn empty, and the
//branch that would run the analysis is skipped. This is the case reported from the cluster share.
TEST(CliRefusal, RgbiWithoutWavefunctionExitsNonZeroAndNamesTheInput)
{
	Cli cli("rgbi");
	const int rc = cli.run({"-rgbi", fixture("RGBI_groups/nh3li.gbw").string()});
	EXPECT_NE(rc, 0);
	const std::string out = cli.output();
	EXPECT_NE(out.find("-rgbi"), std::string::npos) << out;
	EXPECT_NE(out.find("-wfn"), std::string::npos) << out;
}

//NPA sits in the same branch and had the same silence
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

//One misspelling per analysis family. Each of these used to be dropped without a word.
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

//The -nbo/-nrt options are read by the -nbo_native handler itself, from the tokens after its
//wavefunction, so the check above cannot see them: they need their own refusal. A real fixture,
//because the point is that the run dies on the option instead of spending minutes and exiting 0.
TEST(CliRefusal, MisspelledNrtOptionAfterNboNativeIsFatal)
{
	const auto wfn = fixture("RGBI_groups/nh3li.gbw");
	ASSERT_TRUE(std::filesystem::exists(wfn));
	EXPECT_EXIT(parse({"-nbo_native", wfn.string(), "-nrt_maxx", "5"}),
				::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//-wfn refuses a file that is not there; the positional forms of the ELI and topology analyses
//took the name on trust
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

//The refusal above is a trap for whoever adds the next option: a flag in one of the analysis
//families that no digester claims is now fatal. Almost every option is safe by construction - a
//digester claims it, and a claimed flag never reaches the refusal - but the -nbo/-nbo_native
//handlers read their own -nbo_*/-nrt_* options from the tokens after their wavefunction, where no
//digester sees them, and those need an entry in nbo_family_suboptions(). This reads the parser's
//own source and asserts every family flag it mentions is one of the two, so a new option added the
//second way cannot reach a release aborting runs that used to work.
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

//The other half of the check: an option that is real must still parse. -nrt and -nbo_json are
//consumed by the -nbo_native handler, so at the top level no digester claims them - they must not
//be mistaken for typos whichever order they were written in.
TEST(CliRefusal, RealOptionsStillParse)
{
	parse({"-rgbi", "-rgbi_no_sym", "-rgbi_EVs"});
	parse({"-nrt", "-nrt_exhaustive", "-nbo_json", "out.json", "-nbo_threads", "4"});
	parse({"-basin_grid", "3", "-basin_cube"});
	parse({"-rgbi", "-rgbi_theta", "-rgbi_legacy_cutoff"});
	SUCCEED();
}
