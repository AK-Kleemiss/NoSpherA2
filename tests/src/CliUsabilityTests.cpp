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

} // namespace

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
