#include "pch.h"

#include "core/citations.h"

#include <algorithm>
#include <regex>
#include <set>
#include <sstream>
#include <string>

// The table is typed by hand: a blank or malformed DOI still prints and sends the reader nowhere
TEST(Citations, EveryEntryCarriesAResolvableDoi)
{
	// DOI handbook shape 10.<registrant>/<suffix>: 4-9 digits, non-empty suffix without whitespace
	const std::regex doi_shape("^10\\.[0-9]{4,9}/\\S+$");
	ASSERT_FALSE(citations::table().empty()) << "the citation table is empty";
	for (const citations::Reference &r : citations::table())
	{
		ASSERT_NE(r.tag, nullptr);
		ASSERT_NE(r.work, nullptr);
		ASSERT_NE(r.doi, nullptr);
		const std::string tag(r.tag), work(r.work), doi(r.doi);
		EXPECT_FALSE(tag.empty()) << "empty tag next to DOI " << doi;
		EXPECT_FALSE(work.empty()) << "empty reference text next to DOI " << doi;
		EXPECT_FALSE(doi.empty()) << "empty DOI for [" << tag << "] " << work;
		EXPECT_TRUE(std::regex_match(doi, doi_shape)) << "not a DOI: '" << doi << "' for [" << tag << "]";
		// a https:// or doi: prefix already fails the regex; this names the reason
		EXPECT_EQ(doi.rfind("10.", 0), 0u) << "DOI should be bare, no prefix: " << doi;
	}
}

TEST(Citations, FormatsOneLinePerReference)
{
	const citations::Reference bader{citations::Method::QTAIM, "QTAIM",
									 "Bader, Chem. Rev. 91 (1991) 893", "10.1021/cr00005a013"};
	EXPECT_EQ(citations::format(bader),
			  "[QTAIM] Bader, Chem. Rev. 91 (1991) 893, DOI 10.1021/cr00005a013");

	// cite() emits every reference of the method and nothing else
	std::ostringstream os;
	citations::cite(citations::Method::QTAIM, os);
	const citations::Reference book{citations::Method::QTAIM, "QTAIM",
									 "Bader, Atoms in Molecules: A Quantum Theory, Oxford University Press (1990)", "10.1093/oso/9780198551683.001.0001"};
	EXPECT_EQ(os.str(), citations::format(book) + "\n" + citations::format(bader) + "\n") << "QTAIM rests on the book and the review";

	std::ostringstream two;
	citations::cite(citations::Method::HAR, two);
	//str() returns a copy: iterators from two calls point into different temporaries
	const std::string har = two.str();
	EXPECT_EQ(std::count(har.begin(), har.end(), '\n'), 2) << "HAR rests on two papers";
	EXPECT_NE(har.find("10.1107/S0108767308005709"), std::string::npos);
	EXPECT_NE(har.find("10.1107/S2052252514014845"), std::string::npos);
}

// A reader cites inside the caller's open "Reading: <file> ... done!" line, so it queues: printing
// there lands mid-line and shifts every following log line
TEST(Citations, QueueDefersUntilFlushAndThenForgets)
{
	std::ostringstream drain;
	citations::flush(drain); // whatever an earlier test queued, so this one starts empty
	drain.str("");

	citations::queue(citations::Method::PTB);
	citations::queue(citations::Method::Molden);
	EXPECT_TRUE(drain.str().empty()) << "queue() must not write anything";

	std::ostringstream os;
	citations::flush(os);
	const std::string out = os.str();
	EXPECT_NE(out.find("10.1063/5.0137838"), std::string::npos) << "pTB reference missing";
	EXPECT_NE(out.find("10.1023/A:1008193805436"), std::string::npos) << "molden reference missing";
	EXPECT_EQ(out.front(), '[') << "a flushed citation starts its own line";

	std::ostringstream again;
	citations::flush(again);
	EXPECT_TRUE(again.str().empty()) << "flush() must clear the queue";
}

// Two methods may share a DOI (.tsc format and program are one paper); the same method twice is a copy-paste slip
TEST(Citations, NoDuplicateEntries)
{
	std::set<std::string> seen;
	for (const citations::Reference &r : citations::table())
	{
		const std::string key = std::to_string(static_cast<int>(r.method)) + "|" + r.doi;
		EXPECT_TRUE(seen.insert(key).second) << "duplicate entry: " << citations::format(r);
	}
}
