#pragma once
//Per-analysis result files next to the wavefunction: <stem>.rgbi_log, .qtaim_log, .eli_log and
//.nbo_log. The run log keeps everything it always had, byte for byte; a section's result blocks are
//copied into its file as well, the way tee does it, under the program header, a title block and the
//references of that analysis. Progress chatter, and citation lines the header already lists, stay in
//the run log only.
#include "convenience.h"
#include "citations.h"
#include <ctime>
#include <filesystem>
#include <fstream>
#include <vector>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

namespace section_log
{
	//The run log's name, for the pointer in every section file's header
	inline std::string &main_log() { static std::string name = "NoSpherA2.log"; return name; }

	struct open_file
	{
		std::ostream *f;
		std::set<std::string> cited;  //the header's reference lines, as citations::cite prints them
	};
	inline std::map<std::string, open_file> &open_files() { static std::map<std::string, open_file> m; return m; }

	inline std::string rule(const char c) { return std::string(80, c) + "\n"; }

	//Writes to the stream it replaces and, a line at a time and filtered, to the section file
	class teebuf : public std::streambuf
	{
		std::streambuf *main;
		open_file &file;
		std::string line;
		void emit()
		{
			const bool chatter = line.find('\r') != std::string::npos || line.rfind("GridManager:", 0) == 0
				|| line.rfind("Quadrature points sent", 0) == 0;
			if (!chatter && !file.cited.count(line)) *file.f << line << '\n';
			line.clear();
		}
		void put(const char c)
		{
			if (c == '\n') emit();
			else line += c;
		}
	protected:
		int overflow(int c) override
		{
			if (c == traits_type::eof()) return traits_type::not_eof(c);
			put(static_cast<char>(c));
			return main->sputc(static_cast<char>(c));
		}
		std::streamsize xsputn(const char *s, std::streamsize n) override
		{
			for (std::streamsize i = 0; i < n; i++) put(s[i]);
			return main->sputn(s, n);
		}
		int sync() override { file.f->flush(); return main->pubsync(); }
	public:
		teebuf(std::streambuf *m, open_file &f) : main(m), file(f) {}
		~teebuf() override { if (!line.empty()) emit(); }
	};

	inline void heading(std::ostream &f, const std::string &title)
	{
		f << "\n" << rule('-') << "  " << title << "\n" << rule('-');
	}

	//A heading in the open section file of this kind only; the run log does not see it
	inline void heading(const std::string &kind, const std::string &title)
	{
		const auto it = open_files().find(kind);
		if (it != open_files().end()) heading(*it->second.f, title);
	}

	//While alive, std::cout also feeds the open section file of this kind; nothing when none is open
	class tee
	{
		std::streambuf *saved = nullptr;
		teebuf *buf = nullptr;
	public:
		explicit tee(const std::string &kind, const std::string &title = "")
		{
			const auto it = open_files().find(kind);
			if (it == open_files().end()) return;
			std::cout.flush();
			if (!title.empty()) heading(*it->second.f, title);
			buf = new teebuf(std::cout.rdbuf(), it->second);
			saved = std::cout.rdbuf(buf);
		}
		~tee()
		{
			if (!buf) return;
			std::cout.flush();
			std::cout.rdbuf(saved);
			delete buf;
		}
		tee(const tee &) = delete;
		tee &operator=(const tee &) = delete;
	};

	//Opens <wfn stem>.<kind>_log beside the wavefunction and writes its header. Registered sections
	//are what tee finds; an unregistered one (the NBO thread's) is written through stream() only.
	class section
	{
		std::string kind;
		std::ofstream f;
		bool no_date, registered;
		static std::string now()
		{
			const std::time_t t = std::time(nullptr);
			std::tm tm{};
#ifdef _WIN32
			localtime_s(&tm, &t);
#else
			localtime_r(&t, &tm);
#endif
			char s[32];
			std::strftime(s, sizeof s, "%Y-%m-%d %H:%M:%S", &tm);
			return s;
		}
	public:
		section(const std::filesystem::path &wfn, const std::string &k, const std::string &title,
			const std::vector<citations::Method> &refs, const bool nodate, const bool reg = true)
			: kind(k), f(wfn.parent_path() / (wfn.stem().string() + "." + k + "_log")), no_date(nodate), registered(reg)
		{
			f << NoSpherA2_message(no_date);
			if (!no_date) f << build_date;
			f << "\n" << rule('=') << "  " << title << "\n" << rule('=')
				<< "  Wavefunction : " << wfn.filename().string() << "\n";
			if (!no_date) f << "  Started      : " << now() << "\n";
			f << "  Full run log : " << main_log() << " (this file repeats the result sections of it)\n"
				<< "\n  References for this analysis:\n";
			std::set<std::string> cited;
			for (const citations::Method m : refs) {
				std::ostringstream c;
				citations::cite(m, c);
				std::istringstream lines(c.str());
				for (std::string l; std::getline(lines, l);) {
					f << "    " << l << "\n";
					cited.insert(l);
				}
			}
			if (registered) open_files()[kind] = open_file{ &f, cited };
		}
		~section()
		{
			if (registered) open_files().erase(kind);
			f << "\n" << rule('=') << "  End of " << kind << " results";
			if (!no_date) f << " (" << now() << ")";
			f << "\n" << rule('=');
		}
		std::ostream &stream() { return f; }
		section(const section &) = delete;
		section &operator=(const section &) = delete;
	};
}
