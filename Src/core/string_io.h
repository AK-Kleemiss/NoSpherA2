#pragma once
#include "error_check.h"

/** std::getline that drops a CRLF file's trailing \r, which would otherwise spoil comparisons and the last field; every reader uses it. */
inline std::istream& getline_universal(std::istream& is, std::string& line)
{
	std::getline(is, line);
	if (!line.empty() && line.back() == '\r')
		line.pop_back();
	return is;
}

//Line readers for the text formats. At EOF getline leaves the line untouched, so a seek loop
//written as while (line.find(tag) == npos) getline(...) never ends on a truncated file
inline void read_line_or_fail(std::istream& is, std::string& line, const std::string& what, std::ostream& log)
{
	err_checkf(static_cast<bool>(getline_universal(is, line)), "File ends while reading " + what, log);
}
//Advances to the first line containing tag, which is left in line
inline void seek_line(std::istream& is, std::string& line, const std::string& tag, std::ostream& log)
{
	while (line.find(tag) == std::string::npos)
		read_line_or_fail(is, line, tag, log);
}
//Appends every number on line to out; a token that is not a number fails naming the block
template <typename T>
void append_numbers(const std::string& line, std::vector<T>& out, const std::string& what, std::ostream& log)
{
	std::istringstream is(line);
	T v;
	while (is >> v)
		out.push_back(v);
	err_checkf(is.eof(), "Not a number in " + what + ": '" + line + "'", log);
}

//Restores a stream's flags, precision and width on scope exit, including unwinding. std::fixed and
//setprecision are sticky, so every analysis entry point that formats output holds one of these.
struct ostream_format_guard
{
	std::ostream& stream;
	const std::ios_base::fmtflags flags;
	const std::streamsize precision;
	const std::streamsize width;
	explicit ostream_format_guard(std::ostream& s)
		: stream(s), flags(s.flags()), precision(s.precision()), width(s.width()) {}
	ostream_format_guard(const ostream_format_guard&) = delete;
	ostream_format_guard& operator=(const ostream_format_guard&) = delete;
	~ostream_format_guard()
	{
		stream.flags(flags);
		stream.precision(precision);
		stream.width(width);
	}
};

inline void print_centered_text(const std::string& text, int& bar_width, std::ostream& file = std::cout)
{
	const int text_length = static_cast<int>(text.length());
	const int total_padding = bar_width - text_length;
	const int padding_left = total_padding / 2;
	const int padding_right = (total_padding - padding_left) - 1;

	file << "["
		<< std::setw(padding_left) << std::setfill(' ') << ""
		<< text
		<< std::setw(padding_right) << std::setfill(' ') << ""
		<< "]" << std::endl;
}

inline void print_centered_message(const std::string& text, int bar_width, std::ostream& os = std::cout)
{
	const int text_length = static_cast<int>(text.length());
	const int total_padding = bar_width - text_length;
	const int padding_left = total_padding / 2;
	const int padding_right = (total_padding - padding_left) - 1;

	os
		<< std::setw(padding_left) << std::setfill(' ') << ""
		<< text
		<< std::setw(padding_right) << std::setfill(' ') << ""
		<< std::endl;
}
