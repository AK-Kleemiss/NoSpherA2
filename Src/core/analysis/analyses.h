#pragma once
#include <ostream>
#include <vector>

struct options;
//an analysis that reads the wavefunction, prints its tables and ends the run; run_app_impl runs the first
//one asked for and refuse_unread_bonding_options names it, so a new analysis is a row here, its options
//in options_parse.cpp and its help items in help_text.cpp
struct analysis
{
	const char* (*asked)(const options& opt);  //the option that asked for it, nullptr if none did
	bool console;  //tables to the terminal: cout is handed back from the log before run
	int (*run)(options& opt, const char* argv0, std::ostream& log);
};
const std::vector<analysis>& analyses();
