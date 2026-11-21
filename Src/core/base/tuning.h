#pragma once
#include <string>

//Developer knobs (diagnostics, cut-offs, A/B switches), set with -tune NAME[=VALUE]; a bare NAME is "1".
//They were environment variables; an option keeps a variable left in a shell from changing a result unseen.
//The knob's value, or nullptr when it was not given
const char *tuning(const char *name);
//nullptr removes the knob
void set_tuning(const std::string &name, const char *value);
//-tune is a per-run flag, but the knobs are per process, and Olex2 and the in-process tests call run_app repeatedly.
//Puts back what this run's -tune replaced; set_tuning from code (tests/src/main.cpp) stays.
void tuning_end_run();
