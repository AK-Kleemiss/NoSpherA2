#include "pch.h"
#include "tuning.h"

static std::map<std::string, std::string> &tuning_knobs()
{
	static std::map<std::string, std::string> knobs;
	return knobs;
}
const char *tuning(const char *name)
{
	const auto &knobs = tuning_knobs();
	const auto it = knobs.find(name);
	return it == knobs.end() ? nullptr : it->second.c_str();
}
void set_tuning(const std::string &name, const char *value)
{
	if (value)
		tuning_knobs()[name] = value;
	else
		tuning_knobs().erase(name);
}
