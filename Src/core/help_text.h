#pragma once
#include <string>
#include <vector>

//one block of the -h text: an option (flag, args, text in the description column), or prose when flag is empty
struct help_item
{
	std::string topic, sub, flag, args, text;
};
std::vector<help_item>& help_items();
