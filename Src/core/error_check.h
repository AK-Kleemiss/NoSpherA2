#pragma once
#include "pch.h"

void error_check(const bool condition, const std::source_location loc, const std::string& error_mesasge, std::ostream& log_file = std::cout);
void not_implemented(const std::source_location loc, const std::string& error_mesasge, std::ostream& log_file);
//The message is built only on failure: a literal past the SSO length cost a heap allocation per call, also in passing
//checks inside hot loops (get_atom_pos per grid point and atom)
#define err_checkf(condition, error_message, file) ((condition) ? (void)0 : error_check(false, std::source_location::current(), error_message, file))
#define err(error_message, file) error_check(false, std::source_location::current(), error_message, file)
#define err_not_impl_f(error_message, file) not_implemented(std::source_location::current(), error_message, file)
