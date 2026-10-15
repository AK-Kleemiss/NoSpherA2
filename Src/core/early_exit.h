#pragma once
// While nos_in_process is set (the DLL sets it) exit() throws NosEarlyExit, which run_app turns into its
// return code, instead of ending the host process. The exe leaves it unset and exits as before. Included
// last in pch.h, so only Src/ call sites are affected, not the system or OCC headers. A throw out of an
// OpenMP region still ends the process.
#include <cstdlib>
#include <exception>

struct NosEarlyExit
{
	int code;
};

inline bool nos_in_process = false;

// Not [[noreturn]] on purpose: it must return while another exception unwinds.
inline void nos_do_exit(int code)
{
	if (!nos_in_process)
		std::exit(code);
	if (std::uncaught_exceptions() > 0)
	{
		// A destructor is calling exit() while an exception is being propagated.
		// We cannot throw — just return so the destructor finishes cleanly and
		// the original exception continues to unwind.
		return;
	}
	throw NosEarlyExit{ code };
}

// Function-like, so identifiers such as exit_code are unaffected.
#ifdef exit
#  undef exit
#endif
#define exit(code) nos_do_exit(code)
