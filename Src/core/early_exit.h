#pragma once
// In the in-process DLL build exit() throws NosEarlyExit instead of ending the process. Included
// last in pch.h, so only Src/ call sites are affected, not the system or OCC headers.

#ifdef NOSPHERA2_IN_PROCESS

struct NosEarlyExit
{
	int code;
};

// Not [[noreturn]] on purpose: it must return while another exception unwinds.
inline void nos_do_exit(int code)
{
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

#endif // NOSPHERA2_IN_PROCESS
