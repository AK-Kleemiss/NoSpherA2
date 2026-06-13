#pragma once
// early_exit.h — intercepts exit() when NoSpherA2 logic runs in-process
// (e.g. inside the DLL used by the Tests project).
//
// When NOSPHERA2_IN_PROCESS is defined (set only in the DLL build):
//   • NosEarlyExit is a lightweight exception that carries the exit code.
//   • exit() is redefined as a function-like macro that throws instead of
//     terminating the process.  Because this header is included after all
//     system/OCC headers (at the bottom of pch.h), the standard declarations
//     are already visible and only the Src/ call-sites are affected.
//
// When building the stand-alone executable the macro is NOT defined, so all
// exit() calls behave normally.

#ifdef NOSPHERA2_IN_PROCESS

struct NosEarlyExit
{
    int code;
};

// Internal throwing function — named differently so the macro below can
// reference it without recursing.
//
// NOT declared [[noreturn]] intentionally: when called during stack unwinding
// (std::uncaught_exceptions() > 0, i.e. from a destructor while another
// exception is in flight), throwing a second exception would call
// std::terminate().  Instead we return silently and let the already-active
// exception keep propagating.
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

// Redefine exit() for all code compiled after this header.
// Function-like macros only expand when followed by '(', so identifiers such
// as exit_code or exit_fn are unaffected.
#ifdef exit
#  undef exit
#endif
#define exit(code) nos_do_exit(code)

#endif // NOSPHERA2_IN_PROCESS
