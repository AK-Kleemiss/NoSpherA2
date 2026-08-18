//
// Records how the process leaves.
//
// Runs that end with exit code -1 leave no trace anywhere: nothing in the
// program log, nothing on stdout or stderr, and no Windows Error Reporting
// entry. Each hook below writes a symbolised stack to NoSpherA2_exit_trace.txt
// the moment it fires, and *which* hook fires is itself most of the answer:
//
//   exit() called          - somebody in the program asked to leave, and the
//                            stack names them
//   unhandled SEH          - a fault; the code says which kind
//   std::terminate         - an exception escaped, or a noexcept was violated
//   CRT invalid parameter  - the CRT was handed something illegal and would
//                            otherwise die silently in a release build
//   pure virtual call      - a call through a half-destroyed object
//
// If the file stays empty the process did not leave on its own, which points
// outside the program entirely.
//
#ifndef NOSPHERA2_EXIT_DIAGNOSTICS_H
#define NOSPHERA2_EXIT_DIAGNOSTICS_H

#include <cstdio>
#include <cstdlib>
#include <exception>

#ifdef _WIN32
#include <windows.h>
#include <dbghelp.h>
#pragma comment(lib, "dbghelp.lib")

inline void n2_write_trace(const char *why, const char *detail)
{
    FILE *f = nullptr;
    if (fopen_s(&f, "NoSpherA2_exit_trace.txt", "a") != 0 || f == nullptr)
        return;
    fprintf(f, "\n=== %s ===\n", why);
    if (detail != nullptr)
        fprintf(f, "%s\n", detail);

    void *frames[62];
    const USHORT n = CaptureStackBackTrace(0, 62, frames, nullptr);
    const HANDLE proc = GetCurrentProcess();
    static bool symbols_ready = false;
    if (!symbols_ready)
    {
        SymSetOptions(SYMOPT_DEFERRED_LOADS | SYMOPT_LOAD_LINES | SYMOPT_UNDNAME);
        SymInitialize(proc, nullptr, TRUE);
        symbols_ready = true;
    }

    char storage[sizeof(SYMBOL_INFO) + 512];
    SYMBOL_INFO *si = reinterpret_cast<SYMBOL_INFO *>(storage);
    si->SizeOfStruct = sizeof(SYMBOL_INFO);
    si->MaxNameLen = 500;
    for (USHORT i = 0; i < n; ++i)
    {
        const DWORD64 addr = reinterpret_cast<DWORD64>(frames[i]);
        DWORD64 disp = 0;
        if (SymFromAddr(proc, addr, &disp, si))
        {
            IMAGEHLP_LINE64 line;
            DWORD line_disp = 0;
            line.SizeOfStruct = sizeof(IMAGEHLP_LINE64);
            if (SymGetLineFromAddr64(proc, addr, &line_disp, &line))
                fprintf(f, "  %2u  %s + 0x%llx   (%s:%lu)\n", i, si->Name,
                        static_cast<unsigned long long>(disp), line.FileName, line.LineNumber);
            else
                fprintf(f, "  %2u  %s + 0x%llx\n", i, si->Name,
                        static_cast<unsigned long long>(disp));
        }
        else
        {
            fprintf(f, "  %2u  0x%016llx\n", i, static_cast<unsigned long long>(addr));
        }
    }
    fflush(f);
    fclose(f);
}

inline void n2_on_exit_called()
{
    // exit() does not unwind, so this runs on top of the caller's frames and
    // the stack still names whoever asked to leave.
    n2_write_trace("exit() called", nullptr);
}

inline LONG WINAPI n2_on_seh(EXCEPTION_POINTERS *ep)
{
    char detail[160];
    sprintf_s(detail, "code 0x%08lX at 0x%p",
              ep->ExceptionRecord->ExceptionCode, ep->ExceptionRecord->ExceptionAddress);
    n2_write_trace("unhandled SEH exception", detail);
    return EXCEPTION_CONTINUE_SEARCH;
}

inline void n2_on_terminate()
{
    const char *what = "no active exception";
    try
    {
        if (std::current_exception())
            std::rethrow_exception(std::current_exception());
    }
    catch (const std::exception &e) { what = e.what(); }
    catch (...) { what = "non-standard exception"; }
    n2_write_trace("std::terminate", what);
    abort();
}

inline void n2_on_invalid_parameter(const wchar_t *, const wchar_t *, const wchar_t *,
                                    unsigned int line, uintptr_t)
{
    char detail[128];
    sprintf_s(detail, "CRT rejected a parameter, source line %u", line);
    n2_write_trace("CRT invalid parameter", detail);
}

inline void n2_on_purecall()
{
    n2_write_trace("pure virtual call", nullptr);
}

inline void install_exit_diagnostics()
{
    atexit(n2_on_exit_called);
    SetUnhandledExceptionFilter(n2_on_seh);
    std::set_terminate(n2_on_terminate);
    _set_invalid_parameter_handler(n2_on_invalid_parameter);
    _set_purecall_handler(n2_on_purecall);
}
#else
inline void install_exit_diagnostics() {}
#endif

#endif // NOSPHERA2_EXIT_DIAGNOSTICS_H
