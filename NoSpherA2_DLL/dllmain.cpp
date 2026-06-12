// dllmain.cpp : Defines the entry point for the DLL application.
#include "../Src/pch.h"
#include "pch_dll.h"
#include <direct.h>             // _getcwd, _chdir
#include <occ/core/data_directory.h>  // occ::set_data_directory

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
        break;
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

// main() from NoSpherA2.cpp is compiled into this DLL via the Src/*.cpp glob.
// On MSVC a DLL has no CRT startup guard around main(), so calling it as a
// regular function is safe.
extern int main(int argc, char** argv);

// ---------------------------------------------------------------------------
// MSVC enforces two incompatible rules that prevent mixing __try/__except with
// C++ exception handling in the same function:
//
//   C2712 – __try cannot appear in a function that has any object with a
//            non-trivial destructor in scope (even std::streambuf* triggers
//            this because cout.rdbuf() involves C++ machinery).
//   C2713 – Only one form of exception handling (C++ -or- SEH) per function.
//
// The standard workaround is to split into two helpers:
//   nos_cpp_dispatch  – handles C++ exceptions (try/catch only, no __try)
//   nos_seh_dispatch  – handles native SEH failures (__try only, no C++ try;
//                       all locals must be trivially-destructible)
// ---------------------------------------------------------------------------

// C++ exception layer.  Calls main() and translates C++ exceptions to int
// return codes.  savedCout/savedCerr are restored inside the catch blocks so
// that diagnostic messages reach the test runner rather than the (already
// destroyed) log_file buffer that main() may have installed.
static int nos_cpp_dispatch(int argc, char** argv,
                             std::streambuf* savedCout,
                             std::streambuf* savedCerr)
{
    try {
        return main(argc, argv);
    }
    catch (const NosEarlyExit& e) {
        // exit() was called inside NoSpherA2 logic — treat it as a normal return.
        return e.code;
    }
    catch (const std::exception& e) {
        std::cout.rdbuf(savedCout);
        std::cerr.rdbuf(savedCerr);
        std::cerr << "[NoSpherA2_run] Unhandled exception: " << e.what() << "\n";
        return -1;
    }
    catch (...) {
        std::cout.rdbuf(savedCout);
        std::cerr.rdbuf(savedCerr);
        std::cerr << "[NoSpherA2_run] Unhandled unknown exception\n";
        return -2;
    }
}

// SEH layer.  This remains as a narrow last-resort guard for native cleanup
// failures so one broken in-process run does not necessarily kill the whole
// Visual Studio test host.
//
// Requirements for __try/__except:
//   • No C++ try/catch in the same function (C2713).
//   • No locals with non-trivial destructors in scope (C2712).
//   All locals here are raw pointers or int — trivially destructible.
static int nos_seh_dispatch(int argc, char** argv,
                             std::streambuf** pSavedCout,
                             std::streambuf** pSavedCerr)
{
    int result = 0;
    __try {
        result = nos_cpp_dispatch(argc, argv, *pSavedCout, *pSavedCerr);
    }
    __except (EXCEPTION_EXECUTE_HANDLER)
    {
        // Restore streams (they may still be redirected to the now-unwound
        // log_file buffer) before writing the diagnostic.  Always catch native
        // SEH exceptions here so that cout is never left pointing at a
        // destroyed log_file buffer, which would corrupt subsequent tests.
        std::cout.rdbuf(*pSavedCout);
        std::cerr.rdbuf(*pSavedCerr);
        std::cerr << "[NoSpherA2_run] Native exception caught: 0x"
                  << std::hex << GetExceptionCode() << std::dec << "\n";
        result = static_cast<int>(GetExceptionCode());
    }
    return result;
}

// ---------------------------------------------------------------------------
// NoSpherA2_run — exported in-process entry point used by the Tests project.
//
// Saves and restores the process working directory around the call so that
// each test run is isolated to its own temp directory.  Because execution
// stays inside the test host process the VS debugger is already attached:
// just set breakpoints in Src/ and press F5 — no manual attach needed.
//
// Parameters:
//   workDir   – absolute path that becomes the CWD for this run
//   argc/argv – full argument vector (argv[0] is the fake exe name)
// ---------------------------------------------------------------------------
extern "C" DLL_EXPORT int __cdecl NoSpherA2_run(
    const char* workDir,
    int         argc,
    char**      argv)
{
    char savedDir[MAX_PATH] = {};
    _getcwd(savedDir, MAX_PATH);
    _chdir(workDir);

    // The test runner sets OCC_DATA_PATH via SetEnvironmentVariableA (Win32 API).
    // NoSpherA2's ensure_occ_data_path() reads it with _dupenv_s (CRT function).
    // On Windows these two can diverge: SetEnvironmentVariableA updates the
    // Win32 environment block while the CRT maintains its own cached table.
    // If they diverge, _dupenv_s sees stale data and the fallback fires,
    // writing "OCC_DATA_PATH not set or invalid. Using: ..." into the log file.
    //
    // Fix: read OCC_DATA_PATH via GetEnvironmentVariableA (Win32), then sync
    // the value into the CRT table with _putenv_s AND into OCC's own override
    // with occ::set_data_directory().  This ensures all three consumers agree.
    {
        char occBuf[MAX_PATH] = {};
        DWORD occLen = GetEnvironmentVariableA("OCC_DATA_PATH", occBuf, MAX_PATH);
        if (occLen > 0 && occLen < MAX_PATH) {
            occ::set_data_directory(occBuf);          // OCC internal override
            _putenv_s("OCC_DATA_PATH", occBuf);       // sync CRT env table
        }
    }

    // When exit() throws NosEarlyExit, the stack unwinds through main() and
    // destroys log_file — but main() never reaches its explicit
    // std::cout.rdbuf(_coutbuf) restore line.  Save and restore it here so
    // the test host's stdout is always valid after we return.
    std::streambuf* savedCout = std::cout.rdbuf();
    std::streambuf* savedCerr = std::cerr.rdbuf();

    int result = nos_seh_dispatch(argc, argv, &savedCout, &savedCerr);

    // Always restore cout/cerr regardless of how main() exited.
    std::cout.rdbuf(savedCout);
    std::cerr.rdbuf(savedCerr);

    _chdir(savedDir);
    return result;
}
