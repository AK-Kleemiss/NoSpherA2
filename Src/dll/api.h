#pragma once

#if defined(_WIN32)
    #ifdef BUILDING_DLL
        #define NOS_API __declspec(dllexport)
    #else
        #define NOS_API __declspec(dllimport)
    #endif
    #define NOS_CALLCONV __cdecl
#else
    // Export symbols on ELF platforms.
    #if defined(__GNUC__)
        #define NOS_API __attribute__((visibility("default")))
    #else
        #define NOS_API
    #endif
    #define NOS_CALLCONV
#endif

extern "C" {
    NOS_API int NOS_CALLCONV nosphera2_run(int argc, char** argv);
}