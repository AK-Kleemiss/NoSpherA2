#pragma once

#ifdef _WIN32
#ifdef BUILDING_DLL
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT __declspec(dllimport)
#endif
#else
#define DLL_EXPORT
#endif

extern "C" {
    DLL_EXPORT int __cdecl nosphera2_run(int argc, char** argv);
}