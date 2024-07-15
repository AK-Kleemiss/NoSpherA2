#ifndef DLL_HELPER_H
#define DLL_HELPER_H

#include <windows.h>
#include <delayimp.h>
#include <iostream>

FARPROC WINAPI DelayLoadFailureHook(unsigned dliNotify, PDelayLoadInfo pdli);

//extern "C" {
//	PfnDliHook __pfnDliFailureHook2;
//}

#endif // DELAY_LOAD_HELPER_H