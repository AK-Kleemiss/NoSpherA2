#pragma once

#include <windows.h>
#include <delayimp.h>
#include <iostream>

FARPROC WINAPI DelayLoadFailureHook(unsigned dliNotify, PDelayLoadInfo pdli);