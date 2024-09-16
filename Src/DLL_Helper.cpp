#include "DLL_Helper.h"

FARPROC WINAPI DelayLoadFailureHook(unsigned dliNotify, PDelayLoadInfo pdli) {
    if (dliNotify == dliFailLoadLib) {
        std::cerr << "Failed to load DLL: " << pdli->szDll << std::endl;
    }
    else if (dliNotify == dliFailGetProc) {
        std::cerr << "Failed to get procedure: " << pdli->dlp.szProcName << std::endl;
    }
    return nullptr;
}