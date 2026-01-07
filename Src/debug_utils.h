//
// Created by Lucas on 09/12/25.
//
#ifndef NOSPHERA2_DEBUG_UTILS_H
#define NOSPHERA2_DEBUG_UTILS_H
#include <iostream>
#include <cstdlib>

#ifdef _WIN32
    #include <windows.h>
    #include <thread>
    #include <chrono>
#else
    #include <unistd.h>
    #include <thread>
    #include <chrono>
#endif

volatile bool g_debug_attached = false;
inline void wait_for_debugger() {
    if (!std::getenv("DEBUG_WAIT")) {
        return;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "   WAITING FOR DEBUGGER TO ATTACH" << std::endl;
    std::cout << "========================================" << std::endl;

#ifdef _WIN32
    // --- WINDOWS IMPLEMENTATION ---
    std::cout << "Process ID: " << GetCurrentProcessId() << std::endl;
    std::cout << "Action: Attach debugger in your IDE." << std::endl;
    std::cout << "Status: Waiting..." << std::endl;

    // Windows has a built-in API to check if a debugger is attached.
    // We loop until it returns true.
    while (!IsDebuggerPresent()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }

    std::cout << "Debugger detected! Resuming..." << std::endl;
    // Optional: Trigger a breakpoint immediately so you land in the code
    DebugBreak();

#else
    std::cout << "Process ID: " << getpid() << std::endl;
    std::cout << "Action: " << std::endl;
    std::cout << "  1. Open your IDE (CLion/VS Code) and 'Attach to Process'" << std::endl;
    std::cout << "  2. Run command 'expr attached = true'" << std::endl;


    while (!g_debug_attached) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }

    std::cout << "Variable set! Resuming..." << std::endl;
#endif
}
#endif //NOSPHERA2_DEBUG_UTILS_H