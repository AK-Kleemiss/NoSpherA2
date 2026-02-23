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
    #include <fstream>
    #include <string>
    #include <csignal>
#endif

inline volatile bool g_debug_attached = false;

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

    while (!IsDebuggerPresent()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }

    std::cout << "Debugger detected! Resuming..." << std::endl;
    DebugBreak();

#else
    // --- LINUX IMPLEMENTATION ---
    std::cout << "Process ID: " << getpid() << std::endl;
    std::cout << "Action: Attach debugger in CLion." << std::endl;
    std::cout << "Status: Waiting..." << std::endl;

    // Linux equivalent of IsDebuggerPresent()
    auto is_debugger_present = []() -> bool {
        std::ifstream status_file("/proc/self/status");
        std::string line;
        while (std::getline(status_file, line)) {
            if (line.find("TracerPid:") == 0) {
                // Find the first number in the string and convert it
                size_t pos = line.find_first_of("0123456789");
                if (pos != std::string::npos) {
                    return std::stoi(line.substr(pos)) > 0;
                }
            }
        }
        return false;
    };

    // Loop until manual flag is set OR the OS detects the debugger
    while (!g_debug_attached && !is_debugger_present()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }

    std::cout << "Debugger detected! Resuming..." << std::endl;

    // std::raise(SIGTRAP) is the Linux equivalent of Windows DebugBreak()
    // It will automatically pause CLion right here so you can start stepping.
    std::raise(SIGTRAP);
#endif
}
#endif //NOSPHERA2_DEBUG_UTILS_H