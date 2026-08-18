#include "../core/NoSpherA2.h"
#include "debug_utils.h"
#include "exit_diagnostics.h"

int main(int argc, char** argv)
{
    install_exit_diagnostics(); // records how the process leaves
    wait_for_debugger(); // no-op unless DEBUG_WAIT env var is set
    return run_app(argc, argv);
}