#include "../core/NoSpherA2.h"
#include "debug_utils.h"

int main(int argc, char** argv)
{
    wait_for_debugger(); // no-op unless DEBUG_WAIT env var is set
    return run_app(argc, argv);
}