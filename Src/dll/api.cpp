#include "pch.h"
#include "api.h"

int run_app(int argc, char** argv);

extern "C" NOS_API int NOS_CALLCONV nosphera2_run(int argc, char** argv)
{
    return run_app(argc, argv);
}
