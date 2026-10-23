#include "pch.h"
#include "api.h"
#include "core/early_exit.h"

int run_app(int argc, char** argv);

extern "C" NOS_API int NOS_CALLCONV nosphera2_run(int argc, char** argv)
{
	nos_in_process = true;
	return run_app(argc, argv);
}
