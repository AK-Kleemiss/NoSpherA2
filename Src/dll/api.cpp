#include "pch.h"
#include "api.h"
#include "core/early_exit.h"

int run_app(int argc, char** argv);

// set by app/isa_check.c when this CPU lacks the instruction set of the build. Windows fails the load
// instead; ponytail: a static initialiser of the .so that already uses AVX still SIGILLs before this
extern "C" { int nos_isa_refused = 0; }

extern "C" NOS_API int NOS_CALLCONV nosphera2_run(int argc, char** argv)
{
	if (nos_isa_refused)
		return 86;
	nos_in_process = true;
	return run_app(argc, argv);
}

// The NoSpherA2 executable (app/main.cpp) is only a launcher for this: exit() still ends the process,
// exactly as it did when the exe carried the code itself
extern "C" NOS_API int NOS_CALLCONV nosphera2_main(int argc, char** argv)
{
	if (nos_isa_refused)
		return 86;
	return run_app(argc, argv);
}
