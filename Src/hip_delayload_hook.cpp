//A Windows build with AMD kernels delay-loads the HIP runtime (NOSPHERA2_HIP_DELAYLOAD in
//CMakeLists.txt, DelayLoadDLLs in NoSpherA2_gpu.props) so that it can start on a machine
//without ROCm. That alone is not enough: clang registers every kernel object's code with the
//runtime from a static initialiser, so __hipRegisterFatBinary is called before main() and the
//loader raises 0xC06D007E on the missing DLL before a single line of ours has run. This hook
//answers instead. When amdhip64 cannot be loaded, every entry the delay-load helper asks for
//resolves to a stub that reports hipErrorNoDevice, and the registration calls come and go
//without a runtime. Nothing else ever reaches the stubs: each GPU path begins by counting
//devices, and gpu_backend.h answers that call from the presence check, without the runtime,
//when the DLL is absent.
//
//This is compiled into each executable and the DLL rather than into the core library, since
//the linker takes the first definition of __pfnDliFailureHook2 it meets, and an object in a
//static library is only pulled in when something already needs it. It is harmless in a build
//that delay-loads nothing: the hook is only consulted on a failed delay-load, and only reacts
//to the HIP runtime's DLL name.
#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <delayimp.h>
#include <cstring>

namespace {
	//hipErrorNoDevice as a plain number: this file is built into targets that see no HIP
	//headers. All delay-loaded HIP entries return a hipError_t or a handle nobody dereferences.
	int absent_hip_runtime_stub() { return 100; }

	FARPROC WINAPI hip_delay_load_failure(unsigned notification, PDelayLoadInfo info)
	{
		if (_strnicmp(info->szDll, "amdhip64", 8) != 0) return nullptr;
		if (notification == dliFailLoadLib)
			//Any module will do as the stand-in: the helper only uses it for GetProcAddress,
			//which fails and brings the request back here as dliFailGetProc.
			return reinterpret_cast<FARPROC>(GetModuleHandleW(nullptr));
		if (notification == dliFailGetProc)
			return reinterpret_cast<FARPROC>(&absent_hip_runtime_stub);
		return nullptr;
	}
}

extern "C" const PfnDliHook __pfnDliFailureHook2 = hip_delay_load_failure;
#endif
