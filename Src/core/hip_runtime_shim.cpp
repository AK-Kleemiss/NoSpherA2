//The HIP runtime, opened by name instead of linked: a fat binary has to start on a machine
//without ROCm, and linking or delay-loading amdhip64 fails before main() because clang
//registers the kernels from static initialisers. The kernel objects import the entry points
//below by their plain C names; each definition forwards to libamdhip64.so.<major> /
//amdhip64_<major>.dll once opened. All entry points resolve together and a library missing
//one counts as absent: hipErrorNoDevice, a device count of zero, a registration that does
//nothing. hip_runtime_api.h declares the same functions, so a definition whose parameters
//differ does not compile; a new runtime call in a kernel source is an unresolved symbol at
//link time and goes into NOSPHERA2_HIP_RUNTIME_ENTRIES.
#include <hip/hip_runtime_api.h>

#include <cstdio>
#include <cstdlib>
#include <string>

#if defined(_WIN32)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <dlfcn.h>
#endif

//The three the module constructors call, declared in no public header
extern "C" {
void** __hipRegisterFatBinary(const void* data);
void __hipRegisterFunction(void** modules, const void* hostFunction, char* deviceFunction,
	const char* deviceName, unsigned int threadLimit, uint3* tid, uint3* bid, dim3* blockDim,
	dim3* gridDim, int* wSize);
void __hipUnregisterFatBinary(void** modules);
}

//ret, name, parameter list, argument list; nm -u / dumpbin /SYMBOLS on a kernel object lists them
#define NOSPHERA2_HIP_RUNTIME_ENTRIES(F) \
	F(void**, __hipRegisterFatBinary, (const void* data), (data)) \
	F(void, __hipRegisterFunction, (void** modules, const void* hostFunction, char* deviceFunction, const char* deviceName, unsigned int threadLimit, uint3* tid, uint3* bid, dim3* blockDim, dim3* gridDim, int* wSize), (modules, hostFunction, deviceFunction, deviceName, threadLimit, tid, bid, blockDim, gridDim, wSize)) \
	F(void, __hipUnregisterFatBinary, (void** modules), (modules)) \
	F(hipError_t, __hipPushCallConfiguration, (dim3 gridDim, dim3 blockDim, size_t sharedMem, hipStream_t stream), (gridDim, blockDim, sharedMem, stream)) \
	F(hipError_t, __hipPopCallConfiguration, (dim3* gridDim, dim3* blockDim, size_t* sharedMem, hipStream_t* stream), (gridDim, blockDim, sharedMem, stream)) \
	F(hipError_t, hipDeviceSynchronize, (void), ()) \
	F(hipError_t, hipGetDevice, (int* deviceId), (deviceId)) \
	F(hipError_t, hipGetDeviceCount, (int* count), (count)) \
	F(hipError_t, hipGetDeviceProperties, (hipDeviceProp_t* prop, int deviceId), (prop, deviceId)) \
	F(hipError_t, hipGetLastError, (void), ()) \
	F(const char*, hipGetErrorString, (hipError_t hipError), (hipError)) \
	F(hipError_t, hipStreamCreateWithFlags, (hipStream_t* stream, unsigned int flags), (stream, flags)) \
	F(hipError_t, hipStreamDestroy, (hipStream_t stream), (stream)) \
	F(hipError_t, hipStreamSynchronize, (hipStream_t stream), (stream)) \
	F(hipError_t, hipStreamWaitEvent, (hipStream_t stream, hipEvent_t event, unsigned int flags), (stream, event, flags)) \
	F(hipError_t, hipEventCreateWithFlags, (hipEvent_t* event, unsigned flags), (event, flags)) \
	F(hipError_t, hipEventRecord, (hipEvent_t event, hipStream_t stream), (event, stream)) \
	F(hipError_t, hipEventDestroy, (hipEvent_t event), (event)) \
	F(hipError_t, hipMalloc, (void** ptr, size_t size), (ptr, size)) \
	F(hipError_t, hipHostMalloc, (void** ptr, size_t size, unsigned int flags), (ptr, size, flags)) \
	F(hipError_t, hipFree, (void* ptr), (ptr)) \
	F(hipError_t, hipHostFree, (void* ptr), (ptr)) \
	F(hipError_t, hipMemcpy, (void* dst, const void* src, size_t sizeBytes, hipMemcpyKind kind), (dst, src, sizeBytes, kind)) \
	F(hipError_t, hipMemcpyAsync, (void* dst, const void* src, size_t sizeBytes, hipMemcpyKind kind, hipStream_t stream), (dst, src, sizeBytes, kind, stream)) \
	F(hipError_t, hipMemset, (void* dst, int value, size_t sizeBytes), (dst, value, sizeBytes)) \
	F(hipError_t, hipMemGetInfo, (size_t* free, size_t* total), (free, total)) \
	F(hipError_t, hipLaunchKernel, (const void* function_address, dim3 numBlocks, dim3 dimBlocks, void** args, size_t sharedMemBytes, hipStream_t stream), (function_address, numBlocks, dimBlocks, args, sharedMemBytes, stream))

namespace {

#define NOSPHERA2_HIP_STR_(x) #x
#define NOSPHERA2_HIP_STR(x) NOSPHERA2_HIP_STR_(x)

struct hip_runtime {
#define NOSPHERA2_HIP_FIELD(ret, name, params, args) ret (*name) params = nullptr;
	NOSPHERA2_HIP_RUNTIME_ENTRIES(NOSPHERA2_HIP_FIELD)
#undef NOSPHERA2_HIP_FIELD
	bool present = false;
};

#if defined(_WIN32)
using lib_handle = HMODULE;
lib_handle open_lib(const std::string& name) { return LoadLibraryA(name.c_str()); }
void* find_sym(lib_handle h, const char* name) { return reinterpret_cast<void*>(GetProcAddress(h, name)); }
void close_lib(lib_handle h) { FreeLibrary(h); }
const char* const lib_versioned = "amdhip64_" NOSPHERA2_HIP_STR(HIP_VERSION_MAJOR) ".dll";
const char* const lib_plain = "amdhip64.dll";
const char* const lib_subdir = "/bin/";
#else
using lib_handle = void*;
lib_handle open_lib(const std::string& name) { return dlopen(name.c_str(), RTLD_NOW | RTLD_LOCAL); }
void* find_sym(lib_handle h, const char* name) { return dlsym(h, name); }
void close_lib(lib_handle h) { dlclose(h); }
const char* const lib_versioned = "libamdhip64.so." NOSPHERA2_HIP_STR(HIP_VERSION_MAJOR);
const char* const lib_plain = "libamdhip64.so";
const char* const lib_subdir = "/lib/";
#endif

//NOSPHERA2_HIP_RUNTIME names the library file outright. Otherwise the versioned name is
//tried on the loader's own search (PATH, LD_LIBRARY_PATH, the rpath, the usual places),
//then under ROCM_PATH and HIP_PATH, then /opt/rocm, then the unversioned name the same way.
lib_handle open_runtime()
{
	if (const char* file = std::getenv("NOSPHERA2_HIP_RUNTIME")) {
		if (*file) return open_lib(file);
	}
	const char* const names[] = { lib_versioned, lib_plain };
	const char* const roots[] = { std::getenv("ROCM_PATH"), std::getenv("HIP_PATH"),
#if !defined(_WIN32)
		"/opt/rocm",
#endif
		nullptr };
	for (const char* name : names) {
		if (lib_handle h = open_lib(name)) return h;
		for (const char* root : roots) {
			if (!root || !*root) continue;
			if (lib_handle h = open_lib(std::string(root) + lib_subdir + name)) return h;
		}
	}
	return nullptr;
}

//Resolved once, on the first call into any entry point, so the module constructors before
//main() and the device probes later get the same answer
const hip_runtime& runtime()
{
	static const hip_runtime rt = [] {
		hip_runtime r;
		lib_handle h = open_runtime();
		if (!h) return r;
		bool complete = true;
#define NOSPHERA2_HIP_RESOLVE(ret, name, params, args) \
		r.name = reinterpret_cast<ret (*) params>(find_sym(h, NOSPHERA2_HIP_STR(name))); \
		if (!r.name) { \
			std::fprintf(stderr, "NoSpherA2: HIP runtime found but lacks %s; treating it as absent\n", #name); \
			complete = false; \
		}
		NOSPHERA2_HIP_RUNTIME_ENTRIES(NOSPHERA2_HIP_RESOLVE)
#undef NOSPHERA2_HIP_RESOLVE
		if (!complete) {
			close_lib(h);
			return hip_runtime();
		}
		r.present = true;
		return r;
	}();
	return rt;
}

} //namespace

//The plain forwarders: an error code when the runtime is absent.
#define NOSPHERA2_HIP_FORWARD(name, params, args) \
	extern "C" hipError_t name params \
	{ \
		const hip_runtime& rt = runtime(); \
		if (!rt.present) return hipErrorNoDevice; \
		return rt.name args; \
	}

NOSPHERA2_HIP_FORWARD(__hipPushCallConfiguration, (dim3 gridDim, dim3 blockDim, size_t sharedMem, hipStream_t stream), (gridDim, blockDim, sharedMem, stream))
NOSPHERA2_HIP_FORWARD(__hipPopCallConfiguration, (dim3* gridDim, dim3* blockDim, size_t* sharedMem, hipStream_t* stream), (gridDim, blockDim, sharedMem, stream))
NOSPHERA2_HIP_FORWARD(hipDeviceSynchronize, (void), ())
NOSPHERA2_HIP_FORWARD(hipGetDevice, (int* deviceId), (deviceId))
NOSPHERA2_HIP_FORWARD(hipGetDeviceProperties, (hipDeviceProp_t* prop, int deviceId), (prop, deviceId))
NOSPHERA2_HIP_FORWARD(hipGetLastError, (void), ())
NOSPHERA2_HIP_FORWARD(hipStreamCreateWithFlags, (hipStream_t* stream, unsigned int flags), (stream, flags))
NOSPHERA2_HIP_FORWARD(hipStreamDestroy, (hipStream_t stream), (stream))
NOSPHERA2_HIP_FORWARD(hipStreamSynchronize, (hipStream_t stream), (stream))
NOSPHERA2_HIP_FORWARD(hipStreamWaitEvent, (hipStream_t stream, hipEvent_t event, unsigned int flags), (stream, event, flags))
NOSPHERA2_HIP_FORWARD(hipEventCreateWithFlags, (hipEvent_t* event, unsigned flags), (event, flags))
NOSPHERA2_HIP_FORWARD(hipEventRecord, (hipEvent_t event, hipStream_t stream), (event, stream))
NOSPHERA2_HIP_FORWARD(hipEventDestroy, (hipEvent_t event), (event))
NOSPHERA2_HIP_FORWARD(hipMalloc, (void** ptr, size_t size), (ptr, size))
NOSPHERA2_HIP_FORWARD(hipHostMalloc, (void** ptr, size_t size, unsigned int flags), (ptr, size, flags))
NOSPHERA2_HIP_FORWARD(hipFree, (void* ptr), (ptr))
NOSPHERA2_HIP_FORWARD(hipHostFree, (void* ptr), (ptr))
NOSPHERA2_HIP_FORWARD(hipMemcpy, (void* dst, const void* src, size_t sizeBytes, hipMemcpyKind kind), (dst, src, sizeBytes, kind))
NOSPHERA2_HIP_FORWARD(hipMemcpyAsync, (void* dst, const void* src, size_t sizeBytes, hipMemcpyKind kind, hipStream_t stream), (dst, src, sizeBytes, kind, stream))
NOSPHERA2_HIP_FORWARD(hipMemset, (void* dst, int value, size_t sizeBytes), (dst, value, sizeBytes))
NOSPHERA2_HIP_FORWARD(hipMemGetInfo, (size_t* free, size_t* total), (free, total))
NOSPHERA2_HIP_FORWARD(hipLaunchKernel, (const void* function_address, dim3 numBlocks, dim3 dimBlocks, void** args, size_t sharedMemBytes, hipStream_t stream), (function_address, numBlocks, dimBlocks, args, sharedMemBytes, stream))

#undef NOSPHERA2_HIP_FORWARD

//The ones whose absent answer is not an error code.
extern "C" hipError_t hipGetDeviceCount(int* count)
{
	const hip_runtime& rt = runtime();
	if (!rt.present) {
		if (count) *count = 0;
		return hipErrorNoDevice;
	}
	return rt.hipGetDeviceCount(count);
}

extern "C" const char* hipGetErrorString(hipError_t hipError)
{
	const hip_runtime& rt = runtime();
	if (!rt.present) return "HIP runtime (amdhip64) not found";
	return rt.hipGetErrorString(hipError);
}

//Registration without a runtime registers nothing, and the null module handle is ignored
extern "C" void** __hipRegisterFatBinary(const void* data)
{
	const hip_runtime& rt = runtime();
	if (!rt.present) return nullptr;
	return rt.__hipRegisterFatBinary(data);
}

extern "C" void __hipRegisterFunction(void** modules, const void* hostFunction, char* deviceFunction,
	const char* deviceName, unsigned int threadLimit, uint3* tid, uint3* bid, dim3* blockDim,
	dim3* gridDim, int* wSize)
{
	const hip_runtime& rt = runtime();
	if (!rt.present || !modules) return;
	rt.__hipRegisterFunction(modules, hostFunction, deviceFunction, deviceName, threadLimit, tid, bid,
		blockDim, gridDim, wSize);
}

extern "C" void __hipUnregisterFatBinary(void** modules)
{
	const hip_runtime& rt = runtime();
	if (!rt.present || !modules) return;
	rt.__hipUnregisterFatBinary(modules);
}
