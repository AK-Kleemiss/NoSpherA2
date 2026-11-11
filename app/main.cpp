// NoSpherA2 lives in the shared library next to this executable (NoSpherA2.dll, libNoSpherA2.so,
// libNoSpherA2.dylib), which Olex2 also loads in-process; the exe is only its launcher, so the code is
// shipped once. Loaded at run time instead of linked: app/isa_check.c then refuses an unsupported CPU,
// with its message, before any library code or static initialiser runs. Android links it all statically.
#include "debug_utils.h"
#include <cstdio>
#ifdef NOS_STATIC_APP
int run_app(int argc, char** argv);
#elif defined(_WIN32)
#include <string>
#else
#include <dlfcn.h>
#include <string>
#include <unistd.h>
#endif

int main(int argc, char** argv)
{
	wait_for_debugger(); // no-op unless DEBUG_WAIT env var is set
#ifdef NOS_STATIC_APP
	return run_app(argc, argv);
#else
	using nos_main_t = int (*)(int, char**);
	nos_main_t run = nullptr;
#ifdef _WIN32
	std::wstring lib(32768, L'\0');
	lib.resize(GetModuleFileNameW(nullptr, lib.data(), static_cast<DWORD>(lib.size())));
	lib = lib.substr(0, lib.find_last_of(L"\\/") + 1) + L"NoSpherA2.dll";
	// tbb12.dll and libiomp5md.dll are found in the exe's directory, which is the library's
	if (HMODULE h = LoadLibraryW(lib.c_str()))
		run = reinterpret_cast<nos_main_t>(GetProcAddress(h, "nosphera2_main"));
	if (!run)
	{
		fprintf(stderr, "NoSpherA2: cannot load %ls (Windows error %lu)\n", lib.c_str(), GetLastError());
		return 1;
	}
#else
#ifdef __APPLE__
	std::string lib = "@executable_path/libNoSpherA2.dylib";
#else
	std::string lib(4096, '\0');
	const ssize_t n = readlink("/proc/self/exe", lib.data(), lib.size());
	lib.resize(n > 0 ? n : 0);
	lib = lib.substr(0, lib.find_last_of('/') + 1) + "libNoSpherA2.so";
#endif
	if (void* h = dlopen(lib.c_str(), RTLD_NOW | RTLD_LOCAL))
		run = reinterpret_cast<nos_main_t>(dlsym(h, "nosphera2_main"));
	if (!run)
	{
		fprintf(stderr, "NoSpherA2: cannot load %s: %s\n", lib.c_str(), dlerror());
		return 1;
	}
#endif
	return run(argc, argv);
#endif
}
