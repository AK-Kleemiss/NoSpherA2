/* Start-up check that this CPU runs the instruction set NoSpherA2 was built for (NOS_BUILT_ISA from
   nosphera2_isa_check in cmake/NoSpherA2Avx.cmake: 1 SSE4.2, 2 AVX, 3 AVX2 + FMA). Without it an AVX2 build
   on an older CPU dies on its first AVX instruction, often inside a static initialiser, with no output at
   all. It refuses to run, it never picks another code path.
   Exe: exit 86 with the message on stderr and in NoSpherA2.log, which Olex2 (utilities.nos_failure) turns
   into an error for the user.
   Library (NOS_ISA_DLL): it must never end its host (Olex2). On Windows the nonzero C initialiser makes the
   CRT's DllMain return FALSE, so LoadLibrary fails (ERROR_DLL_INIT_FAILED) and the caller falls back to the
   exe, which says why. A shared object cannot fail its dlopen, so nos_isa_refused makes nosphera2_run
   return 86 instead.
   Compiled without SSE3+ and outside LTO, and run before the C++ initialisers. */
#if (defined(_M_X64) || defined(__x86_64__)) && defined(NOS_BUILT_ISA)
#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#include <intrin.h>
#endif

#if NOS_BUILT_ISA >= 3
#define NOS_ISA_NAME "AVX2 + FMA"
#elif NOS_BUILT_ISA == 2
#define NOS_ISA_NAME "AVX"
#else
#define NOS_ISA_NAME "SSE4.2"
#endif
#if NOS_BUILT_ISA >= 2
#define NOS_ISA_HINT "Use the SSE4.2 build of NoSpherA2 that comes with Olex2, or build NoSpherA2 on this machine.\n"
#else
#define NOS_ISA_HINT "NoSpherA2 needs an x86-64 CPU with SSE4.2 (2008 or later).\n"
#endif

static int nos_cpu_ok(void)
{
#ifdef _MSC_VER
	int info[4];
	__cpuid(info, 1);
	/* ECX bit 20 SSE4.2, 27 OSXSAVE, 28 AVX, 12 FMA; XCR0 bits 1|2: the OS saves XMM and YMM */
	if (!(info[2] & (1 << 20)))
		return 0;
	if (NOS_BUILT_ISA >= 2 && (!(info[2] & (1 << 27)) || !(info[2] & (1 << 28)) || (_xgetbv(0) & 6) != 6))
		return 0;
	if (NOS_BUILT_ISA >= 3)
	{
		if (!(info[2] & (1 << 12)))
			return 0;
		__cpuid(info, 0);
		if (info[0] < 7)
			return 0;
		__cpuidex(info, 7, 0);
		if (!(info[1] & (1 << 5)))
			return 0;
	}
	return 1;
#else
	__builtin_cpu_init(); /* needed before main; it includes the OS XSAVE check */
	return __builtin_cpu_supports("sse4.2")
		&& (NOS_BUILT_ISA < 2 || __builtin_cpu_supports("avx"))
		&& (NOS_BUILT_ISA < 3 || (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")));
#endif
}

#ifdef NOS_ISA_DLL
extern int nos_isa_refused; /* Src/dll/api.cpp */
#endif

static int nos_isa_check(void)
{
	static const char msg[] =
		"NoSpherA2 was built for " NOS_ISA_NAME " instructions, which this CPU (or its operating system) does not support.\n"
		NOS_ISA_HINT
		"If that does not help, send this message to florian@olexsys.org.\n";
	if (nos_cpu_ok())
		return 0;
	fputs(msg, stderr);
#ifdef NOS_ISA_DLL
	nos_isa_refused = 1;
	return 1;
#else
	{
		FILE *log = fopen("NoSpherA2.log", "w");
		if (log)
		{
			fputs(msg, log);
			fclose(log);
		}
	}
	exit(86);
#endif
}

#ifdef _MSC_VER
/* .CRT$XIU: the C initialisers, which the CRT runs before every C++ one (.CRT$XC*) */
#pragma section(".CRT$XIU", long, read)
__declspec(allocate(".CRT$XIU")) int (*nos_isa_check_entry)(void) = nos_isa_check;
#else
/* priority 101 comes before every default-priority C++ initialiser of the binary */
__attribute__((constructor(101))) static void nos_isa_check_entry(void) { nos_isa_check(); }
#endif
#endif
