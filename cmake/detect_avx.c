/* Returns 0 if the build host supports AVX (including OS XSAVE support), 1 otherwise.
   Built with NOS_DETECT_AVX2 defined it asks for AVX2 and FMA on top. */
#ifdef _MSC_VER
#include <intrin.h>
int main(void) {
  int info[4];
  __cpuid(info, 1);
  /* CPUID.1:ECX bit 27 = OSXSAVE, bit 28 = AVX */
  if (!(info[2] & (1 << 27)) || !(info[2] & (1 << 28)))
    return 1;
  /* XCR0 bits 1|2: XMM and YMM state enabled by the OS */
  if ((_xgetbv(0) & 6) != 6)
    return 1;
#ifdef NOS_DETECT_AVX2
  /* CPUID.1:ECX bit 12 = FMA, CPUID.(7,0):EBX bit 5 = AVX2 */
  if (!(info[2] & (1 << 12)))
    return 1;
  __cpuid(info, 0);
  if (info[0] < 7)
    return 1;
  __cpuidex(info, 7, 0);
  if (!(info[1] & (1 << 5)))
    return 1;
#endif
  return 0;
}
#else
int main(void) {
#ifdef NOS_DETECT_AVX2
  return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma") ? 0 : 1;
#else
  return __builtin_cpu_supports("avx") ? 0 : 1;
#endif
}
#endif
