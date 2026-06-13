/* Returns 0 if the build host supports AVX (including OS XSAVE support), 1 otherwise. */
#ifdef _MSC_VER
#include <intrin.h>
int main(void) {
  int info[4];
  __cpuid(info, 1);
  /* CPUID.1:ECX bit 27 = OSXSAVE, bit 28 = AVX */
  if (!(info[2] & (1 << 27)) || !(info[2] & (1 << 28)))
    return 1;
  /* XCR0 bits 1|2: XMM and YMM state enabled by the OS */
  return ((_xgetbv(0) & 6) == 6) ? 0 : 1;
}
#else
int main(void) { return __builtin_cpu_supports("avx") ? 0 : 1; }
#endif
