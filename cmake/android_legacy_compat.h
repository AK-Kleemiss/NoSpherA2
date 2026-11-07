// Force-included into every C++ file on Android below API 28 (the legacy APK line reaches API 21):
// bionic declares posix_spawn* and aligned_alloc only from API 28.
#pragma once
#include <errno.h>
#include <spawn.h>
#include <stdlib.h>
#if defined(__ANDROID__) && __ANDROID_API__ < 28
// occ spawns processes only for external programs (xtb, external energy models), which do not
// exist on a phone or tablet: every call reports ENOSYS and occ throws its usual spawn error.
static inline int posix_spawn(pid_t*, const char*, const posix_spawn_file_actions_t*, const posix_spawnattr_t*,
                              char* const[], char* const[]) { return ENOSYS; }
static inline int posix_spawnp(pid_t*, const char*, const posix_spawn_file_actions_t*, const posix_spawnattr_t*,
                               char* const[], char* const[]) { return ENOSYS; }
static inline int posix_spawnattr_init(posix_spawnattr_t*) { return ENOSYS; }
static inline int posix_spawnattr_destroy(posix_spawnattr_t*) { return 0; }
static inline int posix_spawnattr_setflags(posix_spawnattr_t*, short) { return ENOSYS; }
static inline int posix_spawnattr_setsigmask(posix_spawnattr_t*, const sigset_t*) { return ENOSYS; }
static inline int posix_spawn_file_actions_init(posix_spawn_file_actions_t*) { return ENOSYS; }
static inline int posix_spawn_file_actions_destroy(posix_spawn_file_actions_t*) { return 0; }
static inline int posix_spawn_file_actions_adddup2(posix_spawn_file_actions_t*, int, int) { return ENOSYS; }
static inline int posix_spawn_file_actions_addclose(posix_spawn_file_actions_t*, int) { return ENOSYS; }
// pocketfft calls ::aligned_alloc; posix_memalign (API 16) memory is released with free() as well
static inline void* aligned_alloc(size_t align, size_t size) {
  void* p = nullptr;
  return posix_memalign(&p, align < sizeof(void*) ? sizeof(void*) : align, size) ? nullptr : p;
}
#endif
