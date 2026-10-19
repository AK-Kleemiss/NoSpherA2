- NoSpherA2 Windows presets can silently configure as x86 if not run in an x64 VS dev environment.
- Verify `build/release-windows/CMakeFiles/*/CMakeCXXCompiler.cmake` has `CMAKE_CXX_COMPILER_ARCHITECTURE_ID "x64"` and `CMAKE_CXX_SIZEOF_DATA_PTR "8"`.
- OCC builds oneTBB as a shared runtime while keeping OCC libraries static. Package/copy `tbb12.dll`, `libtbb.so*`, or `libtbb.dylib` beside NoSpherA2.
- Do not restore `tbbmalloc_proxy` or `tbbmalloc` linking for NoSpherA2. The allocator-proxy workaround was removed after the OCC DF teardown fix.
- Windows validation: configure with `-DNOSPHERA2_BUILD_TESTS=ON` and run `ctest --preset release-windows` (the Visual Studio solution was removed in October 2026).

- GitHub composite action inputs are stringly-typed in `if:` checks; use `fromJSON(inputs.some_bool)` instead of `if: inputs.some_bool` to avoid `'false'` evaluating truthy and unexpectedly running steps.
