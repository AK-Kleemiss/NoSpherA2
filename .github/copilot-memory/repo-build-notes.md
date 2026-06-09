- NoSpherA2 Windows `windows-clang-cl` can silently configure as x86 if not run in an x64 VS dev environment.
- Verify `build-windows-clang-cl/CMakeFiles/*/CMakeCXXCompiler.cmake` has `CMAKE_CXX_COMPILER_ARCHITECTURE_ID "x64"` and `CMAKE_CXX_SIZEOF_DATA_PTR "8"`.
- `Windows/Build_Dependencies/build_deps.ps1 -Platform x64` correctly imports VsDevCmd x64 and builds OCC dependencies.

- GitHub composite action inputs are stringly-typed in `if:` checks; use `fromJSON(inputs.some_bool)` instead of `if: inputs.some_bool` to avoid `'false'` evaluating truthy and unexpectedly running steps.
