#Run as a script at every build (target NoSpherA2GitCommit in Src/core): writes the commit the
#binary is built from into git_commit.h. configure_file leaves the header untouched while the text
#is unchanged, so a build on the same commit recompiles nothing. "-dirty" means tracked changes
#outside tests/ (whose goldens the test runs rewrite); a tree without git gives "unknown".
set(_c "unknown")
if(GIT_EXECUTABLE)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${SRC}" rev-parse --short=8 HEAD
        OUTPUT_VARIABLE _h OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE _r ERROR_QUIET)
    if(_r EQUAL 0)
        set(_c "${_h}")
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${SRC}" diff --quiet HEAD -- . ":(exclude)tests"
            RESULT_VARIABLE _d ERROR_QUIET)
        if(NOT _d EQUAL 0)
            string(APPEND _c "-dirty")
        endif()
    endif()
endif()
file(CONFIGURE OUTPUT "${OUT}" CONTENT "#pragma once\n#define NOSPHERA2_GIT_COMMIT \"${_c}\"\n" @ONLY)
