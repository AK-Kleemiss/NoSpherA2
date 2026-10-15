#NOSPHERA2_DEPENDENCIES_ONLY: build every dependency target and nothing of NoSpherA2 itself.
#CI builds this tree first and caches it, so the project build that follows only compiles
#NoSpherA2's own sources (.github/workflows/c-cpp_all.yml).
enable_googletests()

add_custom_target(NoSpherA2Dependencies ALL
    DEPENDS
        cint
        tbb
        occ_main
        gtest
        BasisSetConverter
        cargo-build-featomic
)
