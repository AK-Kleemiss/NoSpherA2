# Tests using tests/sucrose_fchk_SF
set_tests_properties(
    TomlIntegrationTests.Fractal
    TomlIntegrationTests.Properties
    TomlIntegrationTests.SucroseSF
    TomlIntegrationTests.SucroseTwin
    PROPERTIES
        RESOURCE_LOCK integration_sucrose_fchk_SF
)

# Tests using tests/sucrose_IAM_SF
set_tests_properties(
    TomlIntegrationTests.SucroseIAM
    TomlIntegrationTests.SucrosePtb
    PROPERTIES
        RESOURCE_LOCK integration_sucrose_IAM_SF
)

# Tests using tests/TFVC
set_tests_properties(
    TomlIntegrationTests.TFVC
    TomlIntegrationTests.TFVCEcp
    PROPERTIES
        RESOURCE_LOCK integration_TFVC
)

# Tests using tests/RGBI_groups (RGBI_NH3Li/RGBI_NH3Li_ANO read nh3li.gbw from there too)
set_tests_properties(
    TomlIntegrationTests.RGBI_Groups_NH3BH3_sym
    TomlIntegrationTests.RGBI_Groups_NH3BH3_sym_ANO
    TomlIntegrationTests.RGBI_NH3Li
    TomlIntegrationTests.RGBI_NH3Li_ANO
    PROPERTIES
        RESOURCE_LOCK integration_RGBI_groups
)

# Tests using tests/SALTED
set_tests_properties(
    TomlIntegrationTests.SALTED
    SALTEDTests.ReadingSALTEDBinaryFile
    BesselTests.AnalyticFourier
    PROPERTIES
        RESOURCE_LOCK integration_SALTED
)

# Tests using tests/P1_test (all four write to the same tests/P1_test/NoSpherA2.log,
# NoSpherA2's log filename is not currently configurable via CLI -- without this lock,
# CTEST_PARALLEL_LEVEL > 1 races these against each other and corrupts the shared log
# with interleaved output from concurrent runs, on every platform)
set_tests_properties(
    TomlIntegrationTests.P1_test_XCW
    TomlIntegrationTests.P1_test_XCW_full
    TomlIntegrationTests.P1_test_XCW_h2
    TomlIntegrationTests.P1_test_XCW_h2_full
    PROPERTIES
        RESOURCE_LOCK integration_P1_test
)
