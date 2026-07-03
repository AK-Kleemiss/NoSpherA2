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

# Tests using tests/RGBI_groups
set_tests_properties(
    TomlIntegrationTests.RGBI_Groups_NH3BH3
    TomlIntegrationTests.RGBI_Groups_NH3BH3_ANO
    TomlIntegrationTests.RGBI_Groups_NH3BH3_sym
    TomlIntegrationTests.RGBI_Groups_NH3BH3_sym_ANO
    TomlIntegrationTests.RGBI_Groups_NH3Li
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
