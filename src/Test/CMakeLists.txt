add_library(
        Test OBJECT

        TestMiscellaneous.cxx MInductanceGeneralTestCases.cxx TestForce.cxx
        TestFieldsPerformance.cxx TestMInductanceZAxis.cxx
        TestMInductanceGeneral.cxx ForceTestCases.cxx
        TestArgumentGeneration.cxx TestMTD.cxx)

target_compile_options(Test PRIVATE -fpic)

target_include_directories(
        Test
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)