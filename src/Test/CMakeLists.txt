add_library(
        Test OBJECT

        TestMiscellaneous.cxx
        MInductanceGeneralTestCases.cxx
        TestMInductanceGeneral.cxx
        ForceTestCases.cxx
        TestForce.cxx
        TestMInductanceZAxis.cxx
        TestArgumentGeneration.cxx
        TestMTD.cxx)

target_compile_options(Test PRIVATE -fpic)

target_include_directories(
        Test
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)