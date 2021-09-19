add_library(
        Test OBJECT

        TestMiscellaneous.cxx BabicTestCases.cxx TestForce.cxx
        TestFieldsPerformance.cxx TestMInductanceZAxis.cxx
        TestMInductanceGeneral.cxx ForceTestCases.cxx
        TestArgumentGeneration.cxx)

target_compile_options(Test PRIVATE -fpic)

target_include_directories(
        Test
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Coil
        PUBLIC ${CMAKE_SOURCE_DIR}/src/OldCoil
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Polynomial
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Tensor
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)