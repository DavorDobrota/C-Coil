add_library(
        Test OBJECT

        TestMiscellaneous.cxx
        TestMInductanceGeneral.cxx
        TestForce.cxx
        TestMInductanceZAxis.cxx
        TestArgumentGeneration.cxx)
set_property(TARGET Test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(Test PRIVATE -fpic)

target_include_directories(
        Test
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)