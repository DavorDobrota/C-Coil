add_library(
        Test OBJECT

        Other/TestMiscellaneous.cxx
        MInductanceAndForce/TestMInductanceGeneral.cxx
        MInductanceAndForce/TestForce.cxx
        MInductanceAndForce/TestMInductanceZAxis.cxx
        Other/TestArgumentGeneration.cxx MInductanceAndForce/TestArrangements.cxx Other/TestFields.cxx)
set_property(TARGET Test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(Test PRIVATE -fpic)

target_include_directories(
        Test
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)