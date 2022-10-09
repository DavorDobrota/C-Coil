add_library(
        Compare OBJECT

        Other/CompMiscellaneous.cxx
        MInductanceAndForce/CompForce.cxx
        MInductanceAndForce/CompMInductanceZAxis.cxx
        MInductanceAndForce/CompMInductanceGeneral.cxx
        Other/CompFieldMTD.cxx)
set_property(TARGET Compare PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(Compare PRIVATE -fpic)

target_include_directories(
        Compare
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)