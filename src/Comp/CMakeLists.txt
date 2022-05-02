add_library(
        Comp OBJECT

        CompMiscellaneous.cxx
        CompForce.cxx
        CompMInductanceZAxis.cxx
        CompMInductanceGeneral.cxx
        CompMTD.cxx)

target_compile_options(Comp PRIVATE -fpic)

target_include_directories(
        Comp
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)