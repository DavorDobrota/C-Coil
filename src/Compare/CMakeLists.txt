add_library(
        Compare OBJECT

        CompareMiscellaneous.cxx
        CompareForce.cxx
        CompareMInductanceZAxis.cxx
        CompareMInductanceGeneral.cxx
        CompareMTD.cxx)

target_compile_options(Compare PRIVATE -fpic)

target_include_directories(
        Compare
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)