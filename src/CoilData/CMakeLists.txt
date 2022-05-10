add_library(CoilData OBJECT CoilData.cxx)

target_compile_options(CoilData PRIVATE -fpic)

target_include_directories(
        CoilData
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)