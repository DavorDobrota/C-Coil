add_library(CoilGroup OBJECT CoilGroup.cxx)

target_compile_options(CoilGroup PRIVATE -fpic)

target_include_directories(
        CoilGroup
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)