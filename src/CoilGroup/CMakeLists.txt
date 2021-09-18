add_library(CoilGroup OBJECT CoilGroup.cxx)

target_compile_options(CoilData PRIVATE -fpic)

target_include_directories(
        CoilGroup
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Coil
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)