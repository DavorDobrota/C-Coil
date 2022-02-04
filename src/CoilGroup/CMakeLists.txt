add_library(CoilGroup OBJECT CoilGroup.cxx)

target_compile_options(CoilGroup PRIVATE -fpic)

target_include_directories(
        CoilGroup
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Coil
        PUBLIC "${CMAKE_SOURCE_DIR}/src/Tensor"
        PUBLIC "${CMAKE_SOURCE_DIR}/src/ThreadPool"
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)