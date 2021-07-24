add_library(Polynomial OBJECT Polynomial.cxx)

target_compile_options(Polynomial PRIVATE -fpic)

target_include_directories(
        Polynomial
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)