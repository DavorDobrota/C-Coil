add_library(ConversionUtil OBJECT ConversionUtil.cxx)

target_compile_options(ConversionUtil PRIVATE -fpic)

target_include_directories(
        ConversionUtil
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Coil
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Vector3
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Polynomial
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)