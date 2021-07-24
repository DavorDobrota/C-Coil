add_library(Coil OBJECT Coil.cxx)

target_compile_options(Coil PRIVATE -fpic)

target_include_directories(
        Coil
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/src/ConversionUtil
        PUBLIC ${CMAKE_SOURCE_DIR}/src/LegendreMatrix
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)