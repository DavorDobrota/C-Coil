add_library(OldCoil OBJECT OldCoil.cxx)

target_compile_options(OldCoil PRIVATE -fpic)

target_include_directories(
        OldCoil
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)