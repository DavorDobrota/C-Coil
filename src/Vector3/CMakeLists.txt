add_library(Vector3 OBJECT CoordVector3.cxx FieldVector3.cxx)

target_compile_options(Vector3 PRIVATE -fpic)

target_include_directories(
        Vector3
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Coil
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)