add_library(Tensor OBJECT CoordVector3.cxx FieldVector3.cxx Matrix3.cxx)

target_compile_options(Tensor PRIVATE -fpic)

target_include_directories(
        Tensor
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)