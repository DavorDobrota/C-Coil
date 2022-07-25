add_library(Tensor OBJECT Vector3.cxx Matrix3.cxx Vector3Array.cxx Matrix3Array.cxx)

target_compile_options(Tensor PRIVATE -fpic)

target_include_directories(
        Tensor
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)