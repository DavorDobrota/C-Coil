add_library(Tensor OBJECT Vector/Vector3.cxx Matrix/Matrix3.cxx Vector/Vector3Array.cxx Matrix/Matrix3Array.cxx)
set_property(TARGET Tensor PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(Tensor PRIVATE -fpic)

target_include_directories(
        Tensor
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)