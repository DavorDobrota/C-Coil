add_library(
        CoilGroupCUDAKernels OBJECT

        FieldAGroup.cu
        FieldBGroup.cu
        FieldEGroup.cu
        FieldGGroup.cu
        MInductanceGroup.cu
        ForceGroup.cu
)

target_compile_options(CoilGroupCUDAKernels PRIVATE "--compiler-options=-fpic")

target_include_directories(
        CoilGroupCUDAKernels
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)
