add_library(
        CoilCUDAKernels OBJECT

        FieldA.cu
        FieldB.cu
        FieldG.cu
        MInductance.cu
        Force.cu
)

target_compile_options(CoilCUDAKernels PRIVATE "--compiler-options=-fpic")

target_include_directories(
        CoilCUDAKernels
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)
