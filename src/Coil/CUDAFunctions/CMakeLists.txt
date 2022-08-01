add_library(
        CUDAFunctions OBJECT

        CoilKernels/FieldA.cu
        CoilKernels/FieldB.cu
        CoilKernels/FieldG.cu
        CoilKernels/MInductance.cu
        CoilKernels/Force.cu
        ErrorCheck/ErrorCheck.cu
        MemoryManagement/GPUMemoryManagement.cu
        CoilGroupKernels/AFieldGroup.cu
        CoilGroupKernels/BFieldGroup.cu
        CoilGroupKernels/GFieldGroup.cu CoilGroupKernels/MInductanceGroup.cu)

target_compile_options(CUDAFunctions PRIVATE "--compiler-options=-fpic")

target_include_directories(
        CUDAFunctions
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)
