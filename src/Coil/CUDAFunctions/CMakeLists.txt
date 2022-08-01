add_library(
        CUDAFunctions OBJECT

        CoilKernels/FieldA.cu
        CoilKernels/FieldB.cu
        CoilKernels/FieldG.cu
        CoilKernels/MInductance.cu
        CoilKernels/Force.cu
        ErrorCheck/ErrorCheck.cu
        MemoryManagement/GPUMemoryManagement.cu
        CoilGroupKernels/FieldAGroup.cu
        CoilGroupKernels/FieldBGroup.cu
        CoilGroupKernels/FieldGGroup.cu CoilGroupKernels/MInductanceGroup.cu CoilGroupKernels/ForceGroup.cu)

target_compile_options(CUDAFunctions PRIVATE "--compiler-options=-fpic")

target_include_directories(
        CUDAFunctions
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)
