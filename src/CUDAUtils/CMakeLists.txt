add_library(
        CUDAUtils OBJECT

        ErrorCheck/ErrorCheck.cu
        MemoryManagement/GPUMemoryManagement.cu
)

target_compile_options(CUDAUtils PRIVATE "--compiler-options=-fpic")

target_include_directories(
        CUDAUtils
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)
