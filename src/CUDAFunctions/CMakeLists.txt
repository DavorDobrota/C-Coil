add_library(CUDAFunctions OBJECT AField.cu BField.cu GField.cu ErrorCheck.cu)

target_compile_options(CUDAFunctions PRIVATE "--compiler-options=-fpic")

target_include_directories(
        CUDAFunctions
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)
