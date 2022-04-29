add_library(CustomMath OBJECT CustomMath.cxx)

target_compile_options(CustomMath PRIVATE -fpic)

target_include_directories(
        CustomMath
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)