add_library(Utils OBJECT Timing.cxx)

target_compile_options(Utils PRIVATE -fpic)

target_include_directories(
        Utils
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)