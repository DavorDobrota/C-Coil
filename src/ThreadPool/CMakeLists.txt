add_library(ThreadPool OBJECT ThreadPool.cxx)
set_property(TARGET ThreadPool PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(ThreadPool PRIVATE -fpic)

target_include_directories(
        ThreadPool
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)