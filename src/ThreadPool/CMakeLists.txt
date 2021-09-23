add_library(ThreadPool OBJECT ThreadPool.cxx)

target_compile_options(ThreadPool PRIVATE -fpic)

target_include_directories(
        ThreadPool
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)