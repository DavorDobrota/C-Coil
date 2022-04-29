add_library(
        Test OBJECT

        BenchmarkMiscellaneous.cxx
        BenchmarkForce.cxx
        BenchmarkFields.cxx
        BenchmarkMInductanceZAxis.cxx
        BenchmarkMInductanceGeneral.cxx
        BenchmarkMTD.cxx)

target_compile_options(Benchmark PRIVATE -fpic)

target_include_directories(
        Benchmark
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)