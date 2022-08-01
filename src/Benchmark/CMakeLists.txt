add_library(
        Benchmark OBJECT

        BenchmarkMiscellaneous.cxx
        BenchmarkForce.cxx
        BenchmarkFields.cxx
        BenchmarkMInductanceZAxis.cxx
        BenchmarkMInductanceGeneral.cxx
        BenchmarkMTD.cxx)
set_property(TARGET Benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(Benchmark PRIVATE -fpic)

target_include_directories(
        Benchmark
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)