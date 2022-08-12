add_library(
        Benchmark OBJECT

        BenchmarkMiscellaneous.cxx
        Force/BenchmarkForce.cxx
        Fields/BenchmarkFields.cxx
        MInductance/BenchmarkMInductanceZAxis.cxx
        MInductance/BenchmarkMInductanceGeneral.cxx
        BenchmarkMTD.cxx Fields/BenchmarkFieldsEveryCoil.cxx Fields/BenchFieldsScaling.cxx)
set_property(TARGET Benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(Benchmark PRIVATE -fpic)

target_include_directories(
        Benchmark
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)