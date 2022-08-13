add_library(
        Benchmark OBJECT

        Other/BenchMathFunctions.cxx
        Coil/MInductanceAndForce/BenchForce.cxx
        Coil/Fields/BenchFields.cxx
        Coil/Fields/BenchFieldsEveryCoil.cxx
        Coil/Fields/BenchFieldsScaling.cxx
        Coil/MInductanceAndForce/BenchMInductanceZAxis.cxx
        Coil/MInductanceAndForce/BenchMInductanceGeneral.cxx
        Coil/MInductanceAndForce/BenchComputeAll.cxx
        CoilGroup/GroupFields.cxx
        CoilGroup/GroupMInductanceAndForce.cxx
    )

set_property(TARGET Benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(Benchmark PRIVATE -fpic)

target_include_directories(
        Benchmark
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)