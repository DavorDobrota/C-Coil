add_library(
        Benchmark OBJECT

        Other/BenchMathFunctions.cxx
        Coil/MInductanceAndForce/BenchCoilForce.cxx
        Coil/Fields/BenchCoilFields.cxx
        Coil/Fields/BenchCoilFieldsEveryCoil.cxx
        Coil/Fields/BenchCoilFieldsScaling.cxx
        Coil/MInductanceAndForce/BenchCoilMInductanceZAxis.cxx
        Coil/MInductanceAndForce/BenchCoilMInductanceGeneral.cxx
        Coil/MInductanceAndForce/BenchCoilComputeAll.cxx
        CoilGroup/BenchGroupFields.cxx
        CoilGroup/BenchGroupMInductanceAndForce.cxx
        CoilGroup/BenchGroupFieldsScaling.cxx)

set_property(TARGET Benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(Benchmark PRIVATE -fpic)

target_include_directories(
        Benchmark
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)