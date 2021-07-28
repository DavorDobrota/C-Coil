add_library(Coil OBJECT Coil.cxx PrecisionArrays.cxx PrecisionFactor.cxx PrecisionArguments.cxx PrecisionGlobalVars.h CoilPairArguments.cxx MainCalculationMethods.cxx PointComputeMethods.cxx CalculateAllMethods.cxx MassComputeMethods.cxx CalculateMInductanceMethods.cxx SelfInductanceMethods.cxx ComputeMInductanceMethods.cxx CalculateForceMethods.cxx ComputeForceMethods.cxx CoilPairArgumentsCPUMethods.cxx CoilPairArgumentsGPUMethods.cxx)

target_compile_options(Coil PRIVATE -fpic)

target_include_directories(
        Coil
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
        # TODO: this can't be standard... fix it
        PUBLIC ${CMAKE_SOURCE_DIR}/src/CoilData
        PUBLIC ${CMAKE_SOURCE_DIR}/src/ConversionUtil
        PUBLIC ${CMAKE_SOURCE_DIR}/src/LegendreMatrix
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)