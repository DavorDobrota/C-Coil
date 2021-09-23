add_library(
        Coil OBJECT

        Coil.cxx
        PrecisionArrays.cxx
        PrecisionFactor.cxx
        PrecisionArguments.cxx
        CoilPairArguments.cxx
        CalculateFieldsMethods.cxx
        ComputePointMethods.cxx
        CalculateAllFieldsMethods.cxx
        ComputeAllPointsMethods.cxx
        CalculateMInductanceMethods.cxx
        ComputeMInductanceMethods.cxx
        CalculateForceMethods.cxx
        ComputeForceMethods.cxx
        CalculateFieldsSlow.cxx
        CalculateFieldsFast.cxx)

target_compile_options(Coil PRIVATE -fpic)

target_include_directories(
        Coil
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${CMAKE_SOURCE_DIR}/src/CoilData
        PUBLIC ${CMAKE_SOURCE_DIR}/src/Tensor
        PUBLIC ${CMAKE_SOURCE_DIR}/src/ConversionUtil
        PUBLIC ${CMAKE_SOURCE_DIR}/src/LegendreMatrix
        PUBLIC ${CMAKE_SOURCE_DIR}/src/ThreadPool
        PUBLIC ${CMAKE_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/src
)