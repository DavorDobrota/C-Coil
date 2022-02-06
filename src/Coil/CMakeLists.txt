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
        PUBLIC ${INCS}
)