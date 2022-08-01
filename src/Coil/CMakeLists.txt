add_library(
        Coil OBJECT

        Coil.cxx
        PrecisionArguments/PrecisionArrays.cxx
        PrecisionArguments/PrecisionFactor.cxx
        PrecisionArguments/PrecisionArguments.cxx
        PrecisionArguments/CoilPairArguments.cxx
        Fields/CalculateFieldsMethods.cxx
        Fields/CalculateAllFieldsMethods.cxx
        MInductance/CalculateMInductanceMethods.cxx
        MInductance/CalculateMInductanceArrangements.cxx
        ForceAndTorque/CalculateForceMethods.cxx
        ForceAndTorque/CalculateForceArrangements.cxx
        Fields/CalculateFieldsSlow.cxx
        Fields/CalculateFieldsFast.cxx
        Utils/SupportFunctions.cxx
        Utils/CalculateAttributes.cxx
)
set_property(TARGET Coil PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(Coil PRIVATE -fpic)

target_include_directories(
        Coil
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)