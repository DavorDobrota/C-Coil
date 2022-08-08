add_library(
        CoilGroup OBJECT

        CoilGroup.cxx
        Fields/CalculateFieldsGPU.cxx
        Fields/CalculateFieldsMTD.cxx MInductance/CalculateMInductanceMTD.cxx MInductance/CalculateMInductanceGPU.cxx ForceAndTorque/CalculateForceMTD.cxx ForceAndTorque/CalculateForceGPU.cxx Utils/GPUArgumentGeneration.cxx Fields/ComputeFields.cxx ForceAndTorque/ComputeForce.cxx MInductance/ComputeMInductance.cxx)
set_property(TARGET CoilGroup PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(CoilGroup PRIVATE -fpic)

target_include_directories(
        CoilGroup
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)