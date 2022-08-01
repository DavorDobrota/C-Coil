add_library(
        CoilGroup OBJECT

        CoilGroup.cxx
        CalculateFieldsGPU.cxx
        ComputeFieldsMethods.cxx
        CalculateMTD.cxx)
set_property(TARGET CoilGroup PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(CoilGroup PRIVATE -fpic)

target_include_directories(
        CoilGroup
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)