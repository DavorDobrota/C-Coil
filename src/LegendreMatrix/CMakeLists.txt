add_library(LegendreMatrix OBJECT LegendreMatrix.cxx)
set_property(TARGET LegendreMatrix PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(LegendreMatrix PRIVATE -fpic)

target_include_directories(
        LegendreMatrix
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC ${INCS}
)