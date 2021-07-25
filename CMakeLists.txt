cmake_minimum_required(VERSION 3.8)

project(General_Coil_Program)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Ofast")

#######################################
add_subdirectory(src/Test)
list(APPEND LIBS "$<TARGET_OBJECTS:Test>")

add_subdirectory(src/Coil)
list(APPEND LIBS "$<TARGET_OBJECTS:Coil>")

add_subdirectory(src/CoilData)
list(APPEND LIBS "$<TARGET_OBJECTS:CoilData>")

add_subdirectory(src/OldCoil)
list(APPEND LIBS "$<TARGET_OBJECTS:OldCoil>")

add_subdirectory(src/LegendreMatrix)
list(APPEND LIBS "$<TARGET_OBJECTS:LegendreMatrix>")

add_subdirectory(src/Polynomial)
list(APPEND LIBS "$<TARGET_OBJECTS:Polynomial>")

add_subdirectory(src/ConversionUtil)
list(APPEND LIBS "$<TARGET_OBJECTS:ConversionUtil>")
#######################################

if(WIN32)
    add_library(hardware_acceleration SHARED IMPORTED)
    set_property(TARGET hardware_acceleration PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/hardware_acceleration.dll")
    set_property(TARGET hardware_acceleration PROPERTY IMPORTED_IMPLIB "${CMAKE_CURRENT_SOURCE_DIR}/lib/hardware_acceleration.lib")
else()
    add_library(hardware_acceleration SHARED IMPORTED)
    set_property(TARGET hardware_acceleration PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/lib/libhardware_acceleration_shared.so")
endif()

add_executable(General_Coil_Program src/main.cxx ${LIBS})
target_include_directories(
        General_Coil_Program
        BEFORE
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src"
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/Coil"
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/OldCoil"
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/Polynomial"
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/Test"
)

if(WIN32)
    target_link_libraries(General_Coil_Program hardware_acceleration)
else()
    target_link_libraries(General_Coil_Program hardware_acceleration pthread)
endif()
