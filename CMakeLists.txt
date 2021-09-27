cmake_minimum_required(VERSION 3.8)

project(General_Coil_Program)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Ofast -msse -msse2 -msse3 -msse4 -msse4a -msse4.1 -msse4.2 -mavx -mavx2 -ffast-math")

option(USE_GPU "Enable GPU acceleration" OFF)

option(OUTPUT_LIB_SHARED "Build a shared library" ON)
option(OUTPUT_EXECUTABLE "Build a test executable" ON)

# TODO: add static lib support

if(USE_GPU)
    add_compile_definitions(USE_GPU)
endif()

#######################################
add_subdirectory(src/Test)
list(APPEND LIBS "$<TARGET_OBJECTS:Test>")

add_subdirectory(src/Coil)
list(APPEND LIBS "$<TARGET_OBJECTS:Coil>")

add_subdirectory(src/CoilData)
list(APPEND LIBS "$<TARGET_OBJECTS:CoilData>")

add_subdirectory(src/CoilGroup)
list(APPEND LIBS "$<TARGET_OBJECTS:CoilGroup>")

add_subdirectory(src/Math)
list(APPEND LIBS "$<TARGET_OBJECTS:CustomMath>")

add_subdirectory(src/LegendreMatrix)
list(APPEND LIBS "$<TARGET_OBJECTS:LegendreMatrix>")

add_subdirectory(src/Tensor)
list(APPEND LIBS "$<TARGET_OBJECTS:Tensor>")

add_subdirectory(src/ThreadPool)
list(APPEND LIBS "$<TARGET_OBJECTS:ThreadPool>")
#######################################

if(OUTPUT_EXECUTABLE)
    add_executable(test_executable src/main.cxx ${LIBS})
    target_include_directories(
            test_executable
            BEFORE
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/Coil"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/CoilGroup"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/Test"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor"
    )
endif()

if(OUTPUT_LIB_SHARED)
    add_library(coil_calculations_shared SHARED ${LIBS})
    target_include_directories(
            coil_calculations_shared
            BEFORE
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/Coil"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/CoilGroup"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/Test"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor"
    )
endif()

if(USE_GPU)
    if(WIN32)
        add_library(hardware_acceleration SHARED IMPORTED)
        set_property(TARGET hardware_acceleration PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/hardware_acceleration.dll")
        set_property(TARGET hardware_acceleration PROPERTY IMPORTED_IMPLIB "${CMAKE_CURRENT_SOURCE_DIR}/lib/hardware_acceleration.lib")
    else()
        add_library(hardware_acceleration SHARED IMPORTED)
        set_property(TARGET hardware_acceleration PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/lib/libhardware_acceleration_shared.so")
    endif()

    if(OUTPUT_EXECUTABLE)
        target_link_libraries(test_executable hardware_acceleration)
    endif()

    if(OUTPUT_LIB_SHARED)
        target_link_libraries(coil_calculations_shared hardware_acceleration)
    endif()
endif()

if(OUTPUT_EXECUTABLE)
    target_link_libraries(test_executable pthread)
endif()

if(OUTPUT_LIB_SHARED)
    target_link_libraries(coil_calculations_shared pthread)
endif()
