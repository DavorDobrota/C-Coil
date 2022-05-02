cmake_minimum_required(VERSION 3.8)

project(General_Coil_Program)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(WIN32)
    set(CMAKE_CXX_FLAGS "/std:c++17 /O2 /arch:AVX2 /arch:AVX /fp:fast")
else()
    set(CMAKE_CXX_FLAGS "-Ofast -mavx2 -ffast-math")
endif()

option(USE_GPU "Enable GPU acceleration" OFF)

option(OUTPUT_LIB_SHARED "Build a shared library" ON)
option(OUTPUT_EXECUTABLE "Build a test executable" ON)

# TODO: add static lib support

if(USE_GPU)
    add_compile_definitions(USE_GPU)
endif()

#######################################
# Add all extern include paths
#######################################
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/extern/CTPL")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/extern/HardwareAcceleration")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/include")
#######################################

#######################################
# Add all local include paths
#######################################
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Test")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Benchmark")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Coil")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/CoilData")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/CoilGroup")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Comp")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/LegendreMatrix")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/ThreadPool")
#######################################

#######################################
# Add all libraries
#######################################
add_subdirectory(src/Benchmark)
list(APPEND LIBS "$<TARGET_OBJECTS:Benchmark>")

add_subdirectory(src/Comp)
list(APPEND LIBS "$<TARGET_OBJECTS:Comp>")

add_subdirectory(src/Test)
list(APPEND LIBS "$<TARGET_OBJECTS:Test>")

add_subdirectory(src/Coil)
list(APPEND LIBS "$<TARGET_OBJECTS:Coil>")

add_subdirectory(src/CoilData)
list(APPEND LIBS "$<TARGET_OBJECTS:CoilData>")

add_subdirectory(src/CoilGroup)
list(APPEND LIBS "$<TARGET_OBJECTS:CoilGroup>")

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
            PUBLIC ${INCS}
    )
endif()

if(OUTPUT_LIB_SHARED)
    add_library(coil_calculations_shared SHARED ${LIBS})
    target_include_directories(
            coil_calculations_shared
            BEFORE
            PUBLIC ${INCS}
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

if(!WIN32)
    if(OUTPUT_EXECUTABLE)
        target_link_libraries(test_executable pthread)
    endif()

    if(OUTPUT_LIB_SHARED)
        target_link_libraries(coil_calculations_shared pthread)
    endif()
endif()
