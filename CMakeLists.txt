cmake_minimum_required(VERSION 3.8)

project(Coil_Evolution)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(WIN32)
    set(CMAKE_CXX_FLAGS "/std:c++17 /O2 /arch:AVX2 /fp:fast")
else()
    set(CMAKE_CXX_FLAGS "-Ofast -mavx2 -ffast-math")
endif()

option(USE_GPU "Enable GPU acceleration" OFF)

option(OUTPUT_LIB_SHARED "Build a shared library" ON)
option(OUTPUT_EXECUTABLE "Build a test executable" ON)

# TODO: add static lib support

if(USE_GPU)
    add_compile_definitions(USE_GPU)

    set(CUDA_NVCC_FLAGS "--use_fast_math --optimise -O3")
    set(CMAKE_CUDA_ARCHITECTURES "52;61;72;75;86")

    project(Coil_Evolution LANGUAGES CXX CUDA)
endif()

#######################################
# Add all extern include paths
#######################################
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/extern/CTPL")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/include")
#######################################

#######################################
# Add all local include paths
#######################################
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Test")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Benchmark")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Coil")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/CoilGroup")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Compare")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/CUDAFunctions")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/LegendreMatrix")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Tensor")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/ThreadPool")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Utils")
#######################################

#######################################
# Add all libraries
#######################################
add_subdirectory(src/Benchmark)
list(APPEND LIBS "$<TARGET_OBJECTS:Benchmark>")

add_subdirectory(src/Compare)
list(APPEND LIBS "$<TARGET_OBJECTS:Compare>")

add_subdirectory(src/Test)
list(APPEND LIBS "$<TARGET_OBJECTS:Test>")

add_subdirectory(src/Coil)
list(APPEND LIBS "$<TARGET_OBJECTS:Coil>")

add_subdirectory(src/CoilGroup)
list(APPEND LIBS "$<TARGET_OBJECTS:CoilGroup>")

if(USE_GPU)
    add_subdirectory(src/CUDAFunctions)
    list(APPEND LIBS "$<TARGET_OBJECTS:CUDAFunctions>")
endif()

add_subdirectory(src/LegendreMatrix)
list(APPEND LIBS "$<TARGET_OBJECTS:LegendreMatrix>")

add_subdirectory(src/Tensor)
list(APPEND LIBS "$<TARGET_OBJECTS:Tensor>")

add_subdirectory(src/Utils)
list(APPEND LIBS "$<TARGET_OBJECTS:Utils>")

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
#    if(USE_GPU)
#        set_target_properties(test_executable PROPERTIES CUDA_ARCHITECTURES "all-major")
#    endif()
    if(!WIN32)
        target_link_libraries(test_executable pthread)
    endif()
endif()

if(OUTPUT_LIB_SHARED)
    add_library(coil_calculations_shared SHARED ${LIBS})
    target_include_directories(
            coil_calculations_shared
            BEFORE
            PUBLIC ${INCS}
    )
#    if(USE_GPU)
#        set_target_properties(coil_calculations_shared PROPERTIES CUDA_ARCHITECTURES "all-major")
#    endif()
    if(!WIN32)
        target_link_libraries(coil_calculations_shared pthread)
    endif()
endif()
