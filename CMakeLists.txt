cmake_minimum_required(VERSION 3.8)

project(C-Coil)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0069 NEW)

#######################################
# Define all build options
#######################################
option(USE_GPU "Enable GPU acceleration" OFF)

set(GPU_INCREMENTS "80" CACHE STRING "GPU increment count")

set(TYPE "float" CACHE STRING "Floating point precision (can be float or double)")

option(OUTPUT_LIB_SHARED "Build a shared library" ON)
option(OUTPUT_EXECUTABLE "Build a test executable" ON)
option(OUTPUT_DOCS "Build the documentation" ON)
# TODO: add static lib support
#######################################

#######################################
# Handle chosen options
#######################################
if(WIN32)
    set(CMAKE_CXX_FLAGS "/std:c++17 /O2 /arch:AVX2 /fp:fast")
else()
    set(CMAKE_CXX_FLAGS "-Ofast -mavx2 -ffast-math -flto=auto")
endif()

if(USE_GPU)
    add_compile_definitions(USE_GPU)

    set(CUDA_NVCC_FLAGS "--use_fast_math --optimise -O3")
    set(CMAKE_CUDA_ARCHITECTURES "61;72;75;86")

    project(C-Coil LANGUAGES CXX CUDA)
endif()

if(${GPU_INCREMENTS} LESS_EQUAL 80 AND ${GPU_INCREMENTS} GREATER_EQUAL 1)
    add_compile_definitions(GPU_INCREMENTS=${GPU_INCREMENTS})
else()
    message(FATAL_ERROR "GPU_INCREMENTS must be set to an integer in the range [1, 80]!")
endif()

add_compile_definitions(TYPE=${TYPE})

if(OUTPUT_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen ALL
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM )
    else()
        message("Doxygen needs to be installed to generate the documentation!")
    endif ()
endif()
#######################################

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
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Coil/CUDAKernels")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/CoilGroup")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/CoilGroup/CUDAKernels")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/CUDAUtils")
list(APPEND INCS "${CMAKE_CURRENT_SOURCE_DIR}/src/Compare")
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

add_subdirectory(src/LegendreMatrix)
list(APPEND LIBS "$<TARGET_OBJECTS:LegendreMatrix>")

add_subdirectory(src/Tensor)
list(APPEND LIBS "$<TARGET_OBJECTS:Tensor>")

add_subdirectory(src/Utils)
list(APPEND LIBS "$<TARGET_OBJECTS:Utils>")

add_subdirectory(src/ThreadPool)
list(APPEND LIBS "$<TARGET_OBJECTS:ThreadPool>")

if(USE_GPU)
    add_subdirectory(src/Coil/CUDAKernels)
    list(APPEND LIBS "$<TARGET_OBJECTS:CoilCUDAKernels>")

    add_subdirectory(src/CoilGroup/CUDAKernels)
    list(APPEND LIBS "$<TARGET_OBJECTS:CoilGroupCUDAKernels>")

    add_subdirectory(src/CUDAUtils)
    list(APPEND LIBS "$<TARGET_OBJECTS:CUDAUtils>")
endif()
#######################################

if(OUTPUT_EXECUTABLE)
    add_executable(test_executable src/main.cxx ${LIBS})
    set_property(TARGET test_executable PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
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
    add_library(c_coil_shared SHARED ${LIBS})
    set_property(TARGET c_coil_shared PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    target_include_directories(
            c_coil_shared
            BEFORE
            PUBLIC ${INCS}
    )
#    if(USE_GPU)
#        set_target_properties(coil_calculations_shared PROPERTIES CUDA_ARCHITECTURES "all-major")
#    endif()
    if(!WIN32)
        target_link_libraries(c_coil_shared pthread)
    endif()
endif()
