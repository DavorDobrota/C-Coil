cmake_minimum_required(VERSION 3.19)
project(General_Coil_Program)

set(CMAKE_CXX_STANDARD 17)

add_executable(General_Coil_Program src/main.cpp src/Coil.cxx include/Coil.h)
include_directories(include)

target_link_libraries(General_Coil_Program -L"C:/Users/unit5/Desktop/C++/FER/General Coil Program/lib")
target_link_libraries(General_Coil_Program hardware_acceleration)