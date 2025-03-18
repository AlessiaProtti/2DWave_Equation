cmake_minimum_required(VERSION 3.30)
project(2DWave_Equation C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
set(CMAKE_C_STANDARD 11)

add_executable(2DWave_Equation Serial/main.c)
target_link_libraries(2DWave_Equation m)