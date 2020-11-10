cmake_minimum_required(VERSION 3.10)
project(power_method C)

set(CMAKE_C_STANDARD 99)

add_executable(power_method main.c test.cpp)
target_link_libraries(power_method  m)