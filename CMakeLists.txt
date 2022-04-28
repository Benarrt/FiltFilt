cmake_minimum_required(VERSION 3.0)
project(filtfilt)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

file(GLOB SRCFILES FilterX.hpp FiltFiltX.hpp SOFiltFilt.hpp)
add_executable(filtfilt "${SRCFILES}" main.cpp)

install(TARGETS filtfilt RUNTIME DESTINATION bin)
