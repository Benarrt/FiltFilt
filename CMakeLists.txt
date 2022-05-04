cmake_minimum_required(VERSION 3.0)
project(filtfilt)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(filtfilt SHARED FiltFilt.cpp)

install(TARGETS filtfilt DESTINATION /usr/lib)
install(FILES include/FiltFilt.h DESTINATION include)
install(FILES include/EigenUtil.h DESTINATION include)
