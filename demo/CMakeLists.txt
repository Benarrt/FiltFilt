cmake_minimum_required(VERSION 3.0)
project(demo)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

find_library(FILTFILT_LIB libfiltfilt.so)
message(STATUS ${FILTFILT_LIB})

add_executable(demo "${SRCFILES}" main.cpp)
install(TARGETS demo RUNTIME DESTINATION bin)
target_link_libraries(demo PRIVATE ${FILTFILT_LIB})