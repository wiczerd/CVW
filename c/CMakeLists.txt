cmake_minimum_required(VERSION 3.13)
project(c C)

set(CMAKE_C_STANDARD 99)

add_executable(CVW main.c)

include_directories(/usr/include/gsl)

#add_library( gsl SHARED IMPORTED )
#target_link_libraries(gsl)

find_package(GSL REQUIRED)

target_link_libraries(CVW GSL::gsl GSL::gslcblas m gomp)