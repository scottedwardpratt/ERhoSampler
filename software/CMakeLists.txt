cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})
project(msu_erhosampler CXX C)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -std=c++17")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/msu_eos/software/include ${GITHOME_MSU}/msu_sampler/software/include)

file(GLOB erhosampler_src ${CMAKE_CURRENT_SOURCE_DIR}/src/**.cc)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(msu_erhosampler STATIC ${erhosampler_src})
target_link_libraries(msu_erhosampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
