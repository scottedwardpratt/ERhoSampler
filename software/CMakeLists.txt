cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

project(newsamplertest CXX)

add_custom_target(extern_msu_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_msu_sampler COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_sampler/software)
add_custom_target(extern_msu_eos COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_eos/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/msu_sampler/software/include ${GITHOME_MSU}/msu_eos/software/include)

link_directories(${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/msu_eos/software/lib ${GITHOME_MSU}/msu_sampler/software/lib)

set(SOURCE_FILES ./newsamplertest.cc)
SET_SOURCE_FILES_PROPERTIES(
      newsamplertest.cc PROPERTIES OBJECT_DEPENDS
    ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a
 )
SET_SOURCE_FILES_PROPERTIES(
      newsamplertest.cc PROPERTIES OBJECT_DEPENDS
    ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a
 )
add_executable(newsamplertest ${SOURCE_FILES})
SET_SOURCE_FILES_PROPERTIES(
      newsamplertest.cc PROPERTIES OBJECT_DEPENDS
    ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a
 )
add_executable(newsamplertest ${SOURCE_FILES})

add_dependencies(newsamplertest extern_msu_commonutils)
add_dependencies(newsamplertest extern_msu_eos)
add_dependencies(newsamplertest extern_msu_sampler)

find_library(${GITHOME_MSU}/msu_sampler/software/lib msu_sampler)
find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/msu_eos/software/lib msu_eos)
target_link_libraries(newsamplertest msu_sampler msu_commonutils msu_eos ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
