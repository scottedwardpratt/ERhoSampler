cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

project(EQCorr CXX)

add_custom_target(extern_msu_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_msu_sampler COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_sampler/software)
add_custom_target(extern_msu_eos COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_eos/software)
add_custom_target(extern_msu_erhosampler COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/ERhoSampler/software)

set(CMAKE_CXX_FLAGS "-O2 -Wall -std=c++17 -mfpmath=sse")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/msu_sampler/software/include ${GITHOME_MSU}/msu_eos/software/include ${GITHOME_MSU}/ERhoSampler/software/include)

link_directories(${GITHOME_MSU}/msu_sampler/software/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/msu_eos/software/lib ${GITHOME_MSU}/ERhoSampler/software/lib)

find_library(${GITHOME_MSU}/msu_sampler/software/lib msu_sampler)
find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/msu_eos/software/lib msu_eos)
find_library(${GITHOME_MSU}/ERhoSampler/software/lib msu_erhosampler)

#######
set(SOURCE_FILES ./corrtest.cc)
SET_SOURCE_FILES_PROPERTIES(corrtest.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a)
SET_SOURCE_FILES_PROPERTIES(corrtest.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a)
SET_SOURCE_FILES_PROPERTIES(corrtest.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a)
SET_SOURCE_FILES_PROPERTIES(corrtest.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/ERhoSampler/software/lib/libmsu_erhosampler.a)
	  
add_executable(corrtest ${SOURCE_FILES})
add_dependencies(corrtest extern_msu_sampler)
add_dependencies(corrtest extern_msu_eos)
add_dependencies(corrtest extern_msu_commonutils)
add_dependencies(corrtest extern_msu_erhosampler)
target_link_libraries(corrtest msu_sampler msu_commonutils msu_eos msu_erhosampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY} ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a ${GITHOME_MSU}/ERhoSampler/software/lib/libmsu_erhosampler.a)
########

########
set(SOURCE_FILES ./hydrocorr.cc)
SET_SOURCE_FILES_PROPERTIES(hydrocorr.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a)
SET_SOURCE_FILES_PROPERTIES(hydrocorr.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a)
SET_SOURCE_FILES_PROPERTIES(hydrocorr.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a)
SET_SOURCE_FILES_PROPERTIES(hydrocorr.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/ERhoSampler/software/lib/libmsu_erhosampler.a)
	
add_executable(hydrocorr ${SOURCE_FILES})
add_dependencies(hydrocorr extern_msu_sampler)
add_dependencies(hydrocorr extern_msu_eos)
add_dependencies(hydrocorr extern_msu_commonutils)
add_dependencies(hydrocorr extern_msu_erhosampler)
target_link_libraries(hydrocorr msu_sampler msu_commonutils msu_eos msu_erhosampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY} ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a ${GITHOME_MSU}/ERhoSampler/software/lib/libmsu_erhosampler.a)
########

########
set(SOURCE_FILES ./decaycorr.cc)
SET_SOURCE_FILES_PROPERTIES(decaycorr.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a)
SET_SOURCE_FILES_PROPERTIES(decaycorr.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a)
SET_SOURCE_FILES_PROPERTIES(decaycorr.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a)
SET_SOURCE_FILES_PROPERTIES(decaycorr.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/ERhoSampler/software/lib/libmsu_erhosampler.a)

add_executable(decaycorr ${SOURCE_FILES})
add_dependencies(decaycorr extern_msu_sampler)
add_dependencies(decaycorr extern_msu_eos)
add_dependencies(decaycorr extern_msu_commonutils)
add_dependencies(decaycorr extern_msu_erhosampler)
target_link_libraries(decaycorr msu_sampler msu_commonutils msu_eos msu_erhosampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY} ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a ${GITHOME_MSU}/ERhoSampler/software/lib/libmsu_erhosampler.a)