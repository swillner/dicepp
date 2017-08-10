cmake_minimum_required(VERSION 2.8)
if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.8)
  cmake_policy(SET CMP0069 NEW) # for INTERPROCEDURAL_OPTIMIZATION
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project(dicepp)

set(DICEPP_VERSION 0.2.0)
if(EXISTS "${CMAKE_SOURCE_DIR}/.git" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/.git")
  find_program(HAVE_GIT git)
  if(HAVE_GIT)
    execute_process(COMMAND git describe --tags --dirty --always OUTPUT_VARIABLE GIT_OUTPUT OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE "^v([0-9]+\\.[0-9]+)\\.(0-)?([0-9]*)((-.+)?)$" "\\1.\\3\\4" DICEPP_VERSION ${GIT_OUTPUT})
  endif()
endif()
message(STATUS "Program version: ${DICEPP_VERSION}")

file(GLOB_RECURSE DICEPP_SOURCES src/*.cpp)
add_executable(dicepp ${DICEPP_SOURCES})

target_include_directories(dicepp PRIVATE include include/climate include/damage lib/cpp-library lib/settingsnode/include lib/settingsnode/lib/yaml-cpp/include)
target_compile_options(dicepp PRIVATE "-std=c++11")
target_compile_definitions(dicepp PRIVATE DICEPP_VERSION="${DICEPP_VERSION}")

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" OR CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
  if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.8)
    message(STATUS "Enabling interprocedural optimization")
    set_property(TARGET dicepp PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
  target_compile_definitions(dicepp PUBLIC NDEBUG)
else()
  target_compile_definitions(dicepp PRIVATE DEBUG)
endif()

option(DICEPP_WITH_NLOPT "NLopt (Optimizer)" ON)
if(DICEPP_WITH_NLOPT)
  find_package(NLOPT REQUIRED)
  message(STATUS "NLopt include directory: ${NLOPT_INCLUDE_DIR}")
  message(STATUS "NLopt library: ${NLOPT_LIBRARY}")
  target_link_libraries(dicepp nlopt)
  target_compile_definitions(dicepp PRIVATE DICEPP_WITH_NLOPT)
endif()

option(DICEPP_WITH_PAGMO "Pagmo (Optimizer)" OFF)
if(DICEPP_WITH_PAGMO)
  mark_as_advanced(Eigen3_DIR nlopt_DIR pagmo_DIR)
  find_package(pagmo REQUIRED)
  target_link_libraries(dicepp Pagmo::pagmo)
  target_compile_definitions(dicepp PRIVATE HAVE_CSTDDEF DICEPP_WITH_PAGMO)
endif()

option(DICEPP_WITH_NETCDF "NetCDF (Output)" ON)
if(DICEPP_WITH_NETCDF)
  find_package(NETCDF REQUIRED)
  message(STATUS "NetCDF include directory: ${NETCDF_INCLUDE_DIR}")
  message(STATUS "NetCDF library: ${NETCDF_LIBRARY}")
  target_link_libraries(dicepp netcdf)

  find_package(NETCDF_CPP4 REQUIRED)
  message(STATUS "NetCDF_c++4 include directory: ${NETCDF_CPP4_INCLUDE_DIR}")
  message(STATUS "NetCDF_c++4 library: ${NETCDF_CPP4_LIBRARY}")
  target_link_libraries(dicepp netcdf_c++4)

  target_compile_definitions(dicepp PRIVATE DICEPP_WITH_NETCDF)
endif()

include(lib/settingsnode/settingsnode.cmake)
target_link_libraries(dicepp settingsnode)
