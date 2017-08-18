cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(cmake/helpers.cmake)

set_default_build_type(Release)

project(dicepp)

file(GLOB_RECURSE DICEPP_SOURCES src/*.cpp)
add_executable(dicepp ${DICEPP_SOURCES})

target_include_directories(dicepp PRIVATE include include/climate include/damage lib/cpp-library)
target_compile_options(dicepp PRIVATE "-std=c++11")

set_advanced_cpp_warnings(dicepp)
set_build_type_specifics(dicepp)
add_git_version(dicepp
  FALLBACK_VERSION 1.0.0)

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

add_cpp_tools(dicepp)
