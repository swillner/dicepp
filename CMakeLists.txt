# Copyright (C) 2019 Sven Willner <sven.willner@pik-potsdam.de>
#
# This file is part of DICE++.
#
# DICE++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# DICE++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with DICE++.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(cmake/helpers.cmake)
include(cmake/libraries.cmake)

set_default_build_type(Release)

project(dicepp)

file(GLOB_RECURSE DICEPP_SOURCES src/*.cpp)
add_executable(dicepp ${DICEPP_SOURCES})

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  target_compile_options(dicepp PUBLIC -fno-math-errno -fno-rtti)
  add_custom_command(
    TARGET dicepp
    POST_BUILD
    COMMAND ${CMAKE_STRIP} dicepp
  )
endif()

target_include_directories(dicepp PRIVATE include include/climate include/damage lib/cpp-library)
target_compile_features(dicepp PUBLIC cxx_std_14)

add_system_paths()
set_advanced_cpp_warnings(dicepp)
set_build_type_specifics(dicepp)
add_git_version(dicepp WITH_DIFF)

option(DICEPP_WITH_BORG "Borg (Optimizer)" OFF)
if(DICEPP_WITH_BORG)
  target_sources(dicepp PUBLIC ${CMAKE_SOURCE_DIR}/lib/borg/borg.c ${CMAKE_SOURCE_DIR}/lib/borg/mt19937ar.c)
  target_include_directories(dicepp PRIVATE lib/borg)
  target_compile_definitions(dicepp PRIVATE DICEPP_WITH_BORG)
endif()

option(DICEPP_WITH_MIDACO "Midaco (Optimizer)" OFF)
if(DICEPP_WITH_MIDACO)
  target_sources(dicepp PUBLIC ${CMAKE_SOURCE_DIR}/lib/midaco/midaco.c)
  target_compile_definitions(dicepp PRIVATE DICEPP_WITH_MIDACO)
endif()

option(DICEPP_WITH_NLOPT "NLopt (Optimizer)" ON)
if(DICEPP_WITH_NLOPT)
  include_nlopt(dicepp ON "v2.6.1")
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
  include_netcdf_cxx4(dicepp ON v4.3.0)
  target_compile_definitions(dicepp PRIVATE DICEPP_WITH_NETCDF)
endif()


include(lib/settingsnode/settingsnode.cmake)
include_settingsnode(dicepp)
include_yaml_cpp(dicepp ON "yaml-cpp-0.6.3")

add_cpp_tools(dicepp STD c++14)
