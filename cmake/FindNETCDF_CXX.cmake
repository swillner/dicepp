# adapted from https://github.com/esa/pagmo2/blob/master/cmake_modules/FindNLOPT.cmake

if(NETCDF_CXX_INCLUDE_DIR AND NETCDF_CXX_LIBRARY)
    set(NETCDF_CXX_FIND_QUIETLY TRUE)
endif()

find_path(NETCDF_CXX_INCLUDE_DIR NAMES ncFile.h)
find_library(NETCDF_CXX_LIBRARY NAMES netcdf_c++4)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(NETCDF_CXX DEFAULT_MSG NETCDF_CXX_INCLUDE_DIR NETCDF_CXX_LIBRARY)

mark_as_advanced(NETCDF_CXX_INCLUDE_DIR NETCDF_CXX_LIBRARY)

# NOTE: this has been adapted from CMake's FindPNG.cmake.
if(NETCDF_CXX_FOUND AND NOT TARGET NETCDF_CXX::netcdf_cxx)
    add_library(NETCDF_CXX::netcdf_cxx UNKNOWN IMPORTED)
    set_target_properties(NETCDF_CXX::netcdf_cxx PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_CXX_INCLUDE_DIR}")
    set_target_properties(NETCDF_CXX::netcdf_cxx PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${NETCDF_CXX_LIBRARY}")
endif()
