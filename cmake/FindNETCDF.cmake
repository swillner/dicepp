# adapted from https://github.com/esa/pagmo2/blob/master/cmake_modules/FindNLOPT.cmake

if(NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
    set(NETCDF_FIND_QUIETLY TRUE)
endif()

find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h)
find_library(NETCDF_LIBRARY NAMES netcdf)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(NETCDF DEFAULT_MSG NETCDF_INCLUDE_DIR NETCDF_LIBRARY)

mark_as_advanced(NETCDF_INCLUDE_DIR NETCDF_LIBRARY)

if(NETCDF_FOUND AND NOT TARGET NETCDF::netcdf)
    add_library(NETCDF::netcdf UNKNOWN IMPORTED)
    set_target_properties(NETCDF::netcdf PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_INCLUDE_DIR}")
    set_target_properties(NETCDF::netcdf PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${NETCDF_LIBRARY}")
endif()
