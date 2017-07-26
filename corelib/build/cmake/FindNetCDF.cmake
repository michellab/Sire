# Find NetCDF library.
#
# Looks for the NetCDF libraries at the default (/usr/local) location 
# or custom location found in the NETCDF_ROOT_DIR environment variable. 
#
# The script defines defines: 
#  NetCDF_FOUND
#  NetCDF_ROOT_DIR
#  NetCDF_INCLUDE_DIR
#  NetCDF_LIBRARY_DIR
#  NetCDF_LIBRARIES      
#

if(NetCDF_INCLUDE_DIR AND NetCDF_LIBRARY_DIR)
    set(NetCDF_FIND_QUIETLY)
endif()

file(TO_CMAKE_PATH "$ENV{NetCDF_ROOT_DIR}" _env_NEFCDF_ROOT_DIR)

find_library(NetCDF_LIBRARIES
    NAMES netcdf
    PATHS "${NetCDF_ROOT_DIR}/lib"
    CACHE STRING "NetCDF libraries")

get_filename_component(NetCDF_LIBRARY_DIR 
    ${NetCDF_LIBRARIES} 
    PATH
    CACHE STRING "NetCDF library path")

find_path(NetCDF_INCLUDE_DIR 
    NAMES netcdf.h 
    PATHS "${NetCDF_ROOT_DIR}/include" "${NetCDF_LIBRARY_DIR}/../include"
    CACHE STRING "NetCDF include directory")    

# if we did not manage to set the root dir at the beginning but found the 
# libs then set the ${NetCDF_LIBRARY_DIR}/.. as root
if(NOT IS_DIRECTORY ${NetCDF_ROOT_DIR})
    if (IS_DIRECTORY "${NetCDF_LIBRARY_DIR}/..") # just double-checking
        get_filename_component(NetCDF_ROOT_DIR 
            "${NetCDF_LIBRARY_DIR}/.." 
            ABSOLUTE)
    endif()   
endif()

if(NOT IS_DIRECTORY ${NetCDF_ROOT_DIR})
    message(STATUS "Could not find NetCDF! Set the NETCDF_ROOT_DIR environment "
    "variable to contain the path of the NetCDF installation if you want to use it.")
endif()

if(NOT IS_DIRECTORY ${NetCDF_LIBRARY_DIR})
    message(STATUS "Could not find NetCDF! Set the NETCDF_ROOT_DIR environment "
    "variable to contain the path of the NetCDF installation if you want to use it.")
    set (NetCDF_FOUND 0)
else()
    set (NetCDF_FOUND 1)
endif()

if ( ${NetCDF_FOUND} )
  if(NOT NetCDF_INCLUDE_DIR)
      message(FATAL_ERROR "Can't find NetCDF includes. Check your NetCDF installation!")
  endif()

  set(NetCDF_ROOT_DIR ${NetCDF_ROOT_DIR} CACHE PATH "NetCDF installation directory")

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(NetCDF DEFAULT_MSG 
                                      NetCDF_ROOT_DIR
                                      NetCDF_LIBRARIES 
                                      NetCDF_LIBRARY_DIR 
                                      NetCDF_INCLUDE_DIR)
 
  mark_as_advanced(NetCDF_INCLUDE_DIR
                   NetCDF_LIBRARIES
                   NetCDF_LIBRARY_DIR)
endif()

