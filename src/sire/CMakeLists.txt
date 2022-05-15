################################
#
#  CMake Project file for sire
#
################################

cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)
cmake_policy(VERSION 3.3)

## Read the version number (so we can set it once in one place)
file (STRINGS "../../version.txt" SIRE_VERSION)

string(REGEX MATCH "([0-9]*).([0-9]*).([0-9]*)" _ ${SIRE_VERSION})

set (S_VERSION_MAJOR ${CMAKE_MATCH_1})
set (S_VERSION_MINOR ${CMAKE_MATCH_2})
set (S_VERSION_PATCH ${CMAKE_MATCH_3})

message( STATUS "Building Sire Python Module ${SIRE_VERSION}")

# use loose loop syntax in CMakeLists files
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

message( STATUS "Building into an existing conda installation..." )
# Where is anaconda installed
set( ANACONDA_BASE "$ENV{HOME}/miniconda3" CACHE PATH
        "Install directory for conda (python3)" )

# Make sure that there is a python3 executable in the bin directory
if ( NOT EXISTS "${ANACONDA_BASE}/bin/python3" AND
      NOT EXISTS "${ANACONDA_BASE}/python.exe"  AND
      NOT EXISTS "${PYTHON_EXECUTABLE}" )
  message( FATAL_ERROR "Cannot find ${ANACONDA_BASE}/bin/python3, ${ANACONDA_BASE}/python.exe "
                        "or ${PYTHON_EXECUTABLE}."
                        "Please set the full path to the conda "
                        "installation directory." )
endif()

set( SIRE_APP "${ANACONDA_BASE}/pkgs/sire-${SIRE_VERSION}" CACHE PATH
        "Install directory for Sire within the conda distribution" )

message( STATUS "Using SIRE_APP = ${SIRE_APP}" )

set( SIRE_INCLUDE_DIR "${SIRE_APP}/include/Sire" )

set( CONDA_INCLUDE_DIR "${ANACONDA_BASE}/include" )
message( STATUS "Using CONDA_INCLUDE_DIR = ${CONDA_INCLUDE_DIR}" )

set ( SIRE_COMP_CMAKE "${SIRE_INCLUDE_DIR}/cmake/SireCompileVariables.cmake" )

# Look for the "SireCompileVariables.cmake" file in the install directory
if (EXISTS "${SIRE_COMP_CMAKE}")
  include( "${SIRE_COMP_CMAKE}" )
else()
  message( FATAL_ERROR "You don't appear to have built corelib into this "
              "Anaconda/Miniconda directory. Please build corelib first. "
              "(there should be a file called ${SIRE_COMP_CMAKE})" )
endif()

# name the project (no language support as just Python modules)
project (SireModule NONE)

message( STATUS "PROJECT COMPLETED" )

# Absolutely can't run cmake in the source directory!
if ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
  message( FATAL_ERROR "You must run CMake in a different directory to the source!" )
endif ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )

if (NOT SIRE_FOUND_PYTHON)
  # Need to find python before we compile the boost python wrappers
  find_package( PythonLibs 3.7 REQUIRED )
endif()

message(STATUS "SIRE_INSTALL_PREFIX ${SIRE_INSTALL_PREFIX}")

set(CMAKE_INSTALL_PREFIX "${SIRE_INSTALL_PREFIX}" CACHE INTERNAL
         "Prefix prepended to install directories" FORCE )

# set the location to install the python modules
include( "build/cmake/install_python.cmake" )
set (SIRE_PYTHON "${PYTHON_SITE_DIR}")

message( STATUS "Installing Python wrappers to ${SIRE_PYTHON}")

# Recurse into the module subdirectories
add_subdirectory (config)
add_subdirectory (utils)
add_subdirectory (qt)
add_subdirectory (error)
add_subdirectory (analysis)
add_subdirectory (base)
add_subdirectory (cas)
add_subdirectory (cluster)
add_subdirectory (ff)
add_subdirectory (id)
add_subdirectory (io)
add_subdirectory (maths)
add_subdirectory (mm)
add_subdirectory (mol)
add_subdirectory (move)
add_subdirectory (search)
add_subdirectory (stream)
add_subdirectory (squire)
add_subdirectory (system)
add_subdirectory (units)
add_subdirectory (vol)

install( FILES __init__.py _load.py _thumbsup.py _pythonize.py
         DESTINATION ${SIRE_PYTHON}/sire
       )

# These commands are used to create an 'uninstall' target
# (this is copied from the vtk.org wiki)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
