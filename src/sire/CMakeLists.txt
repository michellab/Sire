################################
#
#  CMake Project file for sire
#
################################

# project version
set (SIRE_VERSION_MAJOR "2022")
set (SIRE_VERSION_MINOR "2")
set (SIRE_VERSION_PATCH "0")

set (SIRE_VERSION "${SIRE_VERSION_MAJOR}.${SIRE_VERSION_MINOR}.${SIRE_VERSION_PATCH}")

cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)

# use loose loop syntax in CMakeLists files
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

# Option to compile and install into an anaconda or miniconda directory
option( ANACONDA_BUILD "Build into an existing anaconda/miniconda installation" OFF )

if ( ANACONDA_BUILD )
  message( STATUS "Building into an existing Anaconda/Miniconda installation..." )
  # Where is anaconda installed
  set( ANACONDA_BASE "$ENV{HOME}/miniconda3" CACHE PATH
          "Install directory for Anaconda/Miniconda (python3)" )

  # Make sure that there is a python3 executable in the bin directory
  if ( NOT EXISTS "${ANACONDA_BASE}/bin/python3" AND
       NOT EXISTS "${ANACONDA_BASE}/python.exe"  AND
       NOT EXISTS "${PYTHON_EXECUTABLE}" )
    message( FATAL_ERROR "Cannot find ${ANACONDA_BASE}/bin/python3, ${ANACONDA_BASE}/python.exe "
                         "or ${PYTHON_EXECUTABLE}."
                         "Please set the full path to the Anaconda or Miniconda "
                         "installation directory." )
  endif()

  set( SIRE_APP "${ANACONDA_BASE}/pkgs/sire-${SIRE_VERSION}" CACHE PATH
          "Install directory for Sire within the Anaconda/Miniconda distribution" )
#elseif( WIN32 OR WIN64 )
#  # This should be MSYS2/MINGW64, so use the location of cmake to find the location
#  # of the Sire installation
#  message(STATUS "Compiling on windows using MSYS2/MINGW64...")
#  set( MINGW_BASE "C:/msys64/mingw64" CACHE PATH "MSYS MINGW64 root directory" )
#  get_filename_component(MINGW_CMAKE_BIN ${CMAKE_COMMAND} DIRECTORY)
#  get_filename_component(MINGW_CMAKE_ROOT ${MINGW_CMAKE_BIN} DIRECTORY)
#  message(STATUS "Using cmake ${CMAKE_COMMAND} from MING64 base ${MINGW_CMAKE_ROOT}")
#
#  set( MINGW_BASE "${MINGW_CMAKE_ROOT}")
#  set( SIRE_APP "${MINGW_BASE}" )

else()
  #�Option to set the install directory of Sire
  if (NOT DEFINED SIRE_APP)
    set( SIRE_APP "$ENV{HOME}/sire.app" CACHE PATH
            "The full path to the sire.app installed Sire application." )
  endif()
endif()

message( STATUS "Using SIRE_APP = ${SIRE_APP}" )

set( SIRE_INCLUDE_DIR "${SIRE_APP}/include/Sire" )

set( CONDA_INCLUDE_DIR "${ANACONDA_BASE}/include" )
message( STATUS "Using CONDA_INCLUDE_DIR = ${CONDA_INCLUDE_DIR}" )

set ( SIRE_COMP_CMAKE "${SIRE_INCLUDE_DIR}/cmake/SireCompileVariables.cmake" )

# Look for the "SireCompileVariables.cmake" file in the install directory
if (EXISTS "${SIRE_COMP_CMAKE}")
  include( "${SIRE_COMP_CMAKE}" )
else()
  if (ANACONDA_BUILD)
    message( FATAL_ERROR "You don't appear to have built corelib into this "
               "Anaconda/Miniconda directory. Please build corelib first. "
               "(there should be a file called ${SIRE_COMP_CMAKE})" )
  else()
    message( STATUS "You need to provide the full path to sire.app and check that it contains the file ${SIRE_COMP_CMAKE}" )
    message( FATAL_ERROR "Do this by running 'cmake -DSIRE_APP=/path/to/sire.app /path/to/python/source'" )
  endif()
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
find_package( PythonLibs 3.2 REQUIRED )
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
