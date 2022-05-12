################################
#
#  CMake Project file for Sire
#    Python 3 wrappers
#
################################

# require cmake >= 3.3.0 (released in 2018)
cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)
cmake_policy(VERSION 3.3)

# Read the version number (so we can set it once in one place)
file (STRINGS "../version.txt" SIRE_VERSION)

string(REGEX MATCH "([0-9]*).([0-9]*).([0-9]*)" _ ${SIRE_VERSION})

set (S_VERSION_MAJOR ${CMAKE_MATCH_1})
set (S_VERSION_MINOR ${CMAKE_MATCH_2})
set (S_VERSION_PATCH ${CMAKE_MATCH_3})

message( STATUS "Building Sire Wrappers ${SIRE_VERSION}")

############
############ SETTING INITIAL SIRE VARIABLES
############

# use loose loop syntax in CMakeLists files
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

# Get rid of install messages (too verbose!)
set( CMAKE_INSTALL_MESSAGE "NEVER" )

# Absolutely can't run cmake in the source directory!
if ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
  message( FATAL_ERROR "You must run CMake in a different directory to the source!" )
endif ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )

if ( CMAKE_GENERATOR MATCHES "Xcode" )
  message( FATAL_ERROR "You cannot generate the Xcode project for the Python "
            "wrappers. These are auto-generated and should not be hand-edited." )
endif()

# get the number of cores on this machine - useful for
# compiling dependencies
include(ProcessorCount)
ProcessorCount(NCORES)
if (NCORES EQUAL 0)
    set( NCORES 1 )
endif()
message(STATUS "Configuring on a machine with number of cores = ${NCORES}")

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
                       "Please set the full path to the Anaconda or Miniconda "
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
  message( FATAL_ERROR "You don't appear to have built corelib into this"
             "conda directory. Please build corelib first."
             "(there should be a file called ${SIRE_COMP_CMAKE})" )
endif()

# Now set the C++ and C compilers to be used to compile the wrappers
# These have to be the same as used to compile corelib
set(CMAKE_C_COMPILER "${SIRE_C_COMPILER}")
set(CMAKE_CXX_COMPILER "${SIRE_CXX_COMPILER}")

message( STATUS "Setting C compiler to ${SIRE_C_COMPILER} and C++ compiler to ${SIRE_CXX_COMPILER}" )

# conda can set incorrect CXX flags (i.e. -mtune=haswell on MacOS M1!)
# If is safer to completely unset these flags
unset(ENV{CXXFLAGS})
unset(ENV{DEBUG_CXXFLAGS})
unset(ENV{CFLAGS})
unset(ENV{DEBUG_CFLAGS})

# name the project
project (SirePython CXX C)

# Absolutely can't run cmake in the source directory!
if ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
  message( FATAL_ERROR "You must run CMake in a different directory to the source!" )
endif ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )

set(CMAKE_INSTALL_PREFIX "${SIRE_INSTALL_PREFIX}" CACHE INTERNAL
         "Prefix prepended to install directories" FORCE )

message( STATUS "Sire will be compiled and installed to directory ${CMAKE_INSTALL_PREFIX}" )

# We will be managing RPATH ourselves
set( CMAKE_SKIP_RPATH FALSE )
set (CMAKE_SKIP_BUILD_RPATH FALSE)

if (APPLE)
  # make sure that all libraries are named "@rpath/library.dylib". This will allow
  # the library to be found as long as the executable can substitute @rpath with the
  # path to the library
  set (CMAKE_INSTALL_NAME_DIR "@rpath")

  # Add the two passed directories to the @rpaths to search for libraries. @executable_path
  # is the directory containing the executable, so the libraries are included in
  # @exeutable_path/../${SIRE_LIBS}
  set (SIRE_INSTALL_RPATH "@executable_path/../${SIRE_LIBS}")
else()
  set (SIRE_INSTALL_RPATH "\$ORIGIN/../${SIRE_LIBS}:\$ORIGIN/../../../../../${SIRE_LIBS}")
endif()

set (CMAKE_INSTALL_RPATH "${SIRE_INSTALL_RPATH}")
set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# Path to the project's extra cmake files
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/build/cmake)

# Find the Sire libraries
find_package (Sire 0.0.1 REQUIRED)

if (NOT SIRE_FOUND)
   message( FATAL_ERROR "The required version of Sire was not found.")
endif()

message( STATUS "Found Sire version ${SIRE_VERSION}" )
message( STATUS "Sire root = ${SIRE_INSTALL_PREFIX}" )
message( STATUS "Header files in ${SIRE_INCLUDE_DIR}" )
message( STATUS "Libraries in ${SIRE_INSTALL_PREFIX}/${SIRE_LIBS}" )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_strip.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip"
  IMMEDIATE @ONLY)

set ( CMAKE_STRIP "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip" )
set ( SIRE_STRIP_COMMAND "${CMAKE_STRIP} ${SIRE_STRIP_OPTIONS}" )

# Make sure that we are using the same compiler as that used to compile Sire
option (SIRE_CHECK_COMPILER "Turn off to skip the compiler check" ON)

macro( compare_strings RESULT val0 val1)
  if (${val0} STREQUAL ${val1})
    set( ${RESULT} TRUE )
  else()
    set( ${RESULT} FALSE )
  endif()
endmacro()

macro( compare_values RESULT val0 val1)
  if (${val0} EQUAL ${val1})
    set( ${RESULT} TRUE )
  else()
    set( ${RESULT} FALSE )
  endif()
endmacro()

if ( SIRE_CHECK_COMPILER )

  compare_strings( SAME_C_COMPILER ${SIRE_C_COMPILER} ${CMAKE_C_COMPILER} )
  compare_strings( SAME_CXX_COMPILER ${SIRE_CXX_COMPILER} ${CMAKE_CXX_COMPILER} )

  if ( NOT SAME_C_COMPILER OR NOT SAME_CXX_COMPILER )

    if ( NOT SAME_C_COMPILER )
      message( STATUS "Different C compilers! ${SIRE_C_COMPILER} vs. ${CMAKE_C_COMPILER}" )
    endif()

    if ( NOT SAME_CXX_COMPILER )
      message( STATUS "Different C++ compilers! ${SIRE_CXX_COMPILER} vs. ${CMAKE_CXX_COMPILER}" )
    endif()

    message( STATUS "Using different compilers is really not recommended. If you want to "
                    "skip this check, then set SIRE_CHECK_COMPILER to \"OFF\" using ccmake" )

    fatal_message( "You cannot compile the python wrappers with different compilers"
                   "(${CMAKE_C_COMPILER}, ${CMAKE_CXX_COMPILER}) than those used to"
                   "compile the Sire libraries"
                   "(${SIRE_C_COMPILER}, ${SIRE_CXX_COMPILER})."
                   "You need to remove CMakeCache.txt and rerun cmake using the command;\n\n"
                   "CC=${SIRE_C_COMPILER} CXX=${SIRE_CXX_COMPILER} cmake ${CMAKE_SOURCE_DIR}\n" )
  endif()
endif()

# Now add the paths to the header files of the dependencies of Sire
include_directories (BEFORE ${SIRE_INCLUDE_DIR})
include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})

if (TBB_INCLUDE_DIR)
  include_directories (${TBB_INCLUDE_DIR})
  message( STATUS "Using TBB from ${TBB_INCLUDE_DIR}")
endif()

if (SIRE_GSL_INCLUDE_PATH)
  include_directories (${SIRE_GSL_INCLUDE_PATH})
endif()

if(CMAKE_GENERATOR MATCHES "Visual Studio")  # MSBuild
  # Everything we need to do on Windows to allow the wrappers to
  # build using Visual Studio

  # Thanks to this project who have a solution to the tbb12 link
  # issue on Windows
  # https://fossies.org/linux/opencv/cmake/OpenCVFindLAPACK.cmake

  # We may have tbb.lib for 'tbb' target, but not 'tbb12.lib'
  message(STATUS "Adding tbb12 to ignore list")
  set(SIRE_IGNORE_DEFAULTLIB_TBB "tbb12.lib")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIB:${SIRE_IGNORE_DEFAULTLIB_TBB}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:${SIRE_IGNORE_DEFAULTLIB_TBB}")

  message(STATUS "Turning on bigobj")  # needed to link _Mol.pyd
  add_compile_options("/bigobj")
endif()

# Conda can mess up our include paths. This makes sure that every
# directory is included directly, including testing that the
# CONDA include directory is working
include_directories( "${SIRE_APP}/include" "${SIRE_APP}/Library/include" )
include_directories(BEFORE ${CONDA_INCLUDE_DIR})

message(STATUS "Checking if conda/include is a hidden include...")
include(CheckCXXSourceRuns)
get_property(CMAKE_REQUIRED_INCLUDES DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "${CMAKE_REQUIRED_INCLUDES}")
  check_cxx_source_runs("

  #include \"tbb/tbb.h\"

  int main()
  {
    return 0;
  }" IS_CONDA_INCLUDE_WORKING)

if (IS_CONDA_INCLUDE_WORKING)
  message(STATUS "-Iconda/include is used correctly.")
else()
  # Things don't work, so we need to switch off implicit include directories
  unset(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES)
  message(STATUS "Unsetting CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES so that conda/include is not hidden.")
endif()

# The wrappers MUST be built as shared libraries
set (BUILD_SHARED_LIBS ON)

# Uncomment for more verbose compiling/linking
# (or run 'make VERBOSE=1')
# set (CMAKE_VERBOSE_MAKEFILE ON)

# We will use the small compiler flags (as the wrappers are too big and don't need optimising)
message( STATUS "Compiling the python wrappers using a SMALL build." )
set( CMAKE_C_FLAGS ${SIRE_C_FLAGS_SMALL} )
set( CMAKE_CXX_FLAGS ${SIRE_CXX_FLAGS_SMALL} )

if (SIRE_HAS_CPP_LIB)
  message(STATUS "Using libc++ in place of libstdc++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${SIRE_SHARE_LINK_FLAGS}" )
set( CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${SIRE_STATIC_LINK_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${SIRE_EXE_LINK_FLAGS}" )

include (CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG( "-w" SIRE_HAS_NO_WARNINGS_FLAG )

if (SIRE_HAS_NO_WARNINGS_FLAG)
  message(STATUS "Inhibiting all compiler warnings, as the wrappers generate a lot!")
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w" )
else()
  CHECK_CXX_COMPILER_FLAG( "-Wno-unused-variable" SIRE_HAS_NO_UNUSED_VARIABLE )
  CHECK_CXX_COMPILER_FLAG( "-Wno-overloaded-virtual" SIRE_HAS_NO_OVERLOADED_VIRTUAL )

  if (SIRE_HAS_NO_UNUSED_VARIABLE)
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
  endif()

  if (SIRE_HAS_NO_OVERLOADED_VIRTUAL)
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-overloaded-virtual")
  endif()
endif()

if (APPLE)
  # Also on MacOS need to add '-undefined dynamic_lookup' so that we can build
  # modules that can be imported by a python that has statically linked libPython
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup" )
endif()

message( STATUS  "CMAKE_SYSTEM_NAME      == ${CMAKE_SYSTEM_NAME}" )
message( STATUS  "CMAKE_C_COMPILER       == ${CMAKE_C_COMPILER}" )
message( STATUS  "CMAKE_CXX_COMPILER     == ${CMAKE_CXX_COMPILER}" )
message( STATUS  "C compiler flags       == ${CMAKE_C_FLAGS}" )
message( STATUS  "C++ compiler flags     == ${CMAKE_CXX_FLAGS}" )
message( STATUS  "Shared library flags   == ${CMAKE_SHARED_LINKER_FLAGS}" )
message( STATUS  "Static library flags   == ${CMAKE_STATIC_LINKER_FLAGS}" )
message( STATUS  "Executable link flags  == ${CMAKE_EXE_LINKER_FLAGS}" )

# Now lets find Python...
find_package( PythonLibs 3.7 REQUIRED )

if (NOT PYTHON_MODULE_EXTENSION)
  if (WIN32)
    set( PYTHON_MODULE_EXTENSION ".pyd" )
  else()
    set( PYTHON_MODULE_EXTENSION ".so" )
  endif()
endif()

if (NOT PYTHON_LIBRARIES)
  set( PYTHON_LIBRARIES "${PYTHON_LIBRARY}" )
endif()

find_package( PythonInterp REQUIRED )  # needed to get the python version
set( PYTHON_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}" )

# Where we should install our packages
if (WIN32)
  set( PYTHON_SITE_DIR "../../lib/site-packages" )
else()
  set( PYTHON_SITE_DIR "../../lib/python${PYTHON_VERSION}/site-packages" )
endif()

# Don't link against libPython on macOS! Since the conda-forge Python interpreter
# statically links libPython you'll end up with duplicate symbols.
if (APPLE)
  set( PYTHON_LIBRARIES "" )
endif()

# All boost python libraries must be dynamically linked
# (else we could end up with multiple boost::python libraries linked, which
#  causes problems for RTTI, implicit converters and exceptions)
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_STATIC_RUNTIME OFF)
set (Boost_DEBUG ON)  # this gives a nice print out in the cmake list to debug how boost is found and linked
find_package( Boost 1.7 COMPONENTS python REQUIRED )

message( STATUS "Python paths ${PYTHON_LIBRARIES} | ${PYTHON_INCLUDE_DIR} | ${PYTHON_SITE_DIR}" )
include_directories( ${PYTHON_INCLUDE_DIR} )

message( STATUS "boost::python library ${Boost_LIBRARIES}" )

# set the location to install the python modules
set (SIRE_PYTHON "${PYTHON_SITE_DIR}")

# Use the same Qt library as was used to build corelib - this can get
# unset a lot, hence why we keep checking it and updating

# look for Qt5 in the system - this should be found from conda
find_package( Qt5Core REQUIRED )

get_target_property(QT5_LOCATION Qt5::Core LOCATION)
get_filename_component(SIRE_QT5CORE_DIR ${QT5_LOCATION} PATH)

if (QT5_LOCATION MATCHES NOTFOUND)
    message( STATUS "Looking for Qt5Core in directory ${SIRE_QT5CORE_DIR}" )
    set( Qt5Core_DIR "${SIRE_QT5CORE_DIR}/cmake/Qt5Core" )
    find_package( Qt5Core REQUIRED )
endif()

get_target_property(QT5_LOCATION Qt5::Core LOCATION)
message( STATUS "Using Qt5Core in ${QT5_LOCATION}" )

if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
  message( STATUS "Adding useful MSVC compiler options" )
  add_compile_options("/EHsc")  # exceptions are only thrown from C++ throw calls
  add_compile_options("/permissive-")  # disable MSVC's permissive mode - ensures code is standards conforming
endif()

if (NOT DEFINED ${SIRE_VERSION})
  message( STATUS "Strange, SIRE_VERSION has been unset. Resetting it.")
  set (SIRE_VERSION "${S_VERSION_MAJOR}.${S_VERSION_MINOR}.${S_VERSION_PATCH}")
endif()

message( STATUS "SIRE VERSION IS ${SIRE_VERSION}")

# Add the helper library
add_subdirectory (Helpers)

# Recurse into the module subdirectories
add_subdirectory (Qt)
add_subdirectory (Error)
add_subdirectory (Analysis)
add_subdirectory (Base)
add_subdirectory (CAS)
add_subdirectory (Cluster)
add_subdirectory (FF)
add_subdirectory (ID)
add_subdirectory (IO)
add_subdirectory (Maths)
add_subdirectory (MM)
add_subdirectory (Mol)
add_subdirectory (Move)
add_subdirectory (Search)
add_subdirectory (Stream)
add_subdirectory (Squire)
add_subdirectory (System)
add_subdirectory (Units)
add_subdirectory (Vol)

# also compile Sire's python executable
add_subdirectory (python)

# and some additional tools used for e.g. waterswap or nautillus
add_subdirectory (Tools)

# Need to install the dunder for Sire. All C++ wrappers are
# now placed into the sire/legacy folder.
install( FILES __init__.py
         DESTINATION ${SIRE_PYTHON}/sire/legacy
       )

# Finally, write the version header file and finish off the component dependencies file
include( GetGitVersion )

# Create the sire_python_config.h to tell Sire where the
# python libraries will be installed
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/sire_python_config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/sire_python_config.h"
  IMMEDIATE @ONLY)

# These commands are used to create an 'uninstall' target
# (this is copied from the vtk.org wiki)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

if (WIN32)
  # OpenMM.dll is installed in conda/Library/lib when it should be
  # in conda/Library/bin. Check if it is missing and then copy
  # the file across
  file( GLOB files "${ANACONDA_BASE}/Library/lib/OpenMM*.dll" )

  foreach(file ${files})
    get_filename_component(OMMDLL, "${file}")
    set( BINDLL "${ANACONDA_BASE}/Library/bin/${OMMDLL}" )
    if (NOT EXISTS "${BINDLL}")
      message( STATUS "OpenMM DLL in the wrong place? Should be ${BINDLL} but is ${OMMDLL}" )
    endif()
  endforeach()
endif()
