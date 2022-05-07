################################
#
#  CMake Project file for Sire
#    Python 3 wrappers
#
#   (C) Christopher Woods
#
################################

# project version
set (SIRE_VERSION_MAJOR "2022")
set (SIRE_VERSION_MINOR "2")
set (SIRE_VERSION_PATCH "0")

set (SIRE_VERSION "${SIRE_VERSION_MAJOR}.${SIRE_VERSION_MINOR}.${SIRE_VERSION_PATCH}")

############
############ SETTING INITIAL SIRE VARIABLES
############

# require cmake >= 2.8.11.2
cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)
cmake_policy(VERSION 2.8)

# use loose loop syntax in CMakeLists files
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

if(COMMAND CMAKE_POLICY)
  if (POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)  # Anaconda clang is called Clang, not AppleClang
  endif()
  if (POLICY CMP0042)
    cmake_policy(SET CMP0042 OLD)  # managing OS X @rpath manually
  endif()
  if (POLICY CMP0048)
    cmake_policy(SET CMP0048 OLD)  # want to manage versions myself!
  endif()
  if (POLICY CMP0028)
    cmake_policy(SET CMP0028 OLD)  # Qt has colons in names, which new cmake hates
  endif()
  if (POLICY CMP0045)
    cmake_policy(SET CMP0045 OLD)  # Allow non-existant targets with -NOTFOUND
  endif()
  if (POLICY CMP0020)
    cmake_policy(SET CMP0020 NEW)  # Automatically link to qtmain on windows
  endif()
endif()

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

set ( SIRE_COMP_CMAKE "${SIRE_INCLUDE_DIR}/cmake/SireCompileVariables.cmake" )

# Look for the "SireCompileVariables.cmake" file in the install directory
if (EXISTS "${SIRE_COMP_CMAKE}")
  include( "${SIRE_COMP_CMAKE}" )
else()
  if (ANACONDA_BUILD)
    message( FATAL_ERROR "You don't appear to have built corelib into this"
               "Anaconda/Miniconda directory. Please build corelib first."
               "(there should be a file called ${SIRE_COMP_CMAKE})" )
  else()
    message( STATUS "You need to provide the full path to sire.app and check that it contains the file ${SIRE_COMP_CMAKE}" )
    message( FATAL_ERROR "Do this by running 'cmake -DSIRE_APP=/path/to/sire.app /path/to/python/source'" )
  endif()
endif()

if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
  set (SLASH_HYPHEN "/")
else()
  set (SLASH_HYPHEN "-")
endif()

# Now set the C++, C and Fortran compilers to be used to compile the wrappers
# These have to be the same as used to compile corelib
set(CMAKE_C_COMPILER "${SIRE_C_COMPILER}")
set(CMAKE_CXX_COMPILER "${SIRE_CXX_COMPILER}")

message( STATUS "Setting C compiler to ${SIRE_C_COMPILER} and C++ compiler to ${SIRE_CXX_COMPILER}" )

# name the project
project (SirePython CXX C)

message( STATUS "PROJECT COMPLETED" )

set(CMAKE_INSTALL_PREFIX "${SIRE_INSTALL_PREFIX}" CACHE INTERNAL
         "Prefix prepended to install directories" FORCE )

message( STATUS "Sire will be compiled and installed to directory ${CMAKE_INSTALL_PREFIX}" )

set( CMAKE_SKIP_RPATH FALSE )
set (CMAKE_SKIP_BUILD_RPATH FALSE)

if (APPLE)
  # make sure that all libraries are named "@rpath/library.dylib". This will allow
  # the library to be found as long as the executable can substitute @rpath with the
  # path to the library
  set (CMAKE_INSTALL_NAME_DIR "@rpath")

  # Add the two passed directories to the @rpaths to search for libraries. @executable_path
  # is the directory containing the executable, so the libraries are included in
  # @exeutable_path/../${SIRE_LIBS} and @exeutable_path/../${SIRE_BUNDLED_LIBS}
  set (SIRE_INSTALL_RPATH "@executable_path/../${SIRE_LIBS}:@executable_path/../${SIRE_BUNDLED_LIBS}")
else()
  set (SIRE_INSTALL_RPATH "\$ORIGIN/../${SIRE_LIBS}:\$ORIGIN/../${SIRE_BUNDLED_LIBS}:\$ORIGIN/../../../../../${SIRE_LIBS}:\$ORIGIN/../../../../../${SIRE_BUNDLED_LIBS}")
endif()

set (CMAKE_INSTALL_RPATH "${SIRE_INSTALL_RPATH}")
set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# Absolutely can't run cmake in the source directory!
if ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
  message( FATAL_ERROR "You must run CMake in a different directory to the source!" )
endif ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )

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

set ( SIRE_STRIP_COMMAND "${CMAKE_STRIP} ${SIRE_STRIP_OPTIONS}" )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_strip.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip"
  IMMEDIATE @ONLY)

set ( CMAKE_STRIP "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip" )

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
set (TBB_INCLUDE_DIR "/usr/include/tbb")
if (TBB_INCLUDE_DIR)
  include_directories (${TBB_INCLUDE_DIR})
  message( STATUS "Using TBB from ${TBB_INCLUDE_DIR}")
endif()
if (SIRE_GSL_INCLUDE_PATH)
  include_directories (${SIRE_GSL_INCLUDE_PATH})
endif()
if (DEFINED SIRE_APP)
  # For some reason, when using the conda GNU compilers,  CMake stripped
  # "${SIRE_APP}/include" from the directories searched for header files
  # If this happens, uncomment the following line
  # unset(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES)
  include_directories( "${SIRE_APP}/include" "${SIRE_APP}/Library/include" )
endif()

# MUST build shared libraries
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
  # Also need to add '-undefined dynamic_lookup' so that we can build
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

# Option to disable all bundles
option( SIRE_PYTHON_ENABLE_BUNDLING "Allow Sire Python to use bundled libraries and python interpreter" ON )

# Fine-grain control of bundling each dependency
option( SIRE_BUNDLE_PYTHON "Whether Sire should use its bundled copy of Python 3" ON )
option( SIRE_BUNDLE_BOOST_PYTHON "Whether Sire should use its bundled copy of boost::python" ON )

if ( SIRE_PYTHON_ENABLE_BUNDLING )
  if (MSYS)
    message( STATUS "Building a bundled msys Sire (using msys python)" )
  else()
    message( STATUS "Sire Python will be compiled to use its own copy of bundled libraries." )
  endif()

  set( BUNDLE_BUILDDIR "${CMAKE_BINARY_DIR}/build_bundled" )

  # We want to install the bundled libs directly into the sire.app directory
  set( BUNDLE_STAGEDIR "${SIRE_APP}/${SIRE_BUNDLED_DIR}" )

  if (NOT EXISTS "${BUNDLE_BUILDDIR}")
    file(MAKE_DIRECTORY ${BUNDLE_BUILDDIR})
  endif()

  if (NOT EXISTS "${BUNDLE_STAGEDIR}")
   file(MAKE_DIRECTORY ${BUNDLE_STAGEDIR})
  endif()

  if (NOT EXISTS "${BUNDLE_STAGEDIR}/lib")
    file(MAKE_DIRECTORY ${BUNDLE_STAGEDIR}/lib)
  endif()

  if (NOT EXISTS "${BUNDLE_STAGEDIR}/include")
    file(MAKE_DIRECTORY ${BUNDLE_STAGEDIR}/include)
  endif()

  if ( SIRE_BUNDLE_PYTHON )
    include( "bundled/install_python.cmake" )
  endif()

  if (NOT SIRE_FOUND_PYTHON)
    # Need to find python before we compile the boost python wrappers
    find_package( PythonLibs 3.2 REQUIRED )
  endif()

  if ( SIRE_BUNDLE_BOOST_PYTHON )
    include( "bundled/install_boost_python.cmake" )
  endif()

  if ( NOT ${BUNDLE_STAGEDIR} STREQUAL ${CMAKE_INSTALL_PREFIX} )
    # All of the bundled objects are built into the directory ${CMAKE_BINARY_DIR}/bundled
    # This directory must be installed into the sire.app/ directory
    install( DIRECTORY ${BUNDLE_STAGEDIR}
             DESTINATION ${CMAKE_INSTALL_PREFIX}
             USE_SOURCE_PERMISSIONS )
  endif()
endif()

if (NOT ${SIRE_FOUND_PYTHON})
  find_package( PythonLibs 3.2 REQUIRED )
endif()

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

# Don't link against libPython on macOS! Since the conda-forge Python interpreter
# statically links libPython you'll end up with duplicate symbols.
if (APPLE)
  set( PYTHON_LIBRARIES "" )
endif()

if (NOT ${SIRE_FOUND_BOOST_PYTHON})
  # Sire depends on >= boost 1.31 - all libraries must be dynamically linked
  # (else we could end up with multiple boost::python libraries linked, which
  #  causes problems for RTTI and exceptions)
  set (BOOST_ALL_DYN_LINK "YES")
  find_package( Boost 1.31 COMPONENTS python REQUIRED )
  set ( BOOST_PYTHON_LIBRARY "${Boost_PYTHON_LIBRARY}" )
  set ( BOOST_PYTHON_HEADERS "${Boost_INCLUDE_DIR}" )
endif()

set (PYTHON_LIBRARIES "/usr/lib/x86_64-linux-gnu/libpython3.8.so")
set (PYTHON_INCLUDE_DIR "/usr/include/python3.8")
set (BOOST_PYTHON_LIBRARY "/usr/local/lib/libboost_python38.so")
set (BOOST_PYTHON_HEADERS "/usr/local/include/boost")

message( STATUS "Python paths ${PYTHON_LIBRARIES} | ${PYTHON_INCLUDE_DIR} | ${PYTHON_SITE_DIR}" )
include_directories( ${PYTHON_INCLUDE_DIR} )
include_directories( "${BOOST_PYTHON_HEADERS}" )

message( STATUS "boost::python library ${BOOST_PYTHON_LIBRARY} | ${BOOST_PYTHON_HEADERS}" )

# set the location to install the python modules
set (SIRE_PYTHON "${PYTHON_SITE_DIR}")

# Use the same Qt library as was used to build corelib
if (DEFINED Qt5Core_DIR)
    message( STATUS "Looking for Qt5Core in directory ${Qt5Core_DIR}" )
    find_package( Qt5Core )
endif()
get_target_property(QT5_LOCATION Qt5::Core LOCATION)
if (QT5_LOCATION MATCHES NOTFOUND)
    message( STATUS "Looking for Qt5Core in directory ${SIRE_QT5CORE_DIR}" )
    set( Qt5Core_DIR "${SIRE_QT5CORE_DIR}/cmake/Qt5Core" )
    find_package( Qt5Core REQUIRED )
endif()
get_target_property(QT5_LOCATION Qt5::Core LOCATION)
message( STATUS "Using Qt5Core in ${QT5_LOCATION}" )

if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
    add_compile_options("/EHsc")
    add_compile_options("/permissive-")
endif()

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
add_subdirectory (Stream)
add_subdirectory (Squire)
add_subdirectory (System)
add_subdirectory (Units)
add_subdirectory (Vol)

#also compile Sire's python executable
add_subdirectory (python)

add_subdirectory (Tools)

install( FILES __init__.py
         DESTINATION ${SIRE_PYTHON}/Sire
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
