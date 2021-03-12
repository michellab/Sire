################################
#
#  CMake Project file for Sire
#
#   (C) Christopher Woods
#
################################

set (S_VERSION_MAJOR "2021")
set (S_VERSION_MINOR "1")
set (S_VERSION_PATCH "0")

set (SIRE_VERSION "${S_VERSION_MAJOR}.${S_VERSION_MINOR}.${S_VERSION_PATCH}")

############
############ SETTING INITIAL SIRE VARIABLES
############

# require cmake >= 3.0.0
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
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

# Option used to control the number of cores for the compile
set ( BUILD_NCORES "0"
      CACHE STRING
      "The number of cores used to compile Sire. Leave this at 0 to automatically work this out." )

string(REGEX MATCH "[0-9]+" NCORES ${BUILD_NCORES})

# get the number of cores on this machine - useful for
# compiling dependencies
if ( NCORES EQUAL 0 )
  include(ProcessorCount)
  ProcessorCount(NCORES)
  if (NCORES EQUAL 0)
    set( NCORES 1 )
  endif()
endif()

if ( CMAKE_GENERATOR MATCHES "Xcode" )
  message( STATUS "Creating Sire Xcode files. Note that this will "
                  "only work if Sire has already been installed..." )

  # Look for the Sire cmake config file that should have been
  # installed with Sire...
  if (NOT DEFINED SIRE_APP)
    set( SIRE_APP "$ENV{HOME}/sire.app" CACHE PATH
            "The full path to the sire.app installed Sire application." )
  endif()

  set( SIRE_INCLUDE_DIR "${SIRE_APP}/include/Sire" )

  set ( SIRE_COMP_CMAKE "${SIRE_INCLUDE_DIR}/cmake/SireCompileVariables.cmake" )

  # Look for the "SireCompileVariables.cmake" file in the install directory
  if (EXISTS "${SIRE_COMP_CMAKE}")
    include( "${SIRE_COMP_CMAKE}" )
  else()
    # try the anaconda path
    set( CONDA_APP "$ENV{HOME}/sire.app/pkgs/sire-${SIRE_VERSION}" CACHE PATH
          "Install directory for Sire within the Anaconda/Miniconda distribution" )

    set( SIRE_INCLUDE_DIR "${CONDA_APP}/include/Sire" )
    set ( SIRE_COMP_CMAKE "${SIRE_INCLUDE_DIR}/cmake/SireCompileVariables.cmake" )

    if (EXISTS "${SIRE_COMP_CMAKE}")
      include( "${SIRE_COMP_CMAKE}" )
    else()
      message( STATUS "You need to provide the full path to sire.app and check that it contains the file ${SIRE_COMP_CMAKE}" )
      message( FATAL_ERROR "Do this by running 'cmake -G Xcode -DSIRE_APP=/path/to/sire.app /path/to/corelib/source'" )
    endif()
  endif()

  # Now set the C++, C and Fortran compilers to be used to compile the wrappers
  # These have to be the same as used to compile corelib
  set(CMAKE_C_COMPILER "${SIRE_C_COMPILER}")
  set(CMAKE_CXX_COMPILER "${SIRE_CXX_COMPILER}")

  message( STATUS "Setting C compiler to ${SIRE_C_COMPILER} and C++ compiler to ${SIRE_CXX_COMPILER}" )
else()
  message(STATUS "Configuring on a machine with number of cores = ${NCORES}")

  # Option to compile and install into an anaconda or miniconda directory
  option( ANACONDA_BUILD "Build into an existing anaconda/miniconda installation" OFF )

  if ( ANACONDA_BUILD )
    message( STATUS "Building into an existing Anaconda/Miniconda installation..." )
    # Where is anaconda installed
    if (NOT DEFINED ANACONDA_BASE)
      set( ANACONDA_BASE "$ENV{HOME}/miniconda3" CACHE PATH
              "Install directory for Anaconda/Miniconda (python3)" )
    endif()

    # Make sure that there is a python3 executable in the bin directory
    if (DEFINED PYTHON_EXECUTABLE)
        if ( NOT EXISTS "${PYTHON_EXECUTABLE}" )
            message( FATAL_ERROR "Cannot find python3. "
                                 "Please set the full path to the Anaconda or Miniconda "
                                 "installation directory." )
        endif()
    else()
        if ( NOT EXISTS "${ANACONDA_BASE}/bin/python3" AND NOT EXISTS "${ANACONDA_BASE}/python.exe" )
            message( FATAL_ERROR "Cannot find python3. "
                                 "Please set the full path to the Anaconda or Miniconda "
                                 "installation directory." )
        endif()
    endif()

    set( SIRE_INSTALL_PREFIX "${ANACONDA_BASE}/pkgs/sire-${SIRE_VERSION}" CACHE PATH
        "Install directory for Sire within the Anaconda/Miniconda distribution" FORCE)

    set(CMAKE_INSTALL_PREFIX "${SIRE_INSTALL_PREFIX}" CACHE INTERNAL
             "Prefix prepended to install directories" FORCE )

#  elseif( WIN32 OR WIN64 )  # need this instead of MSYS as don't know we are MSYS yet...
#
#    message(STATUS "Compiling on windows using MSYS2/MINGW64...")
#    set( MINGW_BASE "C:/msys64/mingw64" CACHE PATH "MSYS MINGW64 root directory" )
#    get_filename_component(MINGW_CMAKE_BIN ${CMAKE_COMMAND} DIRECTORY)
#    get_filename_component(MINGW_CMAKE_ROOT ${MINGW_CMAKE_BIN} DIRECTORY)
#    message(STATUS "Using cmake ${CMAKE_COMMAND} from MING64 base ${MINGW_CMAKE_ROOT}")
#
#    set( MINGW_BASE "${MINGW_CMAKE_ROOT}")
#    set( MINGW_PYTHON "${MINGW_BASE}/bin/python3.exe" )
#
#    if ( NOT EXISTS "${MINGW_PYTHON}" )
#      message(STATUS "Cannot find mingw python3 executable, ${MINGW_PYTHON}")
#      message(STATUS "Can only install Sire into a MINGW installation with a valid python3.exe.")
#      message(STATUS "This is normally in the file C:/msys64/mingw64/bin/python3.exe.")
#      message(STATUS "Please install MSYS2/MINGW64 with the MINGW python3 package.")
#      message(STATUS "Set MINGW_BASE in cmake to point to your MSYS2/MINGW64 installation directory.")
#      message(FATAL_ERROR "Cannot compile until mingw python3 is installed...")
#    endif()
#
#    set( SIRE_INSTALL_PREFIX "${MINGW_BASE}" )
#    set( CMAKE_INSTALL_PREFIX "${SIRE_INSTALL_PREFIX}" CACHE INTERNAL
#           "Prefix prepended to install directories" FORCE )
  else()
    #�Option to set the install directory of Sire
    if (NOT DEFINED SIRE_APP)
      set(SIRE_APP "$ENV{HOME}/sire.app")
    endif()
    set( SIRE_INSTALL_PREFIX "${SIRE_APP}" CACHE PATH
            "Install directory for Sire. This directory can be moved after compilation" )

    set(CMAKE_INSTALL_PREFIX "${SIRE_INSTALL_PREFIX}" CACHE INTERNAL
             "Prefix prepended to install directories" FORCE )
  endif()
endif()

message( STATUS "Sire will be compiled and installed to directory ${CMAKE_INSTALL_PREFIX}" )
set ( SIRE_APP "${ANACONDA_BASE}" )
message( STATUS "sire.app root directory is ${SIRE_APP}" )

# Absolutely can't run cmake in the source directory!
if ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
  message( STATUS "${CMAKE_BINARY_DIR} | ${CMAKE_SOURCE_DIR}" )
  message( FATAL_ERROR "You must run CMake in a different directory to the source!" )
endif()

# name the project
project (SIRE CXX C)

# Whether or not to build the parts of Sire that depend on
# a graphical environment (e.g. inquire)
option ( SIRE_BUILD_GUI "Build graphical applications" OFF )

# Create a file in which we can save the values of all useful variables.
# This will mean that projects using Sire won't have to set these variables
set (SIRE_VARIABLES_FILE "${CMAKE_BINARY_DIR}/SireCompileVariables.cmake")
file ( REMOVE ${SIRE_VARIABLES_FILE} )

# This function is used to save the names of variable to SireCompileVariables.cmake
function( SAVE_SIRE_VARIABLE _var _value )
  file(APPEND ${SIRE_VARIABLES_FILE}
        "set( ${_var} \"${_value}\" )\n")
endfunction (SAVE_SIRE_VARIABLE)

math (EXPR SIRE_VERSION_NUMBER
      "${S_VERSION_MAJOR}*100000 + ${S_VERSION_MINOR}*100 + ${S_VERSION_PATCH}")

set (SIRE_VERSION_STRING "${S_VERSION_MAJOR}_${S_VERSION_MINOR}_${S_VERSION_PATCH}")

save_sire_variable( "SIRE_APP" "${SIRE_APP}" )
save_sire_variable( "S_VERSION_MAJOR" "${S_VERSION_MAJOR}" )
save_sire_variable( "S_VERSION_MINOR" "${S_VERSION_MINOR}" )
save_sire_variable( "S_VERSION_PATCH" "${S_VERSION_PATCH}" )
save_sire_variable( "SIREE_VERSION" "${SIRE_VERSION}" )
save_sire_variable( "SIRE_VERSION_STRING" "${SIRE_VERSION_STRING}" )
save_sire_variable( "SIRE_APP_DIR" "${CMAKE_INSTALL_PREFIX}" )

save_sire_variable( "SIRE_ANACONDA_BUILD" "${ANACONDA_BUILD}" )
save_sire_variable( "BUILD_NCORES" "${NCORES}" )

if (ANACONDA_BUILD)
  save_sire_variable( "SIRE_ANACONDA_BASE" "${ANACONDA_BASE}" )
endif()

# Ensure that the Sire's directories are searched for header files
# before the system directories - this prevents problems in case
# an older version of Sire is already installed on the system
INCLUDE_DIRECTORIES( BEFORE ${SIRE_SOURCE_DIR} )

# Sire must be built as a set of shared libraries
# (otherwise the modules cannot be loaded in python)
set (BUILD_SHARED_LIBS ON)

# This is where you set your own compiler flags
set (SIRE_USER_CFLAGS "" CACHE STRING "User supplied C compiler flags" )
set (SIRE_USER_CPPFLAGS "" CACHE STRING "User supplied C++ compiler flags" )

save_sire_variable( "SIRE_USER_CFLAGS" "${SIRE_USER_CFLAGS}" )
save_sire_variable( "SIRE_USER_CPPFLAGS" "${SIRE_USER_CPPFLAGS}" )

# Name of executable install dir
set (SIRE_BIN "bin")

# Shared library install dir
set (SIRE_LIBS "lib")

# Where to install all bundled components
set (SIRE_BUNDLED_DIR "bundled" )
set (SIRE_BUNDLED_LIBS "${SIRE_BUNDLED_DIR}/lib")

# Static library install dir
set (SIRE_ARCHIVES "lib/static")
# Include (header) files install dir
set (SIRE_INCLUDES "include/Sire")
# CMake install files install dir
set (SIRE_CMAKEFILES "include/Sire/cmake")
# Sire share directory (for parameters etc.)
set (SIRE_SHARE "share/Sire")
# Where will the test files be saved?
set (SIRE_TEST "test")
set( SIRE_SHARE_EXPORT "${SIRE_SHARE}")

# Add these directories to the RPATH variable for the libraries / executables,
# so that they can find each other when loading
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
  set (SIRE_INSTALL_RPATH "@executable_path/.:@executable_path/../${SIRE_LIBS}:@executable_path/../${SIRE_BUNDLED_LIBS}")
else()
  set (SIRE_INSTALL_RPATH "\$ORIGIN/.:\$ORIGIN/../${SIRE_LIBS}:\$ORIGIN/../${SIRE_BUNDLED_LIBS}")
endif()

set (CMAKE_INSTALL_RPATH "${SIRE_INSTALL_RPATH}")
set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

save_sire_variable( "SIRE_INSTALL_PREFIX" "${CMAKE_INSTALL_PREFIX}" )
save_sire_variable( "SIRE_BIN" "${SIRE_BIN}" )
save_sire_variable( "SIRE_LIBS" "${SIRE_LIBS}" )
save_sire_variable( "SIRE_BUNDLED_DIR" "${SIRE_BUNDLED_DIR}" )
save_sire_variable( "SIRE_BUNDLED_LIBS" "${SIRE_BUNDLED_LIBS}" )
save_sire_variable( "SIRE_ARCHIVES" "${SIRE_ARCHIVES}" )
save_sire_variable( "SIRE_INCLUDES" "${SIRE_INCLUDES}" )
save_sire_variable( "SIRE_SHARE" "${SIRE_SHARE}" )
save_sire_variable( "SIRE_CMAKEFILES" "${SIRE_CMAKEFILES}" )
save_sire_variable( "SIRE_INSTALL_RPATH" "${SIRE_INSTALL_RPATH}" )
save_sire_variable( "SIRE_TEST" "${SIRE_TEST}" )

# Uncomment for more verbose compiling/linking
# (or run 'make VERBOSE=1')
# set (CMAKE_VERBOSE_MAKEFILE ON)

# Path to the project's extra cmake files
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/build/cmake ${CMAKE_SOURCE_DIR}/build/cmake/cmake_cxx11)

# Options to pass to 'strip' when stripping the libraries
# (on linux, aix and mac need to use '-x' to ensure we leave
#  the global symbols)
set ( SIRE_STRIP_OPTIONS "-x" CACHE STRING "Options to pass to 'strip'" )
save_sire_variable( "SIRE_STRIP_OPTIONS" "${SIRE_STRIP_OPTIONS}" )

set ( SIRE_STRIP_COMMAND "${CMAKE_STRIP} ${SIRE_STRIP_OPTIONS}" )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_strip.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip"
  IMMEDIATE @ONLY)

set ( CMAKE_STRIP "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip" )

# Add option to turn on or off symbol hiding
option ( SIRE_SYMBOL_HIDING "Turn on library symbol hiding" ON )
save_sire_variable( "SIRE_SYMBOL_HIDING" "${SIRE_SYMBOL_HIDING}" )

# Add option to turn on or off vectorisation of the code
option ( SIRE_VECTORISE "Turn on code vectorisation" ON )
save_sire_variable( "SIRE_VECTORISE" "${SIRE_VECTORISE}" )

# Add option to disable use of SSE
option ( SIRE_DISABLE_SSE "Prevent use of SSE when compiling" OFF )

# Add option to disable use of AVX
option ( SIRE_DISABLE_AVX "Prevent use of AVX when compiling" OFF )

# Add option to disable use of AVX512F)
option ( SIRE_DISABLE_AVX512F "Prevent use of AVX512F when compiling" ON )

# Add an option to force off Fortran compilation and linking
option ( SIRE_DISABLE_FORTRAN "Turn off Fortran compilation and linking" ON )
save_sire_variable( "SIRE_DISABLE_FORTRAN" "${SIRE_DISABLE_FORTRAN}" )

if (SIRE_DISABLE_FORTRAN)
  message( STATUS "Disabling any code that requires compilation or linking to Fortran" )
endif()

# Add an option to say whether or not to use OpenMM
option( SIRE_USE_OPENMM "Set whether or not to compile in support for OpenMM" ON )
save_sire_variable( "SIRE_USE_OPENMM" "${SIRE_USE_OPENMM}" )

#############
############# PLATFORM TESTS
#############

# Now run some platform tests - the results will be placed in config.h
include (CheckIncludeFiles)
include (CheckSymbolExists)
include (CheckFunctionExists)

check_include_files( unistd.h HAVE_UNISTD_H )
check_function_exists( lseek64 HAVE_LSEEK64 )

check_include_files( malloc/malloc.h HAVE_MALLOC_MALLOC_H )
check_include_files( emmintrin.h HAVE_EMMINTRIN_H )
check_include_files( smmintrin.h HAVE_SMMINTRIN_H )
check_include_files( immintrin.h HAVE_IMMINTRIN_H )
check_include_files( execinfo.h HAVE_EXECINFO_H )
check_include_files( stdint.h HAVE_STDINT_H )
check_function_exists( mstats HAVE_MSTATS )
check_function_exists( mallinfo HAVE_MALLINFO )
check_function_exists( sysctl HAVE_SYSCTL )
check_function_exists( cbrt HAVE_CUBEROOT )
check_symbol_exists( backtrace "execinfo.h" HAVE_BACKTRACE )

#create config.h
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/sire_config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/sire_config.h)

#  We need to get the compiler flags for different options
set( SIRE_SMALL_FLAGS "" CACHE INTERNAL "Compiler flags to produce small code" )
set( SIRE_WARNALL_FLAGS "" CACHE INTERNAL "Compiler flags to turn on all warnings" )
set( SIRE_DEBUG_FLAGS "" CACHE INTERNAL "Compiler flags to produce debugging code" )
set( SIRE_RELEASE_FLAGS "" CACHE INTERNAL "Compiler flags to produce optimised release code" )
set( SIRE_VECTOR_FLAGS "" CACHE INTERNAL "Compiler flags to vectorise code (e.g. SSE)" )
set( SIRE_VISIBILITY_FLAGS "" CACHE INTERNAL "Compiler flags to activate symbol visibility" )
set( SIRE_PLATFORM_FLAGS "" CACHE INTERNAL "Additional flags needed by the compiler" )
set( SIRE_SHARE_LINK_FLAGS "" CACHE INTERNAL "Flags needed to link shared libraries" )
set( SIRE_STATIC_LINK_FLAGS "" CACHE INTERNAL "Flags needed to link static libraries" )
set( SIRE_EXE_LINK_FLAGS "" CACHE INTERNAL "Flags needed to link executables" )

#############
############# DETECTING AND CONFIGURING THE COMPILER
#############

if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
  set (SLASH_HYPHEN "/")
  set (COLON_EQUALS ":")
else()
  set (SLASH_HYPHEN "-")
  set (COLON_EQUALS "=")
endif()

include (CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG( "-stdlib=libc++" SIRE_HAS_CPP_LIB )

include(CheckCXXSourceRuns)

if (SIRE_HAS_CPP_LIB)
  message(STATUS "Checking that we can LINK to libc++")
  set(CMAKE_OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
  set(CMAKE_REQUIRED_FLAGS "-stdlib=libc++")
  check_cxx_source_runs("
           #include <cstdlib>
           #include <iostream>
           int main()
           {
              std::cout << std::abs(-5) << std::endl;
              return 0;
           }" CAN_USE_LIBCPP)
  set(CMAKE_REQUIRED_FLAGS ${CMAKE_OLD_REQUIRED_FLAGS})

  if (CAN_USE_LIBCPP)
    message(STATUS "libc++ is available and can be used")
  else()
    message(STATUS "libc++ is not available. Using the normal alternative.")
    set( SIRE_HAS_CPP_LIB 0 )
  endif()
endif()

if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
  CHECK_CXX_COMPILER_FLAG( "/std:c++11" SIRE_HAS_CPP_11 )
  CHECK_CXX_COMPILER_FLAG( "/std:c++1y" SIRE_HAS_CPP_1Y )
  CHECK_CXX_COMPILER_FLAG( "/std:c++14" SIRE_HAS_CPP_14 )
else()
  CHECK_CXX_COMPILER_FLAG( "-std=c++11" SIRE_HAS_CPP_11 )
  CHECK_CXX_COMPILER_FLAG( "-std=c++1y" SIRE_HAS_CPP_1Y )
  CHECK_CXX_COMPILER_FLAG( "-std=c++14" SIRE_HAS_CPP_14 )
endif()

save_sire_variable( "SIRE_HAS_CPP_LIB" "${SIRE_HAS_CPP_LIB}" )

save_sire_variable( "SIRE_HAS_CPP_11" "${SIRE_HAS_CPP_11}" )
save_sire_variable( "SIRE_HAS_CPP_1Y" "${SIRE_HAS_CPP_1Y}" )
save_sire_variable( "SIRE_HAS_CPP_14" "${SIRE_HAS_CPP_14}" )

option( SIRE_SKIP_LIBC++ "Force Sire to ignore libc++ (e.g. in case of broken clang on linux" OFF )

if (SIRE_SKIP_LIBC++)
  if (SIRE_HAS_CPP_LIB)
    message(STATUS "Deliberately ignoring the presence of libstdc++")
  endif()
else()
  if (SIRE_HAS_CPP_LIB)
    message(STATUS "Using libc++ in place of libstdc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    set(SIRE_CXX_NG_FLAGS "-stdlib=libc++")
  endif()
endif()

set(NEED_UNDEF_STRICT_ANSI FALSE)

if (SIRE_HAS_CPP_14)
  message(STATUS "Compiler supports C++ 2014")
  set(SIRE_CXX_NG_FLAGS "${SIRE_CXX_NG_FLAGS} ${SLASH_HYPHEN}std${COLON_EQUALS}c++14 ${SLASH_HYPHEN}DSIRE_HAS_CPP_14 ${SLASH_HYPHEN}DSIRE_HAS_CPP_1Y ${SLASH_HYPHEN}DSIRE_HAS_CPP_11")

  if (NEED_UNDEF_STRICT_ANSI)
    message( STATUS "Undefining __STRICT_ANSI__ because of compile problems")
    set(SIRE_CXX_NG_FLAGS "${SIRE_CXX_NG_FLAGS} -U__STRICT_ANSI__")
  endif()
elseif (SIRE_HAS_CPP_1Y)
  message(STATUS "Compiler supports C++ 201Y (pre-2014)")
  set(SIRE_CXX_NG_FLAGS "${SIRE_CXX_NG_FLAGS} ${SLASH_HYPHEN}std=c++1y ${SLASH_HYPHEN}DSIRE_HAS_CPP_1Y ${SLASH_HYPHEN}DSIRE_HAS_CPP_11")

  if (NEED_UNDEF_STRICT_ANSI)
    message( STATUS "Undefining __STRICT_ANSI__ because of compile problems")
    set(SIRE_CXX_NG_FLAGS "${SIRE_CXX_NG_FLAGS} -U__STRICT_ANSI__")
  endif()
elseif (SIRE_HAS_CPP_11)
  message(STATUS "Compiler supports C++ 2011")
  set(SIRE_CXX_NG_FLAGS "${SIRE_CXX_NG_FLAGS} ${SLASH_HYPHEN}std${COLON_EQUALS}c++11 ${SLASH_HYPHEN}DSIRE_HAS_CPP_11")

  if (NEED_UNDEF_STRICT_ANSI)
    message( STATUS "Undefining __STRICT_ANSI__ because of compile problems")
    set(SIRE_CXX_NG_FLAGS "${SIRE_CXX_NG_FLAGS} -U__STRICT_ANSI__")
  endif()
else()
  message( FATAL_ERROR "Cannot compile Sire as we now require a modern C++ "
                       "compiler that supports at least C++ 2011 (-std=c++11)")
endif()

save_sire_variable( "SIRE_CXX_NG_FLAGS" "${SIRE_CXX_NG_FLAGS}" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SIRE_CXX_NG_FLAGS}")

# Find which C++ 2011/Y/4 features are supported
include (CheckCXX11Features)

message( STATUS "C++11 supported features: ${CXX11_FEATURE_LIST}" )
set( SIRE_COMPILER "UNKNOWN" )

if (NOT ${CMAKE_C_COMPILER_ID} STREQUAL ${CMAKE_CXX_COMPILER_ID})
  message( WARNING "It is not sensible to compile Sire using different C and C++ compilers. "
                   "You are using ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID}) and "
                   "${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
endif()

string( TOUPPER ${CMAKE_CXX_COMPILER_ID} SIRE_CXX_COMPILER_ID )

# TEMP - TESTING VECTOR DETECTOR
include (FindVector)

# Add compiler specific flags
if ( ${SIRE_CXX_COMPILER_ID} MATCHES "GNU" )

    message( STATUS "Detected and using GCC C and C++ compilers" )
    set( SIRE_COMPILER "GCC" )

    # get the compiler version (via gcc -dumpversion)
    # (I'll assume that we are using the same version
    # of compiler for C and C++...)
    EXEC_PROGRAM( ${CMAKE_CXX_COMPILER}
                  ARGS -dumpversion
                  OUTPUT_VARIABLE exe_output )

    # These two regexps are copied from the Chicken CMake
    # files (which provide a very useful cmake tutorial - thanks :-)
    SET( VERSION_BUILD_REGEX "^([0-9]+)[.]([0-9]+).*$" )
    STRING( REGEX REPLACE ${VERSION_BUILD_REGEX} "\\1"
            GCC_MAJOR_VERSION ${exe_output} )
    STRING( REGEX REPLACE ${VERSION_BUILD_REGEX} "\\2"
            GCC_MINOR_VERSION ${exe_output} )

    message( STATUS "Using GCC version "
                    "${GCC_MAJOR_VERSION}.${GCC_MINOR_VERSION}" )

    # First, we need at least GCC 4.8, as we need C++ 11 support
    if ( GCC_MAJOR_VERSION LESS 5 )
       if ( GCC_MINOR_VERSION LESS 8 )
         message( FATAL_ERROR "Sire requires GCC >= 4.8. Please upgrade "
                        "your version of GCC." )
       endif()
    endif()

    set ( SIRE_SMALL_FLAGS  "-Os" )
    set ( SIRE_WARNALL_FLAGS  "-Wall" )
    set ( SIRE_DEBUG_FLAGS   "-g" )

    if ( GCC_MAJOR_VERSION GREATER 4 )
        # GCC >= 5 still show weird crashes when using O3 which don't
        # show up with other compilers. It is probably my fault, but
        # until I find out why, we have to limit compilation to -O2
        set ( SIRE_RELEASE_FLAGS "-O2 -ffast-math" )
    else()
        set ( SIRE_RELEASE_FLAGS "-O3 -ffast-math" )
    endif()

    # -fomit-frame-pointer breaks backtrace on an apple
    set ( SIRE_RELEASE_FLAGS "${SIRE_RELEASE_FLAGS} -fomit-frame-pointer" )

    set ( SIRE_PLATFORM_FLAGS "-pipe -DSIRE_ALWAYS_INLINE=inline" )
    set ( SIRE_STATIC_LINK_FLAGS "-static" )

    GET_SIRE_VECTOR_FLAGS( "-fopenmp-simd" "-msse2" "-mavx" "-mavx512f" )

    # Now gcc 5 specific options
    if ( GCC_MAJOR_VERSION GREATER 4 )
      if (MSYS)
        message(STATUS "MSYS2 will use builtin OpenMM if available...")
      else()
        # OpenMM with conda uses the old C++ binary API!
        # Tell GCC 5 to respect the old API
        set( SIRE_PLATFORM_FLAGS "${SIRE_PLATFORM_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0" )
      endif()
    endif()

    # support for ppc64le architecture
    if ( CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "ppc64le" )
      set( SIRE_PLATFORM_FLAGS "${SIRE_PLATFORM_FLAGS} -mcpu=power9 -mtune=power9 -mpowerpc64 -mvsx" )
    endif()

    # remove warnings about strict aliasing
    set( SIRE_PLATFORM_FLAGS "${SIRE_PLATFORM_FLAGS} -Wno-strict-aliasing" )
    set ( SIRE_VISIBILITY_FLAGS "-DSIRE_VISIBILITY_AVAILABLE -fvisibility=hidden -fvisibility-inlines-hidden" )

    # Ensure that the libraries don't contain any undefined symbols
    if(APPLE)
        set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} -Wl" )
    else(UNIX)
        CHECK_CXX_COMPILER_FLAG( "-rdynamic" HAVE_RDYNAMIC_FLAG )

        if (HAVE_RDYNAMIC_FLAG)
          set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} -rdynamic -Wl,--no-undefined" )
        else()
          set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} -Wl,--no-undefined" )
        endif()
    endif()

elseif ( ${SIRE_CXX_COMPILER_ID} MATCHES "CLANG" )

    message( STATUS "Detected and using clang C and C++ compilers" )
    set( SIRE_COMPILER "CLANG" )

    set ( SIRE_SMALL_FLAGS  "-Os" )
    set ( SIRE_WARNALL_FLAGS  "-Wall" )
    set ( SIRE_DEBUG_FLAGS   "-g" )

    # an old incompatibility between clang and gcc means that some
    # distributions will fail if use -ffast-math
    set(CMAKE_OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_FLAGS "-O3 -ffast-math")
    CHECK_CXX_SOURCE_RUNS("
           #include <cmath>
           #include <iostream>
           int main()
           {
              std::cout << std::sin(0.0) << std::endl;
              return 0;
           }"
           CAN_USE_FFAST_MATH)
    set(CMAKE_REQUIRED_FLAGS ${CMAKE_OLD_REQUIRED_FLAGS})

    if (CAN_USE_FFAST_MATH)
      set ( SIRE_RELEASE_FLAGS "-O3 -ffast-math" )
    else()
      message("Disabling -ffast-math as incompatible clang and gcc")
      set ( SIRE_RELEASE_FLAGS "-O3" )
    endif()

    # -fomit-frame-pointer breaks backtrace on an apple
    set ( SIRE_RELEASE_FLAGS "${SIRE_RELEASE_FLAGS} -fomit-frame-pointer" )

    set ( SIRE_PLATFORM_FLAGS "-pipe -DSIRE_ALWAYS_INLINE=inline" )
    set ( SIRE_STATIC_LINK_FLAGS "-static" )

    GET_SIRE_VECTOR_FLAGS( "-openmp-simd" "-msse2" "-mavx" "-mavx512f" )

    set ( SIRE_VISIBILITY_FLAGS "-DSIRE_VISIBILITY_AVAILABLE -fvisibility=hidden -fvisibility-inlines-hidden" )

    if (APPLE)
    elseif(UNIX)
        set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} -rdynamic -Wl,--no-undefined" )
    endif()

elseif( ${CMAKE_CXX_COMPILER} MATCHES "icpc" )

    message( STATUS "Compiling with the Intel compiler" )
    set( SIRE_COMPILER "INTEL" )

    set ( SIRE_SMALL_FLAGS  "-Os" )
    set ( SIRE_WARNALL_FLAGS  "-w1" )
    set ( SIRE_RELEASE_FLAGS  "-O3 -finline -finline-functions -ansi-alias -fargument-noalias-global -ip -g -simd" )
    set ( SIRE_DEBUG_FLAGS   "-g" )
    set ( SIRE_VISIBILITY_FLAGS "-DSIRE_NO_VISIBILITY_AVAILABLE" )
    set ( SIRE_SHARE_LINK_FLAGS "-shared" )
    set ( SIRE_STATIC_LINK_FLAGS "-static" )

    set ( SIRE_SMALL_FLAGS  "-Os" )
    set ( SIRE_WARNALL_FLAGS  "-Wall" )
    set ( SIRE_DEBUG_FLAGS   "-g" )

    # -fomit-frame-pointer breaks backtrace on an apple
    set ( SIRE_RELEASE_FLAGS "${SIRE_RELEASE_FLAGS} -fomit-frame-pointer" )

    set ( SIRE_PLATFORM_FLAGS "-pipe -DSIRE_ALWAYS_INLINE=inline" )
    set ( SIRE_STATIC_LINK_FLAGS "-static" )

    CHECK_CXX_COMPILER_FLAG( "-msse2" HAVE_SSE_FLAG )
    CHECK_CXX_COMPILER_FLAG( "-mavx" HAVE_AVX_FLAG )

    GET_SIRE_VECTOR_FLAGS( "-qopenmp-simd" "-msse2" "-mavx" "-xCOMMON-AVX512" )

    set ( SIRE_VISIBILITY_FLAGS "-DSIRE_VISIBILITY_AVAILABLE -fvisibility=hidden -fvisibility-inlines-hidden" )

    if (APPLE)
    elseif(UNIX)
        set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} -rdynamic -Wl,--no-undefined" )
    endif()

elseif (MSVC)

    add_compile_options("/EHsc")
    add_compile_options("/permissive-")
    add_compile_options("/Zc:twoPhase-")
    message( STATUS "Compiling with MSVC" )
    set( SIRE_COMPILER "MSVC" )

    set ( SIRE_SMALL_FLAGS  "/O1" )
    set ( SIRE_WARNALL_FLAGS  "" )
    set ( SIRE_RELEASE_FLAGS  "/O2 /GL /Gw" )
    set ( SIRE_DEBUG_FLAGS   "/Zi" )
    set ( SIRE_VISIBILITY_FLAGS "/DSIRE_NO_VISIBILITY_AVAILABLE" )
    GET_SIRE_VECTOR_FLAGS( "/openmp:experimental" "/arch:SSE2" "/arch:AVX" "/arch:AVX512" )
    set ( SIRE_SHARE_LINK_FLAGS "/LTCG /OPT:REF /OPT:ICF" )
    set ( SIRE_STATIC_LINK_FLAGS "/LTCG /OPT:REF /OPT:ICF" )
    set ( SIRE_EXE_LINK_FLAGS "/LTCG /OPT:REF /OPT:ICF" )

    set ( SIRE_PLATFORM_FLAGS "/DSIRE_ALWAYS_INLINE=__forceinline" )

else()
    message( STATUS "CMAKE_SYSTEM_NAME == ${CMAKE_SYSTEM_NAME}" )
    message( STATUS "CMAKE_C_COMPILER == ${CMAKE_C_COMPILER}" )
    message( STATUS "CMAKE_CXX_COMPILER == ${CMAKE_CXX_COMPILER}" )
    message( STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}" )
    message( STATUS "SIRE_CXX_COMPILER_ID = ${SIRE_CXX_COMPILER_ID}" )
    message( FATAL_ERROR "Could not detect compiler type - cannot set compiler options." )
endif()

save_sire_variable( "SIRE_COMPILER" "${SIRE_COMPILER}" )

################
################ SORTING OUT ALL OF SIRE'S DEPENDENCIES
################

# Tell CMake to run moc when necessary
set(CMAKE_AUTOMOC OFF)

# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if ( CMAKE_GENERATOR MATCHES "Xcode" )
  # Xcode build - use the libraries that are already installed...

else()
  # Option to disable all bundles
  option( SIRE_ENABLE_BUNDLING "Allow Sire to use bundled libraries" ON )

  # Fine-grain control of bundling each dependency
  option( SIRE_BUNDLE_TBB "Whether Sire should use its bundled copy of TBB" ON )
  option( SIRE_BUNDLE_GSL "Whether Sire should use its bundled copy of GSL" ON )
  option( SIRE_BUNDLE_CPUID "Whether Sire should use its bundled copy of libcpuid" ON )
  option( SIRE_BUNDLE_BOOST "Whether Sire should use its bundled copy of Boost" ON )
  option( SIRE_BUNDLE_BLASLAPACK "Whether Sire should use its bundled copy of BLAS and LAPACK" ON )
  option( SIRE_BUNDLE_QT "Whether Sire should use its bundled copy of Qt" ON )
  option( SIRE_BUNDLE_OPENMM "Whether Sire should use its bundled copy of OpenMM" OFF )

  save_sire_variable( "SIRE_ENABLE_BUNDLING" "${SIRE_ENABLE_BUNDLING}" )
  save_sire_variable( "SIRE_BUNDLE_TBB" "${SIRE_BUNDLE_TBB}" )
  save_sire_variable( "SIRE_BUNDLE_CPUID" "${SIRE_BUNDLE_CPUID}" )
  save_sire_variable( "SIRE_BUNDLE_GSL" "${SIRE_BUNDLE_GSL}" )
  save_sire_variable( "SIRE_BUNDLE_BOOST" "${SIRE_BUNDLE_BOOST}" )
  save_sire_variable( "SIRE_BUNDLE_BLASLAPACK" "${SIRE_BUNDLE_BLASLAPACK}" )
  save_sire_variable( "SIRE_BUNDLE_QT" "${SIRE_BUNDLE_QT}" )
  save_sire_variable( "SIRE_BUNDLE_OPENMM" "${SIRE_BUNDLE_OPENMM}" )

  set (SIRE_FOUND_TBB FALSE)
  set (SIRE_FOUND_GSL FALSE)
  set (SIRE_FOUND_BOOST FALSE)
  set (SIRE_FOUND_BLASLAPACK FALSE)
  set (SIRE_FOUND_QT FALSE)
  set (SIRE_FOUND_OPENMM FALSE)
  set (SIRE_FOUND_NETCDF FALSE)
  set (SIRE_FOUND_CPUID FALSE)

  if ( ${SIRE_ENABLE_BUNDLING} )
    message( STATUS "Sire will be compiled to use its own copy of bundled libraries." )

    set( BUNDLE_BUILDDIR "${CMAKE_BINARY_DIR}/build_bundled" )

    if ( MSYS )
      message( STATUS "No bundling - using libraries provided by MSYS2/MINGW64...")
      set( BUNDLE_STAGEDIR "${CMAKE_INSTALL_PREFIX}" )
    else()
      set( BUNDLE_STAGEDIR "${CMAKE_INSTALL_PREFIX}/${SIRE_BUNDLED_DIR}" )

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
    endif()

    if ( ${SIRE_BUNDLE_TBB} )
      include( "src/bundled/install_tbb.cmake" )
    endif()

    if ( ${SIRE_BUNDLE_CPUID} )
      include( "src/bundled/install_cpuid.cmake" )
    endif()

    if ( ${SIRE_BUNDLE_GSL} )
      include( "src/bundled/install_gsl.cmake" )

      if (${SIRE_FOUND_GSL})
        list(APPEND SIREMATHS_EXTRA_LIBRARIES ${GSL_LIBRARY})
        list(APPEND SIREMATHS_EXTRA_LIBRARIES ${GSL_CBLAS_LIBRARY})
      endif()

    endif()

    if ( ${SIRE_BUNDLE_BOOST} )
      include( "src/bundled/install_boost.cmake" )
    endif()

    if ( ${SIRE_BUNDLE_QT} )
      include( "src/bundled/install_qt5.cmake" )
    endif()

    if ( NOT MSYS )
      if ( NOT ${BUNDLE_STAGEDIR} STREQUAL "${CMAKE_INSTALL_PREFIX}/${SIRE_BUNDLED_DIR}" )
        # All of the bundled objects are built into the directory ${CMAKE_BINARY_DIR}/bundled
        # This directory must be installed into the sire.app/ directory
        install( DIRECTORY ${BUNDLE_STAGEDIR}
                 DESTINATION ${CMAKE_INSTALL_PREFIX}
                 USE_SOURCE_PERMISSIONS )
      endif()
    endif()
  endif()
endif()

if ( NOT ${SIRE_FOUND_QT} )
    # look for Qt5 in the system
    find_package( Qt5Core REQUIRED )
endif()

if ( NOT ${SIRE_FOUND_NETCDF} )
    # look for netcdf in the system (ideally from conda)
    if (ANACONDA_BUILD)
        set(NetCDF_ROOT_DIR "${ANACONDA_BASE}")
    endif()

    FIND_PACKAGE( NetCDF )

    if ( ${NetCDF_FOUND} )
        message(STATUS "Found NetCDF in ${NetCDF_INCLUDE_DIR} | ${NetCDF_LIBRARIES}")
        set( SIRE_FOUND_NETCDF 1 )
        save_sire_variable( "SIRE_NetCDF_INCLUDE_DIR" "${NetCDF_INCLUDE_DIR}" )
        save_sire_variable( "SIRE_NetCDF_LIBRARIES" "${NetCDF_LIBRARIES}" )
        save_sire_variable( "SIRE_NetCDF_FOUND" "${NetCDF_FOUND}" )
    else()
        message(STATUS "Cannot find NetCDF - disabling its use!")
        save_sire_variable( "SIRE_NetCDF_FOUND" "${NetCDF_FOUND}" )
        set( SIRE_FOUND_NETCDF 0 )
    endif()
else()
    message(STATUS "Using NetCDF from ${NetCDF_INCLUDE_DIR} | ${NetCDF_LIBRARIES}")
endif()

if ( SIRE_USE_OPENMM )
  if ( NOT ${SIRE_FOUND_OPENMM} )
      # look for OpenMM in the system
      # Now look for OpenMM
      if (ANACONDA_BUILD)
        set(OpenMM_ROOT_DIR "${ANACONDA_BASE}")
        set(OPENMM_ROOT_DIR "${ANACONDA_BASE}")
        FIND_PACKAGE( OpenMM )
      else()
        option (OPENMM_HOME "The location of the OpenMM installation" "/usr/local/openmm")
        set(OpenMM_ROOT_DIR "${OPENMM_HOME}")
        FIND_PACKAGE( OpenMM )
      endif()

      save_sire_variable( "SIRE_OpenMM_FOUND" "${OpenMM_FOUND}" )

      if ( ${OpenMM_FOUND} )
         message( STATUS "Found OpenMM at ${OpenMM_ROOT_DIR}" )
         save_sire_variable( "SIRE_OpenMM_INCLUDE_DIR" "${OpenMM_INCLUDE_DIR}" )
         save_sire_variable( "SIRE_OpenMM_LIBRARIES" "${OpenMM_LIBRARIES}" )
         save_sire_variable( "SIRE_OpenMM_LIBRARY_DIR" "${OpenMM_LIBRARY_DIR}" )
      else()
         message( STATUS "Could not find OpenMM. OpenMM support is disabled." )
      endif()
  else()
    message( STATUS "Disabling use of OpenMM at user request.")
    save_sire_variable( "SIRE_OpenMM_FOUND" "${SIRE_USE_OPENMM}" )
  endif()
endif()

get_target_property(QT5_LOCATION Qt5::Core LOCATION)
get_filename_component(SIRE_QT5CORE_DIR ${QT5_LOCATION} PATH)
message( STATUS "Qt5Core location ${SIRE_QT5CORE_DIR}" )
save_sire_variable( "SIRE_QT5CORE_DIR" "${SIRE_QT5CORE_DIR}" )

if (DEFINED SIRE_APP)
  include_directories( "${SIRE_APP}/include" "${SIRE_APP}/Library/include" )
endif()

if (NOT ${SIRE_FOUND_BOOST})
  # Sire depends on >= boost 1.31 - all libraries must be dynamically linked
  set (BOOST_ALL_DYN_LINK "YES")
  set (Boost_ADDITIONAL_VERSIONS "1.39" "1.39.0")
  FIND_PACKAGE ( Boost 1.31 REQUIRED )

  if (Boost_FOUND)
    message(STATUS "Boost paths ${Boost_LIBRARY_DIRS} | ${Boost_INCLUDE_DIR}" )
    set ( BOOST_INCLUDE_DIRS "${Boost_INCLUDE_DIR}" )
    include_directories( ${Boost_INCLUDE_DIR} )

    #save the path to this include directory so that it can be
    #used by anything compiling against Sire
    save_sire_variable( "SIRE_BOOST_INCLUDE_DIR" "${Boost_INCLUDE_DIR}" )
  else()
    message(FATAL_ERROR "Cannot find the boost libraries.")
  endif()
endif()

if (NOT ${SIRE_FOUND_GSL})
  # Sire depends on >= GSL 1.6
  find_package(GSL 1.6 QUIET)

  message( STATUS "GSL paths ${GSL_LINK_DIRECTORIES} "
                  "${GSL_LIBRARIES} | ${GSL_INCLUDE_DIR}" )

  list(APPEND SIREMATHS_EXTRA_LIBRARIES ${GSL_LIBRARIES})

  include_directories( BEFORE ${GSL_INCLUDE_DIR})

  save_sire_variable( "SIRE_GSL_INCLUDE_PATH" "${GSL_INCLUDE_DIR}" )
endif()

if (NOT ${SIRE_FOUND_BLASLAPACK})
  # Now look for BLAS and LAPACK - if these are not found then
  # a local version will have to be used instead
  set( BLA_STATIC OFF )
  find_package( BLAS )
  find_package( LAPACK )

  if (BLAS_FOUND)
    message( STATUS "Found BLAS library ${BLAS_LIBRARIES}" )
    set( SIRE_USE_BLAS ON )
    save_sire_variable( "SIRE_BLAS_LIBRARIES" "${BLAS_LIBRARIES}" )
  else()
    message( STATUS "Using internal (non-optimal) BLAS functions." )
    save_sire_variable( SIRE_USE_BLAS OFF )
  endif()

  save_sire_variable( "SIRE_USE_BLAS" "${SIRE_USE_BLAS}" )

  if (LAPACK_FOUND)
    message( STATUS "Found LAPACK library ${LAPACK_LIBRARIES}" )
    set( SIRE_USE_LAPACK ON )
    save_sire_variable( "SIRE_LAPACK_LIBRARIES" "${LAPACK_LIBRARIES}" )
  else()
    message( STATUS "Using internal (non-optimal) LAPACK functions." )
    set( SIRE_USE_LAPACK OFF )
  endif()

  save_sire_variable( "SIRE_USE_LAPACK" "${SIRE_USE_LAPACK}" )
endif()

if (NOT ${SIRE_FOUND_CPUID})
  if (NOT ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
    find_library( CPUID_LIBRARY "cpuid" PATHS ${CPUID_LIBRARY_DIR} )
    include_directories (${CPUID_INCLUDE_DIR} )
    set(SIRE_FOUND_CPUID TRUE)
    if (HAVE_STDINT_H)
      set( CPUID_DEFINITIONS "-DHAVE_STDINT_H" )
    endif()
  else()
    message( STATUS "Cannot find libcpuid. Will disable CPU detection code." )
  endif()
endif()

if (NOT ${SIRE_FOUND_TBB})
  find_library( TBB_LIBRARY "tbb" PATHS ${TBB_LIBRARY_DIR} )
  find_library( TBB_MALLOC_LIBRARY "tbbmalloc" PATHS ${TBB_LIBRARY_DIR} )
  include_directories(${TBB_INCLUDE_DIR})
  save_sire_variable( "TBB_LIBRARY" "${TBB_LIBRARY}" )
  save_sire_variable( "TBB_MALLOC_LIBRARY" "${TBB_MALLOC_LIBRARY}" )
  save_sire_variable( "TBB_INCLUDE_DIR" "${TBB_INCLUDE_DIR}" )
endif()

###########
########### PRINTING OUT SETTINGS AND FINALISING COMPILE FLAGS
###########

set( SIRE_DISABLE_WARNINGS ON )

if (SIRE_DISABLE_WARNINGS)
  # Disable some of the unnecessary compiler warnings...
  CHECK_CXX_COMPILER_FLAG( "-Wno-attributes" SIRE_HAS_WNO_ATTRIBUTES )

  if (SIRE_HAS_WNO_ATTRIBUTES)
    list( APPEND SIRE_DISABLE_WARNING_FLAGS "-Wno-attributes" )
  endif()

  set( SIRE_WARNALL_FLAGS "${SIRE_WARNALL_FLAGS} ${SIRE_DISABLE_WARNING_FLAGS}" )

endif()

message( STATUS  "CMAKE_SYSTEM_NAME      == ${CMAKE_SYSTEM_NAME}" )
message( STATUS  "CMAKE_C_COMPILER       == ${CMAKE_C_COMPILER}" )
message( STATUS  "CMAKE_CXX_COMPILER     == ${CMAKE_CXX_COMPILER}" )
message( STATUS  "CMAKE_Fortran_COMPILER == ${CMAKE_Fortran_COMPILER}" )

message ( STATUS "SIRE_CXX_NG_FLAGS       == ${SIRE_CXX_NG_FLAGS}" )
message ( STATUS "SIRE_SMALL_FLAGS       == ${SIRE_SMALL_FLAGS}" )
message ( STATUS "SIRE_WARNALL_FLAGS     == ${SIRE_WARNALL_FLAGS}" )
message ( STATUS "SIRE_DEBUG_FLAGS       == ${SIRE_DEBUG_FLAGS}" )
message ( STATUS "SIRE_RELEASE_FLAGS     == ${SIRE_RELEASE_FLAGS}" )
message ( STATUS "SIRE_VECTOR_FLAGS      == ${SIRE_VECTOR_FLAGS}" )
message ( STATUS "SIRE_VISIBILITY_FLAGS  == ${SIRE_VISIBILITY_FLAGS}" )
message ( STATUS "SIRE_PLATFORM_FLAGS    == ${SIRE_PLATFORM_FLAGS}" )
message ( STATUS "SIRE_SHARE_LINK_FLAGS  == ${SIRE_SHARE_LINK_FLAGS}" )
message ( STATUS "SIRE_STATIC_LINK_FLAGS == ${SIRE_STATIC_LINK_FLAGS}" )
message ( STATUS "SIRE_EXE_LINK_FLAGS    == ${SIRE_EXE_LINK_FLAGS}" )

set( SIRE_COMPILE_FLAGS "${SIRE_WARNALL_FLAGS} ${SIRE_PLATFORM_FLAGS}" )

if ( SIRE_SYMBOL_HIDING )
  set ( SIRE_COMPILE_FLAGS "${SIRE_COMPILE_FLAGS} ${SIRE_VISIBILITY_FLAGS}" )
  #set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} ${SIRE_VISIBILITY_FLAGS}" )
endif()

# Add the small flags here so that they don't accidentally pull in the
# vectorisation flags (this breaks the MultiFloat/etc classes for AVX-512F)
set( SIRE_C_FLAGS_SMALL "${SIRE_COMPILE_FLAGS} ${SIRE_SMALL_FLAGS} ${SIRE_USER_CFLAGS}" )
set( SIRE_CXX_FLAGS_SMALL "${SIRE_COMPILE_FLAGS} ${SIRE_CXX_NG_FLAGS} ${SIRE_SMALL_FLAGS} ${SIRE_USER_CPPFLAGS}" )

if ( SIRE_VECTORISE )
  set ( SIRE_COMPILE_FLAGS "${SIRE_COMPILE_FLAGS} ${SIRE_VECTOR_FLAGS}" )
endif()

set( SIRE_C_FLAGS_RELEASE "${SIRE_COMPILE_FLAGS} ${SIRE_RELEASE_FLAGS} ${SIRE_USER_CFLAGS}" )
set( SIRE_CXX_FLAGS_RELEASE "${SIRE_COMPILE_FLAGS} ${SIRE_CXX_NG_FLAGS} ${SIRE_RELEASE_FLAGS}  ${SIRE_USER_CPPFLAGS}" )

set( SIRE_C_FLAGS_DEBUG "${SIRE_COMPILE_FLAGS} ${SIRE_DEBUG_FLAGS} ${SIRE_USER_CFLAGS}" )
set( SIRE_CXX_FLAGS_DEBUG "${SIRE_COMPILE_FLAGS} ${SIRE_CXX_NG_FLAGS} ${SIRE_DEBUG_FLAGS} ${SIRE_USER_CPPFLAGS}" )

# Copy the Sire compile flags to CMake
set( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${SIRE_C_FLAGS_RELEASE}" )
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${SIRE_CXX_FLAGS_RELEASE}" )
set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${SIRE_C_FLAGS_DEBUG}" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${SIRE_CXX_FLAGS_DEBUG}" )
set( CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} ${SIRE_C_FLAGS_SMALL}" )
set( CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${SIRE_CXX_FLAGS_SMALL}" )

if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    message( STATUS "Compiling Sire using a DEBUG build" )
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEBUG}" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}" )
else()
    message( STATUS "Compiling Sire using a default (RELEASE) build" )
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_RELEASE}" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}" )
endif()

message ( STATUS "C compiler flags       == ${CMAKE_C_FLAGS}" )
message ( STATUS "C++ compiler flags     == ${CMAKE_CXX_FLAGS}" )

# Add this to all shared libraries
set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${SIRE_SHARE_LINK_FLAGS}" )
set( CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${SIRE_STATIC_LINK_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${SIRE_EXE_LINK_FLAGS}" )

message( STATUS "Shared library flags    == ${CMAKE_SHARED_LINKER_FLAGS}" )
message( STATUS "Static library flags    == ${CMAKE_STATIC_LINKER_FLAGS}" )
message( STATUS "Executable link flags   == ${CMAKE_EXE_LINKER_FLAGS}" )

# Now save the compile variables so that they can be re-used by libraries that use Sire
save_sire_variable( "SIRE_C_COMPILER" "${CMAKE_C_COMPILER}" )
save_sire_variable( "SIRE_CXX_COMPILER" "${CMAKE_CXX_COMPILER}" )
save_sire_variable( "SIRE_C_FLAGS_RELEASE" "${CMAKE_C_FLAGS_RELEASE}" )
save_sire_variable( "SIRE_CXX_FLAGS_RELEASE" "${CMAKE_CXX_FLAGS_RELEASE}" )
save_sire_variable( "SIRE_C_FLAGS_SMALL" "${CMAKE_C_FLAGS_MINSIZEREL}" )
save_sire_variable( "SIRE_CXX_FLAGS_SMALL" "${CMAKE_CXX_FLAGS_MINSIZEREL}" )
save_sire_variable( "SIRE_C_FLAGS_DEBUG" "${CMAKE_C_FLAGS_DEBUG}" )
save_sire_variable( "SIRE_CXX_FLAGS_DEBUG" "${CMAKE_CXX_FLAGS_DEBUG}" )
save_sire_variable( "SIRE_SHARE_LINK_FLAGS" "${SIRE_SHARE_LINK_FLAGS}" )
save_sire_variable( "SIRE_STATIC_LINK_FLAGS" "${SIRE_STATIC_LINK_FLAGS}" )
save_sire_variable( "SIRE_EXE_LINK_FLAGS" "${SIRE_EXE_LINK_FLAGS}" )
save_sire_variable( "SIRE_SMALL_FLAGS" "${SIRE_SMALL_FLAGS}" )
save_sire_variable( "SIRE_WARNALL_FLAGS" "${SIRE_WARNALL_FLAGS}" )
save_sire_variable( "SIRE_DEBUG_FLAGS" "${SIRE_DEBUG_FLAGS}" )
save_sire_variable( "SIRE_RELEASE_FLAGS" "${SIRE_RELEASE_FLAGS}" )
save_sire_variable( "SIRE_VECTOR_FLAGS" "${SIRE_VECTOR_FLAGS}" )
save_sire_variable( "SIRE_VISIBILITY_FLAGS" "${SIRE_VISIBILITY_FLAGS}" )
save_sire_variable( "SIRE_PLATFORM_FLAGS" "${SIRE_PLATFORM_FLAGS}" )

#include the build directory so that we can get sire_config.h and sire_version.h
include_directories(${CMAKE_BINARY_DIR})

################
################ TELLING CMAKE WHERE THE SOURCE IS LOCATED
################

# Add the compiler test directory - this allows us to run some compiler
# tests to make sure that Sire will compile properly
add_subdirectory (build/test_compiler)

# Recurse into the src/libs and src/apps subdirectories
add_subdirectory (src/libs)
add_subdirectory (src/apps)

################
################ CONFIG FILES AND INSTALLATION INSTRUCTIONS
################

# Finally, write the version header file and finish off the component dependencies file
include( GetGitVersion )

# install the extra includes file in the Sire/cmake directory
install (FILES ${SIRE_VARIABLES_FILE}
         DESTINATION ${SIRE_CMAKEFILES}
        )

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/sire_config.h
         DESTINATION ${SIRE_INCLUDES})

# install the parameters file (containing forcefield parameters) into
# share/Sire
install (DIRECTORY ${CMAKE_SOURCE_DIR}/parameters
         DESTINATION ${SIRE_SHARE}
         USE_SOURCE_PERMISSIONS
         PATTERN ".svn" EXCLUDE
        )

# install the molecular templates into share/Sire
install (DIRECTORY ${CMAKE_SOURCE_DIR}/templates
         DESTINATION ${SIRE_SHARE}
         USE_SOURCE_PERMISSIONS
         PATTERN ".svn" EXCLUDE
        )

install (DIRECTORY ${CMAKE_SOURCE_DIR}/share/build
         DESTINATION ${SIRE_SHARE}
         USE_SOURCE_PERMISSIONS
         PATTERN ".svn" EXCLUDE
        )

# These commands are used to create an 'uninstall' target
# (this is copied from the vtk.org wiki)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
