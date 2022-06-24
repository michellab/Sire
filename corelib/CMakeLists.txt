################################
#
#  CMake Project file for Sire
#
################################

# require cmake >= 3.3.0 (released in 2018)
cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)
cmake_policy(VERSION 3.3)

# Read the version number (so we can set it once in one place)
file (STRINGS "../version.txt" VERSION_TEXT)

string(REGEX MATCH "([0-9]*).([0-9]*).([0-9]*)" _ ${VERSION_TEXT})

set (S_VERSION_MAJOR ${CMAKE_MATCH_1})
set (S_VERSION_MINOR ${CMAKE_MATCH_2})
set (S_VERSION_PATCH ${CMAKE_MATCH_3})
set (SIRE_VERSION "${S_VERSION_MAJOR}.${S_VERSION_MINOR}.${S_VERSION_PATCH}")

message( STATUS "Building Sire ${SIRE_VERSION}")

############
############ SETTING INITIAL SIRE VARIABLES
############

# use loose loop syntax in CMakeLists files
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

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

message(STATUS "Configuring on a machine with number of cores = ${NCORES}")

message( STATUS "Building into an existing conda installation..." )

# Where is anaconda installed
if (NOT DEFINED ANACONDA_BASE)
  set( ANACONDA_BASE "$ENV{HOME}/miniconda3" CACHE PATH
          "Install directory for conda (python3)" )
endif()

# Make sure that there is a python3 executable in the bin directory
if (DEFINED PYTHON_EXECUTABLE)
    if ( NOT EXISTS "${PYTHON_EXECUTABLE}" )
        message( FATAL_ERROR "Cannot find python3. "
                              "Please set the full path to the conda "
                              "installation directory." )
    endif()
else()
    if ( NOT EXISTS "${ANACONDA_BASE}/bin/python3" AND NOT EXISTS "${ANACONDA_BASE}/python.exe" )
        message( FATAL_ERROR "Cannot find python3. "
                              "Please set the full path to the conda "
                              "installation directory." )
    endif()
endif()

set( SIRE_INSTALL_PREFIX "${ANACONDA_BASE}/pkgs/sire-${SIRE_VERSION}" CACHE PATH
        "Install directory for Sire within the conda" FORCE)

set(CMAKE_INSTALL_PREFIX "${SIRE_INSTALL_PREFIX}" CACHE INTERNAL
          "Prefix prepended to install directories" FORCE )

message( STATUS "Sire will be compiled and installed to directory ${CMAKE_INSTALL_PREFIX}" )
set ( SIRE_APP "${ANACONDA_BASE}" )
message( STATUS "sire.app root directory is ${SIRE_APP}" )

# Absolutely can't run cmake in the source directory!
if ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
  message( STATUS "${CMAKE_BINARY_DIR} | ${CMAKE_SOURCE_DIR}" )
  message( FATAL_ERROR "You must run CMake in a different directory to the source!" )
endif()

# conda can set incorrect CXX flags (i.e. -mtune=haswell on MacOS M1!)
# If is safer to completely unset these flags
unset(ENV{CXXFLAGS})
unset(ENV{DEBUG_CXXFLAGS})
unset(ENV{CFLAGS})
unset(ENV{DEBUG_CFLAGS})

# name the project
project (SIRE CXX C)

# make sure that libraries are not found in any system directories
# This is because everything will be linked to from the conda directory
set (CMAKE_IGNORE_PATH)
foreach (_prefix, /usr /opt /usr/local /usr/local/Frameworks )
  list(APPEND CMAKE_IGNORE_PATH ${_prefix}/include ${_prefix}/lib ${_prefix}/lib64)
endforeach()

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
save_sire_variable( "SIRE_VERSION" "${SIRE_VERSION}" )
save_sire_variable( "SIRE_VERSION_STRING" "${SIRE_VERSION_STRING}" )
save_sire_variable( "SIRE_APP_DIR" "${CMAKE_INSTALL_PREFIX}" )

save_sire_variable( "BUILD_NCORES" "${NCORES}" )
save_sire_variable( "SIRE_ANACONDA_BASE" "${ANACONDA_BASE}" )

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
# so that they can find each other when loading - we will manually handle RPATH
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
  set (SIRE_INSTALL_RPATH "@executable_path/.:@executable_path/../${SIRE_LIBS}")
else()
  set (SIRE_INSTALL_RPATH "\$ORIGIN/.:\$ORIGIN/../${SIRE_LIBS}")
endif()

set (CMAKE_INSTALL_RPATH "${SIRE_INSTALL_RPATH}")
set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

save_sire_variable( "SIRE_INSTALL_PREFIX" "${CMAKE_INSTALL_PREFIX}" )
save_sire_variable( "SIRE_BIN" "${SIRE_BIN}" )
save_sire_variable( "SIRE_LIBS" "${SIRE_LIBS}" )
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
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/build/cmake)

# Options to pass to 'strip' when stripping the libraries
# (on Linux, AIX and MacOS need to use '-x' to ensure we leave
# the global symbols)
set ( SIRE_STRIP_OPTIONS "-x" CACHE STRING "Options to pass to 'strip'" )
save_sire_variable( "SIRE_STRIP_OPTIONS" "${SIRE_STRIP_OPTIONS}" )

set ( SIRE_STRIP_COMMAND "${CMAKE_STRIP} ${SIRE_STRIP_OPTIONS}" )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_strip.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip"
  IMMEDIATE @ONLY)

set ( CMAKE_STRIP "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip" )

# Add option to turn on or off symbol hiding - it is really needed
# to keep the symbol table manageable
option ( SIRE_SYMBOL_HIDING "Turn on library symbol hiding" ON )
save_sire_variable( "SIRE_SYMBOL_HIDING" "${SIRE_SYMBOL_HIDING}" )

# Add option to turn on or off vectorisation of the code
option ( SIRE_VECTORISE "Turn on code vectorisation" ON )
save_sire_variable( "SIRE_VECTORISE" "${SIRE_VECTORISE}" )

# Add option to disable use of SSE
option ( SIRE_DISABLE_SSE "Prevent use of SSE when compiling" OFF )

# Add option to disable use of AVX
option ( SIRE_DISABLE_AVX "Prevent use of AVX when compiling" OFF )

# Add option to disable use of AVX512F - this should almost always
# be disabled as it doesn't improve performance and can be buggy
option ( SIRE_DISABLE_AVX512F "Prevent use of AVX512F when compiling" ON )

# Add an option to say whether or not to use OpenMM - this is installed
# via conda and should be available and usable
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

# Create config.h
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/sire_config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/sire_config.h)

# We need to get the compiler flags for different options
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

# Need to define these symbols for the tests below
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
  CHECK_CXX_COMPILER_FLAG( "/std:c++14" SIRE_HAS_CPP_14 )
else()
  CHECK_CXX_COMPILER_FLAG( "-std=c++14" SIRE_HAS_CPP_14 )
endif()

save_sire_variable( "SIRE_HAS_CPP_LIB" "${SIRE_HAS_CPP_LIB}" )

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
else()
  message( FATAL_ERROR "Cannot compile Sire as we now require a modern C++ "
                       "compiler that supports at least C++ 2014 (-std=c++14)")
endif()

save_sire_variable( "SIRE_CXX_NG_FLAGS" "${SIRE_CXX_NG_FLAGS}" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SIRE_CXX_NG_FLAGS}")

set( SIRE_COMPILER "UNKNOWN" )

if (NOT ${CMAKE_C_COMPILER_ID} STREQUAL ${CMAKE_CXX_COMPILER_ID})
  message( WARNING "It is not sensible to compile Sire using different C and C++ compilers. "
                   "You are using ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID}) and "
                   "${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
endif()

string( TOUPPER ${CMAKE_CXX_COMPILER_ID} SIRE_CXX_COMPILER_ID )

# Brings in the functions used below to find the vectorisation flags
include (FindVector)

# Add compiler specific flags
if ( ${SIRE_CXX_COMPILER_ID} MATCHES "GNU" )

    # This is the default for conda on Linux

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

    # First, we need at least GCC 7, (released 2018). GCC 11 comes with conda
    if ( GCC_MAJOR_VERSION LESS 7 )
      message( FATAL_ERROR "Sire requires GCC >= 7.0. Please upgrade "
                           "your version of GCC." )
    endif()

    set ( SIRE_SMALL_FLAGS  "-Oz" )
    set ( SIRE_WARNALL_FLAGS  "-Wall" )
    set ( SIRE_DEBUG_FLAGS   "-g" )

    set ( SIRE_RELEASE_FLAGS "-O3 -ffast-math" )

    # -fomit-frame-pointer breaks backtrace on an apple
    set ( SIRE_RELEASE_FLAGS "${SIRE_RELEASE_FLAGS} -fomit-frame-pointer" )

    set ( SIRE_PLATFORM_FLAGS "-pipe -DSIRE_ALWAYS_INLINE=inline" )
    set ( SIRE_STATIC_LINK_FLAGS "-static" )

    GET_SIRE_VECTOR_FLAGS( "-fopenmp-simd" "-msse2" "-mavx" "-mavx512f" )

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

    # This is the default for conda on MacOS

    message( STATUS "Detected and using clang C and C++ compilers" )
    set( SIRE_COMPILER "CLANG" )

    set ( SIRE_SMALL_FLAGS  "-Oz" )
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

    # This is not tested as it is not the default with conda.
    # Use at your own risk

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

    # This is the default for conda on Windows

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

# Tell CMake to not run moc when necessary (we don't use moc or the QObject system)
set(CMAKE_AUTOMOC OFF)

# Always useful to look for header files in the current directory
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set (SIRE_FOUND_TBB FALSE)
set (SIRE_FOUND_GSL FALSE)
set (SIRE_FOUND_BOOST FALSE)
set (SIRE_FOUND_BLASLAPACK FALSE)
set (SIRE_FOUND_QT FALSE)
set (SIRE_FOUND_OPENMM FALSE)
set (SIRE_FOUND_NETCDF FALSE)

if ( NOT ${SIRE_FOUND_QT} )
  # look for Qt5 in the system - this should be found from conda
  find_package( Qt5Core REQUIRED )
endif()

get_target_property(QT5_LOCATION Qt5::Core LOCATION)
get_filename_component(SIRE_QT5CORE_DIR ${QT5_LOCATION} PATH)
message( STATUS "Qt5Core location ${SIRE_QT5CORE_DIR}" )
save_sire_variable( "SIRE_QT5CORE_DIR" "${SIRE_QT5CORE_DIR}" )

if ( NOT ${SIRE_FOUND_NETCDF} )
  # look for netcdf in the system from conda
  set(NetCDF_ROOT_DIR "${ANACONDA_BASE}")
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
    set(OpenMM_ROOT_DIR "${ANACONDA_BASE}")
    set(OPENMM_ROOT_DIR "${ANACONDA_BASE}")
    FIND_PACKAGE( OpenMM )

    save_sire_variable( "SIRE_OpenMM_FOUND" "${OpenMM_FOUND}" )

    if ( ${OpenMM_FOUND} )
      message( STATUS "Found OpenMM at ${OpenMM_ROOT_DIR}" )
      save_sire_variable( "SIRE_OpenMM_INCLUDE_DIR" "${OpenMM_INCLUDE_DIR}" )
      save_sire_variable( "SIRE_OpenMM_LIBRARIES" "${OpenMM_LIBRARIES}" )
      save_sire_variable( "SIRE_OpenMM_LIBRARY_DIR" "${OpenMM_LIBRARY_DIR}" )
    else()
      message( STATUS "Could not find OpenMM. OpenMM support is disabled." )
    endif()
  endif()
else()
  message( STATUS "Disabling use of OpenMM at user request.")
  save_sire_variable( "SIRE_OpenMM_FOUND" "${SIRE_USE_OPENMM}" )
endif()

include_directories( "${SIRE_APP}/include" "${SIRE_APP}/Library/include" )

if (NOT ${SIRE_FOUND_BOOST})
  # Sire depends on >= boost 1.7 - all libraries must be dynamically linked
  set (Boost_USE_STATIC_LIBS OFF)
  set (Boost_USE_STATIC_RUNTIME OFF)
  set (Boost_DEBUG ON)  # this gives a nice print out in the cmake list to debug how boost is found and linked
  find_package( Boost 1.7 REQUIRED )

  if (Boost_FOUND)
    message(STATUS "Boost paths ${Boost_LIBRARIES} | ${Boost_INCLUDE_DIR}" )
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
  find_package(GSL 1.6 QUIET REQUIRED)

  message( STATUS "GSL paths ${GSL_LINK_DIRECTORIES} "
                  "${GSL_LIBRARIES} | ${GSL_INCLUDE_DIR}" )

  list(APPEND SIREMATHS_EXTRA_LIBRARIES ${GSL_LIBRARIES})

  include_directories( BEFORE ${GSL_INCLUDE_DIR})

  save_sire_variable( "SIRE_GSL_INCLUDE_PATH" "${GSL_INCLUDE_DIR}" )
endif()

set( SIRE_USE_BLAS OFF )
set( SIRE_USE_LAPACK OFF )
save_sire_variable( "SIRE_USE_BLAS" "${SIRE_USE_BLAS}" )
save_sire_variable( "SIRE_USE_LAPACK" "${SIRE_USE_LAPACK}" )

if (NOT ${SIRE_FOUND_TBB})
  set(TBB_LIBRARY_DIR "${ANACONDA_BASE}/lib")

  message( STATUS "TBB LIBRARY DIR: ${TBB_LIBRARY_DIR}")
  find_library( TBB_LIBRARY "tbb" PATHS ${TBB_LIBRARY_DIR} REQUIRED )
  find_library( TBB_MALLOC_LIBRARY "tbbmalloc" PATHS ${TBB_LIBRARY_DIR} REQUIRED )

  message( STATUS "TBB paths ${TBB_LIBRARY} | ${TBB_MALLOC_LIBRARY} | ${TBB_INCLUDE_DIR}" )

  include_directories(${TBB_INCLUDE_DIR})
  save_sire_variable( "TBB_LIBRARY" "${TBB_LIBRARY}" )
  save_sire_variable( "TBB_MALLOC_LIBRARY" "${TBB_MALLOC_LIBRARY}" )
  save_sire_variable( "TBB_INCLUDE_DIR" "${TBB_INCLUDE_DIR}" )
endif()

# Thanks to this project who have a solution to the tbb12 link
# issue on Windows
# https://fossies.org/linux/opencv/cmake/OpenCVFindLAPACK.cmake
if(CMAKE_GENERATOR MATCHES "Visual Studio")  # MSBuild
  # workaround DEFAULTLIB:tbb12.lib issue
  # We may have tbb.lib for 'tbb' target, but not 'tbb12.lib'
  message(STATUS "Adding tbb12 to ignore list")
  set(SIRE_IGNORE_DEFAULTLIB_TBB "tbb12.lib")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIB:${SIRE_IGNORE_DEFAULTLIB_TBB}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:${SIRE_IGNORE_DEFAULTLIB_TBB}")
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

message ( STATUS "SIRE_CXX_NG_FLAGS      == ${SIRE_CXX_NG_FLAGS}" )
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

message(STATUS "SIRE_COMPILE_FLAGS ${SIRE_COMPILE_FLAGS}")
message(STATUS "SIRE_CXX_NG_FLAGS ${SIRE_CXX_NG_FLAGS}")
message(STATUS "SIRE_SMALL_FLAGS ${SIRE_SMALL_FLAGS}")
message(STATUS "SIRE_USER_CPPFLAGS ${SIRE_USER_CPPFLAGS}")

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
set( CMAKE_C_FLAGS_MINSIZEREL "${SIRE_C_FLAGS_SMALL}" )
set( CMAKE_CXX_FLAGS_MINSIZEREL "${SIRE_CXX_FLAGS_SMALL}" )

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

#�Now save the compile variables so that they can be re-used by libraries that use Sire
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

# Include the build directory so that we can get sire_config.h and sire_version.h
include_directories(${CMAKE_BINARY_DIR})

################
################ TELLING CMAKE WHERE THE SOURCE IS LOCATED
################

if (NOT DEFINED ${SIRE_VERSION})
  message( STATUS "Strange, SIRE_VERSION has been unset. Resetting it.")
  set (SIRE_VERSION "${S_VERSION_MAJOR}.${S_VERSION_MINOR}.${S_VERSION_PATCH}")
endif()

message( STATUS "SIRE VERSION IS ${SIRE_VERSION}")

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
        )

# install the molecular templates into share/Sire
install (DIRECTORY ${CMAKE_SOURCE_DIR}/templates
         DESTINATION ${SIRE_SHARE}
         USE_SOURCE_PERMISSIONS
        )

install (DIRECTORY ${CMAKE_SOURCE_DIR}/share/build
         DESTINATION ${SIRE_SHARE}
         USE_SOURCE_PERMISSIONS
        )

# These commands are used to create an 'uninstall' target
# (this is copied from the vtk.org wiki)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
