# - Find the Sire libraries
# Find the native Sire includes and libraries
#
# Sire uses install(EXPORT) to create native CMake files that
# allow you to link against Sire targets as if they were in 
# your project. The targets are prefixed with "SIRE_", so
# if you want to link against SireBase, just add "SIRE_SireBase"
# to your target_link_libraries.
#
# SIRE_INCLUDE_DIR contains the Sire root directory, and
# SIRE_EXTRA_INCLUDES contains the paths to the header files
# of the dependencies of Sire (GSL, Boost and MPI)
#
#  SIRE_INSTALL_PREFIX - the installation prefix to the installed corelib
#  SIRE_INCLUDE_DIR    - where to find the Sire include files
#  SIRE_EXTRA_INCLUDES - where to find the header files of the dependencies
#  SIRE_FOUND          - True if Sire found.

# I prefer to use all capitals for variables
set (SIRE_FIND_VERSION_MAJOR ${Sire_FIND_VERSION_MAJOR})
set (SIRE_FIND_VERSION_MINOR ${Sire_FIND_VERSION_MINOR})
set (SIRE_FIND_VERSION_PATCH ${Sire_FIND_VERSION_PATCH})
set (SIRE_ROOT ${Sire_ROOT})
set (SIRE_INCLUDEDIR ${Sire_INCLUDEDIR})
set (SIRE_FIND_REQUIRED ${Sire_FIND_REQUIRED})

# If SIRE_INSTALL_PREFIX was defined in the environment, use it.
if (NOT SIRE_INSTALL_PREFIX AND NOT $ENV{SIRE_INSTALL_PREFIX} STREQUAL "")
    set(SIRE_INSTALL_PREFIX $ENV{SIRE_INSTALL_PREFIX})
endif(NOT SIRE_INSTALL_PREFIX AND NOT $ENV{SIRE_INSTALL_PREFIX} STREQUAL "")

# If SIRE_INCLUDEDIR was defined in the environment, use it.
if ( NOT $ENV{SIRE_INCLUDEDIR} STREQUAL "" )
    set(SIRE_INCLUDEDIR $ENV{SIRE_INCLUDEDIR})
endif( NOT $ENV{SIRE_INCLUDEDIR} STREQUAL "" )

if ( SIRE_INSTALL_PREFIX )
    file(TO_CMAKE_PATH ${SIRE_INSTALL_PREFIX} SIRE_INSTALL_PREFIX)
    set(_SIRE_INCLUDE_SEARCH_DIRS 
                 ${SIRE_INSTALL_PREFIX}/include 
                 ${SIRE_INSTALL_PREFIX}
                 ${_SIRE_INCLUDE_SEARCH_DIRS})
endif ( SIRE_INSTALL_PREFIX )

if ( SIRE_INCLUDEDIR )
    file(TO_CMAKE_PATH ${SIRE_INCLUDEDIR} SIRE_INCLUDEDIR)
    set(_SIRE_INCLUDE_SEARCH_DIRS 
                 ${SIRE_INCLUDEDIR} ${_SIRE_INCLUDE_SEARCH_DIRS})
endif ( SIRE_INCLUDEDIR )

message(STATUS "${SIRE_INCLUDE_DIR}")

#Now find the sire_version.h file - this is in the same directory
#as the cmake directory that contains the cmake files necessary
#to import Sire
find_path ( SIRE_INCLUDE_DIR
              NAMES sire_version.h
              HINTS ${_SIRE_INCLUDE_SEARCH_DIRS} ~/local/include
              PATH_SUFFIXES Sire 
          )

set (SIRE_FOUND false)

if (SIRE_INCLUDE_DIR)
    #read the version file to ensure that we have the right version of Sire
    #Read the whole file
    set (SIRE_VERSION 0)
    set (SIRE_LIB_VERSION "")

    if ( NOT EXISTS "${SIRE_INCLUDE_DIR}/sire_version.h" )
        message( FATAL_ERROR "Cannot find ${SIRE_INCLUDE_DIR}/sire_version.h - has "
                             "Sire been removed since you last ran cmake? Is it installed properly?")
    endif()

    file (READ "${SIRE_INCLUDE_DIR}/sire_version.h" _SIRE_VERSION_H_CONTENTS)

    string(REGEX REPLACE ".*#define SIRE_VERSION[ ]+([0-9]+).*" 
                         "\\1" SIRE_VERSION "${_SIRE_VERSION_H_CONTENTS}")

    string(REGEX REPLACE ".*#define SIRE_LIB_VERSION[ ]+([0-9_]+).*"
                         "\\1" SIRE_LIB_VERSION "${_SIRE_VERSION_H_CONTENTS}")

    set (SIRE_LIB_VERSION ${SIRE_LIB_VERSION} CACHE INTERNAL 
                             "The library version string for Sire libraries")
    set (SIRE_VERSION ${SIRE_VERSION} CACHE INTERNAL "The version number for Sire libraries")
    
    if ( NOT "${SIRE_VERSION}" STREQUAL "0")
        math(EXPR SIRE_MAJOR_VERSION "${SIRE_VERSION} / 100000")
        math(EXPR SIRE_MINOR_VERSION "${SIRE_VERSION} / 100 % 1000")
        math(EXPR SIRE_PATCH_VERSION "${SIRE_VERSION} % 100")

        set (SIRE_VERSION_MAJOR ${SIRE_MAJOR_VERSION})
        set (SIRE_VERSION_MINOR ${SIRE_MINOR_VERSION})
        set (SIRE_VERSION_PATCH ${SIRE_PATCH_VERSION})

        set(SIRE_ERROR_REASON  "${SIRE_ERROR_REASON}Sire version: "
                               "${SIRE_MAJOR_VERSION}.${SIRE_MINOR_VERSION}.${SIRE_PATCH_VERSION}\n"
                               "Sire include path: ${SIRE_INCLUDE_DIR}")
    endif (NOT "${SIRE_VERSION}" STREQUAL "0")

else()
    set (SIRE_ERROR_REASON "Could not find the file 'sire_version.h'. Try setting "
                           "SIRE_INCLUDE_DIR to the directory that contains this file, "
                           "or setting SIRE_ROOT to the directory in which Sire is installed.")
endif()

if (SIRE_INCLUDE_DIR)
    SET( SIRE_FOUND TRUE )

    # Check the version of Sire against the requested version.
    if (SIRE_FIND_VERSION AND NOT SIRE_FIND_VERSION_MINOR)
        message(SEND_ERROR "When requesting a specific version of Sire, you must provide "
                           "at least the major and minor version numbers, e.g., 1.2")
    endif (SIRE_FIND_VERSION AND NOT SIRE_FIND_VERSION_MINOR)

    if (SIRE_MAJOR_VERSION LESS "${SIRE_FIND_VERSION_MAJOR}" )
        set( SIRE_FOUND FALSE )
        set( _SIRE_VERSION_AGE "old")
    
    elseif (SIRE_MAJOR_VERSION EQUAL "${SIRE_FIND_VERSION_MAJOR}" )
        if (SIRE_MINOR_VERSION LESS "${SIRE_FIND_VERSION_MINOR}" )
            set( SIRE_FOUND FALSE )
            set( _SIRE_VERSION_AGE "old")
      
        elseif (SIRE_MINOR_VERSION EQUAL "${SIRE_FIND_VERSION_MINOR}" )
            if( SIRE_FIND_VERSION_PATCH AND SIRE_PATCH_VERSION LESS "${SIRE_FIND_VERSION_PATCH}" )
                set( SIRE_FOUND FALSE )
                set( _SIRE_VERSION_AGE "old")
            endif( SIRE_FIND_VERSION_PATCH AND SIRE_PATCH_VERSION LESS "${SIRE_FIND_VERSION_PATCH}" )
        endif( SIRE_MINOR_VERSION LESS "${SIRE_FIND_VERSION_MINOR}" )
    endif( SIRE_MAJOR_VERSION LESS "${SIRE_FIND_VERSION_MAJOR}" )

    if (SIRE_FOUND AND SIRE_FIND_VERSION_EXACT)
        # If the user requested an exact version of Sire, check
        # that. We already know that the Sire version we have is >= the
        # requested version.
        set( _SIRE_VERSION_AGE "new")

        # If the user didn't specify a patchlevel, it's 0.
        if (NOT SIRE_FIND_VERSION_PATCH)
            set(SIRE_FIND_VERSION_PATCH 0)
        endif (NOT SIRE_FIND_VERSION_PATCH)
      
        # We'll set SIRE_FOUND true again if we have an exact version match.
        set(SIRE_FOUND FALSE)
        if (SIRE_MAJOR_VERSION EQUAL "${SIRE_FIND_VERSION_MAJOR}" )
            if (SIRE_MINOR_VERSION EQUAL "${SIRE_FIND_VERSION_MINOR}" )
                if(SIRE_PATCH_VERSION EQUAL "${SIRE_FIND_VERSION_PATCH}" )
                    set( SIRE_FOUND TRUE )
                endif(SIRE_PATCH_VERSION EQUAL "${SIRE_FIND_VERSION_PATCH}" )
            endif( SIRE_MINOR_VERSION EQUAL "${SIRE_FIND_VERSION_MINOR}" )
        endif( SIRE_MAJOR_VERSION EQUAL "${SIRE_FIND_VERSION_MAJOR}" )
    endif (SIRE_FOUND AND SIRE_FIND_VERSION_EXACT)

    if(NOT SIRE_FOUND)
        # State that we found a version of Sire that is too new or too old.
        set(SIRE_ERROR_REASON
            "${SIRE_ERROR_REASON}\nDetected version of Sire is too ${_SIRE_VERSION_AGE}. "
            "Requested version was ${SIRE_FIND_VERSION_MAJOR}.${SIRE_FIND_VERSION_MINOR}")
        if (SIRE_FIND_VERSION_PATCH)
            set(SIRE_ERROR_REASON 
                "${SIRE_ERROR_REASON}.${SIRE_FIND_VERSION_PATCH}")
        endif (SIRE_FIND_VERSION_PATCH)
      
        if (NOT SIRE_FIND_VERSION_EXACT)
            set(SIRE_ERROR_REASON "${SIRE_ERROR_REASON} (or newer)")
        endif (NOT SIRE_FIND_VERSION_EXACT)
     
        set(SIRE_ERROR_REASON "${SIRE_ERROR_REASON}.\n"
                  "Try setting SIRE_ROOT or SIRE_INCLUDE_DIR to point\n"
                  "to an installation of Sire that has the right version.\n")
    endif (NOT SIRE_FOUND)
endif (SIRE_INCLUDE_DIR)

set (Sire_FOUND "${SIRE_FOUND}")

if (SIRE_FOUND)
    if (NOT SIRE_FIND_QUIETLY)
        message(STATUS "Sire version: ${SIRE_MAJOR_VERSION}.${SIRE_MINOR_VERSION}.${SIRE_PATCH_VERSION}")
        message(STATUS "Importing Sire library definitions from "
                       "${SIRE_INCLUDE_DIR}/cmake/SireLibraries.cmake")
    endif (NOT SIRE_FIND_QUIETLY)

    set (SIRE_VERSION "${SIRE_VERSION_MAJOR}.${SIRE_VERSION_MINOR}.${SIRE_VERSION_PATCH}")

    set ( SIRE_LIB_CMAKE "${SIRE_INCLUDE_DIR}/cmake/SireLibraries.cmake" )
    set ( SIRE_COMP_CMAKE "${SIRE_INCLUDE_DIR}/cmake/SireCompileVariables.cmake" )

    if ( NOT EXISTS "${SIRE_LIB_CMAKE}" )
        message(FATAL_ERROR "Cannot find ${SIRE_LIB_CMAKE} - has Sire been installed correctly?")
    endif()

    if ( NOT EXISTS "${SIRE_COMP_CMAKE}" )
        message(FATAL_ERROR "Cannot find ${SIRE_COMP_CMAKE} - has Sire been installed correctly?")
    endif()

    include ("${SIRE_INCLUDE_DIR}/cmake/SireLibraries.cmake")
    include ("${SIRE_INCLUDE_DIR}/cmake/SireCompileVariables.cmake")

else (SIRE_FOUND)
    if (SIRE_FIND_REQUIRED)
        foreach(arg ${SIRE_ERROR_REASON})
            set(_REASON_STRING "${_REASON_STRING} ${arg}")
        endforeach(arg)

        message(SEND_ERROR "Unable to find the requested Sire libraries.\n${_REASON_STRING}")
    endif(SIRE_FIND_REQUIRED)
endif (SIRE_FOUND)
