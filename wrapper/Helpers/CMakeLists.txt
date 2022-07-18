########################################
#
# CMake file for Python Helper Library
#   (SirePython)
#
########################################

# Third Party dependencies of this module
include_directories( ${BOOST_INCLUDE_DIR}
                     ${PYTHON_INCLUDE_DIR} )

# Sire include paths
include_directories( BEFORE ${SIRE_INCLUDE_DIR} )

# Other python wrapping directories
include_directories(${CMAKE_SOURCE_DIR})

# Define the sources in Helpers
set ( WRAPHELPERS_SOURCES

      objectregistry.hpp
      objectregistry.cpp

      release_gil_policy.hpp
      release_gil_policy.cpp

    )

# Create the library that holds all of the class wrappers
add_library (SirePython ${WRAPHELPERS_SOURCES})

set_target_properties (SirePython
                       PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${SIRE_VERSION_MAJOR}
		       )

target_link_libraries (SirePython
                       SIRE_SireCluster
                       SIRE_SireStream
                       SIRE_SireError
                       ${Boost_LIBRARIES}
                       ${PYTHON_LIBRARIES}
                      )

# installation
install( TARGETS SirePython
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )
