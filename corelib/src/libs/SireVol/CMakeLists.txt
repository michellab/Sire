########################################
#
# CMake file for library:SireVol
#
########################################

# Export symbols when building the library
add_definitions( "-DSIREVOL_BUILD" )

# Third Party dependencies of this module

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Define the headers in SireVol
set ( SIREVOL_HEADERS
      aabox.h
      cartesian.h
      combinedspace.h
      combinespaces.h
      coordgroup.h
      errors.h
      grid.h
      gridinfo.h
      patching.h
      periodicbox.h
      space.h
    )

# Define the sources in SireVol
set ( SIREVOL_SOURCES

      register_sirevol.cpp

      aabox.cpp
      cartesian.cpp
      combinedspace.cpp
      combinespaces.cpp
      coordgroup.cpp
      errors.cpp
      grid.cpp
      gridinfo.cpp
      patching.cpp
      periodicbox.cpp
      space.cpp

      ${SIREVOL_HEADERS}
    )

add_library (SireVol ${SIREVOL_SOURCES})

set_target_properties (SireVol PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireVol
                       SireMaths
                       SireBase
                       SireStream
                       )

# installation
install( TARGETS SireVol EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREVOL_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireVol )
