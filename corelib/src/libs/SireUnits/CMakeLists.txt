########################################
#
# CMake file for library:SireUnits
#
########################################

# Third Party dependencies of this module

# Other Sire libraries
include_directories( ${CMAKE_SOURCE_DIR}/src/libs )

set ( SIREUNITS_HEADERS
      angle.h
      convert.h
      dimensions.h
      temperature.h
      units.h
    )

# Define the sources in SireUnits
set ( SIREUNITS_SOURCES

      register_sireunits.cpp

      angle.cpp
      convert.cpp            
    
      ${SIREUNITS_HEADERS}
    )

add_library (SireUnits ${SIREUNITS_SOURCES})

set_target_properties (SireUnits PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireUnits
                       SireStream
                       )

# installation
install( TARGETS SireUnits EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREUNITS_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireUnits )
