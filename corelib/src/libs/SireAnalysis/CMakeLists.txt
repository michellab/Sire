########################################
#
# CMake file for library:SireAnalysis
#
########################################

# Third Party dependencies of this module

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Define the headers in SireAnalysis
set ( SIREANALYSIS_HEADERS
      bennetts.h
      fep.h
      ti.h
      ticomponents.h
    )


# Define the sources in SireAnalysis
set ( SIREANALYSIS_SOURCES

      register_sireanalysis.cpp

      bennetts.cpp
      fep.cpp
      ti.cpp
      ticomponents.cpp

      third_party/regress.cpp

      ${SIREANALYSIS_HEADERS}
    )

add_library (SireAnalysis ${SIREANALYSIS_SOURCES})

set_target_properties (SireAnalysis PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${SIRE_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireAnalysis
                       SireSystem
                       SireCAS
                       SireMaths
                       SireBase
                       SireStream
                       )

# installation
install( TARGETS SireAnalysis EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREANALYSIS_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireAnalysis )

