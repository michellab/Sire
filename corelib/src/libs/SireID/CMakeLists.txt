########################################
#
# CMake file for library:SireID
#
########################################

# Third Party dependencies of this module

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Define the headers in SireID
set ( SIREID_HEADERS
      id.h
      idandset.hpp
      identifier.h
      idorset.hpp
      index.h
      name.h
      number.h
      specify.hpp
    )

# Define the sources in SireID
set ( SIREID_SOURCES

      register_sireid.cpp

      id.cpp
      identifier.cpp
      index.cpp

      ${SIREID_HEADERS}
    )

add_library (SireID ${SIREID_SOURCES})

set_target_properties (SireID PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireID
                       SireStream
                       )

# installation
install( TARGETS SireID EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREID_HEADERS} 
         DESTINATION ${SIRE_INCLUDES}/SireID )

