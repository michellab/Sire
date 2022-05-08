########################################
#
# CMake file for library:SireIO
#
########################################

# Export symbols when building the library
add_definitions( "-DSIRESEARCH_BUILD" )

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# This library uses Intel Threaded Building blocks
include_directories(${TBB_INCLUDE_DIR})

# Define the headers in SireSearch
set ( SIRESEARCH_HEADERS
      ast.h
      grammar.h
      idengine.h
      parser.h
    )

# Define the sources in SireSearch
set ( SIRESEARCH_SOURCES

      register_siresearch.cpp

      ast.cpp
      idengine.cpp
      parser.cpp

      ${SIRESEARCH_HEADERS}
    )

add_library (SireSearch ${SIRESEARCH_SOURCES})

set_target_properties (SireSearch PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireSearch
                       SireMM
                       SireMol
                       SireCAS
                       SireStream
                       )

# installation
install( TARGETS SireSearch EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIRESEARCH_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireSearch )
