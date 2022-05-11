########################################
#
# CMake file for Python wrapper:_Search
#
########################################

# Third Party dependencies of this module
include_directories( ${BOOST_INCLUDE_DIR}
                     ${PYTHON_INCLUDE_DIR} )

# Other Sire libraries
include_directories( ${CMAKE_SOURCE_DIR} )

include_directories( BEFORE ${SIRE_INCLUDE_DIR} )

# Define the sources in SireSearch
set ( WRAPSEARCH_SOURCES

      _Search.main.cpp

      ${PYPP_SOURCES}
    )

# Create the library that holds all of the wrappers
add_library (Search ${WRAPSEARCH_SOURCES})

target_link_libraries (Search
                       SirePython
                       SIRE_SireBase
                       SIRE_SireMol
                       SIRE_SireMM
                       SIRE_SireSearch
                       ${Boost_LIBRARIES}
                       ${PYTHON_LIBRARIES}
                      )


set_target_properties (Search
                       PROPERTIES
                         VERSION ${SIRE_VERSION}
                         SOVERSION ${SIRE_VERSION_MAJOR}
                         OUTPUT_NAME "Search"
                         PREFIX "_"
                         SUFFIX ${PYTHON_MODULE_EXTENSION}
                      )

include( LimitSirePythonExportSymbols )
export_this_symbol_only( "init_Search" "_init_Search" )

# installation
set( INSTALLDIR ${SIRE_PYTHON}/sire/legacy/Search )

install( TARGETS Search
         LIBRARY DESTINATION ${INSTALLDIR}
         ARCHIVE DESTINATION ${INSTALLDIR}
         RUNTIME DESTINATION ${INSTALLDIR}
       )

install( FILES __init__.py
         DESTINATION ${INSTALLDIR}
       )
