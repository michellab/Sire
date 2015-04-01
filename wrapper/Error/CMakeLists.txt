########################################
#
# CMake file for Python wrapper:_Error
#
########################################

# Third Party dependencies of this module
include_directories( ${BOOST_INCLUDE_DIR}
                     ${PYTHON_INCLUDE_DIR} )

# Other Sire libraries
include_directories( BEFORE ${SIRE_INCLUDE_DIR}/SireError
                            ${SIRE_INCLUDE_DIR})

# Define the sources in SireError
set ( ERROR_SOURCES

      wrap_SireError.cpp

      wrap_exceptions.cpp

    )

# Create the library that holds all of the class wrappers
add_library (Error ${ERROR_SOURCES})

set_target_properties (Error 
                       PROPERTIES
                         VERSION ${SIRE_VERSION}
		             SOVERSION ${SIRE_VERSION_MAJOR}
                         OUTPUT_NAME "Error"
                         PREFIX "_"
                         SUFFIX ${PYTHON_MODULE_EXTENSION}
		          )
		       
target_link_libraries (Error
                       SirePython
                       SIRE_SireError
                       ${BOOST_PYTHON_LIBRARY}
                       ${PYTHON_LIBRARIES}
                      )


include( LimitSirePythonExportSymbols )
export_this_symbol_only( "init_Error" "_init_Error" )

# installation
set( INSTALLDIR ${SIRE_PYTHON}/Sire/Error )

install( TARGETS Error
         LIBRARY DESTINATION ${INSTALLDIR}
       )

install( FILES __init__.py
         DESTINATION ${INSTALLDIR} 
       )

