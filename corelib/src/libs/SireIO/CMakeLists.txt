########################################
#
# CMake file for library:SireIO
#
########################################

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# This library uses Intel Threaded Building blocks
include_directories(${TBB_INCLUDE_DIR})

# Define the headers in SireIO
set ( SIREIO_HEADERS
      amber.h
      amberformat.h
      amberprm.h
      amberrst7.h
      cube.h
      errors.h
      flexibilitylibrary.h
      iobase.h
      moleculeparser.h
      pdb.h
      perturbationslibrary.h
      protoms.h
      tinker.h
      trajectorymonitor.h
      zmatrixmaker.h
    )

# Define the sources in SireIO
set ( SIREIO_SOURCES

      register_sireio.cpp

      amber.cpp
      amberformat.cpp
      amberprm.cpp
      amberrst7.cpp
      cube.cpp
      errors.cpp
      flexibilitylibrary.cpp
      iobase.cpp
      moleculeparser.cpp
      pdb.cpp
      perturbationslibrary.cpp
      protoms.cpp
      tinker.cpp
      trajectorymonitor.cpp    
      zmatrixmaker.cpp

      ${SIREIO_HEADERS}
    )

add_library (SireIO ${SIREIO_SOURCES})

set_target_properties (SireIO PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireIO
                       SireMM
                       SireMove
                       SireMol
                       SireStream
                       )

# installation
install( TARGETS SireIO EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREIO_HEADERS} 
         DESTINATION ${SIRE_INCLUDES}/SireIO )
