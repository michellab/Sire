########################################
#
# CMake file for library:Squire
#
########################################

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

set ( SQUIRE_HEADERS

      am1bcc.h
      qmchargecalculator.h
      qmchargeconstraint.h
      qmff.h
      qmmmelecembedpotential.h
      qmmmff.h
      qmmmpotential.h
      qmpotential.h
      qmprogram.h

      gto.h
      hf.h
      latticecharges.h

      orbital.h
      pointcharge.h
      pointdipole.h

      molpro.h
      mopac.h
      sqm.h

      pgto.h
      sgto.h
    )

set ( SQUIRE_SOURCES

      register_squire.cpp

      am1bcc.cpp
      qmchargecalculator.cpp
      qmchargeconstraint.cpp
      qmff.cpp
      qmmmelecembedpotential.cpp
      qmmmff.cpp
      qmmmpotential.cpp
      qmpotential.cpp
      qmprogram.cpp

      gto.cpp
      hf.cpp
      latticecharges.cpp
      
      molpro.cpp
      mopac.cpp
      sqm.cpp

      orbital.cpp
      pointcharge.cpp
      pointdipole.cpp

      pgto.cpp
      sgto.cpp

      ${SQUIRE_HEADERS}
    )

#add the additional include files needed by Squire
include_directories( ${SQUIRE_INCLUDE_DIRS} )
  
add_library (Squire ${SQUIRE_SOURCES})

set_target_properties (Squire PROPERTIES
                         VERSION ${SIRE_VERSION}
                         SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (Squire
                       SireIO
                       SireSystem
                       SireMM
                       SireFF
                       SireMol
                       SireCAS
                       SireBase
                       SireStream
                       SireError
                      )

# installation
install( TARGETS Squire EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SQUIRE_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/Squire )
