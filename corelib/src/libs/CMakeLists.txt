###############################
#
# CMake file to compile the
# Sire libraries
#
###############################

# List of each subdirectory
set ( SUBDIRS SireAnalysis
              SireBase
              SireCAS
              SireCluster
              SireError
              SireFF
              SireID
              SireIO
              SireMaths
              SireMM
              SireMol
              SireMove
              SireSearch
              SireStream
              SireSystem
              SireUnits
              SireVol
              Squire
            )

foreach (DIR ${SUBDIRS})
    add_subdirectory( ${DIR} )
endforeach (DIR)

set ( SIRE_HEADERS
      sireglobal.h
      tostring.h
    )

install( FILES ${SIRE_HEADERS}
         DESTINATION ${SIRE_INCLUDES} )

install( EXPORT SireLibraries
         NAMESPACE SIRE_
         DESTINATION ${SIRE_CMAKEFILES} )
