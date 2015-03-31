########################################
#
# CMake file for library:SireFF
#
########################################

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Define the headers in SireFF
set ( SIREFF_HEADERS
      atomicffparameters.hpp
      energytable.h
      errors.h
      ff.h
      ff3d.h
      ffcomponent.h
      ffid.h
      ffidentifier.h
      ffidx.h
      ffmolgroup.h
      ffname.h
      ffparameters.h
      fieldtable.h
      forcefield.h
      forcefields.h
      forcetable.h
      g1ff.h
      g2ff.h
      inter2b2g3dff.hpp
      inter2b2gff.hpp
      inter2b3dff.hpp
      inter2bff.hpp
      intra2bff.hpp
      intra2b2gff.hpp
      intra2b3dff.hpp
      intra2b2g3dff.hpp
      patch.h
      patches.h
      point.h
      potentialtable.h
      probe.h
    )

set ( SIREFF_DETAIL_HEADERS
      detail/atomiccoords3d.h
      detail/atomicparameters.hpp
      detail/atomicparameters3d.hpp
      detail/ffmolecules.h
      detail/ffmolecules3d.h
    )

# Define the sources in SireFF
set ( SIREFF_SOURCES

      register_sireff.cpp

      energytable.cpp
      errors.cpp
      ff.cpp
      ff3d.cpp
      ffcomponent.cpp
      ffid.cpp
      ffidentifier.cpp
      ffmolgroup.cpp
      ffparameters.cpp
      fieldtable.cpp
      forcefield.cpp
      forcefields.cpp
      forcetable.cpp
      g1ff.cpp
      g2ff.cpp
      patch.cpp
      patches.cpp
      point.cpp
      potentialtable.cpp
      probe.cpp

      detail/atomiccoords3d.cpp
      detail/ffmolecules.cpp
      detail/ffmolecules3d.cpp   

      ${SIREFF_HEADERS}
      ${SIREFF_DETAIL_HEADERS}
    )

add_library (SireFF ${SIREFF_SOURCES})

set_target_properties (SireFF PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${SIRE_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireFF
                       SireMol
                       SireCAS
                       SireMaths
                       SireBase
                       SireStream
                       )

# installation
install( TARGETS SireFF EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREFF_HEADERS} 
         DESTINATION ${SIRE_INCLUDES}/SireFF )

install( FILES ${SIREFF_DETAIL_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireFF/detail )

