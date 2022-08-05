########################################
#
# CMake file for library:SireMM
#
########################################

# Export symbols when building the library
add_definitions( "-DSIREMM_BUILD" )

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# This library uses Intel Threaded Building blocks
include_directories(${TBB_INCLUDE_DIR})

# Define the headers in SireMM
set ( SIREMM_HEADERS
      amberparams.h
      anglerestraint.h
      angle.h
      atomfunctions.h
      atomljs.h
      atompairs.hpp
      bond.h
      clj14group.h
      cljatoms.h
      cljboxes.h
      cljcalculator.h
      cljcomponent.h
      cljdelta.h
      cljextractor.h
      cljfunction.h
      cljgrid.h
      cljgroup.h
      cljparam.h
      cljpotential.h
      cljprobe.h
      cljnbpairs.h
      cljrffunction.h
      cljshiftfunction.h
      cljworkspace.h
      coulombpotential.h
      dihedralrestraint.h
      distancerestraint.h
      errors.h
      fouratomfunctions.h
      gridff.h
      gridff2.h
      gromacsparams.h
      intercljff.h
      intercoulombff.h
      interff.h
      intergroupff.h
      interljff.h
      internalcomponent.h
      internalff.h
      internalparameters.h
      internalperturbation.h
      intersoftcljff.h
      intracljff.h
      intracoulombff.h
      intraff.h
      intragroupff.h
      intraljff.h
      intrasoftcljff.h
      ljfunction.h
      ljpair.h
      ljparameter.h
      ljparameterdb.h
      ljperturbation.h
      ljpotential.h
      mmdetail.h
      multicljcomponent.h
      restraint.h
      restraintcomponent.h
      restraintff.h
      selectorbond.h
      selectormbond.h
      softcljcomponent.h
      softcljpotential.h
      switchingfunction.h
      testff.h
      threeatomfunctions.h
      twoatomfunctions.h
    )

set ( SIREMM_DETAIL_HEADERS
      detail/intrascaledatomicparameters.hpp
    )

# Define the sources in SireMM
set ( SIREMM_SOURCES

      register_siremm.cpp

      amberparams.cpp
      angle.cpp
      anglerestraint.cpp
      atomfunctions.cpp
      atomljs.cpp
      bond.cpp
      clj14group.cpp
      cljatoms.cpp
      cljboxes.cpp
      cljcalculator.cpp
      cljcomponent.cpp
      cljdelta.cpp
      cljextractor.cpp
      cljfunction.cpp
      cljgrid.cpp
      cljgroup.cpp
      cljparam.cpp
      cljpotential.cpp
      cljprobe.cpp
      cljnbpairs.cpp
      cljrffunction.cpp
      cljshiftfunction.cpp
      cljworkspace.cpp
      coulombpotential.cpp
      dihedral.cpp
      dihedralrestraint.cpp
      distancerestraint.cpp
      errors.cpp
      fouratomfunctions.cpp
      gridff.cpp
      gridff2.cpp
      gromacsparams.cpp
      improper.cpp
      intercljff.cpp
      intercoulombff.cpp
      interff.cpp
      intergroupff.cpp
      interljff.cpp
      internalcomponent.cpp
      internalff.cpp
      internalparameters.cpp
      internalperturbation.cpp
      intersoftcljff.cpp
      intracljff.cpp
      intracoulombff.cpp
      intraff.cpp
      intragroupff.cpp
      intraljff.cpp
      intrasoftcljff.cpp
      ljpair.cpp
      ljparameter.cpp
      ljparameterdb.cpp
      ljperturbation.cpp
      ljpotential.cpp
      mmdetail.cpp
      multicljcomponent.cpp
      restraint.cpp
      restraintcomponent.cpp
      restraintff.cpp
      selectorangle.cpp
      selectorbond.cpp
      selectordihedral.cpp
      selectorimproper.cpp
      selectormangle.cpp
      selectormbond.cpp
      selectormdihedral.cpp
      selectormimproper.cpp
      softcljcomponent.cpp
      softcljpotential.cpp
      switchingfunction.cpp
      testff.cpp
      threeatomfunctions.cpp
      twoatomfunctions.cpp

      ${SIREMM_HEADERS}
      ${SIREMM_DETAIL_HEADERS}
    )

add_library (SireMM ${SIREMM_SOURCES})

set_target_properties (SireMM PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireMM
                       SireFF
                       SireMol
                       SireCAS
                       SireBase
                       SireMaths
                       SireUnits
                       SireStream
                       ${TBB_LIBRARY}
                       ${TBB_MALLOC_LIBRARY}
                       )

# installation
install( TARGETS SireMM EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREMM_HEADERS} DESTINATION ${SIRE_INCLUDES}/SireMM )
install( FILES ${SIREMM_DETAIL_HEADERS} DESTINATION ${SIRE_INCLUDES}/SireMM/detail )
