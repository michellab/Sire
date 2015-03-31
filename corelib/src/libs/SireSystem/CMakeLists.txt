########################################
#
# CMake file for library:SireSystem
#
########################################

# Third Party dependencies of this module

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Define the headers in SireSystem
set ( SIRESYSTEM_HEADERS
      anglecomponent.h
      chargeconstraint.h
      closemols.h
      constraint.h
      constraints.h
      checkpoint.h
      delta.h
      dihedralcomponent.h
      distancecomponent.h
      energymonitor.h
      errors.h
      freeenergymonitor.h
      geometrycomponent.h
      idassigner.h
      identityconstraint.h
      moleculeconstraint.h
      monitorcomponent.h
      monitorcomponents.h
      monitorid.h
      monitoridentifier.h
      monitoridx.h
      monitormonitor.h
      monitorname.h
      monitorproperty.h
      perturbationconstraint.h
      polarisecharges.h
      spacewrapper.h
      sysid.h
      sysidentifier.h
      sysidx.h
      sysname.h
      system.h
      systemmonitor.h
      systemmonitors.h
      volmapmonitor.h
    )

# Define the sources in SireSystem
set ( SIRESYSTEM_SOURCES

      register_siresystem.cpp

      anglecomponent.cpp
      chargeconstraint.cpp
      closemols.cpp
      constraint.cpp
      constraints.cpp
      checkpoint.cpp
      delta.cpp
      dihedralcomponent.cpp
      distancecomponent.cpp
      energymonitor.cpp
      errors.cpp
      freeenergymonitor.cpp
      geometrycomponent.cpp
      idassigner.cpp
      identityconstraint.cpp
      moleculeconstraint.cpp
      monitorcomponent.cpp
      monitorcomponents.cpp
      monitorid.cpp
      monitoridentifier.cpp
      monitormonitor.cpp
      monitorproperty.cpp
      perturbationconstraint.cpp
      polarisecharges.cpp
      spacewrapper.cpp
      sysid.cpp
      sysidentifier.cpp
      system.cpp
      systemmonitor.cpp
      systemmonitors.cpp
      volmapmonitor.cpp
    
      ${SIRESYSTEM_HEADERS}
    )

add_library (SireSystem ${SIRESYSTEM_SOURCES})

set_target_properties (SireSystem PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${SIRE_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireSystem
                       SireMM
                       SireFF
                       SireMol
                       SireCAS
                       SireBase
                       SireStream
                       )

# installation
install( TARGETS SireSystem EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIRESYSTEM_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireSystem )

