########################################
#
# CMake file for library:SireBase
#
########################################

# Third Party dependencies of this module

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

if (SIRE_FOUND_CPUID)
  add_definitions( "-DSIRE_FOUND_CPUID" )
  add_definitions( "${CPUID_DEFINITIONS}" )
  include_directories( ${CPUID_INCLUDE_DIR} )
endif()

# Define the header files for SireBase
set ( SIREBASE_HEADERS

      array2d.hpp
      array2d.h
      arrayproperty.hpp
      booleanproperty.h
      chunkedhash.hpp
      chunkedvector.hpp
      combineproperties.h
      countflops.h
      cpuid.h
      errors.h
      findexe.h
      getinstalldir.h
      globalsharedpointer.hpp
      increment.h
      incremint.h
      linktoproperty.h
      lengthproperty.h
      majorminorversion.h
      meminfo.h
      nullproperty.hpp
      numberproperty.h
      packedarray2d.h
      packedarray2d.hpp
      packedarrays.h
      pairmatrix.hpp
      process.h
      properties.h
      property.h
      propertylist.h
      propertymap.h
      quickcopy.hpp
      qvariant_metatype.h
      savestate.h
      shareddatapointer.hpp
      sharedpolypointer.hpp
      sharedpolypointer_cast.hpp
      sparsematrix.hpp
      stringmangler.h
      stringproperty.h
      tempdir.h
      trigarray2d.h
      trigarray2d.hpp

    )

# Define the sources in SireBase
set ( SIREBASE_SOURCES

      register_sirebase.cpp

      array2d.cpp
      arrayproperty.cpp
      booleanproperty.cpp
      chunkedhash.cpp
      chunkedvector.cpp
      combineproperties.cpp
      countflops.cpp
      cpuid.cpp
      errors.cpp
      findexe.cpp
      getinstalldir.cpp
      globalsharedpointer.cpp
      lengthproperty.cpp
      linktoproperty.cpp
      majorminorversion.cpp
      meminfo.cpp
      numberproperty.cpp
      packedarray2d.cpp
      process.cpp
      properties.cpp
      property.cpp
      propertylist.cpp
      propertymap.cpp
      savestate.cpp
      sharedpolypointer.cpp
      stringmangler.cpp 
      stringproperty.cpp    
      tempdir.cpp
      trigarray2d.cpp 

      ${SIREBASE_HEADERS}
    )

add_library (SireBase ${SIREBASE_SOURCES})

set_target_properties (SireBase PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireBase
                       SireUnits
                       SireID
                       SireStream
                       )

if (SIRE_FOUND_CPUID)
  target_link_libraries( SireBase ${CPUID_LIBRARY} )
endif()

# installation
install( TARGETS SireBase EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREBASE_HEADERS} 
         DESTINATION ${SIRE_INCLUDES}/SireBase )
