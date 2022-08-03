########################################
#
# CMake file for library:SireBase
#
########################################

# Export symbols when building the library
add_definitions( "-DSIREBASE_BUILD" )

# Third Party dependencies of this module

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# This library uses Intel Threaded Building blocks
include_directories(${TBB_INCLUDE_DIR})

# Define the header files for SireBase
set ( SIREBASE_HEADERS

      array2d.hpp
      array2d.h
      arrayproperty.hpp
      booleanproperty.h
      chunkedhash.hpp
      chunkedvector.hpp
      combineproperties.h
      convert_property.hpp
      countflops.h
      cpuid.h
      errors.h
      findexe.h
      generalunitproperty.h
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
      range.h
      ranges.h
      refcountdata.h
      sire_process.h
      parallel.h
      podproperty.hpp
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
      slice.h
      sparsematrix.hpp
      stringmangler.h
      stringproperty.h
      tempdir.h
      timeproperty.h
      trigarray2d.h
      trigarray2d.hpp
      unittest.h
      variantproperty.h
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
      generalunitproperty.cpp
      getinstalldir.cpp
      globalsharedpointer.cpp
      lengthproperty.cpp
      linktoproperty.cpp
      majorminorversion.cpp
      meminfo.cpp
      numberproperty.cpp
      packedarray2d.cpp
      parallel.cpp
      process.cpp
      properties.cpp
      property.cpp
      propertylist.cpp
      propertymap.cpp
      range.cpp
      ranges.cpp
      refcountdata.cpp
      savestate.cpp
      sharedpolypointer.cpp
      slice.cpp
      stringmangler.cpp
      stringproperty.cpp
      tempdir.cpp
      timeproperty.cpp
      trigarray2d.cpp
      unittest.cpp
      variantproperty.cpp

      test_sharedptr.cpp

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
                       ${TBB_LIBRARY}
                       ${TBB_MALLOC_LIBRARY}
                       )

# installation
install( TARGETS SireBase EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREBASE_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireBase )
