########################################
#
# CMake file for library:SireStream
#
########################################

# Export symbols when building the library
add_definitions( "-DSIRESTREAM_BUILD" )

include_directories(${CMAKE_SOURCE_DIR}/src/libs)

set ( SIRESTREAM_HEADERS
      datastream.h
      errors.h
      magic_error.h
      md5sum.h
      metatype.h
      registeralternativename.h
      shareddatastream.h
      sharestrings.h
      streamdata.hpp
      streampolypointer.hpp
      version_error.h
      versionid.h
    )

set ( SIRESTREAM_SOURCES

      datastream.cpp
      errors.cpp
      getmagic.cpp
      magic_error.cpp
      md5sum.cpp
      metatype.cpp
      shareddatastream.cpp
      sharestrings.cpp
      streamdata.cpp
      version_error.cpp

      ThirdParty/md5.c
    
      ${SIRESTREAM_HEADERS}
    )


add_library ( SireStream ${SIRESTREAM_SOURCES} )

set_target_properties (SireStream PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

target_link_libraries (SireStream
                       SireError
                      )

# installation
install( TARGETS SireStream  EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIRESTREAM_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireStream )

