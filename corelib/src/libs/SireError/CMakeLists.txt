########################################
#
# CMake file for library:SireError
#
########################################

include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Is backtrace available on this system?
check_symbol_exists( backtrace "execinfo.h" HAVE_BACKTRACE )

if (HAVE_BACKTRACE)
   message( STATUS "Backtrace functionality is available." )
   add_definitions( -D_HAVE_EXECINFO_H_ )
   add_definitions( -DSIRE_ENABLE_BACKTRACE )
endif()

set ( SIREERROR_HEADERS
      errors.h
      exception.h
      getbacktrace.h
      printerror.h
      version_error.h
    )

set ( SIREERROR_SOURCES

      errors.cpp
      exception.cpp
      getbacktrace.cpp
      printerror.cpp
    
      ${SIREERROR_HEADERS}
    )

add_library ( SireError ${SIREERROR_SOURCES} )

set_target_properties (SireError PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${SIRE_VERSION_MAJOR}
                      )

target_link_libraries (SireError
                       Qt5::Core
                      )

# installation
install( TARGETS SireError EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREERROR_HEADERS} 
         DESTINATION ${SIRE_INCLUDES}/SireError )

