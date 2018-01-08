########################################
#
# CMake file for application:sire
#
########################################

# Third Party dependencies of this module
include_directories(${BOOST_INCLUDE_DIR})

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# If we have NetCDF, then include the header files
if ( ${SIRE_FOUND_NETCDF} )
  message( STATUS "Adding in NetCDF support to SireIO using ${NetCDF_ROOT_DIR}" )
  message( STATUS "${NetCDF_INCLUDE_DIR} | ${NetCDF_LIBRARIES}" )
  set ( SIREAPP_NETCDF_LIBRARIES ${NetCDF_LIBRARIES} )
endif()

if ( SIRE_USE_MPI )
  message( STATUS " ** Adding MPI support to sire ** " )
  message( STATUS "MPI headers: ${MPI_INCLUDE_PATH}" )
  message( STATUS "MPI libraries: ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY}" )
  message( STATUS "MPI definitions: ${MPI_DEFINITIONS}" )

  include_directories( BEFORE SYSTEM ${MPI_INCLUDE_PATH} )

  add_definitions( -DSIRE_USE_MPI ${MPI_DEFINITIONS} )
endif()

set ( SIRE_RCCS

      sire.qrc
    )

if ( CMAKE_GENERATOR MATCHES "Xcode" )
  message( STATUS "Xcode project is not compilable..." )
else()
  qt5_add_resources( SIRE_RCC_SRCS ${SIRE_RCCS} )
endif()

set ( SIRE_SOURCES

      main.cpp

      ${SIRE_HEADERS}
    )

# generate the application
add_executable ( sire ${SIRE_SOURCES} ${SIRE_RCC_SRCS} )

# Link to other Sire libraries
target_link_libraries (sire
                       SireIO
                       Squire
                       SireMove
                       SireMol
                       SireCluster
                       SireBase
                       SireStream
                       ${SIREAPP_NETCDF_LIBRARIES}
                      )

# installation
install( TARGETS sire
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

if (APPLE)
  add_custom_command(TARGET sire
      POST_BUILD
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_LIBS}" sire
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_BUNDLED_LIBS}" sire)
endif()
