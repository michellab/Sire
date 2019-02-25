########################################
#
# CMake file for library:SireCluster
#
########################################

# Export symbols when building the library
add_definitions( "-DSIRECLUSTER_BUILD" )

# Third Party dependencies of this module

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

if ( SIRE_USE_MPI )
  message( STATUS " ** Adding MPI support to SireCluster ** " )
  message( STATUS "MPI headers: ${MPI_INCLUDE_PATH}" )
  message( STATUS "MPI libraries: ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY}" )
  message( STATUS "MPI definitions: ${MPI_DEFINITIONS}" )

  include_directories( BEFORE SYSTEM ${MPI_INCLUDE_PATH} )

  add_definitions( -DSIRE_USE_MPI ${MPI_DEFINITIONS} )

  set ( SIRECLUSTER_EXTRA_LIBRARIES ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY} )

  save_sire_variable( "SIRE_MPI_INCLUDE_PATH" "${MPI_INCLUDE_PATH}" )
  save_sire_variable( "SIRE_MPI_DEFINITIONS" "-DSIRE_USE_MPI ${MPI_DEFINITIONS}" )

  set ( SIRECLUSTER_MPI_HEADERS

        mpi/messages.h
        mpi/mpicluster.h
        mpi/mpifrontend.h
        mpi/p2pcomm.h
        mpi/receivequeue.h
        mpi/reply.h
        mpi/reservationmanager.h
        mpi/sendqueue.h

      )

  set ( SIRECLUSTER_MPI_SOURCES

        mpi/messages.cpp
        mpi/mpicluster.cpp
        mpi/mpifrontend.cpp
        mpi/p2pcomm.cpp
        mpi/receivequeue.cpp
        mpi/reply.cpp
        mpi/reservationmanager.cpp
        mpi/sendqueue.cpp

      )

else()
  message( STATUS "SireCluster will be compiled without MPI support." )   
  save_sire_variable( "SIRE_MPI_DEFINITIONS" "" )
endif()

save_sire_variable( "SIRE_USE_MPI" "${MPI_FOUND}" )

# Define the headers in SireCluster
set ( SIRECLUSTER_HEADERS

      backend.h
      cluster.h
      frontend.h
      node.h
      nodes.h
      promise.h      
      workpacket.h
    )

# Define the sources in SireCluster
set ( SIRECLUSTER_SOURCES

      register_sirecluster.cpp

      backend.cpp
      cluster.cpp
      frontend.cpp
      node.cpp
      nodes.cpp
      promise.cpp
      workpacket.cpp
 
      ${SIRECLUSTER_MPI_SOURCES}
      ${SIRECLUSTER_MPI_HEADERS}

      ${SIRECLUSTER_HEADERS}
    )

add_library (SireCluster ${SIRECLUSTER_SOURCES})

set_target_properties (SireCluster PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireCluster
                       SireMaths
                       SireBase
                       SireError
                       SireStream
                       ${SIRECLUSTER_EXTRA_LIBRARIES}
                       )

save_sire_variable( "SIRECLUSTER_EXTRA_LIBRARIES" "${SIRECLUSTER_EXTRA_LIBRARIES}" )

# installation
install( TARGETS SireCluster EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIRECLUSTER_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireCluster )

