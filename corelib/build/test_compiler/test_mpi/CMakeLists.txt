########################################
#
# CMake file for test_mpi
#
########################################

if ( SIRE_USE_MPI )

    include_directories( BEFORE SYSTEM ${MPI_INCLUDE_PATH} )
    add_definitions( ${MPI_DEFINITIONS} )

    add_executable( test_mpi test_mpi.cpp )

    target_link_libraries( test_mpi ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY} )

endif()
