########################################
#
# Set of unit tests to run on a 
# completed Sire installation
#
########################################

set( SIRE_TEST_DIR "${CMAKE_INSTALL_PREFIX}/test" )

# List all of the directories containing unit tests
set (TEST_DIRS "SireBase"
               "SireMM"
               "SireMol"
               "SireQt"
               "SireSystem"
               "SireVol"
               "OpenMM"
               "io"
               "Tools"
     )

# installation
foreach ( TEST_DIR ${TEST_DIRS} )
  install( DIRECTORY ${TEST_DIR} DESTINATION ${SIRE_TEST_DIR}/
           PATTERN .svn EXCLUDE )
endforeach()
