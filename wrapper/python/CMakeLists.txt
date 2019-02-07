########################################
#
# CMake file for sire_python
# (mpi-aware python executable)
#
########################################

# Third Party dependencies of this module
include_directories( ${BOOST_INCLUDE_DIR}
                     ${PYTHON_INCLUDE_DIR} )

# Sire include paths
include_directories( BEFORE ${SIRE_INCLUDE_DIR} )

# Other python wrapping directories
include_directories(${CMAKE_SOURCE_DIR})

if (SIRE_USE_MPI)
  add_definitions( ${SIRE_MPI_DEFINITIONS} )
  include_directories( BEFORE SYSTEM ${SIRE_MPI_INCLUDE_PATH} )
endif()

# Define the sources in Helpers
set ( SIREPYTHON_SOURCES

      main.cpp

    )

# Create the library that holds all of the class wrappers
add_executable (sire_python ${SIREPYTHON_SOURCES})

target_link_libraries (sire_python
                       SirePython
                       SIRE_SireCluster
                       SIRE_Squire
                       SIRE_SireMol
                       SIRE_SireMove
                       SIRE_SireSystem
                       SIRE_SireIO
                       SIRE_SireAnalysis
                       SIRE_SireError
                       ${SIRECLUSTER_EXTRA_LIBRARIES}
                      )

# installation
install( TARGETS sire_python
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

if (APPLE)
  add_custom_command(TARGET sire_python
      POST_BUILD
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_LIBS}" sire_python
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_BUNDLED_LIBS}" sire_python )
endif()
         
install (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts
         DESTINATION ${SIRE_SHARE}
         PATTERN ".svn" EXCLUDE
        )

# get the list of all scripts
file(GLOB EXTRA_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/scripts/*.py")
  
foreach ( EXTRA_SCRIPT ${EXTRA_SCRIPTS} )
  # get the filename and strip off the .py extension
  get_filename_component(PYSCRIPT ${EXTRA_SCRIPT} NAME_WE)
  install( CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\"
                                      -E create_symlink sire_python ${PYSCRIPT}
                                 WORKING_DIRECTORY \"${CMAKE_INSTALL_PREFIX}/${SIRE_BIN}\")" )
endforeach()
