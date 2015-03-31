########################################
#
# CMake file for application:inquire
#
########################################

# Third Party dependencies of this module
include_directories(${BOOST_INCLUDE_DIR})

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

set( CMAKE_AUTOMOC TRUE )

set ( INQUIRE_HEADERS

      qcustomplot.h
    )

set ( INQUIRE_SOURCES

      main.cpp

      qcustomplot.cpp

      ${INQUIRE_HEADERS}
    )

# generate the application
add_executable ( inquire ${INQUIRE_SOURCES} )

# Link to other Sire libraries
target_link_libraries (inquire
                       SireAnalysis
                       Squire
                       SireMM
                       SireSystem
                       SireIO
                       SireStream
                       Qt5::Gui
                       Qt5::Widgets
                       Qt5::PrintSupport
                      )

# installation
install( TARGETS inquire
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

if (APPLE)
  add_custom_command(TARGET inquire
      POST_BUILD
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_LIBS}" inquire
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_BUNDLED_LIBS}" inquire)
endif()
