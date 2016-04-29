########################################
#
# CMake file for application:testmulti
#
########################################

# Third Party dependencies of this module
include_directories(${BOOST_INCLUDE_DIR})

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

set ( TESTMULTI_SOURCES
      main.cpp
    )

set ( SPEEDTEST_SOURCES
      speedtest.cpp
    )

set ( SQRTTEST_SOURCES
      sqrttest.cpp
    )

set ( VECTEST_SOURCES
      vectest.cpp
    )

# generate the application
add_executable ( sire_testmulti ${TESTMULTI_SOURCES} )
add_executable ( sire_multispeedtest ${SPEEDTEST_SOURCES} )
add_executable ( sire_sqrttest ${SQRTTEST_SOURCES} )
add_executable ( sire_vectest ${VECTEST_SOURCES} )

# Link to other Sire libraries
target_link_libraries (sire_testmulti
                       SireMaths
                      )

# installation
install( TARGETS sire_testmulti
         RUNTIME DESTINATION ${SIRE_BIN}
       )

# Link to other Sire libraries
target_link_libraries (sire_multispeedtest
                       SireMaths
                      )

# installation
install( TARGETS sire_multispeedtest
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

# Link to other Sire libraries
target_link_libraries (sire_sqrttest
                       SireMaths
                      )

# installation
install( TARGETS sire_sqrttest
         RUNTIME DESTINATION ${SIRE_BIN}
       )


# Link to other Sire libraries
target_link_libraries (sire_vectest
                       SireMaths
                      )

# installation
install( TARGETS sire_vectest
         RUNTIME DESTINATION ${SIRE_BIN}
       )

if (APPLE)
  add_custom_command(TARGET sire_testmulti
      POST_BUILD 
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_LIBS}" sire_testmulti 
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_BUNDLED_LIBS}" sire_testmulti)
  add_custom_command(TARGET sire_multispeedtest
      POST_BUILD       
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_LIBS}" sire_multispeedtest
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_BUNDLED_LIBS}" sire_multispeedtest)
  add_custom_command(TARGET sire_sqrttest
      POST_BUILD
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_LIBS}" sire_sqrttest
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_BUNDLED_LIBS}" sire_sqrttest )
  add_custom_command(TARGET sire_vectest
      POST_BUILD
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_LIBS}" sire_vectest
      COMMAND ${CMAKE_INSTALL_NAME_TOOL} -add_rpath "@executable_path/../${SIRE_BUNDLED_LIBS}" sire_vectest )
endif() 
