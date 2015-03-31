
# Find the svn program
find_program( SVN_EXECUTABLE
              NAMES svn
              PATHS /usr/bin  /usr/local/bin ~/local/bin )

IF (SVN_EXECUTABLE)
  
   #get info about the repository
   exec_program( ${SVN_EXECUTABLE}
                 ARGS "info ${CMAKE_SOURCE_DIR}"
                 OUTPUT_VARIABLE SVN_XML_VERSION_DATA )

   #convert all returns to spaces
   string( REPLACE "\n" " " SVN_XML_VERSION_DATA "${SVN_XML_VERSION_DATA}" )

   #extract the repository URL  (line is "URL: http://repository.url/sire/branches...")
   string( REGEX REPLACE ".*URL: ([^ ]+).*" "\\1" SVN_REPOSITORY_URL "${SVN_XML_VERSION_DATA}" )

endif()

# Find the svnversion program
find_program( SVNVERSION_EXECUTABLE
              NAMES svnversion
              PATHS /usr/bin  /usr/local/bin ~/local/bin )

if (SVNVERSION_EXECUTABLE)

    #get the version number string for this copy
    exec_program( ${SVNVERSION_EXECUTABLE}
                  ARGS "${CMAKE_SOURCE_DIR}"
                  OUTPUT_VARIABLE SVN_VERSION_NUMBER )

endif()

# now write this version info to a header file
configure_file( ${CMAKE_SOURCE_DIR}/build/cmake/sire_version_template.h
                ${CMAKE_BINARY_DIR}/sire_version.h
                ESCAPE_QUOTES )

# install this file into the Sire include directory
install (FILES ${CMAKE_BINARY_DIR}/sire_version.h 
         DESTINATION ${SIRE_INCLUDES})
