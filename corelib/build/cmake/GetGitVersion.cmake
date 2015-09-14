
# Find the git program
find_program( GIT_EXECUTABLE
              NAMES git
              PATHS /usr/bin  /usr/local/bin ~/local/bin )

IF (GIT_EXECUTABLE)
   #get info about the repository
   execute_process( COMMAND ${GIT_EXECUTABLE} config --get remote.origin.url
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_REMOTE_URL
                    OUTPUT_STRIP_TRAILING_WHITESPACE )

   #get the current commit branch
   execute_process( COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_BRANCH
                    OUTPUT_STRIP_TRAILING_WHITESPACE )

   #get the current hash, author name and commit date and time
   execute_process( COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%H/%cn/%aI
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_COMMIT_INFO
                    OUTPUT_STRIP_TRAILING_WHITESPACE )

   #find out whether the current checkout is clean (has nothing that
   #hasn't yet been committed)
   execute_process( COMMAND ${GIT_EXECUTABLE} status --porcelain
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_UNCLEAN
                    OUTPUT_STRIP_TRAILING_WHITESPACE )

   set( SVN_REPOSITORY_URL "${GIT_REMOTE_URL} | ${GIT_BRANCH}" )
   set( SVN_VERSION_NUMBER "${GIT_COMMIT_INFO}" )

   message( STATUS "REPOSITORY_URL = ${SVN_REPOSITORY_URL}" )
   message( STATUS "VERSION_NUMBER = ${SVN_VERSION_NUMBER}" )

   if ( GIT_UNCLEAN )
     set( SVN_VERSION_NUMBER "${SVN_VERSION_NUMBER} | **UNCLEAN**" )
     message( STATUS "WARNING: Directory contains uncommitted changes!" )
   endif()

else()
   set( SVN_REPOSITORY_URL "UNKNOWN REPOSITORY" )
   set( SVN_VERSION_NUMBER "UNKNOWN VERSION | **UNCLEAN**" )
endif()

# now write this version info to a header file
configure_file( ${CMAKE_SOURCE_DIR}/build/cmake/sire_version_template.h
                ${CMAKE_BINARY_DIR}/sire_version.h
                ESCAPE_QUOTES )

# install this file into the Sire include directory
install (FILES ${CMAKE_BINARY_DIR}/sire_version.h 
         DESTINATION ${SIRE_INCLUDES})
