###
### CMake file to install the Boost headers-only libraries
### in the Sire bundle
###

# First, try to find the include/boost directory in the "bundled" directory. If it exists,
# then we don't need to do anything

if ( MSYS )
  message( STATUS "Using MSYS version of Boost" )

  # Sire depends on >= boost 1.31 - all libraries must be dynamically linked
  set (BOOST_ALL_DYN_LINK "YES")
  set (Boost_ADDITIONAL_VERSIONS "1.39" "1.39.0")
  FIND_PACKAGE ( Boost 1.31 REQUIRED )

  if (Boost_FOUND)
    message(STATUS "Boost paths ${Boost_LIBRARY_DIRS} | ${Boost_INCLUDE_DIR}" )
    set ( BOOST_INCLUDE_DIRS "${Boost_INCLUDE_DIR}" )
    include_directories( ${Boost_INCLUDE_DIR} )

    #save the path to this include directory so that it can be
    #used by anything compiling against Sire
    save_sire_variable( "SIRE_BOOST_INCLUDE_DIR" "${Boost_INCLUDE_DIR}" )
  else()
    message(FATAL_ERROR "Cannot find the boost libraries.")
  endif()

else()
  if (EXISTS "${BUNDLE_STAGEDIR}/include/boost" )
    message( STATUS "Have already unpacked a bundled version of Boost")

  else()
    message( STATUS "Unpacking a bundled version of the Boost header libraries..." )

    set( BOOST_ZIPFILE "${CMAKE_SOURCE_DIR}/src/bundled/boost_headers.tar.gz" )
 
    if (EXISTS "${BOOST_ZIPFILE}")
      message( STATUS "Unzipping ${BOOST_ZIPFILE} to ${BUNDLE_STAGEDIR}/include" )
      execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${BOOST_ZIPFILE}
            WORKING_DIRECTORY ${BUNDLE_STAGEDIR}/include
            OUTPUT_QUIET ERROR_QUIET
        )
      endif()
  endif()

  if (EXISTS "${BUNDLE_STAGEDIR}/include/boost")
    message( STATUS "Using Boost headers in ${BUNDLE_STAGEDIR}/include/boost" )
    set( BOOST_INCLUDE_DIRS "${BUNDLE_STAGEDIR}/include" )
    set( SIRE_FOUND_BOOST TRUE )
  else()
    message( STATUS "Strange? Cannot find the Boost include directory.
                     We cannot unpack it, so will need to rely on the system version..." )
  endif()
endif()
