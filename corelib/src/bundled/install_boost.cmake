###
### CMake file to install the Boost headers-only libraries
### in the Sire bundle
###

# First, try to find the include/boost directory in the "bundled" directory. If it exists,
# then we don't need to do anything

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
          OUTPUT_QUIET
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
