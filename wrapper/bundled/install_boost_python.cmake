###
### CMake file to install boost::python
### in the Sire bundle
###

# I am not using this now as it breaks the finding and linking of boost python 
# in the windows Anaconda build, and we are now moving fully away from the bundled
# dependencies. Please look in wrapper/CMakeLists.txt to see how we are now finding
# and linking to boost python. The main difference is forcing the link to the 
# dynamic library, and not including the fix for unwind_type.hpp below. If this
# fix is still needed then please port it into wrapper/CMakeLists.txt.
#
# I hope this is ok. If there are any problems then please get in 
# touch with me (Christopher Woods) and we can have a chat :-)

unset(Boost_LIBRARIES CACHE)

if ( ANACONDA_BUILD )
  if (MSVC)
    find_package( Boost 1.31 COMPONENTS
      python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED )
    set ( Boost_LIBRARIES "${Boost_LIBRARIES}" )
    set ( BOOST_PYTHON_HEADERS "${Boost_INCLUDE_DIR}" )
    add_definitions("/DBOOST_PYTHON_NO_LIB")
    # unwind_type fails with MSVC >= 15.8.0
    # (https://github.com/boostorg/python/issues/228)
    file(READ "${Boost_INCLUDE_DIR}/boost/python/detail/unwind_type.hpp" UNWIND_TYPE_HPP_DATA)
    if (UNWIND_TYPE_HPP_DATA MATCHES "#ifndef _MSC_VER")
      configure_file("${Boost_INCLUDE_DIR}/boost/python/detail/unwind_type.hpp"
        "${Boost_INCLUDE_DIR}/boost/python/detail/unwind_type.hpp.orig" COPYONLY)
      string(REGEX REPLACE "#ifndef _MSC_VER"
        "#if _MSC_VER >= 1915"
        UNWIND_TYPE_HPP_DATA "${UNWIND_TYPE_HPP_DATA}")
      file(WRITE "${Boost_INCLUDE_DIR}/boost/python/detail/unwind_type.hpp" "${UNWIND_TYPE_HPP_DATA}")
    endif()
  else()
    find_library( Boost_LIBRARIES
                  NAMES boost_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}
                  PATHS ${SIRE_APP}/lib NO_DEFAULT_PATH )
  endif()
elseif ( MSYS )
  message( STATUS "Looking for MSYS version of boost::python..." )
  set (BOOST_ALL_DYN_LINK "YES")
  find_package( Boost 1.31 COMPONENTS python3 REQUIRED )
  set ( Boost_LIBRARIES "${Boost_LIBRARIES}" )
  set ( BOOST_PYTHON_HEADERS "${Boost_INCLUDE_DIR}" )

else()
  find_library( Boost_LIBRARIES
                NAMES boost_python
                PATHS ${BUNDLE_STAGEDIR}/lib NO_DEFAULT_PATH )
endif()

if ( Boost_LIBRARIES )
  message( STATUS "Have already compiled a bundled version of boost::python")
else()
  message( STATUS "Compiling and installing a bundled version of boost::python" )

  set( BOOST_PYTHON_ZIPFILE "${CMAKE_SOURCE_DIR}/bundled/boost_python.tar.gz" )

  if (EXISTS "${BOOST_PYTHON_ZIPFILE}")
    set( BOOST_PYTHON_BUILD_DIR "${BUNDLE_BUILDDIR}/boost_python/build" )

    if (NOT EXISTS "${BOOST_PYTHON_BUILD_DIR}")
      message( STATUS "Unzipping ${BOOST_PYTHON_ZIPFILE} to ${BOOST_PYTHON_BUILD_DIR}" )
      execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf ${BOOST_PYTHON_ZIPFILE}
          WORKING_DIRECTORY ${BUNDLE_BUILDDIR}
          OUTPUT_QUIET ERROR_QUIET
      )
    endif()
  endif()

  message( STATUS "Patience... Configuring boost::python" )
  if (MSYS)
    list( APPEND COMPILE_OPTIONS "-G" )
    list( APPEND COMPILE_OPTIONS "MSYS Makefiles" )
  endif()

  list( APPEND COMPILE_OPTIONS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" )
  list( APPEND COMPILE_OPTIONS "-DPYTHON_LIBRARIES=${PYTHON_LIBRARIES}" )
  list( APPEND COMPILE_OPTIONS "-DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}" )
  list( APPEND COMPILE_OPTIONS "-DCMAKE_INSTALL_PREFIX=${BUNDLE_STAGEDIR}" )
  list( APPEND COMPILE_OPTIONS "-DCMAKE_CXX_FLAGS=${SIRE_CXX_NG_FLAGS} -O2" )
  list( APPEND COMPILE_OPTIONS "../src" )

  message( STATUS "${CMAKE_COMMAND} ${COMPILE_OPTIONS}" )

  execute_process(
          COMMAND ${CMAKE_COMMAND} ${COMPILE_OPTIONS}
          WORKING_DIRECTORY ${BOOST_PYTHON_BUILD_DIR}
          OUTPUT_QUIET ERROR_QUIET
  )

  message( STATUS "Patience... Building boost::python" )
  execute_process(
          COMMAND ${CMAKE_MAKE_PROGRAM} -k -j ${NCORES}
          WORKING_DIRECTORY ${BOOST_PYTHON_BUILD_DIR}
          OUTPUT_QUIET ERROR_QUIET
  )

  message( STATUS "Patience... Installing boost::python" )
  execute_process(
          COMMAND ${CMAKE_MAKE_PROGRAM} install
          WORKING_DIRECTORY ${BOOST_PYTHON_BUILD_DIR}
          OUTPUT_QUIET ERROR_QUIET
  )

  if (APPLE)
    set( Boost_LIBRARIES "${BUNDLE_STAGEDIR}/lib/libboost_python.dylib" )
  else()
    set( Boost_LIBRARIES "${BUNDLE_STAGEDIR}/lib/libboost_python.so" )
  endif()

endif()

if ( MSYS )
  if ( Boost_LIBRARIES )
    message( STATUS "Using boost::python in ${Boost_LIBRARIES} | ${BOOST_PYTHON_HEADERS}" )
  else()
    message( FATAL_ERROR "Cannot find boost::python in MSYS!" )
  endif()
else()
  if ( Boost_LIBRARIES )
    set( BOOST_PYTHON_HEADERS "${SIRE_APP}/include" )

    if ( APPLE )
      message( STATUS "Not linking modules to libPython to prevent double-symbols" )
    else()
      set( Boost_LIBRARIES "${PYTHON_LIBRARIES};${Boost_LIBRARIES}" )
    endif()

    message( STATUS "Using bundled boost::python in ${Boost_LIBRARIES} | ${BOOST_PYTHON_HEADERS}" )
    set( SIRE_FOUND_BOOST_PYTHON TRUE )
  else()
    message( STATUS "Strange? Cannot find the installed boost::python library. We cannot compile it, so will need to rely on the system version..." )
  endif()
endif()
