###
### CMake file to install boost::python
### in the Sire bundle
###

unset(BOOST_PYTHON_LIBRARY CACHE)

if ( ANACONDA_BUILD )
  if (MSVC)
    find_package( Boost 1.31 COMPONENTS
      python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED )
    set ( BOOST_PYTHON_LIBRARY "${Boost_LIBRARIES}" )
    set ( BOOST_PYTHON_HEADERS "${Boost_INCLUDE_DIR}" )
    add_definitions("/DBOOST_PYTHON_NO_LIB")
    #file(READ "${Boost_INCLUDE_DIR}/boost/python/detail/config.hpp" CONFIG_HPP_DATA)
    #if (NOT CONFIG_HPP_DATA MATCHES "BOOST_PYTHON_CONCAT")
    #  configure_file("${Boost_INCLUDE_DIR}/boost/python/detail/config.hpp"
    #    "${Boost_INCLUDE_DIR}/boost/python/detail/config.hpp.orig" COPYONLY)
    #  string(REGEX REPLACE "(#define BOOST_LIB_NAME boost_python##PY_MAJOR_VERSION##PY_MINOR_VERSION)"
    #    "\\1\n#define _BOOST_PYTHON_CONCAT(N, M, m) N ## M ## m\n#define BOOST_PYTHON_CONCAT(N, M, m) _BOOST_PYTHON_CONCAT(N, M, m)\n#define BOOST_LIB_NAME BOOST_PYTHON_CONCAT(boost_python, PY_MAJOR_VERSION, PY_MINOR_VERSION)"
    #    CONFIG_HPP_DATA "${CONFIG_HPP_DATA}")
    #  string(REGEX REPLACE "(#endif  // auto-linking disabled)"
    #    "\\1\n#undef BOOST_PYTHON_CONCAT\n#undef _BOOST_PYTHON_CONCAT"
    #    CONFIG_HPP_DATA "${CONFIG_HPP_DATA}")
    #  file(WRITE "${Boost_INCLUDE_DIR}/boost/python/detail/config.hpp" "${CONFIG_HPP_DATA}")
    #endif()
    message("*** BOOST_PYTHON_LIBRARY ${BOOST_PYTHON_LIBRARY}")
  else()
    find_library( BOOST_PYTHON_LIBRARY
                  NAMES boost_python
                  PATHS ${SIRE_APP}/lib NO_DEFAULT_PATH )
  endif()
elseif ( MSYS )
  message( STATUS "Looking for MSYS version of boost::python..." )
  set (BOOST_ALL_DYN_LINK "YES")
  find_package( Boost 1.31 COMPONENTS python3 REQUIRED )
  set ( BOOST_PYTHON_LIBRARY "${Boost_LIBRARIES}" )
  set ( BOOST_PYTHON_HEADERS "${Boost_INCLUDE_DIR}" )

else()
  find_library( BOOST_PYTHON_LIBRARY 
                NAMES boost_python
                PATHS ${BUNDLE_STAGEDIR}/lib NO_DEFAULT_PATH )
endif()

if ( BOOST_PYTHON_LIBRARY )
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
    set( BOOST_PYTHON_LIBRARY "${BUNDLE_STAGEDIR}/lib/libboost_python.dylib" )
  else()
    set( BOOST_PYTHON_LIBRARY "${BUNDLE_STAGEDIR}/lib/libboost_python.so" )
  endif()

endif()

if ( MSYS )
  if ( BOOST_PYTHON_LIBRARY )
    message( STATUS "Using boost::python in ${BOOST_PYTHON_LIBRARY} | ${BOOST_PYTHON_HEADERS}" )
  else()
    message( FATAL_ERROR "Cannot find boost::python in MSYS!" )
  endif()
else()
  if ( BOOST_PYTHON_LIBRARY )
    set( BOOST_PYTHON_HEADERS "${SIRE_APP}/include" )

    if ( APPLE )
      message( STATUS "Not linking modules to libPython to prevent double-symbols" )
    else()
      set( BOOST_PYTHON_LIBRARY "${PYTHON_LIBRARIES};${BOOST_PYTHON_LIBRARY}" )
    endif()

    message( STATUS "Using bundled boost::python in ${BOOST_PYTHON_LIBRARY} | ${BOOST_PYTHON_HEADERS}" )
    set( SIRE_FOUND_BOOST_PYTHON TRUE )
  else()
    message( STATUS "Strange? Cannot find the installed boost::python library. We cannot compile it, so will need to rely on the system version..." )
  endif()
endif()
