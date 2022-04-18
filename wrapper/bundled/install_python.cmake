###
### CMake file to install the Python 3 interpreter
### in the Sire bundle
###

if (ANACONDA_BUILD)
  if (NOT DEFINED PYTHON_EXECUTABLE)
    # we will just use the python that comes with anaconda
    if (MSYS)
        set (PYTHON_EXECUTABLE "${ANACONDA_BASE}/python" )
    elseif (MSVC)
        set (PYTHON_EXECUTABLE "${ANACONDA_BASE}/python.exe" )
    else()
        set (PYTHON_EXECUTABLE "${ANACONDA_BASE}/bin/python3" )
    endif()
  endif()

  find_package( PythonInterp REQUIRED )

  set( PYTHON_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}" )

  unset(PYTHON_LIBRARY CACHE)

  # sys.abiflags is an empty string for Python >= 3.8
  if (${PYTHON_VERSION} LESS 3.8)
    set(PYTHON_ABIFLAGS "m")
  else()
    set(PYTHON_ABIFLAGS "")
  endif()

  # Find the python library that comes with anaconda
  if (MSYS)
    message( STATUS "Find python${PYTHON_VERSION} in ${ANACONDA_BASE}")
    set( PYTHON_LIBRARY "${ANACONDA_BASE}/lib/python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.lib" )
    set( PYTHON_INCLUDE_DIR "${ANACONDA_BASE}/include" )
    find_package(PythonLibs ${PYTHON_VERSION} REQUIRED)
  elseif (MSVC)
    find_library( PYTHON_LIBRARY
                  NAMES python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.lib
                  PATHS ${ANACONDA_BASE}/libs NO_DEFAULT_PATH )
  else()
    find_library( PYTHON_LIBRARY
                  NAMES python${PYTHON_VERSION}${PYTHON_ABIFLAGS}
                  PATHS ${ANACONDA_BASE}/lib NO_DEFAULT_PATH )
  endif()

  if (NOT PYTHON_LIBRARY)
    message( FATAL_ERROR "Where is the python library that comes with anaconda? "
                         "It cannot be found. Please check that your anaconda "
                         "installation is complete." )
  endif()

elseif (MSYS)
  # use the python that comes with msys
  unset(PYTHON_LIBRARY CACHE)


  if ( PYTHON_EXECUTABLE )
    message( STATUS "Found Python executable ${PYTHON_EXECUTABLE}" )
  else()
    message( FATAL_ERROR "Where is the python executable?" )
  endif()

  find_package( PythonLibs 3.3  )

  if ( PYTHON_LIBRARIES )
    message( STATUS "Python paths ${PYTHON_LIBRARIES} | ${PYTHON_INCLUDE_DIR} | ${PYTHON_VERSION}" )
  else()
    message( FATAL_ERROR "Cannot find the msys installation of Python!" )
  endif()

  find_package( PythonInterp 3.3 )

  set( SIRE_FOUND_PYTHON TRUE )

else()
  # Need to set the version of Python bundled with Sire
  set( PYTHON_VERSION "3.3" )
  set( PYTHON_ABIFLAGS "m" )

  # First, try to find the Python library in the "bundled" directory. If it exists,
  # then we don't need to do anything
  unset(PYTHON_LIBRARY CACHE)

  find_library( PYTHON_LIBRARY
                NAMES python${PYTHON_VERSION}${PYTHON_ABIFLAGS}
                PATHS ${BUNDLE_STAGEDIR}/lib NO_DEFAULT_PATH )

  if ( PYTHON_LIBRARY )
    message( STATUS "Have already compiled a bundled version of Python 3")
  else()
    message( STATUS "Compiling and installing a bundled version of Python 3" )

    set( PYTHON_ZIPFILE "${CMAKE_SOURCE_DIR}/bundled/python3.tar.gz" )

    if (EXISTS "${PYTHON_ZIPFILE}")
      set( PYTHON_BUILD_DIR "${BUNDLE_BUILDDIR}/python3" )

      if (NOT EXISTS "${PYTHON_BUILD_DIR}")
        message( STATUS "Unzipping ${PYTHON_ZIPFILE} to ${PYTHON_BUILD_DIR}" )
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${PYTHON_ZIPFILE}
            WORKING_DIRECTORY ${BUNDLE_BUILDDIR}
        )
      endif()

      list( APPEND PYTHON_OPTIONS "--enable-shared" )
      list( APPEND PYTHON_OPTIONS "--prefix=${BUNDLE_STAGEDIR}" )
      list( APPEND PYTHON_OPTIONS "CC=${CMAKE_C_COMPILER}" )
      list( APPEND PYTHON_OPTIONS "CXX=${CMAKE_CXX_COMPILER}" )

      if (NOT APPLE)
        list( APPEND PYTHON_OPTIONS "LDFLAGS=-Wl,-rpath='$$ORIGIN/../lib'" )
      endif()

      message( STATUS "${PYTHON_OPTIONS}" )

      message( STATUS "Patience... Configuring Python..." )
      execute_process( COMMAND ${PYTHON_BUILD_DIR}/configure ${PYTHON_OPTIONS}
                       WORKING_DIRECTORY ${PYTHON_BUILD_DIR} )

      message( STATUS "Patience... Compiling Python..." )
      execute_process( COMMAND ${CMAKE_MAKE_PROGRAM} -k -j ${NCORES}
                       WORKING_DIRECTORY ${PYTHON_BUILD_DIR} )

      message( STATUS "Patience... Installing Python..." )
      execute_process( COMMAND ${CMAKE_MAKE_PROGRAM} install
                     WORKING_DIRECTORY ${PYTHON_BUILD_DIR} )

      if (APPLE)
        set( PYTHON_LIBRARY "${BUNDLE_STAGEDIR}/lib/libpython${PYTHON_VERSION}${PYTHON_ABIFLAGS}.dylib" )
        execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libpython${PYTHON_VERSION}${PYTHON_ABIFLAGS}.dylib" ${PYTHON_LIBRARY} )
      else()
        set( PYTHON_LIBRARY "${BUNDLE_STAGEDIR}/lib/libpython${PYTHON_VERSION}${PYTHON_ABIFLAGS}.so" )
      endif()
    endif()

    # Now lets unpack and install setuptools and pip, so that we can easily
    # install other python modules
    set( SETUPTOOLS_ZIPFILE "${CMAKE_SOURCE_DIR}/bundled/setuptools.tar.gz" )
    set( PIP_ZIPFILE "${CMAKE_SOURCE_DIR}/bundled/pip.tar.gz" )

    if (EXISTS "${SETUPTOOLS_ZIPFILE}")
      set( SETUPTOOLS_BUILD_DIR "${BUNDLE_BUILDDIR}/setuptools" )

      if (NOT EXISTS "${SETUPTOOLS_BUILD_DIR}")
        message( STATUS "Unzipping ${SETUPTOOLS_ZIPFILE} to ${SETUPTOOLS_BUILD_DIR}" )
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${SETUPTOOLS_ZIPFILE}
            WORKING_DIRECTORY ${BUNDLE_BUILDDIR}
        )

        execute_process(
            COMMAND ${BUNDLE_STAGEDIR}/bin/python3 setup.py build
            WORKING_DIRECTORY ${SETUPTOOLS_BUILD_DIR}
        )

        execute_process(
            COMMAND ${BUNDLE_STAGEDIR}/bin/python3 setup.py install
            WORKING_DIRECTORY ${SETUPTOOLS_BUILD_DIR}
        )
      endif()
    endif()

    if (EXISTS "${PIP_ZIPFILE}")
      set( PIP_BUILD_DIR "${BUNDLE_BUILDDIR}/pip" )

      if (NOT EXISTS "${PIP_BUILD_DIR}")
        message( STATUS "Unzipping ${PIP_ZIPFILE} to ${PIP_BUILD_DIR}" )
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${PIP_ZIPFILE}
            WORKING_DIRECTORY ${BUNDLE_BUILDDIR}
        )

        execute_process(
            COMMAND ${BUNDLE_STAGEDIR}/bin/python3 setup.py build
            WORKING_DIRECTORY ${PIP_BUILD_DIR}
        )

        execute_process(
            COMMAND ${BUNDLE_STAGEDIR}/bin/python3 setup.py install
            WORKING_DIRECTORY ${PIP_BUILD_DIR}
        )
      endif()
    endif()
  endif()
endif()

if ( ANACONDA_BUILD )
  set( PYTHON_LIBRARIES "${PYTHON_LIBRARY}" )
  if (MSVC)
    set( PYTHON_SITE_DIR "../../lib/site-packages" )
  else()
    set( PYTHON_SITE_DIR "../../lib/python${PYTHON_VERSION}/site-packages" )
  endif()

  if (WIN32)
    set( PYTHON_MODULE_EXTENSION ".pyd" )
  else()
    set( PYTHON_MODULE_EXTENSION ".so" )
    get_filename_component(PYTHON_INCLUDE_DIR "${PYTHON_LIBRARY}" DIRECTORY)
    get_filename_component(PYTHON_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}" DIRECTORY)
    set( PYTHON_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}/include")
    if (NOT MSVC)
      set( PYTHON_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}/python${PYTHON_VERSION}${PYTHON_ABIFLAGS}" )
    endif()
  endif()

  message( STATUS "Using anaconda/miniconda python in ${PYTHON_LIBRARIES} | ${PYTHON_INCLUDE_DIR}" )
  message( STATUS "Python modules will be installed to ${PYTHON_SITE_DIR}" )

  set( SIRE_FOUND_PYTHON TRUE )

elseif ( MSYS )
  set( PYTHON_LIBRARIES "${PYTHON_LIBRARY}" )

  # call python to tell it the site-packages directory itself
  execute_process(
     COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils import sysconfig as sc
print(sc.get_python_lib(prefix='', plat_specific=True))"
     OUTPUT_VARIABLE PYTHON_SITE
     OUTPUT_STRIP_TRAILING_WHITESPACE)

  set( PYTHON_SITE_DIR "${SIRE_INSTALL_PREFIX}/${PYTHON_SITE}" )
  set( PYTHON_MODULE_EXTENSION ".pyd" )

  message( STATUS "Using msys python in ${PYTHON_LIBRARIES} | ${PYTHON_INCLUDE_DIR}" )
  message( STATUS "Python modules will be installed to ${PYTHON_SITE_DIR}" )

elseif ( PYTHON_LIBRARY )
  set( PYTHON_LIBRARIES "${PYTHON_LIBRARY}" )
  set( PYTHON_INCLUDE_DIR "${BUNDLE_STAGEDIR}/include/python${PYTHON_VERSION}${PYTHON_ABIFLAGS}")
  set( PYTHON_SITE_DIR "${SIRE_BUNDLED_DIR}/lib/python${PYTHON_VERSION}/site-packages" )
  set( PYTHON_MODULE_EXTENSION ".so" )

  if (APPLE)
    execute_process( COMMAND chmod u+w ${PYTHON_LIBRARY} )
    execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libpython${PYTHON_VERSION}${PYTHON_ABIFLAGS}.dylib" ${PYTHON_LIBRARY} )
  endif()

  message( STATUS "Using bundled python in ${PYTHON_LIBRARIES} | ${PYTHON_INCLUDE_DIR}" )
  message( STATUS "Python modules will be installed to ${PYTHON_SITE_DIR}" )

  # Now lets unpack and install setuptools and pip, so that we can easily
  # install other python modules
  set( SETUPTOOLS_ZIPFILE "${CMAKE_SOURCE_DIR}/bundled/setuptools.tar.gz" )
  set( PIP_ZIPFILE "${CMAKE_SOURCE_DIR}/bundled/pip.tar.gz" )

  if (EXISTS "${SETUPTOOLS_ZIPFILE}")
    set( SETUPTOOLS_BUILD_DIR "${BUNDLE_BUILDDIR}/setuptools" )

    if (NOT EXISTS "${SETUPTOOLS_BUILD_DIR}")
      message( STATUS "Unzipping ${SETUPTOOLS_ZIPFILE} to ${SETUPTOOLS_BUILD_DIR}" )
      execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf ${SETUPTOOLS_ZIPFILE}
          WORKING_DIRECTORY ${BUNDLE_BUILDDIR}
      )

      execute_process(
          COMMAND ${BUNDLE_STAGEDIR}/bin/python3 setup.py build
          WORKING_DIRECTORY ${SETUPTOOLS_BUILD_DIR}
      )

      execute_process(
          COMMAND ${BUNDLE_STAGEDIR}/bin/python3 setup.py install
          WORKING_DIRECTORY ${SETUPTOOLS_BUILD_DIR}
      )
    endif()
  endif()

  if (EXISTS "${PIP_ZIPFILE}")
    set( PIP_BUILD_DIR "${BUNDLE_BUILDDIR}/pip" )

    if (NOT EXISTS "${PIP_BUILD_DIR}")
      message( STATUS "Unzipping ${PIP_ZIPFILE} to ${PIP_BUILD_DIR}" )
      execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf ${PIP_ZIPFILE}
          WORKING_DIRECTORY ${BUNDLE_BUILDDIR}
      )

      execute_process(
          COMMAND ${BUNDLE_STAGEDIR}/bin/python3 setup.py build
          WORKING_DIRECTORY ${PIP_BUILD_DIR}
      )

      execute_process(
          COMMAND ${BUNDLE_STAGEDIR}/bin/python3 setup.py install
          WORKING_DIRECTORY ${PIP_BUILD_DIR}
      )
    endif()
  endif()

  set( SIRE_FOUND_PYTHON TRUE )
else()
  message( STATUS "Strange? Cannot find the installed Python library. We cannot compile it, so will need to rely on the system version..." )
endif()
