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

if (MSVC)
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

set( PYTHON_LIBRARIES "${PYTHON_LIBRARY}" )
if (CMAKE_GENERATOR MATCHES "Visual Studio")  # MSBuild
  set( PYTHON_SITE_DIR "../../lib/site-packages" )
else()
  set( PYTHON_SITE_DIR "../../lib/python${PYTHON_VERSION}/site-packages" )
endif()

if(CMAKE_GENERATOR MATCHES "Visual Studio")  # MSBuild
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
