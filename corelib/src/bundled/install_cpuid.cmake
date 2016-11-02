###
### CMake file to install libcpuid
### in the Sire bundle
###

# First, try to find the cpuid library in the "bundled" directory. If it exists,
# then we don't need to do anything
unset(CPUID_LIBRARY CACHE)
find_library( CPUID_LIBRARY "cpuid" PATHS ${BUNDLE_STAGEDIR}/lib NO_DEFAULT_PATH )

if ( MSYS )
  set( CMAKE_MAKE_PROGRAM "mingw32-make" )
endif()

if ( CPUID_LIBRARY )
  message( STATUS "Have already compiled a bundled version of libcpuid")
else()
  message( STATUS "Compiling and installing a bundled version of libcpuid" )

  set( CPUID_ZIPFILE "${CMAKE_SOURCE_DIR}/src/bundled/libcpuid.tar.gz" )

  if (EXISTS "${CPUID_ZIPFILE}")
    set( CPUID_BUILD_DIR "${BUNDLE_BUILDDIR}/libcpuid" )

    if (NOT EXISTS "${CPUID_BUILD_DIR}")
      message( STATUS "Unzipping ${CPUID_ZIPFILE} to ${CPUID_BUILD_DIR}" )
      execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf ${CPUID_ZIPFILE}
          WORKING_DIRECTORY ${BUNDLE_BUILDDIR} 
          OUTPUT_QUIET ERROR_QUIET
      )
    endif()

    list( APPEND CPUID_OPTIONS "--enable-static=no" )
    list( APPEND CPUID_OPTIONS "--enable-shared=yes" )
    list( APPEND CPUID_OPTIONS "--prefix=${BUNDLE_STAGEDIR}" )
    list( APPEND CPUID_OPTIONS "CC=${CMAKE_C_COMPILER}" )

    if (HAVE_STDINT_H)
      list( APPEND CPUID_OPTIONS "CFLAGS=-DHAVE_STDINT_H" )
    endif()

    message( STATUS "${CPUID_BUILD_DIR}/configure | ${CPUID_OPTIONS}" )

    message( STATUS "Patience... Configuring libcpuid..." )
    if (MSYS)
      execute_process( COMMAND C:/msys64/usr/bin/sh.exe configure ${CPUID_OPTIONS}
                       WORKING_DIRECTORY ${CPUID_BUILD_DIR} )
    else()
      execute_process( COMMAND ${CPUID_BUILD_DIR}/configure ${CPUID_OPTIONS}
                       WORKING_DIRECTORY ${CPUID_BUILD_DIR}
                       OUTPUT_QUIET ERROR_QUIET 
                     )
    endif()

    message( STATUS "Patience... Compiling libcpuid..." )
    execute_process( COMMAND ${CMAKE_MAKE_PROGRAM} -k -j ${NCORES}
                     WORKING_DIRECTORY ${CPUID_BUILD_DIR}
                     OUTPUT_QUIET ERROR_QUIET 
                   )

    message( STATUS "Patience... Installing libcpuid..." )
    execute_process( COMMAND ${CMAKE_MAKE_PROGRAM} install
                     WORKING_DIRECTORY ${CPUID_BUILD_DIR}
                     OUTPUT_QUIET ERROR_QUIET 
                   )

    if (APPLE)
      set( CPUID_LIBRARY "${BUNDLE_STAGEDIR}/lib/libcpuid.dylib" )
      execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libcpuid.10.dylib" ${CPUID_LIBRARY} )
    else()
      set( CPUID_LIBRARY "${BUNDLE_STAGEDIR}/lib/libcpuid.so" )
    endif()

    if (NOT EXISTS "${CPUID_LIBRARY}")
      message( FATAL_ERROR "Cannot find CPUID library ${CPUID_LIBRARY}. Error with compile?")
    endif()

  endif()
endif()

if ( CPUID_LIBRARY )
  message( STATUS "Using libcpuid from ${CPUID_LIBRARY}" )

  set( CPUID_INCLUDE_DIR "${BUNDLE_STAGEDIR}/include")

  if (HAVE_STDINT_H)
    set( CPUID_DEFINITIONS "-DHAVE_STDINT_H" )
  endif()

  set( SIRE_FOUND_CPUID TRUE )
else()
  message( STATUS "Strange? Cannot find the installed libcpuid. We cannot compile it, so will need to rely on the system version..." )
endif()
