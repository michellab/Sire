###
### CMake file to install libcpuid
### in the Sire bundle
###

include (DownloadAndCheckMD5)
# First, try to find the cpuid library in the "sire.app/lib" and "bundled/lib"
# directory. If it exists, then we don't need to do anything.
unset(CPUID_LIBRARY CACHE)

find_library( CPUID_LIBRARY "cpuid" PATHS ${SIRE_APP}/lib NO_DEFAULT_PATH )
if (EXISTS "${SIRE_APP}/include/libcpuid" AND CPUID_LIBRARY )
    message( STATUS "Conda version of libcpuid is already installed.")
    set( CPUID_INCLUDE_DIR "${SIRE_APP}/include" )
else()
    find_library( CPUID_LIBRARY "cpuid" PATHS ${BUNDLE_STAGEDIR}/lib NO_DEFAULT_PATH )

    if ( MSYS )
        set( CMAKE_MAKE_PROGRAM "mingw32-make" )
    endif()

    if ( CPUID_LIBRARY )
        message( STATUS "Have already compiled libcpuid")
    else()
        message( STATUS "Downloading, compiling and installing libcpuid" )

        if(NOT DEFINED CPUID_VERSION)
            set(CPUID_VERSION "0.4.1")
        endif()
        if(NOT DEFINED CPUID_URL)
            set(CPUID_URL "https://github.com/anrieff/libcpuid/archive/v${CPUID_VERSION}.tar.gz")
        endif()
        if(NOT DEFINED CPUID_MD5SUM)
            set(CPUID_MD5SUM "e9ffa7413e14039823ddd7452500ded9")
        endif()
        if(NOT DEFINED CPUID_ZIPFILE)
            set(CPUID_ZIPFILE "${CMAKE_SOURCE_DIR}/src/bundled/libcpuid-${CPUID_VERSION}.tar.gz")
        endif()
        if (NOT EXISTS "${CPUID_ZIPFILE}")
            downloadAndCheckMD5(${CPUID_URL} "${CPUID_ZIPFILE}" ${CPUID_MD5SUM})
        endif()

        if (EXISTS "${CPUID_ZIPFILE}")
            set( CPUID_BUILD_DIR "${BUNDLE_BUILDDIR}/libcpuid-${CPUID_VERSION}" )

            if (NOT EXISTS "${CPUID_BUILD_DIR}")
                message( STATUS "Unzipping ${CPUID_ZIPFILE} to ${CPUID_BUILD_DIR}" )
                execute_process(
                    COMMAND ${CMAKE_COMMAND} -E tar xzf ${CPUID_ZIPFILE}
                    WORKING_DIRECTORY ${BUNDLE_BUILDDIR}
                    OUTPUT_QUIET ERROR_QUIET
                )
            endif()

            if (MSVC)
                set ( MSVC_TOOLSET "" )
                if (CMAKE_GENERATOR MATCHES "^Visual Studio 14 2015.*$")
                    set ( MSVC_TOOLSET "v140" )
                elseif (CMAKE_GENERATOR MATCHES "^Visual Studio 15 2017.*$")
                    set ( MSVC_TOOLSET "v141" )
                elseif (CMAKE_GENERATOR MATCHES "^Visual Studio 16 2019.*$")
                    set ( MSVC_TOOLSET "v142" )
                else()
                    message ("Unknown Visual studio toolset")
                endif()
                if ( NOT MSVC_TOOLSET STREQUAL "" )
                    set ( MSVC_TOOLSET "/p:PlatformToolset=${MSVC_TOOLSET}" )
                endif()
                list( APPEND CPUID_OPTIONS "/m:${NCORES}" )
                list( APPEND CPUID_OPTIONS "/p:Configuration=Release" )
                list( APPEND CPUID_OPTIONS "/p:Platform=${CMAKE_VS_PLATFORM_NAME}" )
                list( APPEND CPUID_OPTIONS "${MSVC_TOOLSET}" )
                list( APPEND CPUID_OPTIONS "libcpuid_vc10.sln" )
                message( STATUS "${CMAKE_MAKE_PROGRAM} | ${CPUID_OPTIONS}" )
            else()
                list( APPEND CPUID_OPTIONS "--enable-static=no" )
                list( APPEND CPUID_OPTIONS "--enable-shared=yes" )
                list( APPEND CPUID_OPTIONS "--prefix=${BUNDLE_STAGEDIR}" )
                list( APPEND CPUID_OPTIONS "--libdir=${BUNDLE_STAGEDIR}/lib" )
                list( APPEND CPUID_OPTIONS "CC=${CMAKE_C_COMPILER}" )

                if (HAVE_STDINT_H)
                    list( APPEND CPUID_OPTIONS "CFLAGS=-DHAVE_STDINT_H" )
                endif()
                message( STATUS "${CPUID_BUILD_DIR}/configure | ${CPUID_OPTIONS}" )
            endif()

            if (MSYS)
                execute_process( COMMAND C:/msys64/usr/bin/sh.exe configure ${CPUID_OPTIONS}
                                WORKING_DIRECTORY ${CPUID_BUILD_DIR} )
            elseif (MSVC)
                message( STATUS "Patience... Compiling libcpuid..." )
                execute_process( COMMAND ${CMAKE_MAKE_PROGRAM} ${CPUID_OPTIONS}
                                WORKING_DIRECTORY ${CPUID_BUILD_DIR}
                            )
            else()
                message( STATUS "Patience... Configuring libcpuid..." )
                execute_process( COMMAND "libtoolize"
                                WORKING_DIRECTORY ${CPUID_BUILD_DIR}
                            )
                execute_process( COMMAND "autoreconf" "--install"
                                WORKING_DIRECTORY ${CPUID_BUILD_DIR}
                            )
                execute_process( COMMAND "${CPUID_BUILD_DIR}/configure" ${CPUID_OPTIONS}
                                WORKING_DIRECTORY ${CPUID_BUILD_DIR}
                            )
                message( STATUS "Patience... Compiling libcpuid..." )
                execute_process( COMMAND "${CMAKE_MAKE_PROGRAM}" -k -j ${NCORES}
                                WORKING_DIRECTORY ${CPUID_BUILD_DIR}
                            )
                message( STATUS "Patience... Installing libcpuid..." )
                execute_process( COMMAND ${CMAKE_MAKE_PROGRAM} install
                                WORKING_DIRECTORY ${CPUID_BUILD_DIR}
                            )
            endif()

            if (APPLE)
                set( CPUID_LIBRARY "${BUNDLE_STAGEDIR}/lib/libcpuid.dylib" )
                execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libcpuid.14.dylib" ${CPUID_LIBRARY} )
            elseif(MSYS)
                unset(CPUID_LIBRARY CACHE)
                find_library( CPUID_LIBRARY "cpuid" PATHS ${BUNDLE_STAGEDIR}/lib NO_DEFAULT_PATH )
            elseif(MSVC)
                unset(CPUID_LIBRARY CACHE)
                find_library( CPUID_LIBRARY "libcpuid.lib" PATHS "${CPUID_BUILD_DIR}/libcpuid/x64/Release" NO_DEFAULT_PATH )
                set( CPUID_INCLUDE_DIR ${CPUID_BUILD_DIR} )
            else()
                set( CPUID_LIBRARY "${BUNDLE_STAGEDIR}/lib/libcpuid.so" )
            endif()

            if (NOT EXISTS "${CPUID_LIBRARY}")
                message( FATAL_ERROR "Cannot find CPUID library ${CPUID_LIBRARY}. Error with compile?")
            endif()
        endif()
    endif()
endif()

if ( CPUID_LIBRARY )
  message( STATUS "Using libcpuid from ${CPUID_LIBRARY}" )

  if (NOT DEFINED CPUID_INCLUDE_DIR)
    set( CPUID_INCLUDE_DIR "${BUNDLE_STAGEDIR}/include")
  endif()

  if (HAVE_STDINT_H)
    set( CPUID_DEFINITIONS "-DHAVE_STDINT_H" )
  endif()

  set( SIRE_FOUND_CPUID TRUE )
else()
  message( STATUS "Strange? Cannot find the installed libcpuid. We cannot compile it, so will need to rely on the system version..." )
endif()
