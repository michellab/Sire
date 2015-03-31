###
### CMake file to install QtCore from Qt 5
### in the Sire bundle
###

# First, try to find the QtCore library in the "bundled" directory. If it exists,
# then we don't need to do anything
set( Qt5Core_DIR "${BUNDLE_STAGEDIR}/lib/cmake/Qt5Core" )
find_package( Qt5Core QUIET )

set ( NEED_BUILD_QT TRUE )

if ( SIRE_BUILD_GUI )
  message( STATUS "Need also to build QtGui to create Sire GUIs" )

  set( Qt5Widgets_DIR "${BUNDLE_STAGEDIR}/lib/cmake/Qt5Widgets" )
  find_package( Qt5Widgets QUIET )

  if ( Qt5Core_FOUND AND Qt5Widgets_FOUND )
    set ( NEED_BUILD_QT FALSE )
    message( STATUS "Have already compiled a bundled version of Qt5Core and Qt5Widgets" )
  endif()
else()
  if ( Qt5Core_FOUND )
    message( STATUS "Have already compiled a bundled version of Qt5Core" )
    set ( NEED_BUILD_QT FALSE )
  endif()
endif()

if ( NEED_BUILD_QT )
  message( STATUS "Compiling and installing a bundled version of QtCore from Qt5" )
  if ( SIRE_BUILD_GUI )
    message( STATUS "Compiling and installing a bundled version of QtWidgets from Qt5" )
  endif()

  message( STATUS "THIS WILL TAKE A VERY LONG TIME" )

  set( QT_ZIPFILE "${CMAKE_SOURCE_DIR}/src/bundled/qtbase.tar.gz" )
  set( QT_PATCHFILE "${CMAKE_SOURCE_DIR}/src/bundled/qtbase_patch.tar.gz" )

  if (EXISTS "${QT_ZIPFILE}")
    set( QT_BUILD_DIR "${BUNDLE_BUILDDIR}/qtbase" )

    if (NOT EXISTS "${QT_BUILD_DIR}")
      message( STATUS "Unzipping ${QT_ZIPFILE} to ${QT_BUILD_DIR}" )
      execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf ${QT_ZIPFILE}
          WORKING_DIRECTORY ${BUNDLE_BUILDDIR}
      )

      if (EXISTS "${QT_PATCHFILE}")
        message( STATUS "Unzipping ${QT_PATCHFILE} to ${QT_BUILD_DIR}" )
        execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf ${QT_PATCHFILE}
          WORKING_DIRECTORY ${BUNDLE_BUILDDIR}
        )
      endif()
    endif()

    list( APPEND QT_OPTIONS "-no-javascript-jit")
    list( APPEND QT_OPTIONS "-no-glib")
    list( APPEND QT_OPTIONS "-no-kms")
    list( APPEND QT_OPTIONS "-no-dbus")
    list( APPEND QT_OPTIONS "-no-pch")
    list( APPEND QT_OPTIONS "-no-icu")
    list( APPEND QT_OPTIONS "-no-iconv")
    list( APPEND QT_OPTIONS "-no-cups")
    list( APPEND QT_OPTIONS "-no-nis")
    list( APPEND QT_OPTIONS "-no-rpath")
    list( APPEND QT_OPTIONS "-no-linuxfb")
    list( APPEND QT_OPTIONS "-no-directfb")
    list( APPEND QT_OPTIONS "-no-eglfs")
    list( APPEND QT_OPTIONS "-no-xcb")
    list( APPEND QT_OPTIONS "-qt-pcre")
    list( APPEND QT_OPTIONS "-no-openssl")
    list( APPEND QT_OPTIONS "-no-gif")
    list( APPEND QT_OPTIONS "-qt-zlib")
    list( APPEND QT_OPTIONS "-no-pkg-config")

    if ( SIRE_BUILD_GUI )
      message( STATUS "Configuring Qt5 to build the graphical and OpenGL libraries..." )
    else()
      list( APPEND QT_OPTIONS "-no-widgets")
      list( APPEND QT_OPTIONS "-no-gui")
      list( APPEND QT_OPTIONS "-no-opengl")
    endif()

    if (SIRE_HAS_CPP_11)
      list( APPEND QT_OPTIONS "-c++11")

      if (NEED_UNDEF_STRICT_ANSI)
        set(ENV{CXXFLAGS} "-U__STRICT_ANSI__")
      endif()
    endif()

    list( APPEND QT_OPTIONS "-nomake;examples")
    list( APPEND QT_OPTIONS "-confirm-license")
    list( APPEND QT_OPTIONS "-opensource")
    list( APPEND QT_OPTIONS "-prefix;${BUNDLE_STAGEDIR}")

    if (APPLE)
      list( APPEND QT_OPTIONS "-no-framework")
    endif()

    message( STATUS "${QT_OPTIONS}" )

    message( STATUS "Patience... Configuring QtCore..." )
    execute_process( COMMAND ${QT_BUILD_DIR}/configure ${QT_OPTIONS}
                     WORKING_DIRECTORY ${QT_BUILD_DIR}
                     RESULT_VARIABLE QT_CONFIGURE_FAILED )

    if (QT_CONFIGURE_FAILED)
      message( FATAL_ERROR "Cannot configure Qt5. Please go to the mailing list for help.")
    endif()

    message( STATUS "Patience... Compiling QtCore..." )
    execute_process( COMMAND ${CMAKE_MAKE_PROGRAM} -k -j ${NCORES}
                     WORKING_DIRECTORY ${QT_BUILD_DIR}
                     RESULT_VARIABLE QT_BUILD_FAILED )

    if (QT_BUILD_FAILED)
      message( FATAL_ERROR "Cannot build Qt5. Please go to the mailing list for help.")
    endif()

    message( STATUS "Patience... Installing QtCore..." )
    execute_process( COMMAND ${CMAKE_MAKE_PROGRAM} install
                     WORKING_DIRECTORY ${QT_BUILD_DIR}
                     RESULT_VARIABLE QT_INSTALL_FAILED )

    if (QT_INSTALL_FAILED)
      message( FATAL_ERROR "Cannot install Qt5. Please go to the mailing list for help.")
    endif()

    set( Qt5Core_DIR "${BUNDLE_STAGEDIR}/lib/cmake/Qt5Core" )
    find_package( Qt5Core )

    if ( Qt5Core_FOUND )
      # need to set the install name so that we can find the library when it is 
      # placed into the bundle directory
      if (APPLE)
        execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libQt5Core.dylib" ${BUNDLE_STAGEDIR}/lib/libQt5Core.dylib )
        execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libQt5Core_debug.dylib" ${BUNDLE_STAGEDIR}/lib/libQt5Core_debug.dylib )
      endif()
    endif()

    if ( SIRE_BUILD_GUI )

      set( Qt5Gui_DIR "${BUNDLE_STAGEDIR}/lib/cmake/Qt5Gui" )
      find_package( Qt5Gui )
      if ( Qt5Gui_FOUND )
        if (APPLE)
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libQt5Gui.dylib" ${BUNDLE_STAGEDIR}/lib/libQt5Gui.dylib )
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libQt5Gui_debug.dylib" ${BUNDLE_STAGEDIR}/lib/libQt5Gui_debug.dylib )
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Core.5.dylib @rpath/libQt5Core.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5Gui.dylib ) 
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Core_debug.5.dylib @rpath/libQt5Core_debug.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5Gui_debug.dylib ) 
        endif()
      endif()

      set( Qt5Widgets_DIR "${BUNDLE_STAGEDIR}/lib/cmake/Qt5Widgets" )
      find_package( Qt5Widgets )
      if ( Qt5Widgets_FOUND )
        if (APPLE)
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libQt5Widgets.dylib" ${BUNDLE_STAGEDIR}/lib/libQt5Widgets.dylib )
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libQt5Widgets_debug.dylib" ${BUNDLE_STAGEDIR}/lib/libQt5Widgets_debug.dylib )
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Core.5.dylib @rpath/libQt5Core.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5Widgets.dylib ) 
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Gui.5.dylib @rpath/libQt5Gui.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5Widgets.dylib ) 
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Core_debug.5.dylib @rpath/libQt5Core_debug.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5Widgets_debug.dylib ) 
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Gui_debug.5.dylib @rpath/libQt5Gui_debug.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5Widgets_debug.dylib ) 
        endif()
      endif()

      set( Qt5PrintSupport_DIR "${BUNDLE_STAGEDIR}/lib/cmake/Qt5PrintSupport" )
      find_package( Qt5PrintSupport )
      if ( Qt5PrintSupport_FOUND )
        if (APPLE)
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libQt5PrintSupport.dylib" ${BUNDLE_STAGEDIR}/lib/libQt5PrintSupport.dylib )
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/libQt5PrintSupport_debug.dylib" ${BUNDLE_STAGEDIR}/lib/libQt5PrintSupport_debug.dylib )
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Core.5.dylib @rpath/libQt5Core.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5PrintSupport.dylib ) 
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Gui.5.dylib @rpath/libQt5Gui.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5PrintSupport.dylib ) 
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Widgets.5.dylib @rpath/libQt5Widgets.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5PrintSupport.dylib ) 
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Core_debug.5.dylib @rpath/libQt5Core_debug.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5PrintSupport_debug.dylib ) 
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Gui_debug.5.dylib @rpath/libQt5Gui_debug.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5PrintSupport_debug.dylib ) 
          execute_process( COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change libQt5Widgets_debug.5.dylib @rpath/libQt5Widgets_debug.5.dylib ${BUNDLE_STAGEDIR}/lib/libQt5PrintSupport_debug.dylib ) 
        endif()
      endif()

    endif()
  endif()
endif()

if ( SIRE_BUILD_GUI )
  if ( Qt5Core_FOUND AND Qt5Widgets_FOUND )
    message( STATUS "Using Qt5 from ${BUNDLE_STAGEDIR}" )
    get_target_property(QtCore_location Qt5::Core LOCATION)
    get_target_property(QtWidgets_location Qt5::Widgets LOCATION)
    message( STATUS "Library: ${QtCore_location}" )
    message( STATUS "Library: ${QtWidgets_location}" )

    message( STATUS "Sire will also use QtGui and QtPrintSupport..." )

    set( Qt5Gui_DIR "${BUNDLE_STAGEDIR}/lib/cmake/Qt5Gui" )
    find_package( Qt5Gui )
    set( Qt5PrintSupport_DIR "${BUNDLE_STAGEDIR}/lib/cmake/Qt5PrintSupport" )
    find_package( Qt5PrintSupport )

    get_target_property(QtGui_location Qt5::Gui LOCATION)
    get_target_property(QtPrint_location Qt5::PrintSupport LOCATION)
    message( STATUS "Library: ${QtGui_location}")
    message( STATUS "Library: ${QtPrint_location}")

    set( SIRE_FOUND_QT TRUE )
  else()
    message( STATUS "Strange? Cannot find the file containing the cmake for Qt5.
                     We cannot compile it, so will need to rely on the system version..." )
  endif()
else()
  if ( Qt5Core_FOUND )
    message( STATUS "Using Qt5Core from ${BUNDLE_STAGEDIR}" )
    get_target_property(QtCore_location Qt5::Core LOCATION)
    message( STATUS "Library: ${QtCore_location}" )
    set( SIRE_FOUND_QT TRUE )
  else()
    message( STATUS "Strange? Cannot find the file containing the cmake for Qt5Core.
                     We cannot compile it, so will need to rely on the system version..." )
  endif()
endif()
