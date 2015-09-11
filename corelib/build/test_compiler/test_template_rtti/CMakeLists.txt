########################################
#
# CMake file for test_template_rtti
#
########################################

include_directories(${CMAKE_SOURCE_DIR}/src/libs)

#create the library
add_library (libFoo libfoo.h libfoo.cpp)

set_target_properties (libFoo PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

target_link_libraries (libFoo
                       SireStream
                       )

add_library (libBar libbar.h libbar.cpp)

set_target_properties (libBar PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

target_link_libraries (libBar
                       libFoo
                      )

#now create the test application
add_executable( test_template_rtti main.cpp )

target_link_libraries( test_template_rtti libBar libFoo )

