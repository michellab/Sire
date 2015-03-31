########################################
#
# CMake file for test_template_rtti
#
########################################

include_directories(${CMAKE_SOURCE_DIR}/src/libs)

#create the library
add_executable ( test_exception_packing main.cpp)

target_link_libraries( test_exception_packing SireStream )

