########################################
#
# CMake file for test_qhash_lookup
#
########################################

include_directories(${CMAKE_SOURCE_DIR}/src/libs)

#create the library
add_executable ( test_qhash_lookup main.cpp )

target_link_libraries( test_qhash_lookup Qt5::Core )
