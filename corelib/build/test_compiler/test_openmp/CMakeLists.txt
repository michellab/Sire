########################################
#
# CMake file for test_stream
#
########################################

include_directories(${CMAKE_SOURCE_DIR}/src/libs)

#create the library
add_executable ( test_openmp main.cpp)

target_link_libraries (test_openmp
                       Qt5::Core
                      )

