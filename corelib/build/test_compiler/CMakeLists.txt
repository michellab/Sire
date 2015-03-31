###############################
#
# CMake file to compile the
# compiler tests
#
###############################

# Add each test
add_subdirectory (test_template_rtti)
add_subdirectory (test_error_packing)
add_subdirectory (test_qhash_lookup)
add_subdirectory (test_stream)
add_subdirectory (test_openmp)
add_subdirectory (test_mpi)
