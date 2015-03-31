###############################
#
# CMake file to compile the
# Sire applications
#
###############################

# Add each application
add_subdirectory (sire)
add_subdirectory (testmulti)
add_subdirectory (test_system)

if (SIRE_BUILD_GUI)
  add_subdirectory (inquire)
endif()
