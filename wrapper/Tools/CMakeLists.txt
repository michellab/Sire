########################################
#
# Set of Tool scripts to simplify
# use of Sire
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        AmberLoader.py
        DCDFile.py
        FDTISingleBound.py
        FDTISingleFree.py
        LSRC.py
        Nautilus.py 
        QuantumToMM.py
        OpenMMMD.py
        Plot.py
        WaterChanger.py
        WaterView.py
        WSRC.py
    )

# If your script needs any default files (e.g. parameter or config files)
# then put those files in a named directory in the Config directory and 
# add the name of that directory here. This directory will be copied
# to the share/Sire/Tools
set ( CONFIGDIRS
        WSRC
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/Sire/Tools )

install( DIRECTORY ${CONFIGDIRS} DESTINATION ${SIRE_SHARE}/Tools
         PATTERN .svn EXCLUDE )

# If this is an anaconda build, then make sure that all of the
# necessary files are copied from the pkgs directory into the
# anaconda lib and bin directories
# (do this here to make sure that it happens after installing
#  everything else)
if (ANACONDA_BUILD)
  # copy the files from the pkg directory to the anaconda lib and bin 
  # directories
  install( CODE "execute_process(COMMAND \"${ANACONDA_BASE}/bin/python\"
                      ${CMAKE_SOURCE_DIR}/build/copy_into_conda.py
                      ${CMAKE_INSTALL_PREFIX}/lib ${ANACONDA_BASE}/lib force
                            WORKING_DIRECTORY \"${ANACONDA_BASE}\")" )

  install( CODE "execute_process(COMMAND \"${ANACONDA_BASE}/bin/python\"
                      ${CMAKE_SOURCE_DIR}/build/copy_into_conda.py
                      ${CMAKE_INSTALL_PREFIX}/bundled/lib ${ANACONDA_BASE}/lib
                            WORKING_DIRECTORY \"${ANACONDA_BASE}\")" )

  install( CODE "execute_process(COMMAND \"${ANACONDA_BASE}/bin/python\"
                      ${CMAKE_SOURCE_DIR}/build/copy_into_conda.py
       	              ${CMAKE_INSTALL_PREFIX}/bin ${ANACONDA_BASE}/bin force
       	                    WORKING_DIRECTORY \"${ANACONDA_BASE}\")" )

endif()
