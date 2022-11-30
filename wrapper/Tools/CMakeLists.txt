########################################
#
# Set of Tool scripts to simplify
# use of Sire
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        ap.py
        AmberLoader.py
        BoreschAnalyticalCorrection.py
        BoreschNumericalCorrection.py
        DCDFile.py
        FreeEnergyAnalysis.py
        LJcutoff.py
        LSRC.py
        Nautilus.py
        QuantumToMM.py
        OpenMMMD.py
        Plot.py
        PSRC.py
        StandardState.py
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
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/legacy/Tools )

install( DIRECTORY ${CONFIGDIRS} DESTINATION ${SIRE_SHARE}/Tools )

# If this is an anaconda build, then make sure that all of the
# necessary files are copied from the pkgs directory into the
# anaconda lib, bin, and share directories
# (do this here to make sure that it happens after installing
#  everything else)
if(CMAKE_GENERATOR MATCHES "Visual Studio")  # MSBuild
  message( STATUS "Copying files from pkgs into Windows conda..." )

  # copy the files from the pkg directory to the anaconda
  # directories (Library/bin for DLLs and Library/lib for libs)
  install( CODE "execute_process(COMMAND \"${ANACONDA_BASE}/python\"
                    ${CMAKE_SOURCE_DIR}/build/copy_into_conda.py
                    ${CMAKE_INSTALL_PREFIX}/bin ${ANACONDA_BASE}/Library/bin force
                          WORKING_DIRECTORY \"${ANACONDA_BASE}\")" )

  install( CODE "execute_process(COMMAND \"${ANACONDA_BASE}/python\"
                          ${CMAKE_SOURCE_DIR}/build/copy_into_conda.py
                          ${CMAKE_INSTALL_PREFIX}/lib ${ANACONDA_BASE}/Library/lib force
                                WORKING_DIRECTORY \"${ANACONDA_BASE}\")" )

  install( CODE "execute_process(COMMAND \"${ANACONDA_BASE}/python\"
                                ${CMAKE_SOURCE_DIR}/build/copy_into_conda.py
                                      ${CMAKE_INSTALL_PREFIX}/share ${ANACONDA_BASE}/share force
                                            WORKING_DIRECTORY \"${ANACONDA_BASE}\")" )
else()
  message( STATUS "Copying files from pkgs into Linux/MacOS conda..." )

  # copy the files from the pkg directory to the anaconda lib and bin
  # directories
  install( CODE "execute_process(COMMAND \"${ANACONDA_BASE}/bin/python\"
                    ${CMAKE_SOURCE_DIR}/build/copy_into_conda.py
                    ${CMAKE_INSTALL_PREFIX}/lib ${ANACONDA_BASE}/lib force
                          WORKING_DIRECTORY \"${ANACONDA_BASE}\")" )

  install( CODE "execute_process(COMMAND \"${ANACONDA_BASE}/bin/python\"
                    ${CMAKE_SOURCE_DIR}/build/copy_into_conda.py
                          ${CMAKE_INSTALL_PREFIX}/bin ${ANACONDA_BASE}/bin force
                                WORKING_DIRECTORY \"${ANACONDA_BASE}\")" )

  install( CODE "execute_process(COMMAND \"${ANACONDA_BASE}/bin/python\"
                    ${CMAKE_SOURCE_DIR}/build/copy_into_conda.py
                          ${CMAKE_INSTALL_PREFIX}/share ${ANACONDA_BASE}/share force
                                WORKING_DIRECTORY \"${ANACONDA_BASE}\")" )
endif()
