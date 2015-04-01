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
        QuantumToMM.py
        OpenMMMD.py
        Plot.py
        WaterChanger.py
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
