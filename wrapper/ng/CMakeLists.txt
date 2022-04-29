########################################
#
# Set of ng (next generation) scripts to 
# simplify use of Sire
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _system.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/Sire/ng )
