########################################
#
# sire.system
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _system.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/system )
