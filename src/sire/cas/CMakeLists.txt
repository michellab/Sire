########################################
#
# sire.cas
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/cas )
