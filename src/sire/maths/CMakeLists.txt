########################################
#
# sire.maths
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _sphere.py
        _vector.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/maths )
