########################################
#
# sire.mol
#
########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _cursor.py
        _element.py
        _trajectory.py
        _view.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/mol )
