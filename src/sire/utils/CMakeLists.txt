###########################################
#
# Set of utility (util) scripts / functions
#
###########################################

# Add your script to this list
set ( SCRIPTS
        __init__.py
        _try_import.py
        _console.py
        _profiler.py
        _simple.py
        _spring_flowers.py
    )

# installation
install( FILES ${SCRIPTS} DESTINATION ${SIRE_PYTHON}/sire/utils )
