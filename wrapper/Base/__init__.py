#############################
##
## The SireBase module
##
## This contains some base classes
## that don't fit anywhere else.
##
## This module may be removed
## as I find homes for some of
## these classes
##
## (C) Christopher Woods
##

from Sire.Base._Base import *

_wrap_functions = []

_base_wrap = wrap

def wrap(value):
    # First, try to wrap the python concrete classes
    if isinstance(value, bool):
        return BooleanProperty(value)

    elif isinstance(value, int) or isinstance(value, float):
        return NumberProperty(value)

    elif isinstance(value, str):
        return StringProperty(value)

    for func in _wrap_functions:
        try:
            return func(value)
        except:
            pass            

    try:
        return _base_wrap(value)
    except:
        pass

    return PropertyList(value)

_original_wrap = wrap

def _add_wrap_function(func):
    _wrap_functions.append(func)
    return _original_wrap

# cludgy quick fix for an anaconda install
_getBundledLibDir = getBundledLibDir
def getBundledLibDir():
    try:
        return _getBundledLibDir()
    except:
        return "%s/lib" % getInstallDir()
