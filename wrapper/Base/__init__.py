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
    for func in _wrap_functions:
        try:
            return func(value)
        except:
            pass            

    return _base_wrap(value)

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
