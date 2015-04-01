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

import Sire.Maths

def wrap(value):
    try:
        return Sire.Base._Base.wrap(value)
    except:
        return Sire.Maths.wrap(value)
