#############################
##
## The SireError module
##
## (C) Christopher Woods
##

import Sire.Qt

from Sire.Error._Error import *

__old_printError = printError

def printError(e):
    __old_printError( str(e) )

