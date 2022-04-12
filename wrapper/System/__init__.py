#############################
##
## The SireSystem module
##
## (C) Christopher Woods
##

import Sire.FF

from Sire.System._System import *

System.__setProperty__ = System.setProperty
System.setProperty = Sire.Base.__set_property__
