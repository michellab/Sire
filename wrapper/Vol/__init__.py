#############################
##
## The SireVol module
##
## This contains all of the classes that provide
## volumes in which the simulation is set up and 
## run (e.g. infinite boxes - 'Cartesian', and 
## periodic boxes - 'PeriodicBox')
##

import Sire.Base
import Sire.Maths
import Sire.Units

# Import all of the classes and functions from the C++ library
from Sire.Vol._Vol import *

# Now define some pure Python functions and classes that are part of 
# this library...


