#############################
##
## The SireFF library.
##
## This contains all of the classes that are used to 
## provide a forcefield (energy+forces) for molecules.
##

import Sire.Mol
import Sire.CAS

# Import all of the classes and functions from the C++ library
from Sire.FF._FF import *

# Now define some pure Python functions and classes that are part of
# the module

