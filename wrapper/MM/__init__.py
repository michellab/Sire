#############################
##
## The SireMM library.
##
## This contains all of the classes that are used to
## provide a molecular mechanics forcefield (partial
## charges, LJ terms, bond, angle, dihedral terms,
## MM parameter database classes etc). It also 
## contains all of the MM forcefields.
##

import Sire.FF
import Sire.CAS

# Import all of the classes and functions from the C++ library
from Sire.MM._MM import *

# Now define some pure Python functions and classes that are part of 
# this library...

###### PROPERTY KLUDGE FIX

__props = [ AtomLJs ]

for __prop in __props:
    Sire.Mol._pvt_property_cludge_fix(__prop)
