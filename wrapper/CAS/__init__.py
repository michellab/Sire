#############################
##
## The SireCAS library.
##
## This module provides a basic CAS (computer algebra system)
## This can be used to provide user-defined functions which 
## can be used to calculate energies and forces, or just
## for your own amusement :-)
##

import Sire.Maths

# Import all of the Qt classes
import Sire.Qt

# Import all of the classes and functions from the C++ library
from Sire.CAS._CAS import *

# Now define some pure Python functions and classes that are part of
# this library...

#enable ** operator for exbase types
ExBase.__pow__ = pow
Expression.__pow__ = pow
