#############################
##
## The SireMaths library.
##
## This module provides lots of maths functions and classes,
## including geometric classes (vector, matrix and quaternion)
## and the algabraic maths engine (MathFunc and derivatives)
##

import Sire.Qt
import Sire.Error
import Sire.Base

# Import all of the classes and functions from the C++ library
from Sire.Maths._Maths import *

# Now define some pure Python functions and classes that are part of
# this library...

# No QVector<float> exposed (would have horrible casting bugs)
MultiFloat.toArray = staticmethod( MultiFloat.toDoubleArray )
