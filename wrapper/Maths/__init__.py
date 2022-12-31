
from .. import Qt as _Qt
from .. import Error as _Error
from .. import Base as _Base

# Import all of the classes and functions from the C++ library
from ._Maths import *

# Now define some pure Python functions and classes that are part of
# this library...

wrap = _Base._add_wrap_function(wrap)

# No QVector<float> exposed (would have horrible casting bugs)
MultiFloat.toArray = staticmethod( MultiFloat.toDoubleArray )
