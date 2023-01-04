
from .. import FF as _FF
from .. import Base as _Base

from ._System import *

System.__setProperty__ = System.setProperty
System.setProperty = _Base.__set_property__
