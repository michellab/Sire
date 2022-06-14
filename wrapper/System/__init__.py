#############################
##
## The SireSystem module
##
## (C) Christopher Woods
##

import Sire.FF

from Sire.System._System import *

def __system_getitem__(system, i):
    try:
        return system.__orig__getitem__(i)
    except Exception as e:
        if e.__class__.__name__ == "ArgumentError":
            return system.__orig__getitem__(Sire.Mol.MolIdx(i))
        else:
            raise e

System.__orig__getitem__ = System.__getitem__
System.__getitem__ = __system_getitem__
