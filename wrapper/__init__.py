"""
.. currentmodule:: Sire

This is the Sire module.

Functions
=========

.. autosummary::
    :toctree: generated/

    load
    save
    atomid
    chainid
    molid
    resid
    segid
    thumbs_up
    disable_thumbs_up
    get_thumbs_up_info

"""

__all__ = [ "load", "save",
            "atomid", "resid", "chainid",
            "segid", "molid" ]

from ._try_import import *

# ensure that the SireQt and SireError libraries are loaded as
# these are vital for the rest of the module
from . import Qt
from . import Error
from . import Config
from . import Base

from . import CAS
from . import FF
from . import ID
from . import IO
from . import MM
from . import Maths
from . import Move
from . import Mol
from . import Stream
from . import System
from . import Units
from . import Vol

qt = Qt
error = Error
config = Config
base = Base
cas = CAS
ff = FF
io = IO
mm = MM
maths = Maths
move = Move
mol = Mol
stream = Stream
system = System
units = Units
vol = Vol


def molid(num: int = None, name: str = None, idx: int = None):
    """Construct an identifer for a Molecule from the passed
       name, number and index.

    Args:
        name (str, optional): The molecule name. Defaults to None.
        num (int, optional): The molecule number. Defaults to None.
        idx (int, optional): The molecule index. Defaults to None.

    Returns:
        MolID : The returned molecule identifier
    """
    ID = None

    if type(num) is str:
        # used in unnamed argument mode
        if name is None:
            name = num
            num = None
        elif type(name) is int:
            (num, name) = (name, num)
        else:
            raise TypeError("The number cannot be a string.")

    if name is not None:
        ID = Mol.MolName(name)

    if num is not None:
        if ID is None:
            ID = Mol.MolNum(num)
        else:
            ID = ID + Mol.MolNum(num)

    if idx is not None:
        if ID is None:
            ID = Mol.MolIdx(idx)
        else:
            ID = ID + Mol.MolIdx(idx)

    if ID is None:
        return Mol.MolIdx()
    else:
        return ID

def atomid(num: int = None, name: str = None, idx: int = None):
    """Construct an identifer for an Atom from the passed
       name, number and index.

    Args:
        name (str, optional): The atom name. Defaults to None.
        num (int, optional): The atom number. Defaults to None.
        idx (int, optional): The atom index. Defaults to None.

    Returns:
        AtomID : The returned atom identifier
    """
    ID = None

    if type(num) is str:
        # used in unnamed argument mode
        if name is None:
            name = num
            num = None
        elif type(name) is int:
            (num, name) = (name, num)
        else:
            raise TypeError("The number cannot be a string.")

    if name is not None:
        ID = Mol.AtomName(name)

    if num is not None:
        if ID is None:
            ID = Mol.AtomNum(num)
        else:
            ID = ID + Mol.AtomNum(num)

    if idx is not None:
        if ID is None:
            ID = Mol.AtomIdx(idx)
        else:
            ID = ID + Mol.AtomIdx(idx)

    if ID is None:
        return Mol.AtomIdx()
    else:
        return ID


def resid(num: int = None, name: str = None, idx: int = None):
    """Construct an identifer for a Residue from the passed
       name, number and index.

    Args:
        name (str, optional): The residue name. Defaults to None.
        number (int, optional): The residue number. Defaults to None.
        index (int, optional): The residue index. Defaults to None.

    Returns:
        ResID : The returned atom identifier
    """
    ID = None

    if type(num) is str:
        # used in unnamed argument mode
        if name is None:
            name = num
            num = None
        elif type(name) is int:
            (num, name) = (name, num)
        else:
            raise TypeError("The number cannot be a string.")

    if name is not None:
        ID = Mol.ResName(name)

    if num is not None:
        if ID is None:
            ID = Mol.ResNum(num)
        else:
            ID = ID + Mol.ResNum(num)

    if idx is not None:
        if ID is None:
            ID = Mol.ResIdx(idx)
        else:
            ID = ID + Mol.ResIdx(idx)

    if ID is None:
        return Mol.ResIdx()
    else:
        return ID


def chainid(idx: int = None, name: str = None):
    """Construct an identifer for a Chain from the passed
       name and index.

    Args:
        name (str, optional): The chain name. Defaults to None.
        index (int, optional): The chain index. Defaults to None.

    Returns:
        ChainID : The returned chain identifier
    """
    ID = None

    if type(idx) is str:
        # used in unnamed argument mode
        if name is None:
            name = idx
            idx = None
        elif type(name) is int:
            (idx, name) = (name, idx)
        else:
            raise TypeError("The index cannot be a string.")

    if name is not None:
        ID = Mol.ChainName(name)

    if idx is not None:
        if ID is None:
            ID = Mol.ChainIdx(idx)
        else:
            ID = ID + Mol.ChainIdx(idx)

    if ID is None:
        return Mol.ChainIdx()
    else:
        return ID


def segid(idx: int = None, name: str = None):
    """Construct an identifer for a Segment from the passed
       name and index.

    Args:
        name (str, optional): The segment name. Defaults to None.
        index (int, optional): The segment index. Defaults to None.

    Returns:
        SegID : The returned chain identifier
    """
    ID = None

    if type(idx) is str:
        # used in unnamed argument mode
        if name is None:
            name = idx
            idx = None
        elif type(name) is int:
            (idx, name) = (name, idx)
        else:
            raise TypeError("The index cannot be a string.")

    if name is not None:
        ID = Mol.SegName(name)

    if idx is not None:
        if ID is None:
            ID = Mol.SegIdx(idx)
        else:
            ID = ID + Mol.SegIdx(idx)

    if ID is None:
        return Mol.SegIdx()
    else:
        return ID



__version__ = Config.__version__

__branch__ = Config.sire_repository_branch
__repository__ = Config.sire_repository_url
__revisionid__ = Config.sire_repository_version[0:7]


def _versionString():
    """Return a nicely formatted string that describes the current Sire version"""
    import Sire.Base as Base

    return """Sire %s [%s|%s, %s]""" % \
              (Base.getReleaseVersion(),
               Base.getRepositoryBranch(),
               Config.sire_repository_version[0:7],
               ["unclean", "clean"][Base.getRepositoryVersionIsClean()])


Config.versionString = _versionString

### Here are the functions and other data that form the public API
### of Sire

from ._load import *
from ._thumbsup import *
