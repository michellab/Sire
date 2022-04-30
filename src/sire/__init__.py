"""
.. currentmodule:: sire

This is the sire module.

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
            "segid", "molid", "thumbs_up",
            "get_thumbs_up_info", "disable_thumbs_up" ]


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

    from .mol import MolName, MolNum, MolIdx

    if name is not None:
        ID = MolName(name)

    if num is not None:
        if ID is None:
            ID = MolNum(num)
        else:
            ID = ID + MolNum(num)

    if idx is not None:
        if ID is None:
            ID = MolIdx(idx)
        else:
            ID = ID + MolIdx(idx)

    if ID is None:
        return MolIdx()
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

    from .mol import AtomName, AtomNum, AtomIdx

    if name is not None:
        ID = AtomName(name)

    if num is not None:
        if ID is None:
            ID = AtomNum(num)
        else:
            ID = ID + AtomNum(num)

    if idx is not None:
        if ID is None:
            ID = AtomIdx(idx)
        else:
            ID = ID + AtomIdx(idx)

    if ID is None:
        return AtomIdx()
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

    from .mol import ResName, ResNum, ResIdx

    if name is not None:
        ID = ResName(name)

    if num is not None:
        if ID is None:
            ID = ResNum(num)
        else:
            ID = ID + ResNum(num)

    if idx is not None:
        if ID is None:
            ID = ResIdx(idx)
        else:
            ID = ID + ResIdx(idx)

    if ID is None:
        return ResIdx()
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

    from .mol import ChainName, ChainIdx

    if name is not None:
        ID = ChainName(name)

    if idx is not None:
        if ID is None:
            ID = ChainIdx(idx)
        else:
            ID = ID + ChainIdx(idx)

    if ID is None:
        return ChainIdx()
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

    from .mol import SegName, SegIdx

    if name is not None:
        ID = SegName(name)

    if idx is not None:
        if ID is None:
            ID = SegIdx(idx)
        else:
            ID = ID + SegIdx(idx)

    if ID is None:
        return SegIdx()
    else:
        return ID


from . import config

__version__ = config.__version__

__branch__ = config.sire_repository_branch
__repository__ = config.sire_repository_url
__revisionid__ = config.sire_repository_version[0:7]

_disable_lazy_import = False

try:
    if _disable_lazy_import:
        raise AssertionError()

    import lazy_import as _lazy_import
    _lazy_import.logging.disable(_lazy_import.logging.DEBUG)
except Exception:
    class _lazy_import:
        """This is not lazy_import, but instead a thin stub that matches the
           API but DOES NOT lazy_import anything. This imports at call time.
        """
        @staticmethod
        def lazy_module(m):
            from importlib import import_module
            return import_module(m, package="sire")

        @staticmethod
        def lazy_function(f):
            module_name, unit_name = f.rsplit('.', 1)
            module = _lazy_import.lazy_module(module_name)
            return getattr(module, unit_name)

        @staticmethod
        def lazy_class(c):
            return _lazy_import.lazy_function(c)


# Lazy import the modules for speed, and also to prevent pythonizing them
# if the users wants to run in legacy mode

analysis = _lazy_import.lazy_module("sire.analysis")
base = _lazy_import.lazy_module("sire.base")
cas = _lazy_import.lazy_module("sire.cas")
cluster = _lazy_import.lazy_module("sire.cluster")
error = _lazy_import.lazy_module("sire.error")
ff = _lazy_import.lazy_module("sire.ff")
id = _lazy_import.lazy_module("sire.id")
io = _lazy_import.lazy_module("sire.io")
maths = _lazy_import.lazy_module("sire.maths")
mm = _lazy_import.lazy_module("sire.mm")
mol = _lazy_import.lazy_module("sire.mol")
move = _lazy_import.lazy_module("sire.move")
qt = _lazy_import.lazy_module("sire.qt")
squire = _lazy_import.lazy_module("sire.squire")
stream = _lazy_import.lazy_module("sire.stream")
units = _lazy_import.lazy_module("sire.units")
vol = _lazy_import.lazy_module("sire.vol")

def _version_string():
    """Return a nicely formatted string that describes the current Sire version"""
    from .base import get_release_version, get_repository_branch, \
        get_repository_version_is_clean

    from .config import sire_repository_version

    return """Sire %s [%s|%s, %s]""" % \
              (get_release_version(),
               get_repository_branch(),
               sire_repository_version[0:7],
               ["unclean", "clean"][get_repository_version_is_clean()])

config.version_string = _version_string

### Here are the functions and other data that form the public API
### of Sire

from ._load import *
from ._thumbsup import *
