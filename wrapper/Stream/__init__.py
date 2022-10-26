"""
.. currentmodule:: sire.legacy.Stream

Classes
=======

.. autosummary::
    :toctree: generated/

    FileHeader
    MD5Sum

Functions
=========

.. autosummary::
    :toctree: generated/

    load
    save

"""

from ._Stream import *

import sys

_pvt_load = load

_pvt_modules = { "SireAnalysis" : "sire.legacy.Analysis",
                 "SireBase"     : "sire.legacy.Base",
                 "SireCAS"      : "sire.legacy.CAS",
                 "SireCluster"  : "sire.legacy.Cluster",
                 "SireError"    : "sire.legacy.Error",
                 "SireFF"       : "sire.legacy.FF",
                 "SireID"       : "sire.legacy.ID",
                 "SireIO"       : "sire.legacy.IO",
                 "SireMM"       : "sire.legacy.MM",
                 "SireMaths"    : "sire.legacy.Maths",
                 "SireMol"      : "sire.legacy.Mol",
                 "SireMove"     : "sire.legacy.Move",
                 "SireSystem"   : "sire.legacy.System",
                 "SireUnits"    : "sire.legacy.Units",
                 "SireVol"      : "sire.legacy.Vol",
                 "Squire"       : "sire.legacy.Squire",
                 "Soiree"       : "sire.legacy.Analysis"  # Soiree was renamed as Analysis
                }

def _pvt_loadLibrary(lib):

    lib = str(lib)

    if lib in _pvt_modules:
        __import__( _pvt_modules[lib] )

def load(data):
    header = getDataHeader(data)

    for lib in header.requiredLibraries():
        _pvt_loadLibrary(lib)

    return _pvt_load(data)

