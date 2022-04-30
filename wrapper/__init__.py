"""
.. currentmodule:: sire.legacy

This is the legacy Sire module. This holds all of the C++ modules
on which the whole of Sire is built. These are designed to be used
by developers, and are not part of the stable public API
(which is sire)

"""

# ensure that the SireQt and SireError libraries are loaded as
# these are vital for the rest of the module
from . import Qt
from . import Error
from . import Config

__version__ = Config.__version__

__branch__ = Config.sire_repository_branch
__repository__ = Config.sire_repository_url
__revisionid__ = Config.sire_repository_version[0:7]


def _versionString():
    """Return a nicely formatted string that describes the current Sire version"""
    from .Base import getReleaseVersion, getRepositoryBranch, \
        getRepositoryVersionIsClean

    from .Config import sire_repository_version

    return """Sire %s [%s|%s, %s]""" % \
              (getReleaseVersion(),
               getRepositoryBranch(),
               sire_repository_version[0:7],
               ["unclean", "clean"][getRepositoryVersionIsClean()])


Config.versionString = _versionString
