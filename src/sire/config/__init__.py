__all__ = [
    "binary_directory",
    "include_directory",
    "library_directory",
    "parameter_directory",
    "share_directory",
    "sire_repository_url",
    "sire_repository_version",
    "sire_repository_branch",
    "version_string",
]

from ..legacy import Config as _Config

from ..legacy.Config import (
    binary_directory,
    include_directory,
    library_directory,
    parameter_directory,
    share_directory,
    sire_repository_url,
    sire_repository_version,
    sire_repository_branch,
)

version_string = _Config.versionString

__version__ = _Config.__version__
