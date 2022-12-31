
__all__ = ["get_coords_array"]

from ..legacy import IO as _IO
from .. import use_new_api as _use_new_api
_use_new_api()


def get_coords_array(mol, units=None, map=None):
    """Return the coordinates of the passed molecule view as a
       numpy array of shape (natoms,3). Specify the length
       units to use, and optionally pass in a map to find
       the coordinates property
    """
    import numpy as np

    if units is None:
        from ..units import angstrom
        units = angstrom

    from ..base import create_map
    map = create_map(map)

    if hasattr(mol, "to_molecule_group"):
        mol = mol.to_molecule_group()

    coords = _IO.getCoordsArray(mol, units, map)

    natoms = int(len(coords)/3)

    return np.reshape(np.asarray(coords, dtype=float), (natoms, 3))


def load_molecules(*args, **kwargs):
    from ..legacy.IO import load_molecules as _load_molecules
    from ..system import System
    return System(_load_molecules(*args, **kwargs))
