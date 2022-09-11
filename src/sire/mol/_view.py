
try:
    import nglview as _nglview
    has_nglview = True
except ImportError:
    has_nglview = False
except AttributeError as e:
    has_nglview = False
    print(f"Failed to import nglview: {e}")


__all__ = ["view"]


if has_nglview:

    @_nglview.register_backend('sire')
    class _SireStructureTrajectory(_nglview.Trajectory, _nglview.Structure):
        def __init__(self, obj=None, sire_obj=None, map=None):
            if sire_obj is not None:
                self._obj = sire_obj.__copy__()
                from ..base import PropertyMap

                if map is None:
                    map = PropertyMap()
                else:
                    map = PropertyMap(map)

                self._map = map
            else:
                self._obj = obj._obj
                self._map = obj._map

            self.ext = 'pdb'
            self.params = {}
            import uuid
            self.id = str(uuid.uuid4())

        def get_structure_string(self):
            from ..legacy.IO import PDB2

            molecules = self._obj

            if molecules.what() != "SireSystem::System":
                from ..legacy.System import System
                from ..legacy.Mol import MoleculeGroup
                s = System()

                try:
                    m = molecules.to_molecule_group()
                except AttributeError:
                    m = MoleculeGroup("all")
                    m.add(molecules)

                s.add(m)
                molecules = s

            pdb2 = PDB2(molecules, map=self._map)

            lines = pdb2.to_lines()

            s = "\n".join(lines)

            return s

        @property
        def n_frames(self):
            return max(1, self._obj.num_frames())

        def get_coordinates(self, index):
            self._obj.load_frame(index)
            coords = self._obj.property("coordinates")
            import numpy as np
            c = np.zeros(shape=(coords.num_atoms(), 3), dtype=float)

            for i, coord in enumerate(coords):
                c[i] = (coord[0], coord[1], coord[2])

            print(c)

            return c

    def view(obj):
        struc_traj = _SireStructureTrajectory(sire_obj=obj)
        view = _nglview.NGLWidget(struc_traj)
        return view

else:
    def view(obj):
        raise ImportError(
            "You need to install nglview to be able to view "
            "molecules. Do this by typing, e.g. "
            "'mamba install nglview'")
