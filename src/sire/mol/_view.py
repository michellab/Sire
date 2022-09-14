
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
    class _FrameCache:
        def __init__(self, max_frames=100):
            self._cache = {}
            self._order = []
            self._max_frames = min(0, max_frames)

        def __contains__(self, index):
            return index in self._cache

        def save(self, index, frame):
            while len(self._order) > self._max_frames:
                idx = self._order[0]
                del self._cache[idx]
                self._order = self._order[1:]

            if index not in self._cache:
                self._order.append(index)

            self._cache[index] = frame

        def is_empty(self):
            return len(self._order) == 0

        def clear(self):
            self._cache = {}
            self._order = {}

        def __getitem__(self, index):
            return self._cache[index]


    @_nglview.register_backend('sire')
    class _SireStructureTrajectory(_nglview.Trajectory, _nglview.Structure):
        def __init__(self, obj=None, map=None):
            if (type(obj) is _SireStructureTrajectory):
                self._traj = obj._traj
                self._map = obj._map
            elif obj is not None:
                from ._trajectory import TrajectoryIterator
                from ..base import PropertyMap

                if map is None:
                    map = PropertyMap()
                else:
                    map = PropertyMap(map)

                self._map = map

                if type(obj) is TrajectoryIterator:
                    self._traj = obj
                else:
                    self._traj = obj.trajectory(self._map)
            else:
                self._traj = None
                self._map = None

            self.ext = 'pdb'
            self.params = {}
            import uuid
            self.id = str(uuid.uuid4())
            self._cache = _FrameCache()

        def __repr__(self):
            return str(self)

        def __str__(self):
            if self._traj is None:
                return "NULL"
            else:
                return str(self._traj)

        def get_structure_string(self):
            from .. import save_to_string
            return "\n".join(save_to_string(self._traj.first(), self.ext))

        @property
        def n_frames(self):
            return max(1, len(self._traj))

        def get_coordinates(self, index):
            if index in self._cache:
                return self._cache[index]

            from ..io import get_coords_array
            from ..units import angstrom

            frame = self._traj[index].current()

            coords = get_coords_array(frame, units=angstrom, map=self._map)

            if coords.size == 0:
                return coords

            if self._cache.is_empty():
                # work out how many frames to cache based on size (no more than 128 MB)
                max_frames = max(1, int(32*1024*1024 / coords.size))
                self._cache = _FrameCache(max_frames=max_frames)

            self._cache.save(index, coords)

            return coords

    def view(obj, map=None):
        struc_traj = _SireStructureTrajectory(obj, map=map)
        view = _nglview.NGLWidget(struc_traj)
        return view

else:
    def view(obj):
        raise ImportError(
            "You need to install nglview to be able to view "
            "molecules. Do this by typing, e.g. "
            "'mamba install nglview'")
