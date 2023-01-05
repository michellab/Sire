__all__ = ["TrajectoryIterator"]


class TrajectoryIterator:
    """
    An iterator that can be used to control which frames of a trajectory
    are accessed or processed.
    """

    def __init__(self, view=None, map=None):
        if view is not None:
            from ..base import create_map

            self._map = create_map(map)
            self._view = view
            self._values = range(0, max(1, self._view.num_frames(self._map)))
            self._times = None
            self._iter = None
            self._frame = None
        else:
            self._view = None
            self._values = []
            self._times = None
            self._iter = None
            self._map = None
            self._frame = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._view is None or self._values is None:
            raise StopIteration()

        if self._iter is None:
            self._iter = self._values.__iter__()

        self._frame = self._iter.__next__()

        return self.current()

    def __len__(self):
        return len(self._values)

    def __getitem__(self, val):
        it = TrajectoryIterator()
        it._view = self._view
        it._map = self._map

        if type(val) is int:
            it._values = [self._values[val]]
        elif type(val) is slice:
            it._values = self._values[val]
        else:
            it._values = [self._values[v] for v in val]

        return it

    def __str__(self):
        if self._view is None:
            return "TrajectoryIterator::null"
        elif len(self._values) <= 1:
            return str(self._view)
        else:
            return f"Trajectory({self._view}, num_frames={len(self._values)})"

    def num_frames(self):
        return len(self._values)

    def current(self):
        """Return the current frame in the trajectory"""
        if self._view is None or self._values is None:
            raise StopIteration()

        if self._frame is None:
            self._frame = self._values[0]

        ret = self._view.clone()

        ret.load_frame(self._frame, map=self._map)

        try:
            mol = ret.molecule()
        except Exception:
            mol = ret[0].molecule()

        time_property = self._map["time"]

        if mol.has_property(time_property):
            time = mol.property(time_property)
            ret.frame_time = lambda: time
        else:
            from ..units import picosecond

            ret.frame_time = lambda: 0 * picosecond

        ret.frame_index = lambda: self._frame

        return ret

    def first(self):
        """Return the first frame in the trajectory"""
        if self._view is None or self._values is None:
            raise StopIteration()

        old_frame = self._frame

        self._frame = self._values[0]

        ret = self.current()

        self._frame = old_frame

        return ret

    def times(self):
        if self._times is not None:
            return self._times

        if self._view is None:
            return {}

        # load the times from the actual underlying trajectory data
        try:
            mol = self._view.molecule()
        except Exception:
            mol = self._view[0].molecule()

        traj = mol.property(self._map["trajectory"])

        self._times = []

        for idx in self._values:
            self._times.append(traj[idx].time())

        return self._times

    def energies(self, obj1=None, forcefield=None, to_pandas=True, map=None):
        if self._view is None:
            return {}

        import numpy as np

        from ..mm import create_forcefield
        from ..legacy.MM import calculate_trajectory_energies
        from .._colname import colname
        from . import _to_molecules
        from ..base import create_map

        map = self._map.merge(create_map(map))

        colnames = []
        forcefields = []

        if obj1 is None:
            for v in self.first():
                colnames.append(colname(v))
                forcefields.append(create_forcefield(v, map=map))
        else:
            if type(obj1) is TrajectoryIterator:
                if obj1.num_frames() != self.num_frames():
                    raise ValueError(
                        "The two trajectories have a different "
                        "number of frames! "
                        f"{self.num_frames()} versus f{obj1.num_frames()}."
                    )

                obj1_mols = _to_molecules(obj1.first())

                for v in self.first():
                    colnames.append(colname(v))
                    forcefields.append(
                        create_forcefield(v, obj1_mols, map=map)
                    )
            else:
                for v in self.first():
                    colnames.append(colname(v))
                    forcefields.append(
                        create_forcefield(v, _to_molecules(obj1), map=map)
                    )

        nframes = len(self)

        times = np.zeros(nframes, dtype=float)
        indexes = np.zeros(nframes, dtype=int)

        t = self.times()

        for i, idx in enumerate(self._values):
            times[i] = t[i].to_default()
            indexes[i] = idx

        time_unit = t[0].get_default().unit_string()
        energy_unit = None

        components = {}

        from ..utils import Console

        import os

        cpu_count = os.cpu_count()

        with Console.progress() as progress:
            task = progress.add_task("Looping through frames", total=nframes)

            num_per_chunk = cpu_count

            i = 0

            import time

            while i < nframes:
                start_time = time.time()
                j = min(i + num_per_chunk, nframes)

                ff_nrgs = calculate_trajectory_energies(
                    forcefields, list(self._values[i:j]), map=map
                )

                if i == 0:
                    for ff_idx in range(0, len(forcefields)):
                        nrg = ff_nrgs[ff_idx][0]

                        if ff_idx == 0:
                            energy_unit = nrg.get_default().unit_string()

                        components[
                            colname(colnames[ff_idx], "total")
                        ] = np.zeros(nframes, dtype=float)

                        for key in nrg.components().keys():
                            components[
                                colname(colnames[ff_idx], key)
                            ] = np.zeros(nframes, dtype=float)

                for idx in range(i, j):
                    for ff_idx in range(0, len(forcefields)):
                        nrg = ff_nrgs[ff_idx][idx - i]
                        components[colname(colnames[ff_idx], "total")][
                            idx
                        ] = nrg.to_default()

                        for key, value in nrg.components().items():
                            try:
                                components[colname(colnames[ff_idx], key)][
                                    idx
                                ] = nrg[key].to_default()
                            except KeyError:
                                k = colname(colnames[ff_idx], key)
                                components[k] = np.zeros(nframes, dtype=float)
                                components[k][idx] = nrg[key].to_default()

                    progress.update(task, completed=idx)

                delta = time.time() - start_time

                if delta > 0.8:
                    # we want about 0.8 seconds between updates
                    num_per_chunk = int(num_per_chunk / 2)
                    if num_per_chunk < cpu_count:
                        num_per_chunk = cpu_count
                elif delta < 0.25:
                    num_per_chunk = num_per_chunk + int(0.5 * num_per_chunk)

                i = j

        data = {}

        data["frame"] = indexes
        data["time"] = times

        colnames = list(components.keys())
        colnames.sort()

        for name in colnames:
            data[name] = components[name]

        if to_pandas:
            import pandas as pd

            df = pd.DataFrame(data)
            df.set_index("frame")

            df.time_unit = lambda: time_unit
            df.energy_unit = lambda: energy_unit

            def pretty_plot(x="time", y=None):
                if y is None:
                    y = colnames

                if x == "time":
                    xlabel = f"Time / {df.time_unit()}"
                elif x == "frame":
                    xlabel = "Frame"
                else:
                    xlabel = x

                ax = df.plot(x=x, y=y)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(f"Energy / {df.energy_unit()}")
                ax.legend(bbox_to_anchor=(1.05, 1.0))

            df.pretty_plot = pretty_plot

            return df

        return data

    def energy(self, obj1=None, to_pandas=True, map=None):
        if self._view is None:
            return {}

        import numpy as np

        from ..mm import create_forcefield
        from ..legacy.MM import calculate_trajectory_energy
        from ..base import create_map

        map = self._map.merge(create_map(map))

        if obj1 is None:
            ff = create_forcefield(self.first(), map=map)
        else:
            if type(obj1) is TrajectoryIterator:
                if obj1.num_frames() != self.num_frames():
                    raise ValueError(
                        "The two trajectories have a different "
                        "number of frames! "
                        f"{self.num_frames()} versus f{obj1.num_frames()}."
                    )

                ff = create_forcefield(self.first(), obj1.first(), map=map)
            else:
                ff = create_forcefield(self.first(), obj1, map=map)

        nframes = len(self)

        times = np.zeros(nframes, dtype=float)
        indexes = np.zeros(nframes, dtype=int)

        t = self.times()

        for i, idx in enumerate(self._values):
            times[i] = t[i].to_default()
            indexes[i] = idx

        time_unit = t[0].get_default().unit_string()
        energy_unit = None

        components = {}

        from ..utils import Console

        import os

        cpu_count = os.cpu_count()

        with Console.progress() as progress:
            task = progress.add_task("Looping through frames", total=nframes)

            num_per_chunk = cpu_count

            i = 0

            import time

            while i < nframes:
                start_time = time.time()
                j = min(i + num_per_chunk, nframes)

                nrgs = calculate_trajectory_energy(
                    ff, list(self._values[i:j]), map
                )

                if i == 0:
                    nrg = nrgs[0]
                    energy_unit = nrg.get_default().unit_string()
                    components["total"] = np.zeros(nframes, dtype=float)
                    for key in nrg.components().keys():
                        components[key] = np.zeros(nframes, dtype=float)

                for idx in range(i, j):
                    nrg = nrgs[idx - i]
                    components["total"][idx] = nrg.to_default()

                    for key, value in nrg.components().items():
                        try:
                            components[key][idx] = nrg[key].to_default()
                        except KeyError:
                            components[key] = np.zeros(nframes, dtype=float)
                            components[key][idx] = nrg[key].to_default()

                    progress.update(task, completed=idx)

                delta = time.time() - start_time

                if delta > 0.8:
                    # we want about 0.8 seconds between updates
                    num_per_chunk = int(num_per_chunk / 2)
                    if num_per_chunk < cpu_count:
                        num_per_chunk = cpu_count
                elif delta < 0.25:
                    num_per_chunk = num_per_chunk + int(0.5 * num_per_chunk)

                i = j

        data = {}

        data["frame"] = indexes
        data["time"] = times

        colnames = list(components.keys())
        colnames.sort()

        for colname in colnames:
            data[colname] = components[colname]

        if to_pandas:
            import pandas as pd

            df = pd.DataFrame(data)
            df.set_index("frame")

            df.time_unit = lambda: time_unit
            df.energy_unit = lambda: energy_unit

            def pretty_plot(x="time", y=None):
                if y is None:
                    y = colnames

                if x == "time":
                    xlabel = f"Time / {df.time_unit()}"
                elif x == "frame":
                    xlabel = "Frame"
                else:
                    xlabel = x

                ax = df.plot(x=x, y=y)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(f"Energy / {df.energy_unit()}")
                ax.legend(bbox_to_anchor=(1.05, 1.0))

            df.pretty_plot = pretty_plot

            return df

        return data

    def _simple_measures(self, to_pandas):
        from .._colname import colname

        if self._view is None:
            return {}

        uses_measures = None

        if hasattr(self._view, "measures"):
            uses_measures = True
        elif hasattr(self._view, "measure"):
            uses_measures = False
        else:
            raise AttributeError(
                f"This view ({self._view}) does not have a `.measure()` "
                "or `.measures()` function, so cannot be measured."
            )

        import numpy as np

        colnames = []
        columns = []

        nframes = len(self)

        times = np.zeros(nframes, dtype=float)
        indexes = np.zeros(nframes, dtype=int)

        time_unit = None
        measure_unit = None

        from ..utils import Console

        if uses_measures:
            for view in self._view:
                colnames.append(colname(view))
                columns.append(np.zeros(nframes, dtype=float))

            with Console.progress() as progress:
                task = progress.add_task(
                    "Looping through frames", total=nframes
                )

                for idx, frame in enumerate(self.__iter__()):
                    for i, measure in enumerate(frame.measures(map=self._map)):
                        columns[i][idx] = measure.to_default()
                        times[idx] = frame.frame_time().to_default()

                        if measure_unit is None:
                            if not measure.is_zero():
                                measure_unit = (
                                    measure.get_default().unit_string()
                                )

                        if time_unit is None:
                            time = frame.frame_time()
                            if not time.is_zero():
                                time_unit = time.get_default().unit_string()

                    indexes[idx] = frame.frame_index()
                    progress.update(task, completed=idx)
        else:
            colnames.append(colname(view))
            column = np.zeros(nframes, dtype=float)

            with Console.progress() as progress:
                task = progress.add_task(
                    "Looping through frames", total=nframes
                )

                for idx, frame in enumerate(self.__iter__()):
                    measure = frame.measure(map=self._map)
                    column[idx] = measure.to_default()
                    times[idx] = frame.frame_time().to_default()

                    if measure_unit is None:
                        if not measure.is_zero():
                            measure_unit = measure.get_default().unit_string()

                    if time_unit is None:
                        time = frame.frame_time()
                        if not time.is_zero():
                            time_unit = time.get_default().unit_string()

                    indexes[idx] = frame.frame_index()
                    progress.update(task, completed=idx)

            columns = [column]

        data = {}

        data["frame"] = indexes
        data["time"] = times

        for i in range(0, len(colnames)):
            data[colnames[i]] = columns[i]

        if to_pandas:
            import pandas as pd

            df = pd.DataFrame(data)
            df.set_index("frame")

            df.time_unit = lambda: time_unit
            df.measure_unit = lambda: measure_unit

            def pretty_plot(x="time", y=None):
                if y is None:
                    y = colnames

                if x == "time":
                    xlabel = f"Time / {df.time_unit()}"
                elif x == "frame":
                    xlabel = "Frame"
                else:
                    xlabel = x

                ax = df.plot(x=x, y=y)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(f"Size / {df.measure_unit()}")
                ax.legend(bbox_to_anchor=(1.05, 1.0))

            df.pretty_plot = pretty_plot

            return df

        return data

    def _custom_measures(self, func, to_pandas):
        if self._view is None:
            return {}

        if not type(func) is dict:
            func = {"custom": func}

        import numpy as np

        colnames = []
        columns = []

        nframes = len(self)

        times = np.zeros(nframes, dtype=float)
        indexes = np.zeros(nframes, dtype=int)

        time_unit = None
        measure_unit = None

        for key in func.keys():
            colnames.append(key)
            columns.append(np.zeros(nframes, dtype=float))

        from ..utils import Console

        with Console.progress() as progress:
            task = progress.add_task("Looping through frames", total=nframes)

            for idx, frame in enumerate(self.__iter__()):
                for i, f in enumerate(func.values()):
                    measure = f(frame)
                    columns[i][idx] = measure.to_default()
                    times[idx] = frame.frame_time().to_default()

                    if measure_unit is None:
                        if not measure.is_zero():
                            measure_unit = measure.get_default().unit_string()

                    if time_unit is None:
                        time = frame.frame_time()
                        if not time.is_zero():
                            time_unit = time.get_default().unit_string()

                indexes[idx] = frame.frame_index()
                progress.update(task, completed=idx)

        data = {}

        data["frame"] = indexes
        data["time"] = times

        for i in range(0, len(colnames)):
            data[colnames[i]] = columns[i]

        if to_pandas:
            import pandas as pd

            df = pd.DataFrame(data)
            df.set_index("frame")

            df.time_unit = lambda: time_unit
            df.measure_unit = lambda: measure_unit

            def pretty_plot(x="time", y=None):
                if y is None:
                    y = colnames

                if x == "time":
                    xlabel = f"Time / {df.time_unit()}"
                elif x == "frame":
                    xlabel = "Frame"
                else:
                    xlabel = x

                ax = df.plot(x=x, y=y)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(f"Size / {df.measure_unit()}")
                ax.legend(bbox_to_anchor=(1.05, 1.0))

            df.pretty_plot = pretty_plot

            return df

        return data

    def measures(self, func=None, to_pandas=True):
        if func is None:
            return self._simple_measures(to_pandas=to_pandas)
        else:
            return self._custom_measures(func=func, to_pandas=to_pandas)

    def apply(self, func, *args, **kwargs):
        """
        Call the passed function on all frames of the trajectory,
        appending the result to a list of results, which
        is returned.

        The function can be either;

        1. a string containing the name of the function to call, or
        2. an actual function (either a normal function or a lambda expression)

        You can optionally pass in positional and keyword arguments
        here that will be passed to the function.

        Args:
            func (str or function): The function to be called, or the name
                                    of the function to be called.

        Returns:
            list: A list of the results, with one result per
                  frame in the trajectory
        """
        result = []

        from ..utils import Console

        nframes = len(self)

        if str(func) == func:
            # we calling a named function
            with Console.progress() as progress:
                task = progress.add_task(
                    "Looping through frames", total=nframes
                )

                for i in range(0, nframes):
                    obj = self.__getitem__(i).current()
                    result.append(getattr(obj, func)(*args, **kwargs))
                    progress.update(task, completed=i + 1)

        else:
            # we have been passed the function to call
            with Console.progress() as progress:
                task = progress.add_task(
                    "Looping through frames", total=nframes
                )

                for i in range(0, nframes):
                    obj = self.__getitem__(i).current()
                    result.append(func(obj, *args, **kwargs))
                    progress.update(task, completed=i + 1)

        return result

    def view(self, *args, **kwargs):
        from ._view import view

        return view(self, *args, **kwargs)
