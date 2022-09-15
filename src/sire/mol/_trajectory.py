
__all__ = ["TrajectoryIterator"]


class TrajectoryIterator:
    def __init__(self, view=None, map=None):
        if view is not None:
            from ..legacy.Base import PropertyMap

            if map is None:
                self._map = PropertyMap()
            else:
                self._map = PropertyMap(map)

            self._view = view
            self._values = range(0, max(1, self._view.num_frames(self._map)))
            self._iter = None
            self._frame = None
        else:
            self._view = None
            self._values = []
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
            ret.frame_time = lambda: 0

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
                task = progress.add_task("Looping through frames", total=nframes)

                for idx, frame in enumerate(self.__iter__()):
                    for i, measure in enumerate(frame.measures(map=self._map)):
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
        else:
            colnames.append(colname(view))
            column = np.zeros(nframes, dtype=float)

            with Console.progress() as progress:
                task = progress.add_task("Looping through frames", total=nframes)

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

    def measures(self, func=None, to_pandas=False):
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
            list: A list of the results, with one result per frame in the trajectory
        """
        result = []

        from ..utils import Console

        nframes = len(self)

        if str(func) == func:
            # we calling a named function
            with Console.progress() as progress:
                task = progress.add_task("Looping through frames", total=nframes)

                for i in range(0, nframes):
                    obj = self.__getitem__(i).current()
                    result.append(getattr(obj, func)(*args, **kwargs))
                    progress.update(task, completed=i+1)

        else:
            # we have been passed the function to call
            with Console.progress() as progress:
                task = progress.add_task("Looping through frames", total=nframes)

                for i in range(0, nframes):
                    obj = self.__getitem__(i).current()
                    result.append(func(obj, *args, **kwargs))
                    progress.update(task, completed=i+1)

        return result

    def view(self, *args, **kwargs):
        from ._view import view
        return view(self, *args, **kwargs)
