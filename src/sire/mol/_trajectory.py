
__all__ = ["TrajectoryIterator"]


class TrajectoryIterator:
    def __init__(self, view=None):
        if view is not None:
            self._view = view
            self._values = range(0, max(1, self._view.num_frames()))
            self._iter = None
        else:
            self._view = None
            self._values = []
            self._iter = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._view is None or self._values is None:
            raise StopIteration()

        if self._iter is None:
            self._iter = self._values.__iter__()

        frame = self._iter.__next__()

        ret = self._view.clone()

        ret.load_frame(frame)

        return ret

    def __getitem__(self, val):
        it = TrajectoryIterator()
        it._view = self._view

        if type(val) is int:
            it._values = [self._values[val]]
        elif type(val) is slice:
            it._values = self._values[val]
        else:
            it._values = [self._values[v] for v in val]

        return it

    def _assert_pandas_loaded(self):
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "You must install pandas to be able to use the measures function! "
                "Install pandas using, e.g. `conda install pandas` or "
                "`mamba install pandas`."
            )

    def measures(self):
        self._assert_pandas_loaded()

        import pandas as pd
        from .._colname import colname

        if self._view is None:
            return pd.DataFrame()

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

        rows = []
        columns = []

        if uses_measures:
            for view in self._view:
                columns.append(colname(view))

            for frame in self.__iter__():
                rows.append([measure.value() for measure in frame.measures()])
        else:
            columns.append(colname(view))

            for frame in self.__iter__():
                rows.append([frame.measure().value()])

        return pd.DataFrame(rows, columns=columns)
