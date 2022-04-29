
__all__ = ["System"]


class System:
    def __init__(self, system=None):
        from Sire.System import System as _System

        if system is None:
            self._system = _System()
        else:
            if _System not in type(system).mro():
                print(type(system).mro())
                raise TypeError(
                    "You can only construct from a Sire.System, "
                    f"not a {type(system)}")

            self._system = system

    def __str__(self):
        return str(self._system)

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, key):
        from Sire.Mol import SelectorMol
        return SelectorMol(self._system)[key]

    def molecules(self, key=None):
        from Sire.Mol import SelectorMol

        if key is None:
            return SelectorMol(self._system)
        else:
            return SelectorMol(self._system, key)

