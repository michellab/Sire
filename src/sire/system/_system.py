
__all__ = ["System"]


class System:
    def __init__(self, system=None):
        from ..legacy.System import System as _System

        if system is None:
            self._system = _System()
        else:
            if _System not in type(system).mro():
                raise TypeError(
                    "You can only construct from a sire.legacy.System.System, "
                    f"not a {type(system)}")

            self._system = system

        self._molecules = None

    def __copy__(self, other):
        if type(other) is not System:
            raise TypeError(f"Cannot copy a {type(other)} to a System")

        self._system = other._system.clone()
        self._molecules = None

    def __deepcopy__(self, other):
        self.__copy__(other)

    def __str__(self):
        return str(self._system)

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, key):
        return self.molecules()[key]

    def clone(self):
        s = System()
        s._system = self._system.clone()
        return s

    def count(self):
        return self.__len__()

    def size(self):
        return self.__len__()

    def __len__(self):
        return len(self.molecules())

    def num_atoms(self):
        return self._system.num_atoms()

    def num_residues(self):
        return self._system.num_residues()

    def num_chains(self):
        return self._system.num_chains()

    def num_segments(self):
        return self._system.num_segments()

    def num_molecules(self):
        return self._system.num_molecules()

    def names(self):
        return self.molecules().names()

    def numbers(self):
        return self.molecules().numbers()

    def molecules(self, key=None):
        if self._molecules is not None:
            if key is None:
                return self._molecules
            else:
                return self._molecules.molecules(key)

        import sire.mol
        self._molecules = sire.mol.SelectorMol(self._system)

        if self._molecules.num_atoms() != self._system.num_atoms():
            # oh dear - this is an edge case where the System does
            # not contain complete molecules. We need to extract
            # the molecules and re-add them
            raise NotImplementedError(
                "sire.system.System does not yet support Systems that hold "
                "partial molecule!. Let us know that you have hit this "
                "bug and we will add support."
            )

        return self.molecules(key)

    def segments(self, key=None):
        if key is None:
            return self.molecules().segments()
        else:
            return self.molecules().segments(key)

    def chains(self, key=None):
        if key is None:
            return self.molecules().chains()
        else:
            return self.molecules().chains(key)

    def residues(self, key=None):
        if key is None:
            return self.molecules().residues()
        else:
            return self.molecules().residues(key)

    def atoms(self, key=None):
        if key is None:
            return self.molecules().atoms()
        else:
            return self.molecules().atoms(key)

    def molecule(self, key):
        return self.molecules().molecule(key)

    def segment(self, key):
        return self.molecules().segment(key)

    def chain(self, key):
        return self.molecules().chain(key)

    def residue(self, key):
        return self.molecules().residue(key)

    def atom(self, key):
        return self.molecules().atom(key)

    def update(self, value):
        self._molecules = None
        self._system.update(value)
