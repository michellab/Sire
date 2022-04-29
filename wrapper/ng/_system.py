
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

        self._molecules = None

    def __str__(self):
        return str(self._system)

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, key):
        return self.molecules()[key]

    def count(self):
        return self.molecules().count()

    def __len__(self):
        return len(self.molecules())

    def num_atoms(self):
        return self._system.nAtoms()

    def num_residues(self):
        return self._system.nResidues()

    def num_chains(self):
        return self._system.nChains()

    def num_segments(self):
        return self._system.nSegments()

    def num_molecules(self):
        return self._system.nMolecules()

    def nAtoms(self):
        return self.num_atoms()

    def nResidues(self):
        return self.num_residues()

    def nChains(self):
        return self.num_chains()

    def nSegments(self):
        return self.num_segments()

    def nMolecules(self):
        return self.num_molecules()

    def molecules(self, key=None):
        if self._molecules is not None:
            if key is None:
                return self._molecules
            else:
                return self._molecules.molecules(key)

        from Sire.Mol import SelectorMol
        self._molecules = SelectorMol(self._system)

        if self._molecules.nAtoms() != self._system.nAtoms():
            # oh dear - this is an edge case where the System does
            # not contain complete molecules. We need to extract
            # the molecules and re-add them
            raise IncompatibleError(
                "Sire.ng.System does not yet support Systems that hold "
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
