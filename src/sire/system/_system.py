
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

    def num_frames(self):
        return self._system.num_frames()

    def load_frame(self, i):
        self._system.load_frame(i)
        self._molecules = None

    def save_frame(self, i=None):
        if i is None:
            self._system.save_frame()
        else:
            self._system.save_frame(i)

        self._molecules = None

    def delete_frame(self, i):
        self._system.delete_frame(i)
        self._molecules = None

    def to_molecule_group(self):
        return self.molecules().to_molecule_group()

    def molecules(self, *args, **kwargs):
        if self._molecules is not None:
            return self._molecules.molecules(*args, **kwargs)

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

        return self.molecules(*args, **kwargs)

    def segments(self, *args, **kwargs):
        return self.molecules().segments(*args, **kwargs)

    def chains(self, *args, **kwargs):
        return self.molecules().chains(*args, **kwargs)

    def residues(self, *args, **kwargs):
        return self.molecules().residues(*args, **kwargs)

    def atoms(self, *args, **kwargs):
        return self.molecules().atoms(*args, **kwargs)

    def bonds(self, *args, **kwargs):
        return self.molecules().bonds(*args, **kwargs)

    def angles(self, *args, **kwargs):
        return self.molecules().angles(*args, **kwargs)

    def dihedrals(self, *args, **kwargs):
        return self.molecules().dihedrals(*args, **kwargs)

    def impropers(self, *args, **kwargs):
        return self.molecules().impropers(*args, **kwargs)

    def molecule(self, *args, **kwargs):
        return self.molecules().molecule(*args, **kwargs)

    def segment(self, *args, **kwargs):
        return self.molecules().segment(*args, **kwargs)

    def chain(self, *args, **kwargs):
        return self.molecules().chain(*args, **kwargs)

    def residue(self, *args, **kwargs):
        return self.molecules().residue(*args, **kwargs)

    def atom(self, *args, **kwargs):
        return self.molecules().atom(*args, **kwargs)

    def bond(self, *args, **kwargs):
        return self.molecules().bond(*args, **kwargs)

    def angle(self, *args, **kwargs):
        return self.molecules().angle(*args, **kwargs)

    def dihedral(self, *args, **kwargs):
        return self.molecules().dihedral(*args, **kwargs)

    def improper(self, *args, **kwargs):
        return self.molecules().improper(*args, **kwargs)

    def trajectory(self):
        return self.molecules().trajectory()

    def energy(self, *args, **kwargs):
        return self.molecules().energy(*args, **kwargs)

    def energies(self, *args, **kwargs):
        return self.molecules().energies(*args, **kwargs)

    def charge(self, *args, **kwargs):
        return self.molecules().charge(*args, **kwargs)

    def mass(self, *args, **kwargs):
        return self.molecules().mass(*args, **kwargs)

    def coordinates(self, *args, **kwargs):
        return self.molecules().coordinates(*args, **kwargs)

    def evaluate(self, *args, **kwargs):
        return self.molecules().evaluate(*args, **kwargs)

    def cursor(self):
        return self.molecules().cursor()

    def view(self, *args, **kwargs):
        return self.molecules().view(*args, **kwargs)

    def update(self, value):
        self._molecules = None
        self._system.update(value)
