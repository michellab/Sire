
__all__ = ["Cursor"]


class _CursorData:
    """This is the shared data class that holds all of the data
       for a Cursor. This is held by all Cursors that are derived
       from the same object, meaning that multiple Cursors
       can share the same molecule editor. Note that the
       Cursor is not thread-safe (unlike the underlying
       system used by Sire)
    """
    def __init__(self, molecule = None):
        if molecule is None:
            self.molecule = None
        else:
            self.molecule = molecule.molecule().edit()

    def update(self, view):
        try:
            return self.molecule[view.index()]
        except Exception:
            return self.molecule


class Cursor:
    """This class provides a cursor that can be used to navigate through
       and edit the properties of Molecules. This makes the whole
       getting and setting of properties more pythonic in writing
       style, while also saving some typing.
    """
    def __init__(self, molecule = None, bond = None,
                 connectivity_property: str="connectivity"):
        """Construct the Cursor to explore and edit the
           properties of the passed MoleculeView.

           Note that you normally don't call this yourself.
           Instead, you would create a Cursor by calling
           the `.cursor()` function on the molecule view
           itself.

           Examples:
               >>> cursor = mol.cursor()
               >>> cursor["cat"] = "meow"
               >>> mol = cursor.commit()
        """
        self._d = _CursorData(molecule)
        self._view = self._d.update(molecule)
        self._connectivity = None
        self._bond = bond
        self._connectivity_property = connectivity_property

    def _update(self):
        if self._bond is not None:
            try:
                self._connectivity = self._d.molecule.getProperty(
                                        self._connectivity_property).edit()
            except Exception:
                # the molecule doesn't have a connectivity. Create one for it
                from Sire.Mol import CovalentBondHunter
                hunter = CovalentBondHunter()

                try:
                    connectivity = hunter(self._d.molecule)
                    self._d.molecule.setProperty(self._connectivity_property,
                                                 connectivity)
                except Exception:
                    pass

        self._view = self._d.update(self._view)

    def __str__(self):
        return f"Cursor({self.type()}:{self.ID()})"

    def __repr__(self):
        return self.__str__()

    def __delitem__(self, key):
        self._update()

        if self._bond is None:
            self._d.molecule.removeProperty(key)
        else:
            self._connectivity.removeProperty(self._bond, key)
            self._d.molecule.setProperty(self._connectivity_property,
                                         self._connectivity.commit())

        self._update()

    def __contains__(self, key):
        self._update()

        if self._bond is None:
            return self._view.hasProperty(key)
        else:
            return self._connectivity.hasProperty(self._bond, key)

    def __getitem__(self, key):
        self._update()

        if self._bond is None:
            return self._view.property(key)
        else:
            return self._connectivity.property(self._bond, key)

    def __setitem__(self, key, value):
        self._update()

        if self._bond is None:
            self._view.setProperty(key, value)
            self._d.molecule = self._view.molecule()
        else:
            self._connectivity.setProperty(self._bond, key, value)
            self._molecule.setProperty(self._connectivity_property,
                                       self._connectivity.commit())

        self._update()

    def set(self, values):
        """Set all of the properties from the passed dictionary of values"""
        self._update()

        if self._bond is None:
            for key in values.keys():
                self._view.setProperty(key, values[key])

            self._d.molecule = self._view.molecule()
        else:
            for key in values.keys():
                self._connectivity.setProperty(self._bond, key, values[key])

                self._molecule.setProperty(self._connectivity_property,
                                           self._connectivity.commit())

        self._update()

    def atoms(self, id=None):
        """Return cursors for all of atoms in this view,
           of, if 'id' is supplied, the atoms in this view
           that match 'id'
        """
        self._update()

        cursors = []

        if id is None:
            atoms = self._view.atoms()
        else:
            atoms = self._view.atoms(id)

        for atom in atoms:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule.atom(atom.index())
            cursors.append(c)

        return cursors

    def residues(self, id=None):
        """Return cursors for all of residues in this view,
           of, if 'id' is supplied, the residues in this view
           that match 'id'
        """
        self._update()

        cursors = []

        if id is None:
            residues = self._view.residues()
        else:
            residues = self._view.residues(id)

        for residue in residues:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule.residue(residue.index())
            cursors.append(c)

        return cursors

    def chains(self, id=None):
        """Return cursors for all of chains in this view,
           of, if 'id' is supplied, the chains in this view
           that match 'id'
        """
        self._update()

        cursors = []

        if id is None:
            chains = self._view.chains()
        else:
            chains = self._view.chains(id)

        for chain in chains:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule.chain(chain.index())
            cursors.append(c)

        return cursors

    def segments(self, id=None):
        """Return cursors for all of segments in this view,
           of, if 'id' is supplied, the segments in this view
           that match 'id'
        """
        self._update()

        cursors = []

        if id is None:
            segments = self._view.segments()
        else:
            segments = self._view.segments(id)

        for segment in segments:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule.segment(segment.index())
            cursors.append(c)

        return cursors

    def atom(self, i):
        """Return the atom in the molecule that matches the passed ID"""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._view.atom(i)
        c._update()

        return c

    def residue(self, i=None):
        """Return the atom in the molecule that matches the passed ID.
           If 'i' is None, then this returns the residue that contains
           this atom (if this is an atom)
        """
        if i is None:
            try:
                c = Cursor()
                c._d = self._d
                c._view = self._view.residue()
                c._update()
                return c
            except Exception:
                raise TypeError(
                    f"There is no residue that contains {self.type()}:{self.ID()}"
                )

        self._update()

        c = Cursor()
        c._d = self._d
        c._view = c._view.residue(i)
        c._update()

        return c

    def chain(self, i=None):
        """Return the chain in the molecule that matches the passed ID.
           If 'i' is None, then this returns the residue that contains
           this atom (if this is an atom)"""
        if i is None:
            try:
                c = Cursor()
                c._d = self._d
                c._view = self._view.chain()
                c._update()
                return c
            except Exception:
                raise TypeError(
                    f"There is no chain that contains {self.type()}:{self.ID()}"
                )

        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._view.chain(i)
        c._update()

        return c

    def segment(self, i=None):
        """Return the segment in the molecule that matches the passed ID.
           If 'i' is None, then this returns the residue that contains
           this atom (if this is an atom)"""
        if i is None:
            try:
                c = Cursor()
                c._d = self._d
                c._view = self._view.segment()
                c._update()
                return c
            except Exception:
                raise TypeError(
                    f"There is no segment that contains {self.type()}:{self.ID()}"
                )

        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._view.segment(i)
        c._update()

        return c

    def molecule(self):
        """Return the molecule"""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule
        return c

    def bond(self, bond, connectivity_property="connectivity"):
        """Return the Cursor for the specified bond. This will
           use the specified connectivity property to locate
           the connectivity that defines the bonds in this molecule
        """
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule
        c._connectivity_property = connectivity_property
        c._bond = bond
        c._update()

        return c

    def parent(self):
        """Return the cursor of the parent object (e.g. parent residue
           of the atom, parent chain of the residue, parent molecule
           of the bond etc. This will return the Cursor for the whole
           molecule if there isn't a suitable parent
        """
        t = self.type()

        c = Cursor()
        c._d = self._d

        c._view = c._d.molecule

        # The parent of an Atom is a Residue (if it is in one), and
        # the parent of a Residue is a Chain (if it is in one).
        # If this fails, then the parent is the Molecule
        try:
            if t == "atom":
                c._view = self._view.residue()
            elif t == "residue":
                c._view = self._view.chain()
        except Exception:
            pass

        c._update()

        return c

    def name(self):
        """Return the name of the current view"""
        self._update()

        if self._bond is not None:
            raise TypeError("A bond does not have a name!")

        return self._view.name()

    def number(self):
        """Return the number of the current view"""
        self._update()

        if self._bond is not None:
            raise TypeError("A bond does not have a number!")

        try:
            return self._view.number()
        except Exception:
            raise TypeError(f"A {self._view.what()} does not have a number!")

    def index(self):
        """Return the index of the current view (e.g. AtomIdx, ResIdx etc)"""
        self._update()

        if self._bond is not None:
            raise TypeError("A bond does not have an index!")

        return self._view.index()

    def ID(self):
        """Return the ID of this view (e.g. AtomIdx, MolNum, BondID)"""
        self._update()

        if self._bond is not None:
            return self._bond

        try:
            return self._view.index()
        except Exception:
            # objects without an index (e.g. molecule) use a number for their ID
            return self._view.number()

    def type(self):
        """Return the type of this Cursor (e.g. 'atom', 'bond',
           'residue', 'chain', 'segment' or 'molecule')
        """
        if self.isBond():
            return "bond"

        w = self._view.what()

        if w.find("Atom") != -1:
            return "atom"
        elif w.find("Res") != -1:
            return "residue"
        elif w.find("Chain") != -1:
            return "chain"
        elif w.find("Seg") != -1:
            return "segment"
        elif w.find("Mol") != -1:
            return "molecule"
        else:
            raise TypeError(f"Cannot identify cursor type {w}")

    def isMolecule(self):
        """Return whether this is pointing to a Molecule"""
        return self.type() == "molecule"

    def isBond(self):
        """Return whether this is pointing to a Bond"""
        self._update()

        return self._bond is not None

    def isAtom(self):
        """Return whether this is pointing to an Atom"""
        return self.type() == "atom"

    def isResidue(self):
        """Return whether this is pointing to a Residue"""
        return self.type() == "residue"

    def isChain(self):
        """Return whether this is pointing to a Chain"""
        return self.type() == "chain"

    def isSegment(self):
        """Return whether this is pointing to a Segment"""
        return self.type() == "segment"

    def commit(self):
        """Commit all of the changes and return the newly
           edited molecule (or MoleculeView)
        """
        self._update()

        if self._connectivity is not None:
            self._molecule.setProperty(self._connectivity_property,
                                       self._connectivity.commit())
            self._update()

        mol = self._d.molecule.commit()

        try:
            return mol[self._view.index()]
        except Exception:
            return mol

    def keys(self):
        self._update()

        if self._bond is None:
            return self._view.propertyKeys()
        else:
            return self._connectivity.propertyKeys(self._bond)

    def values(self):
        self._update()

        try:
            if self._bond is None:
                return self._view.propertyValues()
            else:
                return self._connectivity.propertyValues(self._bond)
        except Exception:
            vals = []

            for key in self.keys():
                vals.append(self.__getitem__(key))

            return vals

    def items(self):
        self._update()

        if self._bond is None:
            keys = self._view.propertyKeys()
            items = []

            for key in keys:
                items.append((key, self._view.property(key)))
        else:
            keys = self._connectivity.propertyKeys(self._bond)
            items = []

            for key in keys:
                items.append((key, self._connectivity.property(self._bond,
                                                               key)))

        return items

    def properties(self):
        """Return the Sire.Base.Properties object for the properties
           of the current view
        """
        from Sire.Base import Properties
        p = Properties()

        for key in self.keys():
            p[key] = self.__getitem__(key)

        return p
