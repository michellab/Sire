
__all__ = ["Cursor"]

class Cursor:
    """This class provides a cursor that can be used to navigate through
       and edit the properties of Molecules. This makes the whole
       getting and setting of properties more pythonic in writing
       style, while also saving some typing.
    """
    def __init__(self, molecule = None):
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
        if molecule is None:
            self._molecule = None
            self._view = None
        else:
            self._molecule = molecule.molecule().edit()

            try:
                self._view = self._molecule[molecule.index()]
            except Exception:
                self._view = self._molecule

        self._connectivity = None
        self._bond = None
        self._connectivity_property = "connectivity"

    def __str__(self):
        if self._bond is None:
            return f"Cursor({self._view})"
        else:
            return f"Cursor(bond:{self._bond})"

    def __repr__(self):
        return self.__str__()

    def __delitem__(self, key):
        if self._bond is None:
            print("HERE")
            self._molecule.removeProperty(key)
            try:
                print("HERE2")
                self._view = self._molecule[self._view.index()]
            except Exception:
                print("HERE3")
                self._view = self._molecule
        else:
            self._connectivity.removeProperty(bond, key)
            self._molecule.setProperty(self._connectivity_property,
                                       self._connectivity.commit())

    def __getitem__(self, key):
        if self._bond is None:
            return self._view.property(key)
        else:
            return self._connectivity.property(bond, key)

    def __setitem__(self, key, value):
        if self._bond is None:
            self._view.setProperty(key, value)
            self._molecule = self._view.molecule()
        else:
            self._connectivity.setProperty(bond, key, value)
            self._molecule.setProperty(self._connectivity_property,
                                       self._connectivity.commit())

    def set(self, values):
        """Set all of the properties from the passed dictionary of values"""
        if self._bond is None:
            for key in values.keys():
                self._view.setProperty(key, values[key])
        else:
            for key in values.keys():
                self._connectivity.setProperty(bond, key, values[key])

        self._molecule.setProperty(self._connectivity_property,
                                   self._connectivity.commit())

    def atom(self, i):
        """Return the atom in the molecule that matches the passed ID"""
        c = Cursor()
        c._molecule = self._molecule
        c._view = self._molecule.atom(i)
        return c

    def residue(self, i):
        """Return the atom in the molecule that matches the passed ID"""
        c = Cursor()
        c._molecule = self._molecule
        c._view = self._molecule.residue(i)
        return c

    def chain(self, i):
        """Return the chain in the molecule that matches the passed ID"""
        c = Cursor()
        c._molecule = self._molecule
        c._view = self._molecule.chain(i)
        return c

    def segment(self, i):
        """Return the segment in the molecule that matches the passed ID"""
        c = Cursor()
        c._molecule = self._molecule
        c._view = self._molecule.segment(i)
        return c

    def molecule(self):
        """Return the molecule"""
        c = Cursor()
        c._molecule = self._molecule
        c._view = self._molecule
        return c

    def bond(self, bond, connectivity_property="connectivity"):
        """Return the Cursor for the specified bond. This will
           use the specified connectivity property to locate
           the connectivity that defines the bonds in this molecule
        """
        c = Cursor()
        c._molecule = self._molecule
        c._view = self._molecule
        c._connectivity = c._molecule.property(connectivity_property).edit()
        c._connectivity_property = connectivity_property
        return c

    def parent(self):
        """Return the cursor of the parent object (e.g. parent residue
           of the atom, parent chain of the residue, parent molecule
           of the bond etc. This will return the Cursor for the whole
           molecule if there isn't a suitable parent
        """
        c = Cursor()
        c._molecule = self._molecule

        try:
            c._view = self._view.parent()
        except Exception:
            c._view = self._molecule

        c._connectivity = None
        c._connectivity_property = None

        return c

    def next(self):
        """Return the cursor to the next logical view (or bond)
           This will go to the next AtomIdx, or next ResIdx,
           or the next bond. This will raise an exception
           (StopIteration) if there is no next view.
        """
        if self._connectivity is None:
            try:
                idx = self._view.index()
                idx += 1

                c = Cursor()
                c._molecule = self._molecule
                c._view = self._molecule[idx]
                return c
            except Exception:
                raise StopIteration()
        else:
            try:
                bonds = self._connectivity.bonds()
                # find index of current bond...
                raise ValueError()

                c = Cursor()
                c._molecule = self._molecule
                c._view = self._molecule
                c._connectivity = self._connectivity
                c._connectivity_property = self._connectivity_property
                c._bond = next_bond
                return c
            except Exception:
                raise StopIteration()

    def prev(self):
        """Return the cursor to the previous logical view (or bond)
           This will go to the previous AtomIdx, or previous ResIdx,
           or the previous bond. This will raise an exception
           (StopIteration) if there is no previous view.
        """
        if self._connectivity is None:
            try:
                idx = self._view.index()
                idx -= 1

                if idx.value() < 0:
                    raise StopIteration()

                c = Cursor()
                c._molecule = self._molecule
                c._view = self._molecule[idx]
                return c
            except Exception:
                raise StopIteration()
        else:
            try:
                bonds = self._connectivity.bonds()
                # find index of current bond...
                raise ValueError()

                c = Cursor()
                c._molecule = self._molecule
                c._view = self._molecule
                c._connectivity = self._connectivity
                c._connectivity_property = self._connectivity_property
                c._bond = next_bond
                return c
            except Exception:
                raise StopIteration()

    def commit(self):
        """Commit all of the changes and return the newly
           edited molecule (or MoleculeView)
        """
        if self._connectivity is not None:
            self._molecule.setProperty(self._connectivity_property,
                                       self._connectivity.commit())

        mol = self._molecule.commit()

        try:
            return mol[self._view.index()]
        except Exception:
            return mol

    def keys(self):
        if self._bond is None:
            return self._view.propertyKeys()
        else:
            return self._connectivity.propertyKeys(bond)

    def values(self):
        try:
            if self._bond is None:
                return self._view.propertyValues()
            else:
                return self._connectivity.propertyValues(bond)
        except Exception:
            vals = []

            for key in self.keys():
                vals.append(self.__getitem__(key))

            return vals

    def items(self):
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
