
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
        self._update()

        if self._bond is None:
            return f"Cursor({self._view})"
        else:
            return f"Cursor(bond:{self._bond})"

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

    def atom(self, i):
        """Return the atom in the molecule that matches the passed ID"""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule.atom(i)

        return c

    def residue(self, i):
        """Return the atom in the molecule that matches the passed ID"""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule.residue(i)
        return c

    def chain(self, i):
        """Return the chain in the molecule that matches the passed ID"""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule.chain(i)
        return c

    def segment(self, i):
        """Return the segment in the molecule that matches the passed ID"""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule.segment(i)
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
        self._update()

        c = Cursor()
        c._d = self._d

        try:
            c._view = self._view.parent()
        except Exception:
            c._view = self._molecule

        c._connectivity = None
        c._connectivity_property = None

        c._update()

        return c

    def next(self):
        """Return the cursor to the next logical view (or bond)
           This will go to the next AtomIdx, or next ResIdx,
           or the next bond. This will raise an exception
           (StopIteration) if there is no next view.
        """
        self._update()

        if self._connectivity is None:
            try:
                idx = self._view.index()
                idx += 1

                c = Cursor()
                c._d = self._d
                c._view = self._d.molecule[idx]
                return c
            except Exception:
                raise StopIteration()
        else:
            try:
                bonds = self._connectivity.bonds()
                # find index of current bond...
                raise ValueError()

                c = Cursor()
                c._d = self._d
                c._connectivity_property = self._connectivity_property
                c._bond = next_bond
                c._update()
                return c
            except Exception:
                raise StopIteration()

    def prev(self):
        """Return the cursor to the previous logical view (or bond)
           This will go to the previous AtomIdx, or previous ResIdx,
           or the previous bond. This will raise an exception
           (StopIteration) if there is no previous view.
        """
        self._update()

        if self._connectivity is None:
            try:
                idx = self._view.index()
                idx -= 1

                if idx.value() < 0:
                    raise StopIteration()

                c = Cursor()
                c._d = self._d
                c._view = self._d.molecule[idx]
                return c
            except Exception:
                raise StopIteration()
        else:
            try:
                bonds = self._connectivity.bonds()
                # find index of current bond...
                raise ValueError()

                c = Cursor()
                c._d = self._d
                c._connectivity_property = self._connectivity_property
                c._bond = next_bond
                c._update()
                return c
            except Exception:
                raise StopIteration()

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
