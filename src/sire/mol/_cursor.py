
from typing import List as _List

__all__ = ["Cursor", "Cursors"]


class _CursorData:
    """This is the shared data class that holds all of the data
       for a Cursor. This is held by all Cursors that are derived
       from the same object, meaning that multiple Cursors
       can share the same molecule editor. Note that the
       Cursor is not thread-safe (unlike the underlying
       system used by Sire)
    """
    def __init__(self, molecule = None,
                 connectivity_property: str="connectivity"):
        if molecule is None:
            self.molecule = None
        else:
            self.molecule = molecule.molecule().edit()

        self.connectivity_property = connectivity_property

        try:
            self.connectivity = self.molecule.property(
                                    self.connectivity_property).edit()
        except Exception:
            # the molecule doesn't have a connectivity. Create one for it
            from ..legacy.Mol import CovalentBondHunter
            hunter = CovalentBondHunter()

            try:
                self.connectivity = hunter(self.molecule)
                self.molecule.set_property(self.connectivity_property,
                                           self.connectivity)
                self.connectivity = self.connectivity.edit()
            except Exception:
                pass

    def remove_bond_property(self, bond, key):
        self.connectivity.remove_property(bond, key)
        self.molecule.set_property(self.connectivity_property,
                                   self.connectivity.commit())

    def set_bond_property(self, bond, key, value):
        if value is None:
            self.remove_bond_property(bond, key)
        else:
            self.connectivity.set_property(bond, key, value)
            self.molecule.set_property(self.connectivity_property,
                                       self.connectivity.commit())

    def set_bond_properties(self, bond, values):
        for key in values.keys():
            value = values[key]

            if value is None:
                self.connectivity.remove_property(bond, key)
            else:
                self.connectivity.set_property(bond, key, value)

        self.molecule.set_property(self.connectivity_property,
                                   self.connectivity.commit())

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
        self._d = _CursorData(molecule=molecule,
                              connectivity_property=connectivity_property)
        self._view = self._d.update(molecule)

        if (molecule is not None) and (bond is None):
            if molecule.what().endswith("Bond"):
                bond = molecule.id()

        self._bond = bond

    def _update(self):
        self._view = self._d.update(self._view)

    def __str__(self):
        if self._d.molecule is None:
            return "Cursor::null"
        elif self._bond is None:
            # This is a view cursor
            try:
                return f"Cursor({self.type()}, {self.name}:{self.number})"
            except Exception:
                return f"Cursor({self.type()}, {self.name})"
        else:
            # This is a bond Cursor
            a0 = self._d.molecule[self._bond.atom0()]
            a1 = self._d.molecule[self._bond.atom1()]
            return f"Cursor({self.type()}, " \
                   f"{a0.name().value()}:{a0.number().value()} => " \
                   f"{a1.name().value()}:{a1.number().value()})"

    def __repr__(self):
        return self.__str__()

    def __delitem__(self, key):
        self._update()

        if self._bond is None:
            if self.is_molecule():
                # remove the property entirely
                self._d.molecule.remove_property(key)
            else:
                # replace the property with a default-constructed value
                self.__setitem__(key, None)
        else:
            self._d.remove_bond_property(self._bond, key)

        self._update()

    def __contains__(self, key):
        self._update()

        if self._bond is None:
            return self._view.has_property(key)
        else:
            return self._d.connectivity.has_property(self._bond, key)

    def __getitem__(self, key):
        self._update()

        if self._bond is None:
            return self._view.property(key)
        else:
            return self._d.connectivity.property(self._bond, key)

    def __setitem__(self, key, value):
        self._update()

        if self._bond is None:
            if value is None:
                # we need to create a default-constructed value
                try:
                    v = self._view.property(key)
                    value = v.__class__()
                except KeyError:
                    # No existing property - assume this is Boolean False
                    value = False

            self._view.set_property(key, value)
            self._d.molecule = self._view.molecule()
        else:
            self._d.set_bond_property(self._bond, key, value)

        self._update()

    def bonds(self):
        """Return cursors for all of the bonds in this
           view or, if 'id' is supplied, the bonds in this
           view that match 'id'
        """
        self._update()

        cursors = []

        from ..mm import SelectorBond
        bonds = SelectorBond(self._view, {"connectivity": self._d.connectivity_property})

        for bond in bonds:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule
            c._bond = bond.id()
            cursors.append(c)

        return Cursors(self, cursors)

    def _views(self, views):
        """Hidden function that allows a set of views, e.g. from
           a Selector_Atom_, SelectorBond etc, to be converted
           to a set of Cursors
        """
        self._update()

        cursors = []

        from sire.mm import Bond

        for view in views:
            c = Cursor()
            c._d = self._d

            if type(view) is Bond:
                # likely a bond
                c._view = self._d.molecule
                c._bond = view.id()
            else:
                # likely an atom, residue, chain or segment
                c._view = self._d.molecule[view.index()]

            cursors.append(c)

        return Cursors(self, cursors)

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

        return Cursors(self, cursors)

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

        return Cursors(self, cursors)

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

        return Cursors(self, cursors)

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

        return Cursors(self, cursors)

    def atom(self, i):
        """Return the atom in the molecule that matches the passed ID"""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._view.atom(i)
        c._update()

        return c

    def residue(self, i=None):
        """Return the residue in the molecule that matches the passed ID.
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
                    f"There is no residue that contains {self.type()}:{self.id()}"
                )

        self._update()

        c = Cursor()
        c._d = self._d
        c._view = c._view.residue(i)
        c._update()

        return c

    def chain(self, i=None):
        """Return the chain in the molecule that matches the passed ID.
           If 'i' is None, then this returns the chain that contains
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
                    f"There is no chain that contains {self.type()}:{self.id()}"
                )

        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._view.chain(i)
        c._update()

        return c

    def segment(self, i=None):
        """Return the segment in the molecule that matches the passed ID.
           If 'i' is None, then this returns the segment that contains
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
                    f"There is no segment that contains {self.type()}:{self.id()}"
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
        c._bond = bond

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

    def get(self, key):
        """Return the property associated with key 'key'"""
        return self.__getitem__(key)

    def set(self, key, value):
        """Set the property associated with key 'key' to the
           passed value
        """
        self.__setitem__(key, value)
        return self

    def delete(self, key):
        """Remove the property associated with the key 'key'"""
        self.__delitem__(key)
        return self

    def get_name(self):
        """Return the name of the current view. Note that this
           returns the name as a simple string (it is not a
           AtomName, ResName etc)"""
        self._update()

        if self._bond is not None:
            raise TypeError("A bond does not have a name!")

        return self._view.name().value()

    def set_name(self, name):
        """Set the name of the object in the current view"""
        self._update()

        if self._bond is not None:
            raise TypeError("A bond does not have a name!")

        # get the type right...
        orig_name = self._view.name()

        self._view.rename(orig_name.__class__(name))
        self._d.molecule = self._view.molecule()

        return self

    def get_number(self):
        """Return the number of the current view. This returns the
           number as a simple number (it is not a AtomNum, ResNum etc)"""
        self._update()

        if self._bond is not None:
            raise TypeError("A bond does not have a number!")

        try:
            return self._view.number().value()
        except Exception:
            raise TypeError(f"A {self._view.what()} does not have a number!")

    def set_number(self, number):
        """Set the number of the object in the current view"""
        self._update()

        if self._bond is not None:
            raise TypeError("A bond does not have a number!")

        try:
            orig_number = self._view.number()
            self._view.renumber(orig_number.__class__(number))
            self._d.molecule = self._view.molecule()
        except AttributeError:
            raise TypeError(f"A {self._view.what()} does not have a number!")

        return self

    def get_index(self):
        """Return the index of the current view. This returns it as
           as simple number (i.e. not as an AtomIdx, ResIdx etc)"""
        self._update()

        if self._bond is not None:
            raise TypeError("A bond does not have an index!")

        return self._view.index().value()

    name = property(get_name, set_name)
    number = property(get_number, set_number)
    index = property(get_index)

    def id(self):
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
        if self.is_bond():
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

    def is_molecule(self):
        """Return whether this is pointing to a Molecule"""
        return self.type() == "molecule"

    def is_bond(self):
        """Return whether this is pointing to a Bond"""
        self._update()

        return self._bond is not None

    def is_atom(self):
        """Return whether this is pointing to an Atom"""
        return self.type() == "atom"

    def is_residue(self):
        """Return whether this is pointing to a Residue"""
        return self.type() == "residue"

    def is_chain(self):
        """Return whether this is pointing to a Chain"""
        return self.type() == "chain"

    def is_segment(self):
        """Return whether this is pointing to a Segment"""
        return self.type() == "segment"

    def evaluate(self, *args, **kwargs):
        """Return a :class:`sire.mol.Evaluator` for the view
           in this cursor
        """
        self._update()

        if self._bond is None:
            return self._view.evaluate(*args, **kwargs)
        else:
            from sire.mm import Bond
            return Bond(self._d.molecule, self._bond)

    def view(self):
        """Return the view underpinning this cursor. This is actually
           the same as committing the changes
        """
        return self.commit()

    def commit(self):
        """Commit all of the changes and return the newly
           edited molecule (or MoleculeView)
        """
        self._update()

        mol = self._d.molecule.commit()

        if self._bond is None:
            try:
                return mol[self._view.index()]
            except Exception:
                return mol
        else:
            from ..mm import Bond
            return Bond(self._d.molecule, self._bond.atom0(), self._bond.atom1())

    def apply(self, func, *args, **kwargs):
        """Apply the passed function (with optional position and keyword
           arguments) to this Cursor. As the function is intended to use
           the Cursor to edit molecules, only this Cursor will be returned.
           This lets you run `.apply(...).commit()` as a single line.

           The function can be either;

           1. a string containing the name of the function to call, or
           2. an actual function (either a normal function or a lambda expression)

           You can optionally pass in positional and keyword arguments
           here that will be passed to the function.

           Args:
               func (str or function): The function to be called, or the name
                                       of the function to be called.

           Returns:
               Cursor: This cursor
        """
        if str(func) == func:
            # we calling a named function
            func = getattr(self, func)

        func(self, *args, **kwargs)
        return self

    def keys(self):
        self._update()

        if self._bond is None:
            return self._view.property_keys()
        else:
            return self._d.connectivity.property_keys(self._bond)

    def values(self):
        self._update()

        try:
            if self._bond is None:
                return self._view.property_values()
            else:
                return self._d.connectivity.property_values(self._bond)
        except Exception:
            vals = []

            for key in self.keys():
                vals.append(self.__getitem__(key))

            return vals

    def items(self):
        self._update()

        if self._bond is None:
            keys = self._view.property_keys()
            items = []

            for key in keys:
                items.append((key, self._view.property(key)))
        else:
            keys = self._d.connectivity.property_keys(self._bond)
            items = []

            for key in keys:
                items.append((key, self._d.connectivity.property(self._bond,
                                                                 key)))

        return items

    def properties(self):
        """Return the sire.base.Properties object for the properties
           of the current view
        """
        from ..base import Properties
        p = Properties()

        if self._bond is None:
            for key in self.keys():
                p[key] = self.__getitem__(key)
        else:
            return self._d.connectivity.properties(self._bond)

        return p


class Cursors:
    """This class holds a list of Cursors. It provides some convenience
       functions that make working with lists of Cursor objects easier.

       This includes being able to commit back to the Cursor that
       created the list, plus being able to apply a function to
       each Cursor in the list
    """
    def __init__(self, parent: Cursor, cursors: _List[Cursor]):
        if type(parent) is not Cursor:
            raise TypeError(f"{parent} must be a Cursor object!")

        for c in cursors:
            if type(c) is not Cursor:
                raise TypeError(f"{c} must be a Cursor object!")

            if c._d is not parent._d:
                raise ValueError(
                    f"The list of cursors must be created from the parent!")

        self._parent = parent
        self._cursors = cursors

    def __getitem__(self, i):
        try:
            idx = int(i)
        except Exception:
            idx = None

        if idx is not None:
            return self._cursors[idx]
        elif type(i) is slice:
                return Cursors(self._parent, self._cursors[i])
        else:
            cs = []

            for idx in i:
                cs.append(self._cursors[idx])

            return Cursors(self._parent, cs)

    def __delitem__(self, i):
        try:
            # if this is an integer, then delete the ith cursor
            self._cursors.__delitem__(i)
            return
        except Exception:
            pass

        # delete this key from all of the cursors
        for c in self._cursors:
            del c[i]

    def __len__(self):
        return len(self._cursors)

    def __str__(self):
        if len(self) == 0:
            return "Cursors::null"
        else:
            lines = []

            n = len(self._cursors)

            if n <= 10:
                for i in range(0, n):
                    lines.append(f"{i+1}: {self._cursors[i]}")
            else:
                for i in range(0, 5):
                    lines.append(f"{i+1}: {self._cursors[i]}")

                lines.append("...")

                for i in range(n-5, n):
                    lines.append(f"{i+1}: {self._cursors[i]}")

            lines = "\n".join(lines)

            return f"Cursors( size={n}\n{lines}\n)"

    def __repr__(self):
        return self.__str__()

    def delete(self, i):
        """Remove either the ith cursor in the list, or i is a string,
           delete that key from all of the cursors
        """
        self.__delitem__(i)
        return self

    def commit(self):
        """Commit all of the changes and return the newly
           edited molecule (or MoleculeView). This commits
           on the parent Cursor that was used to create
           this list, e.g.

           >>> mol.cursor().atoms().commit()

           will commit and return the updated molecule (mol).

           This is equivalent to `self.parent().commit()`
        """
        return self._parent.commit()

    def parent(self):
        """Return the parent cursor"""
        return self._parent

    def apply(self, func, *args, **kwargs):
        """Apply the passed function (with optional position and keyword
           arguments) to all of the cursors in this list of Cursors
           (i.e. everything except the parent). As the function is intended to use
           the Cursor to edit molecules, only this Cursors object will be returned.
           This lets you run `.apply(...).commit()` as a single line.

           The function can be either;

           1. a string containing the name of the function to call, or
           2. an actual function (either a normal function or a lambda expression)

           You can optionally pass in positional and keyword arguments
           here that will be passed to the function.

           Args:
               func (str or function): The function to be called, or the name
                                       of the function to be called.

           Returns:
               Cursors: This list of cursors
        """
        for cursor in self._cursors:
            cursor.apply(func, *args, **kwargs)

        return self
