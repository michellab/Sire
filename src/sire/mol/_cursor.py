
from typing import List as _List

__all__ = ["Cursor", "Cursors", "CursorsM"]


class _CursorData:
    """This is the shared data class that holds all of the data
       for a Cursor. This is held by all Cursors that are derived
       from the same object, meaning that multiple Cursors
       can share the same molecule editor. Note that the
       Cursor is not thread-safe (unlike the underlying
       system used by Sire)
    """
    def __init__(self, molecule = None,
                 connectivity_property: str="connectivity",
                 bond_property: str="bond",
                 angle_property: str="angle",
                 dihedral_property: str="dihedral",
                 improper_property: str="improper"):

        if molecule is None:
            self.molecule = None
        else:
            self.molecule = molecule.molecule().edit()

        from ..base import PropertyMap
        self.map = PropertyMap({"connectivity": connectivity_property,
                                "bond": bond_property,
                                "angle": angle_property,
                                "dihedral": dihedral_property,
                                "improper": improper_property})

        self.connectivity_property = connectivity_property

        try:
            self.connectivity = self.molecule.property(
                                    connectivity_property).edit()
        except Exception:
            # the molecule doesn't have a connectivity. Create one for it
            from ..legacy.Mol import CovalentBondHunter
            hunter = CovalentBondHunter()

            try:
                self.connectivity = hunter(self.molecule)
                self.molecule.set_property(connectivity_property,
                                           self.connectivity)
                self.connectivity = self.connectivity.edit()
            except Exception:
                pass

    def remove_internal_property(self, internal, key):
        self.connectivity.remove_property(internal, key)
        self.molecule.set_property(self.connectivity_property,
                                   self.connectivity.commit())

    def set_internal_property(self, internal, key, value):
        if value is None:
            self.remove_bond_property(internal, key)
        else:
            self.connectivity.set_property(internal, key, value)
            self.molecule.set_property(self.connectivity_property,
                                       self.connectivity.commit())

    def set_internal_properties(self, internal, values):
        for key in values.keys():
            value = values[key]

            if value is None:
                self.connectivity.remove_property(internal, key)
            else:
                self.connectivity.set_property(internal, key, value)

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
    def __init__(self, molecule = None, internal = None,
                 connectivity_property: str="connectivity",
                 bond_property: str="bond",
                 angle_property: str="angle",
                 dihedral_property: str="dihedral",
                 improper_property: str="improper"):
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
                              connectivity_property=connectivity_property,
                              bond_property=bond_property,
                              angle_property=angle_property,
                              dihedral_property=dihedral_property,
                              improper_property=improper_property)
        self._view = self._d.update(molecule)

        if (molecule is not None) and (internal is None):
            w = molecule.what()
            if w.endswith("Bond") or w.endswith("Angle") or \
               w.endswith("Dihedral") or w.endswith("Dihedral"):
                internal = molecule.id()

        self._internal = internal

    def _update(self):
        self._view = self._d.update(self._view)

    def __str__(self):
        if self._d.molecule is None:
            return "Cursor::null"
        elif self.is_internal():
            # This is an Internal Cursor
            a0 = self._d.molecule[self._internal.atom0()]
            a1 = self._d.molecule[self._internal.atom1()]

            if self.is_bond():
                return f"Cursor(bond, " \
                    f"{a0.name().value()}:{a0.number().value()} => " \
                    f"{a1.name().value()}:{a1.number().value()})"
            elif self.is_angle():
                a2 = self._d.molecule[self._internal.atom2()]
                return f"Cursor(angle, " \
                    f"{a0.name().value()}:{a0.number().value()} <= " \
                    f"{a1.name().value()}:{a1.number().value()} => " \
                    f"{a2.name().value()}:{a2.number().value()})"
            else:
                a2 = self._d.molecule[self._internal.atom2()]
                a3 = self._d.molecule[self._internal.atom3()]

                if self.is_dihedral():
                    return f"Cursor(dihedral, " \
                        f"{a0.name().value()}:{a0.number().value()} <= " \
                        f"{a1.name().value()}:{a1.number().value()} == " \
                        f"{a2.name().value()}:{a2.number().value()} => " \
                        f"{a3.name().value()}:{a3.number().value()})"
                else:
                    return f"Cursor(improper, " \
                        f"{a0.name().value()}:{a0.number().value()} => " \
                        f"{a1.name().value()}:{a1.number().value()} == " \
                        f"{a2.name().value()}:{a2.number().value()} <= " \
                        f"{a3.name().value()}:{a3.number().value()})"
        else:
            # This is a view cursor
            try:
                return f"Cursor({self.type()}, {self.name}:{self.number})"
            except Exception:
                return f"Cursor({self.type()}, {self.name})"

    def __repr__(self):
        return self.__str__()

    def __delitem__(self, key):
        self._update()

        if self.is_internal():
            self._d.remove_internal_property(self._internal, key)
        else:
            if self.is_molecule():
                # remove the property entirely
                self._d.molecule.remove_property(key)
            else:
                # replace the property with a default-constructed value
                self.__setitem__(key, None)

        self._update()

    def __contains__(self, key):
        self._update()

        if self.is_internal():
            return self._d.connectivity.has_property(self._internal, key)
        else:
            return self._view.has_property(key)

    def __getitem__(self, key):
        self._update()

        if self.is_internal():
            return self._d.connectivity.property(self._internal, key)
        else:
            return self._view.property(key)

    def __setitem__(self, key, value):
        self._update()

        if self.is_internal():
            self._d.set_internal_property(self._internal, key, value)
        else:
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

        self._update()

    def bonds(self, id=None):
        """Return cursors for all of the bonds in this
           view or, if 'id' is supplied, the bonds in this
           view that match 'id'
        """
        self._update()

        cursors = []

        from ..mm import SelectorBond

        if id is None:
            bonds = SelectorBond(self._view, self._d.map)
        else:
            bonds = SelectorBond(self._view, id, self._d.map)

        for bond in bonds:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule
            c._internal = bond.id()
            cursors.append(c)

        return Cursors(self, cursors, bonds)

    def angles(self, id=None):
        """Return cursors for all of the angles in this
           view or, if 'id' is supplied, the angles in this
           view that match 'id'
        """
        self._update()

        cursors = []

        from ..mm import SelectorAngle

        if id is None:
            angles = SelectorAngle(self._view, self._d.map)
        else:
            angles = SelectorAngle(self._view, id, self._d.map)

        for angle in angles:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule
            c._internal = angle.id()
            cursors.append(c)

        return Cursors(self, cursors, angles)

    def dihedrals(self, id=None):
        """Return cursors for all of the dihedrals in this
           view or, if 'id' is supplied, the dihedrals in this
           view that match 'id'
        """
        self._update()

        cursors = []

        from ..mm import SelectorDihedral

        if id is None:
            dihedrals = SelectorDihedral(self._view, self._d.map)
        else:
            dihedrals = SelectorDihedral(self._view, id, self._d.map)

        for dihedral in dihedrals:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule
            c._internal = dihedral.id()
            cursors.append(c)

        return Cursors(self, cursors, dihedrals)

    def impropers(self, id=None):
        """Return cursors for all of the impropers in this
           view or, if 'id' is supplied, the impropers in this
           view that match 'id'
        """
        self._update()

        cursors = []

        from ..mm import SelectorImproper

        if id is None:
            impropers = SelectorImproper(self._view, self._d.map)
        else:
            impropers = SelectorImproper(self._view, id, self._d.map)

        for improper in impropers:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule
            c._internal = improper.id()
            cursors.append(c)

        return Cursors(self, cursors, impropers)

    def _from_view(self, view):
        """Hidden function that allows a single view of the
           current molecule to be converted to a Cursor
        """
        if not self._d.molecule.is_same_molecule(view):
            raise ValueError(
                f"You cannot create from this view ({view}) as it is from "
                f"a different molecule to {self._d.molecule}.")

        c = Cursor()
        c._d = self._d

        from ..mm import Bond, Angle, Dihedral, Improper
        from . import Molecule

        if type(view) is Molecule:
            c._view = self._d.molecule
        elif type(view) in [Bond, Angle, Dihedral, Improper]:
            c._view = self._d.molecule
            c._internal = view.id()
        else:
            c._view = self._d.molecule[view.index()]

        return c

    def _from_views(self, views):
        """Hidden function that allows a set of views, e.g. from
           a Selector_Atom_, SelectorBond etc, to be converted
           to a set of Cursors
        """
        cursors = []

        from sire.mm import Bond, Angle, Dihedral, Improper

        for view in views:
            c = Cursor()
            c._d = self._d

            if type(view) in [Bond, Angle, Dihedral, Improper]:
                c._view = self._d.molecule
                c._internal = view.id()
            else:
                # likely an atom, residue, chain or segment
                c._view = self._d.molecule[view.index()]

            cursors.append(c)

        return Cursors(self, cursors, views)

    def atoms(self, id=None):
        """Return cursors for all of atoms in this view,
           of, if 'id' is supplied, the atoms in this view
           that match 'id'
        """
        cursors = []

        view = self.commit()

        if id is None:
            atoms = view.atoms()
        else:
            atoms = view.atoms(id)

        for atom in atoms:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule.atom(atom.index())
            cursors.append(c)

        return Cursors(self, cursors, atoms)

    def residues(self, id=None):
        """Return cursors for all of residues in this view,
           of, if 'id' is supplied, the residues in this view
           that match 'id'
        """
        cursors = []

        view = self.commit()

        if id is None:
            residues = view.residues()
        else:
            residues = view.residues(id)

        for residue in residues:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule.residue(residue.index())
            cursors.append(c)

        return Cursors(self, cursors, residues)

    def chains(self, id=None):
        """Return cursors for all of chains in this view,
           of, if 'id' is supplied, the chains in this view
           that match 'id'
        """
        cursors = []

        view = self.commit()

        if id is None:
            chains = view.chains()
        else:
            chains = view.chains(id)

        for chain in chains:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule.chain(chain.index())
            cursors.append(c)

        return Cursors(self, cursors, chains)

    def segments(self, id=None):
        """Return cursors for all of segments in this view,
           of, if 'id' is supplied, the segments in this view
           that match 'id'
        """
        cursors = []

        view = self.commit()

        if id is None:
            segments = view.segments()
        else:
            segments = view.segments(id)

        for segment in segments:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule.segment(segment.index())
            cursors.append(c)

        return Cursors(self, cursors, segments)

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
        c._view = self._view.residue(i)
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

    def bond(self, bond):
        """Return the Cursor for the specified bond."""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule
        c._internal = bond

        # make sure that this works
        c._view.bond(bond)

        return c

    def angle(self, angle):
        """Return the Cursor for the specified angle."""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule
        c._internal = angle

        # make sure that this works
        c._view.angle(angle)

        return c

    def dihedral(self, dihedral):
        """Return the Cursor for the specified dihedral."""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule
        c._internal = dihedral

        # make sure that this works
        c._view.dihedral(dihedral)

        return c

    def improper(self, improper):
        """Return the Cursor for the specified improper."""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule
        c._internal = improper

        # make sure that this works
        c._view.improper(improper)

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

        if self.is_internal():
            raise TypeError("An internal (bond/angle/dihedral/improper) does not have a name!")

        return self._view.name().value()

    def set_name(self, name):
        """Set the name of the object in the current view"""
        self._update()

        if self.is_internal():
            raise TypeError("An internal (bond/angle/dihedral/improper) does not have a name!")

        # get the type right...
        orig_name = self._view.name()

        self._view.rename(orig_name.__class__(name))
        self._d.molecule = self._view.molecule()

        return self

    def get_number(self):
        """Return the number of the current view. This returns the
           number as a simple number (it is not a AtomNum, ResNum etc)"""
        self._update()

        if self.is_internal():
            raise TypeError("An internal (bond/angle/dihedral/improper) does not have a number!")

        try:
            return self._view.number().value()
        except Exception:
            raise TypeError(f"A {self._view.what()} does not have a number!")

    def set_number(self, number):
        """Set the number of the object in the current view"""
        self._update()

        if self.is_internal():
            raise TypeError("An internal (bond/angle/dihedral/improper) does not have a number!")

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

        if self.is_internal():
            raise TypeError("An internal (bond/angle/dihedral/improper) does not have an index!")

        return self._view.index().value()

    name = property(get_name, set_name)
    number = property(get_number, set_number)
    index = property(get_index)

    def id(self):
        """Return the ID of this view (e.g. AtomIdx, MolNum, BondID)"""
        self._update()

        if self.is_internal():
            return self._internal

        try:
            return self._view.index()
        except Exception:
            # objects without an index (e.g. molecule) use a number for their ID
            return self._view.number()

    def type(self):
        """Return the type of this Cursor (e.g. 'atom', 'bond',
           'residue', 'chain', 'segment' or 'molecule')
        """
        if self.is_internal():
            if self.is_bond():
                return "bond"
            elif self.is_angle():
                return "angle"
            elif self.is_dihedral():
                return "dihedral"
            elif self.is_improper():
                return "improper"

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

    def is_internal(self):
        """Return whether or not this is a view of an internal
           (i.e. bond, angle, dihedral or improper)
        """
        return self._internal is not None

    def is_bond(self):
        """Return whether this is pointing to a Bond"""
        from . import BondID
        return BondID in type(self._internal).mro()

    def is_angle(self):
        """Return whether this is pointing to an Angle"""
        from . import AngleID
        return AngleID in type(self._internal).mro()

    def is_dihedral(self):
        """Return whether this is pointing to a Dihedral"""
        from . import DihedralID
        return DihedralID in type(self._internal).mro()

    def is_improper(self):
        """Return whether this is pointing to an Improper"""
        from . import ImproperID
        return ImproperID in type(self._internal).mro()

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

        if self.is_internal():
            if self.is_bond():
                from ..mm import Bond
                return Bond(self._d.molecule, self._internal)
            elif self.is_angle():
                from ..mm import Angle
                return Angle(self._d.molecule, self._internal)
            elif self.is_dihedral():
                from ..mm import Dihedral
                return Dihedral(self._d.molecule, self._internal)
            elif self.is_improper():
                from ..mm import Improper
                return Improper(self._d.molecule, self._internal)

        return self._view.evaluate(*args, **kwargs)

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

        if self.is_internal():
            if self.is_bond():
                from ..mm import Bond
                return Bond(self._d.molecule, self._internal)
            elif self.is_angle():
                from ..mm import Angle
                return Angle(self._d.molecule, self._internal)
            elif self.is_dihedral():
                from ..mm import Dihedral
                return Dihedral(self._d.molecule, self._internal)
            elif self.is_improper():
                from ..mm import Improper
                return Improper(self._d.molecule, self._internal)

        try:
            return mol[self._view.index()]
        except Exception:
            return mol

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

        if self.is_internal():
            return self._d.connectivity.property_keys(self._internal)
        else:
            return self._view.property_keys()

    def values(self):
        self._update()

        try:
            if self.is_internal():
                return self._d.connectivity.property_values(self._internal)
            else:
                return self._view.property_values()
        except Exception:
            vals = []

            for key in self.keys():
                vals.append(self.__getitem__(key))

            return vals

    def items(self):
        self._update()

        if self.is_internal():
            keys = self._d.connectivity.property_keys(self._internal)
            items = []

            for key in keys:
                items.append((key, self._d.connectivity.property(self._internal,
                                                                 key)))
        else:
            keys = self._view.property_keys()
            items = []

            for key in keys:
                items.append((key, self._view.property(key)))

        return items

    def properties(self):
        """Return the sire.base.Properties object for the properties
           of the current view
        """
        from ..base import Properties
        p = Properties()

        if self.is_internal():
            return self._d.connectivity.properties(self._internal)
        else:
            for key in self.keys():
                p[key] = self.__getitem__(key)

        return p

    def translate(self, *args, map=None):
        """Translate all of the atoms operated on by this cursor
           by the passed arguments (these are converted automatically
           to a sr.maths.Vector). Use 'map' to specify the property
           map to use to find the coordinates property
        """
        from ..maths import Vector
        delta = Vector(*args)

        if map is None:
            from ..base import PropertyMap
            map = PropertyMap()

        view = self.commit()
        view = view.move().translate(delta, map=map).commit()

        self._d.molecule = view.molecule().edit()
        self._update()

    def rotate(self, angle=None, axis=None, center=None,
               quaternion=None, matrix=None,
               map=None):
        """Rotate all of the atoms operated on by this cursor
           by the passed arguments. Use 'map' to specify the
           property map to use to find the coordinates property.

           There are many ways to specify the rotation, hence
           the number of named arguments:

           angle: (float or angle)
                The angle to rotate by - this is interpreted as
                degrees if you pass in a float. Otherwise use
                sire.units.degrees or sire.units.radians to specify
                the angle unit. This is superseded by the
                quaternion or matrix arguments.

            axis: sire.maths.Vector (or anything that can convert to a Vector)
                The vector about which to rotate. If this is not
                specified, and no other rotation specification is
                used, then the rotation is about the z axis.
                This is superseded by the quaternion or
                matrix arguments.

            center: sire.maths.Vector (or anything that can convert to a Vector)
                The center for the rotation. If this isn't passed then
                the center of mass of the atoms operated on by this
                cursor is used.

            quaternion: sire.maths.Quaternion
                The Quaternion description of the rotation. Note that,
                if you pass this, then the angle, axis and matrix
                arguments will be ignored.

            matrix: sire.maths.Matrix
                The 3x3 rotation matrix that describes the rotation.
                Note that, if you pass this, then the angle and axis
                arguments will be ignored. This is superseded by
                the quaternion argument.

            map: None, dict or sire.base.PropertyMap
                The property map used to find the coordinates property
        """
        from ..maths import create_quaternion, Vector

        quaternion = create_quaternion(angle=angle, axis=axis,
                                       matrix=matrix, quaternion=quaternion)

        view = self.commit()

        if map is None:
            from ..base import PropertyMap
            map = PropertyMap()

        if center is None:
            center = view.evaluate().center_of_mass(map=map)
        else:
            center = Vector(center)

        view = view.move().rotate(quaternion, center, map=map).commit()

        self._d.molecule = view.molecule().edit()
        self._update()


class Cursors:
    """This class holds a list of Cursors. It provides some convenience
       functions that make working with lists of Cursor objects easier.

       This includes being able to commit back to the Cursor that
       created the list, plus being able to apply a function to
       each Cursor in the list
    """
    def __init__(self, parent: Cursor, cursors: _List[Cursor], view):
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
        self._view = view

    def __getitem__(self, i):
        try:
            # do the simple thing if we are asking for the ith cursor
            idx = int(i)
            return self._cursors[idx]
        except Exception:
            pass

        # otherwise we need to search the parent and
        # get the cursors as needed
        view = self._view[i]

        if view.is_selector():
            return self._parent._from_views(view)
        else:
            return self._parent._from_view(view)

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

    def view(self):
        """Return the view underpinning this cursor. This is
           the same as calling '.commit()'
        """
        ret = self._view.clone()

        from . import Molecules
        ret.update(Molecules(self.commit()))
        return ret

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

    def atoms(self, id=None):
        """Return cursors for all of atoms in this view,
           of, if 'id' is supplied, the atoms in this view
           that match 'id'
        """
        if id is None:
            return self._parent._from_views(self._view.atoms())
        else:
            return self._parent._from_views(self._view.atoms(id))

    def residues(self, id=None):
        """Return cursors for all of residues in this view,
           of, if 'id' is supplied, the residues in this view
           that match 'id'
        """
        if id is None:
            return self._parent._from_views(self._view.residues())
        else:
            return self._parent._from_views(self._view.residues(id))

    def chains(self, id=None):
        """Return cursors for all of chains in this view,
           of, if 'id' is supplied, the chains in this view
           that match 'id'
        """
        if id is None:
            return self._parent._from_views(self._view.chains())
        else:
            return self._parent._from_views(self._view.chains(id))

    def segments(self, id=None):
        """Return cursors for all of segments in this view,
           of, if 'id' is supplied, the segments in this view
           that match 'id'
        """
        if id is None:
            return self._parent._from_views(self._view.segments())
        else:
            return self._parent._from_views(self._view.segments(id))

    def bonds(self, id=None):
        """Return cursors for all of the bonds in this
           view or, if 'id' is supplied, the bonds in this
           view that match 'id'
        """
        if id is None:
            return self._parent._from_views(self._view.bonds())
        else:
            return self._parent._from_views(self._view.bonds(id))

    def angles(self, id=None):
        """Return cursors for all of the angles in this
           view or, if 'id' is supplied, the angles in this
           view that match 'id'
        """
        if id is None:
            return self._parent._from_views(self._view.angles())
        else:
            return self._parent._from_views(self._view.angles(id))

    def dihedrals(self, id=None):
        """Return cursors for all of the dihedrals in this
           view or, if 'id' is supplied, the dihedrals in this
           view that match 'id'
        """
        if id is None:
            return self._parent._from_views(self._view.dihedrals())
        else:
            return self._parent._from_views(self._view.dihedrals(id))

    def impropers(self, id=None):
        """Return cursors for all of the impropers in this
           view or, if 'id' is supplied, the impropers in this
           view that match 'id'
        """
        if id is None:
            return self._parent._from_views(self._view.impropers())
        else:
            return self._parent._from_views(self._view.impropers(id))

    def atom(self, i):
        """Return the atom in the molecule that matches the passed ID"""
        return self._parent._from_view(self._view.atom(i))

    def residue(self, i=None):
        """Return the residue in the molecule that matches the passed ID.
           If 'i' is None, then this returns the residue that contains
           this atom (if this is an atom)
        """
        if i is None:
            return self._parent._from_view(self._view.residue())
        else:
            return self._parent._from_view(self._view.residue(i))

    def chain(self, i=None):
        """Return the chain in the molecule that matches the passed ID.
           If 'i' is None, then this returns the chain that contains
           this atom (if this is an atom)"""
        if i is None:
            return self._parent._from_view(self._view.chain())
        else:
            return self._parent._from_view(self._view.chain(i))

    def segment(self, i=None):
        """Return the segment in the molecule that matches the passed ID.
           If 'i' is None, then this returns the segment that contains
           this atom (if this is an atom)"""
        if i is None:
            return self._parent._from_view(self._view.segment())
        else:
            return self._parent._from_view(self._view.segment(i))

    def molecule(self):
        """Return the molecule"""
        return self._parent.molecule()

    def bond(self, bond):
        """Return the Cursor for the specified bond."""
        if bond is None:
            return self._parent._from_view(self._view.bond())
        else:
            return self._parent._from_view(self._view.bond(bond))

    def angle(self, angle):
        """Return the Cursor for the specified angle."""
        if angle is None:
            return self._parent._from_view(self._view.angle())
        else:
            return self._parent._from_view(self._view.angle(angle))

    def dihedral(self, dihedral):
        """Return the Cursor for the specified dihedral."""
        if dihedral is None:
            return self._parent._from_view(self._view.dihedral())
        else:
            return self._parent._from_view(self._view.dihedral(dihedral))

    def improper(self, improper):
        """Return the Cursor for the specified improper."""
        if improper is None:
            return self._parent._from_view(self._view.improper())
        else:
            return self._parent._from_view(self._view.improper(improper))

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

    def translate(self, *args, map=None):
        """Translate all of the atoms operated on by these cursors
           by the passed arguments (these are converted automatically
           to a sr.maths.Vector). Use 'map' to specify the property
           map to use to find the coordinates property
        """
        for cursor in self._cursors:
            cursor.translate(*args, map=map)

    def rotate(self, angle=None, axis=None, center=None,
               quaternion=None, matrix=None,
               map=None):
        """Rotate all of the atoms operated on by this cursor
           by the passed arguments. Use 'map' to specify the
           property map to use to find the coordinates property.

           There are many ways to specify the rotation, hence
           the number of named arguments:

           angle: (float or angle)
                The angle to rotate by - this is interpreted as
                degrees if you pass in a float. Otherwise use
                sire.units.degrees or sire.units.radians to specify
                the angle unit. This is superseded by the
                quaternion or matrix arguments.

            axis: sire.maths.Vector (or anything that can convert to a Vector)
                The vector about which to rotate. If this is not
                specified, and no other rotation specification is
                used, then the rotation is about the z axis.
                This is superseded by the quaternion or
                matrix arguments.

            center: sire.maths.Vector (or anything that can convert to a Vector)
                The center for the rotation. If this isn't passed then
                the center of mass of the atoms operated on by this
                cursor is used.

            quaternion: sire.maths.Quaternion
                The Quaternion description of the rotation. Note that,
                if you pass this, then the angle, axis and matrix
                arguments will be ignored.

            matrix: sire.maths.Matrix
                The 3x3 rotation matrix that describes the rotation.
                Note that, if you pass this, then the angle and axis
                arguments will be ignored. This is superseded by
                the quaternion argument.

            map: None, dict or sire.base.PropertyMap
                The property map used to find the coordinates property
        """
        from ..maths import create_quaternion

        quaternion = create_quaternion(angle=angle, axis=axis,
                                       matrix=matrix, quaternion=quaternion)

        for cursor in self._cursors:
            cursor.rotate(quaternion=quaternion, center=center, map=map)


class CursorsM:
    """This class holds a list of Cursor/Cursors that operate across
       multiple molecules. This allows you to perform editing
       operations across many molecules at the same time.
    """
    def __init__(self, parent=None):
        self._parent = None
        self._cursors = []

        self._molcursors = {}

        if parent is not None:
            self._parent = parent.clone()

            for child in parent:
                child_mol = child.molecule()
                molnum = child_mol.number()

                if molnum not in self._molcursors:
                    self._molcursors[molnum] = child_mol.cursor()

                self._cursors.append(
                        self._molcursors[molnum]._from_view(child))

    def _from_view(self, view):
        """Internal function that constructs from a single view"""
        molnum = view.molecule().number()
        molcursor = self._molcursors[molnum]

        c = Cursor()
        c._d = molcursor._d

        from ..mm import Bond, Angle, Dihedral, Improper
        from . import Molecule

        if type(view) is Molecule:
            c._view = molcursor._d.molecule
        elif type(view) in [Bond, Angle, Dihedral, Improper]:
            c._view = molcursor._d.molecule
            c._internal = view.id()
        else:
            c._view = molcursor._d.molecule[view.index()]

        c._update()

        return c

    def _from_views(self, views):
        """Internal function to construct from a set of views"""
        ret = CursorsM()
        ret._parent = views.clone()
        ret._molcursors = self._molcursors

        for view in views:
            molnum = view.molecule().number()
            ret._cursors.append(ret._molcursors[molnum]._from_view(view))

        return ret

    def __getitem__(self, i):
        try:
            # try the simplest case - the ith cursor
            idx = int(i)
            return self._cursors[idx]
        except Exception:
            pass

        return self._from_views(self._parent[i])

    def __len__(self):
        return len(self._cursors)

    def __str__(self):
        if len(self) == 0:
            return "CursorsM::null"
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

            return f"CursorsM( size={n}\n{lines}\n)"

    def __repr__(self):
        return self.__str__()

    def commit(self):
        """Commit all of the changes and return the newly
           edited multi-molecule view.
        """
        from . import Molecules

        updated = Molecules()
        updated.reserve(len(self._cursors))

        for cursor in self._cursors:
            updated.add(cursor.commit())

        self._parent.update(updated)
        return self._parent

    def atoms(self, id=None):
        """Return cursors for all of atoms in this view,
           or, if 'id' is supplied, the atoms in this view
           that match 'id'
        """
        if id is None:
            return self._from_views(self._parent.atoms())
        else:
            return self._from_views(self._parent.atoms(id))

    def residues(self, id=None):
        """Return cursors for all of residues in this view,
           or, if 'id' is supplied, the residues in this view
           that match 'id'
        """
        if id is None:
            return self._from_views(self._parent.residues())
        else:
            return self._from_views(self._parent.residues(id))

    def chains(self, id=None):
        """Return cursors for all of chains in this view,
           or, if 'id' is supplied, the chains in this view
           that match 'id'
        """
        if id is None:
            return self._from_views(self._parent.chains())
        else:
            return self._from_views(self._parent.chains(id))

    def segments(self, id=None):
        """Return cursors for all of segments in this view,
           or, if 'id' is supplied, the segments in this view
           that match 'id'
        """
        if id is None:
            return self._from_views(self._parent.segments())
        else:
            return self._from_views(self._parent.segments(id))

    def molecules(self, id=None):
        """Return cursors for all of the molecules in this view,
           or, if 'id' is supplied, the molecules in this view
           that match 'id'
        """
        if id is None:
            return self._from_views(self._parent.molecules())
        else:
            return self._from_views(self._parent.molecules(id))

    def bonds(self, id=None):
        """Return cursors for all of the bonds in this
           view or, if 'id' is supplied, the bonds in this
           view that match 'id'
        """
        if id is None:
            return self._from_views(self._parent.bonds())
        else:
            return self._from_views(self._parent.bonds(id))

    def angles(self, id=None):
        """Return cursors for all of the angles in this
           view or, if 'id' is supplied, the angles in this
           view that match 'id'
        """
        if id is None:
            return self._from_views(self._parent.angles())
        else:
            return self._from_views(self._parent.angles(id))

    def dihedrals(self, id=None):
        """Return cursors for all of the dihedrals in this
           view or, if 'id' is supplied, the dihedrals in this
           view that match 'id'
        """
        if id is None:
            return self._from_views(self._parent.dihedrals())
        else:
            return self._from_views(self._parent.dihedrals(id))

    def impropers(self, id=None):
        """Return cursors for all of the impropers in this
           view or, if 'id' is supplied, the impropers in this
           view that match 'id'
        """
        if id is None:
            return self._from_views(self._parent.impropers())
        else:
            return self._from_views(self._parent.impropers(id))

    def atom(self, i):
        """Return the atom in the molecule that matches the passed ID"""
        return self._from_view(self._parent.atom(i))

    def residue(self, i=None):
        """Return the residue in the molecule that matches the passed ID.
           If 'i' is None, then this returns the residue that contains
           this atom (if this is an atom)
        """
        if i is None:
            return self._from_view(self._parent.residue())
        else:
            return self._from_view(self._parent.residue(i))

    def chain(self, i=None):
        """Return the chain in the molecule that matches the passed ID.
           If 'i' is None, then this returns the chain that contains
           this atom (if this is an atom)"""
        if i is None:
            return self._from_view(self._parent.chain())
        else:
            return self._from_view(self._parent.chain(i))

    def segment(self, i=None):
        """Return the segment in the molecule that matches the passed ID.
           If 'i' is None, then this returns the segment that contains
           this atom (if this is an atom)"""
        if i is None:
            return self._from_view(self._parent.segment())
        else:
            return self._from_view(self._parent.segment(i))

    def molecule(self, i=None):
        """Return the molecule"""
        if i is None:
            mols = self._parent.molecules()
            if len(mols) != 1:
                raise ValueError(
                    f"There is more than one molecule in this view ({len(mols)}) ."
                    "You need to specify which one you want."
                )

            return self._from_view(mols[0])
        else:
            return self._from_view(self._parent.molecule(i))

    def bond(self, bond):
        """Return the Cursor for the specified bond."""
        if bond is None:
            return self._from_view(self._parent.bond())
        else:
            return self._from_view(self._parent.bond(bond))

    def angle(self, angle):
        """Return the Cursor for the specified angle."""
        if angle is None:
            return self._from_view(self._parent.angle())
        else:
            return self._from_view(self._parent.angle(angle))

    def dihedral(self, dihedral):
        """Return the Cursor for the specified dihedral."""
        if dihedral is None:
            return self._from_view(self._parent.dihedral())
        else:
            return self._from_view(self._parent.dihedral(dihedral))

    def improper(self, improper):
        """Return the Cursor for the specified improper."""
        if improper is None:
            return self._from_view(self._parent.improper())
        else:
            return self._from_view(self._parent.improper(improper))

    def apply(self, func, *args, **kwargs):
        """Apply the passed function (with optional position and keyword
           arguments) to all of the cursors in this list of Cursors.
           As the function is intended to use
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

    def translate(self, *args, map=None):
        """Translate all of the atoms operated on by these cursors
           by the passed arguments (these are converted automatically
           to a sr.maths.Vector). Use 'map' to specify the property
           map to use to find the coordinates property
        """
        for cursor in self._cursors:
            cursor.translate(*args, map=map)

    def rotate(self, angle=None, axis=None, center=None,
               quaternion=None, matrix=None,
               map=None):
        """Rotate all of the atoms operated on by this cursor
           by the passed arguments. Use 'map' to specify the
           property map to use to find the coordinates property.

           There are many ways to specify the rotation, hence
           the number of named arguments:

           angle: (float or angle)
                The angle to rotate by - this is interpreted as
                degrees if you pass in a float. Otherwise use
                sire.units.degrees or sire.units.radians to specify
                the angle unit. This is superseded by the
                quaternion or matrix arguments.

            axis: sire.maths.Vector (or anything that can convert to a Vector)
                The vector about which to rotate. If this is not
                specified, and no other rotation specification is
                used, then the rotation is about the z axis.
                This is superseded by the quaternion or
                matrix arguments.

            center: sire.maths.Vector (or anything that can convert to a Vector)
                The center for the rotation. If this isn't passed then
                the center of mass of the atoms operated on by this
                cursor is used.

            quaternion: sire.maths.Quaternion
                The Quaternion description of the rotation. Note that,
                if you pass this, then the angle, axis and matrix
                arguments will be ignored.

            matrix: sire.maths.Matrix
                The 3x3 rotation matrix that describes the rotation.
                Note that, if you pass this, then the angle and axis
                arguments will be ignored. This is superseded by
                the quaternion argument.

            map: None, dict or sire.base.PropertyMap
                The property map used to find the coordinates property
        """
        from ..maths import create_quaternion

        quaternion = create_quaternion(angle=angle, axis=axis,
                                       matrix=matrix, quaternion=quaternion)

        for cursor in self._cursors:
            cursor.rotate(quaternion=quaternion, center=center, map=map)
