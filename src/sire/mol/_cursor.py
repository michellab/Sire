
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
    def __init__(self, molecule = None, map=None):
        if molecule is None:
            self.molecule = None
            self.map = None
            self.connectivity = None
            self.connectivity_property = None
            return
        else:
            self.molecule = molecule.molecule().edit()

        from ..base import PropertyMap

        if map is None:
            self.map = PropertyMap()
        else:
            self.map = PropertyMap(map)

        self.connectivity_property = self.map["connectivity"].source()

        if self.connectivity_property is None:
            # we cannot support value-based connectivity properties
            self.map.set("connectivity", "connectivity")
            self.connectivity_property = self.map["connectivity"].source()

        try:
            self.connectivity = self.molecule.property(
                                    self.connectivity_property).edit()
        except Exception:
            # the molecule doesn't have a connectivity. Create one for it
            from ..legacy.Mol import CovalentBondHunter
            hunter = CovalentBondHunter()

            try:
                connectivity = hunter(self.molecule)
                self.molecule.set_property(self.connectivity_property,
                                           connectivity)
                self.connectivity = connectivity.edit()
            except Exception as e:
                from ..utils import Console
                Console.warning(
                    f"Cannot auto-generate a connectivity for {self.molecule}. "
                    f"The error is:\n\n{e}")

    def number(self):
        """Return the molnum number of the molecule being edited
           by this cursor
        """
        return self.molecule.number()

    def merge(self, map):
        """Return a property map that is the combination
           of self.map and the passed map. The properties
           set in the passed map have precedence.
        """
        if map is None:
            return self.map
        else:
            return self.map.merge(map)

    def remove_internal_property(self, internal, key):
        self.connectivity.remove_property(internal, key)
        self.molecule.set_property(self.connectivity_property,
                                   self.connectivity.commit())

    def set_internal_property(self, internal, key, value):
        if value is None:
            self.remove_internal_property(internal, key)
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


def _process_move_options(view, anchor, weighting, map):
    """Internal function used to process the passed move options
       and return a property map with those options set.

       The default is to have no anchors and to have
       the AbsFromNumber weighting (this is defined
       in weightfunction.h in the C++ layer).

       Here we convert the weighting (as a string option)
       into a WeightFunction, plus we convert the anchor
       (as a view and then selection into that view) into
       an AtomSelection object
    """
    if anchor is None and weighting is None:
        return map

    from ..base import PropertyMap

    if map is None:
        map = PropertyMap()
    else:
        map = PropertyMap(map)

    if anchor is not None:
        try:
            selection = view[anchor].selection()
            map.set("anchors", selection)
        except Exception as e:
            from ..utils import Console
            Console.warning(
                f"Unable to find anchors '{anchor}'.\n\n{e}")

    if weighting is not None:
        from ..legacy.Mol import AbsFromMass, RelFromMass, \
                                 AbsFromNumber, RelFromNumber

        weightfuncs = {
            "relative_mass": RelFromMass,
            "absolute_mass": AbsFromMass,
            "relative_number": RelFromNumber,
            "absolute_number": AbsFromNumber
        }

        try:
            weighting = weightfuncs[weighting]()
        except Exception:
            raise ValueError(
                f"Unsupported weighting: {weighting}. Supported values "
                f"are {weightfuncs.keys()}.")

        map.set("weight function", weighting)

    return map


class Cursor:
    """This class provides a cursor that can be used to navigate through
       and edit the properties of Molecules. This makes the whole
       getting and setting of properties more pythonic in writing
       style, while also saving some typing.
    """
    def __init__(self, molecule = None, internal = None, map=None):
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
        self._d = _CursorData(molecule=molecule, map=map)
        self._view = self._d.update(molecule)

        if (molecule is not None) and (internal is None):
            w = molecule.what()
            if w.endswith("Bond") or w.endswith("Angle") or \
               w.endswith("Dihedral") or w.endswith("Dihedral"):
                internal = molecule.id()

        self._internal = internal
        self._add_extra_functions()

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
        """Delete the property with specified key
        """
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
        """Return whether a property with the specified key is
           contained in this view.
        """
        self._update()

        if self.is_internal():
            return self._d.connectivity.has_property(self._internal, key)
        else:
            return self._view.has_property(key)

    def __call__(self, key):
        """Return a cursor that represents the sub-view of this
           cursor, indexed by key. For example,
           cursor("element C") would return a cursor for all
           of the carbon atoms in this view.
        """
        view = self.view()[key]

        if view.is_selector():
            return self._from_views(view)
        else:
            return self._from_view(view)

    def __getitem__(self, key):
        """Return the property that matches the passed key, OR the
           sub-view that matches the key. This will only look for
           the sub-view if there is no matching property. Use
           the __call__ function if you only want to search for
           sub-views.
        """
        if type(key) is not str:
            return self.__call__(key)

        try:
            return self.get(key)
        except Exception as e:
            property_error = e

        # We can't find the property so try the sub-view
        try:
            return self.__call__(key)
        except Exception:
            pass

        # We can't find the sub-view, but since this searches the
        # property first, we will raise the property exception
        raise property_error

    def __setitem__(self, key, value):
        """Set the property with key 'key' to the passed 'value'
        """
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

    def _add_bond_functions(self):
        """Internal function used to add member functions that are
           specific only to cursors that operate on bonds
        """
        def get_length(map=None):
            """Return the length of the bond being edited by this cursor"""
            map = self._d.merge(map)
            return self.view().length(map=map)

        def set_length(value, anchor=None, weighting=None,
                       auto_align=True, map=None):
            """Set the length of the bond being edited by this cursor to
               'value'. This should be either a length unit, or a float
               (in which case it is converted into a value with default
               length units)

                value: float or length
                    The length to which to set this bond.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            from ..units import length

            value = length(value)

            view = self._d.molecule.commit()

            map = _process_move_options(view=view, anchor=anchor,
                                        weighting=weighting, map=map)
            map = self._d.merge(map)

            moved = view.move().set(self._internal, value, map)

            if auto_align and (anchor is None):
                moved.align(view, map)

            self._d.molecule = moved.commit().molecule().edit()
            self._update()

            return self

        def change_length(delta, anchor=None, weighting=None,
                          auto_align=True, map=None):
            """Change the length of the bond being edited by this cursor by
               'delta'. This should be either a length unit, or a float
               (in which case it is converted into a value with default
               length units)

                delta: float or length
                    The length by which to change this bond.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            from ..units import length

            delta = length(delta)

            view = self._d.molecule.commit()

            map = _process_move_options(view=view, anchor=anchor,
                                        weighting=weighting, map=map)
            map = self._d.merge(map)

            moved = view.move().change(self._internal, delta, map)

            if auto_align and (anchor is None):
                moved.align(view, map)

            self._d.molecule = moved.commit().molecule().edit()
            self._update()

            return self

        self.length = get_length
        self.set_length = set_length
        self.change_length = change_length
        self.measure = get_length
        self.set_measure = set_length
        self.change_measure = change_length

    def _add_angle_functions(self):
        """Internal function used to add member functions that are
           specific only to cursors that operate on angles, dihedrals
           or impropers
        """
        def get_size(map=None):
            """Return the size of the internal being edited by this cursor"""
            map = self._d.merge(map)
            return self.view().size(map=map)

        def set_size(value, anchor=None, weighting=None,
                     auto_align=True, move_all=True, map=None):
            """Set the size of the internal being edited by this cursor to
               'delta'. This should be either an angle unit, or a float
               (in which case it is converted into a value with default
                angle units)

                value: float or angle
                    The angle to which to set this internal.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            from ..units import angle

            value = angle(value)

            view = self._d.molecule.commit()

            map = _process_move_options(view=view,
                                        anchor=anchor,
                                        weighting=weighting,
                                        map=map)

            map = self._d.merge(map)

            if move_all and self.is_dihedral():
                moved = view.move().set_all(self._internal, value, map)
            else:
                moved = view.move().set(self._internal, value, map)

            if auto_align and (anchor is None):
                moved.align(view, map)

            self._d.molecule = moved.commit().molecule().edit()
            self._update()

            return self

        def change_size(delta, anchor=None, weighting=None,
                        auto_align=True, move_all=True, map=None):
            """Change the size of the internal being edited by this cursor by
               'delta'. This should be either an angle unit, or a float
               (in which case it is converted into a value with default
                angle units)

                delta: float or angle
                    The angle by which this internal is changed.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            from ..units import angle

            delta = angle(delta)
            view = self._d.molecule.commit()

            map = _process_move_options(view=view, anchor=anchor,
                                        weighting=weighting, map=map)

            map = self._d.merge(map)

            if move_all and self.is_dihedral():
                from . import BondID
                center_bond = BondID(self._internal.atom1(), self._internal.atom2())
                moved = view.move().change(center_bond, delta, map)
            else:
                moved = view.move().change(self._internal, delta, map)

            if auto_align and (anchor is None):
                moved.align(view, map)

            self._d.molecule = moved.commit().molecule().edit()
            self._update()

        self.size = get_size
        self.set_size = set_size
        self.change_size = change_size
        self.measure = get_size
        self.set_measure = set_size
        self.change_measure = change_size

    def _add_extra_functions(self):
        """Internal function that adds additional functions to this
           cursor depending on what type of object is being edited.
        """
        if self.is_bond():
            self._add_bond_functions()
        elif self.is_internal():
            self._add_angle_functions()

    def is_same_editor(self, other):
        """Return whether this Cursor is using the same editor to edit
           the molecule as 'other'. This returns true if the underlying
           editor for both cursors is the same, i.e. changes made by
           one cursor would be seen and be editable by the other cursor.
        """
        try:
            # Other is a Cursor
            return self._d is other._d
        except AttributeError:
            pass

        try:
            # Other is a Cursors
            return self._d is other._parent._d
        except AttributeError:
            # Other is something else (a CursorsM?)
            return False

    def bonds(self, *args, **kwargs):
        """Return cursors for all of the bonds in this
           view or, if 'id' is supplied, the bonds in this
           view that match 'id'
        """
        self._update()

        cursors = []

        bonds = self.view().bonds(*args, **kwargs)

        for bond in bonds:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule
            c._internal = bond.id()
            c._add_extra_functions()
            cursors.append(c)

        return Cursors(self, cursors, bonds)

    def angles(self, *args, **kwargs):
        """Return cursors for all of the angles in this
           view or, if 'id' is supplied, the angles in this
           view that match 'id'
        """
        self._update()

        cursors = []

        angles = self.view().angles(*args, **kwargs)

        for angle in angles:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule
            c._internal = angle.id()
            c._add_extra_functions()
            cursors.append(c)

        return Cursors(self, cursors, angles)

    def dihedrals(self, *args, **kwargs):
        """Return cursors for all of the dihedrals in this
           view or, if 'id' is supplied, the dihedrals in this
           view that match 'id'
        """
        self._update()

        cursors = []

        dihedrals = self.view().dihedrals(*args, **kwargs)

        for dihedral in dihedrals:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule
            c._internal = dihedral.id()
            c._add_extra_functions()
            cursors.append(c)

        return Cursors(self, cursors, dihedrals)

    def impropers(self, *args, **kwargs):
        """Return cursors for all of the impropers in this
           view or, if 'id' is supplied, the impropers in this
           view that match 'id'
        """
        self._update()

        cursors = []

        impropers = self.view().impropers(*args, **kwargs)

        for improper in impropers:
            c = Cursor()
            c._d = self._d
            c._view = self._d.molecule
            c._internal = improper.id()
            c._add_extra_functions()
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

        c._add_extra_functions()

        return c

    def _from_views(self, views):
        """Hidden function that allows a set of views, e.g. from
           a Selector_Atom_, SelectorBond etc, to be converted
           to a set of Cursors
        """
        cursors = []

        from ..mm import Bond, Angle, Dihedral, Improper

        for view in views:
            c = Cursor()
            c._d = self._d

            if type(view) in [Bond, Angle, Dihedral, Improper]:
                c._view = self._d.molecule
                c._internal = view.id()
            else:
                # likely an atom, residue, chain or segment
                c._view = self._d.molecule[view.index()]

            c._add_extra_functions()

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
            c._add_extra_functions()
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
            c._add_extra_functions()
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
            c._add_extra_functions()
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
            c._add_extra_functions()
            cursors.append(c)

        return Cursors(self, cursors, segments)

    def atom(self, i):
        """Return the atom in the molecule that matches the passed ID"""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._view.atom(i)
        c._add_extra_functions()
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
                c._add_extra_functions()
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
        c._add_extra_functions()
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
                c._add_extra_functions()
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
        c._add_extra_functions()
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
                c._add_extra_functions()
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
        c._add_extra_functions()
        c._update()

        return c

    def molecule(self):
        """Return the molecule"""
        self._update()

        c = Cursor()
        c._d = self._d
        c._view = self._d.molecule
        c._add_extra_functions()
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

        c._add_extra_functions()

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

        c._add_extra_functions()

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

        c._add_extra_functions()

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

        c._add_extra_functions()

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

        c._add_extra_functions()
        c._update()

        return c

    def get(self, key):
        """Return the property associated with key 'key'"""
        self._update()

        if self.is_internal():
            return self._d.connectivity.property(self._internal, key)
        else:
            return self._view.property(key)

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

    def num_frames(self):
        """Return the number of trajectory frames contained by the molecule"""
        return self._d.molecule.num_frames()

    def load_frame(self, *args, **kwargs):
        """Call the `load_frame` function on the contained view, passing
           the arguments directly. This is equivalent to calling
           `load_frame` directly on the contained view.
        """
        self._d.molecule.load_frame(*args, **kwargs)
        self._update()
        return self

    def save_frame(self, *args, **kwargs):
        """Call the `save_frame` function on the contained view, passing
           the arguments directly. This is equivalent to calling
           `save_frame` directly on the contained view.
        """
        self._d.molecule.save_frame(*args, **kwargs)
        self._update()
        return self

    def delete_frame(self, *args, **kwargs):
        """Call the `delete_frame` function on the contained view, passing
           the arguments directly. This is equivalent to calling
           `delete_frame` directly on the contained view.
        """
        self._d.molecule.delete_frame(*args, **kwargs)
        self._update()
        return self

    def translate(self, *args, map=None):
        """Translate all of the atoms operated on by this cursor
           by the passed arguments (these are converted automatically
           to a sr.maths.Vector). Use 'map' to specify the property
           map to use to find the coordinates property
        """
        from ..maths import Vector
        delta = Vector(*args)

        map = self._d.merge(map)
        view = self.commit()
        view = view.move().translate(delta, map=map).commit()

        self._d.molecule = view.molecule().edit()
        self._update()

        return self

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

        map = self._d.merge(map)

        if center is None:
            center = view.evaluate().center_of_mass(map=map)
        else:
            center = Vector(center)

        view = view.move().rotate(quaternion, center, map=map).commit()

        self._d.molecule = view.molecule().edit()
        self._update()

        return self


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

        self._add_extra_functions()

    def __call__(self, i):
        """Return the sub-view(s) of this cursor that match the index 'i'.
           Note that this will not look at the properties of the cursors.
        """
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

    def __getitem__(self, key):
        """Either find the property that matches the passed 'key' or
           look for the sub-view(s) of the cursor that match the passed key.
           This will only look for sub-view(s) if there are no matching
           properties. Use the __call__ operator to skip the property search.
        """
        if type(key) is not str:
            return self.__call__(key)

        try:
            return self.get(key)
        except Exception as e:
            property_error = e

        # We couldn't find the property, so instead
        # try to find the matching sub-view
        try:
            return self.__call__(key)
        except Exception:
            pass

        # This interface is primarily for properties,
        # so raise the missing property error
        raise property_error

    def __setitem__(self, key, value):
        for cursor in self._cursors:
            cursor.set(key, value)

    def __delitem__(self, key):
        for cursor in self._cursors:
            cursor.delete(key)

    def __contains__(self, key):
        for cursor in self._cursors:
            if key in cursor:
                return True

        return False

    def _update(self):
        """Internal function used to ensure that all
           child cursors are up to date"""
        for cursor in self._cursors:
            cursor._update()

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

    def _add_bond_functions(self):
        """Internal function that adds functions to this cursor
           that are useful when editing bonds
        """
        def get_lengths(map=None):
            """Return the lengths of all bonds being edited by these
               cursors
            """
            map = self._parent._d.merge(map)

            lengths = []

            for cursor in self._cursors:
                lengths.append(cursor.length(map=map))

            return lengths

        def set_lengths(values, anchor=None, weighting=None,
                        auto_align=True, map=None):
            """Set the lengths of the bonds being edited by this cursor to the specified
               values. Note that there should be the same number of values
               as there are cursors, and they should all be floats or
               lengths.

                values: list[float] or list[length]
                    The lengths to which to set these bonds.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            if len(values) != len(self._cursors):
                raise ValueError(
                    f"The number of length values ({len(values)}) does "
                    f"not equal the number of cursors ({len(self._cursors)}).")

            from ..units import length
            values = [length(value) for value in values]

            molecule = self._parent._d.molecule.commit()

            map = _process_move_options(view=molecule, anchor=anchor,
                                        weighting=weighting, map=map)
            map = self._parent._d.merge(map)

            moved = molecule.move()

            for value, cursor in zip(values, self._cursors):
                moved.set(cursor._internal, value, map)

            if auto_align and (anchor is None):
                moved.align(molecule, map)

            self._parent._d.molecule = moved.commit().molecule().edit()
            self._update()
            return self

        def set_length(value, anchor=None, weighting=None,
                       auto_align=True, map=None):
            """Set all bonds edited by this cursor to the supplied length.

                value: float or length
                    The length to which to set these bonds.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            values = [value for _ in range(0, len(self._cursors))]
            set_lengths(values, anchor=anchor, weighting=weighting,
                        auto_align=auto_align, map=map)
            return self

        def change_lengths(deltas, anchor=None, weighting=None,
                           auto_align=True, map=None):
            """Change the bonds being edited by this cursor by the specified
               values. Note that there should be the same number of values
               as there are cursors, and they should all be floats or
               lengths.

                deltas: list[float] or list[length]
                    The lengths by which to change these bonds.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            if len(deltas) != len(self._cursors):
                raise ValueError(
                    f"The number of length values ({len(deltas)}) does "
                    f"not equal the number of cursors ({len(self._cursors)}).")

            from ..units import length
            deltas = [length(delta) for delta in deltas]

            molecule = self._parent._d.molecule.commit()

            map = _process_move_options(view=molecule, anchor=anchor,
                                        weighting=weighting, map=map)
            map = self._parent._d.merge(map)

            moved = molecule.move()

            for delta, cursor in zip(deltas, self._cursors):
                moved.change(cursor._internal, delta, map)

            if auto_align and (anchor is None):
                moved.align(molecule, map)

            self._parent._d.molecule = moved.commit().molecule().edit()
            self._update()
            return self

        def change_length(delta, anchor=None, weighting=None,
                          auto_align=True, map=None):
            """Change all bonds edited by this cursor by the supplied length.

                delta: float or length
                    The length by which to change these bonds.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.            """
            deltas = [delta for _ in range(0, len(self._cursors))]
            change_lengths(deltas, anchor=anchor, weighting=weighting,
                           auto_align=auto_align, map=map)

        self.lengths = get_lengths
        self.set_length = set_length
        self.set_lengths = set_lengths
        self.change_length = change_length
        self.change_lengths = change_lengths

        self.measures = get_lengths
        self.set_measure = set_length
        self.set_measures = set_lengths
        self.change_measure = set_length
        self.change_measures = set_lengths

    def _add_angle_functions(self):
        """Internal function that adds functions to this cursor
           that are useful when editing angles, dihedrals or impropers
        """
        def get_sizes(map=None):
            """Return the sizes of all internals being edited by these
               cursors
            """
            map = self._parent._d.merge(map)

            sizes = []

            for cursor in self._cursors:
                sizes.append(cursor.size(map=map))

            return sizes

        def set_sizes(values, anchor=None, weighting=None,
                      auto_align=True, move_all=True, map=None):
            """Set the internals being edited by this cursor to the specified
               values. Note that there should be the same number of values
               as there are cursors, and they should all be floats or
               angles.

                values: list[float] or list[angle]
                    The angles to which to set these internals.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            if len(values) != len(self._cursors):
                raise ValueError(
                    f"The number of angle values ({len(values)}) does "
                    f"not equal the number of cursors ({len(self._cursors)}).")

            from ..units import angle
            values = [angle(value) for value in values]

            view = self._parent._d.molecule.commit()

            map = _process_move_options(view=view, anchor=anchor,
                                        weighting=weighting, map=map)

            map = self._parent._d.merge(map)

            moved = view.move()

            if move_all and self._cursors[0].is_dihedral():
                for value, cursor in zip(values, self._cursors):
                    moved.set_all(cursor._internal, value, map)
            else:
                for value, cursor in zip(values, self._cursors):
                    moved.set(cursor._internal, value, map)

            if auto_align and (anchor is None):
                moved.align(view, map)

            self._parent._d.molecule = moved.commit().molecule().edit()
            self._update()
            return self

        def set_size(value, anchor=None, weighting=None,
                     auto_align=True, move_all=True, map=None):
            """Set all internals edited by this cursor to the supplied angle.

                value: float or angle
                    The angle to which to set these internals.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            values = [value for _ in range(0, len(self._cursors))]
            set_sizes(values, anchor=anchor, weighting=weighting,
                      auto_align=auto_align, move_all=move_all, map=map)
            return self

        def change_sizes(deltas, anchor=None, weighting=None,
                         auto_align=True, move_all=True, map=None):
            """Change the internals being edited by this cursor by the specified
               values. Note that there should be the same number of values
               as there are cursors, and they should all be floats or
               angles.

                deltas: list[float] or list[angle]
                    The angles by which to change these internals.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            if len(deltas) != len(self._cursors):
                raise ValueError(
                    f"The number of angle values ({len(deltas)}) does "
                    f"not equal the number of cursors ({len(self._cursors)}).")

            from ..units import angle
            deltas = [angle(delta) for delta in deltas]

            view = self._parent._d.molecule.commit()

            map = _process_move_options(view=view, anchor=anchor,
                                        weighting=weighting, map=map)

            map = self._parent._d.merge(map)

            moved = view.move()

            if move_all and self._cursors[0].is_dihedral():
                from . import BondID
                for delta, cursor in zip(deltas, self._cursors):
                    bond = BondID(cursor._internal.atom0(),
                                  cursor._internal.atom1())
                    moved.change(bond, delta, map)
            else:
                for delta, cursor in zip(deltas, self._cursors):
                    moved.change(cursor._internal, delta, map)

            if auto_align and (anchor is None):
                moved.align(view, map)

            self._parent._d.molecule = moved.commit().molecule().edit()
            self._update()
            return self

        def change_size(delta, anchor=None, weighting=None,
                        auto_align=True, move_all=True, map=None):
            """Change all internals edited by this cursor by the supplied angle.

                delta: float or angle
                    The angle by which to change these internals.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            deltas = [delta for _ in range(0, len(self._cursors))]
            change_sizes(deltas, anchor=anchor, weighting=weighting,
                         auto_align=auto_align, move_all=move_all, map=map)
            return self

        self.sizes = get_sizes
        self.set_size = set_size
        self.set_sizes = set_sizes
        self.change_size = change_size
        self.change_sizes = change_sizes

        self.measures = get_sizes
        self.set_measure = set_size
        self.set_measures = set_sizes
        self.change_measure = set_size
        self.change_measures = set_sizes

    def _add_extra_functions(self):
        """Internal function used to add extra member functions to this
           cursor depending on the type of object being edited
        """
        if hasattr(self._cursors[0], "set_length"):
            self._add_bond_functions()
        elif hasattr(self._cursors[0], "set_size"):
            self._add_angle_functions()

    def is_same_editor(self, other):
        """Return whether this is using the same editor to edit
           the molecule as 'other'. This returns true if the underlying
           editor for both cursors is the same, i.e. changes made by
           one cursor would be seen and be editable by the other cursor.
        """
        return self._parent.is_same_editor(other)

    def get(self, key):
        """Return the property associated with key 'key'"""
        values = []

        for cursor in self._cursors:
            values.append(cursor.get(key))

        return values

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

    def bonds(self, *args, **kwargs):
        """Return cursors for all of the bonds in this
           view or, if 'id' is supplied, the bonds in this
           view that match 'id'
        """
        return self._parent._from_views(self._view.bonds(*args, **kwargs))

    def angles(self, *args, **kwargs):
        """Return cursors for all of the angles in this
           view or, if 'id' is supplied, the angles in this
           view that match 'id'
        """
        return self._parent._from_views(self._view.angles(*args, **kwargs))

    def dihedrals(self, *args, **kwargs):
        """Return cursors for all of the dihedrals in this
           view or, if 'id' is supplied, the dihedrals in this
           view that match 'id'
        """
        return self._parent._from_views(self._view.dihedrals(*args, **kwargs))

    def impropers(self, *args, **kwargs):
        """Return cursors for all of the impropers in this
           view or, if 'id' is supplied, the impropers in this
           view that match 'id'
        """
        return self._parent._from_views(self._view.impropers(*args, **kwargs))

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

    def invert(self):
        """Return the inverse view of this cursor (i.e. all views that
           are not selected - same as view.invert())"""
        return self._parent._from_views(self._view.invert())

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

    def num_frames(self):
        """Return the number of frames in the trajectory held by this molecule"""
        return self._parent.num_frames()

    def load_frame(self, *args, **kwargs):
        """Call `load_frame` with these arguments on all contained cursors"""
        self._parent.load_frame(*args, **kwargs)
        self._update()

        return self

    def save_frame(self, *args, **kwargs):
        """Call 'save_frame' with these arguments on all contained cursors"""
        self._parent.save_frame(*args, **kwargs)
        self._update()

        return self

    def delete_frame(self, *args, **kwargs):
        """Call 'delete_frame' with these arguments on all contained cursors"""
        self._parent.delete_frame(*args, **kwargs)
        self._update()

        return self

    def translate(self, *args, map=None):
        """Translate all of the atoms operated on by these cursors
           by the passed arguments (these are converted automatically
           to a sr.maths.Vector). Use 'map' to specify the property
           map to use to find the coordinates property
        """
        for cursor in self._cursors:
            cursor.translate(*args, map=map)

        return self

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

        return self


class CursorsM:
    """This class holds a list of Cursor/Cursors that operate across
       multiple molecules. This allows you to perform editing
       operations across many molecules at the same time.
    """
    def __init__(self, parent=None, map=None):
        self._parent = None
        self._cursors = []

        self._molcursors = {}

        if parent is not None:
            self._parent = parent.clone()

            for child in parent:
                child_mol = child.molecule()
                molnum = child_mol.number()

                if molnum not in self._molcursors:
                    self._molcursors[molnum] = child_mol.cursor(map=map)

                self._cursors.append(
                        self._molcursors[molnum]._from_view(child))

            self._add_extra_functions()

    def _add_bond_functions(self):
        """Internal function used to add functions for editing bonds"""

        def get_lengths(map=None):
            """Return the lengths of all the bonds edited by this cursor."""
            map = self._cursors[0]._d.merge(map)

            lengths = []

            for cursor in self._cursors:
                lengths.append(cursor.length(map=map))

            return lengths

        def set_lengths(values, anchor=None, weighting=None,
                        auto_align=True, map=None):
            """Set the lengths of the bonds being edited by this cursor to the specified
               values. Note that there should be the same number of values
               as there are cursors, and they should all be floats or
               lengths.

                values: list[float] or list[length]
                    The lengths to which to set these bonds.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            if len(values) != len(self._cursors):
                raise ValueError(
                    f"The number of length values ({len(values)}) does "
                    f"not equal the number of cursors ({len(self._cursors)}).")

            from ..units import length
            values = [length(value) for value in values]

            movers = {}
            molecules = {}
            maps = {}

            for molnum, view in self._molcursors.items():
                molecules[molnum] = view._d.molecule.commit()
                movers[molnum] = molecules[molnum].move()
                maps[molnum] = self._cursors[0]._d.merge(
                                _process_move_options(view=molecules[molnum],
                                                      anchor=anchor,
                                                      weighting=weighting,
                                                      map=map))

            for value, cursor in zip(values, self._cursors):
                molnum = cursor._d.number()
                movers[molnum].set(cursor._internal, value,
                                   maps[molnum])

            for molnum, mover in movers.items():
                if auto_align and (anchor is None):
                    mover.align(molecules[molnum], maps[molnum])

                self._molcursors[molnum]._d.molecule = mover.commit().edit()
                self._molcursors[molnum]._update()

            self._update()
            return self

        def set_length(value, anchor=None, weighting=None,
                       auto_align=True, map=None):
            """Set the lengths of all of the bonds edited by this
               cursor to the passed value.

                value: float or length
                    The length to which to set these bonds.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            values = [value for _ in self._cursors]
            set_lengths(values, anchor=anchor, weighting=weighting,
                        auto_align=auto_align, map=map)
            return self

        def change_lengths(deltas, anchor=None, weighting=None,
                           auto_align=True, map=None):
            """Change the lengths of the bonds being edited by this cursor by the specified
               deltas. Note that there should be the same number of values
               as there are cursors, and they should all be floats or
               lengths.

                deltas: list[float] or list[length]
                    The lengths by which to change these bonds.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            if len(deltas) != len(self._cursors):
                raise ValueError(
                    f"The number of length values ({len(deltas)}) does "
                    f"not equal the number of cursors ({len(self._cursors)}).")

            from ..units import length
            deltas = [length(delta) for delta in deltas]

            movers = {}
            molecules = {}
            maps = {}

            for molnum, view in self._molcursors.items():
                molecules[molnum] = view._d.molecule.commit()
                movers[molnum] = molecules[molnum].move()
                maps[molnum] = self._cursors[0]._d.merge(
                        _process_move_options(view=molecules[molnum],
                                              anchor=anchor,
                                              weighting=weighting,
                                              map=map))

            for delta, cursor in zip(deltas, self._cursors):
                molnum = cursor._d.number()
                movers[molnum].change(cursor._internal,
                                      delta, maps[molnum])

            for molnum, mover in movers.items():
                if auto_align and (anchor is None):
                    mover.align(molecules[molnum], maps[molnum])

                self._molcursors[molnum]._d.molecule = mover.commit().edit()
                self._molcursors[molnum]._update()

            self._update()
            return self

        def change_length(delta, anchor=None, weighting=None,
                          auto_align=True, map=None):
            """Change the lengths of all of the bonds edited by this
               cursor by the passed value.

                delta: float or length
                    The length by which to change these bonds.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            deltas = [delta for _ in self._cursors]
            change_lengths(deltas, anchor=anchor, weighting=weighting,
                           auto_align=auto_align, map=map)
            return self

        self.lengths = get_lengths
        self.set_length = set_length
        self.set_lengths = set_lengths
        self.change_length = change_length
        self.change_lengths = change_lengths

        self.measures = get_lengths
        self.set_measure = set_length
        self.set_measures = set_lengths
        self.change_measure = change_length
        self.change_measures = change_lengths

    def _add_angle_functions(self):
        """Internal function used to add functions for editing
           angles, dihedrals or impropers
        """

        def get_sizes(map=None):
            """Return the sizes of all the internals edited by this cursor."""
            map = self._cursors[0]._d.merge(map)

            sizes = []

            for cursor in self._cursors:
                sizes.append(cursor.size(map=map))

            return sizes

        def set_sizes(values, anchor=None, weighting=None,
                      auto_align=True, move_all=True, map=None):
            """Set the sizes of the internals being edited by this cursor to the specified
               values. Note that there should be the same number of values
               as there are cursors, and they should all be floats or
               angles.

                values: list[float] or list[angle]
                    The angles to which to set these internals.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            if len(values) != len(self._cursors):
                raise ValueError(
                    f"The number of angle values ({len(values)}) does "
                    f"not equal the number of cursors ({len(self._cursors)}).")

            from ..units import angle
            values = [angle(value) for value in values]

            molecules = {}
            movers = {}
            maps = {}

            map = self._cursors[0]._d.merge(map)

            for molnum, view in self._molcursors.items():
                molecules[molnum] = view._d.molecule.commit()
                movers[molnum] = molecules[molnum].move()
                maps[molnum] = self._cursors[0]._d.merge(
                                _process_move_options(view=molecules[molnum],
                                                      anchor=anchor,
                                                      weighting=weighting,
                                                      map=map))

            for value, cursor in zip(values, self._cursors):
                molnum = cursor._d.number()
                movers[molnum].set(cursor._internal, value,
                                   maps[molnum])

            for molnum, mover in movers.items():
                if auto_align and (anchor is None):
                    mover.align(molecules[molnum], map[molnum])

                self._molcursors[molnum]._d.molecule = mover.commit().edit()
                self._molcursors[molnum]._update()

            self._update()

            return self

        def set_size(value, anchor=None, weighting=None,
                     auto_align=True, move_all=True, map=None):
            """Set the sizes of all of the internals edited by this
               cursor to the passed value.

                value: float or angle
                    The angle to which to set these internals.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            values = [value for _ in self._cursors]
            set_sizes(values, anchor=anchor, weighting=weighting,
                      auto_align=auto_align, move_all=move_all, map=map)
            return self

        def change_sizes(deltas, anchor=None, weighting=None,
                         auto_align=True, move_all=True, map=None):
            """Change the sizes of the internals being edited by this cursor by the specified
               values. Note that there should be the same number of values
               as there are cursors, and they should all be floats or
               angles.

                deltas: list[float] or list[angle]
                    The angles by which to change these internals.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            if len(deltas) != len(self._cursors):
                raise ValueError(
                    f"The number of angle values ({len(deltas)}) does "
                    f"not equal the number of cursors ({len(self._cursors)}).")

            from ..units import angle
            deltas = [angle(delta) for delta in deltas]

            molecules = {}
            movers = {}
            maps = {}

            for molnum, view in self._molcursors.items():
                molecules[molnum] = view._d.molecule.commit()
                movers[molnum] = molecules[molnum].move()
                maps[molnum] = self._cursors[0]._d.merge(
                                _process_move_options(view=molecules[molnum],
                                                      anchor=anchor,
                                                      weighting=weighting,
                                                      map=map))

            if move_all and self._cursors[0].is_dihedral():
                from . import BondID
                for delta, cursor in zip(deltas, self._cursors):
                    molnum = cursor._d.number()
                    bond = BondID(cursor._internal.atom0(),
                                  cursor._internal.atom1())
                    movers[molnum].change(bond, delta, maps[molnum])
            else:
                for delta, cursor in zip(deltas, self._cursors):
                    molnum = cursor._d.number()
                    movers[molnum].change(cursor._internal, delta,
                                          maps[molnum])

            for molnum, mover in movers.items():
                if auto_align and (anchor is None):
                    mover.align(molecules[molnum], maps[molnum])

                self._molcursors[molnum]._d.molecule = mover.commit().edit()
                self._molcursors[molnum]._update()

            self._update()

            return self

        def change_size(value, anchor=None, weighting=None,
                        auto_align=True, move_all=True, map=None):
            """Change the sizes of all of the internals edited by this
               cursor by the passed value.

                delta: float or angle
                    The angle by which to change these internals.

                anchor: string or ID
                    The search or ID to identify atoms in the view that
                    are anchored, and which cannot be moved when the
                    internal is set.

                weighting: string
                    The weighting function used to distribute the move
                    across the atoms in the view. This is either;
                    'absolute_number', 'relative_number',
                    'absolute_mass' or 'relative_mass'. It defaults
                    to 'absolute_number'

                auto_align: bool
                    Whether or not to align the molecule against itself
                    after the move

                move_all: bool
                    Option only used for dihedrals. Whether or not to
                    move all atoms around the dihedral or just the
                    specified dihedral

                map: sire.base.PropertyMap or dict
                    Map of property keys that are passed through to
                    control which properties are used for the move,
                    plus fine-grained control of the weight function
                    or anchors.
            """
            values = [value for _ in self._cursors]
            change_sizes(values, anchor=anchor, weighting=weighting,
                         auto_align=auto_align, move_all=move_all, map=map)
            return self

        self.sizes = get_sizes
        self.set_size = set_size
        self.set_sizes = set_sizes
        self.change_size = change_size
        self.change_sizes = change_sizes

        self.measures = get_sizes
        self.set_measure = set_size
        self.set_measures = set_sizes
        self.change_measure = change_size
        self.change_measures = change_sizes

    def _add_extra_functions(self):
        """Internal function used to add extra member functions depending
           on the type of object being edited
        """
        if hasattr(self._cursors[0], "set_length"):
            self._add_bond_functions()
        elif hasattr(self._cursors[0], "set_size"):
            self._add_angle_functions()

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

        c._add_extra_functions()
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

        ret._add_extra_functions()

        return ret

    def __call__(self, i):
        try:
            # try the simplest case - the ith cursor
            idx = int(i)
            return self._cursors[idx]
        except Exception:
            pass

        return self._from_views(self._parent[i])

    def __getitem__(self, key):
        if type(key) is not str:
            return self.__call__(key)

        try:
            return self.get(key)
        except Exception as e:
            property_error = e

        try:
            return self.__call__(key)
        except Exception:
            pass

        raise property_error

    def __setitem__(self, key, value):
        for cursor in self._cursors:
            cursor.set(key, value)

    def __delitem__(self, key):
        for cursor in self._cursors:
            cursor.delete(key)

    def __contains__(self, key):
        for cursor in self._cursors:
            if key in cursor:
                return True

        return False

    def _update(self):
        """Ensure that all cursors are up to date"""
        for cursor in self._cursors:
            cursor._update()

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

    def get(self, key):
        """Return the property associated with key 'key'"""
        values = []

        for cursor in self._cursors:
            values.append(cursor.get(key))

        return values

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

    def view(self):
        """Return the view that underlies this Cursor. Note that
           this may not be updated to reflect the edits made.
           Use .commit() to get the up-to-date version of this view.
        """
        return self._parent.clone()

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

    def bonds(self, *args, **kwargs):
        """Return cursors for all of the bonds in this
           view or, if 'id' is supplied, the bonds in this
           view that match 'id'
        """
        return self._from_views(self._parent.bonds(*args, **kwargs))

    def angles(self, *args, **kwargs):
        """Return cursors for all of the angles in this
           view or, if 'id' is supplied, the angles in this
           view that match 'id'
        """
        return self._from_views(self._parent.angles(*args, **kwargs))

    def dihedrals(self, *args, **kwargs):
        """Return cursors for all of the dihedrals in this
           view or, if 'id' is supplied, the dihedrals in this
           view that match 'id'
        """
        return self._from_views(self._parent.dihedrals(*args, **kwargs))

    def impropers(self, *args, **kwargs):
        """Return cursors for all of the impropers in this
           view or, if 'id' is supplied, the impropers in this
           view that match 'id'
        """
        return self._from_views(self._parent.impropers(*args, **kwargs))

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

    def invert(self):
        """Return the inverse view of this cursor (i.e. all views that
           are not selected - same as view.invert())"""
        return self._from_views(self._parent.invert())

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

    def num_frames(self):
        """Return the number of frames in the trajectories held by
           this cursor
        """
        num = None

        for cursor in self._molcursors.values():
            if num is None:
                num = cursor.num_frames()
            else:
                n = cursor.num_frames()

                if n < num:
                    num = n

        return num

    def load_frame(self, *args, **kwargs):
        """Call the 'load_frame' function with these arguments on all
           contained cursors"""
        for cursor in self._molcursors.values():
            cursor.load_frame(*args, **kwargs)

        self._update()

        return self

    def save_frame(self, *args, **kwargs):
        """Call the 'save_frame' function with these arguments on all
           contained cursors"""
        for cursor in self._molcursors.values():
            cursor.save_frame(*args, **kwargs)

        self._update()

        return self

    def delete_frame(self, *args, **kwargs):
        """Call the 'delete_frame' function with these arguments on all
           contained cursors"""
        for cursor in self._molcursors.values():
            cursor.delete_frame(*args, **kwargs)

        self._update()

        return self

    def translate(self, *args, map=None):
        """Translate all of the atoms operated on by these cursors
           by the passed arguments (these are converted automatically
           to a sr.maths.Vector). Use 'map' to specify the property
           map to use to find the coordinates property
        """
        for cursor in self._cursors:
            cursor.translate(*args, map=map)

        return self

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

        return self
