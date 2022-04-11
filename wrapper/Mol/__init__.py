"""
.. currentmodule:: Sire.Mol

Classes
=======

.. autosummary::
    :toctree: generated/

    Atom
    Cursor
    Molecule
    Residue
    Segment


Functions
=========

.. autosummary::
    :toctree: generated/

    get_molview

"""

from calendar import c
from importlib.util import resolve_name
from typing import ChainMap
import Sire.Maths
import Sire.Base
import Sire.ID
import Sire.Qt
import Sire.CAS
import Sire.Vol
import Sire.Units

from Sire.Mol._Mol import *


def __get_property__(molview, key):
    property_type = molview.propertyType(key).replace("::","_")

    return getattr(molview, "_get_property_%s" % property_type)(key)


def __get_metadata__(molview, *args):

    if len(args) == 1:
        metakey = args[0]
        property_type = molview.metadataType(metakey).replace("::","_")
        return getattr(molview, "_get_metadata_%s" % property_type)(metakey)

    elif len(args) == 2:
         (key, metakey) = args
         property_type = molview.metadataType(key, metakey).replace("::","_")
         return getattr(molview, "_get_metadata_%s" % property_type)(key, metakey)

    else:
        raise AttributeError( "Only molview.metadata(metakey) or molview.metadata(key, metakey) are valid!" )

_typename_mapping = {"SireMol_Velocity3D" : "SireMaths_Vector3D_SireUnits_Dimension_Velocity_"}


def __get_typename__(obj):
    try:
        typename = obj.typeName().replace("::","_")
        return (_typename_mapping.get(typename, typename), obj)
    except:
        if isinstance(obj, float):
            return ("double", obj)
        elif isinstance(obj, int):
            return ("qint64", obj)
        elif isinstance(obj, str):
            return ("QString", obj)
        else:
            return ("QVariant", Sire.Qt.QVariant(obj))


def _match_to_type(typename, property):
    """Match the passed type of the property to the typename
       of the AtomProperty, CGProperty etc that is used to
       hold that type.

       This is useful to, e.g. allow a AtomStringArrayProperty
       to be set on a per-atom basis from DoubleArrayProperty
       values.
    """
    if typename.endswith("StringArrayProperty"):
        return Sire.Base.StringArrayProperty(property)
    elif typename.endswith("DoubleArrayProperty"):
        return Sire.Base.DoubleArrayProperty(property)
    elif typename.endswith("IntegerArrayProperty"):
        return Sire.Base.IntegerArrayProperty(property)
    elif typename.endswith("PropertyList"):
        return Sire.Base.PropertyList(property)
    else:
        return property


def _set_property(molview, key, property):
    if molview.hasProperty(key):
        # get the type of the existing property
        typename = molview.propertyType(key)
        property = _match_to_type(typename, property)

    (typename, property) = __get_typename__(property)

    return getattr(molview, "_set_property_%s" % typename)(key, property)


def __set_property__(molview, key, property):
    try:
        return _set_property(molview, key, property)
    except Exception as e:
        if e.__class__.__name__ == "ArgumentError":
            return _set_property(molview, key, Sire.Base.wrap(property))
        else:
            raise e


def __set_bond_property__(connectivity, bond, key, property):
    try:
        return connectivity.__setProperty__(bond, key, property)
    except Exception as e:
        if e.__class__.__name__ == "ArgumentError":
            return connectivity.__setProperty__(bond, key,
                                                Sire.Base.wrap(property))
        else:
            raise e


def __set_metadata__(molview, *args):

    if len(args) == 2:
        metakey = args[0]
        property = args[1]

        (typename, property) = __get_typename__(property)

        return getattr(molview, "_set_metadata_%s" % typename)(metakey, property)

    elif len(args) == 3:
         (key, metakey, property) = args

         (typename, property) = __get_typename__(property)

         return getattr(molview, "_set_metadata_%s" % typename)(key, metakey, property)

    else:
        raise AttributeError( "Only molview.setMetadata(metakey, property) " + \
                              "or molview.setMetadata(key, metakey, property) are valid!" )


Atom.property = __get_property__
AtomEditorBase.setProperty = __set_property__
Atom.metadata = __get_metadata__
AtomEditorBase.setMetadata = __set_metadata__

CutGroup.property = __get_property__
CGEditorBase.setProperty = __set_property__
CutGroup.metadata = __get_metadata__
CGEditorBase.setMetadata = __set_metadata__

Residue.property = __get_property__
ResEditorBase.setProperty = __set_property__
Residue.metadata = __get_metadata__
ResEditorBase.setMetadata = __set_metadata__

Chain.property = __get_property__
ChainEditorBase.setProperty = __set_property__
Chain.metadata = __get_metadata__
ChainEditorBase.setMetadata = __set_metadata__

Segment.property = __get_property__
SegEditorBase.setProperty = __set_property__
Segment.metadata = __get_metadata__
SegEditorBase.setMetadata = __set_metadata__

ConnectivityEditor.__setProperty__ = ConnectivityEditor.setProperty
ConnectivityEditor.setProperty = __set_bond_property__

MolEditor.__setProperty__ = MolEditor.setProperty
MolEditor.setProperty = Sire.Base.__set_property__


def get_molview(mol):
    """Convert the passed molecule into the most appropriate view,
       e.g. a PartialMolecule containing all atoms will be returned
       as a Molecule, a PartialMolecule with a single atom will be
       returned as an Atom etc."""
    return mol


Atom._indexType = AtomIdx
Atom._nameType = AtomName
Atom._numberType = AtomNum

CutGroup._indexType = CGIdx
CutGroup._nameType = CGName
CutGroup._childType = Atom

Residue._indexType = ResIdx
Residue._nameType = ResName
Residue._numberType = ResNum
Residue._childType = Atom

Chain._indexType = ChainIdx
Chain._nameType = ChainName
Chain._childType = Residue

Segment._indexType = SegIdx
Segment._nameType = SegName
Segment._childType = Atom

Molecule._indexType = MolIdx
Molecule._nameType = MolName
Molecule._numberType = MolNum
Molecule._childType = Atom


def __getitem__view__(view, i):
    """Return the object at the specified index."""
    try:
        idx = view._childType._indexType(int(i))
    except Exception:
        idx = None

    if idx is not None:
        return view.__orig_getitem__(idx)

    try:
        return view.__orig_getitem__(i)
    except Exception:
        pass

    idx = view._childType._nameType(str(i))

    try:
        return view.__orig_getitem__(i)
    except Exception:
        pass

    idx = str(i)

    if idx.find(":") != -1:
        (name, number) = idx.split(":")[0:2]
        Name = view._childType._nameType
        Number = view._childType._numberType
        idx = Name(name) + Number(number)
        return view.__orig_getitem__(idx)

    # nothing worked - just use the original 'i' so that we
    # get the right error message
    if type(i) is str:
        i = view._childType._nameType(i)

    return view.__orig_getitem__(i)


def _fix_atom(view, i):
    try:
        idx = int(i)
    except Exception:
        idx = None

    if idx is not None:
        return view._orig_atom(AtomIdx(idx))

    try:
        return view._orig_atom(AtomName(str(i)))
    except Exception:
        pass

    if type(i) is str:
        i = AtomName(i)

    return view._orig_atom(i)


def _fix_residue(view, i):
    try:
        idx = int(i)
    except Exception:
        idx = None

    if idx is not None:
        return view._orig_residue(ResIdx(idx))

    try:
        return view._orig_residue(ResName(str(i)))
    except Exception:
        pass

    if type(i) is str:
        i = ResName(i)

    return view._orig_residue(i)


def _fix_getitem(view):

    if not hasattr(view, "__orig_getitem__"):
        view.__orig_getitem__ = view.__getitem__

    view.__getitem__ = __getitem__view__

    if hasattr(view, "atom"):
        if not hasattr(view, "_orig_atom"):
            view._orig_atom = view.atom

        view.atom = _fix_atom

    if hasattr(view, "residue"):
        if not hasattr(view, "_orig_residue"):
            view._orig_residue = view.residue

        view.residue = _fix_residue


_fix_getitem(Atom)
_fix_getitem(CutGroup)
_fix_getitem(Residue)
_fix_getitem(Chain)
_fix_getitem(Segment)
_fix_getitem(MolEditor)
_fix_getitem(Molecule)

class IncompatibleError(Exception):
    pass


# Python automatically converts ViewsOfMol into a list. Need
# to add a Molecule.join function to convert a list back into
# a single molecule
@staticmethod
def _molecule_join( views ):
    """Join the passed views of a molecule into a single molecule"""
    if len(views) == 0:
        return Molecule()
    elif len(views) == 1:
        return views[0]
    else:
        atoms = views[0].selection()
        molnum = views[0].molecule().number()

        for i in range(1,len(views)):
            if views[i].molecule().number() != molnum:
                raise IncompatibleError( \
                    "Cannot join different molecules together! %s vs. %s" \
                        % (molnum, views[i].number()) )

            atoms = atoms.unite(views[i].selection())

        return get_molview( PartialMolecule(views[i], atoms) )

Molecule.join = _molecule_join

##########
########## CLUDGY WORKAROUND
##########

# python wrappers can't distinguish between AtomProperty
# typedefs, and the full template classes,
#  (e.g. AtomLJs.array gives an lvalue error
#   as it wants a AtomProperty<LJParameter>)
#
#  I can fix this by accessing the arrays first
#  via the following code

__p = Sire.Base.Properties()

def _pvt_property_cludge_fix(C):
   __p.setProperty("c", C())
   t = __p.property("c").array()

__props = [ AtomCharges, AtomElements,
            AtomStringArrayProperty,
            AtomPropertyList,
            AtomDoubleArrayProperty,
            AtomIntegerArrayProperty ]

for __prop in __props:
    _pvt_property_cludge_fix(__prop)

##########
########## END OF CLUDGY WORKAROUND
##########

#### Here are some extra classes / functions defined as part of the
#### public API

class Cursor:
    """This class provides a cursor that can be used to navigate through
       and edit the properties of Molecules. This makes the whole
       getting and setting of properties more pythonic in writing
       style, while also saving some typing.
    """
    def __init__(self, molecule: MoleculeView = None):
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
        p = Sire.Base.Properties()

        for key in self.keys():
            p[key] = self.__getitem__(key)

        return p


def _cursor(view):
    """Return a Cursor that can be used to edit the properties
       of this view
    """
    return Cursor(view)


Atom.cursor = _cursor
Residue.cursor = _cursor
Chain.cursor = _cursor
Segment.cursor = _cursor
Molecule.cursor = _cursor
