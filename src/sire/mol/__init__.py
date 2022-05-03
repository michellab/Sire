"""
.. currentmodule:: sire.mol

"""

from ..legacy import Mol as _Mol
from .. import use_new_api as _use_new_api
_use_new_api()

from ..legacy import Base as _Base

from ..legacy.Mol import AtomName, AtomNum, AtomIdx, AtomID, \
                         ResName, ResNum, ResIdx, ResID, \
                         ChainName, ChainIdx, ChainID, \
                         SegName, SegIdx, SegID, \
                         CGName, CGIdx, CGID, \
                         MolName, MolNum, MolIdx, MolID, \
                         Atom, Selector_Atom_, SelectorM_Atom_, \
                         CutGroup, Selector_CutGroup_, SelectorM_CutGroup_, \
                         Residue, Selector_Residue_, SelectorM_Residue_, \
                         Chain, Selector_Chain_, SelectorM_Chain_, \
                         Segment, Selector_Segment_, SelectorM_Segment_, \
                         Molecule, SelectorMol, \
                         MoleculeView

# Here I will define some functions that make accessing
# things from moleculeviews more convenient
def __is_molecule_class(obj):
    mro = type(obj).mro()

    return Molecule in mro or \
            SelectorMol in mro


def __is_atom_class(obj):
    mro = type(obj).mro()

    return Atom in mro or \
            Selector_Atom_ in mro or \
            SelectorM_Atom_ in mro


def __is_residue_class(obj):
    mro = type(obj).mro()

    return Residue in mro or \
            Selector_Residue_ in mro or \
            SelectorM_Residue_ in mro


def __is_chain_class(obj):
    mro = type(obj).mro()

    return Chain in mro or \
            Selector_Chain_ in mro or \
            SelectorM_Chain_ in mro


def __is_segment_class(obj):
    mro = type(obj).mro()

    return Segment in mro or \
            Selector_Segment_ in mro or \
            SelectorM_Segment_ in mro


def __is_cutgroup_class(obj):
    mro = type(obj).mro()

    return CutGroup in mro or \
            Selector_CutGroup_ in mro or \
            SelectorM_CutGroup_ in mro


def __is_selector_class(obj):
    return obj.what().find("SireMol::Selector") != -1


def __is_list_class(obj):
    if type(obj) is list:
        return True
    else:
        try:
            return obj.what().find("::Selector") != -1
        except Exception:
            return False


def __from_select_result(obj):
    """Convert the passed SelectResult from a search into the
       most appropriate MoleculeView-derived class
    """
    views = []

    molnums = obj.mol_nums()

    if len(molnums) == 0:
        raise KeyError("Nothing matched the search.")

    elif len(molnums) == 1:
        v = obj.views(molnums[0])

        if __is_list_class(v) and len(v) == 1:
            return v[0]
        else:
            return v

    else:
        typ = obj.get_common_type()

        if typ == Molecule.typename():
            return SelectorMol(obj)
        elif typ == Atom.typename():
            return SelectorM_Atom_(obj)
        elif typ == Residue.typename():
            return SelectorM_Residue_(obj)
        elif typ == Chain.typename():
            return SelectorM_Chain_(obj)
        elif typ == Segment.typename():
            return SelectorM_Segment_(obj)
        elif typ == CutGroup.typename():
            return SelectorM_CutGroup_(obj)
        else:
            print(Atom.typename())
            print(f"Unrecognised type: {typ}. Returning as atoms.")
            return SelectorM_Atom_(obj)


def __fixed__getitem__(obj, key):
    if type(key) is int:
        if __is_selector_class(obj):
            return obj.__orig__getitem__(key)
        elif __is_chain_class(obj):
            return obj.residue(key)
        else:
            return obj.atom(key)
    elif type(key) is str:
        # is this a search object - if so, then return whatever is
        # most relevant from the search
        try:
            return __from_select_result(obj.search(key))
        except SyntaxError:
            pass
    elif AtomID in type(key).mro():
        return obj.atoms(key, auto_reduce=True)
    elif ResID in type(key).mro():
        return obj.residues(key, auto_reduce=True)
    elif ChainID in type(key).mro():
        return obj.chains(key, auto_reduce=True)
    elif SegID in type(key).mro():
        return obj.segments(key, auto_reduce=True)

    if __is_selector_class(obj):
        return obj.__orig__getitem__(key)
    elif __is_chain_class(obj):
        return obj.residues(key, auto_reduce=True)
    else:
        return obj.atoms(key, auto_reduce=True)


def __fixed__atoms__(obj, idx=None, auto_reduce=False):
    if idx is None:
        result = obj.__orig__atoms()
    elif type(idx) is range:
        result = obj.__orig__atoms(list(idx))
    else:
        result = obj.__orig__atoms(idx)

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__residues__(obj, idx=None, auto_reduce=False):
    if idx is None:
        result = obj.__orig__residues()
    elif type(idx) is range:
        result = obj.__orig__residues(list(idx))
    else:
        result = obj.__orig__residues(idx)

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__chains__(obj, idx=None, auto_reduce=False):
    if idx is None:
        result = obj.__orig__chains()
    elif type(idx) is range:
        result = obj.__orig__chains(list(idx))
    else:
        result = obj.__orig__chains(idx)

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__segments__(obj, idx=None, auto_reduce=False):
    if idx is None:
        result = obj.__orig__segments()
    elif type(idx) is range:
        result = obj.__orig__segments(list(idx))
    else:
        result = obj.__orig__segments(idx)

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fix_getitem(C):
    if not hasattr(C, "__orig__getitem__"):
        C.__orig__getitem__ = C.__getitem__

    if not hasattr(C, "__orig__atoms"):
        C.__orig__atoms = C.atoms

    if not hasattr(C, "__orig__residues"):
        C.__orig__residues = C.residues

    if not hasattr(C, "__orig__chains"):
        C.__orig__chains = C.chains

    if not hasattr(C, "__orig__segments"):
        C.__orig__segments = C.segments

    C.__getitem__ = __fixed__getitem__
    C.atoms = __fixed__atoms__
    C.residues = __fixed__residues__
    C.chains = __fixed__chains__
    C.segments = __fixed__segments__


for C in [Atom, CutGroup, Residue, Chain, Segment, Molecule,
          Selector_Atom_, Selector_Residue_,
          Selector_Chain_, Selector_Segment_,
          Selector_CutGroup_,
          SelectorMol, SelectorM_Atom_, SelectorM_Residue_,
          SelectorM_Chain_, SelectorM_Segment_,
          SelectorM_CutGroup_]:
    __fix_getitem(C)


MoleculeView.coordinates = lambda x : x.property("coordinates")
Atom.element = lambda x : x.property("element")
Atom.x = lambda x : x.property("coordinates").x()
Atom.y = lambda x : x.property("coordinates").y()
Atom.z = lambda x : x.property("coordinates").z()

def _add_evals(obj):
    obj.mass = lambda x : x.evaluate().mass()
    obj.charge = lambda x : x.evaluate().charge()

for C in [MoleculeView, SelectorMol, SelectorM_Atom_,
          SelectorM_Residue_, SelectorM_Chain_,
          SelectorM_CutGroup_, SelectorM_Segment_]:
    _add_evals(C)

def _get_atom_mass(x):
    if x.has_property("mass"):
        return x.property("mass")
    elif x.has_property("element"):
        return x.property("element").mass()
    else:
        return 0

Atom.mass = _get_atom_mass

def _get_atom_charge(x):
    if x.has_property("charge"):
        return x.property("charge")
    elif x.has_property("formal_charge"):
        return x.property("formal_charge")
    else:
        return 0

Atom.charge = _get_atom_charge

Molecule.connectivity = lambda x : x.property("connectivity")


## Fixing the property functions
def __get_property__(molview, key):
    property_type = molview.property_type(key).replace("::","_")

    return getattr(molview, "_get_property_%s" % property_type)(key)


def __get_metadata__(molview, *args):

    if len(args) == 1:
        metakey = args[0]
        property_type = molview.metadata_type(metakey).replace("::","_")
        return getattr(molview, "_get_metadata_%s" % property_type)(metakey)

    elif len(args) == 2:
         (key, metakey) = args
         property_type = molview.metadata_type(key, metakey).replace("::","_")
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
            return ("QVariant", _Qt.QVariant(obj))


def _match_to_type(typename, property):
    """Match the passed type of the property to the typename
       of the AtomProperty, CGProperty etc that is used to
       hold that type.

       This is useful to, e.g. allow a AtomStringArrayProperty
       to be set on a per-atom basis from DoubleArrayProperty
       values.
    """
    if typename.endswith("StringArrayProperty"):
        return _Base.StringArrayProperty(property)
    elif typename.endswith("DoubleArrayProperty"):
        return _Base.DoubleArrayProperty(property)
    elif typename.endswith("IntegerArrayProperty"):
        return _Base.IntegerArrayProperty(property)
    elif typename.endswith("PropertyList"):
        return _Base.PropertyList(property)
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
        if e.__class__.__name__ == "ArgumentError" or e.__class__.__name__ == "AttributeError":
            return _set_property(molview, key, _Base.wrap(property))
        else:
            raise e


def __set_bond_property__(connectivity, bond, key, property):
    try:
        return connectivity.__setProperty__(bond, key, property)
    except Exception as e:
        if e.__class__.__name__ == "ArgumentError":
            return connectivity.__setProperty__(bond, key,
                                                _Base.wrap(property))
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


_Mol.Atom.property = __get_property__
_Mol.AtomEditorBase.set_property = __set_property__
_Mol.Atom.metadata = __get_metadata__
_Mol.AtomEditorBase.set_metadata = __set_metadata__

_Mol.CutGroup.property = __get_property__
_Mol.CGEditorBase.set_property = __set_property__
_Mol.CutGroup.metadata = __get_metadata__
_Mol.CGEditorBase.set_metadata = __set_metadata__

_Mol.Residue.property = __get_property__
_Mol.ResEditorBase.set_property = __set_property__
_Mol.Residue.metadata = __get_metadata__
_Mol.ResEditorBase.set_metadata = __set_metadata__

_Mol.Chain.property = __get_property__
_Mol.ChainEditorBase.set_property = __set_property__
_Mol.Chain.metadata = __get_metadata__
_Mol.ChainEditorBase.set_metadata = __set_metadata__

_Mol.Segment.property = __get_property__
_Mol.SegEditorBase.set_property = __set_property__
_Mol.Segment.metadata = __get_metadata__
_Mol.SegEditorBase.set_metadata = __set_metadata__

_Mol.ConnectivityEditor.__set_property__ = _Mol.ConnectivityEditor.set_property
_Mol.ConnectivityEditor.set_property = __set_bond_property__

_Mol.MolEditor.__set_property__ = _Mol.MolEditor.set_property
_Mol.MolEditor.set_property = _Base.__set_property__


#### Here are some extra classes / functions defined as part of the
#### public API

from ._cursor import *

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
