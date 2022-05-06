"""
.. currentmodule:: sire.mol

"""

from ..legacy import Mol as _Mol
from .. import use_new_api as _use_new_api
_use_new_api()

from ..legacy import Search as _Search
_Search.install_search_parser()

from ..legacy import Base as _Base

from ..legacy.Mol import AtomName, AtomNum, AtomIdx, AtomID, \
                         ResName, ResNum, ResIdx, ResID, \
                         ChainName, ChainIdx, ChainID, \
                         SegName, SegIdx, SegID, \
                         CGName, CGIdx, CGID, \
                         MolName, MolNum, MolIdx, MolID, \
                         BondID, AngleID, DihedralID, ImproperID, \
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


def __is_bond_class(obj):
    mro = type(obj).mro()

    from sire.mm import Bond, SelectorBond

    return Bond in mro or \
        SelectorBond in mro


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
    t = obj.what()
    return t.find("SireMol::Selector") != -1 or \
                t.find("SireMM::Selector") != -1


def __is_list_class(obj):
    if type(obj) is list:
        return True
    else:
        try:
            return obj.what().find("::Selector") != -1
        except Exception:
            return False


def __fix_obj(obj):
    """This is needed because MolViewPtr objects that hold Selector_T_ types
       do not convert properly.
    """
    w = obj.what()

    if w == Selector_Atom_.typename():
        return obj.atoms()
    elif w == Selector_Residue_.typename():
        return obj.residues()
    elif w == Selector_Chain_.typename():
        return obj.chains()
    elif w == Selector_Segment_.typename():
        return obj.segments()
    elif w == Selector_CutGroup_.typename():
        return obj.cutgroups()
    else:
        return obj


def __from_select_result(obj):
    """Convert the passed SelectResult from a search into the
       most appropriate MoleculeView-derived class
    """
    if obj.list_count() == 0:
        raise KeyError("Nothing matched the search.")

    typ = obj.get_common_type()

    if obj.list_count() == 1:
        obj = __fix_obj(obj.list_at(0))

        if obj.what() != typ:
            if typ == Molecule.typename():
                return obj.molecule()
            elif typ == Segment.typename():
                return obj.segments(auto_reduce=True)
            elif typ == Chain.typename():
                return obj.chains(auto_reduce=True)
            elif typ == Residue.typename():
                return obj.residues(auto_reduce=True)
            elif typ == Atom.typename():
                return obj.atoms(auto_reduce=True)

        return obj

    from ..mm import SelectorBond

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
    elif typ == SelectorBond.typename():
        # Eventually want a SelectorMBond
        return [obj.list_at(i) for i in range(0, obj.list_count())]
    else:
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
    elif BondID in type(key).mro():
        return obj.bonds(key, auto_reduce=True)

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


def __fixed__bonds__(obj, idx=None, idx1=None, auto_reduce=False):
    if idx is None and idx1 is not None:
        idx = idx1
        idx1 = None

    if idx is None:
        from ..mm import SelectorBond
        result = SelectorBond(obj)
    elif idx1 is None:
        from ..mm import SelectorBond
        result = SelectorBond(obj.atoms(idx))
    else:
        from ..mm import SelectorBond
        result = SelectorBond(obj.atoms(idx), obj.atoms(idx1))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__bond__(obj, idx=None, idx1=None):
    bonds = __fixed__bonds__(obj, idx, idx1, auto_reduce=False)

    if len(bonds) == 0:
        raise KeyError("There is no matching bond in this view.")
    elif len(bonds) > 1:
        raise KeyError(
            f"More than one bond matchs. Number of matches is {len(bonds)}.")

    return bonds[0]


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


def __fixed__molecules__(obj, idx=None, auto_reduce=False):
    if idx is None:
        result = obj.__orig__molecules()
    elif type(idx) is range:
        result = obj.__orig__molecules(list(idx))
    else:
        result = obj.__orig__molecules(idx)

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

    C.count = C.__len__
    C.size = C.__len__

    if hasattr(C, "molecules"):
        if not hasattr(C, "__orig__molecules"):
            C.__orig__molecules = C.molecules

        C.molecules = __fixed__molecules__
    else:
        # currently don't have this on multi-molecule containers
        C.bonds = __fixed__bonds__
        C.bond = __fixed__bond__


Residue.__len__ = Residue.nAtoms
Chain.__len__ = Chain.nResidues
Segment.__len__ = Segment.nAtoms
CutGroup.__len__ = CutGroup.nAtoms
Molecule.__len__ = Molecule.nAtoms


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
