"""
.. currentmodule:: sire.mol

"""

from multiprocessing import reduction
from tkinter import E
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
                         BondID, AngleID, DihedralID, ImproperID, \
                         Atom, Selector_Atom_, SelectorM_Atom_, \
                         CutGroup, Selector_CutGroup_, SelectorM_CutGroup_, \
                         Residue, Selector_Residue_, SelectorM_Residue_, \
                         Chain, Selector_Chain_, SelectorM_Chain_, \
                         Segment, Selector_Segment_, SelectorM_Segment_, \
                         Molecule, SelectorMol, \
                         MoleculeView, Select, \
                         BondType, Stereoscopy, \
                         AtomCoords, AtomMasses, AtomCharges

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


def __is_angle_class(obj):
    mro = type(obj).mro()

    from sire.mm import Angle, SelectorAngle

    return Angle in mro or \
        SelectorAngle in mro


def __is_dihedral_class(obj):
    mro = type(obj).mro()

    from sire.mm import Dihedral, SelectorDihedral

    return Dihedral in mro or \
        SelectorDihedral in mro


def __is_improper_class(obj):
    mro = type(obj).mro()

    from sire.mm import Improper, SelectorImproper

    return Improper in mro or \
        SelectorImproper in mro


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
    try:
        t = obj.what()
        return t.find("SireMol::Selector") != -1 or \
                    t.find("SireMM::Selector") != -1
    except Exception:
        return False


def __is_internal_class(obj):
    from ..mm import Bond, Angle, Dihedral, Improper
    return type(obj) in [Bond, Angle, Dihedral, Improper]


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

    if w == Molecule.typename():
        return obj
    elif w == Selector_Atom_.typename():
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
    if hasattr(obj, "listCount") and not hasattr(obj, "list_count"):
        # Sometimes the SelectResult hasn't been converted, i.e. because
        # it has come from an old api or mixed version of Sire, or
        # BioSimSpace has done something weird...
        raise SystemError(
            "Something has gone wrong with sire. Despite being loaded "
            "with the new or mixed API, it is being passed the object "
            f"'{obj}' of type {type(obj)} which only has the old API active. "
            "Has Sire been loaded with support for old module names?")

    if obj.list_count() == 0:
        raise KeyError("Nothing matched the search.")

    typ = obj.get_common_type()

    if obj.list_count() == 1:
        obj = __fix_obj(obj.list_at(0))

        if obj.what() in ["SireMM::SelectorBond",
                          "SireMM::SelectorAngle",
                          "SireMM::SelectorDihedral",
                          "SireMM::SelectorImproper"]:
            if obj.count() == 1:
                obj = obj[0]
        elif obj.what() != typ:
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
        from ..mm import SelectorBond, SelectorMBond
        if SelectorBond in type(obj.list_at(0)).mro():
            return SelectorMBond(obj)

        from ..mm import SelectorAngle, SelectorMAngle
        if SelectorAngle in type(obj.list_at(0)).mro():
            return SelectorMAngle(obj)

        from ..mm import SelectorDihedral, SelectorMDihedral
        if SelectorDihedral in type(obj.list_at(0)).mro():
            return SelectorMDihedral(obj)

        from ..mm import SelectorImproper, SelectorMImproper
        if SelectorImproper in type(obj.list_at(0)).mro():
            return SelectorMImproper(obj)

        # return this as a raw list
        return obj.to_list()


def __select_call__(obj, molecules, map={}):
    """Search for the desired objects in the passed molecules,
       optionally passing in a property map to identify the properties
    """
    from ..system import System
    if type(molecules) is System:
        molecules = molecules._system

    return __from_select_result(obj.__orig_call__(molecules, map))


if not hasattr(Select, "__orig_call__"):
    Select.__orig_call__ = Select.__call__
    Select.__call__ = __select_call__


def __fixed__getitem__(obj, key):
    if type(key) is int:
        if __is_selector_class(obj):
            return obj.__orig__getitem__(key)
        elif __is_chain_class(obj):
            return obj.residue(key)
        elif __is_internal_class(obj):
            return obj.atom(obj.id()[key])
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
    elif AngleID in type(key).mro():
        return obj.angles(key, auto_reduce=True)
    elif DihedralID in type(key).mro():
        return obj.dihedrals(key, auto_reduce=True)
    elif ImproperID in type(key).mro():
        return obj.impropers(key, auto_reduce=True)

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


def __fixed__bonds__(obj, idx=None, idx1=None, auto_reduce=False, map=None):
    if map is None:
        from ..base import PropertyMap
        map = PropertyMap()

    if idx is None and idx1 is not None:
        idx = idx1
        idx1 = None

    if hasattr(obj, "molecules"):
        # this is a multi-molecule container
        from ..mm import SelectorMBond
        C = SelectorMBond
        def _fromBondID(obj, bondid):
            return SelectorMBond(obj.to_select_result(), bondid, map=map)
    else:
        from ..mm import SelectorBond
        C = SelectorBond
        def _fromBondID(obj, bondid):
            return SelectorBond(obj, bondid, map=map)

    if idx is None:
        result = C(obj)
    elif idx1 is None:
        if BondID in type(idx).mro():
            result = _fromBondID(obj, idx)
        else:
            result = C(obj.atoms(idx))
    else:
        result = C(obj.atoms(idx), obj.atoms(idx1))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__angles__(obj, idx=None, idx1=None, idx2=None, auto_reduce=False,
                      map=None):
    if map is None:
        from ..base import PropertyMap
        map = PropertyMap()

    if idx1 is None and idx2 is not None:
        idx1 = idx2
        idx2 = None

    if idx is None and idx1 is not None:
        idx = idx1
        idx1 = None

    if hasattr(obj, "molecules"):
        # this is a multi-molecule container
        from ..mm import SelectorMAngle
        C = SelectorMAngle
        def _fromAngleID(obj, angid):
            return SelectorMAngle(obj.to_select_result(), angid, map=map)
    else:
        from ..mm import SelectorAngle
        C = SelectorAngle
        def _fromAngleID(obj, angid):
            return SelectorAngle(obj, angid, map=map)

    if idx is None:
        result = C(obj)
    elif idx1 is None:
        if AngleID in type(idx).mro():
            result = _fromAngleID(obj, idx)
        else:
            result = C(obj.atoms(idx))
    elif idx2 is None:
        result = C(obj.atoms(idx), obj.atoms(idx1))
    else:
        result = C(obj.atoms(idx), obj.atoms(idx1), obj.atoms(idx2))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__dihedrals__(obj, idx=None, idx1=None,
                         idx2=None, idx3=None, auto_reduce=False, map=None):
    if map is None:
        from ..base import PropertyMap
        map = PropertyMap()

    if idx2 is None and idx3 is not None:
        idx2 = idx3
        idx3 = None

    if idx1 is None and idx2 is not None:
        idx1 = idx2
        idx2 = None

    if idx is None and idx1 is not None:
        idx = idx1
        idx1 = None

    if hasattr(obj, "molecules"):
        # this is a multi-molecule container
        from ..mm import SelectorMDihedral
        C = SelectorMDihedral
        def _fromDihedralID(obj, dihid):
            return SelectorMDihedral(obj.to_select_result(), dihid, map=map)
    else:
        from ..mm import SelectorDihedral
        C = SelectorDihedral
        def _fromDihedralID(obj, dihid):
            return SelectorDihedral(obj, dihid, map=map)

    if idx is None:
        result = C(obj)
    elif idx1 is None:
        if DihedralID in type(idx).mro():
            result = _fromDihedralID(obj, idx)
        else:
            result = C(obj.atoms(idx))
    elif idx2 is None:
        result = C(obj.atoms(idx), obj.atoms(idx1))
    elif idx3 is None:
        result = C(obj.atoms(idx), obj.atoms(idx1), obj.atoms(idx2))
    else:
        result = C(obj.atoms(idx), obj.atoms(idx1),
                   obj.atoms(idx2), obj.atoms(idx3))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__impropers__(obj, idx=None, idx1=None,
                         idx2=None, idx3=None, auto_reduce=False, map=None):
    if map is None:
        from ..base import PropertyMap
        map = PropertyMap()

    if idx2 is None and idx3 is not None:
        idx2 = idx3
        idx3 = None

    if idx1 is None and idx2 is not None:
        idx1 = idx2
        idx2 = None

    if idx is None and idx1 is not None:
        idx = idx1
        idx1 = None

    if hasattr(obj, "molecules"):
        # this is a multi-molecule container
        from ..mm import SelectorMImproper
        C = SelectorMImproper
        def _fromImproperID(obj, impid):
            return SelectorMImproper(obj.to_select_result(), impid, map=map)
    else:
        from ..mm import SelectorImproper
        C = SelectorImproper
        def _fromImproperID(obj, impid):
            return SelectorImproper(obj, impid, map=map)

    if idx is None:
        result = C(obj)
    elif idx1 is None:
        if ImproperID in type(idx).mro():
            result = _fromImproperID(obj, idx)
        else:
            result = C(obj.atoms(idx))
    elif idx2 is None:
        result = C(obj.atoms(idx), obj.atoms(idx1))
    elif idx3 is None:
        result = C(obj.atoms(idx), obj.atoms(idx1), obj.atoms(idx2))
    else:
        result = C(obj.atoms(idx), obj.atoms(idx1),
                   obj.atoms(idx2), obj.atoms(idx3))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__bond__(obj, idx=None, idx1=None, map=None):
    bonds = __fixed__bonds__(obj, idx, idx1, auto_reduce=False, map=map)

    if len(bonds) == 0:
        raise KeyError("There is no matching bond in this view.")
    elif len(bonds) > 1:
        raise KeyError(
            f"More than one bond matches. Number of matches is {len(bonds)}.")

    return bonds[0]


def __fixed__angle__(obj, idx=None, idx1=None, idx2=None, map=None):
    angles = __fixed__angles__(obj, idx, idx1, idx2,
                               auto_reduce=False, map=map)

    if len(angles) == 0:
        raise KeyError("There is no matching angle in this view.")
    elif len(angles) > 1:
        raise KeyError(
            f"More than one angle matches. Number of matches is {len(angles)}.")

    return angles[0]


def __fixed__dihedral__(obj, idx=None, idx1=None,
                        idx2=None, idx3=None, map=None):
    dihedrals = __fixed__dihedrals__(obj, idx, idx1, idx2, idx3,
                                     auto_reduce=False, map=map)

    if len(dihedrals) == 0:
        raise KeyError("There is no matching dihedral in this view.")
    elif len(dihedrals) > 1:
        raise KeyError(
            f"More than one dihedral matches. Number of matches is {len(dihedrals)}.")

    return dihedrals[0]


def __fixed__improper__(obj, idx=None, idx1=None,
                        idx2=None, idx3=None, map=None):
    impropers = __fixed__impropers__(obj, idx, idx1, idx2, idx3,
                                     auto_reduce=False, map=map)

    if len(impropers) == 0:
        raise KeyError("There is no matching improper in this view.")
    elif len(impropers) > 1:
        raise KeyError(
            f"More than one improper matches. Number of matches is {len(impropers)}.")

    return impropers[0]


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

    if hasattr(C, "measure"):
        # make sure we use the right `size` function
        C.size = C.measure
    else:
        C.size = C.__len__

    if hasattr(C, "molecules"):
        if not hasattr(C, "__orig__molecules"):
            C.__orig__molecules = C.molecules

        C.molecules = __fixed__molecules__

    C.bonds = __fixed__bonds__
    C.bond = __fixed__bond__
    C.angles = __fixed__angles__
    C.angle = __fixed__angle__
    C.dihedrals = __fixed__dihedrals__
    C.dihedral = __fixed__dihedral__
    C.impropers = __fixed__impropers__
    C.improper = __fixed__improper__

try:
    Residue.__len__ = Residue.nAtoms
    Chain.__len__ = Chain.nResidues
    Segment.__len__ = Segment.nAtoms
    CutGroup.__len__ = CutGroup.nAtoms
    Molecule.__len__ = Molecule.nAtoms
except AttributeError:
    Residue.__len__ = Residue.num_atoms
    Chain.__len__ = Chain.num_residues
    Segment.__len__ = Segment.num_atoms
    CutGroup.__len__ = CutGroup.num_atoms
    Molecule.__len__ = Molecule.num_atoms

for C in [Atom, CutGroup, Residue, Chain, Segment, Molecule,
          Selector_Atom_, Selector_Residue_,
          Selector_Chain_, Selector_Segment_,
          Selector_CutGroup_,
          SelectorMol, SelectorM_Atom_, SelectorM_Residue_,
          SelectorM_Chain_, SelectorM_Segment_,
          SelectorM_CutGroup_]:
    __fix_getitem(C)

Atom.element = lambda x : x.property("element")
Atom.lj = lambda x : x.property("LJ")
Atom.coordinates = lambda x : x.property("coordinates")
Atom.coords = Atom.coordinates
Atom.x = lambda x : x.coordinates().x()
Atom.y = lambda x : x.coordinates().y()
Atom.z = lambda x : x.coordinates().z()


def __atomcoords__str__(obj):
    n = len(obj)

    if n == 0:
        return "AtomCoords::empty"

    parts = []

    if n <= 10:
        for i in range(0, 10):
            parts.append(f"{i}: {obj[i]}")
    else:
        for i in range(0, 5):
            parts.append(f"{i}: {obj[i]}")

        parts.append("...")

        for i in range(n-5, n):
            parts.append(f"{i}: {obj[i]}")

    joined = "\n".join(parts)

    return f"AtomCoords( size={n}\n{joined}\n)"


AtomCoords.__str__ = __atomcoords__str__
AtomCoords.__repr__ = __atomcoords__str__


def _add_evals(obj):
    obj.mass = lambda x : x.evaluate().mass()
    obj.charge = lambda x : x.evaluate().charge()
    obj.coordinates = lambda x : x.evaluate().center_of_mass()
    obj.coords = obj.coordinates
    obj.x = lambda x : x.coordinates().x()
    obj.y = lambda x : x.coordinates().y()
    obj.z = lambda x : x.coordinates().z()


def _get_property(x, key):
    try:
        return x.__orig__property(key)
    except Exception as e:
        saved_exception = e

    mol = x.molecule()

    prop = mol.property(key)

    import sire

    if issubclass(prop.__class__, sire.legacy.Mol.AtomProp):
        vals = []
        for atom in x.atoms():
            vals.append(atom.property(key))

        return vals
    else:
        raise saved_exception


def _apply(objs, func, *args, **kwargs):
    """
    Call the passed function on all views in the container,
    appending the result to a list of results, which
    is returned.

    The function can be either;

    1. a string containing the name of the function to call, or
    2. an actual function (either a normal function or a lambda expression)

    You can optionally pass in positional and keyword arguments
    here that will be passed to the function.

    Args:
        objs (self): The container itself (this is self)
        func (str or function): The function to be called, or the name
                                of the function to be called.

    Returns:
        list: A list of the results, with one result per view in the container.
    """
    result = []

    if str(func) == func:
        # we calling a named function
        for obj in objs:
            result.append(getattr(obj, func)(*args, **kwargs))
    else:
        # we have been passed the function to call
        for obj in objs:
            result.append(func(obj, *args, **kwargs))

    return result


def _apply_reduce(objs, func, reduction_func=None, *args, **kwargs):
    """
    Call the passed function on all views in the container,
    reducing the result into a single value via the 'reduce' function.

    This is equivalent to calling

    ```
    reduce(reduction_func, objs.apply(func, *args, **kwargs))

    The function can be either;

    1. a string containing the name of the function to call, or
    2. an actual function (either a normal function or a lambda expression)

    The reduction function should be a function that can be passed
    to `reduce`. If this isn't passed, then it is assumed to
    be operator.add.

    You can optionally pass in positional and keyword arguments
    here that will be passed to the applied function.

    Args:
        objs (self): The container itself (this is self)
        func (str or function): The function to be called, or the name
                                of the function to be called.
        reduction_func: The function used to reduce the result. This
                        is operator.add by default.

    Returns:
        result: The reduced result
    """
    if reduction_func is None:
        from operator import add
        reduction_func = add

    from functools import reduce
    return reduce(reduction_func, objs.apply(func, *args, **kwargs))


def _add_apply_func(obj):
    if hasattr(obj, "apply"):
        return

    obj.apply = _apply
    obj.apply_reduce = _apply_reduce


def _add_property_func(obj):
    if hasattr(obj, "__orig__property"):
        return

    if hasattr(obj, "property"):
        obj.__orig__property = obj.property

    obj.property = _get_property


for C in [MoleculeView, SelectorMol, SelectorM_Atom_,
          SelectorM_Residue_, SelectorM_Chain_,
          SelectorM_CutGroup_, SelectorM_Segment_]:
    _add_evals(C)
    _add_property_func(C)
    _add_apply_func(C)

for C in [Residue, Chain, Segment]:
    _add_property_func(C)


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


def _cursors(views):
    """Return the Cursors object that contains cursors for all
       of the views in this collection. Note that the `parent`
       cursor of this list will be the molecule
    """
    cursor = views.molecule().cursor()
    return cursor._views(views)


Selector_Atom_.cursor = _cursors
Selector_Residue_.cursor = _cursors
Selector_Chain_.cursor = _cursors
Selector_Segment_.cursor = _cursors

from ._element import *


def _trajectory(obj, map=None):
    from ._trajectory import TrajectoryIterator
    return TrajectoryIterator(obj, map=map)


MoleculeView.trajectory = _trajectory
SelectorM_Atom_.trajectory = _trajectory
SelectorM_Residue_.trajectory = _trajectory
SelectorM_Chain_.trajectory = _trajectory
SelectorM_Segment_.trajectory = _trajectory
SelectorM_CutGroup_.trajectory = _trajectory
SelectorMol.trajectory = _trajectory


def _to_molecules(obj):
    if obj is None:
        return None

    if hasattr(obj, "to_molecule_group"):
        return obj.to_molecule_group().molecules()
    else:
        from ..legacy.Mol import Molecules
        m = Molecules(obj)
        return m


def _energy(obj, obj1=None, map=None):
    from ..mm import calculate_energy

    if map is None:
        if obj1 is None:
            return calculate_energy(obj)
        else:
            return calculate_energy(obj, _to_molecules(obj1))
    elif obj1 is None:
        return calculate_energy(obj, map=map)
    else:
        return calculate_energy(obj, _to_molecules(obj1), map=map)


def _energies(obj, obj1, map=None):
    return obj.apply("energy", obj1, map=map)


def _atom_energy(obj, obj1=None, map=None):
    # An individual atom has a zero energy
    if obj1 is None:
        from ..units import GeneralUnit
        return GeneralUnit(0)
    elif map is None:
        from ..mm import calculate_energy
        return calculate_energy(obj, _to_molecules(obj1))
    else:
        from ..mm import calculate_energy
        return calculate_energy(obj, _to_molecules(obj1), map=map)


def _total_energy(obj, obj1=None, map=None):
    if hasattr(obj, "to_molecule_group"):
        mols = obj.to_molecule_group()
    else:
        from ..legacy.Mol import MoleculeGroup
        mols = MoleculeGroup("all")
        mols.add(obj)

    from ..mm import calculate_energy

    if map is None:
        if obj1 is None:
            return calculate_energy(mols.molecules())
        else:
            return calculate_energy(mols.molecules(),
                                    _to_molecules(obj1))
    elif obj1 is None:
        return calculate_energy(mols.molecules(), map=map)
    else:
        return calculate_energy(mols.molecules(),
                                _to_molecules(obj1), map=map)


Atom.energy = _atom_energy
Residue.energy =  _energy
Chain.energy = _energy
Segment.energy = _energy
CutGroup.energy = _energy
Molecule.energy = _energy

SelectorMol.energy = _total_energy
Selector_Atom_.energy = _total_energy
Selector_Residue_.energy = _total_energy
Selector_Chain_.energy = _total_energy
Selector_Segment_.energy = _total_energy
Selector_CutGroup_.energy = _total_energy

SelectorM_Atom_.energy = _total_energy
SelectorM_Residue_.energy = _total_energy
SelectorM_Chain_.energy = _total_energy
SelectorM_Segment_.energy = _total_energy
SelectorM_CutGroup_.energy = _total_energy

SelectorMol.energies = _energies
Selector_Atom_.energies = _energies
Selector_Residue_.energies = _energies
Selector_Chain_.energies = _energies
Selector_Segment_.energies = _energies
Selector_CutGroup_.energies = _energies

SelectorM_Atom_.energies = _energies
SelectorM_Residue_.energies = _energies
SelectorM_Chain_.energies = _energies
SelectorM_Segment_.energies = _energies
SelectorM_CutGroup_.energies = _energies

from ._view import view as _viewfunc

MoleculeView.view = _viewfunc
SelectorMol.view = _viewfunc
Selector_Atom_.view = _viewfunc
Selector_Residue_.view = _viewfunc
Selector_Chain_.view = _viewfunc
Selector_Segment_.view = _viewfunc
Selector_CutGroup_.view = _viewfunc
SelectorM_Atom_.view = _viewfunc
SelectorM_Residue_.view = _viewfunc
SelectorM_Chain_.view = _viewfunc
SelectorM_Segment_.view = _viewfunc
SelectorM_CutGroup_.view = _viewfunc
