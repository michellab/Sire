"""
.. currentmodule:: Sire.Mol

Classes
=======

.. autosummary::
    :toctree: generated/

    AbsFromMass
    AbsFromNumber
    AngleID
    AnglePerturbation
    Atom
    AtomBeads
    AtomCharges
    AtomCoords
    AtomCutting
    AtomDoubleArrayProperty
    AtomEditor
    AtomElements
    AtomEnergies
    AtomFloatProperty
    AtomForces
    AtomID
    AtomIDMatcher
    AtomIdx
    AtomIdxMatcher
    AtomIntegerArrayProperty
    AtomIntProperty
    AtomMasses
    AtomMatcher
    AtomMatchInverter
    AtomMCSMatcher
    AtomMultiMatcher
    AtomName
    AtomNameMatcher
    AtomNum
    AtomPolarisabilities
    AtomPropertyList
    AtomRadicals
    AtomRadii
    AtomResultMatcher
    AtomSelection
    AtomStringArrayProperty
    AtomStringProperty
    AtomStructureEditor
    AtomVariantProperty
    AtomVelocities
    Bead
    BeadEditor
    BeadFloatProperty
    BeadID
    BeadIdx
    Beading
    BeadIntProperty
    BeadNum
    Beads
    BeadStringProperty
    BeadVariantProperty
    BondHunter
    BondID
    BondPerturbation
    BondType
    Chain
    ChainAtomID
    ChainEditor
    ChainFloatProperty
    ChainID
    ChainIdx
    ChainIntProperty
    ChainName
    ChainResID
    ChainStringProperty
    ChainStructureEditor
    ChainsWithAtoms
    ChainsWithRes
    ChainVariantProperty
    ChargePerturbation
    ChemicalBondHunter
    Connectivity
    ConnectivityEditor
    CovalentBondHunter
    CuttingFunction
    DihedralID
    DihedralPerturbation
    Element
    Evaluator
    Force3D
    GeometryPerturbation
    GeometryPerturbations
    ImproperID
    MGID
    MGIdx
    MGName
    MGNum
    Molecule
    MoleculeBeading
    MoleculeGroup
    MoleculeGroups
    MoleculeInfo
    Molecules
    MoleculeView
    MolEditor
    MolID
    MolIdx
    MolInfo
    MolName
    MolNum
    MolResNum
    MolStructureEditor
    MoverBase
    PartialMolecule
    Perturbation
    Perturbations
    PerturbationSymbols
    Radical
    RelFromMass
    RelFromNumber
    ResAtomID
    ResEditor
    ResFloatProperty
    ResID
    Residue
    ResidueBeading
    ResidueCutting
    ResIdx
    ResIdxAtomCoordMatcher
    ResIdxAtomMCSMatcher
    ResIdxAtomNameMatcher
    ResIntProperty
    ResName
    ResNum
    ResNumAtomNameMatcher
    ResStringProperty
    ResStructureEditor
    ResVariantProperty
    ResWithAtoms
    SegAtomID
    SegChainID
    SegEditor
    SegFloatProperty
    SegID
    SegIdx
    SegIntProperty
    Segment
    SegName
    SegResID
    SegStringProperty
    SegStructureEditor
    SegsWithAtoms
    SegVariantProperty
    Select
    SelectResult
    SelectResultMover
    SpecifyMol
    Stereoscopy
    UserBeading
    Velocity3D
    VolumeMap
    WeightFunction
    Within

Functions
=========

.. autosummary::
    :toctree: generated/

    getAlignment

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
        if e.__class__.__name__ == "ArgumentError" or e.__class__.__name__ == "AttributeError":
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

# Here I will define some functions that make accessing
# things from moleculeviews more convenient

def __fixed__getitem__(obj, key):
    if type(key) is slice:
        if obj.what() == "SireMol::Chain":
            return [residue for residue in obj.residues(key)]
        else:
            return [atom for atom in obj.atoms(key)]
    else:
        return obj.__orig__getitem__(key)


def __fixed__atoms__(obj, idx=None):
    if idx is None:
        return obj.__orig__atoms()
    elif type(idx) is range:
        return obj.__orig__atoms(list(idx))
    else:
        return obj.__orig__atoms(idx)


def __fixed__residues__(obj, idx=None):
    if idx is None:
        return obj.__orig__residues()
    elif type(idx) is range:
        return obj.__orig__residues(list(idx))
    else:
        return obj.__orig__residues(idx)


def __fixed__chains__(obj, idx=None):
    if idx is None:
        return obj.__orig__chains()
    elif type(idx) is range:
        return obj.__orig__chains(list(idx))
    else:
        return obj.__orig__chains(idx)


def __fixed__segments__(obj, idx=None):
    if idx is None:
        return obj.__orig__segments()
    elif type(idx) is range:
        return obj.__orig__segments(list(idx))
    else:
        return obj.__orig__segments(idx)


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


for C in [Atom, CutGroup, Residue, Chain, Segment, Molecule]:
    __fix_getitem(C)


MoleculeView.coordinates = lambda x : x.property("coordinates")
Atom.element = lambda x : x.property("element")

MoleculeView.mass = lambda x : x.evaluate().mass()
MoleculeView.charge = lambda x : x.evaluate().charge()

def _get_atom_mass(x):
    if x.hasProperty("mass"):
        return x.property("mass")
    elif x.hasProperty("element"):
        return x.property("element").mass()
    else:
        return 0

Atom.mass = _get_atom_mass

def _get_atom_charge(x):
    if x.hasProperty("charge"):
        return x.property("charge")
    elif x.hasProperty("formal_charge"):
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
