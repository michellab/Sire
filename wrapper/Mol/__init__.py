"""
.. currentmodule:: sire.legacy.Mol

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

from .. import Maths as _Maths
from .. import Base as _Base
from .. import ID as _ID
from .. import Qt as _Qt
from .. import CAS as _CAS
from .. import Vol as _Vol
from .. import Units as _Units

from ._Mol import *

from ..Search._Search import install_search_parser as _install_search_parser
_install_search_parser()


def __get_property__(molview, key):
    if hasattr(molview, "property_type"):
        property_type = molview.property_type(key).replace("::","_")
    else:
        property_type = molview.propertyType(key).replace("::","_")

    return getattr(molview, "_get_property_%s" % property_type)(key)


def __get_metadata__(molview, *args):

    if len(args) == 1:
        metakey = args[0]
        if hasattr(molview, "metadata_type"):
            property_type = molview.metadata_type(metakey).replace("::","_")
        else:
            property_type = molview.metadataType(metakey).replace("::","_")
        return getattr(molview, "_get_metadata_%s" % property_type)(metakey)

    elif len(args) == 2:
        (key, metakey) = args
        if hasattr(molview, "metadata_type"):
            property_type = molview.metadata_type(key, metakey).replace("::","_")
        else:
            property_type = molview.metadataType(key, metakey).replace("::","_")
        return getattr(molview, "_get_metadata_%s" % property_type)(key, metakey)

    else:
        raise AttributeError( "Only molview.metadata(metakey) or molview.metadata(key, metakey) are valid!" )

_typename_mapping = {"SireMol_Velocity3D" : "SireMaths_Vector3D_SireUnits_Dimension_Velocity_"}


def __get_typename__(obj):
    try:
        if hasattr(obj, "what"):
            typename = obj.what().replace("::","_")
        elif hasattr(obj, "typename"):
            typename = obj.typename().replace("::","_")
        else:
            typename = obj.typeName().replace("::","_")

        return (_typename_mapping.get(typename, typename), obj)
    except Exception as e:
        if isinstance(obj, float):
            return ("double", obj)
        elif isinstance(obj, int):
            return ("qint64", obj)
        elif isinstance(obj, str):
            return ("QString", obj)
        else:
            raise TypeError(f"Unable to wrap type {type(obj)}: {obj} : {e}")


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
    elif typename.endswith("Coords"):
        return _Maths.Vector(property)
    else:
        return property


def _set_property(molview, key, property):
    if hasattr(molview, "has_property"):
        if molview.has_property(key):
            typename = molview.property_type(key)
            property = _match_to_type(typename, property)
    else:
        if molview.hasProperty(key):
            # get the type of the existing property
            typename = molview.propertyType(key)
            property = _match_to_type(typename, property)

    (typename, property) = __get_typename__(property)

    try:
        return getattr(molview, "_set_property_%s" % typename)(key, property)
    except AttributeError as e:
        # no matching function - see if we have the generic 'PropertyProperty'
        if hasattr(molview, "_set_property_SireBase_PropertyPtr"):
            return molview._set_property_SireBase_PropertyPtr(key, _Base.wrap(property))

        raise e


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
MolEditor.setProperty = _Base.__set_property__


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

__p = _Base.Properties()

def _pvt_property_cludge_fix(C):
    __p.setProperty("c", C())
    try:
        t = __p.property("c").array()
    except Exception as e:
        # this catches cases when the underlying PackedArray2D class
        # is not wrapped - but this class isn't really needed any more
        # print(f"WARNING: Problem with {C} : {e}")
        pass


__props = [ AtomCharges, AtomElements,
            AtomStringArrayProperty,
            AtomPropertyList,
            AtomDoubleArrayProperty,
            AtomIntegerArrayProperty,
            AtomEnergies, AtomFloatProperty,
            AtomForces, AtomIntProperty,
            AtomMasses, AtomPolarisabilities,
            AtomRadicals, AtomRadii,
            AtomStringProperty,
            AtomVariantProperty,
            AtomVelocities
          ]

for __prop in __props:
    _pvt_property_cludge_fix(__prop)

##########
########## END OF CLUDGY WORKAROUND
##########
