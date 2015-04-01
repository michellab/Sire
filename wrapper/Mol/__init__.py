#############################
##
## The SireMol module
##
## (C) Christopher Woods
##

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

def __get_typename__(obj):
    try:
        return (obj.typeName().replace("::","_"), obj)
    except:
        if isinstance(obj, float):
            return ("double", obj)
        elif isinstance(obj, int):
            return ("qint64", obj)
        elif isinstance(obj, str):
            return ("QString", obj)
        else:
            return ("QVariant", Sire.Qt.QVariant(obj)) 

def __set_property__(molview, key, property):
    (typename, property) = __get_typename__(property)

    return getattr(molview, "_set_property_%s" % typename)(key, property)     

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

__props = [ AtomCharges, AtomElements ]

for __prop in __props:
    _pvt_property_cludge_fix(__prop)

##########
########## END OF CLUDGY WORKAROUND
##########
