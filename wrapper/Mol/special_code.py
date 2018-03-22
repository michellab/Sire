###############################################
#
# This file contains special code to help
# with the wrapping of SireVol classes
#
#

import sys
import pickle

from pyplusplus.module_builder import call_policies

sys.path.append("../AutoGenerate")
from scanheaders import *

atomprops = pickle.load( open("atomprops.data", "rb") )
chainprops = pickle.load( open("chainprops.data", "rb") )
cgprops = pickle.load( open("cgprops.data", "rb") )
resprops = pickle.load( open("resprops.data", "rb") )
segprops = pickle.load( open("segprops.data", "rb") )
beadprops = pickle.load( open("beadprops.data", "rb") )

active_headers = pickle.load( open("active_headers.data", "rb") )

return_const = "bp::return_value_policy<bp::copy_const_reference>()"
return_self = "bp::return_self< >()"

def fix_MolView(c, molview, props):
   #now add in all of the header files
   for header in props.dependencies():
       c.add_declaration_code( "#include %s" % header )

   #add accessor functions for all of the view properties
   for property in props.properties():
       p = property[0]
       prop = property[1].replace("::","_").replace("<","_").replace(">","_")

       c.add_registration_code( "def( \"_get_property_%s\", &%s::property< %s >, %s)" \
                                      % (prop, molview, p, return_const) )
       c.add_registration_code( "def( \"_get_metadata_%s\", get_Metadata_%s_function1, %s)" \
                                      % (prop, prop, return_const) )
       c.add_registration_code( "def( \"_get_metadata_%s\", &get_Metadata_%s_function2, %s)" \
                                      % (prop, prop, return_const) )

       c.add_declaration_code( """ const %s& get_Metadata_%s_function1(const %s &atom,
                                   const QString &metakey){ return atom.metadata< %s >(metakey); }""" \
                                      % (p, prop, molview, p) )
 
       c.add_declaration_code( """ const %s& get_Metadata_%s_function2(const %s &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< %s >(key, metakey); }""" \
                                      % (p, prop, molview, p) )

def fix_Atom(c):
   fix_MolView(c, "SireMol::Atom", atomprops)

def fix_Chain(c):
   fix_MolView(c, "SireMol::Chain", chainprops)

def fix_CutGroup(c):
   fix_MolView(c, "SireMol::CutGroup", cgprops)

def fix_Residue(c):
   fix_MolView(c, "SireMol::Residue", resprops)

def fix_Segment(c):
   fix_MolView(c, "SireMol::Segment", segprops)

def fix_Bead(c):
   fix_MolView(c, "SireMol::Bead", beadprops)

def fix_MolEditorBase(c):
   c.decls( "removeProperty" ).call_policies = call_policies.return_self()
   c.decls( "removeMetadata" ).call_policies = call_policies.return_self()   

   c.add_registration_code( """def( \"setProperty\",
                               &SireMol::MolEditorBase::setProperty<SireBase::Property>, %s )""" \
                               % (return_self ) )
   
   c.add_registration_code( "def( \"setMetadata\", &set_Metadata_function1, %s)" \
                               % (return_self) )

   c.add_registration_code( "def( \"setMetadata\", &set_Metadata_function2, %s)" \
                               % (return_self) )

   c.add_declaration_code( """SireMol::MolEditorBase& set_Metadata_function1(
                              SireMol::MolEditorBase &molview,
                              const QString &metakey, const SireBase::Property &p)
                              { return molview.setMetadata<SireBase::Property>(metakey, p); }""" )

   c.add_declaration_code( """SireMol::MolEditorBase& set_Metadata_function2(
                              SireMol::MolEditorBase &molview,
                              const QString &key, const QString &metakey, 
                              const SireBase::Property &p)
                              { return molview.setMetadata<SireBase::Property>(key, metakey, p); }""" )

def fix_MolViewEditorBase(c, molview, props):
   c.decls( "removeProperty" ).call_policies = call_policies.return_self()
   c.decls( "removeMetadata" ).call_policies = call_policies.return_self()

   #now add the code to set properties and metadata
   for header in props.dependencies():
       c.add_declaration_code( "#include %s" % header )

   #add accessor functions for all of the atom properties
   for property in props.properties():
       p = property[0]
       p_rep = p.replace("::","_").replace("<","_").replace(">","_")
       prop = property[1].replace("::","_").replace("<","_").replace(">","_")

       c.add_registration_code( """def( \"_set_property_%s\", 
                                   &%s::setProperty< %s >, %s )""" \
                                   % (p_rep, molview, p, return_self ) )

       c.add_registration_code( "def( \"_set_metadata_%s\", &set_Metadata_%s_function1, %s)" \
                                   % (p_rep, prop, return_self) )

       c.add_registration_code( "def( \"_set_metadata_%s\", &set_Metadata_%s_function2, %s)" \
                                   % (p_rep, prop, return_self) )

       c.add_declaration_code( """%s& set_Metadata_%s_function1(
                                  %s &molview,
                                   const QString &metakey, const %s &p)
                                   { return molview.setMetadata< %s >(metakey, p); }""" \
                                      % (molview, prop, molview, p, p) )

       c.add_declaration_code( """%s& set_Metadata_%s_function2(
                                  %s &molview,
                                   const QString &key, const QString &metakey, const %s &p)
                                   { return molview.setMetadata< %s >(key, metakey, p); }""" \
                                      % (molview, prop, molview, p, p) )

def fix_AtomEditorBase(c):
    fix_MolViewEditorBase(c, "SireMol::AtomEditorBase", atomprops)

def fix_ChainEditorBase(c):
    fix_MolViewEditorBase(c, "SireMol::ChainEditorBase", chainprops)

def fix_CGEditorBase(c):
    fix_MolViewEditorBase(c, "SireMol::CGEditorBase", cgprops)

def fix_ResEditorBase(c):
    fix_MolViewEditorBase(c, "SireMol::ResEditorBase", resprops)

def fix_SegEditorBase(c):
    fix_MolViewEditorBase(c, "SireMol::SegEditorBase", segprops)

def fix_BeadEditorBase(c):
    fix_MolViewEditorBase(c, "SireMol::BeadEditorBase", beadprops)

def fix_AtomEditor(c):
   c.decls( "rename" ).call_policies = call_policies.return_self()
   c.decls( "renumber" ).call_policies = call_policies.return_self()

def fix_AtomStructureEditor(c):
   fix_AtomEditor(c)

   c.decls( "reindex" ).call_policies = call_policies.return_self()
   c.decls( "reparent" ).call_policies = call_policies.return_self()   

def fix_AtomSelection(c):
   c.decls( "selectAll" ).call_policies = call_policies.return_self()
   c.decls( "deselectAll" ).call_policies = call_policies.return_self()
   c.decls( "selectNone" ).call_policies = call_policies.return_self()
   c.decls( "select" ).call_policies = call_policies.return_self()
   c.decls( "deselect" ).call_policies = call_policies.return_self()
   c.decls( "selectOnly" ).call_policies = call_policies.return_self()
   c.decls( "invert" ).call_policies = call_policies.return_self()
   c.decls( "intersect" ).call_policies = call_policies.return_self()
   c.decls( "unite" ).call_policies = call_policies.return_self()
   c.decls( "subtract" ).call_policies = call_policies.return_self()
   c.decls( "mask" ).call_policies = call_policies.return_self()

def fix_CGEditor(c):
   c.decls( "rename" ).call_policies = call_policies.return_self()

def fix_CGStructureEditor(c):
   fix_CGEditor(c)
   
   c.decls( "reindex" ).call_policies = call_policies.return_self()
   c.decls( "remove" ).call_policies = call_policies.return_self()
   c.decls( "transfer" ).call_policies = call_policies.return_self()
   c.decls( "transferAll" ).call_policies = call_policies.return_self()

fix_ChainEditor = fix_CGEditor
fix_ChainStructureEditor = fix_CGStructureEditor

fix_SegEditor = fix_CGEditor
fix_SegStructureEditor = fix_CGStructureEditor

def fix_ResEditor(c):
    c.decls( "renumber" ).call_policies = call_policies.return_self()
    c.decls( "rename" ).call_policies = call_policies.return_self()

def fix_ResStructureEditor(c):
    fix_ResEditor(c)

    c.decls( "reindex" ).call_policies = call_policies.return_self()
    c.decls( "reparent" ).call_policies = call_policies.return_self()
    c.decls( "remove" ).call_policies = call_policies.return_self()
    c.decls( "transfer" ).call_policies = call_policies.return_self()
    c.decls( "transferAll" ).call_policies = call_policies.return_self()

def fix_MolEditor(c):
    c.decls( "renumber" ).call_policies = call_policies.return_self()
    c.decls( "rename" ).call_policies = call_policies.return_self()

def fix_MolStructureEditor(c):
    fix_MolEditor(c)

    c.decls( "remove" ).call_policies = call_policies.return_self()
    c.decls( "removeAllAtoms" ).call_policies = call_policies.return_self()
    c.decls( "removeAllCutGroups" ).call_policies = call_policies.return_self()
    c.decls( "removeAllResidues" ).call_policies = call_policies.return_self()
    c.decls( "removeAllChains" ).call_policies = call_policies.return_self()
    c.decls( "removeAllSegments" ).call_policies = call_policies.return_self()

def fix_ConnectivityEditor(c):
    c.decls( "connect" ).call_policies = call_policies.return_self()
    c.decls( "disconnect" ).call_policies = call_policies.return_self()
    c.decls( "disconnectAll" ).call_policies = call_policies.return_self()

def fix_MGNum(c):
    c.add_declaration_code( "#include \"mgid.h\"" )
    c.add_declaration_code( "#include \"mgidx.h\"" )
    c.add_declaration_code( "#include \"mgname.h\"" )
    c.add_declaration_code( "#include \"mgnum.h\"" )
    c.add_declaration_code( "#include \"moleculegroups.h\"")

fix_MGIdx = fix_MGNum
fix_MGName = fix_MGNum

def fix_MolNum(c):
    c.add_declaration_code( "#include \"molid.h\"" )
    c.add_declaration_code( "#include \"molidx.h\"" )
    c.add_declaration_code( "#include \"molnum.h\"" )
    c.add_declaration_code( "#include \"molname.h\"" )
    c.add_declaration_code( "#include \"moleculegroup.h\"" )
    c.add_declaration_code( "#include \"moleculegroups.h\"" )
    c.add_declaration_code( "#include \"mover.hpp\"" )

fix_MolName = fix_MolNum
fix_MolIdx = fix_MolNum

def fix_MolInfo(c):
    c.add_declaration_code( "#include \"moleculeinfodata.h\"" )
    c.add_declaration_code( "#include \"atomselection.h\"" )

def fix_Mover(c):
    c.decls("mapInto").call_policies = call_policies.return_self()
    c.decls("transform").call_policies = call_policies.return_self()
    c.decls("translate").call_policies = call_policies.return_self()
    c.decls("rotate").call_policies = call_policies.return_self()
    c.decls("transform").call_policies = call_policies.return_self()
    c.decls("changeFrame").call_policies = call_policies.return_self()
    c.decls("change").call_policies = call_policies.return_self()
    c.decls("set").call_policies = call_policies.return_self()
    c.decls("setAll").call_policies = call_policies.return_self()
    c.decls("alignTo").call_policies = call_policies.return_self()
    c.decls("align").call_policies = call_policies.return_self()

    #also include all of the header files included in mover.cpp
    for header in active_headers["mover.h"].dependencies():
        c.add_declaration_code( "#include %s" % header )    

def fix_MolViewProperty(c):
    c.add_declaration_code( "#include \"SireMaths/vector.h\"" )
    c.add_declaration_code( "#include \"SireMol/moleculeview.h\"" )
    c.decls("set").call_policies = call_policies.return_self()

def fix_AtomCoords(c):
    fix_MolViewProperty(c)

    c.add_declaration_code("#include \"SireMaths/quaternion.h\"")
    c.add_declaration_code("#include \"SireMaths/matrix.h\"")
    c.add_declaration_code("#include \"SireVol/aabox.h\"")
    c.add_declaration_code("#include \"SireMaths/axisset.h\"")

def fix_CGAtomIdx(c):
    c.add_declaration_code("#include \"cgidx.h\"")
    c.add_declaration_code("#include \"SireID/index.h\"")

def fix_CGIdx(c):
    c.add_declaration_code("#include \"SireID/index.h\"")
    c.add_declaration_code("#include \"cgatomidx.h\"")
    c.add_registration_code("def( other<SireID::Index>() + self )")

def fix_AtomID(c):
    #also include all of the header files included in atomid.cpp
    for header in active_headers["atomid.h"].dependencies():
        c.add_declaration_code( "#include %s" % header )
    
def fix_CGID(c):
    #also include all of the header files included in cgid.cpp
    for header in active_headers["cgid.h"].dependencies():
        c.add_declaration_code( "#include %s" % header )

def fix_ChainID(c):
    #also include all of the header files included in chainid.cpp
    for header in active_headers["chainid.h"].dependencies():
        c.add_declaration_code( "#include %s" % header )

def fix_ResID(c):
    #also include all of the header files included in resid.cpp
    for header in active_headers["resid.h"].dependencies():
        c.add_declaration_code( "#include %s" % header )

def fix_SegID(c):
    #also include all of the header files included in segid.cpp
    for header in active_headers["segid.h"].dependencies():
        c.add_declaration_code( "#include %s" % header )

def fix_BeadID(c):
    #also include all of the header files included in segid.cpp
    for header in active_headers["beadid.h"].dependencies():
        c.add_declaration_code( "#include %s" % header )


def fix_PerturbationSymbols(c):
    c.mem_funs("lambda").rename("Lambda")

special_code = { "SireMol::Atom" : fix_Atom,
                 "SireMol::Editor<SireMol::AtomEditor, SireMol::Atom>" : fix_AtomEditorBase,
                 "SireMol::AtomEditor" : fix_AtomEditor,
                 "SireMol::AtomSelection" : fix_AtomSelection,
                 "SireMol::AtomStructureEditor" : fix_AtomStructureEditor,
                 "SireMol::Mover<SireMol::Atom>" : fix_Mover,
                 "SireMol::Mover<SireMol::Selector<SireMol::Atom> >" : fix_Mover,

                 "SireMol::AtomIdx" : fix_AtomID,
                 "SireMol::AtomName" : fix_AtomID,
                 "SireMol::AtomNum" : fix_AtomID,
                 "SireMol::CGIdx" : fix_CGID,
                 "SireMol::CGName" : fix_CGID,
                 "SireMol::ChainIdx" : fix_ChainID,
                 "SireMol::ChainName" : fix_ChainID,
                 "SireMol::ResIdx" : fix_ResID,
                 "SireMol::ResName" : fix_ResID,
                 "SireMol::ResNum" : fix_ResID,
                 "SireMol::SegIdx" : fix_SegID,
                 "SireMol::SegName" : fix_SegID,
                 "SireMol::BeadIdx" : fix_BeadID,
                 "SireMol::BeadNum" : fix_BeadID,

                 "SireMol::Bead" : fix_Bead,
                 "SireMol::Editor<SireMol::BeadEditor, SireMol::Bead>" : fix_BeadEditorBase,
                 "SireMol::Mover<SireMol::Bead>" : fix_Mover,
                 
                 "SireMol::Mover<SireMol::Beads>" : fix_Mover,

                 "SireMol::CutGroup" : fix_CutGroup,
                 "SireMol::Editor<SireMol::CGEditor, SireMol::CutGroup>" : fix_CGEditorBase,
                 "SireMol::CGEditor" : fix_CGEditor,
                 "SireMol::CGStructureEditor" : fix_CGStructureEditor,
                 "SireMol::Mover<SireMol::CutGroup>" : fix_Mover,
                 "SireMol::Mover<SireMol::Selector<SireMol::CutGroup> >" : fix_Mover,

                 "SireMol::Chain" : fix_Chain,
                 "SireMol::Editor<SireMol::ChainEditor, SireMol::Chain>" : fix_ChainEditorBase,
                 "SireMol::ChainEditor" : fix_ChainEditor,
                 "SireMol::ChainStructureEditor" : fix_ChainStructureEditor,
                 "SireMol::Mover<SireMol::Chain>" : fix_Mover,
                 "SireMol::Mover<SireMol::Selector<SireMol::Chain> >" : fix_Mover,

                 "SireMol::Residue" : fix_Residue,
                 "SireMol::Editor<SireMol::ResEditor, SireMol::Residue>" : fix_ResEditorBase,
                 "SireMol::ResEditor" : fix_ResEditor,
                 "SireMol::ResStructureEditor" : fix_ResStructureEditor,
                 "SireMol::Mover<SireMol::Residue>" : fix_Mover,
                 "SireMol::Mover<SireMol::Selector<SireMol::Residue> >" : fix_Mover,

                 "SireMol::Segment" : fix_Segment,
                 "SireMol::Editor<SireMol::SegEditor, SireMol::Segment>" : fix_SegEditorBase,
                 "SireMol::SegEditor" : fix_SegEditor,
                 "SireMol::SegStructureEditor" : fix_SegStructureEditor,
                 "SireMol::Mover<SireMol::Segment>" : fix_Mover,
                 "SireMol::Mover<SireMol::Selector<SireMol::Segment> >" : fix_Mover,

                 "SireMol::MolEditor" : fix_MolEditor,
                 "SireMol::Editor<SireMol::MolEditor, SireMol::Molecule>" : fix_MolEditorBase,
                 "SireMol::MolStructureEditor" : fix_MolStructureEditor,
                 "SireMol::Mover<SireMol::Molecule>" : fix_Mover,
                 "SireMol::Mover<SireMol::PartialMolecule>" : fix_Mover,

                 "AtomStringProperty" : fix_MolViewProperty,
                 "AtomIntProperty" : fix_MolViewProperty,
                 "AtomFloatProperty" : fix_MolViewProperty,
                 "AtomVariantProperty" : fix_MolViewProperty,
                 "BeadStringProperty" : fix_MolViewProperty,
                 "BeadIntProperty" : fix_MolViewProperty,
                 "BeadFloatProperty" : fix_MolViewProperty,
                 "BeadVariantProperty" : fix_MolViewProperty,
                 "CGStringProperty" : fix_MolViewProperty,
                 "CGIntProperty" : fix_MolViewProperty,
                 "CGFloatProperty" : fix_MolViewProperty,
                 "CGVariantProperty" : fix_MolViewProperty,
                 "ResStringProperty" : fix_MolViewProperty,
                 "ResIntProperty" : fix_MolViewProperty,
                 "ResFloatProperty" : fix_MolViewProperty,
                 "ResVariantProperty" : fix_MolViewProperty,
                 "ChainStringProperty" : fix_MolViewProperty,
                 "ChainIntProperty" : fix_MolViewProperty,
                 "ChainFloatProperty" : fix_MolViewProperty,
                 "ChainVariantProperty" : fix_MolViewProperty,
                 "SegStringProperty" : fix_MolViewProperty,
                 "SegIntProperty" : fix_MolViewProperty,
                 "SegFloatProperty" : fix_MolViewProperty,
                 "SegVariantProperty" : fix_MolViewProperty,

                 "AtomBeads" : fix_MolViewProperty,
                 "AtomCoords" : fix_AtomCoords,
                 "AtomCharges" : fix_MolViewProperty,
                 "AtomElements" : fix_MolViewProperty,
                 "AtomEnergies" : fix_MolViewProperty,
                 "AtomForces" : fix_MolViewProperty,                
                 "AtomMasses"  : fix_MolViewProperty,
                 "AtomVelocities" : fix_MolViewProperty,
                 "AtomPolarisabilities" : fix_MolViewProperty,
                 "AtomRadii" : fix_MolViewProperty,

                 "SireMol::ConnectivityEditor" : fix_ConnectivityEditor,
                 "SireMol::MGName" : fix_MGName,
                 "SireMol::MGIdx" : fix_MGIdx,
                 "SireMol::MGNum" : fix_MGNum,
                 "SireMol::MolNum" : fix_MolNum,
                 "SireMol::MolName" : fix_MolName,
                 "SireMol::MolIdx" : fix_MolIdx,
                 "SireMol::MolInfo" : fix_MolInfo, 
                 "SireMol::MoleculeInfo" : fix_MolInfo,

                 "SireMol::PerturbationSymbols" : fix_PerturbationSymbols,

                 "SireMol::CGIdx" : fix_CGIdx,
                 "SireMol::CGAtomIdx" : fix_CGAtomIdx }

implicitly_convertible = [ ("SireMol::AtomID", "SireMol::AtomIdentifier"),
                           ("SireMol::CGID", "SireMol::CGIdentifier"),
                           ("SireMol::ChainID", "SireMol::ChainIdentifier"),
                           ("SireMol::ResID", "SireMol::ResIdentifier"),
                           ("SireMol::SegID", "SireMol::SegIdentifier"),
                           ("SireMol::MolID", "SireMol::MolIdentifier"),
                           ("SireMol::MGID", "SireMol::MGIdentifier"),
                           ("SireMol::MoleculeView", "SireMol::MoleculeData"),
                           ("SireMol::MoleculeView", "SireMol::PartialMolecule"),
                           ("SireMol::MoleculeInfoData", "SireMol::MoleculeInfo"),
                           ("SireMol::MoleculeInfo", "SireMol::MoleculeInfoData") ]

def fixMB(mb):
    mb.add_declaration_code("#include \"SireMol/moleculedata.h\"")
    mb.add_declaration_code("#include \"SireMol/moleculeview.h\"")
    mb.add_declaration_code("#include \"SireMol/partialmolecule.h\"")
    mb.add_declaration_code("#include \"SireMol/mover.hpp\"")
    mb.add_declaration_code("#include \"SireMol/mgidentifier.h\"")
    mb.add_declaration_code("#include \"SireMol/moleculeinfo.h\"")
