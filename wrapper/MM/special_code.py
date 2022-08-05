###############################################
#
# This file contains special code to help
# with the wrapping of SireMM classes
#
#

import sys
import pickle

active_headers = pickle.load( open("active_headers.data", "rb") )
mol_headers = pickle.load( open("../Mol/active_headers.data", "rb") )

from pyplusplus.module_builder import call_policies

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
    c.decls("align").call_policies = call_policies.return_self()

    #also include all of the header files included in mover.cpp
    for header in mol_headers["mover.h"].dependencies():
        c.add_declaration_code( "#include %s" % header )

def fix_MolViewProperty(c):
    c.decls("set").call_policies = call_policies.return_self()

def fix_AtomLJs(c):
    fix_MolViewProperty(c)
    c.add_declaration_code("#include \"SireMol/moleculeview.h\"")

def fix_AtomFunctions(c):
    c.add_declaration_code("#include \"SireMol/moleculedata.h\"")

def fix_CLJFunction(c):
    c.add_declaration_code("#include \"SireMM/cljshiftfunction.h\"")

special_code = { "AtomLJs" : fix_AtomLJs,
                 "SireMM::CLJFunction" : fix_CLJFunction,
                 "SireMM::FourAtomFunctions" : fix_AtomFunctions,
                 "SireMM::ThreeAtomFunctions" : fix_AtomFunctions,
                 "SireMM::TwoAtomFunctions" : fix_AtomFunctions,
                 "SireMol::Mover<SireMM::Bond>" : fix_Mover,
                 "SireMol::Mover<SireMM::SelectorBond>" : fix_Mover,
                 "SireMol::Mover<SireMM::Angle>" : fix_Mover,
                 "SireMol::Mover<SireMM::SelectorAngle>" : fix_Mover,
                 "SireMol::Mover<SireMM::Dihedral>" : fix_Mover,
                 "SireMol::Mover<SireMM::SelectorDihedral>" : fix_Mover,
                 "SireMol::Mover<SireMM::Improper>" : fix_Mover,
                 "SireMol::Mover<SireMM::SelectorImproper>" : fix_Mover,
               }
