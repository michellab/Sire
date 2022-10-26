###############################################
#
# This file contains special code to help
# with the wrapping of SireSystem classes
#
#

import sys
import pickle

from pyplusplus.module_builder import call_policies

sys.path.append("../AutoGenerate")
from scanheaders import *

def fix_MonitorID(c):
    c.add_declaration_code( "#include \"systemmonitors.h\"" )
    c.add_declaration_code( "#include \"monitorname.h\"" )

def fix_System(c):
    c.add_declaration_code( "#include \"SireBase/slice.h\"" )

fix_MonitorName = fix_MonitorID
fix_MonitorIdx = fix_MonitorID

special_code = { "SireSystem::MonitorID" : fix_MonitorID,
                 "SireSystem::MonitorName" : fix_MonitorName,
                 "SireSystem::MonitorIdx" : fix_MonitorIdx,
                 "SireSystem::System" : fix_System
               }

implicitly_convertible = [ ("SireMol::MoleculeGroup", "SireSystem::AssignerGroup"),
                           ("SireSystem::IDAssigner", "SireSystem::AssignerGroup") ]


def fixMB(mb):   
    mb.add_declaration_code("#include \"SireSystem/freeenergymonitor.h\"")
    mb.add_declaration_code("#include \"SireMol/moleculegroup.h\"")

