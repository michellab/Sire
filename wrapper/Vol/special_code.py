###############################################
#
# This file contains special code to help
# with the wrapping of SireVol classes
#
#

from pyplusplus.module_builder import call_policies

def fix_CoordGroupEditor(c):
   c.decls( "translate" ).call_policies = call_policies.return_self()
   c.decls( "rotate" ).call_policies = call_policies.return_self()
   c.decls( "transform" ).call_policies = call_policies.return_self()
   c.decls( "setCoordinates" ).call_policies = call_policies.return_self()
   c.decls( "mapInto" ).call_policies = call_policies.return_self()
   c.decls( "changeFrame" ).call_policies = call_policies.return_self()

special_code = { "SireVol::CoordGroupEditor" : fix_CoordGroupEditor }

def fixMB(mb):
   mb.add_declaration_code("#include <QVector>")
   mb.add_declaration_code("#include \"SireMaths/vector.h\"")
   mb.add_declaration_code("#include \"SireVol/coordgroup.h\"")

implicitly_convertible = [ ("QVector<SireMaths::Vector>","SireVol::CoordGroup") ]

