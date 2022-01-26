###############################################
#
# This file contains special code to help
# with the wrapping of SireFF classes
#
#

from pyplusplus.module_builder import call_policies

def fix_FFID(c):
   c.add_declaration_code("#include \"forcefields.h\"")
   c.add_declaration_code("#include \"ffidx.h\"")

def fix_Table(c):
   for f in c.decls( "getTable" ):
       if str(f.return_type).find("const") == -1:
           f.exclude()

def fix_ForceFields(c):
   c.add_declaration_code("#include \"forcetable.h\"")
   c.add_declaration_code("#include \"energytable.h\"")

def fix_G2FF(c):
   c.include_files.insert(0, "SireMol/molecule.h" )
   c.include_files.insert(0, "SireMol/partialmolecule.h" )

special_code = { "SireFF::ForceTable" : fix_Table,
                 "SireFF::FieldTable" : fix_Table,
                 "SireFF::PotentialTable" : fix_Table,
                 "SireFF::EnergyTable" : fix_Table,
                 "SireFF::FFID" : fix_FFID,
                 "SireFF::FFIdx" : fix_FFID,
                 "SireFF::FFName" : fix_FFID,
                 "SireFF::ForceFields" : fix_ForceFields,
                 "SireFF::G2FF" : fix_G2FF }

implicitly_convertible = [ ("SireMaths::Vector", "SireFF::PointRef"),
                           ("SireMol::Atom", "SireFF::PointRef"),
                           ("SireFF::Point", "SireFF::PointRef"),
                           ("SireFF::PointRef", "SireFF::PointPtr") ]

def fixMB(mb):
    mb.add_declaration_code("#include \"SireFF/point.h\"")
    mb.add_declaration_code("#include \"SireMol/molecules.h\"")
    mb.add_declaration_code("#include \"SireMol/atom.h\"")
    mb.add_declaration_code("#include \"SireMaths/vector.h\"")

