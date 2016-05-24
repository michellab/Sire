
import re

from pyplusplus.module_builder import call_policies

def findGlobals():
    #read in the information about this module
    lines = open("module_info", "r").readlines()
    root = lines[2].split()[1]
    sourcedir = lines[1].split()[1]

    lines = open("%s/%s/constants.h" % (root,sourcedir), "r").readlines()
    
    FILE = open("_Maths_global_variables.pyman.hpp", "w")
    
    print("#ifndef _Maths_global_variables_hpp", file=FILE)
    print("#define _Maths_global_variables_hpp", file=FILE)
    print("\nvoid register_man_global_variables();\n", file=FILE)
    print("#endif", file=FILE)
     
    FILE.close()

    FILE = open("_Maths_global_variables.pyman.cpp", "w")
    
    print("\n#include \"_Maths_global_variables.pyman.hpp\"", file=FILE)
    print("#include <boost/python.hpp>", file=FILE)
    print("#include \"SireMaths/constants.h\"", file=FILE)
    print("\nusing namespace boost::python;", file=FILE)
    print("using namespace SireMaths;\n", file=FILE)
    
    print("void register_man_global_variables()", file=FILE)
    print("{", file=FILE)
    
    
    for line in lines:
        match = re.search(r"const (double|int)\s+(\w+)\s*=", line)
    
        if match:
            name = match.group(2)
            print("    scope().attr(\"%s\") = %s;\n" % (name,name), file=FILE)

    print("}\n", file=FILE)

def fix_Array2D(c):
   for o in c.operators("()"):
       if o.call_policies is None:
           o.exclude()

   c.add_declaration_code( "#include \"SireBase/array2d.hpp\"" )

def fix_Multi(c):

   try:
       c.decls( "multiplyAdd" ).call_policies = call_policies.return_self()
   except:
       pass

   c.add_declaration_code("#include \"multifloat.h\"")
   c.add_declaration_code("#include \"multiint.h\"")
   c.add_declaration_code("#include \"multidouble.h\"")
   c.add_declaration_code("#include \"multivector.h\"")
   c.add_declaration_code("#include \"multiquaternion.h\"")

special_code = { "SireBase::Array2D<SireBase::PropPtr<SireMaths::Accumulator> >" : fix_Array2D,
                 "SireMaths::MultiFloat" : fix_Multi,
                 "SireMaths::MultiFixed" : fix_Multi,
                 "SireMaths::MultiDouble" : fix_Multi,
                 "SireMaths::MultiUInt" : fix_Multi,
                 "SireMaths::MultiVector" : fix_Multi,
                 "SireMaths::MultiQuaternion" : fix_Multi }

def fixMB(mb):
   mb.add_declaration_code("#include \"_Maths_global_variables.pyman.hpp\"")
   mb.add_registration_code("register_man_global_variables();")

   #add all of the global physical constants to the module
   findGlobals()
