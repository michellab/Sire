
def fix_Array2D(c):
   for o in c.operators("()"):
       if o.call_policies is None:
           o.exclude()

special_code = { "SireBase::Array2D<double>" : fix_Array2D }

implicitly_convertible = [ ("QString", "SireBase::PropertyName"),
                           ("SireBase::Property", "SireBase::PropertyName"),
                           ("QHash<QString,SireBase::PropertyName>", "SireBase::PropertyMap"),
                           ("SireBase::DoubleArrayProperty", "SireBase::StringArrayProperty"),
                           ("SireBase::DoubleArrayProperty", "SireBase::IntegerArrayProperty"),
                           ("SireBase::DoubleArrayProperty", "SireBase::PropertyList"),
                           ("SireBase::IntegerArrayProperty", "SireBase::StringArrayProperty"),
                           ("SireBase::IntegerArrayProperty", "SireBase::DoubleArrayProperty"),
                           ("SireBase::IntegerArrayProperty", "SireBase::PropertyList"),
                           ("SireBase::StringArrayProperty", "SireBase::PropertyList") ]

def fixMB(mb):
    mb.add_declaration_code( "#include \"SireBase/propertymap.h\"" )
    mb.add_declaration_code( "#include \"SireBase/stringproperty.h\"" )
    mb.add_declaration_code( "#include \"SireBase/numberproperty.h\"" )
    mb.add_declaration_code( "#include \"SireBase/lengthproperty.h\"" )
    mb.add_declaration_code( "#include \"SireBase/propertylist.h\"" )
    mb.add_declaration_code( "#include <QString>" )
    mb.add_declaration_code( "#include \"SireBase/slice.h\"" )
    mb.add_declaration_code( "void autoconvert_Slice();")
    mb.add_definition_code( "autoconvert_Slice();")

