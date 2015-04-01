
def fix_Array2D(c):
   for o in c.operators("()"):
       if o.call_policies is None:
           o.exclude()

special_code = { "SireBase::Array2D<double>" : fix_Array2D }

implicitly_convertible = [ ("QString", "SireBase::PropertyName"),
                           ("QString", "SireBase::VariantProperty"),
                           ("double", "SireBase::VariantProperty"),
                           ("QString", "SireBase::PropertyName"),
                           ("SireBase::Property", "SireBase::PropertyName"),
                           ("QString", "SireBase::StringProperty"),
                           ("SireBase::StringProperty", "QString"),
                           ("double", "SireBase::NumberProperty"),
                           ("int", "SireBase::NumberProperty"),
                           ("SireBase::LengthProperty", "SireUnits::Dimension::Length"),
                           ("QHash<QString,SireBase::PropertyName>", "SireBase::PropertyMap") ] 

def fixMB(mb):
    mb.add_declaration_code( "#include \"SireBase/propertymap.h\"" )
    mb.add_declaration_code( "#include \"SireBase/stringproperty.h\"" )
    mb.add_declaration_code( "#include \"SireBase/numberproperty.h\"" )
    mb.add_declaration_code( "#include \"SireBase/lengthproperty.h\"" )
    mb.add_declaration_code( "#include <QString>" )

