// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Property.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireError/getbacktrace.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "generalunitproperty.h"

#include "property.h"

#include "propertylist.h"

#include <QDebug>

#include <QMutex>

#include "property.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Property_class(){

    { //::SireBase::Property
        typedef bp::class_< SireBase::Property, boost::noncopyable > Property_exposer_t;
        Property_exposer_t Property_exposer = Property_exposer_t( "Property", "This is the base class of all properties that may be attached to a\nmolecule. Properties are used to assign extra information to a molecule,\nwhich may then be carried by the molecule throughout its passage\nthrough the simulation. Examples of properties may include the file\nfrom which the molecule was read, the charge parameters on the atoms,\nthe PDB code etc.\n\nProperties form a polymorphic hierarchy which are implicitly shared\nvia SireBase::SharedPolyPointer smart pointers.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Property_scope( Property_exposer );
        { //::SireBase::Property::asABoolean
        
            typedef bool ( ::SireBase::Property::*asABoolean_function_type)(  ) const;
            asABoolean_function_type asABoolean_function_value( &::SireBase::Property::asABoolean );
            
            Property_exposer.def( 
                "asABoolean"
                , asABoolean_function_value
                , bp::release_gil_policy()
                , "Return this property converted to a bool. This throws an invalid\ncast if this is not possible" );
        
        }
        { //::SireBase::Property::asADouble
        
            typedef double ( ::SireBase::Property::*asADouble_function_type)(  ) const;
            asADouble_function_type asADouble_function_value( &::SireBase::Property::asADouble );
            
            Property_exposer.def( 
                "asADouble"
                , asADouble_function_value
                , bp::release_gil_policy()
                , "Return this property converted to a double. This throws an invalid\ncast if this is not possible" );
        
        }
        { //::SireBase::Property::asAString
        
            typedef ::QString ( ::SireBase::Property::*asAString_function_type)(  ) const;
            asAString_function_type asAString_function_value( &::SireBase::Property::asAString );
            
            Property_exposer.def( 
                "asAString"
                , asAString_function_value
                , bp::release_gil_policy()
                , "Return this property converted to a string. This throws an invalid\ncast if this is not possible" );
        
        }
        { //::SireBase::Property::asAUnit
        
            typedef ::SireUnits::Dimension::GeneralUnit ( ::SireBase::Property::*asAUnit_function_type)(  ) const;
            asAUnit_function_type asAUnit_function_value( &::SireBase::Property::asAUnit );
            
            Property_exposer.def( 
                "asAUnit"
                , asAUnit_function_value
                , bp::release_gil_policy()
                , "Return this property converted to a unit" );
        
        }
        { //::SireBase::Property::asAnArray
        
            typedef ::SireBase::PropertyList ( ::SireBase::Property::*asAnArray_function_type)(  ) const;
            asAnArray_function_type asAnArray_function_value( &::SireBase::Property::asAnArray );
            
            Property_exposer.def( 
                "asAnArray"
                , asAnArray_function_value
                , bp::release_gil_policy()
                , "Return this property converted to an array property. By default, this\nautomatically puts this property into a PropertyList and returns that" );
        
        }
        { //::SireBase::Property::asAnInteger
        
            typedef int ( ::SireBase::Property::*asAnInteger_function_type)(  ) const;
            asAnInteger_function_type asAnInteger_function_value( &::SireBase::Property::asAnInteger );
            
            Property_exposer.def( 
                "asAnInteger"
                , asAnInteger_function_value
                , bp::release_gil_policy()
                , "Return this property converted to an integer. This throws an invalid\ncast if this is not possible" );
        
        }
        { //::SireBase::Property::copy
        
            typedef void ( ::SireBase::Property::*copy_function_type)( ::SireBase::Property const & ) ;
            copy_function_type copy_function_value( &::SireBase::Property::copy );
            
            Property_exposer.def( 
                "copy"
                , copy_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Property::equals
        
            typedef bool ( ::SireBase::Property::*equals_function_type)( ::SireBase::Property const & ) const;
            equals_function_type equals_function_value( &::SireBase::Property::equals );
            
            Property_exposer.def( 
                "equals"
                , equals_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Property::isABoolean
        
            typedef bool ( ::SireBase::Property::*isABoolean_function_type)(  ) const;
            isABoolean_function_type isABoolean_function_value( &::SireBase::Property::isABoolean );
            
            Property_exposer.def( 
                "isABoolean"
                , isABoolean_function_value
                , bp::release_gil_policy()
                , "Return whether or not this property holds a bool (or can convert\nto a bool)" );
        
        }
        { //::SireBase::Property::isADouble
        
            typedef bool ( ::SireBase::Property::*isADouble_function_type)(  ) const;
            isADouble_function_type isADouble_function_value( &::SireBase::Property::isADouble );
            
            Property_exposer.def( 
                "isADouble"
                , isADouble_function_value
                , bp::release_gil_policy()
                , "Return whether or not this property holds a double (or can convert\nto a double)" );
        
        }
        { //::SireBase::Property::isAString
        
            typedef bool ( ::SireBase::Property::*isAString_function_type)(  ) const;
            isAString_function_type isAString_function_value( &::SireBase::Property::isAString );
            
            Property_exposer.def( 
                "isAString"
                , isAString_function_value
                , bp::release_gil_policy()
                , "Return whether or not this property holds a string (or can convert\nto a string)" );
        
        }
        { //::SireBase::Property::isAUnit
        
            typedef bool ( ::SireBase::Property::*isAUnit_function_type)(  ) const;
            isAUnit_function_type isAUnit_function_value( &::SireBase::Property::isAUnit );
            
            Property_exposer.def( 
                "isAUnit"
                , isAUnit_function_value
                , bp::release_gil_policy()
                , "Return whether or not this property holds a unit (or can convert\nto a unit)\n" );
        
        }
        { //::SireBase::Property::isAnArray
        
            typedef bool ( ::SireBase::Property::*isAnArray_function_type)(  ) const;
            isAnArray_function_type isAnArray_function_value( &::SireBase::Property::isAnArray );
            
            Property_exposer.def( 
                "isAnArray"
                , isAnArray_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is an array property (or can convert to an\narray property)" );
        
        }
        { //::SireBase::Property::isAnInteger
        
            typedef bool ( ::SireBase::Property::*isAnInteger_function_type)(  ) const;
            isAnInteger_function_type isAnInteger_function_value( &::SireBase::Property::isAnInteger );
            
            Property_exposer.def( 
                "isAnInteger"
                , isAnInteger_function_value
                , bp::release_gil_policy()
                , "Return whether or not this property holds an integer (or can convert\nto an integer)" );
        
        }
        { //::SireBase::Property::load
        
            typedef void ( ::SireBase::Property::*load_function_type)( ::QDataStream & ) ;
            load_function_type load_function_value( &::SireBase::Property::load );
            
            Property_exposer.def( 
                "load"
                , load_function_value
                , ( bp::arg("ds") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Property::null
        
            typedef ::SireBase::NullProperty const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireBase::Property::null );
            
            Property_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the global null property" );
        
        }
        { //::SireBase::Property::save
        
            typedef void ( ::SireBase::Property::*save_function_type)( ::QDataStream & ) const;
            save_function_type save_function_value( &::SireBase::Property::save );
            
            Property_exposer.def( 
                "save"
                , save_function_value
                , ( bp::arg("ds") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Property::toString
        
            typedef ::QString ( ::SireBase::Property::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::Property::toString );
            
            Property_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Default toString() function for properties - it would\nhelp if all properties output something more sensible" );
        
        }
        { //::SireBase::Property::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::Property::typeName );
            
            Property_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Property::what
        
            typedef char const * ( ::SireBase::Property::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireBase::Property::what );
            
            Property_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Property_exposer.staticmethod( "null" );
        Property_exposer.staticmethod( "typeName" );
        Property_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::Property >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Property_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::Property >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Property_exposer.def_pickle(sire_pickle_suite< ::SireBase::Property >());
        Property_exposer.def( "__str__", &__str__< ::SireBase::Property > );
        Property_exposer.def( "__repr__", &__str__< ::SireBase::Property > );
    }

}
