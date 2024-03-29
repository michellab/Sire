// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Number.pypp.hpp"

namespace bp = boost::python;

#include "number.h"

#include "number.h"

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireID::Number&){ return "SireID::Number";}

#include "Helpers/release_gil_policy.hpp"

void register_Number_class(){

    { //::SireID::Number
        typedef bp::class_< SireID::Number, boost::noncopyable > Number_exposer_t;
        Number_exposer_t Number_exposer = Number_exposer_t( "Number", "This is the base class of all Number ID objects. A Number\nis used to provide an object with an identifying number.\nThis could be the number of a residue in a molecule, a\nuser-supplied number of a CutGroup in a molecule, or\nperhaps the automatic unique ID numbers of molecules,\nforcefields or molecule groups that are assigned by the\nprogram. The key point of a Number is to provide an ID\nthat can be queried and compared rapidly, and that\ndoes not change as the object is moved between different\ncontainers. Generally an object should keep its number\nthroughout its lifetime.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Number_scope( Number_exposer );
        { //::SireID::Number::hash
        
            typedef ::uint ( ::SireID::Number::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::Number::hash );
            
            Number_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::Number::isNull
        
            typedef bool ( ::SireID::Number::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::Number::isNull );
            
            Number_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::Number::null
        
            typedef ::qint32 ( *null_function_type )(  );
            null_function_type null_function_value( &::SireID::Number::null );
            
            Number_exposer.def( 
                "null"
                , null_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::Number::value
        
            typedef ::qint32 ( ::SireID::Number::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireID::Number::value );
            
            Number_exposer.def( 
                "value"
                , value_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Number_exposer.staticmethod( "null" );
        Number_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::Number >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Number_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::Number >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Number_exposer.def_pickle(sire_pickle_suite< ::SireID::Number >());
        Number_exposer.def( "__str__", &pvt_get_name);
        Number_exposer.def( "__repr__", &pvt_get_name);
        Number_exposer.def( "__hash__", &::SireID::Number::hash );
    }

}
