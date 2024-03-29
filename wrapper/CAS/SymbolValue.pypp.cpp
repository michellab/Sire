// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SymbolValue.pypp.hpp"

namespace bp = boost::python;

#include "symbolvalue.h"

#include "symbolvalue.h"

#include "SireCAS/values.h"

SireCAS::SymbolValue __copy__(const SireCAS::SymbolValue &other){ return SireCAS::SymbolValue(other); }

const char* pvt_get_name(const SireCAS::SymbolValue&){ return "SireCAS::SymbolValue";}

#include "Helpers/release_gil_policy.hpp"

void register_SymbolValue_class(){

    { //::SireCAS::SymbolValue
        typedef bp::class_< SireCAS::SymbolValue > SymbolValue_exposer_t;
        SymbolValue_exposer_t SymbolValue_exposer = SymbolValue_exposer_t( "SymbolValue", "Small class that holds a SymbolID number and an associated value", bp::init< SireCAS::SymbolID, double >(( bp::arg("id"), bp::arg("val") ), "") );
        bp::scope SymbolValue_scope( SymbolValue_exposer );
        { //::SireCAS::SymbolValue::ID
        
            typedef ::SireCAS::SymbolID ( ::SireCAS::SymbolValue::*ID_function_type)(  ) const;
            ID_function_type ID_function_value( &::SireCAS::SymbolValue::ID );
            
            SymbolValue_exposer.def( 
                "ID"
                , ID_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::SymbolValue::value
        
            typedef double ( ::SireCAS::SymbolValue::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireCAS::SymbolValue::value );
            
            SymbolValue_exposer.def( 
                "value"
                , value_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SymbolValue_exposer.def( bp::self + bp::self );
        SymbolValue_exposer.def( bp::self + bp::other< SireCAS::Values >() );
        SymbolValue_exposer.def( "__copy__", &__copy__);
        SymbolValue_exposer.def( "__deepcopy__", &__copy__);
        SymbolValue_exposer.def( "clone", &__copy__);
        SymbolValue_exposer.def( "__str__", &pvt_get_name);
        SymbolValue_exposer.def( "__repr__", &pvt_get_name);
    }

}
