// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Sech.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "hyperbolicfuncs.h"

#include "identities.h"

#include "trigfuncs.h"

#include "hyperbolicfuncs.h"

SireCAS::Sech __copy__(const SireCAS::Sech &other){ return SireCAS::Sech(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Sech_class(){

    { //::SireCAS::Sech
        typedef bp::class_< SireCAS::Sech, bp::bases< SireCAS::SingleFunc, SireCAS::ExBase > > Sech_exposer_t;
        Sech_exposer_t Sech_exposer = Sech_exposer_t( "Sech", "Hyperbolic secant", bp::init< >("Null constructor") );
        bp::scope Sech_scope( Sech_exposer );
        Sech_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") ), "Construct cos(expression)") );
        Sech_exposer.def( bp::init< SireCAS::Sech const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::Sech::evaluate
        
            typedef double ( ::SireCAS::Sech::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Sech::evaluate );
            
            Sech_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "Evaluate this function" );
        
        }
        { //::SireCAS::Sech::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Sech::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Sech::evaluate );
            
            Sech_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "Complex evaluation" );
        
        }
        Sech_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::Sech::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::Sech::typeName );
            
            Sech_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Sech::what
        
            typedef char const * ( ::SireCAS::Sech::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::Sech::what );
            
            Sech_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Sech_exposer.staticmethod( "typeName" );
        Sech_exposer.def( "__copy__", &__copy__);
        Sech_exposer.def( "__deepcopy__", &__copy__);
        Sech_exposer.def( "clone", &__copy__);
        Sech_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::Sech >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Sech_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::Sech >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Sech_exposer.def_pickle(sire_pickle_suite< ::SireCAS::Sech >());
        Sech_exposer.def( "__str__", &__str__< ::SireCAS::Sech > );
        Sech_exposer.def( "__repr__", &__str__< ::SireCAS::Sech > );
        Sech_exposer.def( "__hash__", &::SireCAS::Sech::hash );
    }

}
