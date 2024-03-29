// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AlwaysFalse.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireError/errors.h"

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "complexvalues.h"

#include "conditional.h"

#include "expressions.h"

#include "functions.h"

#include "identities.h"

#include "symbols.h"

#include "values.h"

#include "conditional.h"

SireCAS::AlwaysFalse __copy__(const SireCAS::AlwaysFalse &other){ return SireCAS::AlwaysFalse(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_AlwaysFalse_class(){

    { //::SireCAS::AlwaysFalse
        typedef bp::class_< SireCAS::AlwaysFalse, bp::bases< SireCAS::Condition, SireCAS::ExBase > > AlwaysFalse_exposer_t;
        AlwaysFalse_exposer_t AlwaysFalse_exposer = AlwaysFalse_exposer_t( "AlwaysFalse", "This is an overloaded conditional that is always false", bp::init< >("Constructor") );
        bp::scope AlwaysFalse_scope( AlwaysFalse_exposer );
        AlwaysFalse_exposer.def( bp::init< SireCAS::AlwaysFalse const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::AlwaysFalse::alwaysFalse
        
            typedef bool ( ::SireCAS::AlwaysFalse::*alwaysFalse_function_type)(  ) const;
            alwaysFalse_function_type alwaysFalse_function_value( &::SireCAS::AlwaysFalse::alwaysFalse );
            
            AlwaysFalse_exposer.def( 
                "alwaysFalse"
                , alwaysFalse_function_value
                , bp::release_gil_policy()
                , "AlwaysFalse is always false" );
        
        }
        { //::SireCAS::AlwaysFalse::alwaysTrue
        
            typedef bool ( ::SireCAS::AlwaysFalse::*alwaysTrue_function_type)(  ) const;
            alwaysTrue_function_type alwaysTrue_function_value( &::SireCAS::AlwaysFalse::alwaysTrue );
            
            AlwaysFalse_exposer.def( 
                "alwaysTrue"
                , alwaysTrue_function_value
                , bp::release_gil_policy()
                , "AlwaysFalse is never true" );
        
        }
        { //::SireCAS::AlwaysFalse::children
        
            typedef ::SireCAS::Expressions ( ::SireCAS::AlwaysFalse::*children_function_type)(  ) const;
            children_function_type children_function_value( &::SireCAS::AlwaysFalse::children );
            
            AlwaysFalse_exposer.def( 
                "children"
                , children_function_value
                , bp::release_gil_policy()
                , "False has no children" );
        
        }
        { //::SireCAS::AlwaysFalse::evaluate
        
            typedef double ( ::SireCAS::AlwaysFalse::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::AlwaysFalse::evaluate );
            
            AlwaysFalse_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "False is never true" );
        
        }
        { //::SireCAS::AlwaysFalse::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::AlwaysFalse::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::AlwaysFalse::evaluate );
            
            AlwaysFalse_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "False is never true" );
        
        }
        { //::SireCAS::AlwaysFalse::evaluateCondition
        
            typedef bool ( ::SireCAS::AlwaysFalse::*evaluateCondition_function_type)( ::SireCAS::Values const & ) const;
            evaluateCondition_function_type evaluateCondition_function_value( &::SireCAS::AlwaysFalse::evaluateCondition );
            
            AlwaysFalse_exposer.def( 
                "evaluateCondition"
                , evaluateCondition_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "False is never true" );
        
        }
        { //::SireCAS::AlwaysFalse::evaluateCondition
        
            typedef bool ( ::SireCAS::AlwaysFalse::*evaluateCondition_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluateCondition_function_type evaluateCondition_function_value( &::SireCAS::AlwaysFalse::evaluateCondition );
            
            AlwaysFalse_exposer.def( 
                "evaluateCondition"
                , evaluateCondition_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "False is never true" );
        
        }
        { //::SireCAS::AlwaysFalse::expand
        
            typedef ::QList< SireCAS::Factor > ( ::SireCAS::AlwaysFalse::*expand_function_type)( ::SireCAS::Symbol const & ) const;
            expand_function_type expand_function_value( &::SireCAS::AlwaysFalse::expand );
            
            AlwaysFalse_exposer.def( 
                "expand"
                , expand_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "False cannot be expanded" );
        
        }
        { //::SireCAS::AlwaysFalse::functions
        
            typedef ::SireCAS::Functions ( ::SireCAS::AlwaysFalse::*functions_function_type)(  ) const;
            functions_function_type functions_function_value( &::SireCAS::AlwaysFalse::functions );
            
            AlwaysFalse_exposer.def( 
                "functions"
                , functions_function_value
                , bp::release_gil_policy()
                , "There are no functions in false" );
        
        }
        { //::SireCAS::AlwaysFalse::hash
        
            typedef ::uint ( ::SireCAS::AlwaysFalse::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireCAS::AlwaysFalse::hash );
            
            AlwaysFalse_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "Hash false" );
        
        }
        { //::SireCAS::AlwaysFalse::isComplex
        
            typedef bool ( ::SireCAS::AlwaysFalse::*isComplex_function_type)(  ) const;
            isComplex_function_type isComplex_function_value( &::SireCAS::AlwaysFalse::isComplex );
            
            AlwaysFalse_exposer.def( 
                "isComplex"
                , isComplex_function_value
                , bp::release_gil_policy()
                , "False is never complex" );
        
        }
        { //::SireCAS::AlwaysFalse::isCompound
        
            typedef bool ( ::SireCAS::AlwaysFalse::*isCompound_function_type)(  ) const;
            isCompound_function_type isCompound_function_value( &::SireCAS::AlwaysFalse::isCompound );
            
            AlwaysFalse_exposer.def( 
                "isCompound"
                , isCompound_function_value
                , bp::release_gil_policy()
                , "False is always simple" );
        
        }
        { //::SireCAS::AlwaysFalse::isConstant
        
            typedef bool ( ::SireCAS::AlwaysFalse::*isConstant_function_type)(  ) const;
            isConstant_function_type isConstant_function_value( &::SireCAS::AlwaysFalse::isConstant );
            
            AlwaysFalse_exposer.def( 
                "isConstant"
                , isConstant_function_value
                , bp::release_gil_policy()
                , "Truth is always constant" );
        
        }
        { //::SireCAS::AlwaysFalse::isFunction
        
            typedef bool ( ::SireCAS::AlwaysFalse::*isFunction_function_type)( ::SireCAS::Symbol const & ) const;
            isFunction_function_type isFunction_function_value( &::SireCAS::AlwaysFalse::isFunction );
            
            AlwaysFalse_exposer.def( 
                "isFunction"
                , isFunction_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "This is not a function of anything" );
        
        }
        { //::SireCAS::AlwaysFalse::isNull
        
            typedef bool ( ::SireCAS::AlwaysFalse::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireCAS::AlwaysFalse::isNull );
            
            AlwaysFalse_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "False is never empty" );
        
        }
        { //::SireCAS::AlwaysFalse::operator=
        
            typedef ::SireCAS::AlwaysFalse & ( ::SireCAS::AlwaysFalse::*assign_function_type)( ::SireCAS::AlwaysFalse const & ) ;
            assign_function_type assign_function_value( &::SireCAS::AlwaysFalse::operator= );
            
            AlwaysFalse_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AlwaysFalse_exposer.def( bp::self == bp::self );
        AlwaysFalse_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::AlwaysFalse::simplify
        
            typedef ::SireCAS::Expression ( ::SireCAS::AlwaysFalse::*simplify_function_type)( int ) const;
            simplify_function_type simplify_function_value( &::SireCAS::AlwaysFalse::simplify );
            
            AlwaysFalse_exposer.def( 
                "simplify"
                , simplify_function_value
                , ( bp::arg("options") )
                , bp::release_gil_policy()
                , "This cannot be further simplified" );
        
        }
        { //::SireCAS::AlwaysFalse::substitute
        
            typedef ::SireCAS::Expression ( ::SireCAS::AlwaysFalse::*substitute_function_type)( ::SireCAS::Identities const & ) const;
            substitute_function_type substitute_function_value( &::SireCAS::AlwaysFalse::substitute );
            
            AlwaysFalse_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") )
                , bp::release_gil_policy()
                , "There is no substituting false" );
        
        }
        { //::SireCAS::AlwaysFalse::symbols
        
            typedef ::SireCAS::Symbols ( ::SireCAS::AlwaysFalse::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::AlwaysFalse::symbols );
            
            AlwaysFalse_exposer.def( 
                "symbols"
                , symbols_function_value
                , bp::release_gil_policy()
                , "There are no symbols in false" );
        
        }
        { //::SireCAS::AlwaysFalse::toString
        
            typedef ::QString ( ::SireCAS::AlwaysFalse::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireCAS::AlwaysFalse::toString );
            
            AlwaysFalse_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of truth" );
        
        }
        { //::SireCAS::AlwaysFalse::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::AlwaysFalse::typeName );
            
            AlwaysFalse_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::AlwaysFalse::what
        
            typedef char const * ( ::SireCAS::AlwaysFalse::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::AlwaysFalse::what );
            
            AlwaysFalse_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AlwaysFalse_exposer.staticmethod( "typeName" );
        AlwaysFalse_exposer.def( "__copy__", &__copy__);
        AlwaysFalse_exposer.def( "__deepcopy__", &__copy__);
        AlwaysFalse_exposer.def( "clone", &__copy__);
        AlwaysFalse_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::AlwaysFalse >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AlwaysFalse_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::AlwaysFalse >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AlwaysFalse_exposer.def_pickle(sire_pickle_suite< ::SireCAS::AlwaysFalse >());
        AlwaysFalse_exposer.def( "__str__", &__str__< ::SireCAS::AlwaysFalse > );
        AlwaysFalse_exposer.def( "__repr__", &__str__< ::SireCAS::AlwaysFalse > );
        AlwaysFalse_exposer.def( "__hash__", &::SireCAS::AlwaysFalse::hash );
    }

}
