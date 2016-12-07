// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PowerFunction.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireError/errors.h"

#include "SireMaths/complex.h"

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "exp.h"

#include "identities.h"

#include "integrationconstant.h"

#include "power.h"

#include "powerconstant.h"

#include "values.h"

#include <QDebug>

#include "power.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_PowerFunction_class(){

    { //::SireCAS::PowerFunction
        typedef bp::class_< SireCAS::PowerFunction, bp::bases< SireCAS::ExBase >, boost::noncopyable > PowerFunction_exposer_t;
        PowerFunction_exposer_t PowerFunction_exposer = PowerFunction_exposer_t( "PowerFunction", "\nThis is the base class of all power expressions, e.g. x^y (all of the form core^power). There are several sub-classes that depend on exactly what is being raised to which power, e.g. Exp is e^y, Power is x^y, PowerConstant is c^y and ConstantPower is x^c (with ConstantPower further derived into RationalPower and RealPower based on whether the constant is rational). All of these can be constructed transparently by creating a Power and then calling reduce on the resulting object.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope PowerFunction_scope( PowerFunction_exposer );
        { //::SireCAS::PowerFunction::children
        
            typedef ::SireCAS::Expressions ( ::SireCAS::PowerFunction::*children_function_type)(  ) const;
            children_function_type children_function_value( &::SireCAS::PowerFunction::children );
            
            PowerFunction_exposer.def( 
                "children"
                , children_function_value
                , "Return the child expressions of this Power - this contains the core() and the power()" );
        
        }
        { //::SireCAS::PowerFunction::core
        
            typedef ::SireCAS::Expression ( ::SireCAS::PowerFunction::*core_function_type)(  ) const;
            core_function_type core_function_value( &::SireCAS::PowerFunction::core );
            
            PowerFunction_exposer.def( 
                "core"
                , core_function_value
                , "" );
        
        }
        { //::SireCAS::PowerFunction::differentiate
        
            typedef ::SireCAS::Expression ( ::SireCAS::PowerFunction::*differentiate_function_type)( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireCAS::PowerFunction::differentiate );
            
            PowerFunction_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") )
                , "Return the differential of this expression with respect to symbol" );
        
        }
        { //::SireCAS::PowerFunction::expand
        
            typedef ::QList< SireCAS::Factor > ( ::SireCAS::PowerFunction::*expand_function_type)( ::SireCAS::Symbol const & ) const;
            expand_function_type expand_function_value( &::SireCAS::PowerFunction::expand );
            
            PowerFunction_exposer.def( 
                "expand"
                , expand_function_value
                , ( bp::arg("symbol") )
                , "" );
        
        }
        { //::SireCAS::PowerFunction::functions
        
            typedef ::SireCAS::Functions ( ::SireCAS::PowerFunction::*functions_function_type)(  ) const;
            functions_function_type functions_function_value( &::SireCAS::PowerFunction::functions );
            
            PowerFunction_exposer.def( 
                "functions"
                , functions_function_value
                , "" );
        
        }
        { //::SireCAS::PowerFunction::integrate
        
            typedef ::SireCAS::Expression ( ::SireCAS::PowerFunction::*integrate_function_type)( ::SireCAS::Symbol const & ) const;
            integrate_function_type integrate_function_value( &::SireCAS::PowerFunction::integrate );
            
            PowerFunction_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("symbol") )
                , "Return the integral of this power with respect to symbol" );
        
        }
        { //::SireCAS::PowerFunction::isCompound
        
            typedef bool ( ::SireCAS::PowerFunction::*isCompound_function_type)(  ) const;
            isCompound_function_type isCompound_function_value( &::SireCAS::PowerFunction::isCompound );
            
            PowerFunction_exposer.def( 
                "isCompound"
                , isCompound_function_value
                , "" );
        
        }
        { //::SireCAS::PowerFunction::isConstant
        
            typedef bool ( ::SireCAS::PowerFunction::*isConstant_function_type)(  ) const;
            isConstant_function_type isConstant_function_value( &::SireCAS::PowerFunction::isConstant );
            
            PowerFunction_exposer.def( 
                "isConstant"
                , isConstant_function_value
                , "Return whether or not this is a constant" );
        
        }
        { //::SireCAS::PowerFunction::isFunction
        
            typedef bool ( ::SireCAS::PowerFunction::*isFunction_function_type)( ::SireCAS::Symbol const & ) const;
            isFunction_function_type isFunction_function_value( &::SireCAS::PowerFunction::isFunction );
            
            PowerFunction_exposer.def( 
                "isFunction"
                , isFunction_function_value
                , ( bp::arg("symbol") )
                , "Return whether this is a function of symbol" );
        
        }
        { //::SireCAS::PowerFunction::power
        
            typedef ::SireCAS::Expression ( ::SireCAS::PowerFunction::*power_function_type)(  ) const;
            power_function_type power_function_value( &::SireCAS::PowerFunction::power );
            
            PowerFunction_exposer.def( 
                "power"
                , power_function_value
                , "" );
        
        }
        { //::SireCAS::PowerFunction::reduce
        
            typedef ::SireCAS::Expression ( ::SireCAS::PowerFunction::*reduce_function_type)(  ) const;
            reduce_function_type reduce_function_value( &::SireCAS::PowerFunction::reduce );
            
            PowerFunction_exposer.def( 
                "reduce"
                , reduce_function_value
                , "Reduce this Power to a simplified expression (if possible)" );
        
        }
        { //::SireCAS::PowerFunction::substitute
        
            typedef ::SireCAS::Expression ( ::SireCAS::PowerFunction::*substitute_function_type)( ::SireCAS::Identities const & ) const;
            substitute_function_type substitute_function_value( &::SireCAS::PowerFunction::substitute );
            
            PowerFunction_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") )
                , "Return this expression with the supplied substitutions" );
        
        }
        { //::SireCAS::PowerFunction::symbols
        
            typedef ::SireCAS::Symbols ( ::SireCAS::PowerFunction::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::PowerFunction::symbols );
            
            PowerFunction_exposer.def( 
                "symbols"
                , symbols_function_value
                , "" );
        
        }
        { //::SireCAS::PowerFunction::toOpenMMString
        
            typedef ::QString ( ::SireCAS::PowerFunction::*toOpenMMString_function_type)(  ) const;
            toOpenMMString_function_type toOpenMMString_function_value( &::SireCAS::PowerFunction::toOpenMMString );
            
            PowerFunction_exposer.def( 
                "toOpenMMString"
                , toOpenMMString_function_value
                , "Return a string representation of this power in the OpenMM syntax" );
        
        }
        { //::SireCAS::PowerFunction::toString
        
            typedef ::QString ( ::SireCAS::PowerFunction::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireCAS::PowerFunction::toString );
            
            PowerFunction_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this power" );
        
        }
        { //::SireCAS::PowerFunction::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::PowerFunction::typeName );
            
            PowerFunction_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        PowerFunction_exposer.staticmethod( "typeName" );
        PowerFunction_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::PowerFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PowerFunction_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::PowerFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PowerFunction_exposer.def( "__str__", &__str__< ::SireCAS::PowerFunction > );
        PowerFunction_exposer.def( "__repr__", &__str__< ::SireCAS::PowerFunction > );
        PowerFunction_exposer.def( "__hash__", &::SireCAS::PowerFunction::hash );
    }

}
