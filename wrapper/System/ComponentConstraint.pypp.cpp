// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "ComponentConstraint.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/errors.h"

#include "constraint.h"

#include "delta.h"

#include "system.h"

#include <QDebug>

#include "constraint.h"

SireSystem::ComponentConstraint __copy__(const SireSystem::ComponentConstraint &other){ return SireSystem::ComponentConstraint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ComponentConstraint_class(){

    { //::SireSystem::ComponentConstraint
        typedef bp::class_< SireSystem::ComponentConstraint, bp::bases< SireSystem::Constraint, SireBase::Property > > ComponentConstraint_exposer_t;
        ComponentConstraint_exposer_t ComponentConstraint_exposer = ComponentConstraint_exposer_t( "ComponentConstraint", "This constraint is used to constrain the value of a\ncomponent of the system to a specific value, or the result\nof an expression based on other components in the system.\n\nYou can use this constraint, for example, to constrain\nthe value of lambda_forwards to equal Min( 1, lambda+delta_lambda )\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope ComponentConstraint_scope( ComponentConstraint_exposer );
        ComponentConstraint_exposer.def( bp::init< SireCAS::Symbol const &, SireCAS::Expression const & >(( bp::arg("component"), bp::arg("expression") ), "Construct to constrain the component with symbol component\nto have the value resulting from the expression expression") );
        ComponentConstraint_exposer.def( bp::init< SireSystem::ComponentConstraint const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::ComponentConstraint::component
        
            typedef ::SireCAS::Symbol const & ( ::SireSystem::ComponentConstraint::*component_function_type)(  ) const;
            component_function_type component_function_value( &::SireSystem::ComponentConstraint::component );
            
            ComponentConstraint_exposer.def( 
                "component"
                , component_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the symbol representing the component being constrained" );
        
        }
        { //::SireSystem::ComponentConstraint::expression
        
            typedef ::SireCAS::Expression const & ( ::SireSystem::ComponentConstraint::*expression_function_type)(  ) const;
            expression_function_type expression_function_value( &::SireSystem::ComponentConstraint::expression );
            
            ComponentConstraint_exposer.def( 
                "expression"
                , expression_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the expression used to evaluate the constraint" );
        
        }
        ComponentConstraint_exposer.def( bp::self != bp::self );
        { //::SireSystem::ComponentConstraint::operator=
        
            typedef ::SireSystem::ComponentConstraint & ( ::SireSystem::ComponentConstraint::*assign_function_type)( ::SireSystem::ComponentConstraint const & ) ;
            assign_function_type assign_function_value( &::SireSystem::ComponentConstraint::operator= );
            
            ComponentConstraint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ComponentConstraint_exposer.def( bp::self == bp::self );
        { //::SireSystem::ComponentConstraint::toString
        
            typedef ::QString ( ::SireSystem::ComponentConstraint::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireSystem::ComponentConstraint::toString );
            
            ComponentConstraint_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of the constraint" );
        
        }
        { //::SireSystem::ComponentConstraint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::ComponentConstraint::typeName );
            
            ComponentConstraint_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        ComponentConstraint_exposer.staticmethod( "typeName" );
        ComponentConstraint_exposer.def( "__copy__", &__copy__);
        ComponentConstraint_exposer.def( "__deepcopy__", &__copy__);
        ComponentConstraint_exposer.def( "clone", &__copy__);
        ComponentConstraint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::ComponentConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ComponentConstraint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::ComponentConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ComponentConstraint_exposer.def( "__str__", &__str__< ::SireSystem::ComponentConstraint > );
        ComponentConstraint_exposer.def( "__repr__", &__str__< ::SireSystem::ComponentConstraint > );
    }

}
