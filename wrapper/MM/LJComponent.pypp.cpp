// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "LJComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "cljcomponent.h"

#include "cljcomponent.h"

SireMM::LJComponent __copy__(const SireMM::LJComponent &other){ return SireMM::LJComponent(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_LJComponent_class(){

    { //::SireMM::LJComponent
        typedef bp::class_< SireMM::LJComponent, bp::bases< SireFF::FFComponent, SireCAS::Symbol, SireCAS::ExBase > > LJComponent_exposer_t;
        LJComponent_exposer_t LJComponent_exposer = LJComponent_exposer_t( "LJComponent", "", bp::init< bp::optional< SireFF::FFName const & > >(( bp::arg("ffname")=SireFF::FFName() ), "Constructor") );
        bp::scope LJComponent_scope( LJComponent_exposer );
        LJComponent_exposer.def( bp::init< SireFF::FFName const &, QString const & >(( bp::arg("ffname"), bp::arg("suffix") ), "Construct using the name of the forcefield, and the passed suffix") );
        LJComponent_exposer.def( bp::init< SireCAS::Symbol const & >(( bp::arg("symbol") ), "Construct from a symbol\nThrow: SireError::incompatible_error\n") );
        LJComponent_exposer.def( bp::init< SireMM::LJComponent const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::LJComponent::changeEnergy
        
            typedef void ( ::SireMM::LJComponent::*changeEnergy_function_type)( ::SireFF::FF &,::SireMM::LJEnergy const & ) const;
            changeEnergy_function_type changeEnergy_function_value( &::SireMM::LJComponent::changeEnergy );
            
            LJComponent_exposer.def( 
                "changeEnergy"
                , changeEnergy_function_value
                , ( bp::arg("ff"), bp::arg("ljnrg") )
                , "Change the LJ component of the energy in the forcefield ff\nby delta" );
        
        }
        { //::SireMM::LJComponent::setEnergy
        
            typedef void ( ::SireMM::LJComponent::*setEnergy_function_type)( ::SireFF::FF &,::SireMM::LJEnergy const & ) const;
            setEnergy_function_type setEnergy_function_value( &::SireMM::LJComponent::setEnergy );
            
            LJComponent_exposer.def( 
                "setEnergy"
                , setEnergy_function_value
                , ( bp::arg("ff"), bp::arg("ljnrg") )
                , "Set the LJ component of the energy in the forcefield ff\nto equal to the passed LJEnergy" );
        
        }
        { //::SireMM::LJComponent::symbols
        
            typedef ::SireCAS::Symbols ( ::SireMM::LJComponent::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireMM::LJComponent::symbols );
            
            LJComponent_exposer.def( 
                "symbols"
                , symbols_function_value
                , "" );
        
        }
        { //::SireMM::LJComponent::total
        
            typedef ::SireMM::LJComponent const & ( ::SireMM::LJComponent::*total_function_type)(  ) const;
            total_function_type total_function_value( &::SireMM::LJComponent::total );
            
            LJComponent_exposer.def( 
                "total"
                , total_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireMM::LJComponent::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::LJComponent::typeName );
            
            LJComponent_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::LJComponent::what
        
            typedef char const * ( ::SireMM::LJComponent::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::LJComponent::what );
            
            LJComponent_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        LJComponent_exposer.staticmethod( "typeName" );
        LJComponent_exposer.def( "__copy__", &__copy__);
        LJComponent_exposer.def( "__deepcopy__", &__copy__);
        LJComponent_exposer.def( "clone", &__copy__);
        LJComponent_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::LJComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        LJComponent_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::LJComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        LJComponent_exposer.def( "__str__", &__str__< ::SireMM::LJComponent > );
        LJComponent_exposer.def( "__repr__", &__str__< ::SireMM::LJComponent > );
        LJComponent_exposer.def( "__hash__", &::SireMM::LJComponent::hash );
    }

}
