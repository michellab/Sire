// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "NullDeleter.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "moldeleter.h"

#include "uniformsampler.h"

#include "moldeleter.h"

SireMove::NullDeleter __copy__(const SireMove::NullDeleter &other){ return SireMove::NullDeleter(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_NullDeleter_class(){

    { //::SireMove::NullDeleter
        typedef bp::class_< SireMove::NullDeleter, bp::bases< SireMove::MolDeleter, SireBase::Property > > NullDeleter_exposer_t;
        NullDeleter_exposer_t NullDeleter_exposer = NullDeleter_exposer_t( "NullDeleter", "This is a null deleter - this deletes nothing", bp::init< >("Constructor") );
        bp::scope NullDeleter_scope( NullDeleter_exposer );
        NullDeleter_exposer.def( bp::init< SireMove::NullDeleter const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::NullDeleter::deleteFrom
        
            typedef ::boost::tuples::tuple< SireMol::Molecule, double, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireMove::NullDeleter::*deleteFrom_function_type)( ::SireSystem::System & ) ;
            deleteFrom_function_type deleteFrom_function_value( &::SireMove::NullDeleter::deleteFrom );
            
            NullDeleter_exposer.def( 
                "deleteFrom"
                , deleteFrom_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Delete a molecule from the system - well this does nothing too" );
        
        }
        { //::SireMove::NullDeleter::generator
        
            typedef ::SireMaths::RanGenerator const & ( ::SireMove::NullDeleter::*generator_function_type)(  ) const;
            generator_function_type generator_function_value( &::SireMove::NullDeleter::generator );
            
            NullDeleter_exposer.def( 
                "generator"
                , generator_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the generator used to select random molecules\n(this just returns the global generator)" );
        
        }
        NullDeleter_exposer.def( bp::self != bp::self );
        { //::SireMove::NullDeleter::operator=
        
            typedef ::SireMove::NullDeleter & ( ::SireMove::NullDeleter::*assign_function_type)( ::SireMove::NullDeleter const & ) ;
            assign_function_type assign_function_value( &::SireMove::NullDeleter::operator= );
            
            NullDeleter_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        NullDeleter_exposer.def( bp::self == bp::self );
        { //::SireMove::NullDeleter::setGenerator
        
            typedef void ( ::SireMove::NullDeleter::*setGenerator_function_type)( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::NullDeleter::setGenerator );
            
            NullDeleter_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("generator") )
                , bp::release_gil_policy()
                , "Set the generator used to select random molecules\n(this does nothing)" );
        
        }
        { //::SireMove::NullDeleter::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::NullDeleter::typeName );
            
            NullDeleter_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        NullDeleter_exposer.staticmethod( "typeName" );
        NullDeleter_exposer.def( "__copy__", &__copy__);
        NullDeleter_exposer.def( "__deepcopy__", &__copy__);
        NullDeleter_exposer.def( "clone", &__copy__);
        NullDeleter_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::NullDeleter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullDeleter_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::NullDeleter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullDeleter_exposer.def_pickle(sire_pickle_suite< ::SireMove::NullDeleter >());
        NullDeleter_exposer.def( "__str__", &__str__< ::SireMove::NullDeleter > );
        NullDeleter_exposer.def( "__repr__", &__str__< ::SireMove::NullDeleter > );
    }

}
