// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "NullQM.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMol/molecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "latticecharges.h"

#include "qmprogram.h"

#include <QMutex>

#include "qmprogram.h"

Squire::NullQM __copy__(const Squire::NullQM &other){ return Squire::NullQM(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_NullQM_class(){

    { //::Squire::NullQM
        typedef bp::class_< Squire::NullQM, bp::bases< Squire::QMProgram, SireBase::Property > > NullQM_exposer_t;
        NullQM_exposer_t NullQM_exposer = NullQM_exposer_t( "NullQM", "This is the null QM program that returns zero energy and force", bp::init< >("Constructor") );
        bp::scope NullQM_scope( NullQM_exposer );
        NullQM_exposer.def( bp::init< Squire::NullQM const & >(( bp::arg("other") ), "Copy constructor") );
        NullQM_exposer.def( bp::self != bp::self );
        { //::Squire::NullQM::operator=
        
            typedef ::Squire::NullQM & ( ::Squire::NullQM::*assign_function_type)( ::Squire::NullQM const & ) ;
            assign_function_type assign_function_value( &::Squire::NullQM::operator= );
            
            NullQM_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        NullQM_exposer.def( bp::self == bp::self );
        { //::Squire::NullQM::supportsLatticeCharges
        
            typedef bool ( ::Squire::NullQM::*supportsLatticeCharges_function_type)(  ) const;
            supportsLatticeCharges_function_type supportsLatticeCharges_function_value( &::Squire::NullQM::supportsLatticeCharges );
            
            NullQM_exposer.def( 
                "supportsLatticeCharges"
                , supportsLatticeCharges_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::Squire::NullQM::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::NullQM::typeName );
            
            NullQM_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        NullQM_exposer.staticmethod( "typeName" );
        NullQM_exposer.def( "__copy__", &__copy__);
        NullQM_exposer.def( "__deepcopy__", &__copy__);
        NullQM_exposer.def( "clone", &__copy__);
        NullQM_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::NullQM >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullQM_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::NullQM >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullQM_exposer.def_pickle(sire_pickle_suite< ::Squire::NullQM >());
        NullQM_exposer.def( "__str__", &__str__< ::Squire::NullQM > );
        NullQM_exposer.def( "__repr__", &__str__< ::Squire::NullQM > );
    }

}
