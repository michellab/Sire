// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "NullInserter.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/quaternion.h"

#include "SireMol/core.h"

#include "SireMol/molecule.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "molinserter.h"

#include "molinserter.h"

SireMove::NullInserter __copy__(const SireMove::NullInserter &other){ return SireMove::NullInserter(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_NullInserter_class(){

    { //::SireMove::NullInserter
        typedef bp::class_< SireMove::NullInserter, bp::bases< SireMove::MolInserter, SireBase::Property > > NullInserter_exposer_t;
        NullInserter_exposer_t NullInserter_exposer = NullInserter_exposer_t( "NullInserter", "This is the null inserter - this doesnt insert anything\ninto anything\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope NullInserter_scope( NullInserter_exposer );
        NullInserter_exposer.def( bp::init< SireMove::NullInserter const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::NullInserter::insert
        
            typedef double ( ::SireMove::NullInserter::*insert_function_type)( ::SireMol::Molecule const &,::SireSystem::System &,::SireVol::Space const & ) ;
            insert_function_type insert_function_value( &::SireMove::NullInserter::insert );
            
            NullInserter_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("molecule"), bp::arg("system"), bp::arg("space") )
                , bp::release_gil_policy()
                , "This does nothing" );
        
        }
        { //::SireMove::NullInserter::insert
        
            typedef double ( ::SireMove::NullInserter::*insert_function_type)( ::SireMol::PartialMolecule const &,::SireSystem::System &,::SireVol::Space const & ) ;
            insert_function_type insert_function_value( &::SireMove::NullInserter::insert );
            
            NullInserter_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("molecule"), bp::arg("system"), bp::arg("space") )
                , bp::release_gil_policy()
                , "This does nothing" );
        
        }
        NullInserter_exposer.def( bp::self != bp::self );
        { //::SireMove::NullInserter::operator=
        
            typedef ::SireMove::NullInserter & ( ::SireMove::NullInserter::*assign_function_type)( ::SireMove::NullInserter const & ) ;
            assign_function_type assign_function_value( &::SireMove::NullInserter::operator= );
            
            NullInserter_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        NullInserter_exposer.def( bp::self == bp::self );
        { //::SireMove::NullInserter::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::NullInserter::typeName );
            
            NullInserter_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        NullInserter_exposer.staticmethod( "typeName" );
        NullInserter_exposer.def( "__copy__", &__copy__);
        NullInserter_exposer.def( "__deepcopy__", &__copy__);
        NullInserter_exposer.def( "clone", &__copy__);
        NullInserter_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::NullInserter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullInserter_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::NullInserter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullInserter_exposer.def_pickle(sire_pickle_suite< ::SireMove::NullInserter >());
        NullInserter_exposer.def( "__str__", &__str__< ::SireMove::NullInserter > );
        NullInserter_exposer.def( "__repr__", &__str__< ::SireMove::NullInserter > );
    }

}
