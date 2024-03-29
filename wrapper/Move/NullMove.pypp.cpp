// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "NullMove.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/rangenerator.h"

#include "SireMol/core.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "ensemble.h"

#include "move.h"

#include <QMutex>

#include "move.h"

SireMove::NullMove __copy__(const SireMove::NullMove &other){ return SireMove::NullMove(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_NullMove_class(){

    { //::SireMove::NullMove
        typedef bp::class_< SireMove::NullMove, bp::bases< SireMove::Move, SireBase::Property > > NullMove_exposer_t;
        NullMove_exposer_t NullMove_exposer = NullMove_exposer_t( "NullMove", "This is a null move - it doesnt change the system at all", bp::init< >("Constructor") );
        bp::scope NullMove_scope( NullMove_exposer );
        NullMove_exposer.def( bp::init< SireMove::NullMove const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::NullMove::clearStatistics
        
            typedef void ( ::SireMove::NullMove::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireMove::NullMove::clearStatistics );
            
            NullMove_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , bp::release_gil_policy()
                , "There are no statistics to clear" );
        
        }
        { //::SireMove::NullMove::ensemble
        
            typedef ::SireMove::Ensemble ( ::SireMove::NullMove::*ensemble_function_type)(  ) const;
            ensemble_function_type ensemble_function_value( &::SireMove::NullMove::ensemble );
            
            NullMove_exposer.def( 
                "ensemble"
                , ensemble_function_value
                , bp::release_gil_policy()
                , "NullMove doesnt change anything (so must be NVE)" );
        
        }
        { //::SireMove::NullMove::move
        
            typedef void ( ::SireMove::NullMove::*move_function_type)( ::SireSystem::System &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::NullMove::move );
            
            NullMove_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats") )
                , bp::release_gil_policy()
                , "NullMove doesnt perform any moves - no matter how hard you try" );
        
        }
        { //::SireMove::NullMove::nMoves
        
            typedef int ( ::SireMove::NullMove::*nMoves_function_type)(  ) const;
            nMoves_function_type nMoves_function_value( &::SireMove::NullMove::nMoves );
            
            NullMove_exposer.def( 
                "nMoves"
                , nMoves_function_value
                , bp::release_gil_policy()
                , "There have been and never will be any NullMove events" );
        
        }
        NullMove_exposer.def( bp::self != bp::self );
        { //::SireMove::NullMove::operator=
        
            typedef ::SireMove::NullMove & ( ::SireMove::NullMove::*assign_function_type)( ::SireMove::NullMove const & ) ;
            assign_function_type assign_function_value( &::SireMove::NullMove::operator= );
            
            NullMove_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        NullMove_exposer.def( bp::self == bp::self );
        { //::SireMove::NullMove::setGenerator
        
            typedef void ( ::SireMove::NullMove::*setGenerator_function_type)( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::NullMove::setGenerator );
            
            NullMove_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("rangenerator") )
                , bp::release_gil_policy()
                , "The NullMove does not use a random number generator" );
        
        }
        { //::SireMove::NullMove::toString
        
            typedef ::QString ( ::SireMove::NullMove::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::NullMove::toString );
            
            NullMove_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation" );
        
        }
        { //::SireMove::NullMove::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::NullMove::typeName );
            
            NullMove_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        NullMove_exposer.staticmethod( "typeName" );
        NullMove_exposer.def( "__copy__", &__copy__);
        NullMove_exposer.def( "__deepcopy__", &__copy__);
        NullMove_exposer.def( "clone", &__copy__);
        NullMove_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::NullMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullMove_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::NullMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullMove_exposer.def_pickle(sire_pickle_suite< ::SireMove::NullMove >());
        NullMove_exposer.def( "__str__", &__str__< ::SireMove::NullMove > );
        NullMove_exposer.def( "__repr__", &__str__< ::SireMove::NullMove > );
    }

}
