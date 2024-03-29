// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SupraSubMove.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "suprasubmove.h"

#include "suprasubsystem.h"

#include "suprasubmove.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_SupraSubMove_class(){

    { //::SireMove::SupraSubMove
        typedef bp::class_< SireMove::SupraSubMove, bp::bases< SireBase::Property >, boost::noncopyable > SupraSubMove_exposer_t;
        SupraSubMove_exposer_t SupraSubMove_exposer = SupraSubMove_exposer_t( "SupraSubMove", "This is the base class of the controller for the sub-moves\nthat are performed on the sub-systems of each SupraSystem\nas part of a SupraMove\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope SupraSubMove_scope( SupraSubMove_exposer );
        { //::SireMove::SupraSubMove::clearStatistics
        
            typedef void ( ::SireMove::SupraSubMove::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireMove::SupraSubMove::clearStatistics );
            
            SupraSubMove_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , bp::release_gil_policy()
                , "Clear the move statistics" );
        
        }
        { //::SireMove::SupraSubMove::move
        
            typedef void ( ::SireMove::SupraSubMove::*move_function_type)( ::SireMove::SupraSubSystem &,int,int,bool ) ;
            move_function_type move_function_value( &::SireMove::SupraSubMove::move );
            
            SupraSubMove_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("n_supra_moves"), bp::arg("n_supra_moves_per_block"), bp::arg("record_stats")=(bool)(true) )
                , "" );
        
        }
        { //::SireMove::SupraSubMove::null
        
            typedef ::SireMove::NullSupraSubMove const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMove::SupraSubMove::null );
            
            SupraSubMove_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the global null SupraSubMove" );
        
        }
        { //::SireMove::SupraSubMove::toString
        
            typedef ::QString ( ::SireMove::SupraSubMove::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::SupraSubMove::toString );
            
            SupraSubMove_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::SupraSubMove::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::SupraSubMove::typeName );
            
            SupraSubMove_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SupraSubMove_exposer.staticmethod( "null" );
        SupraSubMove_exposer.staticmethod( "typeName" );
        SupraSubMove_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::SupraSubMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SupraSubMove_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::SupraSubMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SupraSubMove_exposer.def_pickle(sire_pickle_suite< ::SireMove::SupraSubMove >());
        SupraSubMove_exposer.def( "__str__", &__str__< ::SireMove::SupraSubMove > );
        SupraSubMove_exposer.def( "__repr__", &__str__< ::SireMove::SupraSubMove > );
    }

}
