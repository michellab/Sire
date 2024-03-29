// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "TitrationMove.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "titrationmove.h"

#include "titrator.h"

#include <QDebug>

#include "titrationmove.h"

SireMove::TitrationMove __copy__(const SireMove::TitrationMove &other){ return SireMove::TitrationMove(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_TitrationMove_class(){

    { //::SireMove::TitrationMove
        typedef bp::class_< SireMove::TitrationMove, bp::bases< SireMove::MonteCarlo, SireMove::Move, SireBase::Property > > TitrationMove_exposer_t;
        TitrationMove_exposer_t TitrationMove_exposer = TitrationMove_exposer_t( "TitrationMove", "This class performs a Monte Carlo titration move. This moves\na charge from one place to another by swapping the coordinates\nof once molecule with another, e.g. swapping a charge with a water.\nThis allows ions to move quickly through a simulation box, and for\nions to equilibrate between boxes (e.g. during a WSRC calcualtion)\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope TitrationMove_scope( TitrationMove_exposer );
        TitrationMove_exposer.def( bp::init< SireMove::TitrationMove const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::TitrationMove::move
        
            typedef void ( ::SireMove::TitrationMove::*move_function_type)( ::SireSystem::System &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::TitrationMove::move );
            
            TitrationMove_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Actually perform the move" );
        
        }
        TitrationMove_exposer.def( bp::self != bp::self );
        { //::SireMove::TitrationMove::operator=
        
            typedef ::SireMove::TitrationMove & ( ::SireMove::TitrationMove::*assign_function_type)( ::SireMove::TitrationMove const & ) ;
            assign_function_type assign_function_value( &::SireMove::TitrationMove::operator= );
            
            TitrationMove_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        TitrationMove_exposer.def( bp::self == bp::self );
        { //::SireMove::TitrationMove::toString
        
            typedef ::QString ( ::SireMove::TitrationMove::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::TitrationMove::toString );
            
            TitrationMove_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of the move" );
        
        }
        { //::SireMove::TitrationMove::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::TitrationMove::typeName );
            
            TitrationMove_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::TitrationMove::what
        
            typedef char const * ( ::SireMove::TitrationMove::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMove::TitrationMove::what );
            
            TitrationMove_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        TitrationMove_exposer.staticmethod( "typeName" );
        TitrationMove_exposer.def( "__copy__", &__copy__);
        TitrationMove_exposer.def( "__deepcopy__", &__copy__);
        TitrationMove_exposer.def( "clone", &__copy__);
        TitrationMove_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::TitrationMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TitrationMove_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::TitrationMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TitrationMove_exposer.def_pickle(sire_pickle_suite< ::SireMove::TitrationMove >());
        TitrationMove_exposer.def( "__str__", &__str__< ::SireMove::TitrationMove > );
        TitrationMove_exposer.def( "__repr__", &__str__< ::SireMove::TitrationMove > );
    }

}
