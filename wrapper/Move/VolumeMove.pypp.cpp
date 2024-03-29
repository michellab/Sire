// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "VolumeMove.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/savestate.h"

#include "SireError/errors.h"

#include "SireMol/mgname.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "volumemove.h"

#include "volumemove.h"

SireMove::VolumeMove __copy__(const SireMove::VolumeMove &other){ return SireMove::VolumeMove(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_VolumeMove_class(){

    { //::SireMove::VolumeMove
        typedef bp::class_< SireMove::VolumeMove, bp::bases< SireMove::MonteCarlo, SireMove::Move, SireBase::Property > > VolumeMove_exposer_t;
        VolumeMove_exposer_t VolumeMove_exposer = VolumeMove_exposer_t( "VolumeMove", "This is a Monte Carlo volume move. This is used to allow\nthe pressure to be kept constant\n\nAuthor: Christopher Woods\n", bp::init< bp::optional< SireBase::PropertyMap const & > >(( bp::arg("map")=SireBase::PropertyMap() ), "Null constructor") );
        bp::scope VolumeMove_scope( VolumeMove_exposer );
        VolumeMove_exposer.def( bp::init< SireMol::MGID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("mgid"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a volume move that can be used to generate the ensemble\nfor a temperature of 25 C, pressure of 1 atm, and with a maximum\nchange of 100 A^3 by moving the molecules in the\nmolecule groups that match the ID mgid\nusing a ScaleVolumeFromCenter centered on the origin") );
        VolumeMove_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a volume move that can be used to generate the ensemble\nfor a temperature of 25 C, pressure of 1 atm, and with a maximum\nchange of 100 A^3 by moving the molecules in molgroup\nusing a ScaleVolumeFromCenter centered on the origin") );
        VolumeMove_exposer.def( bp::init< SireMove::VolumeChanger const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("volchanger"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a volume move that can be used to generate the ensemble\nfor a temperature of 25 C, pressure of 1 atm, and with a maximum\nchange of 100 A^3 using the passed volume changer") );
        VolumeMove_exposer.def( bp::init< SireMove::VolumeMove const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::VolumeMove::groupID
        
            typedef ::SireMol::MGID const & ( ::SireMove::VolumeMove::*groupID_function_type)(  ) const;
            groupID_function_type groupID_function_value( &::SireMove::VolumeMove::groupID );
            
            VolumeMove_exposer.def( 
                "groupID"
                , groupID_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the ID that matches the molecule groups that\nwill be affected by this move" );
        
        }
        { //::SireMove::VolumeMove::maximumVolumeChange
        
            typedef ::SireUnits::Dimension::Volume const & ( ::SireMove::VolumeMove::*maximumVolumeChange_function_type)(  ) const;
            maximumVolumeChange_function_type maximumVolumeChange_function_value( &::SireMove::VolumeMove::maximumVolumeChange );
            
            VolumeMove_exposer.def( 
                "maximumVolumeChange"
                , maximumVolumeChange_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the maximum change of volume attempted by a move" );
        
        }
        { //::SireMove::VolumeMove::move
        
            typedef void ( ::SireMove::VolumeMove::*move_function_type)( ::SireSystem::System &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::VolumeMove::move );
            
            VolumeMove_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Perform nmoves volume moves on the passed system, optionally\nrecording simulation statistics if record_stats is true" );
        
        }
        VolumeMove_exposer.def( bp::self != bp::self );
        { //::SireMove::VolumeMove::operator=
        
            typedef ::SireMove::VolumeMove & ( ::SireMove::VolumeMove::*assign_function_type)( ::SireMove::VolumeMove const & ) ;
            assign_function_type assign_function_value( &::SireMove::VolumeMove::operator= );
            
            VolumeMove_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        VolumeMove_exposer.def( bp::self == bp::self );
        { //::SireMove::VolumeMove::setGenerator
        
            typedef void ( ::SireMove::VolumeMove::*setGenerator_function_type)( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::VolumeMove::setGenerator );
            
            VolumeMove_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("rangenerator") )
                , bp::release_gil_policy()
                , "Set the random number generator used by this move" );
        
        }
        { //::SireMove::VolumeMove::setMaximumVolumeChange
        
            typedef void ( ::SireMove::VolumeMove::*setMaximumVolumeChange_function_type)( ::SireUnits::Dimension::Volume const & ) ;
            setMaximumVolumeChange_function_type setMaximumVolumeChange_function_value( &::SireMove::VolumeMove::setMaximumVolumeChange );
            
            VolumeMove_exposer.def( 
                "setMaximumVolumeChange"
                , setMaximumVolumeChange_function_value
                , ( bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the maximum change in volume" );
        
        }
        { //::SireMove::VolumeMove::setVolumeChanger
        
            typedef void ( ::SireMove::VolumeMove::*setVolumeChanger_function_type)( ::SireMove::VolumeChanger const & ) ;
            setVolumeChanger_function_type setVolumeChanger_function_value( &::SireMove::VolumeMove::setVolumeChanger );
            
            VolumeMove_exposer.def( 
                "setVolumeChanger"
                , setVolumeChanger_function_value
                , ( bp::arg("volchanger") )
                , bp::release_gil_policy()
                , "Set the volume changer used to change the volume to volchanger" );
        
        }
        { //::SireMove::VolumeMove::setVolumeChanger
        
            typedef void ( ::SireMove::VolumeMove::*setVolumeChanger_function_type)( ::SireMol::MoleculeGroup const & ) ;
            setVolumeChanger_function_type setVolumeChanger_function_value( &::SireMove::VolumeMove::setVolumeChanger );
            
            VolumeMove_exposer.def( 
                "setVolumeChanger"
                , setVolumeChanger_function_value
                , ( bp::arg("molgroup") )
                , bp::release_gil_policy()
                , "Set the volume changer used to change the volume to a\nScaleVolumeFromCenter that scales the molecules in molgroup\nfrom the center of a box centered at (0,0,0)" );
        
        }
        { //::SireMove::VolumeMove::toString
        
            typedef ::QString ( ::SireMove::VolumeMove::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::VolumeMove::toString );
            
            VolumeMove_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this move" );
        
        }
        { //::SireMove::VolumeMove::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::VolumeMove::typeName );
            
            VolumeMove_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::VolumeMove::volumeChanger
        
            typedef ::SireMove::VolumeChanger const & ( ::SireMove::VolumeMove::*volumeChanger_function_type)(  ) const;
            volumeChanger_function_type volumeChanger_function_value( &::SireMove::VolumeMove::volumeChanger );
            
            VolumeMove_exposer.def( 
                "volumeChanger"
                , volumeChanger_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the volume changer used to change the volume" );
        
        }
        VolumeMove_exposer.staticmethod( "typeName" );
        VolumeMove_exposer.def( "__copy__", &__copy__);
        VolumeMove_exposer.def( "__deepcopy__", &__copy__);
        VolumeMove_exposer.def( "clone", &__copy__);
        VolumeMove_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::VolumeMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        VolumeMove_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::VolumeMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        VolumeMove_exposer.def_pickle(sire_pickle_suite< ::SireMove::VolumeMove >());
        VolumeMove_exposer.def( "__str__", &__str__< ::SireMove::VolumeMove > );
        VolumeMove_exposer.def( "__repr__", &__str__< ::SireMove::VolumeMove > );
    }

}
