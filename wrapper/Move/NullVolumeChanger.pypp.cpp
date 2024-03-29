// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "NullVolumeChanger.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/moleculegroup.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "math.h"

#include "sire_config.h"

#include "volumechanger.h"

#include <cmath>

#include "volumechanger.h"

SireMove::NullVolumeChanger __copy__(const SireMove::NullVolumeChanger &other){ return SireMove::NullVolumeChanger(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_NullVolumeChanger_class(){

    { //::SireMove::NullVolumeChanger
        typedef bp::class_< SireMove::NullVolumeChanger, bp::bases< SireMove::VolumeChanger, SireBase::Property > > NullVolumeChanger_exposer_t;
        NullVolumeChanger_exposer_t NullVolumeChanger_exposer = NullVolumeChanger_exposer_t( "NullVolumeChanger", "This is a null volume changer that does nothing", bp::init< >("Constructor") );
        bp::scope NullVolumeChanger_scope( NullVolumeChanger_exposer );
        NullVolumeChanger_exposer.def( bp::init< SireMove::NullVolumeChanger const & >(( bp::arg("other") ), "Copy constructor") );
        NullVolumeChanger_exposer.def( bp::self != bp::self );
        { //::SireMove::NullVolumeChanger::operator=
        
            typedef ::SireMove::NullVolumeChanger & ( ::SireMove::NullVolumeChanger::*assign_function_type)( ::SireMove::NullVolumeChanger const & ) ;
            assign_function_type assign_function_value( &::SireMove::NullVolumeChanger::operator= );
            
            NullVolumeChanger_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        NullVolumeChanger_exposer.def( bp::self == bp::self );
        { //::SireMove::NullVolumeChanger::setVolume
        
            typedef int ( ::SireMove::NullVolumeChanger::*setVolume_function_type)( ::SireSystem::System &,::SireUnits::Dimension::Volume const &,::SireBase::PropertyMap const & ) const;
            setVolume_function_type setVolume_function_value( &::SireMove::NullVolumeChanger::setVolume );
            
            NullVolumeChanger_exposer.def( 
                "setVolume"
                , setVolume_function_value
                , ( bp::arg("system"), bp::arg("volume"), bp::arg("map")=SireBase::PropertyMap() )
                , "The null volume changer doesnt change anything" );
        
        }
        { //::SireMove::NullVolumeChanger::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::NullVolumeChanger::typeName );
            
            NullVolumeChanger_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        NullVolumeChanger_exposer.staticmethod( "typeName" );
        NullVolumeChanger_exposer.def( "__copy__", &__copy__);
        NullVolumeChanger_exposer.def( "__deepcopy__", &__copy__);
        NullVolumeChanger_exposer.def( "clone", &__copy__);
        NullVolumeChanger_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::NullVolumeChanger >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullVolumeChanger_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::NullVolumeChanger >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullVolumeChanger_exposer.def_pickle(sire_pickle_suite< ::SireMove::NullVolumeChanger >());
        NullVolumeChanger_exposer.def( "__str__", &__str__< ::SireMove::NullVolumeChanger > );
        NullVolumeChanger_exposer.def( "__repr__", &__str__< ::SireMove::NullVolumeChanger > );
    }

}
