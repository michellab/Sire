// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Trajectory.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/space.h"

#include "trajectory.h"

#include "trajectory.h"

SireMol::Trajectory __copy__(const SireMol::Trajectory &other){ return SireMol::Trajectory(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Trajectory_class(){

    { //::SireMol::Trajectory
        typedef bp::class_< SireMol::Trajectory, bp::bases< SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property > > Trajectory_exposer_t;
        Trajectory_exposer_t Trajectory_exposer = Trajectory_exposer_t( "Trajectory", "This is a molecular property that holds the handle to the\ntrajectory data for that molecule. In addition to the\nhandle, this also holds the index of the first atom\nin the underlying trajectory data (trajectory data is a\nvector of coordinates in atom index order for each molecule)\n", bp::init< >("") );
        bp::scope Trajectory_scope( Trajectory_exposer );
        Trajectory_exposer.def( bp::init< SireMol::TrajectoryData const &, qint64 >(( bp::arg("data"), bp::arg("index") ), "") );
        Trajectory_exposer.def( bp::init< SireMol::Trajectory const & >(( bp::arg("other") ), "") );
        { //::SireMol::Trajectory::appendFrame
        
            typedef void ( ::SireMol::Trajectory::*appendFrame_function_type)( ::SireMol::AtomCoords const & ) ;
            appendFrame_function_type appendFrame_function_value( &::SireMol::Trajectory::appendFrame );
            
            Trajectory_exposer.def( 
                "appendFrame"
                , appendFrame_function_value
                , ( bp::arg("coords") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::appendFrame
        
            typedef void ( ::SireMol::Trajectory::*appendFrame_function_type)( ::SireMol::AtomVelocities const & ) ;
            appendFrame_function_type appendFrame_function_value( &::SireMol::Trajectory::appendFrame );
            
            Trajectory_exposer.def( 
                "appendFrame"
                , appendFrame_function_value
                , ( bp::arg("velocities") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::appendFrame
        
            typedef void ( ::SireMol::Trajectory::*appendFrame_function_type)( ::SireMol::AtomForces const & ) ;
            appendFrame_function_type appendFrame_function_value( &::SireMol::Trajectory::appendFrame );
            
            Trajectory_exposer.def( 
                "appendFrame"
                , appendFrame_function_value
                , ( bp::arg("forces") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::assertContainsCoordinates
        
            typedef void ( ::SireMol::Trajectory::*assertContainsCoordinates_function_type)(  ) const;
            assertContainsCoordinates_function_type assertContainsCoordinates_function_value( &::SireMol::Trajectory::assertContainsCoordinates );
            
            Trajectory_exposer.def( 
                "assertContainsCoordinates"
                , assertContainsCoordinates_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::assertContainsForces
        
            typedef void ( ::SireMol::Trajectory::*assertContainsForces_function_type)(  ) const;
            assertContainsForces_function_type assertContainsForces_function_value( &::SireMol::Trajectory::assertContainsForces );
            
            Trajectory_exposer.def( 
                "assertContainsForces"
                , assertContainsForces_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::assertContainsSpace
        
            typedef void ( ::SireMol::Trajectory::*assertContainsSpace_function_type)(  ) const;
            assertContainsSpace_function_type assertContainsSpace_function_value( &::SireMol::Trajectory::assertContainsSpace );
            
            Trajectory_exposer.def( 
                "assertContainsSpace"
                , assertContainsSpace_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::assertContainsVelocities
        
            typedef void ( ::SireMol::Trajectory::*assertContainsVelocities_function_type)(  ) const;
            assertContainsVelocities_function_type assertContainsVelocities_function_value( &::SireMol::Trajectory::assertContainsVelocities );
            
            Trajectory_exposer.def( 
                "assertContainsVelocities"
                , assertContainsVelocities_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::assertSupportsCoordinates
        
            typedef void ( ::SireMol::Trajectory::*assertSupportsCoordinates_function_type)(  ) const;
            assertSupportsCoordinates_function_type assertSupportsCoordinates_function_value( &::SireMol::Trajectory::assertSupportsCoordinates );
            
            Trajectory_exposer.def( 
                "assertSupportsCoordinates"
                , assertSupportsCoordinates_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::assertSupportsForces
        
            typedef void ( ::SireMol::Trajectory::*assertSupportsForces_function_type)(  ) const;
            assertSupportsForces_function_type assertSupportsForces_function_value( &::SireMol::Trajectory::assertSupportsForces );
            
            Trajectory_exposer.def( 
                "assertSupportsForces"
                , assertSupportsForces_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::assertSupportsSpace
        
            typedef void ( ::SireMol::Trajectory::*assertSupportsSpace_function_type)(  ) const;
            assertSupportsSpace_function_type assertSupportsSpace_function_value( &::SireMol::Trajectory::assertSupportsSpace );
            
            Trajectory_exposer.def( 
                "assertSupportsSpace"
                , assertSupportsSpace_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::assertSupportsVelocities
        
            typedef void ( ::SireMol::Trajectory::*assertSupportsVelocities_function_type)(  ) const;
            assertSupportsVelocities_function_type assertSupportsVelocities_function_value( &::SireMol::Trajectory::assertSupportsVelocities );
            
            Trajectory_exposer.def( 
                "assertSupportsVelocities"
                , assertSupportsVelocities_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::containsCoordinates
        
            typedef bool ( ::SireMol::Trajectory::*containsCoordinates_function_type)(  ) const;
            containsCoordinates_function_type containsCoordinates_function_value( &::SireMol::Trajectory::containsCoordinates );
            
            Trajectory_exposer.def( 
                "containsCoordinates"
                , containsCoordinates_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::containsForces
        
            typedef bool ( ::SireMol::Trajectory::*containsForces_function_type)(  ) const;
            containsForces_function_type containsForces_function_value( &::SireMol::Trajectory::containsForces );
            
            Trajectory_exposer.def( 
                "containsForces"
                , containsForces_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::containsSpace
        
            typedef bool ( ::SireMol::Trajectory::*containsSpace_function_type)(  ) const;
            containsSpace_function_type containsSpace_function_value( &::SireMol::Trajectory::containsSpace );
            
            Trajectory_exposer.def( 
                "containsSpace"
                , containsSpace_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::containsVelocities
        
            typedef bool ( ::SireMol::Trajectory::*containsVelocities_function_type)(  ) const;
            containsVelocities_function_type containsVelocities_function_value( &::SireMol::Trajectory::containsVelocities );
            
            Trajectory_exposer.def( 
                "containsVelocities"
                , containsVelocities_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::deleteFrame
        
            typedef void ( ::SireMol::Trajectory::*deleteFrame_function_type)( int ) ;
            deleteFrame_function_type deleteFrame_function_value( &::SireMol::Trajectory::deleteFrame );
            
            Trajectory_exposer.def( 
                "deleteFrame"
                , deleteFrame_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::getFrame
        
            typedef ::SireMol::AtomCoords ( ::SireMol::Trajectory::*getFrame_function_type)( int,::SireMol::AtomCoords const & ) const;
            getFrame_function_type getFrame_function_value( &::SireMol::Trajectory::getFrame );
            
            Trajectory_exposer.def( 
                "getFrame"
                , getFrame_function_value
                , ( bp::arg("i"), bp::arg("coords") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::getFrame
        
            typedef ::SireMol::AtomVelocities ( ::SireMol::Trajectory::*getFrame_function_type)( int,::SireMol::AtomVelocities const & ) const;
            getFrame_function_type getFrame_function_value( &::SireMol::Trajectory::getFrame );
            
            Trajectory_exposer.def( 
                "getFrame"
                , getFrame_function_value
                , ( bp::arg("i"), bp::arg("velocities") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::getFrame
        
            typedef ::SireMol::AtomForces ( ::SireMol::Trajectory::*getFrame_function_type)( int,::SireMol::AtomForces const & ) const;
            getFrame_function_type getFrame_function_value( &::SireMol::Trajectory::getFrame );
            
            Trajectory_exposer.def( 
                "getFrame"
                , getFrame_function_value
                , ( bp::arg("i"), bp::arg("forces") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::getSpace
        
            typedef ::SireVol::Space const & ( ::SireMol::Trajectory::*getSpace_function_type)( int ) const;
            getSpace_function_type getSpace_function_value( &::SireMol::Trajectory::getSpace );
            
            Trajectory_exposer.def( 
                "getSpace"
                , getSpace_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::Trajectory::insertFrame
        
            typedef void ( ::SireMol::Trajectory::*insertFrame_function_type)( int,::SireMol::AtomCoords const & ) ;
            insertFrame_function_type insertFrame_function_value( &::SireMol::Trajectory::insertFrame );
            
            Trajectory_exposer.def( 
                "insertFrame"
                , insertFrame_function_value
                , ( bp::arg("i"), bp::arg("coords") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::insertFrame
        
            typedef void ( ::SireMol::Trajectory::*insertFrame_function_type)( int,::SireMol::AtomVelocities const & ) ;
            insertFrame_function_type insertFrame_function_value( &::SireMol::Trajectory::insertFrame );
            
            Trajectory_exposer.def( 
                "insertFrame"
                , insertFrame_function_value
                , ( bp::arg("i"), bp::arg("velocities") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::insertFrame
        
            typedef void ( ::SireMol::Trajectory::*insertFrame_function_type)( int,::SireMol::AtomForces const & ) ;
            insertFrame_function_type insertFrame_function_value( &::SireMol::Trajectory::insertFrame );
            
            Trajectory_exposer.def( 
                "insertFrame"
                , insertFrame_function_value
                , ( bp::arg("i"), bp::arg("forces") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::isCompatibleWith
        
            typedef bool ( ::SireMol::Trajectory::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::Trajectory::isCompatibleWith );
            
            Trajectory_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::isEmpty
        
            typedef bool ( ::SireMol::Trajectory::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Trajectory::isEmpty );
            
            Trajectory_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::nFrames
        
            typedef int ( ::SireMol::Trajectory::*nFrames_function_type)(  ) const;
            nFrames_function_type nFrames_function_value( &::SireMol::Trajectory::nFrames );
            
            Trajectory_exposer.def( 
                "nFrames"
                , nFrames_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Trajectory_exposer.def( bp::self != bp::self );
        { //::SireMol::Trajectory::operator=
        
            typedef ::SireMol::Trajectory & ( ::SireMol::Trajectory::*assign_function_type)( ::SireMol::Trajectory const & ) ;
            assign_function_type assign_function_value( &::SireMol::Trajectory::operator= );
            
            Trajectory_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Trajectory_exposer.def( bp::self == bp::self );
        { //::SireMol::Trajectory::setFrame
        
            typedef void ( ::SireMol::Trajectory::*setFrame_function_type)( int,::SireMol::AtomCoords const & ) ;
            setFrame_function_type setFrame_function_value( &::SireMol::Trajectory::setFrame );
            
            Trajectory_exposer.def( 
                "setFrame"
                , setFrame_function_value
                , ( bp::arg("i"), bp::arg("coords") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::setFrame
        
            typedef void ( ::SireMol::Trajectory::*setFrame_function_type)( int,::SireMol::AtomVelocities const & ) ;
            setFrame_function_type setFrame_function_value( &::SireMol::Trajectory::setFrame );
            
            Trajectory_exposer.def( 
                "setFrame"
                , setFrame_function_value
                , ( bp::arg("i"), bp::arg("velocities") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::setFrame
        
            typedef void ( ::SireMol::Trajectory::*setFrame_function_type)( int,::SireMol::AtomForces const & ) ;
            setFrame_function_type setFrame_function_value( &::SireMol::Trajectory::setFrame );
            
            Trajectory_exposer.def( 
                "setFrame"
                , setFrame_function_value
                , ( bp::arg("i"), bp::arg("forces") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::supportsCoordinates
        
            typedef bool ( ::SireMol::Trajectory::*supportsCoordinates_function_type)(  ) const;
            supportsCoordinates_function_type supportsCoordinates_function_value( &::SireMol::Trajectory::supportsCoordinates );
            
            Trajectory_exposer.def( 
                "supportsCoordinates"
                , supportsCoordinates_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::supportsForces
        
            typedef bool ( ::SireMol::Trajectory::*supportsForces_function_type)(  ) const;
            supportsForces_function_type supportsForces_function_value( &::SireMol::Trajectory::supportsForces );
            
            Trajectory_exposer.def( 
                "supportsForces"
                , supportsForces_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::supportsSpace
        
            typedef bool ( ::SireMol::Trajectory::*supportsSpace_function_type)(  ) const;
            supportsSpace_function_type supportsSpace_function_value( &::SireMol::Trajectory::supportsSpace );
            
            Trajectory_exposer.def( 
                "supportsSpace"
                , supportsSpace_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::supportsVelocities
        
            typedef bool ( ::SireMol::Trajectory::*supportsVelocities_function_type)(  ) const;
            supportsVelocities_function_type supportsVelocities_function_value( &::SireMol::Trajectory::supportsVelocities );
            
            Trajectory_exposer.def( 
                "supportsVelocities"
                , supportsVelocities_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::timestep
        
            typedef ::SireUnits::Dimension::Time ( ::SireMol::Trajectory::*timestep_function_type)(  ) const;
            timestep_function_type timestep_function_value( &::SireMol::Trajectory::timestep );
            
            Trajectory_exposer.def( 
                "timestep"
                , timestep_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Trajectory::typeName );
            
            Trajectory_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Trajectory::what
        
            typedef char const * ( ::SireMol::Trajectory::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::Trajectory::what );
            
            Trajectory_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Trajectory_exposer.staticmethod( "typeName" );
        Trajectory_exposer.def( "__copy__", &__copy__);
        Trajectory_exposer.def( "__deepcopy__", &__copy__);
        Trajectory_exposer.def( "clone", &__copy__);
        Trajectory_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Trajectory >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Trajectory_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Trajectory >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Trajectory_exposer.def_pickle(sire_pickle_suite< ::SireMol::Trajectory >());
        Trajectory_exposer.def( "__str__", &__str__< ::SireMol::Trajectory > );
        Trajectory_exposer.def( "__repr__", &__str__< ::SireMol::Trajectory > );
    }

}