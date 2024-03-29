// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "MolecularDynamics.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/rangenerator.h"

#include "SireMol/moleculegroup.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "moleculardynamics.h"

#include "velocityverlet.h"

#include "moleculardynamics.h"

SireMove::MolecularDynamics __copy__(const SireMove::MolecularDynamics &other){ return SireMove::MolecularDynamics(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_MolecularDynamics_class(){

    { //::SireMove::MolecularDynamics
        typedef bp::class_< SireMove::MolecularDynamics, bp::bases< SireMove::Dynamics, SireMove::Move, SireBase::Property > > MolecularDynamics_exposer_t;
        MolecularDynamics_exposer_t MolecularDynamics_exposer = MolecularDynamics_exposer_t( "MolecularDynamics", "This class implements a molecular dynamics move.\n\nAuthor: Christopher Woods\n", bp::init< bp::optional< SireBase::PropertyMap const & > >(( bp::arg("map")=SireBase::PropertyMap() ), "Constructor") );
        bp::scope MolecularDynamics_scope( MolecularDynamics_exposer );
        MolecularDynamics_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to perform moves on the molecules in the group molgroup. This\ndefaults to an all-atom velocity-verlet integrator") );
        MolecularDynamics_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireMove::Integrator const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("integrator"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a move for the passed molecule group, integrated\nusing the supplied integrator") );
        MolecularDynamics_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireUnits::Dimension::Time, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("timestep"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a move for the passed molecule group, integrated with\nthe passed timestep") );
        MolecularDynamics_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireMove::Integrator const &, SireUnits::Dimension::Time, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("integrator"), bp::arg("timestep"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a move for the passed molecule group, integrated\nusing the passed integrator using the passed timestep") );
        MolecularDynamics_exposer.def( bp::init< SireMove::MolecularDynamics const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::MolecularDynamics::clearStatistics
        
            typedef void ( ::SireMove::MolecularDynamics::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireMove::MolecularDynamics::clearStatistics );
            
            MolecularDynamics_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , bp::release_gil_policy()
                , "Completely clear any move statistics - this clears all existing\nvelocities" );
        
        }
        { //::SireMove::MolecularDynamics::coordinatesProperty
        
            typedef ::SireBase::PropertyName ( ::SireMove::MolecularDynamics::*coordinatesProperty_function_type)(  ) const;
            coordinatesProperty_function_type coordinatesProperty_function_value( &::SireMove::MolecularDynamics::coordinatesProperty );
            
            MolecularDynamics_exposer.def( 
                "coordinatesProperty"
                , coordinatesProperty_function_value
                , bp::release_gil_policy()
                , "Return the property used to find the molecular coordinates" );
        
        }
        { //::SireMove::MolecularDynamics::elementsProperty
        
            typedef ::SireBase::PropertyName ( ::SireMove::MolecularDynamics::*elementsProperty_function_type)(  ) const;
            elementsProperty_function_type elementsProperty_function_value( &::SireMove::MolecularDynamics::elementsProperty );
            
            MolecularDynamics_exposer.def( 
                "elementsProperty"
                , elementsProperty_function_value
                , bp::release_gil_policy()
                , "Return the property used to find the atomic elements" );
        
        }
        { //::SireMove::MolecularDynamics::integrator
        
            typedef ::SireMove::Integrator const & ( ::SireMove::MolecularDynamics::*integrator_function_type)(  ) const;
            integrator_function_type integrator_function_value( &::SireMove::MolecularDynamics::integrator );
            
            MolecularDynamics_exposer.def( 
                "integrator"
                , integrator_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the integrator used to advance the coordinates\nfrom one timestep to the next" );
        
        }
        { //::SireMove::MolecularDynamics::kineticEnergy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::MolecularDynamics::*kineticEnergy_function_type)(  ) const;
            kineticEnergy_function_type kineticEnergy_function_value( &::SireMove::MolecularDynamics::kineticEnergy );
            
            MolecularDynamics_exposer.def( 
                "kineticEnergy"
                , kineticEnergy_function_value
                , bp::release_gil_policy()
                , "Return the kinetic energy of the system at the last move." );
        
        }
        { //::SireMove::MolecularDynamics::massesProperty
        
            typedef ::SireBase::PropertyName ( ::SireMove::MolecularDynamics::*massesProperty_function_type)(  ) const;
            massesProperty_function_type massesProperty_function_value( &::SireMove::MolecularDynamics::massesProperty );
            
            MolecularDynamics_exposer.def( 
                "massesProperty"
                , massesProperty_function_value
                , bp::release_gil_policy()
                , "Return the property used to find the molecular masses" );
        
        }
        { //::SireMove::MolecularDynamics::moleculeGroup
        
            typedef ::SireMol::MoleculeGroup const & ( ::SireMove::MolecularDynamics::*moleculeGroup_function_type)(  ) const;
            moleculeGroup_function_type moleculeGroup_function_value( &::SireMove::MolecularDynamics::moleculeGroup );
            
            MolecularDynamics_exposer.def( 
                "moleculeGroup"
                , moleculeGroup_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the molecule group on which this move operates" );
        
        }
        { //::SireMove::MolecularDynamics::move
        
            typedef void ( ::SireMove::MolecularDynamics::*move_function_type)( ::SireSystem::System &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::MolecularDynamics::move );
            
            MolecularDynamics_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Perform this move on the System system - perform the move\nnmoves times, optionally recording simulation statistics\nif record_stats is true" );
        
        }
        { //::SireMove::MolecularDynamics::nMoves
        
            typedef int ( ::SireMove::MolecularDynamics::*nMoves_function_type)(  ) const;
            nMoves_function_type nMoves_function_value( &::SireMove::MolecularDynamics::nMoves );
            
            MolecularDynamics_exposer.def( 
                "nMoves"
                , nMoves_function_value
                , bp::release_gil_policy()
                , "Return the number of moves completed using this object" );
        
        }
        MolecularDynamics_exposer.def( bp::self != bp::self );
        { //::SireMove::MolecularDynamics::operator=
        
            typedef ::SireMove::MolecularDynamics & ( ::SireMove::MolecularDynamics::*assign_function_type)( ::SireMove::MolecularDynamics const & ) ;
            assign_function_type assign_function_value( &::SireMove::MolecularDynamics::operator= );
            
            MolecularDynamics_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        MolecularDynamics_exposer.def( bp::self == bp::self );
        { //::SireMove::MolecularDynamics::regenerateVelocities
        
            typedef void ( ::SireMove::MolecularDynamics::*regenerateVelocities_function_type)( ::SireSystem::System const &,::SireMove::VelocityGenerator const & ) ;
            regenerateVelocities_function_type regenerateVelocities_function_value( &::SireMove::MolecularDynamics::regenerateVelocities );
            
            MolecularDynamics_exposer.def( 
                "regenerateVelocities"
                , regenerateVelocities_function_value
                , ( bp::arg("system"), bp::arg("generator") )
                , bp::release_gil_policy()
                , "Regenerate all of the velocities using the passed velocity generator" );
        
        }
        { //::SireMove::MolecularDynamics::setCoordinatesProperty
        
            typedef void ( ::SireMove::MolecularDynamics::*setCoordinatesProperty_function_type)( ::SireBase::PropertyName const & ) ;
            setCoordinatesProperty_function_type setCoordinatesProperty_function_value( &::SireMove::MolecularDynamics::setCoordinatesProperty );
            
            MolecularDynamics_exposer.def( 
                "setCoordinatesProperty"
                , setCoordinatesProperty_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the property used to find the coordinates of the molecules" );
        
        }
        { //::SireMove::MolecularDynamics::setElementsProperty
        
            typedef void ( ::SireMove::MolecularDynamics::*setElementsProperty_function_type)( ::SireBase::PropertyName const & ) ;
            setElementsProperty_function_type setElementsProperty_function_value( &::SireMove::MolecularDynamics::setElementsProperty );
            
            MolecularDynamics_exposer.def( 
                "setElementsProperty"
                , setElementsProperty_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the property used to find the elements of the atoms" );
        
        }
        { //::SireMove::MolecularDynamics::setGenerator
        
            typedef void ( ::SireMove::MolecularDynamics::*setGenerator_function_type)( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::MolecularDynamics::setGenerator );
            
            MolecularDynamics_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("generator") )
                , bp::release_gil_policy()
                , "Set the random number generator used by this move\n(this move may be completely deterministic, so may not\nuse a generator)" );
        
        }
        { //::SireMove::MolecularDynamics::setIntegrator
        
            typedef void ( ::SireMove::MolecularDynamics::*setIntegrator_function_type)( ::SireMove::Integrator const & ) ;
            setIntegrator_function_type setIntegrator_function_value( &::SireMove::MolecularDynamics::setIntegrator );
            
            MolecularDynamics_exposer.def( 
                "setIntegrator"
                , setIntegrator_function_value
                , ( bp::arg("integrator") )
                , bp::release_gil_policy()
                , "Set the integrator to be used to advance the coordinates from\none timestep to the next." );
        
        }
        { //::SireMove::MolecularDynamics::setMassesProperty
        
            typedef void ( ::SireMove::MolecularDynamics::*setMassesProperty_function_type)( ::SireBase::PropertyName const & ) ;
            setMassesProperty_function_type setMassesProperty_function_value( &::SireMove::MolecularDynamics::setMassesProperty );
            
            MolecularDynamics_exposer.def( 
                "setMassesProperty"
                , setMassesProperty_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the property used to find the molecular masses" );
        
        }
        { //::SireMove::MolecularDynamics::setMoleculeGroup
        
            typedef void ( ::SireMove::MolecularDynamics::*setMoleculeGroup_function_type)( ::SireMol::MoleculeGroup const & ) ;
            setMoleculeGroup_function_type setMoleculeGroup_function_value( &::SireMove::MolecularDynamics::setMoleculeGroup );
            
            MolecularDynamics_exposer.def( 
                "setMoleculeGroup"
                , setMoleculeGroup_function_value
                , ( bp::arg("molgroup") )
                , bp::release_gil_policy()
                , "Set the molecule group containing the molecules to be moved" );
        
        }
        { //::SireMove::MolecularDynamics::setMoleculeGroup
        
            typedef void ( ::SireMove::MolecularDynamics::*setMoleculeGroup_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) ;
            setMoleculeGroup_function_type setMoleculeGroup_function_value( &::SireMove::MolecularDynamics::setMoleculeGroup );
            
            MolecularDynamics_exposer.def( 
                "setMoleculeGroup"
                , setMoleculeGroup_function_value
                , ( bp::arg("molgroup"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Set the molecule group containing the molecules to be moved" );
        
        }
        { //::SireMove::MolecularDynamics::setSpaceProperty
        
            typedef void ( ::SireMove::MolecularDynamics::*setSpaceProperty_function_type)( ::SireBase::PropertyName const & ) ;
            setSpaceProperty_function_type setSpaceProperty_function_value( &::SireMove::MolecularDynamics::setSpaceProperty );
            
            MolecularDynamics_exposer.def( 
                "setSpaceProperty"
                , setSpaceProperty_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the property used to find the system space" );
        
        }
        { //::SireMove::MolecularDynamics::setTimeStep
        
            typedef void ( ::SireMove::MolecularDynamics::*setTimeStep_function_type)( ::SireUnits::Dimension::Time const & ) ;
            setTimeStep_function_type setTimeStep_function_value( &::SireMove::MolecularDynamics::setTimeStep );
            
            MolecularDynamics_exposer.def( 
                "setTimeStep"
                , setTimeStep_function_value
                , ( bp::arg("timestep") )
                , bp::release_gil_policy()
                , "Set the timestep for the dynamics integration" );
        
        }
        { //::SireMove::MolecularDynamics::setVelocitiesProperty
        
            typedef void ( ::SireMove::MolecularDynamics::*setVelocitiesProperty_function_type)( ::SireBase::PropertyName const & ) ;
            setVelocitiesProperty_function_type setVelocitiesProperty_function_value( &::SireMove::MolecularDynamics::setVelocitiesProperty );
            
            MolecularDynamics_exposer.def( 
                "setVelocitiesProperty"
                , setVelocitiesProperty_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the property used to find the molecular velocities" );
        
        }
        { //::SireMove::MolecularDynamics::setVelocityGeneratorProperty
        
            typedef void ( ::SireMove::MolecularDynamics::*setVelocityGeneratorProperty_function_type)( ::SireBase::PropertyName const & ) ;
            setVelocityGeneratorProperty_function_type setVelocityGeneratorProperty_function_value( &::SireMove::MolecularDynamics::setVelocityGeneratorProperty );
            
            MolecularDynamics_exposer.def( 
                "setVelocityGeneratorProperty"
                , setVelocityGeneratorProperty_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the property used to find the generator used to\ngenerate velocities when they are missing" );
        
        }
        { //::SireMove::MolecularDynamics::spaceProperty
        
            typedef ::SireBase::PropertyName ( ::SireMove::MolecularDynamics::*spaceProperty_function_type)(  ) const;
            spaceProperty_function_type spaceProperty_function_value( &::SireMove::MolecularDynamics::spaceProperty );
            
            MolecularDynamics_exposer.def( 
                "spaceProperty"
                , spaceProperty_function_value
                , bp::release_gil_policy()
                , "Return the property used to find the system space" );
        
        }
        { //::SireMove::MolecularDynamics::temperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::SireMove::MolecularDynamics::*temperature_function_type)(  ) const;
            temperature_function_type temperature_function_value( &::SireMove::MolecularDynamics::temperature );
            
            MolecularDynamics_exposer.def( 
                "temperature"
                , temperature_function_value
                , bp::release_gil_policy()
                , "Return the temperature of the system at the last move" );
        
        }
        { //::SireMove::MolecularDynamics::timeStep
        
            typedef ::SireUnits::Dimension::Time ( ::SireMove::MolecularDynamics::*timeStep_function_type)(  ) const;
            timeStep_function_type timeStep_function_value( &::SireMove::MolecularDynamics::timeStep );
            
            MolecularDynamics_exposer.def( 
                "timeStep"
                , timeStep_function_value
                , bp::release_gil_policy()
                , "Return the timestep for the integration" );
        
        }
        { //::SireMove::MolecularDynamics::toString
        
            typedef ::QString ( ::SireMove::MolecularDynamics::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::MolecularDynamics::toString );
            
            MolecularDynamics_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this move" );
        
        }
        { //::SireMove::MolecularDynamics::totalTime
        
            typedef ::SireUnits::Dimension::Time ( ::SireMove::MolecularDynamics::*totalTime_function_type)(  ) const;
            totalTime_function_type totalTime_function_value( &::SireMove::MolecularDynamics::totalTime );
            
            MolecularDynamics_exposer.def( 
                "totalTime"
                , totalTime_function_value
                , bp::release_gil_policy()
                , "Return the total amount of time simulated using these moves" );
        
        }
        { //::SireMove::MolecularDynamics::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::MolecularDynamics::typeName );
            
            MolecularDynamics_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::MolecularDynamics::velocitiesProperty
        
            typedef ::SireBase::PropertyName ( ::SireMove::MolecularDynamics::*velocitiesProperty_function_type)(  ) const;
            velocitiesProperty_function_type velocitiesProperty_function_value( &::SireMove::MolecularDynamics::velocitiesProperty );
            
            MolecularDynamics_exposer.def( 
                "velocitiesProperty"
                , velocitiesProperty_function_value
                , bp::release_gil_policy()
                , "Return the property used to find the molecular velocities" );
        
        }
        { //::SireMove::MolecularDynamics::velocityGeneratorProperty
        
            typedef ::SireBase::PropertyName ( ::SireMove::MolecularDynamics::*velocityGeneratorProperty_function_type)(  ) const;
            velocityGeneratorProperty_function_type velocityGeneratorProperty_function_value( &::SireMove::MolecularDynamics::velocityGeneratorProperty );
            
            MolecularDynamics_exposer.def( 
                "velocityGeneratorProperty"
                , velocityGeneratorProperty_function_value
                , bp::release_gil_policy()
                , "Return the property used to find the generator for\nmissing velocities" );
        
        }
        MolecularDynamics_exposer.staticmethod( "typeName" );
        MolecularDynamics_exposer.def( "__copy__", &__copy__);
        MolecularDynamics_exposer.def( "__deepcopy__", &__copy__);
        MolecularDynamics_exposer.def( "clone", &__copy__);
        MolecularDynamics_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::MolecularDynamics >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolecularDynamics_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::MolecularDynamics >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolecularDynamics_exposer.def_pickle(sire_pickle_suite< ::SireMove::MolecularDynamics >());
        MolecularDynamics_exposer.def( "__str__", &__str__< ::SireMove::MolecularDynamics > );
        MolecularDynamics_exposer.def( "__repr__", &__str__< ::SireMove::MolecularDynamics > );
    }

}
