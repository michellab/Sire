// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "OpenMMMDIntegrator.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/variantproperty.h"

#include "SireFF/forcetable.h"

#include "SireIO/amber.h"

#include "SireMM/atomljs.h"

#include "SireMM/internalff.h"

#include "SireMaths/constants.h"

#include "SireMaths/rangenerator.h"

#include "SireMaths/vector.h"

#include "SireMol/amberparameters.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atommasses.h"

#include "SireMol/bondid.h"

#include "SireMol/connectivity.h"

#include "SireMol/mgname.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleditor.h"

#include "SireMol/partialmolecule.h"

#include "SireMove/flexibility.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/convert.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "ensemble.h"

#include "openmmmdintegrator.h"

#include <QDebug>

#include <QTime>

#include <iostream>

#include "openmmmdintegrator.h"

SireMove::OpenMMMDIntegrator __copy__(const SireMove::OpenMMMDIntegrator &other){ return SireMove::OpenMMMDIntegrator(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_OpenMMMDIntegrator_class(){

    { //::SireMove::OpenMMMDIntegrator
        typedef bp::class_< SireMove::OpenMMMDIntegrator, bp::bases< SireMove::Integrator, SireBase::Property > > OpenMMMDIntegrator_exposer_t;
        OpenMMMDIntegrator_exposer_t OpenMMMDIntegrator_exposer = OpenMMMDIntegrator_exposer_t( "OpenMMMDIntegrator", "This class implements a pure MD integrator using OpenMM.\nNo free energy methods are supported.\n\nAuthor: Julien Michel and Gaetano Calabro\n", bp::init< bp::optional< bool > >(( bp::arg("frequent_save_velocities")=(bool)(false) ), "Constructor") );
        bp::scope OpenMMMDIntegrator_scope( OpenMMMDIntegrator_exposer );
        OpenMMMDIntegrator_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< bool > >(( bp::arg("molecule_group"), bp::arg("frequent_save_velocities")=(bool)(false) ), "Constructor using the passed molecule group") );
        OpenMMMDIntegrator_exposer.def( bp::init< SireMove::OpenMMMDIntegrator const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::OpenMMMDIntegrator::createWorkspace
        
            typedef ::SireMove::IntegratorWorkspacePtr ( ::SireMove::OpenMMMDIntegrator::*createWorkspace_function_type)( ::SireBase::PropertyMap const & ) const;
            createWorkspace_function_type createWorkspace_function_value( &::SireMove::OpenMMMDIntegrator::createWorkspace );
            
            OpenMMMDIntegrator_exposer.def( 
                "createWorkspace"
                , createWorkspace_function_value
                , ( bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Create an empty workspace" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::createWorkspace
        
            typedef ::SireMove::IntegratorWorkspacePtr ( ::SireMove::OpenMMMDIntegrator::*createWorkspace_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) const;
            createWorkspace_function_type createWorkspace_function_value( &::SireMove::OpenMMMDIntegrator::createWorkspace );
            
            OpenMMMDIntegrator_exposer.def( 
                "createWorkspace"
                , createWorkspace_function_value
                , ( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::release_gil_policy()
                , "Create a workspace for this integrator for the molecule group molgroup" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::ensemble
        
            typedef ::SireMove::Ensemble ( ::SireMove::OpenMMMDIntegrator::*ensemble_function_type)(  ) const;
            ensemble_function_type ensemble_function_value( &::SireMove::OpenMMMDIntegrator::ensemble );
            
            OpenMMMDIntegrator_exposer.def( 
                "ensemble"
                , ensemble_function_value
                , bp::release_gil_policy()
                , "Return the ensemble of this integrator" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::equilibrateSystem
        
            typedef ::SireSystem::System ( ::SireMove::OpenMMMDIntegrator::*equilibrateSystem_function_type)( ::SireSystem::System &,::SireUnits::Dimension::Time,int ) ;
            equilibrateSystem_function_type equilibrateSystem_function_value( &::SireMove::OpenMMMDIntegrator::equilibrateSystem );
            
            OpenMMMDIntegrator_exposer.def( 
                "equilibrateSystem"
                , equilibrateSystem_function_value
                , ( bp::arg("system"), bp::arg("equib_time_step"), bp::arg("equib_steps") )
                , bp::release_gil_policy()
                , "\n annealLambda will equilibrate the system to the current alchemical lambda\n value of the system\n Par:am system                Sire System including molegroup, forcefield\n                              positions etc\n Par:am timestep              Default = 0.005. Time step used of the\n equilibration to the desired lambda\n Par:am annealingSteps        Default = 1000. Number of steps used for the\n annealing\n Return:                      Sire system with updated coordinates and\n velocities.\n" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getAndersen
        
            typedef bool ( ::SireMove::OpenMMMDIntegrator::*getAndersen_function_type)(  ) ;
            getAndersen_function_type getAndersen_function_value( &::SireMove::OpenMMMDIntegrator::getAndersen );
            
            OpenMMMDIntegrator_exposer.def( 
                "getAndersen"
                , getAndersen_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getAndersenFrequency
        
            typedef double ( ::SireMove::OpenMMMDIntegrator::*getAndersenFrequency_function_type)(  ) ;
            getAndersenFrequency_function_type getAndersenFrequency_function_value( &::SireMove::OpenMMMDIntegrator::getAndersenFrequency );
            
            OpenMMMDIntegrator_exposer.def( 
                "getAndersenFrequency"
                , getAndersenFrequency_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getBufferFrequency
        
            typedef int ( ::SireMove::OpenMMMDIntegrator::*getBufferFrequency_function_type)(  ) ;
            getBufferFrequency_function_type getBufferFrequency_function_value( &::SireMove::OpenMMMDIntegrator::getBufferFrequency );
            
            OpenMMMDIntegrator_exposer.def( 
                "getBufferFrequency"
                , getBufferFrequency_function_value
                , bp::release_gil_policy()
                , "Get the frequency of buffering coordinates" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getCMMremovalFrequency
        
            typedef int ( ::SireMove::OpenMMMDIntegrator::*getCMMremovalFrequency_function_type)(  ) ;
            getCMMremovalFrequency_function_type getCMMremovalFrequency_function_value( &::SireMove::OpenMMMDIntegrator::getCMMremovalFrequency );
            
            OpenMMMDIntegrator_exposer.def( 
                "getCMMremovalFrequency"
                , getCMMremovalFrequency_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getConstraintType
        
            typedef ::QString ( ::SireMove::OpenMMMDIntegrator::*getConstraintType_function_type)(  ) ;
            getConstraintType_function_type getConstraintType_function_value( &::SireMove::OpenMMMDIntegrator::getConstraintType );
            
            OpenMMMDIntegrator_exposer.def( 
                "getConstraintType"
                , getConstraintType_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getCutoffDistance
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::OpenMMMDIntegrator::*getCutoffDistance_function_type)(  ) ;
            getCutoffDistance_function_type getCutoffDistance_function_value( &::SireMove::OpenMMMDIntegrator::getCutoffDistance );
            
            OpenMMMDIntegrator_exposer.def( 
                "getCutoffDistance"
                , getCutoffDistance_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getCutoffType
        
            typedef ::QString ( ::SireMove::OpenMMMDIntegrator::*getCutoffType_function_type)(  ) ;
            getCutoffType_function_type getCutoffType_function_value( &::SireMove::OpenMMMDIntegrator::getCutoffType );
            
            OpenMMMDIntegrator_exposer.def( 
                "getCutoffType"
                , getCutoffType_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getDeviceIndex
        
            typedef ::QString ( ::SireMove::OpenMMMDIntegrator::*getDeviceIndex_function_type)(  ) ;
            getDeviceIndex_function_type getDeviceIndex_function_value( &::SireMove::OpenMMMDIntegrator::getDeviceIndex );
            
            OpenMMMDIntegrator_exposer.def( 
                "getDeviceIndex"
                , getDeviceIndex_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getFieldDielectric
        
            typedef double ( ::SireMove::OpenMMMDIntegrator::*getFieldDielectric_function_type)(  ) ;
            getFieldDielectric_function_type getFieldDielectric_function_value( &::SireMove::OpenMMMDIntegrator::getFieldDielectric );
            
            OpenMMMDIntegrator_exposer.def( 
                "getFieldDielectric"
                , getFieldDielectric_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getFriction
        
            typedef ::SireUnits::Dimension::Time ( ::SireMove::OpenMMMDIntegrator::*getFriction_function_type)(  ) ;
            getFriction_function_type getFriction_function_value( &::SireMove::OpenMMMDIntegrator::getFriction );
            
            OpenMMMDIntegrator_exposer.def( 
                "getFriction"
                , getFriction_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getIntegrationTolerance
        
            typedef double ( ::SireMove::OpenMMMDIntegrator::*getIntegrationTolerance_function_type)(  ) ;
            getIntegrationTolerance_function_type getIntegrationTolerance_function_value( &::SireMove::OpenMMMDIntegrator::getIntegrationTolerance );
            
            OpenMMMDIntegrator_exposer.def( 
                "getIntegrationTolerance"
                , getIntegrationTolerance_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getIntegrator
        
            typedef ::QString ( ::SireMove::OpenMMMDIntegrator::*getIntegrator_function_type)(  ) ;
            getIntegrator_function_type getIntegrator_function_value( &::SireMove::OpenMMMDIntegrator::getIntegrator );
            
            OpenMMMDIntegrator_exposer.def( 
                "getIntegrator"
                , getIntegrator_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getKineticEnergy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::OpenMMMDIntegrator::*getKineticEnergy_function_type)(  ) ;
            getKineticEnergy_function_type getKineticEnergy_function_value( &::SireMove::OpenMMMDIntegrator::getKineticEnergy );
            
            OpenMMMDIntegrator_exposer.def( 
                "getKineticEnergy"
                , getKineticEnergy_function_value
                , bp::release_gil_policy()
                , "\n <Returns the kinetic energy of the OpenMM system>\n minimizeEnergy will find the nearest local potential energy minimum,\n given the current Sire::System. It calls the\n LocalEnergyMinimizer :: minimize() function of OpenMM.\n Par:am system                Sire System including molegroup, forcefield\n                              positions etc\n Return:                      Kinetic energy computed with OpenMM.\n" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getLJDispersion
        
            typedef bool ( ::SireMove::OpenMMMDIntegrator::*getLJDispersion_function_type)(  ) ;
            getLJDispersion_function_type getLJDispersion_function_value( &::SireMove::OpenMMMDIntegrator::getLJDispersion );
            
            OpenMMMDIntegrator_exposer.def( 
                "getLJDispersion"
                , getLJDispersion_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getMCBarostat
        
            typedef bool ( ::SireMove::OpenMMMDIntegrator::*getMCBarostat_function_type)(  ) ;
            getMCBarostat_function_type getMCBarostat_function_value( &::SireMove::OpenMMMDIntegrator::getMCBarostat );
            
            OpenMMMDIntegrator_exposer.def( 
                "getMCBarostat"
                , getMCBarostat_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getMCBarostatFrequency
        
            typedef int ( ::SireMove::OpenMMMDIntegrator::*getMCBarostatFrequency_function_type)(  ) ;
            getMCBarostatFrequency_function_type getMCBarostatFrequency_function_value( &::SireMove::OpenMMMDIntegrator::getMCBarostatFrequency );
            
            OpenMMMDIntegrator_exposer.def( 
                "getMCBarostatFrequency"
                , getMCBarostatFrequency_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getPlatform
        
            typedef ::QString ( ::SireMove::OpenMMMDIntegrator::*getPlatform_function_type)(  ) ;
            getPlatform_function_type getPlatform_function_value( &::SireMove::OpenMMMDIntegrator::getPlatform );
            
            OpenMMMDIntegrator_exposer.def( 
                "getPlatform"
                , getPlatform_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getPotentialEnergy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::OpenMMMDIntegrator::*getPotentialEnergy_function_type)( ::SireSystem::System const & ) ;
            getPotentialEnergy_function_type getPotentialEnergy_function_value( &::SireMove::OpenMMMDIntegrator::getPotentialEnergy );
            
            OpenMMMDIntegrator_exposer.def( 
                "getPotentialEnergy"
                , getPotentialEnergy_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getPrecision
        
            typedef ::QString ( ::SireMove::OpenMMMDIntegrator::*getPrecision_function_type)(  ) ;
            getPrecision_function_type getPrecision_function_value( &::SireMove::OpenMMMDIntegrator::getPrecision );
            
            OpenMMMDIntegrator_exposer.def( 
                "getPrecision"
                , getPrecision_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getPressure
        
            typedef ::SireUnits::Dimension::Pressure ( ::SireMove::OpenMMMDIntegrator::*getPressure_function_type)(  ) ;
            getPressure_function_type getPressure_function_value( &::SireMove::OpenMMMDIntegrator::getPressure );
            
            OpenMMMDIntegrator_exposer.def( 
                "getPressure"
                , getPressure_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getRestraint
        
            typedef bool ( ::SireMove::OpenMMMDIntegrator::*getRestraint_function_type)(  ) ;
            getRestraint_function_type getRestraint_function_value( &::SireMove::OpenMMMDIntegrator::getRestraint );
            
            OpenMMMDIntegrator_exposer.def( 
                "getRestraint"
                , getRestraint_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getTemperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::SireMove::OpenMMMDIntegrator::*getTemperature_function_type)(  ) ;
            getTemperature_function_type getTemperature_function_value( &::SireMove::OpenMMMDIntegrator::getTemperature );
            
            OpenMMMDIntegrator_exposer.def( 
                "getTemperature"
                , getTemperature_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getTimetoSkip
        
            typedef ::SireUnits::Dimension::Time ( ::SireMove::OpenMMMDIntegrator::*getTimetoSkip_function_type)(  ) ;
            getTimetoSkip_function_type getTimetoSkip_function_value( &::SireMove::OpenMMMDIntegrator::getTimetoSkip );
            
            OpenMMMDIntegrator_exposer.def( 
                "getTimetoSkip"
                , getTimetoSkip_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::getToleranceEwaldPME
        
            typedef double ( ::SireMove::OpenMMMDIntegrator::*getToleranceEwaldPME_function_type)(  ) ;
            getToleranceEwaldPME_function_type getToleranceEwaldPME_function_value( &::SireMove::OpenMMMDIntegrator::getToleranceEwaldPME );
            
            OpenMMMDIntegrator_exposer.def( 
                "getToleranceEwaldPME"
                , getToleranceEwaldPME_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::initialise
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*initialise_function_type)(  ) ;
            initialise_function_type initialise_function_value( &::SireMove::OpenMMMDIntegrator::initialise );
            
            OpenMMMDIntegrator_exposer.def( 
                "initialise"
                , initialise_function_value
                , bp::release_gil_policy()
                , "Integrate the coordinates of the atoms in the molecules in molgroup\nusing the forces in forcetable, using the optionally supplied\nproperty map to find the necessary molecular properties\nThrow: SireMol::missing_molecule\nThrow: SireBase::missing_property\nThrow: SireError:invalid_cast\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::integrate
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*integrate_function_type)( ::SireMove::IntegratorWorkspace &,::SireCAS::Symbol const &,::SireUnits::Dimension::Time,int,bool ) ;
            integrate_function_type integrate_function_value( &::SireMove::OpenMMMDIntegrator::integrate );
            
            OpenMMMDIntegrator_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("workspace"), bp::arg("nrg_component"), bp::arg("timestep"), bp::arg("nmoves"), bp::arg("record_stats") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::isTimeReversible
        
            typedef bool ( ::SireMove::OpenMMMDIntegrator::*isTimeReversible_function_type)(  ) const;
            isTimeReversible_function_type isTimeReversible_function_value( &::SireMove::OpenMMMDIntegrator::isTimeReversible );
            
            OpenMMMDIntegrator_exposer.def( 
                "isTimeReversible"
                , isTimeReversible_function_value
                , bp::release_gil_policy()
                , "Return whether or not this integrator is time-reversible" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::minimiseEnergy
        
            typedef ::SireSystem::System ( ::SireMove::OpenMMMDIntegrator::*minimiseEnergy_function_type)( ::SireSystem::System &,double,int ) ;
            minimiseEnergy_function_type minimiseEnergy_function_value( &::SireMove::OpenMMMDIntegrator::minimiseEnergy );
            
            OpenMMMDIntegrator_exposer.def( 
                "minimiseEnergy"
                , minimiseEnergy_function_value
                , ( bp::arg("system"), bp::arg("tolerance"), bp::arg("max_iteration") )
                , bp::release_gil_policy()
                , "\n <Runs an energy Minimisation on the current system.>\n minimizeEnergy will find the nearest local potential energy minimum,\n given the current Sire::System. It calls the\n LocalEnergyMinimizer :: minimize() function of OpenMM.\n Par:am system                Sire System including molegroup, forcefield\n                              positions etc\n Par:am tolerance             Default = 1. This specifies how precisely the\n energy minimum must be located. Minimisation will be halted once the\n root-mean-square value of all force components reaches this tolerance.\n Par:am max_iteration         Default = 1000. this specifies the number of\n iterations are run for the minimisation. If max_iteration = 0, the\n iteration will run until convergence.\n\n Return:                      Sire System, with the updated energy\n minimised coordinates.\n" );
        
        }
        OpenMMMDIntegrator_exposer.def( bp::self != bp::self );
        { //::SireMove::OpenMMMDIntegrator::operator=
        
            typedef ::SireMove::OpenMMMDIntegrator & ( ::SireMove::OpenMMMDIntegrator::*assign_function_type)( ::SireMove::OpenMMMDIntegrator const & ) ;
            assign_function_type assign_function_value( &::SireMove::OpenMMMDIntegrator::operator= );
            
            OpenMMMDIntegrator_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        OpenMMMDIntegrator_exposer.def( bp::self == bp::self );
        { //::SireMove::OpenMMMDIntegrator::setAndersen
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setAndersen_function_type)( bool ) ;
            setAndersen_function_type setAndersen_function_value( &::SireMove::OpenMMMDIntegrator::setAndersen );
            
            OpenMMMDIntegrator_exposer.def( 
                "setAndersen"
                , setAndersen_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set Andersen thermostat" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setAndersenFrequency
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setAndersenFrequency_function_type)( double ) ;
            setAndersenFrequency_function_type setAndersenFrequency_function_value( &::SireMove::OpenMMMDIntegrator::setAndersenFrequency );
            
            OpenMMMDIntegrator_exposer.def( 
                "setAndersenFrequency"
                , setAndersenFrequency_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Andersen Thermostat frequency collision" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setBufferFrequency
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setBufferFrequency_function_type)( int ) ;
            setBufferFrequency_function_type setBufferFrequency_function_value( &::SireMove::OpenMMMDIntegrator::setBufferFrequency );
            
            OpenMMMDIntegrator_exposer.def( 
                "setBufferFrequency"
                , setBufferFrequency_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Center of Mass motion removal frequency" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setCMMremovalFrequency
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setCMMremovalFrequency_function_type)( int ) ;
            setCMMremovalFrequency_function_type setCMMremovalFrequency_function_value( &::SireMove::OpenMMMDIntegrator::setCMMremovalFrequency );
            
            OpenMMMDIntegrator_exposer.def( 
                "setCMMremovalFrequency"
                , setCMMremovalFrequency_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Center of Mass motion removal frequency" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setConstraintType
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setConstraintType_function_type)( ::QString ) ;
            setConstraintType_function_type setConstraintType_function_value( &::SireMove::OpenMMMDIntegrator::setConstraintType );
            
            OpenMMMDIntegrator_exposer.def( 
                "setConstraintType"
                , setConstraintType_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Constraint type: none, hbonds, allbonds, hangles" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setCutoffDistance
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setCutoffDistance_function_type)( ::SireUnits::Dimension::Length ) ;
            setCutoffDistance_function_type setCutoffDistance_function_value( &::SireMove::OpenMMMDIntegrator::setCutoffDistance );
            
            OpenMMMDIntegrator_exposer.def( 
                "setCutoffDistance"
                , setCutoffDistance_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the cutoff distance in A" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setCutoffType
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setCutoffType_function_type)( ::QString ) ;
            setCutoffType_function_type setCutoffType_function_value( &::SireMove::OpenMMMDIntegrator::setCutoffType );
            
            OpenMMMDIntegrator_exposer.def( 
                "setCutoffType"
                , setCutoffType_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the cufott type: nocutoff, cutoffnonperiodic, cutoffperiodic" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setDeviceIndex
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setDeviceIndex_function_type)( ::QString ) ;
            setDeviceIndex_function_type setDeviceIndex_function_value( &::SireMove::OpenMMMDIntegrator::setDeviceIndex );
            
            OpenMMMDIntegrator_exposer.def( 
                "setDeviceIndex"
                , setDeviceIndex_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the OpenMM Platform: CUDA, OpenCL, CPU" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setFieldDielectric
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setFieldDielectric_function_type)( double ) ;
            setFieldDielectric_function_type setFieldDielectric_function_value( &::SireMove::OpenMMMDIntegrator::setFieldDielectric );
            
            OpenMMMDIntegrator_exposer.def( 
                "setFieldDielectric"
                , setFieldDielectric_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the dielectric constant" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setFriction
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setFriction_function_type)( ::SireUnits::Dimension::Time ) ;
            setFriction_function_type setFriction_function_value( &::SireMove::OpenMMMDIntegrator::setFriction );
            
            OpenMMMDIntegrator_exposer.def( 
                "setFriction"
                , setFriction_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setIntegrationTolerance
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setIntegrationTolerance_function_type)( double ) ;
            setIntegrationTolerance_function_type setIntegrationTolerance_function_value( &::SireMove::OpenMMMDIntegrator::setIntegrationTolerance );
            
            OpenMMMDIntegrator_exposer.def( 
                "setIntegrationTolerance"
                , setIntegrationTolerance_function_value
                , ( bp::arg("tollerance") )
                , bp::release_gil_policy()
                , "Set the integration tolerance" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setIntegrator
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setIntegrator_function_type)( ::QString ) ;
            setIntegrator_function_type setIntegrator_function_value( &::SireMove::OpenMMMDIntegrator::setIntegrator );
            
            OpenMMMDIntegrator_exposer.def( 
                "setIntegrator"
                , setIntegrator_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setLJDispersion
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setLJDispersion_function_type)( bool ) ;
            setLJDispersion_function_type setLJDispersion_function_value( &::SireMove::OpenMMMDIntegrator::setLJDispersion );
            
            OpenMMMDIntegrator_exposer.def( 
                "setLJDispersion"
                , setLJDispersion_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Retraint mode" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setMCBarostat
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setMCBarostat_function_type)( bool ) ;
            setMCBarostat_function_type setMCBarostat_function_value( &::SireMove::OpenMMMDIntegrator::setMCBarostat );
            
            OpenMMMDIntegrator_exposer.def( 
                "setMCBarostat"
                , setMCBarostat_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set Monte Carlo Barostat onoff" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setMCBarostatFrequency
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setMCBarostatFrequency_function_type)( int ) ;
            setMCBarostatFrequency_function_type setMCBarostatFrequency_function_value( &::SireMove::OpenMMMDIntegrator::setMCBarostatFrequency );
            
            OpenMMMDIntegrator_exposer.def( 
                "setMCBarostatFrequency"
                , setMCBarostatFrequency_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Monte Carlo Barostat frequency in time speps" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setPlatform
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setPlatform_function_type)( ::QString ) ;
            setPlatform_function_type setPlatform_function_value( &::SireMove::OpenMMMDIntegrator::setPlatform );
            
            OpenMMMDIntegrator_exposer.def( 
                "setPlatform"
                , setPlatform_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the OpenMM Platform: CUDA, OpenCL, CPU" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setPrecision
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setPrecision_function_type)( ::QString ) ;
            setPrecision_function_type setPrecision_function_value( &::SireMove::OpenMMMDIntegrator::setPrecision );
            
            OpenMMMDIntegrator_exposer.def( 
                "setPrecision"
                , setPrecision_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the OpenMM Precision" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setPressure
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setPressure_function_type)( ::SireUnits::Dimension::Pressure ) ;
            setPressure_function_type setPressure_function_value( &::SireMove::OpenMMMDIntegrator::setPressure );
            
            OpenMMMDIntegrator_exposer.def( 
                "setPressure"
                , setPressure_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Pressure" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setReinitialiseContext
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setReinitialiseContext_function_type)( bool ) ;
            setReinitialiseContext_function_type setReinitialiseContext_function_value( &::SireMove::OpenMMMDIntegrator::setReinitialiseContext );
            
            OpenMMMDIntegrator_exposer.def( 
                "setReinitialiseContext"
                , setReinitialiseContext_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the flag to reinitialise the context" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setRestraint
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setRestraint_function_type)( bool ) ;
            setRestraint_function_type setRestraint_function_value( &::SireMove::OpenMMMDIntegrator::setRestraint );
            
            OpenMMMDIntegrator_exposer.def( 
                "setRestraint"
                , setRestraint_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Retraint mode" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setTemperature
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setTemperature_function_type)( ::SireUnits::Dimension::Temperature ) ;
            setTemperature_function_type setTemperature_function_value( &::SireMove::OpenMMMDIntegrator::setTemperature );
            
            OpenMMMDIntegrator_exposer.def( 
                "setTemperature"
                , setTemperature_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Temperature" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setTimetoSkip
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setTimetoSkip_function_type)( ::SireUnits::Dimension::Time ) ;
            setTimetoSkip_function_type setTimetoSkip_function_value( &::SireMove::OpenMMMDIntegrator::setTimetoSkip );
            
            OpenMMMDIntegrator_exposer.def( 
                "setTimetoSkip"
                , setTimetoSkip_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Get total time to skip" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::setToleranceEwaldPME
        
            typedef void ( ::SireMove::OpenMMMDIntegrator::*setToleranceEwaldPME_function_type)( double ) ;
            setToleranceEwaldPME_function_type setToleranceEwaldPME_function_value( &::SireMove::OpenMMMDIntegrator::setToleranceEwaldPME );
            
            OpenMMMDIntegrator_exposer.def( 
                "setToleranceEwaldPME"
                , setToleranceEwaldPME_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::toString
        
            typedef ::QString ( ::SireMove::OpenMMMDIntegrator::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::OpenMMMDIntegrator::toString );
            
            OpenMMMDIntegrator_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this integrator" );
        
        }
        { //::SireMove::OpenMMMDIntegrator::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::OpenMMMDIntegrator::typeName );
            
            OpenMMMDIntegrator_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        OpenMMMDIntegrator_exposer.staticmethod( "typeName" );
        OpenMMMDIntegrator_exposer.def( "__copy__", &__copy__);
        OpenMMMDIntegrator_exposer.def( "__deepcopy__", &__copy__);
        OpenMMMDIntegrator_exposer.def( "clone", &__copy__);
        OpenMMMDIntegrator_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::OpenMMMDIntegrator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        OpenMMMDIntegrator_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::OpenMMMDIntegrator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        OpenMMMDIntegrator_exposer.def_pickle(sire_pickle_suite< ::SireMove::OpenMMMDIntegrator >());
        OpenMMMDIntegrator_exposer.def( "__str__", &__str__< ::SireMove::OpenMMMDIntegrator > );
        OpenMMMDIntegrator_exposer.def( "__repr__", &__str__< ::SireMove::OpenMMMDIntegrator > );
    }

}
