// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SQM.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/findexe.h"

#include "SireBase/sire_process.h"

#include "SireBase/tempdir.h"

#include "SireError/errors.h"

#include "SireFF/potentialtable.h"

#include "SireMM/cljprobe.h"

#include "SireMaths/vector.h"

#include "SireMol/element.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "SireVol/grid.h"

#include "latticecharges.h"

#include "qmpotential.h"

#include "sqm.h"

#include "tostring.h"

#include <QDebug>

#include <QTime>

#include "sqm.h"

Squire::SQM __copy__(const Squire::SQM &other){ return Squire::SQM(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_SQM_class(){

    { //::Squire::SQM
        typedef bp::class_< Squire::SQM, bp::bases< Squire::QMProgram, SireBase::Property > > SQM_exposer_t;
        SQM_exposer_t SQM_exposer = SQM_exposer_t( "SQM", "This is a wrapper that allows SQM to be used to calculate\nQM and QMMM energies (SQM is the semiempirical QM program\nthat comes free with AmberTools)\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope SQM_scope( SQM_exposer );
        SQM_exposer.def( bp::init< QString const & >(( bp::arg("SQM") ), "Copy constructor") );
        SQM_exposer.def( bp::init< Squire::SQM const & >(( bp::arg("other") ), "Copy constructor") );
        { //::Squire::SQM::energyTemplate
        
            typedef ::QString const & ( ::Squire::SQM::*energyTemplate_function_type)(  ) const;
            energyTemplate_function_type energyTemplate_function_value( &::Squire::SQM::energyTemplate );
            
            SQM_exposer.def( 
                "energyTemplate"
                , energyTemplate_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the template for the command file to be used to get SQM\nto calculate the energy." );
        
        }
        { //::Squire::SQM::environment
        
            typedef ::QHash< QString, QString > const & ( ::Squire::SQM::*environment_function_type)(  ) const;
            environment_function_type environment_function_value( &::Squire::SQM::environment );
            
            SQM_exposer.def( 
                "environment"
                , environment_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return all of the environmental variables that are to be set explicitly\nwhen SQM is run. This does not include any environmental variables\nthat have not been explicitly set, but do have values" );
        
        }
        { //::Squire::SQM::environment
        
            typedef ::QString ( ::Squire::SQM::*environment_function_type)( ::QString const & ) const;
            environment_function_type environment_function_value( &::Squire::SQM::environment );
            
            SQM_exposer.def( 
                "environment"
                , environment_function_value
                , ( bp::arg("variable") )
                , bp::release_gil_policy()
                , "Return the value of the explicitly set environmental variable variable.\nA null string is returned if this variable has not been set\nexplicitly (this does not mean the variable doesnt exist - merely\nthat a specific value has not been set)" );
        
        }
        { //::Squire::SQM::executable
        
            typedef ::QString ( ::Squire::SQM::*executable_function_type)(  ) const;
            executable_function_type executable_function_value( &::Squire::SQM::executable );
            
            SQM_exposer.def( 
                "executable"
                , executable_function_value
                , bp::release_gil_policy()
                , "Return the executable (full path and also arguments) to be used. This\nis null if the executable is searched for in the path" );
        
        }
        { //::Squire::SQM::expectedNumberOfQMAtoms
        
            typedef int ( ::Squire::SQM::*expectedNumberOfQMAtoms_function_type)(  ) const;
            expectedNumberOfQMAtoms_function_type expectedNumberOfQMAtoms_function_value( &::Squire::SQM::expectedNumberOfQMAtoms );
            
            SQM_exposer.def( 
                "expectedNumberOfQMAtoms"
                , expectedNumberOfQMAtoms_function_value
                , bp::release_gil_policy()
                , "Return the maximum number of expected QM atoms. This returns -1 if\nwe dont expect any QM atoms" );
        
        }
        { //::Squire::SQM::forceTemplate
        
            typedef ::QString const & ( ::Squire::SQM::*forceTemplate_function_type)(  ) const;
            forceTemplate_function_type forceTemplate_function_value( &::Squire::SQM::forceTemplate );
            
            SQM_exposer.def( 
                "forceTemplate"
                , forceTemplate_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the template for the command file to be used to get SQM\nto calculate the forces." );
        
        }
        { //::Squire::SQM::maximumNumberOfSQMInputLines
        
            typedef int ( ::Squire::SQM::*maximumNumberOfSQMInputLines_function_type)(  ) const;
            maximumNumberOfSQMInputLines_function_type maximumNumberOfSQMInputLines_function_value( &::Squire::SQM::maximumNumberOfSQMInputLines );
            
            SQM_exposer.def( 
                "maximumNumberOfSQMInputLines"
                , maximumNumberOfSQMInputLines_function_value
                , bp::release_gil_policy()
                , "Return the maximum number of supported SQM input lines. This returns\n-1 if SQM doesnt have a file size limit" );
        
        }
        { //::Squire::SQM::maximumRunTime
        
            typedef int ( ::Squire::SQM::*maximumRunTime_function_type)(  ) const;
            maximumRunTime_function_type maximumRunTime_function_value( &::Squire::SQM::maximumRunTime );
            
            SQM_exposer.def( 
                "maximumRunTime"
                , maximumRunTime_function_value
                , bp::release_gil_policy()
                , "Return the maximum runtime allowed for a SQM job, in milliseconds" );
        
        }
        { //::Squire::SQM::method
        
            typedef ::QString const & ( ::Squire::SQM::*method_function_type)(  ) const;
            method_function_type method_function_value( &::Squire::SQM::method );
            
            SQM_exposer.def( 
                "method"
                , method_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the QM method to be used to calculate the energy or\nforce (e.g. AM1, PM3, AM1d etc.). This will substitute for\n@QM_METHOD@ in the command file templates" );
        
        }
        { //::Squire::SQM::numberOfMMAtomsLimit
        
            typedef int ( ::Squire::SQM::*numberOfMMAtomsLimit_function_type)(  ) const;
            numberOfMMAtomsLimit_function_type numberOfMMAtomsLimit_function_value( &::Squire::SQM::numberOfMMAtomsLimit );
            
            SQM_exposer.def( 
                "numberOfMMAtomsLimit"
                , numberOfMMAtomsLimit_function_value
                , bp::release_gil_policy()
                , "Return the maximum number of MM atoms supported by SQM. This returns\n-1 if there is no limit on the number of atoms" );
        
        }
        { //::Squire::SQM::numberOfMMAtomsLimit
        
            typedef int ( ::Squire::SQM::*numberOfMMAtomsLimit_function_type)( int ) const;
            numberOfMMAtomsLimit_function_type numberOfMMAtomsLimit_function_value( &::Squire::SQM::numberOfMMAtomsLimit );
            
            SQM_exposer.def( 
                "numberOfMMAtomsLimit"
                , numberOfMMAtomsLimit_function_value
                , ( bp::arg("num_qm_atoms") )
                , bp::release_gil_policy()
                , "Return the maximum number of MM atoms supported by SQM if there\nare num_qm_atoms QM atoms. This returns\n-1 if there is no limit on the number of atoms" );
        
        }
        SQM_exposer.def( bp::self != bp::self );
        { //::Squire::SQM::operator=
        
            typedef ::Squire::SQM & ( ::Squire::SQM::*assign_function_type)( ::Squire::SQM const & ) ;
            assign_function_type assign_function_value( &::Squire::SQM::operator= );
            
            SQM_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        SQM_exposer.def( bp::self == bp::self );
        { //::Squire::SQM::setEnergyTemplate
        
            typedef void ( ::Squire::SQM::*setEnergyTemplate_function_type)( ::QString const & ) ;
            setEnergyTemplate_function_type setEnergyTemplate_function_value( &::Squire::SQM::setEnergyTemplate );
            
            SQM_exposer.def( 
                "setEnergyTemplate"
                , setEnergyTemplate_function_value
                , ( bp::arg("energy_template") )
                , bp::release_gil_policy()
                , "Set the template for the command file to be used to get\nSQM to calculate an energy. The following tags will\n" );
        
        }
        { //::Squire::SQM::setEnvironment
        
            typedef void ( ::Squire::SQM::*setEnvironment_function_type)( ::QString const &,::QString const & ) ;
            setEnvironment_function_type setEnvironment_function_value( &::Squire::SQM::setEnvironment );
            
            SQM_exposer.def( 
                "setEnvironment"
                , setEnvironment_function_value
                , ( bp::arg("variable"), bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the environmental variable variable to have the value value\nwhen the SQM executable is run. This replaces any existing\nvalue of this environmental variable" );
        
        }
        { //::Squire::SQM::setExecutable
        
            typedef void ( ::Squire::SQM::*setExecutable_function_type)( ::QString const & ) ;
            setExecutable_function_type setExecutable_function_value( &::Squire::SQM::setExecutable );
            
            SQM_exposer.def( 
                "setExecutable"
                , setExecutable_function_value
                , ( bp::arg("SQM_exe") )
                , bp::release_gil_policy()
                , "Set the SQM executable (full path and also arguments) to be used" );
        
        }
        { //::Squire::SQM::setExpectedNumberOfQMAtoms
        
            typedef void ( ::Squire::SQM::*setExpectedNumberOfQMAtoms_function_type)( int ) ;
            setExpectedNumberOfQMAtoms_function_type setExpectedNumberOfQMAtoms_function_value( &::Squire::SQM::setExpectedNumberOfQMAtoms );
            
            SQM_exposer.def( 
                "setExpectedNumberOfQMAtoms"
                , setExpectedNumberOfQMAtoms_function_value
                , ( bp::arg("natoms") )
                , bp::release_gil_policy()
                , "Set the maximum number of expected QM atoms. This is used, together with\nthe maximum number of lines in a SQM input file, to work out the maximum\nnumber of supported MM atoms" );
        
        }
        { //::Squire::SQM::setForceTemplate
        
            typedef void ( ::Squire::SQM::*setForceTemplate_function_type)( ::QString const & ) ;
            setForceTemplate_function_type setForceTemplate_function_value( &::Squire::SQM::setForceTemplate );
            
            SQM_exposer.def( 
                "setForceTemplate"
                , setForceTemplate_function_value
                , ( bp::arg("force_template") )
                , bp::release_gil_policy()
                , "Set the template for the command file to be used to get\nSQM to calculate the forces. The following tags will\n" );
        
        }
        { //::Squire::SQM::setMaximumNumberOfSQMInputLines
        
            typedef void ( ::Squire::SQM::*setMaximumNumberOfSQMInputLines_function_type)( int ) ;
            setMaximumNumberOfSQMInputLines_function_type setMaximumNumberOfSQMInputLines_function_value( &::Squire::SQM::setMaximumNumberOfSQMInputLines );
            
            SQM_exposer.def( 
                "setMaximumNumberOfSQMInputLines"
                , setMaximumNumberOfSQMInputLines_function_value
                , ( bp::arg("numlines") )
                , bp::release_gil_policy()
                , "Set the maximum number of lines that can be parsed from an SQM input file.\nCurrently, SQM has a hard-coded limit of 1000 lines" );
        
        }
        { //::Squire::SQM::setMaximumRunTime
        
            typedef void ( ::Squire::SQM::*setMaximumRunTime_function_type)( int ) ;
            setMaximumRunTime_function_type setMaximumRunTime_function_value( &::Squire::SQM::setMaximumRunTime );
            
            SQM_exposer.def( 
                "setMaximumRunTime"
                , setMaximumRunTime_function_value
                , ( bp::arg("max_runtime") )
                , bp::release_gil_policy()
                , "Set the maximum allowed runtime for the SQM job - this is used\nto detect hangs - if the SQM job takes longer than this\ntime then it is killed and an exception raised. The maximum\nruntime is measured in milliseconds" );
        
        }
        { //::Squire::SQM::setMethod
        
            typedef void ( ::Squire::SQM::*setMethod_function_type)( ::QString const & ) ;
            setMethod_function_type setMethod_function_value( &::Squire::SQM::setMethod );
            
            SQM_exposer.def( 
                "setMethod"
                , setMethod_function_value
                , ( bp::arg("method") )
                , bp::release_gil_policy()
                , "Set the QM method to be used to calculate the energy or\nforce (e.g. AM1, PM3, AM1d etc. See the AmberTools documentation\nfor SQM to find the supported methods and the string used to\nspecify that method). This will substitute for\n@QM_METHOD@ in the command file templates, and should be the same\nstring used in SQM as specified in the SQM documentation" );
        
        }
        { //::Squire::SQM::setTotalCharge
        
            typedef void ( ::Squire::SQM::*setTotalCharge_function_type)( int ) ;
            setTotalCharge_function_type setTotalCharge_function_value( &::Squire::SQM::setTotalCharge );
            
            SQM_exposer.def( 
                "setTotalCharge"
                , setTotalCharge_function_value
                , ( bp::arg("charge") )
                , bp::release_gil_policy()
                , "Set the total charge of the system (in unit charges)" );
        
        }
        { //::Squire::SQM::supportsLatticeCharges
        
            typedef bool ( ::Squire::SQM::*supportsLatticeCharges_function_type)(  ) const;
            supportsLatticeCharges_function_type supportsLatticeCharges_function_value( &::Squire::SQM::supportsLatticeCharges );
            
            SQM_exposer.def( 
                "supportsLatticeCharges"
                , supportsLatticeCharges_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::Squire::SQM::toString
        
            typedef ::QString ( ::Squire::SQM::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::Squire::SQM::toString );
            
            SQM_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::Squire::SQM::totalCharge
        
            typedef int ( ::Squire::SQM::*totalCharge_function_type)(  ) const;
            totalCharge_function_type totalCharge_function_value( &::Squire::SQM::totalCharge );
            
            SQM_exposer.def( 
                "totalCharge"
                , totalCharge_function_value
                , bp::release_gil_policy()
                , "Return the total charge of the system" );
        
        }
        { //::Squire::SQM::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::SQM::typeName );
            
            SQM_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SQM_exposer.staticmethod( "typeName" );
        SQM_exposer.def( "__copy__", &__copy__);
        SQM_exposer.def( "__deepcopy__", &__copy__);
        SQM_exposer.def( "clone", &__copy__);
        SQM_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::SQM >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SQM_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::SQM >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SQM_exposer.def_pickle(sire_pickle_suite< ::Squire::SQM >());
        SQM_exposer.def( "__str__", &__str__< ::Squire::SQM > );
        SQM_exposer.def( "__repr__", &__str__< ::Squire::SQM > );
    }

}
