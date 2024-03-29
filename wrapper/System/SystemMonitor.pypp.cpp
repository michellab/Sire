// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SystemMonitor.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "system.h"

#include "systemmonitor.h"

#include <QMutex>

#include "systemmonitor.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_SystemMonitor_class(){

    { //::SireSystem::SystemMonitor
        typedef bp::class_< SireSystem::SystemMonitor, bp::bases< SireBase::Property >, boost::noncopyable > SystemMonitor_exposer_t;
        SystemMonitor_exposer_t SystemMonitor_exposer = SystemMonitor_exposer_t( "SystemMonitor", "This is the virtual base class of all system monitors. A system\nmonitor is an object that monitors a system during a simulation,\ne.g. collecting the average energy, saving a radial distribution\nfunction etc.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope SystemMonitor_scope( SystemMonitor_exposer );
        { //::SireSystem::SystemMonitor::clearStatistics
        
            typedef void ( ::SireSystem::SystemMonitor::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireSystem::SystemMonitor::clearStatistics );
            
            SystemMonitor_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireSystem::SystemMonitor::monitor
        
            typedef void ( ::SireSystem::SystemMonitor::*monitor_function_type)( ::SireSystem::System & ) ;
            monitor_function_type monitor_function_value( &::SireSystem::SystemMonitor::monitor );
            
            SystemMonitor_exposer.def( 
                "monitor"
                , monitor_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireSystem::SystemMonitor::null
        
            typedef ::SireSystem::NullMonitor const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireSystem::SystemMonitor::null );
            
            SystemMonitor_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireSystem::SystemMonitor::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::SystemMonitor::typeName );
            
            SystemMonitor_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SystemMonitor_exposer.staticmethod( "null" );
        SystemMonitor_exposer.staticmethod( "typeName" );
        SystemMonitor_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::SystemMonitor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SystemMonitor_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::SystemMonitor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SystemMonitor_exposer.def_pickle(sire_pickle_suite< ::SireSystem::SystemMonitor >());
        SystemMonitor_exposer.def( "__str__", &__str__< ::SireSystem::SystemMonitor > );
        SystemMonitor_exposer.def( "__repr__", &__str__< ::SireSystem::SystemMonitor > );
    }

}
