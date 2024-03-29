// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "BennettsRatios.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/maths.h"

#include "SireStream/registeralternativename.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "bennetts.h"

#include "tostring.h"

#include "bennetts.h"

SireAnalysis::BennettsRatios __copy__(const SireAnalysis::BennettsRatios &other){ return SireAnalysis::BennettsRatios(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_BennettsRatios_class(){

    { //::SireAnalysis::BennettsRatios
        typedef bp::class_< SireAnalysis::BennettsRatios, bp::bases< SireBase::Property > > BennettsRatios_exposer_t;
        BennettsRatios_exposer_t BennettsRatios_exposer = BennettsRatios_exposer_t( "BennettsRatios", "This class is used to hold a set of Bennets acceptance ratios\nfor a single iteration\n\nAuthor: Christopher Woods\n", bp::init< >("Construct an empty set of deltas") );
        bp::scope BennettsRatios_scope( BennettsRatios_exposer );
        BennettsRatios_exposer.def( bp::init< QList< double > const &, QMap< double, SireMaths::BennettsFreeEnergyAverage > const &, QMap< double, SireMaths::BennettsFreeEnergyAverage > const & >(( bp::arg("windows"), bp::arg("forwards_ratios"), bp::arg("backwards_ratios") ), "Construct the ratios as the ratios between each window and the windows above\n(forwards_ratios) and windows below (backwards_ratios)") );
        BennettsRatios_exposer.def( bp::init< SireAnalysis::BennettsRatios const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireAnalysis::BennettsRatios::backwardsData
        
            typedef ::QMap< double, SireMaths::BennettsFreeEnergyAverage > ( ::SireAnalysis::BennettsRatios::*backwardsData_function_type)(  ) const;
            backwardsData_function_type backwardsData_function_value( &::SireAnalysis::BennettsRatios::backwardsData );
            
            BennettsRatios_exposer.def( 
                "backwardsData"
                , backwardsData_function_value
                , bp::release_gil_policy()
                , "Return the raw data for the backwards ratios" );
        
        }
        { //::SireAnalysis::BennettsRatios::backwardsRatios
        
            typedef ::QMap< double, SireMaths::BennettsFreeEnergyAverage > ( ::SireAnalysis::BennettsRatios::*backwardsRatios_function_type)(  ) const;
            backwardsRatios_function_type backwardsRatios_function_value( &::SireAnalysis::BennettsRatios::backwardsRatios );
            
            BennettsRatios_exposer.def( 
                "backwardsRatios"
                , backwardsRatios_function_value
                , bp::release_gil_policy()
                , "Return the raw data for the backwards ratios" );
        
        }
        { //::SireAnalysis::BennettsRatios::constants
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::BennettsRatios::*constants_function_type)(  ) const;
            constants_function_type constants_function_value( &::SireAnalysis::BennettsRatios::constants );
            
            BennettsRatios_exposer.def( 
                "constants"
                , constants_function_value
                , bp::release_gil_policy()
                , "Return the constants for each set of Bennetts acceptance ratios. This\nreturns the lambda value of the from window, together with the constant\nused for the numerator from this window to the next window (which must\nbe the same as the constant used for the denominator for the next window\nback to this window)" );
        
        }
        { //::SireAnalysis::BennettsRatios::denominators
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::BennettsRatios::*denominators_function_type)(  ) const;
            denominators_function_type denominators_function_value( &::SireAnalysis::BennettsRatios::denominators );
            
            BennettsRatios_exposer.def( 
                "denominators"
                , denominators_function_value
                , bp::release_gil_policy()
                , "Return the denominators for the Bennetts acceptance ratio. This returns the\nlambda value of the previous window, together with the Bennetts ratio for the\nenergy difference from this window to the previous window" );
        
        }
        { //::SireAnalysis::BennettsRatios::forwardsData
        
            typedef ::QMap< double, SireMaths::BennettsFreeEnergyAverage > ( ::SireAnalysis::BennettsRatios::*forwardsData_function_type)(  ) const;
            forwardsData_function_type forwardsData_function_value( &::SireAnalysis::BennettsRatios::forwardsData );
            
            BennettsRatios_exposer.def( 
                "forwardsData"
                , forwardsData_function_value
                , bp::release_gil_policy()
                , "Return the raw data for the fowards ratios" );
        
        }
        { //::SireAnalysis::BennettsRatios::forwardsRatios
        
            typedef ::QMap< double, SireMaths::BennettsFreeEnergyAverage > ( ::SireAnalysis::BennettsRatios::*forwardsRatios_function_type)(  ) const;
            forwardsRatios_function_type forwardsRatios_function_value( &::SireAnalysis::BennettsRatios::forwardsRatios );
            
            BennettsRatios_exposer.def( 
                "forwardsRatios"
                , forwardsRatios_function_value
                , bp::release_gil_policy()
                , "Return the raw data for the fowards ratios" );
        
        }
        { //::SireAnalysis::BennettsRatios::integrate
        
            typedef ::SireAnalysis::PMF ( ::SireAnalysis::BennettsRatios::*integrate_function_type)(  ) const;
            integrate_function_type integrate_function_value( &::SireAnalysis::BennettsRatios::integrate );
            
            BennettsRatios_exposer.def( 
                "integrate"
                , integrate_function_value
                , bp::release_gil_policy()
                , "Integrate (sum) the deltas across the windows to return the PMF" );
        
        }
        { //::SireAnalysis::BennettsRatios::isEmpty
        
            typedef bool ( ::SireAnalysis::BennettsRatios::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireAnalysis::BennettsRatios::isEmpty );
            
            BennettsRatios_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireAnalysis::BennettsRatios::lambdaValues
        
            typedef ::QList< double > ( ::SireAnalysis::BennettsRatios::*lambdaValues_function_type)(  ) const;
            lambdaValues_function_type lambdaValues_function_value( &::SireAnalysis::BennettsRatios::lambdaValues );
            
            BennettsRatios_exposer.def( 
                "lambdaValues"
                , lambdaValues_function_value
                , bp::release_gil_policy()
                , "Return the lambda values for all of the windows" );
        
        }
        { //::SireAnalysis::BennettsRatios::merge
        
            typedef ::SireAnalysis::BennettsRatios ( *merge_function_type )( ::QList< SireAnalysis::BennettsRatios > const & );
            merge_function_type merge_function_value( &::SireAnalysis::BennettsRatios::merge );
            
            BennettsRatios_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("deltas") )
                , bp::release_gil_policy()
                , "Merge together all of the passed BennettsRatios into a single object" );
        
        }
        { //::SireAnalysis::BennettsRatios::nLambdaValues
        
            typedef int ( ::SireAnalysis::BennettsRatios::*nLambdaValues_function_type)(  ) const;
            nLambdaValues_function_type nLambdaValues_function_value( &::SireAnalysis::BennettsRatios::nLambdaValues );
            
            BennettsRatios_exposer.def( 
                "nLambdaValues"
                , nLambdaValues_function_value
                , bp::release_gil_policy()
                , "Return the number of lambda values (windows)" );
        
        }
        { //::SireAnalysis::BennettsRatios::nSamples
        
            typedef ::qint64 ( ::SireAnalysis::BennettsRatios::*nSamples_function_type)(  ) const;
            nSamples_function_type nSamples_function_value( &::SireAnalysis::BennettsRatios::nSamples );
            
            BennettsRatios_exposer.def( 
                "nSamples"
                , nSamples_function_value
                , bp::release_gil_policy()
                , "Return the total number of samples in the deltas" );
        
        }
        { //::SireAnalysis::BennettsRatios::nWindows
        
            typedef int ( ::SireAnalysis::BennettsRatios::*nWindows_function_type)(  ) const;
            nWindows_function_type nWindows_function_value( &::SireAnalysis::BennettsRatios::nWindows );
            
            BennettsRatios_exposer.def( 
                "nWindows"
                , nWindows_function_value
                , bp::release_gil_policy()
                , "Return the number of windows" );
        
        }
        { //::SireAnalysis::BennettsRatios::numerators
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::BennettsRatios::*numerators_function_type)(  ) const;
            numerators_function_type numerators_function_value( &::SireAnalysis::BennettsRatios::numerators );
            
            BennettsRatios_exposer.def( 
                "numerators"
                , numerators_function_value
                , bp::release_gil_policy()
                , "Return the numerators for the Bennetts acceptance ratio. This returns the\nlambda value of the from window, together with the Bennetts ratio for the\nenergy difference from this window to the next window" );
        
        }
        BennettsRatios_exposer.def( bp::self != bp::self );
        BennettsRatios_exposer.def( bp::self + bp::self );
        { //::SireAnalysis::BennettsRatios::operator=
        
            typedef ::SireAnalysis::BennettsRatios & ( ::SireAnalysis::BennettsRatios::*assign_function_type)( ::SireAnalysis::BennettsRatios const & ) ;
            assign_function_type assign_function_value( &::SireAnalysis::BennettsRatios::operator= );
            
            BennettsRatios_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        BennettsRatios_exposer.def( bp::self == bp::self );
        { //::SireAnalysis::BennettsRatios::sum
        
            typedef ::SireAnalysis::PMF ( ::SireAnalysis::BennettsRatios::*sum_function_type)(  ) const;
            sum_function_type sum_function_value( &::SireAnalysis::BennettsRatios::sum );
            
            BennettsRatios_exposer.def( 
                "sum"
                , sum_function_value
                , bp::release_gil_policy()
                , "Integrate (sum) the deltas across the windows to return the PMF" );
        
        }
        { //::SireAnalysis::BennettsRatios::temperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::SireAnalysis::BennettsRatios::*temperature_function_type)(  ) const;
            temperature_function_type temperature_function_value( &::SireAnalysis::BennettsRatios::temperature );
            
            BennettsRatios_exposer.def( 
                "temperature"
                , temperature_function_value
                , bp::release_gil_policy()
                , "Return the temperature at which the Bennetts deltas were all collected" );
        
        }
        { //::SireAnalysis::BennettsRatios::toString
        
            typedef ::QString ( ::SireAnalysis::BennettsRatios::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireAnalysis::BennettsRatios::toString );
            
            BennettsRatios_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireAnalysis::BennettsRatios::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireAnalysis::BennettsRatios::typeName );
            
            BennettsRatios_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireAnalysis::BennettsRatios::values
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::BennettsRatios::*values_function_type)(  ) const;
            values_function_type values_function_value( &::SireAnalysis::BennettsRatios::values );
            
            BennettsRatios_exposer.def( 
                "values"
                , values_function_value
                , bp::release_gil_policy()
                , "Return the values between windows. This returns the value of lambda of the from\nwindow, and the difference in free energy between this and the next window" );
        
        }
        { //::SireAnalysis::BennettsRatios::what
        
            typedef char const * ( ::SireAnalysis::BennettsRatios::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireAnalysis::BennettsRatios::what );
            
            BennettsRatios_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireAnalysis::BennettsRatios::windows
        
            typedef ::QList< double > ( ::SireAnalysis::BennettsRatios::*windows_function_type)(  ) const;
            windows_function_type windows_function_value( &::SireAnalysis::BennettsRatios::windows );
            
            BennettsRatios_exposer.def( 
                "windows"
                , windows_function_value
                , bp::release_gil_policy()
                , "Return the values of all of the windows" );
        
        }
        BennettsRatios_exposer.staticmethod( "merge" );
        BennettsRatios_exposer.staticmethod( "typeName" );
        BennettsRatios_exposer.def( "__copy__", &__copy__);
        BennettsRatios_exposer.def( "__deepcopy__", &__copy__);
        BennettsRatios_exposer.def( "clone", &__copy__);
        BennettsRatios_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireAnalysis::BennettsRatios >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BennettsRatios_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireAnalysis::BennettsRatios >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BennettsRatios_exposer.def_pickle(sire_pickle_suite< ::SireAnalysis::BennettsRatios >());
        BennettsRatios_exposer.def( "__str__", &__str__< ::SireAnalysis::BennettsRatios > );
        BennettsRatios_exposer.def( "__repr__", &__str__< ::SireAnalysis::BennettsRatios > );
    }

}
