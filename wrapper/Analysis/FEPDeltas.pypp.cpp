// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "FEPDeltas.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/maths.h"

#include "SireStream/registeralternativename.h"

#include "SireStream/shareddatastream.h"

#include "fep.h"

#include "tostring.h"

#include "fep.h"

SireAnalysis::FEPDeltas __copy__(const SireAnalysis::FEPDeltas &other){ return SireAnalysis::FEPDeltas(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_FEPDeltas_class(){

    { //::SireAnalysis::FEPDeltas
        typedef bp::class_< SireAnalysis::FEPDeltas, bp::bases< SireBase::Property > > FEPDeltas_exposer_t;
        FEPDeltas_exposer_t FEPDeltas_exposer = FEPDeltas_exposer_t( "FEPDeltas", "This class is used to hold the set of FEP deltas from a single\niteration\n\nAuthor: Christopher Woods\n", bp::init< >("Construct an empty set of deltas") );
        bp::scope FEPDeltas_scope( FEPDeltas_exposer );
        FEPDeltas_exposer.def( bp::init< QList< double > const &, QMap< double, SireMaths::FreeEnergyAverage > const & >(( bp::arg("windows"), bp::arg("deltas") ), "Construct the deltas as the deltas between each window and the window above") );
        FEPDeltas_exposer.def( bp::init< QList< double > const &, QMap< double, SireMaths::FreeEnergyAverage > const &, QMap< double, SireMaths::FreeEnergyAverage > const & >(( bp::arg("windows"), bp::arg("forwards_deltas"), bp::arg("backwards_deltas") ), "Construct the deltas as the deltas between each window and the windows above\n(forwards_deltas) and windows below (backwards deltas)") );
        FEPDeltas_exposer.def( bp::init< SireAnalysis::FEPDeltas const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireAnalysis::FEPDeltas::backwardsData
        
            typedef ::QMap< double, SireMaths::FreeEnergyAverage > ( ::SireAnalysis::FEPDeltas::*backwardsData_function_type)(  ) const;
            backwardsData_function_type backwardsData_function_value( &::SireAnalysis::FEPDeltas::backwardsData );
            
            FEPDeltas_exposer.def( 
                "backwardsData"
                , backwardsData_function_value
                , bp::release_gil_policy()
                , "Return the raw data for the backwards deltas" );
        
        }
        { //::SireAnalysis::FEPDeltas::backwardsDeltas
        
            typedef ::QMap< double, SireMaths::FreeEnergyAverage > ( ::SireAnalysis::FEPDeltas::*backwardsDeltas_function_type)(  ) const;
            backwardsDeltas_function_type backwardsDeltas_function_value( &::SireAnalysis::FEPDeltas::backwardsDeltas );
            
            FEPDeltas_exposer.def( 
                "backwardsDeltas"
                , backwardsDeltas_function_value
                , bp::release_gil_policy()
                , "Return the raw data for the backwards deltas" );
        
        }
        { //::SireAnalysis::FEPDeltas::backwardsValues
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::FEPDeltas::*backwardsValues_function_type)(  ) const;
            backwardsValues_function_type backwardsValues_function_value( &::SireAnalysis::FEPDeltas::backwardsValues );
            
            FEPDeltas_exposer.def( 
                "backwardsValues"
                , backwardsValues_function_value
                , bp::release_gil_policy()
                , "Return the backwards deltas. This returns the lambda value of the from window,\ntogether with the free energy delta" );
        
        }
        { //::SireAnalysis::FEPDeltas::forwardsData
        
            typedef ::QMap< double, SireMaths::FreeEnergyAverage > ( ::SireAnalysis::FEPDeltas::*forwardsData_function_type)(  ) const;
            forwardsData_function_type forwardsData_function_value( &::SireAnalysis::FEPDeltas::forwardsData );
            
            FEPDeltas_exposer.def( 
                "forwardsData"
                , forwardsData_function_value
                , bp::release_gil_policy()
                , "Return the raw data for the fowards deltas" );
        
        }
        { //::SireAnalysis::FEPDeltas::forwardsDeltas
        
            typedef ::QMap< double, SireMaths::FreeEnergyAverage > ( ::SireAnalysis::FEPDeltas::*forwardsDeltas_function_type)(  ) const;
            forwardsDeltas_function_type forwardsDeltas_function_value( &::SireAnalysis::FEPDeltas::forwardsDeltas );
            
            FEPDeltas_exposer.def( 
                "forwardsDeltas"
                , forwardsDeltas_function_value
                , bp::release_gil_policy()
                , "Return the raw data for the fowards deltas" );
        
        }
        { //::SireAnalysis::FEPDeltas::forwardsValues
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::FEPDeltas::*forwardsValues_function_type)(  ) const;
            forwardsValues_function_type forwardsValues_function_value( &::SireAnalysis::FEPDeltas::forwardsValues );
            
            FEPDeltas_exposer.def( 
                "forwardsValues"
                , forwardsValues_function_value
                , bp::release_gil_policy()
                , "Return the forwards deltas. This returns the lambda value of the from window,\ntogether with the value of the free energy delta" );
        
        }
        { //::SireAnalysis::FEPDeltas::integrate
        
            typedef ::SireAnalysis::PMF ( ::SireAnalysis::FEPDeltas::*integrate_function_type)(  ) const;
            integrate_function_type integrate_function_value( &::SireAnalysis::FEPDeltas::integrate );
            
            FEPDeltas_exposer.def( 
                "integrate"
                , integrate_function_value
                , bp::release_gil_policy()
                , "Integrate (sum) the deltas across the windows to return the PMF" );
        
        }
        { //::SireAnalysis::FEPDeltas::isEmpty
        
            typedef bool ( ::SireAnalysis::FEPDeltas::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireAnalysis::FEPDeltas::isEmpty );
            
            FEPDeltas_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireAnalysis::FEPDeltas::lambdaValues
        
            typedef ::QList< double > ( ::SireAnalysis::FEPDeltas::*lambdaValues_function_type)(  ) const;
            lambdaValues_function_type lambdaValues_function_value( &::SireAnalysis::FEPDeltas::lambdaValues );
            
            FEPDeltas_exposer.def( 
                "lambdaValues"
                , lambdaValues_function_value
                , bp::release_gil_policy()
                , "Return the lambda values for all of the windows" );
        
        }
        { //::SireAnalysis::FEPDeltas::merge
        
            typedef ::SireAnalysis::FEPDeltas ( *merge_function_type )( ::QList< SireAnalysis::FEPDeltas > const & );
            merge_function_type merge_function_value( &::SireAnalysis::FEPDeltas::merge );
            
            FEPDeltas_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("deltas") )
                , bp::release_gil_policy()
                , "Merge together all of the passed FEPDeltas into a single object" );
        
        }
        { //::SireAnalysis::FEPDeltas::nLambdaValues
        
            typedef int ( ::SireAnalysis::FEPDeltas::*nLambdaValues_function_type)(  ) const;
            nLambdaValues_function_type nLambdaValues_function_value( &::SireAnalysis::FEPDeltas::nLambdaValues );
            
            FEPDeltas_exposer.def( 
                "nLambdaValues"
                , nLambdaValues_function_value
                , bp::release_gil_policy()
                , "Return the number of lambda values (windows)" );
        
        }
        { //::SireAnalysis::FEPDeltas::nSamples
        
            typedef ::qint64 ( ::SireAnalysis::FEPDeltas::*nSamples_function_type)(  ) const;
            nSamples_function_type nSamples_function_value( &::SireAnalysis::FEPDeltas::nSamples );
            
            FEPDeltas_exposer.def( 
                "nSamples"
                , nSamples_function_value
                , bp::release_gil_policy()
                , "Return the total number of samples in the deltas" );
        
        }
        { //::SireAnalysis::FEPDeltas::nWindows
        
            typedef int ( ::SireAnalysis::FEPDeltas::*nWindows_function_type)(  ) const;
            nWindows_function_type nWindows_function_value( &::SireAnalysis::FEPDeltas::nWindows );
            
            FEPDeltas_exposer.def( 
                "nWindows"
                , nWindows_function_value
                , bp::release_gil_policy()
                , "Return the number of windows" );
        
        }
        FEPDeltas_exposer.def( bp::self != bp::self );
        FEPDeltas_exposer.def( bp::self + bp::self );
        { //::SireAnalysis::FEPDeltas::operator=
        
            typedef ::SireAnalysis::FEPDeltas & ( ::SireAnalysis::FEPDeltas::*assign_function_type)( ::SireAnalysis::FEPDeltas const & ) ;
            assign_function_type assign_function_value( &::SireAnalysis::FEPDeltas::operator= );
            
            FEPDeltas_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        FEPDeltas_exposer.def( bp::self == bp::self );
        { //::SireAnalysis::FEPDeltas::sum
        
            typedef ::SireAnalysis::PMF ( ::SireAnalysis::FEPDeltas::*sum_function_type)(  ) const;
            sum_function_type sum_function_value( &::SireAnalysis::FEPDeltas::sum );
            
            FEPDeltas_exposer.def( 
                "sum"
                , sum_function_value
                , bp::release_gil_policy()
                , "Integrate (sum) the deltas across the windows to return the PMF" );
        
        }
        { //::SireAnalysis::FEPDeltas::sumBackwards
        
            typedef ::SireAnalysis::PMF ( ::SireAnalysis::FEPDeltas::*sumBackwards_function_type)(  ) const;
            sumBackwards_function_type sumBackwards_function_value( &::SireAnalysis::FEPDeltas::sumBackwards );
            
            FEPDeltas_exposer.def( 
                "sumBackwards"
                , sumBackwards_function_value
                , bp::release_gil_policy()
                , "Integrate (sum) the backwards deltas across the windows to return the PMF" );
        
        }
        { //::SireAnalysis::FEPDeltas::sumBackwardsTaylor
        
            typedef ::SireAnalysis::PMF ( ::SireAnalysis::FEPDeltas::*sumBackwardsTaylor_function_type)(  ) const;
            sumBackwardsTaylor_function_type sumBackwardsTaylor_function_value( &::SireAnalysis::FEPDeltas::sumBackwardsTaylor );
            
            FEPDeltas_exposer.def( 
                "sumBackwardsTaylor"
                , sumBackwardsTaylor_function_value
                , bp::release_gil_policy()
                , "Integrate (sum) the backwards Taylor expansions across the windows to return the PMF" );
        
        }
        { //::SireAnalysis::FEPDeltas::sumForwards
        
            typedef ::SireAnalysis::PMF ( ::SireAnalysis::FEPDeltas::*sumForwards_function_type)(  ) const;
            sumForwards_function_type sumForwards_function_value( &::SireAnalysis::FEPDeltas::sumForwards );
            
            FEPDeltas_exposer.def( 
                "sumForwards"
                , sumForwards_function_value
                , bp::release_gil_policy()
                , "Integrate (sum) the forwards deltas across the windows to return the PMF" );
        
        }
        { //::SireAnalysis::FEPDeltas::sumForwardsTaylor
        
            typedef ::SireAnalysis::PMF ( ::SireAnalysis::FEPDeltas::*sumForwardsTaylor_function_type)(  ) const;
            sumForwardsTaylor_function_type sumForwardsTaylor_function_value( &::SireAnalysis::FEPDeltas::sumForwardsTaylor );
            
            FEPDeltas_exposer.def( 
                "sumForwardsTaylor"
                , sumForwardsTaylor_function_value
                , bp::release_gil_policy()
                , "Integrate (sum) the forwards taylor expansions across the windows to return the PMF" );
        
        }
        { //::SireAnalysis::FEPDeltas::temperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::SireAnalysis::FEPDeltas::*temperature_function_type)(  ) const;
            temperature_function_type temperature_function_value( &::SireAnalysis::FEPDeltas::temperature );
            
            FEPDeltas_exposer.def( 
                "temperature"
                , temperature_function_value
                , bp::release_gil_policy()
                , "Return the temperature at which the FEP deltas were all collected" );
        
        }
        { //::SireAnalysis::FEPDeltas::toString
        
            typedef ::QString ( ::SireAnalysis::FEPDeltas::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireAnalysis::FEPDeltas::toString );
            
            FEPDeltas_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireAnalysis::FEPDeltas::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireAnalysis::FEPDeltas::typeName );
            
            FEPDeltas_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireAnalysis::FEPDeltas::values
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::FEPDeltas::*values_function_type)(  ) const;
            values_function_type values_function_value( &::SireAnalysis::FEPDeltas::values );
            
            FEPDeltas_exposer.def( 
                "values"
                , values_function_value
                , bp::release_gil_policy()
                , "Return the values between windows. This returns the average of the\nforwards and backwards values" );
        
        }
        { //::SireAnalysis::FEPDeltas::what
        
            typedef char const * ( ::SireAnalysis::FEPDeltas::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireAnalysis::FEPDeltas::what );
            
            FEPDeltas_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireAnalysis::FEPDeltas::windows
        
            typedef ::QList< double > ( ::SireAnalysis::FEPDeltas::*windows_function_type)(  ) const;
            windows_function_type windows_function_value( &::SireAnalysis::FEPDeltas::windows );
            
            FEPDeltas_exposer.def( 
                "windows"
                , windows_function_value
                , bp::release_gil_policy()
                , "Return the values of all of the windows" );
        
        }
        FEPDeltas_exposer.staticmethod( "merge" );
        FEPDeltas_exposer.staticmethod( "typeName" );
        FEPDeltas_exposer.def( "__copy__", &__copy__);
        FEPDeltas_exposer.def( "__deepcopy__", &__copy__);
        FEPDeltas_exposer.def( "clone", &__copy__);
        FEPDeltas_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireAnalysis::FEPDeltas >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FEPDeltas_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireAnalysis::FEPDeltas >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FEPDeltas_exposer.def_pickle(sire_pickle_suite< ::SireAnalysis::FEPDeltas >());
        FEPDeltas_exposer.def( "__str__", &__str__< ::SireAnalysis::FEPDeltas > );
        FEPDeltas_exposer.def( "__repr__", &__str__< ::SireAnalysis::FEPDeltas > );
    }

}
