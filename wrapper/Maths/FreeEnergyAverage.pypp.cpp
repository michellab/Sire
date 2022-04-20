// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "FreeEnergyAverage.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/histogram.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "freeenergyaverage.h"

#include <QDebug>

#include "freeenergyaverage.h"

SireMaths::FreeEnergyAverage __copy__(const SireMaths::FreeEnergyAverage &other){ return SireMaths::FreeEnergyAverage(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_FreeEnergyAverage_class(){

    { //::SireMaths::FreeEnergyAverage
        typedef bp::class_< SireMaths::FreeEnergyAverage, bp::bases< SireMaths::ExpAverage, SireMaths::Accumulator, SireBase::Property > > FreeEnergyAverage_exposer_t;
        FreeEnergyAverage_exposer_t FreeEnergyAverage_exposer = FreeEnergyAverage_exposer_t( "FreeEnergyAverage", "This class is used to accumulate the free energy\naverage. As well as calculating the average, it also\nrecords a histogram of values that can be used for\nerror analysis\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor - this defaults to accumulating the average\nat room temperature (25 C) and collects statistics about the\nfree energy using a histogram of bin width 0.5 kcal mol-1") );
        bp::scope FreeEnergyAverage_scope( FreeEnergyAverage_exposer );
        FreeEnergyAverage_exposer.def( bp::init< bool >(( bp::arg("forwards_free_energy") ), "Constructor - this defaults to accumulating the average\nat room temperature (25 C) and collects statistics about the\nfree energy using a histogram of bin width 0.5 kcal mol-1, specifying\nwhether or not this is a forwards free energy") );
        FreeEnergyAverage_exposer.def( bp::init< SireUnits::Dimension::Temperature const &, bp::optional< bool > >(( bp::arg("temperature"), bp::arg("forwards_free_energy")=(bool)(true) ), "Construct an accumulator to accumulate the free energy average\nat the specified temperature, and to collect statistics about\nthe free energy using a histogram of bin width 0.5 kcal mol-1") );
        FreeEnergyAverage_exposer.def( bp::init< SireUnits::Dimension::MolarEnergy const &, bp::optional< bool > >(( bp::arg("binwidth"), bp::arg("forwards_free_energy")=(bool)(true) ), "Constructor - this defaults to accumulating the average\nat room temperature (25 C) and collects statistics about the\nfree energy using a histogram of the passed bin width. If the binwidth\nis zero, then a histogram of energies is not collected") );
        FreeEnergyAverage_exposer.def( bp::init< SireUnits::Dimension::Temperature const &, SireUnits::Dimension::MolarEnergy const &, bp::optional< bool > >(( bp::arg("temperature"), bp::arg("binwidth"), bp::arg("forwards_free_energy")=(bool)(true) ), "Construct an accumulator to accumulate the free energy average\nat the specified temperature, and to collect statistics about\nthe free energy using a histogram of passed bin width") );
        FreeEnergyAverage_exposer.def( bp::init< SireMaths::FreeEnergyAverage const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMaths::FreeEnergyAverage::accumulate
        
            typedef void ( ::SireMaths::FreeEnergyAverage::*accumulate_function_type)( double ) ;
            accumulate_function_type accumulate_function_value( &::SireMaths::FreeEnergyAverage::accumulate );
            
            FreeEnergyAverage_exposer.def( 
                "accumulate"
                , accumulate_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Accumulate the passed energy difference onto the free energy average" );
        
        }
        { //::SireMaths::FreeEnergyAverage::average
        
            typedef double ( ::SireMaths::FreeEnergyAverage::*average_function_type)(  ) const;
            average_function_type average_function_value( &::SireMaths::FreeEnergyAverage::average );
            
            FreeEnergyAverage_exposer.def( 
                "average"
                , average_function_value
                , bp::release_gil_policy()
                , "Return the average free energy. Note that if this is a backwards free energy,\nthen this will return the negative (so that it is easy to combine backwards\nand forwards values)" );
        
        }
        { //::SireMaths::FreeEnergyAverage::average2
        
            typedef double ( ::SireMaths::FreeEnergyAverage::*average2_function_type)(  ) const;
            average2_function_type average2_function_value( &::SireMaths::FreeEnergyAverage::average2 );
            
            FreeEnergyAverage_exposer.def( 
                "average2"
                , average2_function_value
                , bp::release_gil_policy()
                , "Return the square of the average free energy. Note that if this is a backwards\nfree energy, then this will return the negative (so that it is easy to combine\nbackwards and forwards values)" );
        
        }
        { //::SireMaths::FreeEnergyAverage::clear
        
            typedef void ( ::SireMaths::FreeEnergyAverage::*clear_function_type)(  ) ;
            clear_function_type clear_function_value( &::SireMaths::FreeEnergyAverage::clear );
            
            FreeEnergyAverage_exposer.def( 
                "clear"
                , clear_function_value
                , bp::release_gil_policy()
                , "Clear all data from the accumulator" );
        
        }
        { //::SireMaths::FreeEnergyAverage::fepFreeEnergy
        
            typedef double ( ::SireMaths::FreeEnergyAverage::*fepFreeEnergy_function_type)(  ) const;
            fepFreeEnergy_function_type fepFreeEnergy_function_value( &::SireMaths::FreeEnergyAverage::fepFreeEnergy );
            
            FreeEnergyAverage_exposer.def( 
                "fepFreeEnergy"
                , fepFreeEnergy_function_value
                , bp::release_gil_policy()
                , "Return the average free energy. Note that if this is a backwards free energy,\nthen this will return the negative (so that it is easy to combine backwards\nand forwards values)" );
        
        }
        { //::SireMaths::FreeEnergyAverage::histogram
        
            typedef ::SireMaths::Histogram const & ( ::SireMaths::FreeEnergyAverage::*histogram_function_type)(  ) const;
            histogram_function_type histogram_function_value( &::SireMaths::FreeEnergyAverage::histogram );
            
            FreeEnergyAverage_exposer.def( 
                "histogram"
                , histogram_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the histogram of energies" );
        
        }
        { //::SireMaths::FreeEnergyAverage::isBackwardsFreeEnergy
        
            typedef bool ( ::SireMaths::FreeEnergyAverage::*isBackwardsFreeEnergy_function_type)(  ) const;
            isBackwardsFreeEnergy_function_type isBackwardsFreeEnergy_function_value( &::SireMaths::FreeEnergyAverage::isBackwardsFreeEnergy );
            
            FreeEnergyAverage_exposer.def( 
                "isBackwardsFreeEnergy"
                , isBackwardsFreeEnergy_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is a backwards free energy" );
        
        }
        { //::SireMaths::FreeEnergyAverage::isForwardsFreeEnergy
        
            typedef bool ( ::SireMaths::FreeEnergyAverage::*isForwardsFreeEnergy_function_type)(  ) const;
            isForwardsFreeEnergy_function_type isForwardsFreeEnergy_function_value( &::SireMaths::FreeEnergyAverage::isForwardsFreeEnergy );
            
            FreeEnergyAverage_exposer.def( 
                "isForwardsFreeEnergy"
                , isForwardsFreeEnergy_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is a forwards free energy" );
        
        }
        FreeEnergyAverage_exposer.def( bp::self != bp::self );
        FreeEnergyAverage_exposer.def( bp::self + bp::self );
        { //::SireMaths::FreeEnergyAverage::operator=
        
            typedef ::SireMaths::FreeEnergyAverage & ( ::SireMaths::FreeEnergyAverage::*assign_function_type)( ::SireMaths::FreeEnergyAverage const & ) ;
            assign_function_type assign_function_value( &::SireMaths::FreeEnergyAverage::operator= );
            
            FreeEnergyAverage_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        FreeEnergyAverage_exposer.def( bp::self == bp::self );
        { //::SireMaths::FreeEnergyAverage::taylorExpansion
        
            typedef double ( ::SireMaths::FreeEnergyAverage::*taylorExpansion_function_type)(  ) const;
            taylorExpansion_function_type taylorExpansion_function_value( &::SireMaths::FreeEnergyAverage::taylorExpansion );
            
            FreeEnergyAverage_exposer.def( 
                "taylorExpansion"
                , taylorExpansion_function_value
                , bp::release_gil_policy()
                , "Return the Taylor series expansion estimate the difference in free energy" );
        
        }
        { //::SireMaths::FreeEnergyAverage::temperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::SireMaths::FreeEnergyAverage::*temperature_function_type)(  ) const;
            temperature_function_type temperature_function_value( &::SireMaths::FreeEnergyAverage::temperature );
            
            FreeEnergyAverage_exposer.def( 
                "temperature"
                , temperature_function_value
                , bp::release_gil_policy()
                , "Return the temperature at which the free energy average\nis being accumulated" );
        
        }
        { //::SireMaths::FreeEnergyAverage::toString
        
            typedef ::QString ( ::SireMaths::FreeEnergyAverage::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::FreeEnergyAverage::toString );
            
            FreeEnergyAverage_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::FreeEnergyAverage::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::FreeEnergyAverage::typeName );
            
            FreeEnergyAverage_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        FreeEnergyAverage_exposer.staticmethod( "typeName" );
        FreeEnergyAverage_exposer.def( "__copy__", &__copy__);
        FreeEnergyAverage_exposer.def( "__deepcopy__", &__copy__);
        FreeEnergyAverage_exposer.def( "clone", &__copy__);
        FreeEnergyAverage_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::FreeEnergyAverage >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FreeEnergyAverage_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::FreeEnergyAverage >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FreeEnergyAverage_exposer.def_pickle(sire_pickle_suite< ::SireMaths::FreeEnergyAverage >());
        FreeEnergyAverage_exposer.def( "__str__", &__str__< ::SireMaths::FreeEnergyAverage > );
        FreeEnergyAverage_exposer.def( "__repr__", &__str__< ::SireMaths::FreeEnergyAverage > );
    }

}
