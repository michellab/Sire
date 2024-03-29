// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "HistogramValue.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "histogram.h"

#include "histogram.h"

SireMaths::HistogramValue __copy__(const SireMaths::HistogramValue &other){ return SireMaths::HistogramValue(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_HistogramValue_class(){

    { //::SireMaths::HistogramValue
        typedef bp::class_< SireMaths::HistogramValue, bp::bases< SireMaths::HistogramBin > > HistogramValue_exposer_t;
        HistogramValue_exposer_t HistogramValue_exposer = HistogramValue_exposer_t( "HistogramValue", "This class represents a single histogram bin with its associated value", bp::init< >("Null constructor") );
        bp::scope HistogramValue_scope( HistogramValue_exposer );
        HistogramValue_exposer.def( bp::init< SireMaths::HistogramBin const &, double >(( bp::arg("bin"), bp::arg("value") ), "Construct the value for the bin bin equal to value") );
        HistogramValue_exposer.def( bp::init< SireMaths::HistogramValue const & >(( bp::arg("other") ), "Copy constructor") );
        HistogramValue_exposer.def( bp::self != bp::self );
        { //::SireMaths::HistogramValue::operator=
        
            typedef ::SireMaths::HistogramValue & ( ::SireMaths::HistogramValue::*assign_function_type)( ::SireMaths::HistogramValue const & ) ;
            assign_function_type assign_function_value( &::SireMaths::HistogramValue::operator= );
            
            HistogramValue_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        HistogramValue_exposer.def( bp::self == bp::self );
        { //::SireMaths::HistogramValue::toString
        
            typedef ::QString ( ::SireMaths::HistogramValue::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::HistogramValue::toString );
            
            HistogramValue_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation" );
        
        }
        { //::SireMaths::HistogramValue::value
        
            typedef double ( ::SireMaths::HistogramValue::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireMaths::HistogramValue::value );
            
            HistogramValue_exposer.def( 
                "value"
                , value_function_value
                , bp::release_gil_policy()
                , "Return the value of the bin" );
        
        }
        HistogramValue_exposer.def( "__copy__", &__copy__);
        HistogramValue_exposer.def( "__deepcopy__", &__copy__);
        HistogramValue_exposer.def( "clone", &__copy__);
        HistogramValue_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::HistogramValue >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        HistogramValue_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::HistogramValue >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        HistogramValue_exposer.def_pickle(sire_pickle_suite< ::SireMaths::HistogramValue >());
        HistogramValue_exposer.def( "__str__", &__str__< ::SireMaths::HistogramValue > );
        HistogramValue_exposer.def( "__repr__", &__str__< ::SireMaths::HistogramValue > );
    }

}
