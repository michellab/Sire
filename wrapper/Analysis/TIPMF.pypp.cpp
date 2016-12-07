// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "TIPMF.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/registeralternativename.h"

#include "SireStream/shareddatastream.h"

#include "ti.h"

#include "tostring.h"

#include <cmath>

#include "ti.h"

SireAnalysis::TIPMF __copy__(const SireAnalysis::TIPMF &other){ return SireAnalysis::TIPMF(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_TIPMF_class(){

    { //::SireAnalysis::TIPMF
        typedef bp::class_< SireAnalysis::TIPMF, bp::bases< SireAnalysis::PMF, SireBase::Property > > TIPMF_exposer_t;
        TIPMF_exposer_t TIPMF_exposer = TIPMF_exposer_t( "TIPMF", "This class contains the complete potential of mean force\nthat has been created by integrating the TI gradients\n\nAuthor: Christopher Woods\n", bp::init< >("Construct a PMF that will use 10 polynomials to fit\nand integrate the gradients between 0 and 1") );
        bp::scope TIPMF_scope( TIPMF_exposer );
        TIPMF_exposer.def( bp::init< int >(( bp::arg("order") ), "Construct a PMF that will use the passed number of polynomials\nto fit and integrate the gradients between 0 and 1") );
        TIPMF_exposer.def( bp::init< double, double >(( bp::arg("min_x"), bp::arg("max_x") ), "Construct a PMF that will use 10 polynomials to fit and integrate\nthe gradients in the passed range") );
        TIPMF_exposer.def( bp::init< double, double, int >(( bp::arg("min_x"), bp::arg("max_x"), bp::arg("order") ), "Construct a PMF that will use the passed number of polynomials to fit and\nintegrate the gradients in the passed range") );
        TIPMF_exposer.def( bp::init< SireAnalysis::TIPMF const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireAnalysis::TIPMF::dropEndPoints
        
            typedef ::SireAnalysis::TIPMF ( ::SireAnalysis::TIPMF::*dropEndPoints_function_type)(  ) const;
            dropEndPoints_function_type dropEndPoints_function_value( &::SireAnalysis::TIPMF::dropEndPoints );
            
            TIPMF_exposer.def( 
                "dropEndPoints"
                , dropEndPoints_function_value
                , "Return a copy of the PMF where the gradients at the end points\n(the first and last gradients) have been removed. This can be used\nto estimate the effect of end-point error" );
        
        }
        { //::SireAnalysis::TIPMF::gradients
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::TIPMF::*gradients_function_type)(  ) const;
            gradients_function_type gradients_function_value( &::SireAnalysis::TIPMF::gradients );
            
            TIPMF_exposer.def( 
                "gradients"
                , gradients_function_value
                , "Return the raw gradients used to calculate the PMF" );
        
        }
        { //::SireAnalysis::TIPMF::integral
        
            typedef double ( ::SireAnalysis::TIPMF::*integral_function_type)(  ) const;
            integral_function_type integral_function_value( &::SireAnalysis::TIPMF::integral );
            
            TIPMF_exposer.def( 
                "integral"
                , integral_function_value
                , "Return the free energy calculated using integration of the\npolynomial fitted to the gradients" );
        
        }
        TIPMF_exposer.def( bp::self != bp::self );
        { //::SireAnalysis::TIPMF::operator=
        
            typedef ::SireAnalysis::TIPMF & ( ::SireAnalysis::TIPMF::*assign_function_type)( ::SireAnalysis::TIPMF const & ) ;
            assign_function_type assign_function_value( &::SireAnalysis::TIPMF::operator= );
            
            TIPMF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        TIPMF_exposer.def( bp::self == bp::self );
        { //::SireAnalysis::TIPMF::order
        
            typedef int ( ::SireAnalysis::TIPMF::*order_function_type)(  ) const;
            order_function_type order_function_value( &::SireAnalysis::TIPMF::order );
            
            TIPMF_exposer.def( 
                "order"
                , order_function_value
                , "Return the order (number of polynomials) used to integrate\nthe gradients to get the PMF" );
        
        }
        { //::SireAnalysis::TIPMF::quadrature
        
            typedef double ( ::SireAnalysis::TIPMF::*quadrature_function_type)(  ) const;
            quadrature_function_type quadrature_function_value( &::SireAnalysis::TIPMF::quadrature );
            
            TIPMF_exposer.def( 
                "quadrature"
                , quadrature_function_value
                , "Return the free energy calculated using trapezium quadrature\nfrom the raw gradients" );
        
        }
        { //::SireAnalysis::TIPMF::rangeMax
        
            typedef double ( ::SireAnalysis::TIPMF::*rangeMax_function_type)(  ) const;
            rangeMax_function_type rangeMax_function_value( &::SireAnalysis::TIPMF::rangeMax );
            
            TIPMF_exposer.def( 
                "rangeMax"
                , rangeMax_function_value
                , "Return the maximum value of the range of integration" );
        
        }
        { //::SireAnalysis::TIPMF::rangeMin
        
            typedef double ( ::SireAnalysis::TIPMF::*rangeMin_function_type)(  ) const;
            rangeMin_function_type rangeMin_function_value( &::SireAnalysis::TIPMF::rangeMin );
            
            TIPMF_exposer.def( 
                "rangeMin"
                , rangeMin_function_value
                , "Return the minimum value of the range of integration" );
        
        }
        { //::SireAnalysis::TIPMF::setGradients
        
            typedef void ( ::SireAnalysis::TIPMF::*setGradients_function_type)( ::QVector< SireAnalysis::DataPoint > const & ) ;
            setGradients_function_type setGradients_function_value( &::SireAnalysis::TIPMF::setGradients );
            
            TIPMF_exposer.def( 
                "setGradients"
                , setGradients_function_value
                , ( bp::arg("gradients") )
                , "Set the raw gradients to be integrated" );
        
        }
        { //::SireAnalysis::TIPMF::setOrder
        
            typedef void ( ::SireAnalysis::TIPMF::*setOrder_function_type)( ::qint32 ) ;
            setOrder_function_type setOrder_function_value( &::SireAnalysis::TIPMF::setOrder );
            
            TIPMF_exposer.def( 
                "setOrder"
                , setOrder_function_value
                , ( bp::arg("order") )
                , "Set the order (number of polynomials) to fit the gradients for\nPMF integration" );
        
        }
        { //::SireAnalysis::TIPMF::setRange
        
            typedef void ( ::SireAnalysis::TIPMF::*setRange_function_type)( double,double ) ;
            setRange_function_type setRange_function_value( &::SireAnalysis::TIPMF::setRange );
            
            TIPMF_exposer.def( 
                "setRange"
                , setRange_function_value
                , ( bp::arg("min_x"), bp::arg("max_x") )
                , "Set the range of integration" );
        
        }
        { //::SireAnalysis::TIPMF::smoothedGradients
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::TIPMF::*smoothedGradients_function_type)(  ) const;
            smoothedGradients_function_type smoothedGradients_function_value( &::SireAnalysis::TIPMF::smoothedGradients );
            
            TIPMF_exposer.def( 
                "smoothedGradients"
                , smoothedGradients_function_value
                , "Return the smoothed (fitted) gradients used to calculate the PMF" );
        
        }
        { //::SireAnalysis::TIPMF::toString
        
            typedef ::QString ( ::SireAnalysis::TIPMF::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireAnalysis::TIPMF::toString );
            
            TIPMF_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireAnalysis::TIPMF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireAnalysis::TIPMF::typeName );
            
            TIPMF_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireAnalysis::TIPMF::what
        
            typedef char const * ( ::SireAnalysis::TIPMF::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireAnalysis::TIPMF::what );
            
            TIPMF_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        TIPMF_exposer.staticmethod( "typeName" );
        TIPMF_exposer.def( "__copy__", &__copy__);
        TIPMF_exposer.def( "__deepcopy__", &__copy__);
        TIPMF_exposer.def( "clone", &__copy__);
        TIPMF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireAnalysis::TIPMF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TIPMF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireAnalysis::TIPMF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TIPMF_exposer.def( "__str__", &__str__< ::SireAnalysis::TIPMF > );
        TIPMF_exposer.def( "__repr__", &__str__< ::SireAnalysis::TIPMF > );
    }

}
