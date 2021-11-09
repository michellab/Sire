// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CLJBoxDistance.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMM/cljboxes.h"

#include "SireMM/cljdelta.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/aabox.h"

#include "SireVol/periodicbox.h"

#include "SireVol/space.h"

#include "SireVol/triclinicbox.h"

#include "cljboxes.h"

#include "tostring.h"

#include <QDebug>

#include <QElapsedTimer>

#include "cljboxes.h"

SireMM::CLJBoxDistance __copy__(const SireMM::CLJBoxDistance &other){ return SireMM::CLJBoxDistance(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CLJBoxDistance_class(){

    { //::SireMM::CLJBoxDistance
        typedef bp::class_< SireMM::CLJBoxDistance > CLJBoxDistance_exposer_t;
        CLJBoxDistance_exposer_t CLJBoxDistance_exposer = CLJBoxDistance_exposer_t( "CLJBoxDistance", "This simple class holds the minimum distance between the two\nCLJBoxes at specified indicies\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope CLJBoxDistance_scope( CLJBoxDistance_exposer );
        CLJBoxDistance_exposer.def( bp::init< quint32, quint32, float >(( bp::arg("box0"), bp::arg("box1"), bp::arg("distance") ), "Construct saying that the minimum distance between the box with index box0\nand the box with index box1 is distance") );
        CLJBoxDistance_exposer.def( bp::init< SireMM::CLJBoxDistance const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::CLJBoxDistance::box0
        
            typedef ::quint32 ( ::SireMM::CLJBoxDistance::*box0_function_type)(  ) const;
            box0_function_type box0_function_value( &::SireMM::CLJBoxDistance::box0 );
            
            CLJBoxDistance_exposer.def( 
                "box0"
                , box0_function_value
                , "" );
        
        }
        { //::SireMM::CLJBoxDistance::box1
        
            typedef ::quint32 ( ::SireMM::CLJBoxDistance::*box1_function_type)(  ) const;
            box1_function_type box1_function_value( &::SireMM::CLJBoxDistance::box1 );
            
            CLJBoxDistance_exposer.def( 
                "box1"
                , box1_function_value
                , "" );
        
        }
        { //::SireMM::CLJBoxDistance::distance
        
            typedef float ( ::SireMM::CLJBoxDistance::*distance_function_type)(  ) const;
            distance_function_type distance_function_value( &::SireMM::CLJBoxDistance::distance );
            
            CLJBoxDistance_exposer.def( 
                "distance"
                , distance_function_value
                , "" );
        
        }
        CLJBoxDistance_exposer.def( bp::self != bp::self );
        CLJBoxDistance_exposer.def( bp::self < bp::self );
        { //::SireMM::CLJBoxDistance::operator=
        
            typedef ::SireMM::CLJBoxDistance & ( ::SireMM::CLJBoxDistance::*assign_function_type)( ::SireMM::CLJBoxDistance const & ) ;
            assign_function_type assign_function_value( &::SireMM::CLJBoxDistance::operator= );
            
            CLJBoxDistance_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CLJBoxDistance_exposer.def( bp::self == bp::self );
        CLJBoxDistance_exposer.def( bp::self > bp::self );
        { //::SireMM::CLJBoxDistance::toString
        
            typedef ::QString ( ::SireMM::CLJBoxDistance::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::CLJBoxDistance::toString );
            
            CLJBoxDistance_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMM::CLJBoxDistance::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::CLJBoxDistance::typeName );
            
            CLJBoxDistance_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::CLJBoxDistance::what
        
            typedef char const * ( ::SireMM::CLJBoxDistance::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::CLJBoxDistance::what );
            
            CLJBoxDistance_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        CLJBoxDistance_exposer.staticmethod( "typeName" );
        CLJBoxDistance_exposer.def( "__copy__", &__copy__);
        CLJBoxDistance_exposer.def( "__deepcopy__", &__copy__);
        CLJBoxDistance_exposer.def( "clone", &__copy__);
        CLJBoxDistance_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJBoxDistance >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJBoxDistance_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJBoxDistance >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJBoxDistance_exposer.def( "__str__", &__str__< ::SireMM::CLJBoxDistance > );
        CLJBoxDistance_exposer.def( "__repr__", &__str__< ::SireMM::CLJBoxDistance > );
    }

}
