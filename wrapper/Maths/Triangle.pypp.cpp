// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Triangle.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "triangle.h"

#include <QDataStream>

#include "triangle.h"

SireMaths::Triangle __copy__(const SireMaths::Triangle &other){ return SireMaths::Triangle(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Triangle_class(){

    { //::SireMaths::Triangle
        typedef bp::class_< SireMaths::Triangle > Triangle_exposer_t;
        Triangle_exposer_t Triangle_exposer = Triangle_exposer_t( "Triangle", "\nThis class represents a triangle in three-dimensional space. (or three points)\n\nAuthor: Christopher Woods\n", bp::init< >("Create a zero triangle") );
        bp::scope Triangle_scope( Triangle_exposer );
        Triangle_exposer.def( bp::init< SireMaths::Vector const &, SireMaths::Vector const &, SireMaths::Vector const & >(( bp::arg("point0"), bp::arg("point1"), bp::arg("point2") ), "Create a triangle from point 0 -> 1 -> 2") );
        Triangle_exposer.def( bp::init< SireMaths::Triangle const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMaths::Triangle::angle
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMaths::Triangle::*angle_function_type)(  ) const;
            angle_function_type angle_function_value( &::SireMaths::Triangle::angle );
            
            Triangle_exposer.def( 
                "angle"
                , angle_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::angle0
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMaths::Triangle::*angle0_function_type)(  ) const;
            angle0_function_type angle0_function_value( &::SireMaths::Triangle::angle0 );
            
            Triangle_exposer.def( 
                "angle0"
                , angle0_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::angle1
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMaths::Triangle::*angle1_function_type)(  ) const;
            angle1_function_type angle1_function_value( &::SireMaths::Triangle::angle1 );
            
            Triangle_exposer.def( 
                "angle1"
                , angle1_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::angle2
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMaths::Triangle::*angle2_function_type)(  ) const;
            angle2_function_type angle2_function_value( &::SireMaths::Triangle::angle2 );
            
            Triangle_exposer.def( 
                "angle2"
                , angle2_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::at
        
            typedef ::SireMaths::Vector const & ( ::SireMaths::Triangle::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMaths::Triangle::at );
            
            Triangle_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMaths::Triangle::count
        
            typedef int ( ::SireMaths::Triangle::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMaths::Triangle::count );
            
            Triangle_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::line0
        
            typedef ::SireMaths::Line ( ::SireMaths::Triangle::*line0_function_type)(  ) const;
            line0_function_type line0_function_value( &::SireMaths::Triangle::line0 );
            
            Triangle_exposer.def( 
                "line0"
                , line0_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::line1
        
            typedef ::SireMaths::Line ( ::SireMaths::Triangle::*line1_function_type)(  ) const;
            line1_function_type line1_function_value( &::SireMaths::Triangle::line1 );
            
            Triangle_exposer.def( 
                "line1"
                , line1_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::line2
        
            typedef ::SireMaths::Line ( ::SireMaths::Triangle::*line2_function_type)(  ) const;
            line2_function_type line2_function_value( &::SireMaths::Triangle::line2 );
            
            Triangle_exposer.def( 
                "line2"
                , line2_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::operator[]
        
            typedef ::SireMaths::Vector const & ( ::SireMaths::Triangle::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMaths::Triangle::operator[] );
            
            Triangle_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMaths::Triangle::point
        
            typedef ::SireMaths::Vector const & ( ::SireMaths::Triangle::*point_function_type)( int ) const;
            point_function_type point_function_value( &::SireMaths::Triangle::point );
            
            Triangle_exposer.def( 
                "point"
                , point_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMaths::Triangle::toString
        
            typedef ::QString ( ::SireMaths::Triangle::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::Triangle::toString );
            
            Triangle_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of the triangle" );
        
        }
        { //::SireMaths::Triangle::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::Triangle::typeName );
            
            Triangle_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::vector
        
            typedef ::SireMaths::Vector ( ::SireMaths::Triangle::*vector_function_type)(  ) const;
            vector_function_type vector_function_value( &::SireMaths::Triangle::vector );
            
            Triangle_exposer.def( 
                "vector"
                , vector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::vector0
        
            typedef ::SireMaths::Vector ( ::SireMaths::Triangle::*vector0_function_type)(  ) const;
            vector0_function_type vector0_function_value( &::SireMaths::Triangle::vector0 );
            
            Triangle_exposer.def( 
                "vector0"
                , vector0_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::vector1
        
            typedef ::SireMaths::Vector ( ::SireMaths::Triangle::*vector1_function_type)(  ) const;
            vector1_function_type vector1_function_value( &::SireMaths::Triangle::vector1 );
            
            Triangle_exposer.def( 
                "vector1"
                , vector1_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::vector2
        
            typedef ::SireMaths::Vector ( ::SireMaths::Triangle::*vector2_function_type)(  ) const;
            vector2_function_type vector2_function_value( &::SireMaths::Triangle::vector2 );
            
            Triangle_exposer.def( 
                "vector2"
                , vector2_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::Triangle::what
        
            typedef char const * ( ::SireMaths::Triangle::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::Triangle::what );
            
            Triangle_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Triangle_exposer.staticmethod( "typeName" );
        Triangle_exposer.def( "__copy__", &__copy__);
        Triangle_exposer.def( "__deepcopy__", &__copy__);
        Triangle_exposer.def( "clone", &__copy__);
        Triangle_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::Triangle >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Triangle_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::Triangle >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Triangle_exposer.def_pickle(sire_pickle_suite< ::SireMaths::Triangle >());
        Triangle_exposer.def( "__str__", &__str__< ::SireMaths::Triangle > );
        Triangle_exposer.def( "__repr__", &__str__< ::SireMaths::Triangle > );
        Triangle_exposer.def( "__len__", &__len_count< ::SireMaths::Triangle > );
    }

}
