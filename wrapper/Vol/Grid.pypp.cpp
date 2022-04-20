// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Grid.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/rotate.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "grid.h"

#include "grid.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Grid_class(){

    { //::SireVol::Grid
        typedef bp::class_< SireVol::Grid, bp::bases< SireBase::Property >, boost::noncopyable > Grid_exposer_t;
        Grid_exposer_t Grid_exposer = Grid_exposer_t( "Grid", "This is the base class of all grids. A grid is a set of points\nin space which are laid out in some manner (e.g. regularly spaced,\nor randomly distributed).\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Grid_scope( Grid_exposer );
        { //::SireVol::Grid::aaBox
        
            typedef ::SireVol::AABox const & ( ::SireVol::Grid::*aaBox_function_type)(  ) const;
            aaBox_function_type aaBox_function_value( &::SireVol::Grid::aaBox );
            
            Grid_exposer.def( 
                "aaBox"
                , aaBox_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireVol::Grid::center
        
            typedef ::SireMaths::Vector ( ::SireVol::Grid::*center_function_type)(  ) const;
            center_function_type center_function_value( &::SireVol::Grid::center );
            
            Grid_exposer.def( 
                "center"
                , center_function_value
                , bp::release_gil_policy()
                , "Return the center of the grid" );
        
        }
        { //::SireVol::Grid::count
        
            typedef int ( ::SireVol::Grid::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireVol::Grid::count );
            
            Grid_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of points in the grid" );
        
        }
        { //::SireVol::Grid::hasWeights
        
            typedef bool ( ::SireVol::Grid::*hasWeights_function_type)(  ) const;
            hasWeights_function_type hasWeights_function_value( &::SireVol::Grid::hasWeights );
            
            Grid_exposer.def( 
                "hasWeights"
                , hasWeights_function_value
                , bp::release_gil_policy()
                , "Return whether or not the grid points have weights\n(i.e. they are not equally weighted)" );
        
        }
        { //::SireVol::Grid::isEmpty
        
            typedef bool ( ::SireVol::Grid::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireVol::Grid::isEmpty );
            
            Grid_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether this grid is empty" );
        
        }
        { //::SireVol::Grid::maxCoords
        
            typedef ::SireMaths::Vector ( ::SireVol::Grid::*maxCoords_function_type)(  ) const;
            maxCoords_function_type maxCoords_function_value( &::SireVol::Grid::maxCoords );
            
            Grid_exposer.def( 
                "maxCoords"
                , maxCoords_function_value
                , bp::release_gil_policy()
                , "Return the maximum coordinates of the grid" );
        
        }
        { //::SireVol::Grid::minCoords
        
            typedef ::SireMaths::Vector ( ::SireVol::Grid::*minCoords_function_type)(  ) const;
            minCoords_function_type minCoords_function_value( &::SireVol::Grid::minCoords );
            
            Grid_exposer.def( 
                "minCoords"
                , minCoords_function_value
                , bp::release_gil_policy()
                , "Return the minimum coordinates of the grid" );
        
        }
        { //::SireVol::Grid::nPoints
        
            typedef int ( ::SireVol::Grid::*nPoints_function_type)(  ) const;
            nPoints_function_type nPoints_function_value( &::SireVol::Grid::nPoints );
            
            Grid_exposer.def( 
                "nPoints"
                , nPoints_function_value
                , bp::release_gil_policy()
                , "Return the number of points in the grid" );
        
        }
        { //::SireVol::Grid::null
        
            typedef ::SireVol::NullGrid const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireVol::Grid::null );
            
            Grid_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireVol::Grid::points
        
            typedef ::QVector< SireMaths::Vector > const & ( ::SireVol::Grid::*points_function_type)(  ) const;
            points_function_type points_function_value( &::SireVol::Grid::points );
            
            Grid_exposer.def( 
                "points"
                , points_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the array of grid points" );
        
        }
        { //::SireVol::Grid::recenter
        
            typedef ::SireVol::GridPtr ( ::SireVol::Grid::*recenter_function_type)( ::SireMaths::Vector const & ) const;
            recenter_function_type recenter_function_value( &::SireVol::Grid::recenter );
            
            Grid_exposer.def( 
                "recenter"
                , recenter_function_value
                , ( bp::arg("center") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Grid::rotate
        
            typedef ::SireVol::GridPtr ( ::SireVol::Grid::*rotate_function_type)( ::SireMaths::Matrix const &,::SireMaths::Vector const & ) const;
            rotate_function_type rotate_function_value( &::SireVol::Grid::rotate );
            
            Grid_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("rotmat"), bp::arg("center")=SireMaths::Vector(0) )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Grid::rotate
        
            typedef ::SireVol::GridPtr ( ::SireVol::Grid::*rotate_function_type)( ::SireMaths::Quaternion const &,::SireMaths::Vector const & ) const;
            rotate_function_type rotate_function_value( &::SireVol::Grid::rotate );
            
            Grid_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("quat"), bp::arg("center")=SireMaths::Vector(0) )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Grid::scale
        
            typedef ::SireVol::GridPtr ( ::SireVol::Grid::*scale_function_type)( double ) const;
            scale_function_type scale_function_value( &::SireVol::Grid::scale );
            
            Grid_exposer.def( 
                "scale"
                , scale_function_value
                , ( bp::arg("scalefactor") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Grid::translate
        
            typedef ::SireVol::GridPtr ( ::SireVol::Grid::*translate_function_type)( ::SireMaths::Vector const & ) const;
            translate_function_type translate_function_value( &::SireVol::Grid::translate );
            
            Grid_exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("delta") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Grid::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireVol::Grid::typeName );
            
            Grid_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Grid::weights
        
            typedef ::QVector< double > const & ( ::SireVol::Grid::*weights_function_type)(  ) const;
            weights_function_type weights_function_value( &::SireVol::Grid::weights );
            
            Grid_exposer.def( 
                "weights"
                , weights_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the array of grid weights - this will be empty\nif all of the points are equally weighted" );
        
        }
        Grid_exposer.staticmethod( "null" );
        Grid_exposer.staticmethod( "typeName" );
        Grid_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireVol::Grid >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Grid_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireVol::Grid >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Grid_exposer.def_pickle(sire_pickle_suite< ::SireVol::Grid >());
        Grid_exposer.def( "__str__", &__str__< ::SireVol::Grid > );
        Grid_exposer.def( "__repr__", &__str__< ::SireVol::Grid > );
        Grid_exposer.def( "__len__", &__len_count< ::SireVol::Grid > );
    }

}
