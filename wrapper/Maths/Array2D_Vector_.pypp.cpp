// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Array2D_Vector_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/array2d.hpp"

#include "SireBase/trigarray2d.hpp"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "matrix.h"

#include "nmatrix.h"

#include "nvector.h"

#include "sire_blas.h"

#include "sire_lapack.h"

#include "sire_linpack.h"

#include "trigmatrix.h"

#include "vector.h"

#include "nmatrix.h"

SireBase::Array2D<SireMaths::Vector> __copy__(const SireBase::Array2D<SireMaths::Vector> &other){ return SireBase::Array2D<SireMaths::Vector>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Array2D_Vector__class(){

    { //::SireBase::Array2D< SireMaths::Vector >
        typedef bp::class_< SireBase::Array2D< SireMaths::Vector >, bp::bases< SireBase::Array2DBase > > Array2D_Vector__exposer_t;
        Array2D_Vector__exposer_t Array2D_Vector__exposer = Array2D_Vector__exposer_t( "Array2D_Vector_", "", bp::init< >("") );
        bp::scope Array2D_Vector__scope( Array2D_Vector__exposer );
        Array2D_Vector__exposer.def( bp::init< int, int >(( bp::arg("nrows"), bp::arg("ncolumns") ), "") );
        Array2D_Vector__exposer.def( bp::init< int, int, SireMaths::Vector const & >(( bp::arg("nrows"), bp::arg("ncolumns"), bp::arg("default_value") ), "") );
        Array2D_Vector__exposer.def( bp::init< SireBase::Array2D< SireMaths::Vector > const & >(( bp::arg("other") ), "") );
        { //::SireBase::Array2D< SireMaths::Vector >::at
        
            typedef SireBase::Array2D< SireMaths::Vector > exported_class_t;
            typedef ::SireMaths::Vector const & ( ::SireBase::Array2D< SireMaths::Vector >::*at_function_type)( int,int ) const;
            at_function_type at_function_value( &::SireBase::Array2D< SireMaths::Vector >::at );
            
            Array2D_Vector__exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::Vector >::get
        
            typedef SireBase::Array2D< SireMaths::Vector > exported_class_t;
            typedef ::SireMaths::Vector const & ( ::SireBase::Array2D< SireMaths::Vector >::*get_function_type)( int,int ) const;
            get_function_type get_function_value( &::SireBase::Array2D< SireMaths::Vector >::get );
            
            Array2D_Vector__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        Array2D_Vector__exposer.def( bp::self != bp::self );
        { //::SireBase::Array2D< SireMaths::Vector >::operator()
        
            typedef SireBase::Array2D< SireMaths::Vector > exported_class_t;
            typedef ::SireMaths::Vector const & ( ::SireBase::Array2D< SireMaths::Vector >::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireBase::Array2D< SireMaths::Vector >::operator() );
            
            Array2D_Vector__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::Vector >::operator=
        
            typedef SireBase::Array2D< SireMaths::Vector > exported_class_t;
            typedef ::SireBase::Array2D< SireMaths::Vector > & ( ::SireBase::Array2D< SireMaths::Vector >::*assign_function_type)( ::SireBase::Array2D< SireMaths::Vector > const & ) ;
            assign_function_type assign_function_value( &::SireBase::Array2D< SireMaths::Vector >::operator= );
            
            Array2D_Vector__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Array2D_Vector__exposer.def( bp::self == bp::self );
        { //::SireBase::Array2D< SireMaths::Vector >::redimension
        
            typedef SireBase::Array2D< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::Vector >::*redimension_function_type)( int,int ) ;
            redimension_function_type redimension_function_value( &::SireBase::Array2D< SireMaths::Vector >::redimension );
            
            Array2D_Vector__exposer.def( 
                "redimension"
                , redimension_function_value
                , ( bp::arg("nrows"), bp::arg("ncolumns") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::Vector >::set
        
            typedef SireBase::Array2D< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::Vector >::*set_function_type)( int,int,::SireMaths::Vector const & ) ;
            set_function_type set_function_value( &::SireBase::Array2D< SireMaths::Vector >::set );
            
            Array2D_Vector__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::Vector >::setAll
        
            typedef SireBase::Array2D< SireMaths::Vector > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::Vector >::*setAll_function_type)( ::SireMaths::Vector const & ) ;
            setAll_function_type setAll_function_value( &::SireBase::Array2D< SireMaths::Vector >::setAll );
            
            Array2D_Vector__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::Vector >::toString
        
            typedef SireBase::Array2D< SireMaths::Vector > exported_class_t;
            typedef ::QString ( ::SireBase::Array2D< SireMaths::Vector >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::Array2D< SireMaths::Vector >::toString );
            
            Array2D_Vector__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::Vector >::transpose
        
            typedef SireBase::Array2D< SireMaths::Vector > exported_class_t;
            typedef ::SireBase::Array2D< SireMaths::Vector > ( ::SireBase::Array2D< SireMaths::Vector >::*transpose_function_type)(  ) const;
            transpose_function_type transpose_function_value( &::SireBase::Array2D< SireMaths::Vector >::transpose );
            
            Array2D_Vector__exposer.def( 
                "transpose"
                , transpose_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Array2D_Vector__exposer.def( "__copy__", &__copy__);
        Array2D_Vector__exposer.def( "__deepcopy__", &__copy__);
        Array2D_Vector__exposer.def( "clone", &__copy__);
        Array2D_Vector__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::Array2D<SireMaths::Vector> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Array2D_Vector__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::Array2D<SireMaths::Vector> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Array2D_Vector__exposer.def_pickle(sire_pickle_suite< ::SireBase::Array2D<SireMaths::Vector> >());
        Array2D_Vector__exposer.def( "__str__", &__str__< ::SireBase::Array2D<SireMaths::Vector> > );
        Array2D_Vector__exposer.def( "__repr__", &__str__< ::SireBase::Array2D<SireMaths::Vector> > );
    }

}
