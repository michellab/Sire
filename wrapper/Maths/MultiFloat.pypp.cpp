// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MultiFloat.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/unittest.h"

#include "SireError/errors.h"

#include "multidouble.h"

#include "multifloat.h"

#include "multiint.h"

#include "multiuint.h"

#include "sincos.h"

#include <QDebug>

#include <QStringList>

#include "multifloat.h"

#include "multifloat.h"

#include "multiint.h"

#include "multidouble.h"

#include "multivector.h"

#include "multiquaternion.h"

SireMaths::MultiFloat __copy__(const SireMaths::MultiFloat &other){ return SireMaths::MultiFloat(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_MultiFloat_class(){

    { //::SireMaths::MultiFloat
        typedef bp::class_< SireMaths::MultiFloat > MultiFloat_exposer_t;
        MultiFloat_exposer_t MultiFloat_exposer = MultiFloat_exposer_t( "MultiFloat", "This class provides a vectorised float. This represents\na single vector of floats on the compiled machine, e.g.\n4 floats if we use SSE2, 8 floats for AVX\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope MultiFloat_scope( MultiFloat_exposer );
        MultiFloat_exposer.def( bp::init< float >(( bp::arg("value") ), "Construct from the passed array - this must be the same size as the vector") );
        MultiFloat_exposer.def( bp::init< double >(( bp::arg("value") ), "Construct from the passed array - this must be the same size as the vector") );
        MultiFloat_exposer.def( bp::init< int >(( bp::arg("value") ), "Construct from a MultiInt") );
        MultiFloat_exposer.def( bp::init< float const *, int >(( bp::arg("array"), bp::arg("size") ), "Construct from the passed array. If size is greater than MultiFloat::size()\nthen an error will be raised. If size is less than MultiFloat::size() then\nthis float will be padded with zeroes") );
        MultiFloat_exposer.def( bp::init< float const *, SireMaths::MultiInt const & >(( bp::arg("array"), bp::arg("indicies") ), "Construct from the passed array, taking the values of each element\nof the vector from the index in the associated MultiInt, e.g.\nMultiFloat[i] = array[ MultiInt[i] ]\n") );
        MultiFloat_exposer.def( bp::init< QVector< float > const & >(( bp::arg("array") ), "Construct from the passed array - this must be the same size as the vector") );
        MultiFloat_exposer.def( bp::init< QVector< double > const & >(( bp::arg("array") ), "Construct from the passed array - this must be the same size as the vector") );
        MultiFloat_exposer.def( bp::init< SireMaths::MultiDouble const & >(( bp::arg("other") ), "Construct from a MultiInt") );
        MultiFloat_exposer.def( bp::init< SireMaths::MultiFloat const & >(( bp::arg("other") ), "Construct from a MultiInt") );
        MultiFloat_exposer.def( bp::init< SireMaths::MultiInt const & >(( bp::arg("other") ), "Construct from a MultiInt") );
        { //::SireMaths::MultiFloat::abs
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*abs_function_type)(  ) const;
            abs_function_type abs_function_value( &::SireMaths::MultiFloat::abs );
            
            MultiFloat_exposer.def( 
                "abs"
                , abs_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::at
        
            typedef float ( ::SireMaths::MultiFloat::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMaths::MultiFloat::at );
            
            MultiFloat_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::compareEqual
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*compareEqual_function_type)( ::SireMaths::MultiFloat const & ) const;
            compareEqual_function_type compareEqual_function_value( &::SireMaths::MultiFloat::compareEqual );
            
            MultiFloat_exposer.def( 
                "compareEqual"
                , compareEqual_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::compareGreater
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*compareGreater_function_type)( ::SireMaths::MultiFloat const & ) const;
            compareGreater_function_type compareGreater_function_value( &::SireMaths::MultiFloat::compareGreater );
            
            MultiFloat_exposer.def( 
                "compareGreater"
                , compareGreater_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::compareGreaterEqual
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*compareGreaterEqual_function_type)( ::SireMaths::MultiFloat const & ) const;
            compareGreaterEqual_function_type compareGreaterEqual_function_value( &::SireMaths::MultiFloat::compareGreaterEqual );
            
            MultiFloat_exposer.def( 
                "compareGreaterEqual"
                , compareGreaterEqual_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::compareLess
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*compareLess_function_type)( ::SireMaths::MultiFloat const & ) const;
            compareLess_function_type compareLess_function_value( &::SireMaths::MultiFloat::compareLess );
            
            MultiFloat_exposer.def( 
                "compareLess"
                , compareLess_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::compareLessEqual
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*compareLessEqual_function_type)( ::SireMaths::MultiFloat const & ) const;
            compareLessEqual_function_type compareLessEqual_function_value( &::SireMaths::MultiFloat::compareLessEqual );
            
            MultiFloat_exposer.def( 
                "compareLessEqual"
                , compareLessEqual_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::compareNotEqual
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*compareNotEqual_function_type)( ::SireMaths::MultiFloat const & ) const;
            compareNotEqual_function_type compareNotEqual_function_value( &::SireMaths::MultiFloat::compareNotEqual );
            
            MultiFloat_exposer.def( 
                "compareNotEqual"
                , compareNotEqual_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::count
        
            typedef int ( *count_function_type )(  );
            count_function_type count_function_value( &::SireMaths::MultiFloat::count );
            
            MultiFloat_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::doubleSum
        
            typedef double ( ::SireMaths::MultiFloat::*doubleSum_function_type)(  ) const;
            doubleSum_function_type doubleSum_function_value( &::SireMaths::MultiFloat::doubleSum );
            
            MultiFloat_exposer.def( 
                "doubleSum"
                , doubleSum_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::fromArray
        
            typedef ::QVector< SireMaths::MultiFloat > ( *fromArray_function_type )( ::QVector< double > const & );
            fromArray_function_type fromArray_function_value( &::SireMaths::MultiFloat::fromArray );
            
            MultiFloat_exposer.def( 
                "fromArray"
                , fromArray_function_value
                , ( bp::arg("array") )
                , bp::release_gil_policy()
                , "Create an array of MultiFloats from the passed array of doubles. This\nwill pad the end of the array with zeroes if necessary" );
        
        }
        { //::SireMaths::MultiFloat::fromArray
        
            typedef ::QVector< SireMaths::MultiFloat > ( *fromArray_function_type )( ::QVector< float > const & );
            fromArray_function_type fromArray_function_value( &::SireMaths::MultiFloat::fromArray );
            
            MultiFloat_exposer.def( 
                "fromArray"
                , fromArray_function_value
                , ( bp::arg("array") )
                , bp::release_gil_policy()
                , "Create an array of MultiFloats from the passed array of floats. This will\npad the end of the array with zeroes if necessary" );
        
        }
        { //::SireMaths::MultiFloat::fromArray
        
            typedef ::QVector< SireMaths::MultiFloat > ( *fromArray_function_type )( double const *,int );
            fromArray_function_type fromArray_function_value( &::SireMaths::MultiFloat::fromArray );
            
            MultiFloat_exposer.def( 
                "fromArray"
                , fromArray_function_value
                , ( bp::arg("array"), bp::arg("size") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::fromArray
        
            typedef ::QVector< SireMaths::MultiFloat > ( *fromArray_function_type )( float const *,int );
            fromArray_function_type fromArray_function_value( &::SireMaths::MultiFloat::fromArray );
            
            MultiFloat_exposer.def( 
                "fromArray"
                , fromArray_function_value
                , ( bp::arg("array"), bp::arg("size") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::get
        
            typedef float ( ::SireMaths::MultiFloat::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMaths::MultiFloat::get );
            
            MultiFloat_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return the ith value in the multifloat" );
        
        }
        { //::SireMaths::MultiFloat::getitem
        
            typedef float ( ::SireMaths::MultiFloat::*getitem_function_type)( int ) const;
            getitem_function_type getitem_function_value( &::SireMaths::MultiFloat::getitem );
            
            MultiFloat_exposer.def( 
                "getitem"
                , getitem_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::hasBinaryOne
        
            typedef bool ( ::SireMaths::MultiFloat::*hasBinaryOne_function_type)(  ) const;
            hasBinaryOne_function_type hasBinaryOne_function_value( &::SireMaths::MultiFloat::hasBinaryOne );
            
            MultiFloat_exposer.def( 
                "hasBinaryOne"
                , hasBinaryOne_function_value
                , bp::release_gil_policy()
                , "Return whether or not at least one of the elements of this vector\nis binary one (the float is equal to 0xFFFFFFFF)" );
        
        }
        { //::SireMaths::MultiFloat::hasBinaryZero
        
            typedef bool ( ::SireMaths::MultiFloat::*hasBinaryZero_function_type)(  ) const;
            hasBinaryZero_function_type hasBinaryZero_function_value( &::SireMaths::MultiFloat::hasBinaryZero );
            
            MultiFloat_exposer.def( 
                "hasBinaryZero"
                , hasBinaryZero_function_value
                , bp::release_gil_policy()
                , "Return whether or not at least one of the elements of this vector\nis binary zero (the float is equal to 0x00000000)" );
        
        }
        { //::SireMaths::MultiFloat::isAligned
        
            typedef bool ( ::SireMaths::MultiFloat::*isAligned_function_type)(  ) const;
            isAligned_function_type isAligned_function_value( &::SireMaths::MultiFloat::isAligned );
            
            MultiFloat_exposer.def( 
                "isAligned"
                , isAligned_function_value
                , bp::release_gil_policy()
                , "Return whether or not this MultiFloat is correctly aligned. If it is not,\nthen any SSE operations will fail" );
        
        }
        { //::SireMaths::MultiFloat::isBinaryOne
        
            typedef bool ( ::SireMaths::MultiFloat::*isBinaryOne_function_type)(  ) const;
            isBinaryOne_function_type isBinaryOne_function_value( &::SireMaths::MultiFloat::isBinaryOne );
            
            MultiFloat_exposer.def( 
                "isBinaryOne"
                , isBinaryOne_function_value
                , bp::release_gil_policy()
                , "Return whether all of the elements of this MultiFloat are\nequal to 0xFFFFFFFF (e.g. every bit in the entire vector is 1)" );
        
        }
        { //::SireMaths::MultiFloat::isBinaryZero
        
            typedef bool ( ::SireMaths::MultiFloat::*isBinaryZero_function_type)(  ) const;
            isBinaryZero_function_type isBinaryZero_function_value( &::SireMaths::MultiFloat::isBinaryZero );
            
            MultiFloat_exposer.def( 
                "isBinaryZero"
                , isBinaryZero_function_value
                , bp::release_gil_policy()
                , "Return whether all of the elements of this MultiFloat are\nequal to 0x00000000 (e.g. every bit in the entire vector is 0)" );
        
        }
        { //::SireMaths::MultiFloat::isNotBinaryOne
        
            typedef bool ( ::SireMaths::MultiFloat::*isNotBinaryOne_function_type)(  ) const;
            isNotBinaryOne_function_type isNotBinaryOne_function_value( &::SireMaths::MultiFloat::isNotBinaryOne );
            
            MultiFloat_exposer.def( 
                "isNotBinaryOne"
                , isNotBinaryOne_function_value
                , bp::release_gil_policy()
                , "Return whether all of the elements of this MultiFloat are\nnot equal to 0xFFFFFFFF (e.g. at least one bit in the entire vector is 0)" );
        
        }
        { //::SireMaths::MultiFloat::isNotBinaryZero
        
            typedef bool ( ::SireMaths::MultiFloat::*isNotBinaryZero_function_type)(  ) const;
            isNotBinaryZero_function_type isNotBinaryZero_function_value( &::SireMaths::MultiFloat::isNotBinaryZero );
            
            MultiFloat_exposer.def( 
                "isNotBinaryZero"
                , isNotBinaryZero_function_value
                , bp::release_gil_policy()
                , "Return whether all of the elements of this MultiFloat are\nnot equal to 0x00000000 (e.g. at least one bit in the entire vector is 1)" );
        
        }
        { //::SireMaths::MultiFloat::logicalAnd
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*logicalAnd_function_type)( ::SireMaths::MultiFloat const & ) const;
            logicalAnd_function_type logicalAnd_function_value( &::SireMaths::MultiFloat::logicalAnd );
            
            MultiFloat_exposer.def( 
                "logicalAnd"
                , logicalAnd_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::logicalAnd
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*logicalAnd_function_type)( ::SireMaths::MultiUInt const & ) const;
            logicalAnd_function_type logicalAnd_function_value( &::SireMaths::MultiFloat::logicalAnd );
            
            MultiFloat_exposer.def( 
                "logicalAnd"
                , logicalAnd_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::logicalAnd
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*logicalAnd_function_type)( ::SireMaths::MultiInt const & ) const;
            logicalAnd_function_type logicalAnd_function_value( &::SireMaths::MultiFloat::logicalAnd );
            
            MultiFloat_exposer.def( 
                "logicalAnd"
                , logicalAnd_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::logicalAndNot
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*logicalAndNot_function_type)( ::SireMaths::MultiFloat const & ) const;
            logicalAndNot_function_type logicalAndNot_function_value( &::SireMaths::MultiFloat::logicalAndNot );
            
            MultiFloat_exposer.def( 
                "logicalAndNot"
                , logicalAndNot_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::logicalAndNot
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*logicalAndNot_function_type)( ::SireMaths::MultiInt const & ) const;
            logicalAndNot_function_type logicalAndNot_function_value( &::SireMaths::MultiFloat::logicalAndNot );
            
            MultiFloat_exposer.def( 
                "logicalAndNot"
                , logicalAndNot_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::logicalAndNot
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*logicalAndNot_function_type)( ::SireMaths::MultiUInt const & ) const;
            logicalAndNot_function_type logicalAndNot_function_value( &::SireMaths::MultiFloat::logicalAndNot );
            
            MultiFloat_exposer.def( 
                "logicalAndNot"
                , logicalAndNot_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::logicalNot
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*logicalNot_function_type)(  ) const;
            logicalNot_function_type logicalNot_function_value( &::SireMaths::MultiFloat::logicalNot );
            
            MultiFloat_exposer.def( 
                "logicalNot"
                , logicalNot_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::logicalOr
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*logicalOr_function_type)( ::SireMaths::MultiFloat const & ) const;
            logicalOr_function_type logicalOr_function_value( &::SireMaths::MultiFloat::logicalOr );
            
            MultiFloat_exposer.def( 
                "logicalOr"
                , logicalOr_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::logicalXor
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*logicalXor_function_type)( ::SireMaths::MultiFloat const & ) const;
            logicalXor_function_type logicalXor_function_value( &::SireMaths::MultiFloat::logicalXor );
            
            MultiFloat_exposer.def( 
                "logicalXor"
                , logicalXor_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::max
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*max_function_type)( ::SireMaths::MultiFloat const & ) const;
            max_function_type max_function_value( &::SireMaths::MultiFloat::max );
            
            MultiFloat_exposer.def( 
                "max"
                , max_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::max
        
            typedef float ( ::SireMaths::MultiFloat::*max_function_type)(  ) const;
            max_function_type max_function_value( &::SireMaths::MultiFloat::max );
            
            MultiFloat_exposer.def( 
                "max"
                , max_function_value
                , bp::release_gil_policy()
                , "Return the maximum value in the vector" );
        
        }
        { //::SireMaths::MultiFloat::min
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*min_function_type)( ::SireMaths::MultiFloat const & ) const;
            min_function_type min_function_value( &::SireMaths::MultiFloat::min );
            
            MultiFloat_exposer.def( 
                "min"
                , min_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::min
        
            typedef float ( ::SireMaths::MultiFloat::*min_function_type)(  ) const;
            min_function_type min_function_value( &::SireMaths::MultiFloat::min );
            
            MultiFloat_exposer.def( 
                "min"
                , min_function_value
                , bp::release_gil_policy()
                , "Return the minimum value in the vector" );
        
        }
        { //::SireMaths::MultiFloat::multiplyAdd
        
            typedef ::SireMaths::MultiFloat & ( ::SireMaths::MultiFloat::*multiplyAdd_function_type)( ::SireMaths::MultiFloat const &,::SireMaths::MultiFloat const & ) ;
            multiplyAdd_function_type multiplyAdd_function_value( &::SireMaths::MultiFloat::multiplyAdd );
            
            MultiFloat_exposer.def( 
                "multiplyAdd"
                , multiplyAdd_function_value
                , ( bp::arg("val0"), bp::arg("val1") )
                , bp::return_self< >()
                , "" );
        
        }
        MultiFloat_exposer.def( !bp::self );
        MultiFloat_exposer.def( bp::self != bp::self );
        MultiFloat_exposer.def( bp::self & bp::self );
        MultiFloat_exposer.def( bp::self * bp::self );
        MultiFloat_exposer.def( bp::self + bp::self );
        MultiFloat_exposer.def( -bp::self );
        MultiFloat_exposer.def( bp::self - bp::self );
        MultiFloat_exposer.def( bp::self / bp::self );
        MultiFloat_exposer.def( bp::self < bp::self );
        MultiFloat_exposer.def( bp::self <= bp::self );
        { //::SireMaths::MultiFloat::operator=
        
            typedef ::SireMaths::MultiFloat & ( ::SireMaths::MultiFloat::*assign_function_type)( ::SireMaths::MultiFloat const & ) ;
            assign_function_type assign_function_value( &::SireMaths::MultiFloat::operator= );
            
            MultiFloat_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::operator=
        
            typedef ::SireMaths::MultiFloat & ( ::SireMaths::MultiFloat::*assign_function_type)( ::SireMaths::MultiDouble const & ) ;
            assign_function_type assign_function_value( &::SireMaths::MultiFloat::operator= );
            
            MultiFloat_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::operator=
        
            typedef ::SireMaths::MultiFloat & ( ::SireMaths::MultiFloat::*assign_function_type)( ::SireMaths::MultiInt const & ) ;
            assign_function_type assign_function_value( &::SireMaths::MultiFloat::operator= );
            
            MultiFloat_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::operator=
        
            typedef ::SireMaths::MultiFloat & ( ::SireMaths::MultiFloat::*assign_function_type)( float ) ;
            assign_function_type assign_function_value( &::SireMaths::MultiFloat::operator= );
            
            MultiFloat_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        MultiFloat_exposer.def( bp::self == bp::self );
        MultiFloat_exposer.def( bp::self > bp::self );
        MultiFloat_exposer.def( bp::self >= bp::self );
        { //::SireMaths::MultiFloat::operator[]
        
            typedef float ( ::SireMaths::MultiFloat::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMaths::MultiFloat::operator[] );
            
            MultiFloat_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        MultiFloat_exposer.def( bp::self ^ bp::self );
        MultiFloat_exposer.def( bp::self | bp::self );
        { //::SireMaths::MultiFloat::quickSet
        
            typedef void ( ::SireMaths::MultiFloat::*quickSet_function_type)( int,float ) ;
            quickSet_function_type quickSet_function_value( &::SireMaths::MultiFloat::quickSet );
            
            MultiFloat_exposer.def( 
                "quickSet"
                , quickSet_function_value
                , ( bp::arg("i"), bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::reciprocal
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*reciprocal_function_type)(  ) const;
            reciprocal_function_type reciprocal_function_value( &::SireMaths::MultiFloat::reciprocal );
            
            MultiFloat_exposer.def( 
                "reciprocal"
                , reciprocal_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::reciprocal_approx
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*reciprocal_approx_function_type)(  ) const;
            reciprocal_approx_function_type reciprocal_approx_function_value( &::SireMaths::MultiFloat::reciprocal_approx );
            
            MultiFloat_exposer.def( 
                "reciprocal_approx"
                , reciprocal_approx_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::reciprocal_approx_nr
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*reciprocal_approx_nr_function_type)(  ) const;
            reciprocal_approx_nr_function_type reciprocal_approx_nr_function_value( &::SireMaths::MultiFloat::reciprocal_approx_nr );
            
            MultiFloat_exposer.def( 
                "reciprocal_approx_nr"
                , reciprocal_approx_nr_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::rotate
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*rotate_function_type)(  ) const;
            rotate_function_type rotate_function_value( &::SireMaths::MultiFloat::rotate );
            
            MultiFloat_exposer.def( 
                "rotate"
                , rotate_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::rsqrt
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*rsqrt_function_type)(  ) const;
            rsqrt_function_type rsqrt_function_value( &::SireMaths::MultiFloat::rsqrt );
            
            MultiFloat_exposer.def( 
                "rsqrt"
                , rsqrt_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::rsqrt_approx
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*rsqrt_approx_function_type)(  ) const;
            rsqrt_approx_function_type rsqrt_approx_function_value( &::SireMaths::MultiFloat::rsqrt_approx );
            
            MultiFloat_exposer.def( 
                "rsqrt_approx"
                , rsqrt_approx_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::rsqrt_approx_nr
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*rsqrt_approx_nr_function_type)(  ) const;
            rsqrt_approx_nr_function_type rsqrt_approx_nr_function_value( &::SireMaths::MultiFloat::rsqrt_approx_nr );
            
            MultiFloat_exposer.def( 
                "rsqrt_approx_nr"
                , rsqrt_approx_nr_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::set
        
            typedef void ( ::SireMaths::MultiFloat::*set_function_type)( int,float ) ;
            set_function_type set_function_value( &::SireMaths::MultiFloat::set );
            
            MultiFloat_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the ith value of the multifloat to value" );
        
        }
        { //::SireMaths::MultiFloat::size
        
            typedef int ( *size_function_type )(  );
            size_function_type size_function_value( &::SireMaths::MultiFloat::size );
            
            MultiFloat_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::sqrt
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*sqrt_function_type)(  ) const;
            sqrt_function_type sqrt_function_value( &::SireMaths::MultiFloat::sqrt );
            
            MultiFloat_exposer.def( 
                "sqrt"
                , sqrt_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::sqrt_approx
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*sqrt_approx_function_type)(  ) const;
            sqrt_approx_function_type sqrt_approx_function_value( &::SireMaths::MultiFloat::sqrt_approx );
            
            MultiFloat_exposer.def( 
                "sqrt_approx"
                , sqrt_approx_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::sqrt_approx_nr
        
            typedef ::SireMaths::MultiFloat ( ::SireMaths::MultiFloat::*sqrt_approx_nr_function_type)(  ) const;
            sqrt_approx_nr_function_type sqrt_approx_nr_function_value( &::SireMaths::MultiFloat::sqrt_approx_nr );
            
            MultiFloat_exposer.def( 
                "sqrt_approx_nr"
                , sqrt_approx_nr_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::sum
        
            typedef float ( ::SireMaths::MultiFloat::*sum_function_type)(  ) const;
            sum_function_type sum_function_value( &::SireMaths::MultiFloat::sum );
            
            MultiFloat_exposer.def( 
                "sum"
                , sum_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::swap
        
            typedef void ( *swap_function_type )( ::SireMaths::MultiFloat &,int,::SireMaths::MultiFloat &,int );
            swap_function_type swap_function_value( &::SireMaths::MultiFloat::swap );
            
            MultiFloat_exposer.def( 
                "swap"
                , swap_function_value
                , ( bp::arg("f0"), bp::arg("idx0"), bp::arg("f1"), bp::arg("idx1") )
                , bp::release_gil_policy()
                , "Swap the values of the value at index idx0 in f0 with the value at index idx in f1" );
        
        }
        { //::SireMaths::MultiFloat::toArray
        
            typedef ::QVector< float > ( *toArray_function_type )( ::QVector< SireMaths::MultiFloat > const & );
            toArray_function_type toArray_function_value( &::SireMaths::MultiFloat::toArray );
            
            MultiFloat_exposer.def( 
                "toArray"
                , toArray_function_value
                , ( bp::arg("array") )
                , bp::release_gil_policy()
                , "Return the passed MultiFloat converted back into a normal array" );
        
        }
        { //::SireMaths::MultiFloat::toBinaryString
        
            typedef ::QString ( ::SireMaths::MultiFloat::*toBinaryString_function_type)(  ) const;
            toBinaryString_function_type toBinaryString_function_value( &::SireMaths::MultiFloat::toBinaryString );
            
            MultiFloat_exposer.def( 
                "toBinaryString"
                , toBinaryString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::toDoubleArray
        
            typedef ::QVector< double > ( *toDoubleArray_function_type )( ::QVector< SireMaths::MultiFloat > const & );
            toDoubleArray_function_type toDoubleArray_function_value( &::SireMaths::MultiFloat::toDoubleArray );
            
            MultiFloat_exposer.def( 
                "toDoubleArray"
                , toDoubleArray_function_value
                , ( bp::arg("array") )
                , bp::release_gil_policy()
                , "Return the passed MultiFloat converted back into a normal array of doubles" );
        
        }
        { //::SireMaths::MultiFloat::toString
        
            typedef ::QString ( ::SireMaths::MultiFloat::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::MultiFloat::toString );
            
            MultiFloat_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::MultiFloat::typeName );
            
            MultiFloat_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::MultiFloat::what
        
            typedef char const * ( ::SireMaths::MultiFloat::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::MultiFloat::what );
            
            MultiFloat_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        MultiFloat_exposer.staticmethod( "count" );
        MultiFloat_exposer.staticmethod( "fromArray" );
        MultiFloat_exposer.staticmethod( "size" );
        MultiFloat_exposer.staticmethod( "swap" );
        MultiFloat_exposer.staticmethod( "toArray" );
        MultiFloat_exposer.staticmethod( "toDoubleArray" );
        MultiFloat_exposer.staticmethod( "typeName" );
        MultiFloat_exposer.def( "__copy__", &__copy__);
        MultiFloat_exposer.def( "__deepcopy__", &__copy__);
        MultiFloat_exposer.def( "clone", &__copy__);
        MultiFloat_exposer.def( "__str__", &__str__< ::SireMaths::MultiFloat > );
        MultiFloat_exposer.def( "__repr__", &__str__< ::SireMaths::MultiFloat > );
        MultiFloat_exposer.def( "__len__", &__len_size< ::SireMaths::MultiFloat > );
        MultiFloat_exposer.def( "__getitem__", &::SireMaths::MultiFloat::getitem );
    }

}
