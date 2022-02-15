// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "NVector.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "nmatrix.h"

#include "nvector.h"

#include "vector.h"

#include "nvector.h"

SireMaths::NVector __copy__(const SireMaths::NVector &other){ return SireMaths::NVector(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_NVector_class(){

    { //::SireMaths::NVector
        typedef bp::class_< SireMaths::NVector > NVector_exposer_t;
        NVector_exposer_t NVector_exposer = NVector_exposer_t( "NVector", "This is a dense, double column Vector\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope NVector_scope( NVector_exposer );
        NVector_exposer.def( bp::init< int >(( bp::arg("nrows") ), "Construct a column vector with nrows rows - the values\nof the rows are not initialised") );
        NVector_exposer.def( bp::init< int, double >(( bp::arg("nrows"), bp::arg("initial_value") ), "Construct a column vector with nrows rows, all\ninitialised to the value initial_value") );
        NVector_exposer.def( bp::init< SireMaths::Vector const & >(( bp::arg("vector") ), "Construct from the passed 3-vector") );
        NVector_exposer.def( bp::init< QVector< double > const & >(( bp::arg("vector") ), "Construct from the passed vector") );
        NVector_exposer.def( bp::init< SireMaths::NVector const & >(( bp::arg("vector") ), "Construct from the passed 3-vector") );
        { //::SireMaths::NVector::assertNColumns
        
            typedef void ( ::SireMaths::NVector::*assertNColumns_function_type)( int ) const;
            assertNColumns_function_type assertNColumns_function_value( &::SireMaths::NVector::assertNColumns );
            
            NVector_exposer.def( 
                "assertNColumns"
                , assertNColumns_function_value
                , ( bp::arg("ncolumns") )
                , "Assert that this column vector has ncolumns columns - note\nthat this is a column vector, so only has 1 column\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::NVector::assertNRows
        
            typedef void ( ::SireMaths::NVector::*assertNRows_function_type)( int ) const;
            assertNRows_function_type assertNRows_function_value( &::SireMaths::NVector::assertNRows );
            
            NVector_exposer.def( 
                "assertNRows"
                , assertNRows_function_value
                , ( bp::arg("nrows") )
                , "Assert that this column vector has nrows rows\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::NVector::assertValidIndex
        
            typedef void ( ::SireMaths::NVector::*assertValidIndex_function_type)( int ) const;
            assertValidIndex_function_type assertValidIndex_function_value( &::SireMaths::NVector::assertValidIndex );
            
            NVector_exposer.def( 
                "assertValidIndex"
                , assertValidIndex_function_value
                , ( bp::arg("i") )
                , "Assert that the index i is valid\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::NVector::assertValidIndex
        
            typedef void ( ::SireMaths::NVector::*assertValidIndex_function_type)( int,int ) const;
            assertValidIndex_function_type assertValidIndex_function_value( &::SireMaths::NVector::assertValidIndex );
            
            NVector_exposer.def( 
                "assertValidIndex"
                , assertValidIndex_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , "Assert that the index [i,j] is valid\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::NVector::count
        
            typedef int ( ::SireMaths::NVector::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMaths::NVector::count );
            
            NVector_exposer.def( 
                "count"
                , count_function_value
                , "Return the number of values in this vector" );
        
        }
        { //::SireMaths::NVector::cross
        
            typedef ::SireMaths::NVector ( ::SireMaths::NVector::*cross_function_type)( ::SireMaths::NVector const & ) const;
            cross_function_type cross_function_value( &::SireMaths::NVector::cross );
            
            NVector_exposer.def( 
                "cross"
                , cross_function_value
                , ( bp::arg("other") )
                , "Return the cross product of this vector with other\nThrow: SireError::incompatible_error\nThrow: SireMaths::domain_error\n" );
        
        }
        { //::SireMaths::NVector::dot
        
            typedef double ( ::SireMaths::NVector::*dot_function_type)( ::SireMaths::NVector const & ) const;
            dot_function_type dot_function_value( &::SireMaths::NVector::dot );
            
            NVector_exposer.def( 
                "dot"
                , dot_function_value
                , ( bp::arg("other") )
                , "Return the dot product of this vector with other\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::NVector::isZero
        
            typedef bool ( ::SireMaths::NVector::*isZero_function_type)(  ) const;
            isZero_function_type isZero_function_value( &::SireMaths::NVector::isZero );
            
            NVector_exposer.def( 
                "isZero"
                , isZero_function_value
                , "Return whether or not this is a zero vector" );
        
        }
        { //::SireMaths::NVector::length
        
            typedef double ( ::SireMaths::NVector::*length_function_type)(  ) const;
            length_function_type length_function_value( &::SireMaths::NVector::length );
            
            NVector_exposer.def( 
                "length"
                , length_function_value
                , "Return the length of this vector" );
        
        }
        { //::SireMaths::NVector::length2
        
            typedef double ( ::SireMaths::NVector::*length2_function_type)(  ) const;
            length2_function_type length2_function_value( &::SireMaths::NVector::length2 );
            
            NVector_exposer.def( 
                "length2"
                , length2_function_value
                , "Return the length squared of this vector" );
        
        }
        { //::SireMaths::NVector::nColumns
        
            typedef int ( ::SireMaths::NVector::*nColumns_function_type)(  ) const;
            nColumns_function_type nColumns_function_value( &::SireMaths::NVector::nColumns );
            
            NVector_exposer.def( 
                "nColumns"
                , nColumns_function_value
                , "Return the number of columns in this column vector (0 or 1)" );
        
        }
        { //::SireMaths::NVector::nRows
        
            typedef int ( ::SireMaths::NVector::*nRows_function_type)(  ) const;
            nRows_function_type nRows_function_value( &::SireMaths::NVector::nRows );
            
            NVector_exposer.def( 
                "nRows"
                , nRows_function_value
                , "Return the number of rows in this column vector" );
        
        }
        { //::SireMaths::NVector::normalise
        
            typedef ::SireMaths::NVector ( ::SireMaths::NVector::*normalise_function_type)(  ) const;
            normalise_function_type normalise_function_value( &::SireMaths::NVector::normalise );
            
            NVector_exposer.def( 
                "normalise"
                , normalise_function_value
                , "Normalise this vector\nThrow: SireMaths::domain_error\n" );
        
        }
        NVector_exposer.def( bp::self != bp::self );
        { //::SireMaths::NVector::operator()
        
            typedef double const & ( ::SireMaths::NVector::*__call___function_type)( int ) const;
            __call___function_type __call___function_value( &::SireMaths::NVector::operator() );
            
            NVector_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMaths::NVector::operator()
        
            typedef double const & ( ::SireMaths::NVector::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireMaths::NVector::operator() );
            
            NVector_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        NVector_exposer.def( bp::self * bp::other< double >() );
        NVector_exposer.def( bp::self + bp::self );
        NVector_exposer.def( -bp::self );
        NVector_exposer.def( bp::self - bp::self );
        NVector_exposer.def( bp::self / bp::other< double >() );
        { //::SireMaths::NVector::operator=
        
            typedef ::SireMaths::NVector & ( ::SireMaths::NVector::*assign_function_type)( ::SireMaths::NVector const & ) ;
            assign_function_type assign_function_value( &::SireMaths::NVector::operator= );
            
            NVector_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        NVector_exposer.def( bp::self == bp::self );
        { //::SireMaths::NVector::operator[]
        
            typedef double & ( ::SireMaths::NVector::*__getitem___function_type)( int ) ;
            __getitem___function_type __getitem___function_value( &::SireMaths::NVector::operator[] );
            
            NVector_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_non_const_reference >()
                , "" );
        
        }
        { //::SireMaths::NVector::operator[]
        
            typedef double const & ( ::SireMaths::NVector::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMaths::NVector::operator[] );
            
            NVector_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMaths::NVector::set
        
            typedef void ( ::SireMaths::NVector::*set_function_type)( int,double ) ;
            set_function_type set_function_value( &::SireMaths::NVector::set );
            
            NVector_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("value") )
                , "Set the value at [i,0] to value\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::NVector::set
        
            typedef void ( ::SireMaths::NVector::*set_function_type)( int,int,double ) ;
            set_function_type set_function_value( &::SireMaths::NVector::set );
            
            NVector_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("value") )
                , "Set the value at [i,j] to value\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::NVector::setAll
        
            typedef void ( ::SireMaths::NVector::*setAll_function_type)( double ) ;
            setAll_function_type setAll_function_value( &::SireMaths::NVector::setAll );
            
            NVector_exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("value") )
                , "Set all values in this vector to value" );
        
        }
        { //::SireMaths::NVector::size
        
            typedef int ( ::SireMaths::NVector::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMaths::NVector::size );
            
            NVector_exposer.def( 
                "size"
                , size_function_value
                , "Return the number of values in this vector" );
        
        }
        { //::SireMaths::NVector::sum
        
            typedef double ( ::SireMaths::NVector::*sum_function_type)(  ) const;
            sum_function_type sum_function_value( &::SireMaths::NVector::sum );
            
            NVector_exposer.def( 
                "sum"
                , sum_function_value
                , "Return the sum of the elements of this vector" );
        
        }
        { //::SireMaths::NVector::toString
        
            typedef ::QString ( ::SireMaths::NVector::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::NVector::toString );
            
            NVector_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this vector" );
        
        }
        { //::SireMaths::NVector::transpose
        
            typedef ::SireMaths::NMatrix ( ::SireMaths::NVector::*transpose_function_type)(  ) const;
            transpose_function_type transpose_function_value( &::SireMaths::NVector::transpose );
            
            NVector_exposer.def( 
                "transpose"
                , transpose_function_value
                , "Return the transpose of this column vector (a row vector, which\nis implemented in this code as a NMatrix)" );
        
        }
        { //::SireMaths::NVector::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::NVector::typeName );
            
            NVector_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMaths::NVector::what
        
            typedef char const * ( ::SireMaths::NVector::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::NVector::what );
            
            NVector_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        NVector_exposer.staticmethod( "typeName" );
        NVector_exposer.def( bp::other< double >() * bp::self );
        NVector_exposer.def( "__copy__", &__copy__);
        NVector_exposer.def( "__deepcopy__", &__copy__);
        NVector_exposer.def( "clone", &__copy__);
        NVector_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::NVector >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NVector_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::NVector >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NVector_exposer.def_pickle(sire_pickle_suite< ::SireMaths::NVector >());
        NVector_exposer.def( "__str__", &__str__< ::SireMaths::NVector > );
        NVector_exposer.def( "__repr__", &__str__< ::SireMaths::NVector > );
        NVector_exposer.def( "__len__", &__len_size< ::SireMaths::NVector > );
    }

}
