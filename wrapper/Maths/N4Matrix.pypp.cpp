// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "N4Matrix.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/array2d.hpp"

#include "SireError/errors.h"

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "n4matrix.h"

#include "nmatrix.h"

#include "nvector.h"

#include "n4matrix.h"

SireMaths::N4Matrix __copy__(const SireMaths::N4Matrix &other){ return SireMaths::N4Matrix(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_N4Matrix_class(){

    { //::SireMaths::N4Matrix
        typedef bp::class_< SireMaths::N4Matrix > N4Matrix_exposer_t;
        N4Matrix_exposer_t N4Matrix_exposer = N4Matrix_exposer_t( "N4Matrix", "This is a dense, double, general NMLK 4-dimensional matrix.\nThe data is stored as a column-major 2D matrix of column-major\n2D matricies (so each 2D sub-matrix is suitable for\nuse with Fortran BLAS or LAPACK functions). This is\ndesigned for high speed.\n\nThe data is implicitly shared (copy on write), so\ncopying a matrix is very fast.\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope N4Matrix_scope( N4Matrix_exposer );
        N4Matrix_exposer.def( bp::init< int, int, int, int >(( bp::arg("nbigrows"), bp::arg("nbigcolumns"), bp::arg("nrows"), bp::arg("columns") ), "Construct a matrix with nbigrows big rows, nbigcolumns big columns,\nnrows rows and ncolumns columns. The values in the matrix are not initialised") );
        N4Matrix_exposer.def( bp::init< int, int, int, int, double >(( bp::arg("nbigrows"), bp::arg("nbigcolumn"), bp::arg("nrows"), bp::arg("ncolumns"), bp::arg("initial_value") ), "Construct a matrix with nbigrows big rows, nbigcolumns big columns,\nnrows rows and ncolumns columns. The values in the matrix are\ninitialised to be equal to initial_value") );
        N4Matrix_exposer.def( bp::init< SireMaths::NMatrix const & >(( bp::arg("matrix") ), "Construct from the passed matrix - this creates a matrix\nof dimension [1, 1, matrix.nRows(), matrix.nColumns()]") );
        N4Matrix_exposer.def( bp::init< SireBase::Array2D< SireMaths::NMatrix > const & >(( bp::arg("matrix") ), "Construct from the passed Array or Matricies") );
        N4Matrix_exposer.def( bp::init< QVector< QVector< QVector< QVector< double > > > > const & >(( bp::arg("matrix") ), "Construct from the passed vector of vector of vector of vectors...") );
        N4Matrix_exposer.def( bp::init< SireMaths::N4Matrix const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMaths::N4Matrix::add
        
            typedef void ( ::SireMaths::N4Matrix::*add_function_type)( int,int,::SireMaths::NMatrix const & ) ;
            add_function_type add_function_value( &::SireMaths::N4Matrix::add );
            
            N4Matrix_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("matrix") )
                , bp::release_gil_policy()
                , "Add the contents of matrix to the sub-matrix view at [i,j]\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::N4Matrix::assertNBigColumns
        
            typedef void ( ::SireMaths::N4Matrix::*assertNBigColumns_function_type)( int ) const;
            assertNBigColumns_function_type assertNBigColumns_function_value( &::SireMaths::N4Matrix::assertNBigColumns );
            
            N4Matrix_exposer.def( 
                "assertNBigColumns"
                , assertNBigColumns_function_value
                , ( bp::arg("nbigcolumns") )
                , bp::release_gil_policy()
                , "Assert that this matrix has nbigcolumns big columns\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::N4Matrix::assertNBigRows
        
            typedef void ( ::SireMaths::N4Matrix::*assertNBigRows_function_type)( int ) const;
            assertNBigRows_function_type assertNBigRows_function_value( &::SireMaths::N4Matrix::assertNBigRows );
            
            N4Matrix_exposer.def( 
                "assertNBigRows"
                , assertNBigRows_function_value
                , ( bp::arg("nbigrows") )
                , bp::release_gil_policy()
                , "Assert that this matrix has nbigrows big rows\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::N4Matrix::assertNColumns
        
            typedef void ( ::SireMaths::N4Matrix::*assertNColumns_function_type)( int ) const;
            assertNColumns_function_type assertNColumns_function_value( &::SireMaths::N4Matrix::assertNColumns );
            
            N4Matrix_exposer.def( 
                "assertNColumns"
                , assertNColumns_function_value
                , ( bp::arg("ncolumns") )
                , bp::release_gil_policy()
                , "Assert that this matrix has ncolumns columns\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::N4Matrix::assertNRows
        
            typedef void ( ::SireMaths::N4Matrix::*assertNRows_function_type)( int ) const;
            assertNRows_function_type assertNRows_function_value( &::SireMaths::N4Matrix::assertNRows );
            
            N4Matrix_exposer.def( 
                "assertNRows"
                , assertNRows_function_value
                , ( bp::arg("nrows") )
                , bp::release_gil_policy()
                , "Assert that this matrix has nrows rows\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::N4Matrix::assertValidBigColumn
        
            typedef void ( ::SireMaths::N4Matrix::*assertValidBigColumn_function_type)( int ) const;
            assertValidBigColumn_function_type assertValidBigColumn_function_value( &::SireMaths::N4Matrix::assertValidBigColumn );
            
            N4Matrix_exposer.def( 
                "assertValidBigColumn"
                , assertValidBigColumn_function_value
                , ( bp::arg("j") )
                , bp::release_gil_policy()
                , "Assert that there is an jth big column\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::N4Matrix::assertValidBigRow
        
            typedef void ( ::SireMaths::N4Matrix::*assertValidBigRow_function_type)( int ) const;
            assertValidBigRow_function_type assertValidBigRow_function_value( &::SireMaths::N4Matrix::assertValidBigRow );
            
            N4Matrix_exposer.def( 
                "assertValidBigRow"
                , assertValidBigRow_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Assert that there is an ith big row\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::N4Matrix::assertValidColumn
        
            typedef void ( ::SireMaths::N4Matrix::*assertValidColumn_function_type)( int ) const;
            assertValidColumn_function_type assertValidColumn_function_value( &::SireMaths::N4Matrix::assertValidColumn );
            
            N4Matrix_exposer.def( 
                "assertValidColumn"
                , assertValidColumn_function_value
                , ( bp::arg("l") )
                , bp::release_gil_policy()
                , "Assert that there is an lth column\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::N4Matrix::assertValidIndex
        
            typedef void ( ::SireMaths::N4Matrix::*assertValidIndex_function_type)( int,int,int,int ) const;
            assertValidIndex_function_type assertValidIndex_function_value( &::SireMaths::N4Matrix::assertValidIndex );
            
            N4Matrix_exposer.def( 
                "assertValidIndex"
                , assertValidIndex_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("k"), bp::arg("l") )
                , bp::release_gil_policy()
                , "Assert that the index [i,j,k,l] is valid for this matrix\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::N4Matrix::assertValidRow
        
            typedef void ( ::SireMaths::N4Matrix::*assertValidRow_function_type)( int ) const;
            assertValidRow_function_type assertValidRow_function_value( &::SireMaths::N4Matrix::assertValidRow );
            
            N4Matrix_exposer.def( 
                "assertValidRow"
                , assertValidRow_function_value
                , ( bp::arg("k") )
                , bp::release_gil_policy()
                , "Assert that there is an kth row\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::N4Matrix::checkedOffset
        
            typedef int ( ::SireMaths::N4Matrix::*checkedOffset_function_type)( int,int,int,int ) const;
            checkedOffset_function_type checkedOffset_function_value( &::SireMaths::N4Matrix::checkedOffset );
            
            N4Matrix_exposer.def( 
                "checkedOffset"
                , checkedOffset_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("k"), bp::arg("l") )
                , bp::release_gil_policy()
                , "Calculate the offset in the 1D array of the value\nat index [i,j,k,l]\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::N4Matrix::memory
        
            typedef ::QVector< double > ( ::SireMaths::N4Matrix::*memory_function_type)(  ) const;
            memory_function_type memory_function_value( &::SireMaths::N4Matrix::memory );
            
            N4Matrix_exposer.def( 
                "memory"
                , memory_function_value
                , bp::release_gil_policy()
                , "Return the raw QVector memory used by this matrix" );
        
        }
        { //::SireMaths::N4Matrix::nBigColumns
        
            typedef int ( ::SireMaths::N4Matrix::*nBigColumns_function_type)(  ) const;
            nBigColumns_function_type nBigColumns_function_value( &::SireMaths::N4Matrix::nBigColumns );
            
            N4Matrix_exposer.def( 
                "nBigColumns"
                , nBigColumns_function_value
                , bp::release_gil_policy()
                , "Return the number of big columns in this matrix" );
        
        }
        { //::SireMaths::N4Matrix::nBigRows
        
            typedef int ( ::SireMaths::N4Matrix::*nBigRows_function_type)(  ) const;
            nBigRows_function_type nBigRows_function_value( &::SireMaths::N4Matrix::nBigRows );
            
            N4Matrix_exposer.def( 
                "nBigRows"
                , nBigRows_function_value
                , bp::release_gil_policy()
                , "Return the number of big rows in this matrix" );
        
        }
        { //::SireMaths::N4Matrix::nColumns
        
            typedef int ( ::SireMaths::N4Matrix::*nColumns_function_type)(  ) const;
            nColumns_function_type nColumns_function_value( &::SireMaths::N4Matrix::nColumns );
            
            N4Matrix_exposer.def( 
                "nColumns"
                , nColumns_function_value
                , bp::release_gil_policy()
                , "Return the number of columns in this matrix" );
        
        }
        { //::SireMaths::N4Matrix::nRows
        
            typedef int ( ::SireMaths::N4Matrix::*nRows_function_type)(  ) const;
            nRows_function_type nRows_function_value( &::SireMaths::N4Matrix::nRows );
            
            N4Matrix_exposer.def( 
                "nRows"
                , nRows_function_value
                , bp::release_gil_policy()
                , "Return the number of rows in this matrix" );
        
        }
        { //::SireMaths::N4Matrix::offset
        
            typedef int ( ::SireMaths::N4Matrix::*offset_function_type)( int,int,int,int ) const;
            offset_function_type offset_function_value( &::SireMaths::N4Matrix::offset );
            
            N4Matrix_exposer.def( 
                "offset"
                , offset_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("k"), bp::arg("l") )
                , bp::release_gil_policy()
                , "" );
        
        }
        N4Matrix_exposer.def( bp::self != bp::self );
        { //::SireMaths::N4Matrix::operator()
        
            typedef double const & ( ::SireMaths::N4Matrix::*__call___function_type)( int,int,int,int ) const;
            __call___function_type __call___function_value( &::SireMaths::N4Matrix::operator() );
            
            N4Matrix_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("k"), bp::arg("l") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMaths::N4Matrix::operator()
        
            typedef ::SireMaths::NMatrix ( ::SireMaths::N4Matrix::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireMaths::N4Matrix::operator() );
            
            N4Matrix_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , "" );
        
        }
        N4Matrix_exposer.def( bp::self * bp::other< double >() );
        N4Matrix_exposer.def( bp::self + bp::self );
        N4Matrix_exposer.def( -bp::self );
        N4Matrix_exposer.def( bp::self - bp::self );
        N4Matrix_exposer.def( bp::self / bp::other< double >() );
        { //::SireMaths::N4Matrix::operator=
        
            typedef ::SireMaths::N4Matrix & ( ::SireMaths::N4Matrix::*assign_function_type)( ::SireMaths::N4Matrix const & ) ;
            assign_function_type assign_function_value( &::SireMaths::N4Matrix::operator= );
            
            N4Matrix_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        N4Matrix_exposer.def( bp::self == bp::self );
        { //::SireMaths::N4Matrix::redimension
        
            typedef void ( ::SireMaths::N4Matrix::*redimension_function_type)( int,int,int,int ) ;
            redimension_function_type redimension_function_value( &::SireMaths::N4Matrix::redimension );
            
            N4Matrix_exposer.def( 
                "redimension"
                , redimension_function_value
                , ( bp::arg("nbigrows"), bp::arg("nbigcolumns"), bp::arg("nrows"), bp::arg("ncolumns") )
                , bp::release_gil_policy()
                , "Redimension this matrix to have nbigrows big rows,\nnbigcolumns big columns, nrows rows and ncolumns\ncolumns. The contents of this matrix are undefined after\nthis redimension. This function will only reallocate\nmemory if there is not enough memory allocated to store\nthe new matrix. Use this function if you want to use\nthe same piece of memory over and over again for lots\nof different size matricies - just create a matrix with\nthe maximum dimension, then call this redimension function\nwhenever you want to change. It is very fast, as it just\nupdates the internal record of the size of the matrix" );
        
        }
        { //::SireMaths::N4Matrix::set
        
            typedef void ( ::SireMaths::N4Matrix::*set_function_type)( int,int,::SireMaths::NMatrix const & ) ;
            set_function_type set_function_value( &::SireMaths::N4Matrix::set );
            
            N4Matrix_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("matrix") )
                , bp::release_gil_policy()
                , "Set the view at [i,j] equal to matrix\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::N4Matrix::set
        
            typedef void ( ::SireMaths::N4Matrix::*set_function_type)( int,int,int,int,double ) ;
            set_function_type set_function_value( &::SireMaths::N4Matrix::set );
            
            N4Matrix_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("k"), bp::arg("l"), bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the value at [i,j,k,l] equal to value\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::N4Matrix::setAll
        
            typedef void ( ::SireMaths::N4Matrix::*setAll_function_type)( double ) ;
            setAll_function_type setAll_function_value( &::SireMaths::N4Matrix::setAll );
            
            N4Matrix_exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Set all entries in the matrix to the value value" );
        
        }
        { //::SireMaths::N4Matrix::subtract
        
            typedef void ( ::SireMaths::N4Matrix::*subtract_function_type)( int,int,::SireMaths::NMatrix const & ) ;
            subtract_function_type subtract_function_value( &::SireMaths::N4Matrix::subtract );
            
            N4Matrix_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("matrix") )
                , bp::release_gil_policy()
                , "Subtract the contents of matrix from the sub-matrix view at [i,j]\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::N4Matrix::toString
        
            typedef ::QString ( ::SireMaths::N4Matrix::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::N4Matrix::toString );
            
            N4Matrix_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this matrix" );
        
        }
        { //::SireMaths::N4Matrix::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::N4Matrix::typeName );
            
            N4Matrix_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::N4Matrix::view
        
            typedef ::SireMaths::NMatrix ( ::SireMaths::N4Matrix::*view_function_type)( int,int ) const;
            view_function_type view_function_value( &::SireMaths::N4Matrix::view );
            
            N4Matrix_exposer.def( 
                "view"
                , view_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "Return the sub-matrix view at [i,j,k,l]\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMaths::N4Matrix::what
        
            typedef char const * ( ::SireMaths::N4Matrix::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::N4Matrix::what );
            
            N4Matrix_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        N4Matrix_exposer.staticmethod( "typeName" );
        N4Matrix_exposer.def( "__copy__", &__copy__);
        N4Matrix_exposer.def( "__deepcopy__", &__copy__);
        N4Matrix_exposer.def( "clone", &__copy__);
        N4Matrix_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::N4Matrix >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        N4Matrix_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::N4Matrix >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        N4Matrix_exposer.def_pickle(sire_pickle_suite< ::SireMaths::N4Matrix >());
        N4Matrix_exposer.def( "__str__", &__str__< ::SireMaths::N4Matrix > );
        N4Matrix_exposer.def( "__repr__", &__str__< ::SireMaths::N4Matrix > );
    }

}
