/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#include "n4matrix.h"
#include "nmatrix.h"
#include "nvector.h"

#include "SireBase/array2d.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"
#include "SireMaths/errors.h"

using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<N4Matrix> r_n4matrix(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const N4Matrix &matrix)
{
    writeHeader(ds, r_n4matrix, 1);
    
    SharedDataStream sds(ds);
    
    sds << matrix.array << matrix.nbigrows << matrix.nbigcolumns 
        << matrix.nrows << matrix.ncolumns;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, N4Matrix &matrix)
{
    VersionID v = readHeader(ds, r_n4matrix);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> matrix.array >> matrix.nbigrows >> matrix.nbigcolumns 
            >> matrix.nrows >> matrix.nbigcolumns;
    }
    else
        throw version_error(v, "1", r_n4matrix, CODELOC);
        
    return ds;
}

/** Null constructor */
N4Matrix::N4Matrix() : nbigrows(0), nbigcolumns(0), nrows(0), ncolumns(0)
{}

/** Construct a matrix with 'nbigrows' big rows, 'nbigcolumns' big columns,
    'nrows' rows and 'ncolumns' columns. The values in the matrix are not initialised */
N4Matrix::N4Matrix(int nbr, int nbc, int nr, int nc)
         : nbigrows(nbr), nbigcolumns(nbc), nrows(nr), ncolumns(nc)
{
    if (nbr <= 0 or nbc <= 0 or nr <= 0 or nc <= 0)
    {
        nbigrows = 0;
        nbigcolumns = 0;
        nrows = 0;
        ncolumns = 0;
    }
    else
    {
        array = QVector<double>(nbigrows*nbigcolumns*nrows*ncolumns);
        array.squeeze();
    }
}

/** Construct a matrix with 'nbigrows' big rows, 'nbigcolumns' big columns,
    'nrows' rows and 'ncolumns' columns. The values in the matrix are
    initialised to be equal to 'initial_value' */
N4Matrix::N4Matrix(int nbr, int nbc, int nr, int nc, double initial_value)
         : nbigrows(nbr), nbigcolumns(nbc), nrows(nr), ncolumns(nc)
{
    if (nbr <= 0 or nbc <= 0 or nr <= 0 or nc <= 0)
    {
        nbigrows = 0;
        nbigcolumns = 0;
        nrows = 0;
        ncolumns = 0;
    }
    else
    {
        array = QVector<double>(nbigrows*nbigcolumns*nrows*ncolumns, initial_value);
        array.squeeze();
    }
}

/** Construct from the passed matrix - this creates a matrix
    of dimension [1, 1, matrix.nRows(), matrix.nColumns()] */
N4Matrix::N4Matrix(const NMatrix &matrix)
         : nbigrows(0), nbigcolumns(0), nrows(0), ncolumns(0) 
{
    if (matrix.nRows() != 0)
    {
        nbigrows = 1;
        nbigcolumns = 1;
    
        if (matrix.isTransposed())
        {
            NMatrix mat_t = matrix.transpose().fullTranspose();
        
            array = mat_t.memory();
            
            nrows = mat_t.nRows();
            ncolumns = mat_t.nColumns();
        }
        else
        {
            nrows = matrix.nRows();
            ncolumns = matrix.nColumns();
            array = matrix.memory();
        }
    }
}

/** Construct from the passed Array or Matricies */
N4Matrix::N4Matrix(const SireBase::Array2D<NMatrix> &array4d)
         : nbigrows(array4d.nRows()), nbigcolumns(array4d.nColumns())
{
    int sz = array4d.nRows() * array4d.nColumns();

    if (sz > 0)
    {
        //get the size of the largest 2D Matrix within the array
        int max_nrow = 0;
        int max_ncolumn = 0;
        
        const NMatrix *array4d_data = array4d.constData();
        
        for (uint i=0; i<array4d.nRows(); ++i)
        {
            for (uint j=0; j<array4d.nColumns(); ++j)
            {
                int idx = array4d.map(i,j);
            
                max_nrow = qMax(max_nrow, array4d_data[idx].nRows());
                max_ncolumn = qMax(max_ncolumn, array4d_data[idx].nColumns());
            }
        }

        sz *= (max_nrow*max_ncolumn);
    
        array = QVector<double>(sz, 0);
        array.squeeze();
        
        nbigrows = array4d.nRows();
        nbigcolumns = array4d.nColumns();
        nrows = max_nrow;
        ncolumns = max_ncolumn;
        
        double *data = array.data();
        
        for (int i=0; i<nbigrows; ++i)
        {
            for (int j=0; j<nbigcolumns; ++j)
            {
                int idx = array4d.map(i,j);

                const NMatrix &mat = array4d_data[idx];
                
                const double *mat_array = mat.constData();
                
                for (int k=0; k<mat.nRows(); ++k)
                {
                    for (int l=0; l<mat.nColumns(); ++l)
                    {
                        data[ offset(i,j,k,l) ] = mat_array[ mat.offset(k,l) ];
                    }
                }
            }
        }
    }
}

/** Construct from the passed vector of vector of vector of vectors... */
N4Matrix::N4Matrix(const QVector< QVector< QVector< QVector<double> > > > &array4d)
         : nbigrows(0), nbigcolumns(0), nrows(0), ncolumns(0)
{
    const QVector< QVector< QVector<double> > > *array4d_data = array4d.constData();
    
    nrows = array4d.count();
    ncolumns = 0;
    
    for (int i=0; i<nrows; ++i)
    {
        ncolumns = qMax(ncolumns, array4d_data[i].count());
    }
    
    if (ncolumns > 0)
    {
        Array2D<NMatrix> arrays(nrows, ncolumns);
        
        for (int i=0; i<nrows; ++i)
        {
            for (int j=0; j<array4d_data[i].count(); ++j)
            {
                arrays(i,j) = NMatrix(array4d_data[i][j]);
            }
        }

        this->operator=( N4Matrix(arrays) );
    }
}

/** Copy constructor */
N4Matrix::N4Matrix(const N4Matrix &other)
         : array(other.array), nbigrows(other.nbigrows),
           nbigcolumns(other.nbigcolumns), nrows(other.nrows),
           ncolumns(other.ncolumns)
{}

/** Destructor */
N4Matrix::~N4Matrix()
{}

const char* N4Matrix::typeName()
{
    return QMetaType::typeName( qMetaTypeId<N4Matrix>() );
}

const char* N4Matrix::what() const
{
    return N4Matrix::typeName();
}

/** Copy assignment operator */
N4Matrix& N4Matrix::operator=(const N4Matrix &other)
{
    if (this != &other)
    {
        array = other.array;
        nbigrows = other.nbigrows;
        nbigcolumns = other.nbigcolumns;
        nrows = other.nrows;
        ncolumns = other.ncolumns;
    }

    return *this;
}

/** Comparison operator */
bool N4Matrix::operator==(const N4Matrix &other) const
{
    return this == &other or
           (nbigrows == other.nbigrows and nbigcolumns == other.nbigcolumns and
            nrows == other.nrows and ncolumns == other.ncolumns and
            array == other.array);
}

/** Comparison operator */
bool N4Matrix::operator!=(const N4Matrix &other) const
{
    return not this->operator==(other);
}

/** Assert that the index [i,j,k,l] is valid for this matrix

    \throw SireError::invalid_index
*/
void N4Matrix::assertValidIndex(int i, int j, int k, int l) const
{
    if ( i < 0 or i >= nbigrows or j < 0 or j >= nbigcolumns or
         k < 0 or k >= nrows or l < 0 or l >= ncolumns )
    {
        throw SireError::invalid_index( QObject::tr(
            "The index [%1,%2,%3,%4] is invalid for the matrix with "
            "dimension [%5,%6,%7,%8].")
                .arg(i).arg(j).arg(k).arg(l)
                .arg(nbigrows).arg(nbigcolumns)
                .arg(nrows).arg(ncolumns), CODELOC );
    }
}

/** Return a modifiable reference to the value at [i,j,k,l]

    \throw SireError::invalid_index
*/
double& N4Matrix::operator()(int i, int j, int k, int l)
{
    this->assertValidIndex(i,j,k,l);
    
    return array.data()[ this->offset(i,j,k,l) ];
}

/** Return a reference to the value at [i,j,k,l]

    \throw SireError::invalid_index
*/
const double& N4Matrix::operator()(int i, int j, int k, int l) const
{
    this->assertValidIndex(i,j,k,l);
    
    return array.data()[ this->offset(i,j,k,l) ];
}

/** Return the sub-matrix view at [i,j] 

    \throw SireError::invalid_index
*/
NMatrix N4Matrix::operator()(int i, int j) const
{
    return this->view(i,j);
}

/** Assert that this matrix has 'nbigrows' big rows

    \throw SireError::incompatible_error
*/
void N4Matrix::assertNBigRows(int nr) const
{
    if (nr != nbigrows)
        throw SireError::incompatible_error( QObject::tr(
                "The number of big rows in this matrix (dimension "
                "[%1,%2,%3,%4]) is not "
                "equal to %5.")
                    .arg(nbigrows).arg(nbigcolumns)
                    .arg(nrows).arg(ncolumns).arg(nr), CODELOC );
}

/** Assert that this matrix has 'nbigcolumns' big columns

    \throw SireError::incompatible_error
*/
void N4Matrix::assertNBigColumns(int nc) const
{
    if (nc != ncolumns)
        throw SireError::incompatible_error( QObject::tr(
                "The number of big columns in this matrix (dimension "
                "[%1,%2,%3,%4]) is not "
                "equal to %5.")
                    .arg(nbigrows).arg(nbigcolumns)
                    .arg(nrows).arg(ncolumns).arg(nc), CODELOC );
}

/** Assert that this matrix has 'nrows' rows

    \throw SireError::incompatible_error
*/
void N4Matrix::assertNRows(int nr) const
{
    if (nr != nrows)
        throw SireError::incompatible_error( QObject::tr(
                "The number of rows in this matrix (dimension "
                "[%1,%2,%3,%4]) is not "
                "equal to %5.")
                    .arg(nbigrows).arg(nbigcolumns)
                    .arg(nrows).arg(ncolumns).arg(nr), CODELOC );
}

/** Assert that this matrix has 'ncolumns' columns

    \throw SireError::incompatible_error
*/
void N4Matrix::assertNColumns(int nc) const
{
    if (nc != ncolumns)
        throw SireError::incompatible_error( QObject::tr(
                "The number of columns in this matrix (dimension "
                "[%1,%2,%3,%4]) is not "
                "equal to %5.")
                    .arg(nbigrows).arg(nbigcolumns)
                    .arg(nrows).arg(ncolumns).arg(nc), CODELOC );
}

/** Matrix addition 

    \throw SireError::incompatible_error
*/
N4Matrix& N4Matrix::operator+=(const N4Matrix &other)
{
    assertNBigRows(other.nBigRows());
    assertNBigColumns(other.nBigColumns());
    assertNRows(other.nRows());
    assertNColumns(other.nColumns());
    
    double *data = array.data();
    const double *other_data = other.array.constData();
    const int sz = array.count();
    
    for (int i=0; i<sz; ++i)
    {
        data[i] += other_data[i];
    }
    
    return *this;
}

/** Matrix subtraction 

    \throw SireError::incompatible_error
*/
N4Matrix& N4Matrix::operator-=(const N4Matrix &other)
{
    assertNBigRows(other.nBigRows());
    assertNBigColumns(other.nBigColumns());
    assertNRows(other.nRows());
    assertNColumns(other.nColumns());
    
    double *data = array.data();
    const double *other_data = other.array.constData();
    const int sz = array.count();
    
    for (int i=0; i<sz; ++i)
    {
        data[i] -= other_data[i];
    }
    
    return *this;
}

/** Multiply all elements of this matrix by 'scale' */
N4Matrix& N4Matrix::operator*=(double scale)
{
    double *data = array.data();
    const int sz = array.count();

    if (scale == 0)
    {
        for (int i=0; i<sz; ++i)
        {
            data[i] = 0;
        }
    }
    else
    {
        for (int i=0; i<sz; ++i)
        {
            data[i] *= sz;
        }
    }
    
    return *this;
}

/** Divide all elements of this matrix by 'scale' */
N4Matrix& N4Matrix::operator/=(double scale)
{
    if (scale == 0)
        throw SireMaths::domain_error( QObject::tr(
            "This code does not support dividing a matrix by zero!"), CODELOC );
            
    return this->operator*=( 1/scale );
}

/** Return the negative of this matrix */
N4Matrix N4Matrix::operator-() const
{
    N4Matrix ret(*this);

    const int sz = array.count();
    double *ret_data = ret.array.data();
    
    for (int i=0; i<sz; ++i)
    {
        ret_data[i] = -ret_data[i];
    }
    
    return ret;
}

/** Matrix addition 

    \throw SireError::incompatible_error
*/
N4Matrix N4Matrix::operator+(const N4Matrix &other) const
{
    N4Matrix ret(*this);
    ret += other;
    
    return ret;
}

/** Matrix subtraction 

    \throw SireError::incompatible_error
*/
N4Matrix N4Matrix::operator-(const N4Matrix &other) const
{
    N4Matrix ret(*this);
    ret -= other;
    
    return ret;
}

/** Multiply all elements of this matrix by 'scale' */
N4Matrix N4Matrix::operator*(double scale) const
{
    N4Matrix ret(*this);
    
    ret *= scale;
    
    return ret;
}

/** Divide all elements of this matrix by 'scale' */
N4Matrix N4Matrix::operator/(double scale) const
{
    N4Matrix ret(*this);
    
    ret /= scale;
    
    return ret;
}

/** Return the number of big rows in this matrix */
int N4Matrix::nBigRows() const
{
    return nbigrows;
}

/** Return the number of big columns in this matrix */
int N4Matrix::nBigColumns() const
{
    return nbigcolumns;
}

/** Return the number of rows in this matrix */
int N4Matrix::nRows() const
{
    return nrows;
}

/** Return the number of columns in this matrix */
int N4Matrix::nColumns() const
{
    return ncolumns;
}

/** Redimension this matrix to have 'nbigrows' big rows, 
    'nbigcolumns' big columns, 'nrows' rows and 'ncolumns' 
    columns. The contents of this matrix are undefined after
    this redimension. This function will only reallocate
    memory if there is not enough memory allocated to store
    the new matrix. Use this function if you want to use
    the same piece of memory over and over again for lots
    of different size matricies - just create a matrix with
    the maximum dimension, then call this 'redimension' function
    whenever you want to change. It is very fast, as it just
    updates the internal record of the size of the matrix */
void N4Matrix::redimension(int nbr, int nbc, int nr, int nc)
{
    const int sz = nbr * nbc * nr * nc;
    
    if (sz <= 0)
    {
        nbigrows = 0;
        nbigcolumns = 0;
        nrows = 0;
        ncolumns = 0;
    }
    else
    {
        if (sz > array.count())
            array.resize(sz);
            
        nbigrows = nbr;
        nbigcolumns = nbc;
        nrows = nr;
        ncolumns = nc;
    }
}

/** Assert that there is an ith big row! 

    \throw SireError::invalid_index
*/
void N4Matrix::assertValidBigRow(int i) const
{
    if (i < 0 or i >= nbigrows)
        throw SireError::invalid_index( QObject::tr(
                "The matrix with dimension [%1,%2,%3,%4] does not have "
                "a big row with index %5.")
                    .arg(nbigrows).arg(nbigcolumns)
                    .arg(nrows).arg(ncolumns).arg(i), CODELOC );
}

/** Assert that there is an jth big column! 

    \throw SireError::invalid_index
*/
void N4Matrix::assertValidBigColumn(int j) const
{
    if (j < 0 or j >= nbigcolumns)
        throw SireError::invalid_index( QObject::tr(
                "The matrix with dimension [%1,%2,%3,%4] does not have "
                "a big column with index %5.")
                    .arg(nbigrows).arg(nbigcolumns)
                    .arg(nrows).arg(ncolumns).arg(j), CODELOC );
}

/** Assert that there is an kth row! 

    \throw SireError::invalid_index
*/
void N4Matrix::assertValidRow(int k) const
{
    if (k < 0 or k >= nrows)
        throw SireError::invalid_index( QObject::tr(
                "The matrix with dimension [%1,%2,%3,%4] does not have "
                "a row with index %5.")
                    .arg(nbigrows).arg(nbigcolumns)
                    .arg(nrows).arg(ncolumns).arg(k), CODELOC );
}

/** Assert that there is an lth column! 

    \throw SireError::invalid_index
*/
void N4Matrix::assertValidColumn(int l) const
{
    if (l < 0 or l >= ncolumns)
        throw SireError::invalid_index( QObject::tr(
                "The matrix with dimension [%1,%2,%3,%4] does not have "
                "a big row with index %5.")
                    .arg(nbigrows).arg(nbigcolumns)
                    .arg(nrows).arg(ncolumns).arg(l), CODELOC );
}

/** Return the sub-matrix view at [i,j,k,l] 

    \throw SireError::invalid_index
*/
NMatrix N4Matrix::view(int i, int j) const
{
    const double *ptr = array.constData() + checkedOffset(i,j,0,0);
    
    NMatrix mat(nrows, ncolumns);
    
    memcpy( mat.data(), ptr, nrows*ncolumns*sizeof(double) );
    
    return mat;
}

/** Set the view at [i,j] equal to 'matrix'

    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void N4Matrix::set(int i, int j, const NMatrix &matrix)
{
    assertNRows(matrix.nRows());
    assertNColumns(matrix.nColumns());
    
    double *ptr = array.data() + checkedOffset(i,j,0,0);
    
    if (matrix.isTransposed())
    {
        NMatrix mat_c = matrix.transpose().fullTranspose();
        
        memcpy( ptr, mat_c.constData(), nrows*ncolumns*sizeof(double) );
    }
    else
    {
        memcpy( ptr, matrix.constData(), nrows*ncolumns*sizeof(double) );
    }
}

/** Add the contents of 'matrix' to the sub-matrix view at [i,j]
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void N4Matrix::add(int i, int j, const NMatrix &matrix)
{
    assertNRows(matrix.nRows());
    assertNColumns(matrix.nColumns());
    
    assertValidIndex(i,j,0,0);
    
    if (matrix.isTransposed())
    {
        this->add( i, j, matrix.transpose().fullTranspose() );
        return;
    }
    
    double *ptr = array.data() + checkedOffset(i,j,0,0);
    
    const double *data = matrix.constData();
    const int sz = nrows * ncolumns;
    
    for (int i=0; i<sz; ++i)
    {
        ptr[i] += data[i];
    }
}

/** Subtract the contents of 'matrix' from the sub-matrix view at [i,j]
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void N4Matrix::subtract(int i, int j, const NMatrix &matrix)
{
    assertNRows(matrix.nRows());
    assertNColumns(matrix.nColumns());
    
    assertValidIndex(i,j,0,0);
    
    if (matrix.isTransposed())
    {
        this->subtract( i, j, matrix.transpose().fullTranspose() );
        return;
    }
    
    double *ptr = array.data() + checkedOffset(i,j,0,0);
    
    const double *data = matrix.constData();
    const int sz = nrows * ncolumns;
    
    for (int i=0; i<sz; ++i)
    {
        ptr[i] -= data[i];
    }
}

/** Set the value at [i,j,k,l] equal to 'value'

    \throw SireError::invalid_index
*/
void N4Matrix::set(int i, int j, int k, int l, double value)
{
    array[ checkedOffset(i,j,k,l) ] = value;
}

/** Set all entries in the matrix to the value 'value' */
void N4Matrix::setAll(double value)
{
    double *d = array.data();
    int sz = array.count();
    
    for (int i=0; i<sz; ++i)
    {
        d[i] = value;
    }
}

/** Return a raw pointer to the data of this matrix. The data is
    stored in column-major order (same as Fortran - not same as C++ or C).
    To be safe, use the 'offset' function to get the offset of 
    the value at [i,j] in this array */
double* N4Matrix::data()
{
    return array.data();
}

/** Return a raw pointer to the data of this matrix. The data is
    stored in column-major order (same as Fortran - not same as C++ or C).
    To be safe, use the 'offset' function to get the offset of 
    the value at [i,j] in this array */
const double* N4Matrix::data() const
{
    return array.constData();
}

/** Return a raw pointer to the data of this matrix. The data is
    stored in column-major order (same as Fortran - not same as C++ or C).
    To be safe, use the 'offset' function to get the offset of 
    the value at [i,j] in this array */
const double* N4Matrix::constData() const
{
    return array.constData();
}

/** Return the raw QVector memory used by this matrix */
QVector<double> N4Matrix::memory() const
{
    return array;
}

/** Calculate the offset in the 1D array of the value
    at index [i,j,k,l]
    
    \throw SireError::invalid_index
*/
int N4Matrix::checkedOffset(int i, int j, int k, int l) const
{
    this->assertValidIndex(i,j,k,l);
    return this->offset(i,j,k,l);
}

/** Return a string representation of this matrix */
QString N4Matrix::toString() const
{
    return QObject::tr( "N4Matrix( %1 by %2 by %3 by %4 )" )
                .arg(nbigrows).arg(nbigcolumns)
                .arg(nrows).arg(ncolumns);
}
